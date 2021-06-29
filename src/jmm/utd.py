import logging
import numpy as np
import scipy.special

def _F(x):
    '''Evaluate the Kouyoumjian transition function for (a vector of)
    argument(s) `x`. This is implemented in terms of the modified
    negative Fresnel integral. See Appendix B of "Introduction to
    Uniform Geometrical Theory of Diffraction" by McNamara, Pistorius,
    and Malherbe (or equation (4.72) on pp. 184 of the same
    reference).

    '''
    try:
        sqrt_x = np.sqrt(x)
    except:
        import pdb; pdb.set_trace()
        pass
    fm = scipy.special.modfresnelm(sqrt_x)[0]
    return 2*1j*sqrt_x*np.exp(1j*x)*fm

def _N(beta, n, sign=1):
    ints = np.outer(np.ones_like(beta), [-1, 0, 1])
    abs_diff = abs(2*np.pi*n*ints - (beta + sign*np.pi).reshape(-1, 1))
    argmin = np.argmin(abs_diff, axis=1).reshape(-1, 1)
    return np.take_along_axis(ints, argmin, axis=1).ravel()

def _a(beta, n, sign=1):
    return 2*np.cos((2*n*np.pi*_N(beta, n, sign) - beta)/2)**2

def _Di(k, n, beta0, phi, phip, Li, sign1=1, sign2=1):
    if sign1 not in {1, -1}:
        raise ValueError('sign1 should be +1 or -1')
    if sign2 not in {1, -1}:
        raise ValueError('sign2 should be +1 or -1')

    beta = phi + sign2*phip

    # See the appendix of Tsingos and Funkhauser or (33) and (34) of
    # Kouyoumjian and Pathak for an explanation of this. This special
    # case handles the case when the argument of the cotangent below
    # is 0.
    N = _N(beta, n, sign1)
    eps = beta - sign1*(2*np.pi*n*N - np.pi)
    tol = np.finfo(np.float64).resolution
    mask = abs(eps) <= tol

    Di = np.empty(phi.size, dtype=np.complex128)

    if mask.sum() > 0:
        eps_masked = eps[mask]

        sgn = np.empty(mask.sum(), dtype=np.float64)
        sgn[eps_masked > 0] = 1
        sgn[eps_masked <= 0] = -1

        Di[mask] = np.sqrt(2*np.pi*k*Li[mask])*sgn
        Di[mask] -= 2*k*Li[mask]*eps_masked*np.exp(1j*np.pi/4)
        Di[mask] *= n*np.exp(1j*np.pi/4)

    # Now we compute the rest of the term the normal way
    tmp1 = -np.exp(-1j*np.pi/4)/(2*n*np.sqrt(2*np.pi*k)*np.sin(beta0[~mask]))
    tmp2 = 1/np.tan((np.pi + sign1*beta[~mask])/(2*n))
    tmp3 = _F(k*Li[~mask]*_a(beta[~mask], n, sign=sign1))
    Di[~mask] = tmp1*tmp2*tmp3

    return Di

def D_from_utd_params(k, n, beta0, phi, phip, Li, refl_coef=1):
    '''Compute the diffraction coefficient for a straight, sound-hard
    wedge with planar facets.

    Parameters
    ----------
    k : float [1/m]
        The wavenumber of the incident wave.
    n : float [dimensionless]
        Satisfies alpha = (2 - n)*pi, where alpha is the wedge angle.
    beta0 : float [rad]
        The angle that the incident and outgoing rays make with the
        diffracting edge.
    phi : float [rad]
        The angle that the plane spanned by the edge and the
        diffracted ray makes with the o-face.
    phip : float [rad]
        The angle that the plane spanned by the edge and the
        incident ray makes with the o-face.
    Li : float [dimensionless]
        A spreading coefficient used internally to calculate the
        diffraction coefficient.

    '''

    if k <= 0:
        raise ValueError('k should be positive')
    if not 0 <= n <= 2:
        raise ValueError('n should be in the range [0, 2]')

    D1 = _Di(k, n, beta0, phi, phip, Li, sign1=1, sign2=-1)
    D2 = _Di(k, n, beta0, phi, phip, Li, sign1=-1, sign2=-1)
    D3 = _Di(k, n, beta0, phi, phip, Li, sign1=1, sign2=1)
    D4 = _Di(k, n, beta0, phi, phip, Li, sign1=-1, sign2=1)

    return D1 + D2 + refl_coef*(D3 + D4)

def D_from_geometry(k, alpha, no, e, s, sp, t, hess, refl_coef=1):
    '''Compute the diffraction coefficient for a straight, sound-hard
    wedge with planar facets from a description of the local
    geometry. This computes the inputs to `D` and then evaluates it.

    Parameters
    ----------
    k : float [1/m]
        The wavenumber of the incident wave.
    alpha : float [rad]
        The angle of the wedge at the point of diffraction.
    no : array_like
        The surface normal for the "o-face" (the face on the far side
        of the point of incident from the incident ray).
    e : array_like
        The unit edge tangent vector, oriented so that `to =
        np.cross(no, e)` points into the "o-face".
    s : array_like
        The unit tangent vector of the diffracted ray at the point of
        diffraction.
    sp : array_like
        The unit tangent vector of the incident ray at the point of
        diffraction.
    t : array_like
        The distance from the point of diffraction.
    hess : array_like
        The Hessian of the eikonal at each point.

    '''

    # Since `s` will be undefined on the diffracting edge, we want to
    # mask these values out to avoid bad values while we compute
    # `D`. The main purpose of this is less confusing debugging.
    mask = np.logical_not(np.isnan(s).any(1))

    # We also want to mask out nodes where s or sp are exactly aligned
    # with e. This should only happen far away from the diffracting
    # edge where s or sp get aligned with e coincidentally; i.e., in a
    # disconnected component... (TODO: verify this...)
    mask[(s == e).all(1)] = False
    mask[(s == -e).all(1)] = False
    mask[(sp == e).all(1)] = False
    mask[(sp == -e).all(1)] = False

    s_masked = s[mask]
    sp_masked = sp[mask]
    hess_masked = hess[mask]
    t_masked = t[mask]

    if not np.isfinite(sp_masked.ravel()).all():
        raise RuntimeError('missing some entries of sp after masking')

    if not np.isfinite(hess_masked.ravel()).all():
        raise RuntimeError('missing some entries of hess after masking')

    if not np.isfinite(t_masked.ravel()).all():
        raise RuntimeError('missing some entries of t after masking')

    if (t_masked < 0).any():
        raise RuntimeError('passed negative distances (t)')

    n = 2 - alpha/np.pi
    beta0 = np.arccos(np.clip(s_masked@e, -1, 1))

    to = np.cross(no, e)
    to /= np.linalg.norm(to)

    # NOTE: when we compute phi and phip below, we just clamp them to
    # the desired range. There might be cases where these angles fall
    # outside [0, n*pi), but this should only happen as a result of
    # numerical error, so what we're doing here should be OK.

    tol = np.finfo(np.float64).resolution

    log = logging.getLogger('utd.py')

    st = s_masked - np.outer(s_masked@e, e)
    st /= np.sqrt(np.sum(st**2, axis=1)).reshape(-1, 1)
    phi = np.pi - (np.pi - np.arccos(st@to))*np.sign(st@no)
    if np.isnan(phi).any():
        raise RuntimeError('computed bad angles (phi)')

    if phi.min() < -tol:
        log.warn('out of range angle (phi: by %g)', abs(phi.min()))
    if 2*np.pi - alpha + tol < phi.max():
        log.warn('out of range angle (phi: by %g)',
                 abs(2*np.pi - alpha - phi.max()))
    phi = np.clip(phi, 0, n*np.pi)

    spt = sp_masked - np.outer(sp_masked@e, e)
    spt /= np.sqrt(np.sum(spt**2, axis=1)).reshape(-1, 1)
    phip = np.pi - (np.pi - np.arccos(-spt@to))*np.sign(-spt@no)
    if np.isnan(phip).any():
        raise RuntimeError('computed bad angles (phip)')

    if phip.min() < -tol:
        log.warn('out of range angle (phip: by %g)', abs(phip.min()))
    if 2*np.pi - alpha + tol < phip.max():
        log.warn('out of range angle (phip: by %g)',
                 abs(2*np.pi - alpha - phip.max()))
    phip = np.clip(phip, 0, n*np.pi)

    # Compute the radius of curvature in the "edge-fixed plane of
    # incidence"
    qe = e - (sp_masked.T*(sp_masked@e)).T
    qe /= np.sqrt(np.sum(qe**2, axis=1)).reshape(-1, 1)
    kappae = np.maximum(
        np.finfo(np.float64).resolution,
        np.einsum('ijk,ij,ik->i', hess_masked, qe, qe)
    )
    rhoe = 1/kappae

    q1 = no - (sp_masked.T*(sp_masked@no)).T
    q1 /= np.sqrt(np.sum(q1**2, axis=1)).reshape(-1, 1)
    kappa1 = np.maximum(
        np.finfo(np.float64).resolution,
        np.einsum('ijk,ij,ik->i', hess_masked, q1, q1)
    )
    rho1 = 1/kappa1

    q2 = np.cross(q1, to)
    q2 /= np.sqrt(np.sum(q2**2, axis=1)).reshape(-1, 1)
    kappa2 = np.maximum(
        np.finfo(np.float64).resolution,
        np.einsum('ijk,ij,ik->i', hess_masked, q2, q2)
    )
    rho2 = 1/kappa2

    Li = t_masked*(rhoe + t_masked)*rho1*rho2
    Li /= rhoe*(rho1 + t_masked)*(rho2 * t_masked)
    Li *= np.sin(beta0)**2

    if (Li < 0).any():
        raise RuntimeError('computed negative spreading factors (Li)')

    D = np.empty_like(t, dtype=np.complex128)
    D[~mask] = 1
    D[mask] = D_from_utd_params(k, n, beta0, phi, phip, Li, refl_coef)

    return D
