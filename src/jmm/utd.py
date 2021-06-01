import numpy as np
import scipy.special

def F(x):
    '''Evaluate the Kouyoumjian transition function for (a vector of)
    argument(s) `x`. This is implemented in terms of the modified
    negative Fresnel integral. See Appendix B of "Introduction to
    Uniform Geometrical Theory of Diffraction" by McNamara, Pistorius,
    and Malherbe (or equation (4.72) on pp. 184 of the same
    reference).

    '''
    sqrt_x = np.sqrt(x)
    return 2*1j*sqrt_x*np.exp(1j*x)*scipy.special.fm(sqrt_x)

def N(beta, n, sign=1):
    if not all(abs(beta) < 2*np.pi):
        raise ValueError('beta should be in the range [0, 2pi)')
    if not all(0 <= n <= 2):
        raise ValueError('n should be in the range [0, 2]')
    if sign not in {1, -1}:
        raise ValueError('sign should be +1 or -1')
    ints = np.outer(np.ones_like(beta), [-1, 0, 1])
    return ints[np.argmin(abs(2*np.pi*n*I - beta - sign*np.pi), axis=1)]

def a(beta, n, sign=1):
    if sign not in {1, -1}:
        raise ValueError('sign should be +1 or -1')
    return 2*np.cos((2*n*np.pi*N(beta, n, sign) - beta)/2)**2

def Di(om, n, beta0, phi, phip, sign1=1, sign2=1):
    if om <= 0:
        raise ValueError('om should be positive')
    if not all(0 <= n <= 2):
        raise ValueError('n should be in the range [0, 2]')
    if sign1 not in {1, -1}:
        raise ValueError('sign1 should be +1 or -1')
    if sign2 not in {1, -1}:
        raise ValueError('sign2 should be +1 or -1')
    beta = phi + sign2*phip
    tmp1 = -np.exp(-1j*np.pi/4)/(2*n*np.sqrt(2*np.pi*om)*np.sin(beta0))
    tmp2 = np.cot((np.pi + sign1*beta)/(2*n))
    tmp3 = F(om*Li*a(beta, n, sign=sign1))
    return tmp1*tmp2*tmp3

def D(om, n, beta0, phi, phip, refl_coef):
    '''Compute the diffraction coefficient for a straight, sound-hard
    wedge with planar facets.

    Parameters
    ----------
    om : float [Hz]
        The frequency of the incident wave.
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
    refl_coef : float [dimensionless]
        The reflection coefficient at the point of diffraction.

    '''
    D1 = Di(om, n, beta0, phi, phip, sign1=1, sign2=-1)
    D2 = Di(om, n, beta0, phi, phip, sign1=-1, sign2=-1)
    D3 = Di(om, n, beta0, phi, phip, sign1=1, sign2=1)
    D4 = Di(om, n, beta0, phi, phip, sign1=-1, sign2=1)
    return D1 + D2 + refl_coef*(D3 + D4)

def D_from_geometry(om, alpha, no, e, s, sp, refl_coef):
    '''Compute the diffraction coefficient for a straight, sound-hard
    wedge with planar facets from a description of the local
    geometry. This computes the inputs to `D` and then evaluates it.

    Parameters
    ----------
    om : float [Hz]
        The frequency of the incident wave.
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

    '''
    if s@e < 0:
        raise RuntimeError('s and e point in opposite directions')
    if sp@e < 0:
        raise RuntimeError('sp and e point in opposite directions')

    n = 2 - alpha/np.pi
    beta0 = np.arccos(np.clip(np.dot(e, s), -1, 1))
    to = np.cross(no, e)

    st = s - (e@s)e
    st /= np.linalg.norm(st)
    phi = np.arctan2(to@st, n@st)

    spt = sp - (e@sp)e
    spt /= np.linalg.norm(spt)
    phip = np.arctan2(to@spt, n@spt)

    return D(om, n, beta0, phi, phi, refl_coef)
