import datetime
import logging
import numpy as np
import os
import pickle

from abc import ABC
from cached_property import threaded_cached_property
from pytransform3d.rotations import matrix_from_axis_angle

import jmm.bb
import jmm.defs
import jmm.eik
import jmm.jet
import jmm.mesh
import jmm.slerp
import jmm.utd
import jmm.util
import jmm.utri

import colorcet as cc
import pyvista as pv
from jmm.plot import *

def _get_kappa1_and_kappa2(H):
    return np.linalg.svd(H, compute_uv=False)[1::-1]

def _get_rho1_and_rho2(H):
    kappa1, kappa2 = _get_kappa1_and_kappa2(H)
    rho1 = 0 if kappa1 == 0 else 1/kappa1
    rho2 = 0 if kappa2 == 0 else 1/kappa2
    return rho1, rho2

class Logger(object):
    def __init__(self):
        self.log = logging.getLogger(type(self).__name__)

class Domain(Logger):
    def __init__(self, verts, cells, refl_coef):
        super().__init__()

        # For now, just assume that the reflection coefficient is
        # provided as a single scalar. TODO: Later, we'll allow it to
        # be specified at each boundary vertex.
        self.refl_coef = refl_coef

        self.mesh = jmm.mesh.Mesh3(verts, cells)

        self.log.info('domain has %d vertices and %d cells',
                      self.num_verts, self.num_cells)

    def __reduce__(self):
        return (self.__class__, (self.verts, self.cells, self.refl_coef))

    @property
    def verts(self):
        return self.mesh.verts

    @property
    def num_verts(self):
        return self.mesh.num_verts

    @property
    def cells(self):
        return self.mesh.cells

    @property
    def num_cells(self):
        return self.mesh.num_cells

    @property
    def boundary_faces(self):
        return np.array(list(self.mesh.get_boundary_faces()))

    @property
    def h(self):
        return self.mesh.min_edge_length

class Field(ABC, Logger):
    speed_of_sound = 343
    minimum_magnitude = 1e-3 # == -60 dB

    def __init__(self, domain, index, omega, ftype, bd_inds, parent=None):
        super().__init__()

        self.domain = domain
        self.index = index
        self.omega = omega

        self._ftype = ftype

        self.bd_inds = bd_inds
        self.parent = parent

        orig = None
        if self.parent is not None:
            orig = self.parent.eik

        # TODO: instantiate these lazily to save memory!
        self.eik = jmm.eik.Eik3.from_mesh_and_ftype(domain.mesh, ftype, orig)

        self.solved = False
        self._scattered_fields = []

    def __reduce__(self):
        args = (self.domain, self.index, self.omega, self._ftype,
                self.bd_inds, self.parent)
        return (self.__class__, args)

    @property
    def parent_labels(self):
        s = '%s%d' % (type(self).__name__[0], self.index)
        if self.parent is not None:
            s += ' -> ' + self.parent.parent_labels
        return s

    @property
    def k(self):
        '''The wavenumber of the field (equal to `self.omega` divided by the
        speed of sound, obtained from the `Field` class).

        '''
        return self.omega/Field.speed_of_sound

    @property
    def h(self):
        return self.domain.h

    def solve(self):
        if self.solved:
            return

        self.solved = True

        jmm.util.tic()
        self.eik.solve()
        self.log.info('solved eikonal equation [%1.2fs]', jmm.util.toc())

    @property
    def is_solved(self):
        return self.eik.is_solved

    def _get_reflection_BCs(self, faces):
        max_mag_bd = np.nanmax(abs(self.amplitude[np.unique(faces)]))
        if max_mag_bd < Field.minimum_magnitude:
            return

        # Get the surface normal and reflection matrix for the
        # reflector
        nu = self.domain.mesh.get_face_normal(*faces[0])
        refl = np.eye(nu.size) - 2*np.outer(nu, nu)

        # Get a basis for the tangent plane of the reflector
        Q = np.linalg.svd(np.eye(nu.size) - np.outer(nu, nu))[0][:, :2]

        # Set up transform matrices for computing the reflected Hessian
        X_in = np.empty((3, 3))
        X_in[:, 1:] = Q
        X_out = np.empty((3, 3))
        X_out[:, 1:] = Q

        # Traverse the faces and pull out the BCs for the reflection
        bd_faces, bd_T, bd_grad_T, bd_hess_T, bd_t_in, bd_amplitude = \
            [], [], [], [], [], []
        for lf in faces:
            t_in = self.eik.grad_T[lf]

            t_in_is_nan = np.isnan(t_in).any(1)
            if sum(t_in_is_nan) > 1:
                continue

            # Get the eikonal gradients at the face vertices. Skip
            # this face if any of the gradients graze the surface
            # or seem to be emitted from the surface. If we have
            # t_in leaving a face, something unphysical is
            # happening.
            if (np.dot(t_in, nu) > -self.domain.mesh.eps).any():
                continue

            # Reflect the ray directions over the reflector to get the
            # BCs for the eikonal gradient for the reflection
            t_out = np.dot(t_in, refl)

            # Reflect the Hessian
            hess = np.empty((3, 3, 3))
            for i, l in enumerate(lf):
                X_in[:, 0] = t_in[i]
                X_out[:, 0] = t_out[i]
                Y = X_in@np.linalg.inv(X_out)
                hess[i] = Y.T@self.eik.hess[lf[i]]@Y

            bd_faces.append(lf)
            bd_T.append(self.eik.T[lf])
            bd_grad_T.append(t_out)
            bd_hess_T.append(hess)
            bd_t_in.append(t_in)
            bd_amplitude.append(self.amplitude[lf])

        if not bd_faces or len(bd_faces) == 1:
            return

        if np.nanmax(abs(np.array(bd_amplitude))) < Field.minimum_magnitude:
            return

        return np.array(bd_faces), np.array(bd_T), np.array(bd_grad_T), \
            np.array(bd_hess_T), np.array(bd_t_in), np.array(bd_amplitude)

    def _impute_diffracted_bd_grad_T(self, le, t_in):
        '''When an edge-diffracted field is generated by another
        edge-diffracted field which is incident upon it, there will be
        one node on the new edge which will be incident on the old
        edge. This node will be missing gradient data, which *is*
        correct, but we need to construct the Hermite polynomial
        approximating T on `le`. To do this, we take `t_in` from the
        originating field, and rotate it around the original
        diffracting edge until it's aligned with the new edge. We can
        then use this `grad_T` to compute the directional derivative
        that we need for the Hermite polynomial approximating `T`.

        '''
        # Order the indices of `le` so that `l0` indexes the node
        # that's on the diffracting edge (i.e., is missing gradient
        # data), and `l1` indexes the node that *isn't* on the
        # diffracting edge (i.e., isn't).
        l0, l1 = le[np.argsort(np.isfinite(self.eik.grad_T[le]).all(1))]

        # Find the diffracting edge which `l0` is incident on. We'll
        # use this as the rotation axis.
        inc_diff_edges = [
            le for le in self.domain.mesh.get_inc_diff_edges(l0)
            if all(self.eik.has_BCs(l) for l in le)
        ]
        if not len(inc_diff_edges) == 1:
            raise RuntimeError('entered weird state while imputing grad_T')
        le_diff = inc_diff_edges[0]
        assert any(self.eik.grad_T[_].all() for _ in le_diff)

        # Compute the rotation axis from the edge endpoints.
        q = -np.subtract(*self.domain.mesh.verts[le_diff])
        q /= np.linalg.norm(q)

        # Get the `t_in` vector at `l0` and make sure it's finite
        t_in_0 = self.eik.t_in[l0]
        if not np.isfinite(t_in_0).all():
            raise RuntimeError('entered weird state while imputing grad_T')

        # Compute the unit vector pointing from `l0` to `l1` (i.e., in
        # the direction of the new diffracting edge
        u = self.domain.mesh.verts[l1] - self.domain.mesh.verts[l0]
        u /= np.linalg.norm(u)

        # Project `l0`'s `t_in` vector into the orthogonal complement
        # of the rotation axis so that we can compute the amount that
        # we need to rotate the `t_in` vector to bring it into
        # alignment with the new diffracting edge
        v = np.dot(np.eye(3) - np.outer(q, q), t_in_0)
        v_norm = np.linalg.norm(v)
        if v_norm == 0:
            raise RuntimeError('bad t_in vector while imputing...')
        v /= v_norm

        # Rotate the `t_in` vector about `q` until it's aligned with
        # the edge `le`
        theta = np.arccos(np.dot(u, v))
        rot = matrix_from_axis_angle([*q, -theta])
        grad_T_imputed = rot@t_in_0

        return grad_T_imputed

    def _get_diffraction_BCs(self, edges):
        max_mag_bd = np.nanmax(abs(self.amplitude[np.unique(edges)]))
        if max_mag_bd < Field.minimum_magnitude:
            return

        # Get the tangent vector for the diffracting edge
        e = -np.subtract(*self.domain.verts[edges[0]])
        e /= np.linalg.norm(e)

        bd_edges, bd_T_bb, bd_rho1, bd_t_in, bd_amplitude = \
            [], [], [], [], []

        # Traverse the faces and pull out the BCs for the
        # edge-diffracted field
        for le in edges:
            T = self.eik.T[le]

            # Get the eikonal gradients (the `t_in` vectors) at the
            # edge, and skip this edge if any of them graze the edge
            grad_T = self.eik.grad_T[le]
            mask = np.isnan(grad_T).any(1)
            if (abs(grad_T[~mask]@e - 1) < self.domain.mesh.eps).any():
                continue

            num_finite = np.sum(~mask)
            if num_finite == 0:
                raise RuntimeError('bad t_in vectors')
            elif num_finite == 1:
                grad_T[mask] = self._impute_diffracted_bd_grad_T(le, grad_T)

            # Compute the directional derivative along `le` from
            # `grad_T` and then compute the cubic Hermite polynomial
            # approximating T along the edge
            dx = -np.subtract(*self.domain.verts[le])
            DT = grad_T@dx
            x = np.array([0, np.linalg.norm(dx)])
            T_bb = jmm.bb.Bb31.from_1d_data(T, DT, x)

            t_in = self.eik.grad_T[le]
            t_in[mask] = np.nan

            rho1 = []
            bad_rho1 = False
            for i, t in enumerate(t_in):
                if np.isnan(t).any():
                    rho1.append(0)
                else:
                    q1 = e - (t@e)*t
                    if np.linalg.norm(q1) < 1e-14:
                        bad_rho1 = True
                        break
                    hess = self.eik.hess[le[i]]
                    rho1.append(1/np.einsum('ij,i,j', hess, q1, q1)
                                if np.isfinite(hess).all() else 0)
            if bad_rho1:
                continue

            bd_edges.append(le)
            bd_T_bb.append(T_bb)
            bd_rho1.append(rho1)
            bd_t_in.append(t_in)
            bd_amplitude.append(self.amplitude[le])

        if not bd_edges:
            return

        if np.nanmax(abs(np.array(bd_amplitude))) < Field.minimum_magnitude:
            return

        return np.array(bd_edges), np.array(bd_T_bb), np.array(bd_rho1), \
            np.array(bd_t_in), np.array(bd_amplitude)

    def _init_scattered_fields(self):
        '''Set up the scattered fields. To do this, we find each reflector and
        diffractor in `self.domain` incident on `self.mask`. If `self`
        is a `ReflectedField` or `DiffractedField`, we also skip the
        reflector or diffractor corresponding to *this* field. This
        makes sense for a constant speed of sound but is too
        restrictive for a varying speed of sound, so this will need to
        be fixed later.

        '''
        if not self.solved:
            self.solve()

        mesh, fields = self.domain.mesh, []

        for index, faces in enumerate(mesh.reflectors):
            if isinstance(self, ReflectedField) and index == self.index:
                continue

            # Don't reflect from reflectors which contain the parent
            # diffractor
            if isinstance(self, DiffractedField) and \
               np.in1d(np.unique(self.bd_inds), np.unique(faces)).all():
                continue

            BCs = self._get_reflection_BCs(faces)
            if not BCs:
                continue

            fields.append(
                ReflectedField(
                    self.domain, index, self.omega, *BCs, parent=self))

        for index, edges in enumerate(mesh.diffractors):
            if isinstance(self, DiffractedField) and index == self.index:
                continue

            # Don't diffract from edges which are contained in the
            # parent reflector
            if isinstance(self, ReflectedField) and \
               np.in1d(np.unique(edges), np.unique(self.bd_inds)).all():
                continue

            BCs = self._get_diffraction_BCs(edges)
            if not BCs:
                continue
            fields.append(
                DiffractedField(
                    self.domain, index, self.omega, *BCs, parent=self))

        self._scattered_fields = fields

    @property
    def scattered_fields(self):
        '''Iterate over the fields of all orders scattered by this field,
        returning them in nondecreasing order of their minimum eikonal
        value.

        '''

        skip = dict()
        skip_tol = np.finfo(self.dtype).resolution

        fields = [self]

        while fields:
            field = fields.pop(0)

            lmin, tmin = field.time_argmin_and_min
            key = (type(field), field.index, lmin)
            if key in skip and abs(skip[key] - tmin) < skip_tol:
                continue

            skip[key] = tmin

            # DEBUG: dump the field before solving it for debugging
            pickled_field_path = 'field.pickle'
            if os.path.exists(pickled_field_path):
                os.remove(pickled_field_path)
            with open(pickled_field_path, 'wb') as f:
                pickle.dump(field, f)
                self.log.info('saved pickled field to %s', pickled_field_path)

            # Solve the next scattered field...
            field.solve()

            # ... and return it as the next iterate
            yield field

            # Initialize the next field's scattered fields and merge
            # them into `fields`
            if not field._scattered_fields:
                field._init_scattered_fields()
            fields.extend(field._scattered_fields)
            fields = sorted(fields, key=lambda _: _.eik.T.min())

    @property
    def dtype(self):
        return np.float64

    @property
    def _boundary_mask(self):
        mask = np.zeros(self.domain.num_verts, dtype=np.bool_)
        mask[np.unique(self.bd_inds)] = True
        return mask

    @property
    def r(self):
        return self.eik.T

    @property
    def time(self):
        time = np.empty(self.domain.num_verts)
        time[...] = self.r[...]/Field.speed_of_sound
        if not np.isfinite(time).all():
            raise RuntimeError('found bad eikonal values while computing time')
        return time

    @property
    def time_argmin_and_min(self):
        time = np.empty(self.domain.num_verts)
        try:
            time[...] = self.r[...]/Field.speed_of_sound
        except:
            import pdb; pdb.set_trace()
            pass
        l = np.nanargmin(time)
        if not np.isfinite(time[l]):
            raise RuntimeError("minimum time isn't finite")
        return l, time[l]

    @property
    def direction(self):
        direction = np.empty((self.domain.num_verts, 3))
        direction[...] = self.eik.grad_T[...]
        if not np.isfinite(direction[~self._boundary_mask]).all():
            raise RuntimeError('found bad arrival directions')
        return direction

    @property
    def _scale_tol(self):
        return self.h/2

    @property
    def _scale(self):
        scale = np.ones(self.domain.num_verts)

        dot = np.sum(self.eik.t_out*self.eik.grad_T, axis=1)
        dot = np.clip(dot, -1, 1)

        arc_length = np.arccos(dot)

        # We compute the scale so that it equals 1 if the arc length
        # is less than h, and otherwise apply a Gaussian window so
        # that as arc_length increases from h to 2*h, it drops 60 dB
        # (i.e. scale decays from 1 to 1e-3)

        log = np.log(Field.minimum_magnitude)
        log *= ((arc_length - self._scale_tol)/self._scale_tol)**2

        min_log = np.log(np.finfo(np.float64).eps)
        mask = log < min_log

        scale[mask] = 0
        scale[~mask] = np.exp(log[~mask])
        scale[self._boundary_mask | (arc_length < self._scale_tol)] = 1

        return scale

    @property
    def magnitude_mask(self):
        return abs(self.amplitude) >= Field.minimum_magnitude

    def transport_scalars(self, values, skip_filled=False):
        # Verify that transporting is going to complete successfully
        for l in range(self.domain.num_verts):
            if self.eik.has_par(l):
                par = self.eik.get_par(l)
                for l_ in par.l:
                    if self.eik.has_BCs(l_) and not np.isfinite(values[l_]):
                        raise ValueError('bad value will propagate')

        return self.eik.transport_scalars(values, skip_filled)

    def transport_curvature(self, kappa, skip_filled=False):
        # TODO: error handling...

        return self.eik.transport_curvature(kappa, skip_filled)

    def _fill_in_missing_amplitude_in(self, amplitude_in):
        # TODO: filling the amplitude in when a boundary is incident
        # on a diffracting edge is really challenging... probably can
        # only handel this correctly if we work with the geometric
        # spreading instead. so, for now, just set the NAN values to 0
        # so we can transport. anyway, as h -> 0 this should do
        # something somewhat reasonable.

        nan_bd_inds = np.unique(self.bd_inds)
        nan_bd_inds = nan_bd_inds[np.isnan(amplitude_in[nan_bd_inds])]

        amplitude_in[nan_bd_inds] = 0 # :-(

        return amplitude_in

class PointSourceField(Field):
    def __init__(self, domain, src_index, omega):
        index = src_index # re-use index of point source for field index
        ftype = jmm.defs.Ftype.PointSource
        bd_inds = np.array([src_index])
        super().__init__(domain, index, omega, ftype, bd_inds)

        jet = jmm.jet.Jet3.make_point_source()
        self.eik.add_pt_src_BCs(src_index, jet)

    def __reduce__(self):
        args = (self.domain, self.bd_inds[0], self.omega)
        return (self.__class__, args)

    @threaded_cached_property
    def amplitude(self):
        if not self.is_solved:
            raise RuntimeError('solve before accessing amplitude')

        mask = np.zeros(self.domain.num_verts, dtype=np.bool_)
        mask[self.bd_inds] = True

        A = np.empty(self.domain.num_verts, dtype=np.complex128)
        A[mask] = np.nan
        A[~mask] = self._scale[~mask]/self.r[~mask]

        return A

class ReflectedField(Field):
    def __init__(self, domain, index, omega, bd_faces, bd_T,
                 bd_grad_T, bd_hess_T, bd_t_in, bd_amplitude,
                 parent=None):
        num_faces = bd_faces.shape[0]
        if bd_T.shape[0] != num_faces or bd_grad_T.shape[0] != num_faces:
            raise ValueError('boundary faces and BCs must have the same shape')

        if bd_faces.shape[0] == 0:
            raise ValueError('no boundary faces were passed!')

        ftype = jmm.defs.Ftype.Reflection

        super().__init__(domain, index, omega, ftype, bd_faces, parent=parent)

        self.bd_T = bd_T
        self.bd_grad_T = bd_grad_T
        self.bd_hess = bd_hess_T
        self.bd_t_in = bd_t_in
        self.bd_amplitude = bd_amplitude

        for lf, T, grad_T, hess_T, t_in in zip(
                bd_faces, bd_T, bd_grad_T, bd_hess_T, bd_t_in):
            jets = [jmm.jet.Jet3(t, *dt) for t, dt in zip(T, grad_T)]
            self.eik.add_refl_BCs(*lf, *jets, hess_T, t_in)

    def __reduce__(self):
        args = (self.domain, self.index, self.omega, self.bd_inds,
                self.bd_T, self.bd_grad_T, self.bd_hess, self.bd_t_in,
                self.bd_amplitude, self.parent)
        return (self.__class__, args)

    @threaded_cached_property
    def reflector_face_normal(self):
        return self.domain.mesh.get_face_normal(*self.bd_inds[0])

    @property
    def _scale(self):
        nu = self.reflector_face_normal
        refl = np.eye(3) - 2*np.outer(nu, nu)

        t_in_refl = self.eik.t_in@refl
        t_out = self.eik.t_out
        dot = np.clip(np.sum(t_in_refl*t_out, axis=1), -1, 1)
        arc_length = np.arccos(dot)

        log = np.log(Field.minimum_magnitude)
        log *= ((arc_length - self._scale_tol)/self._scale_tol)**2

        min_log = np.log(np.finfo(np.float64).eps)
        mask = log < min_log

        scale = np.empty(self.domain.num_verts)
        scale[mask] = 0
        scale[~mask] = np.exp(log[~mask])
        scale[self._boundary_mask | (arc_length < self._scale_tol)] = 1

        return scale*super()._scale

    def _fill_in_missing_kappa1_in_and_kappa2_in(self, kappa1_in, kappa2_in,
                                                 nan_bd_inds):
        for l0 in range(self.domain.num_verts):
            if not self.eik.has_par(l0):
                continue
            l, b = self.eik.get_par(l0).get_active()
            l_nan = [_ for _ in l if _ in nan_bd_inds]
            if not l_nan:
                continue
            if len(l_nan) != 1:
                import pdb; pdb.set_trace()
                pass
            l_nan = l_nan[0]
            le = [_ for _ in self.domain.mesh.get_inc_diff_edges(l_nan)
                  if all(self.parent.eik.has_BCs(m) for m in _)]
            if len(le) != 1:
                import pdb; pdb.set_trace()
                pass
            le = le[0]
            xb = b@self.domain.verts[l]
            utri = jmm.utri.Utri.from_eik_without_l(self.parent.eik, xb, le)
            H = utri.approx_hess()
            k1, k2 = _get_kappa1_and_kappa2(H)
            if not np.isfinite(k1) or not np.isfinite(k2):
                raise RuntimeError('bad curvature!')
            kappa1_in[l0] = k1
            kappa2_in[l0] = k2

        return kappa1_in, kappa2_in

    @threaded_cached_property
    def amplitude(self):
        if not self.is_solved:
            raise RuntimeError('solve before accessing amplitude')

        refl_coef = self.domain.refl_coef # TODO: make this varying

        rho1, rho2 = np.array([_get_rho1_and_rho2(H) for H in self.eik.hess]).T

        num_verts = self.domain.num_verts

        # We initially compute and transport kappa1 in case any of the
        # curvatures are zero (can't take linear combinations in that
        # case).

        kappa1_in = np.empty(num_verts)
        kappa1_in[...] = np.nan

        kappa2_in = np.empty(num_verts)
        kappa2_in[...] = np.nan

        amplitude_in = np.empty(num_verts, dtype=self.bd_amplitude.dtype)
        amplitude_in[...] = np.nan

        nan_bd_inds = set()

        for lf, hess, amplitude in zip(
                self.bd_inds, self.bd_hess, self.bd_amplitude):
            k1f, k2f = np.array([_get_kappa1_and_kappa2(_) for _ in hess]).T
            for l, k1, k2, A in zip(lf, k1f, k2f, amplitude):
                kappa1_in[l] = k1
                kappa2_in[l] = k2
                amplitude_in[l] = A

                # Keep track of the boundary inds that have a NaN
                # curvature value. This will speed up "filling in"
                # curvature values in the next step
                if np.isnan(k1) or np.isnan(k2):
                    nan_bd_inds.add(l)

        if self.parent is None and nan_bd_inds:
            raise RuntimeError('entered bad state...')

        if nan_bd_inds:
            kappa1_in, kappa2_in = \
                self._fill_in_missing_kappa1_in_and_kappa2_in(
                    kappa1_in, kappa2_in, nan_bd_inds)
            amplitude_in = \
                self._fill_in_missing_amplitude_in(amplitude_in)

        kappa1_in = self.transport_curvature(kappa1_in)
        kappa2_in = self.transport_curvature(kappa2_in)
        amplitude_in = self.transport_scalars(amplitude_in)

        rho1_in = np.empty_like(kappa1_in)
        rho1_in[kappa1_in == 0] = np.inf
        rho1_in[kappa1_in != 0] = 1/kappa1_in[kappa1_in != 0]

        rho2_in = np.empty_like(kappa2_in)
        rho2_in[kappa2_in == 0] = np.inf
        rho2_in[kappa2_in != 0] = 1/kappa2_in[kappa2_in != 0]

        numer = rho1_in*rho2_in
        denom = rho1*rho2
        if (np.isinf(numer) & np.isinf(denom)).any():
            raise RuntimeError('found inconsistent spreading factor')

        mask = denom == 0

        sqrt_J = np.sqrt(numer[~mask]/denom[~mask])

        amplitude = np.empty_like(amplitude_in)
        amplitude[mask] = np.nan
        amplitude[~mask] = \
            amplitude_in[~mask]*refl_coef*sqrt_J*self._scale[~mask]

        return amplitude

class DiffractedField(Field):
    def __init__(self, domain, index, omega, bd_edges, bd_T_bb,
                 bd_rho1, bd_t_in, bd_amplitude, parent=None):
        num_edges = bd_edges.shape[0]
        if num_edges == 0:
            raise ValueError('no boundary edges were passed!')

        if bd_T_bb.shape[0] != num_edges or \
           bd_rho1.shape[0] != num_edges or \
           bd_t_in.shape[0] != num_edges or \
           bd_amplitude.shape[0] != num_edges:
            raise ValueError('BCs and edges must be compatible')

        ftype = jmm.defs.Ftype.EdgeDiffraction
        super().__init__(domain, index, omega, ftype, bd_edges, parent=parent)

        self.bd_T_bb = bd_T_bb
        self.bd_rho1 = bd_rho1
        self.bd_t_in = bd_t_in
        self.bd_amplitude = bd_amplitude

        for le, T_bb, rho1, t_in in zip(bd_edges, bd_T_bb, bd_rho1, bd_t_in):
            self.eik.add_diff_edge_BCs(le, T_bb, rho1, t_in)

    def __reduce__(self):
        args = (self.domain, self.index, self.omega, self.bd_inds,
                self.bd_T_bb, self.bd_rho1, self.bd_t_in,
                self.bd_amplitude, self.parent)
        return (self.__class__, args)

    @property
    def diffractor_tangent_vector(self):
        '''The tangent vector of the diffracting edge'''
        x0, x1 = self.domain.mesh.verts[self.bd_inds[0]]
        e = x1 - x0
        return e/np.linalg.norm(e)

    @property
    def _o_face_normal(self):
        l0, l1 = self.bd_inds[0]

        x0 = self.domain.mesh.verts[l0]
        x1 = self.domain.mesh.verts[l1]

        e = x1 - x0
        e /= np.linalg.norm(e)

        faces = self.domain.boundary_faces
        mask = (faces == l0).any(1) & (faces == l1).any(1)
        lf1, lf2 = faces[mask]

        l2 = next(_ for _ in lf1 if _ != l0 and _ != l1)
        x2 = self.domain.mesh.verts[l2]
        t = x2 - x0
        t /= np.linalg.norm(x2 - x0)
        no = np.cross(e, t)
        no /= np.linalg.norm(no)

        n = self.domain.mesh.get_face_normal(*lf1)
        if n@no > 0:
            return no

        l2 = next(_ for _ in lf2 if _ != l0 and _ != l1)
        x2 = self.domain.mesh.verts[l2]
        t = x2 - x0
        t /= np.linalg.norm(x2 - x0)
        no = np.cross(e, t)
        no /= np.linalg.norm(no)

        n = self.domain.mesh.get_face_normal(*lf2)
        if n@no > 0:
            return no

        raise RuntimeError('entered a weird state')

    @property
    def _scale(self):
        t = self.diffractor_tangent_vector

        arc_length_in = np.arccos(np.clip(self.eik.t_in@t, -1, 1))
        arc_length_out = np.arccos(np.clip(self.eik.t_out@t, -1, 1))
        arc_length_diff = arc_length_in - arc_length_out

        log = np.log(Field.minimum_magnitude)
        log *= ((arc_length_diff - self._scale_tol)/self._scale_tol)**2

        min_log = np.log(np.finfo(np.float64).eps)
        mask = log < min_log

        scale = np.empty(self.domain.num_verts)
        scale[mask] = 0
        scale[~mask] = np.exp(log[~mask])
        scale[self._boundary_mask | (arc_length_diff < self._scale_tol)] = 1

        return scale*super()._scale

    @property
    def _wedge_angle(self):
        return self.domain.mesh.get_edge_ext_angle(*self.bd_inds[0])

    def _fill_in_missing_kappa1_in(self, kappa1_in):
        '''Fill in missing kappa1_in values here: for nodes that have a single
        parent with kappa1 == inf, we compute kappa1 = 1/r, where r is
        the update distance (we treat the parent as a point source)

        '''
        for l in range(self.domain.num_verts):
            if not self.eik.has_par(l):
                continue
            par = self.eik.get_par(l)
            if par.num_active != 1:
                continue
            l0 = par.get_active()[0][0]
            if np.isinf(kappa1_in[l0]):
                dx = self.domain.verts[l] - self.domain.verts[l0]
                kappa1_in[l] = 1/np.linalg.norm(dx)
        return kappa1_in

    @threaded_cached_property
    def amplitude(self):
        if not self.is_solved:
            raise RuntimeError('solve before accessing amplitude')

        refl_coef = self.domain.refl_coef  # TODO: make this varying

        rho1, rho2 = np.array([_get_rho1_and_rho2(H) for H in self.eik.hess]).T
        if (rho1 < 0).any():
            raise RuntimeError('computed negative curvatures (rho1_in)')

        if (rho1 < 0).any():
            raise RuntimeError('computed negative curvatures (rho1)')

        num_verts = self.domain.num_verts


        # Transport the curvature instead of the radius of curvature
        # so that we can linearly interpolate when the curvature is
        # zero

        kappa1_in = np.empty(num_verts, dtype=np.float64)
        kappa1_in[...] = np.nan

        amplitude_in = np.empty(num_verts, dtype=np.complex128)
        amplitude_in[...] = np.nan

        for le, bd_rho1, bd_amplitude in zip(
                self.bd_inds, self.bd_rho1, self.bd_amplitude):
            for l, r1, A in zip(le, bd_rho1, bd_amplitude):
                kappa1_in[l] = 1/r1 if r1 != 0 else np.inf
                amplitude_in[l] = A

        kappa1_in = self._fill_in_missing_kappa1_in(kappa1_in)
        amplitude_in = self._fill_in_missing_amplitude_in(amplitude_in)

        kappa1_in = self.transport_curvature(kappa1_in, skip_filled=True)
        amplitude_in = self.transport_scalars(amplitude_in)

        rho1_in = np.empty_like(kappa1_in)
        rho1_in[kappa1_in == 0] = np.inf
        rho1_in[kappa1_in != 0] = 1/kappa1_in[kappa1_in != 0]

        numer = rho1_in
        denom = rho1*rho2
        if (np.isinf(numer) & np.isinf(denom)).any():
            raise RuntimeError('found inconsistent spreading factor')

        mask = denom != 0

        sqrt_J = np.sqrt(numer[mask]/denom[mask])

        bd_inds = np.unique(self.bd_inds)

        T_in = np.empty(num_verts, dtype=np.float64)
        T_in[bd_inds] = self.eik.T[bd_inds]
        T_in = self.transport_scalars(T_in)

        D = jmm.utd.D_from_geometry(self.k, self._wedge_angle,
                                    self._o_face_normal,
                                    self.diffractor_tangent_vector,
                                    self.eik.t_out, self.eik.t_in,
                                    T_in, self.eik.hess, refl_coef)

        amplitude = np.empty_like(amplitude_in)
        amplitude[~mask] = np.nan
        amplitude[mask] = amplitude_in[mask]*D[mask]*sqrt_J*self._scale[mask]

        return amplitude

class MultipleArrivals(Logger):
    def __init__(self, domain, root_field, num_arrivals):
        super().__init__()

        self.domain = domain
        self.root_field = root_field
        self.num_arrivals = num_arrivals

        num_verts = self.domain.num_verts

        shape = (num_verts, num_arrivals + 1)

        self._time = np.empty(shape, dtype=np.float64)
        self._time[...] = np.inf

        self._amplitude = np.empty(shape, dtype=np.complex128)
        self._amplitude[...] = np.nan

        self._direction = np.empty((*shape, 3), dtype=np.float64)
        self._direction[...] = np.nan

        self._fields = []

    def traverse(self):
        self.log.info('traversing scattered fields')

        num_accepted = 0

        for field in self.root_field.scattered_fields:
            self._fields.append(field)

            time = field.time

            min_time = time.min()
            min_time_str = str(datetime.timedelta(seconds=min_time))

            max_mag_dB = np.maximum(
                np.finfo(np.float64).eps, np.nanmax(abs(field.amplitude)))
            max_mag_dB = 20*np.log10(max_mag_dB)

            perc_above_thresh = \
                100*field.magnitude_mask.sum()/self.domain.num_verts

            self.log.info('accepted %s (%d fields)', type(field).__name__,
                          num_accepted + 1)
            self.log.info('    labels: %s', field.parent_labels)
            self.log.info('    earliest arrival: %s', min_time_str)
            self.log.info('    max |A|: %f dB', max_mag_dB)
            self.log.info('    above amp thresh: %1.2f%%', perc_above_thresh)

            if not np.isfinite(min_time):
                raise RuntimeError('bad field')

            self._time[:, -1] = time
            self._time[~field.magnitude_mask, -1] = np.inf

            self._amplitude[:, -1] = field.amplitude
            self._direction[:, -1, :] = field.direction

            perm = np.argsort(self._time, axis=1)
            self._time = np.take_along_axis(self._time, perm, axis=1)
            self._amplitude = np.take_along_axis(self._amplitude, perm, axis=1)
            for i in range(self._direction.shape[-1]):
                self._direction[:, :, i] = \
                    np.take_along_axis(self._direction[:, :, i], perm, axis=1)

            finite_time = np.isfinite(self.time)

            perc_finite = 100*finite_time.sum()/np.prod(finite_time.shape)
            self.log.info('computed %1.1f%% of all arrivals', perc_finite)

            if finite_time.all():
                break

            num_accepted += 1

    @property
    def time(self):
        return self._time[:, :-1]

    @property
    def amplitude(self):
        return self._amplitude[:, :-1]

    @property
    def direction(self):
        return self._direction[:, :-1, :]
