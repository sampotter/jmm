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
import jmm.util

import colorcet as cc
import pyvista as pv
from jmm.plot import *

class Logger(object):
    def __init__(self):
        self.log = logging.getLogger(type(self).__name__)

class Domain(Logger):
    def __init__(self, extended_verts, extended_cells, num_dom_verts,
                 num_dom_cells):
        super().__init__()

        verts = extended_verts[:num_dom_verts]
        cells = extended_cells[:num_dom_cells]
        self.mesh = jmm.mesh.Mesh3(verts, cells)

        self.log.info('domain has %d vertices and %d cells',
                      num_dom_verts, num_dom_cells)

        # Create extended mesh. When we do this, we need to insert the
        # boundary faces and diffracting edges into the new mesh to
        # ensure consistency between the computed solutions.
        self.extended_mesh = jmm.mesh.Mesh3(extended_verts, extended_cells)
        for le in self.mesh.get_diff_edges():
            self.extended_mesh.set_boundary_edge(*le, True)
        for lf in self.mesh.get_boundary_faces():
            self.extended_mesh.set_boundary_face(*lf, True)

        self.log.info('extended domain has %d vertices and %d cells',
                      extended_verts.shape[0], extended_cells.shape[0])

        self.h = self.extended_mesh.min_edge_length

    def __reduce__(self):
        extended_verts = self.extended_mesh.verts
        extended_cells = self.extended_mesh.cells
        num_dom_verts = self.mesh.num_verts
        num_dom_cells = self.mesh.num_cells
        args = (extended_verts, extended_cells, num_dom_verts, num_dom_cells)
        return (self.__class__, args)

    @property
    def num_verts(self):
        return self.mesh.num_verts

    def get_active_reflectors(self, mask):
        for faces in self.mesh.reflectors:
            active_faces = faces[mask[faces].any(1)]
            if active_faces.size > 0:
                yield active_faces

    def get_active_diffractors(self, mask):
        for edges in self.mesh.diffractors:
            active_edges = edges[mask[edges].any(1)]
            if active_edges.size > 0:
                yield active_edges

class Field(ABC, Logger):
    speed_of_sound = 343
    num_mask_threshold_bins = 513

    def __init__(self, domain, index, ftype, bd_inds, bd_T, bd_grad_T,
                 parent=None):
        super().__init__()

        self.domain = domain
        self.index = index
        self.bd_inds = bd_inds
        self.bd_T = bd_T
        self.bd_grad_T = bd_grad_T
        self.parent = parent

        # TODO: instantiate these lazily to save memory!
        self.eik = jmm.eik.Eik3.from_mesh_and_ftype(self.domain.mesh, ftype)
        self.extended_eik = jmm.eik.Eik3.from_mesh_and_ftype(
            self.domain.extended_mesh, ftype)

        self.solved = False
        self._scattered_fields = []

    def __reduce__(self):
        args = (self.domain, self.bd_inds, self.bd_T, self.bd_grad_T)
        kwargs = {'parent': self.parent}
        return (self.__class__, args, kwargs)

    @property
    def h(self):
        return self.domain.h

    def solve(self):
        if self.solved:
            return
        self.solved = True

        jmm.util.tic()
        self.eik.solve()
        self.log.info(
            'solved eikonal equation on domain [%1.2fs]',
            jmm.util.toc())

        jmm.util.tic()
        self.extended_eik.solve()
        self.log.info(
            'solved eikonal equation on extended domain [%1.2fs]',
            jmm.util.toc())

    def _get_reflection_BCs(self, faces):
        # Get the surface normal and reflection matrix for the
        # reflector
        nu = self.domain.mesh.get_face_normal(*faces[0])
        refl = np.eye(nu.size) - 2*np.outer(nu, nu)

        # Traverse the faces and pull out the BCs for the reflection
        bd_faces, bd_T, bd_grad_T, bd_t_in = [], [], [], []
        for lf in faces:
            # Get the eikonal gradients at the face vertices. Skip
            # this face if any of the gradients graze the surface or
            # seem to be emitted from the surface. If we have t_in
            # leaving a face, something unphysical is happening.
            t_in = self.eik.grad_T[lf]
            if (np.dot(t_in, nu) > -self.domain.mesh.eps).all():
                continue

            # Reflect the ray directions over the reflector to get the
            # BCs for the eikonal gradient for the reflection
            t_out = np.dot(t_in, refl)

            bd_faces.append(lf)
            bd_T.append(self.eik.T[lf])
            bd_grad_T.append(t_out)
            bd_t_in.append(t_in)

        if not bd_faces:
            return

        return np.array(bd_faces), np.array(bd_T), np.array(bd_grad_T), \
            np.array(bd_t_in)

    def _get_diffraction_BCs(self, edges):
        # Get the tangent vector for the diffracting edge
        t_e = -np.subtract(*self.domain.mesh.verts[edges[0]])
        t_e /= np.linalg.norm(t_e)

        # Traverse the faces and pull out the BCs for the
        # edge-diffracted field
        bd_edges, bd_T, bd_grad_T = [], [], []
        for le in edges:
            # Get the eikonal gradients (the `t_in` vectors) at the
            # edge, and skip this edge if any of them graze the edge
            t_in = self.eik.grad_T[le]
            if (abs(np.dot(t_in, t_e) - 1) < self.domain.mesh.eps).any():
                continue

            num_finite = np.isfinite(t_in).all(1).sum()
            if num_finite == 0:
                raise RuntimeError('bad t_in vectors')
            if num_finite == 1 and not isinstance(self, DiffractedField):
                raise RuntimeError('bad t_in vectors and should be diff field')
            if num_finite == 1:
                t_in_imputed = self._impute_diffracted_t_in(le, t_in)
                t_in[np.isnan(t_in).all(1)] = t_in_imputed

            bd_edges.append(le)
            bd_T.append(self.eik.T[le])
            bd_grad_T.append(t_in)

        if not bd_edges:
            return

        return np.array(bd_edges), np.array(bd_T), np.array(bd_grad_T)

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

        fields = []

        active_reflectors = self.domain.get_active_reflectors(self.mask)
        for index, active_faces in enumerate(active_reflectors):
            if isinstance(self, ReflectedField) and index == self.index:
                continue

            BCs = self._get_reflection_BCs(active_faces)
            if not BCs:
                continue

            fields.append(ReflectedField(self.domain, index, *BCs, parent=self))

        active_diffractors = self.domain.get_active_diffractors(self.mask)
        for index, active_edges in enumerate(active_diffractors):
            if isinstance(self, DiffractedField) and index == self.index:
                continue

            BCs = self._get_diffraction_BCs(active_edges)
            if not BCs:
                continue

            fields.append(DiffractedField(self.domain, index, *BCs, parent=self))

        self._scattered_fields = fields

    @property
    def scattered_fields(self):
        '''Iterate over the fields of all orders scattered by this field,
        returning them in nondecreasing order of their minimum eikonal
        value.

        '''

        fields = [self]

        while fields:
            field = fields.pop(0)

            # DEBUG: dump the field before solving it for debugging
            pickled_field_path = 'field.pickle'
            if os.path.exists(pickled_field_path):
                os.remove(pickled_field_path)
            with open(pickled_field_path, 'wb') as f:
                pickle.dump(field, f)

            # Solve the next scattered field and return it as the next iterate
            field.solve()
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
    def _error_thresh(self):
        return np.maximum(np.finfo(self.dtype).resolution, self.h**4)

    @threaded_cached_property
    def _phi_shadow(self):
        T = self.eik.T
        extended_T = self.extended_eik.T[:self.domain.num_verts]

        abs_diff = abs(T - extended_T)

        # A good number of points will agree to machine precision. We
        # want to mask those out first.
        #
        # TODO: we could do this more robustly by doing three-way Otsu
        # here instead of using this manual mask
        error_mask = abs_diff > self._error_thresh

        # Get the base-10 exponent of the absolute domain errors
        log_abs_diff = np.log10(np.maximum(1e-16, abs_diff[error_mask]))

        finite_mask = np.isfinite(log_abs_diff)

        # Compute a threshold in order to separate the distribution of
        # exponents using Otsu's method
        # p = jmm.util.otsu(
        #     log_abs_diff[finite_mask], bins=Field.num_mask_threshold_bins)
        # self._shadow_mask_thresh = 10**p
        self._shadow_mask_thresh = self.h**2

        # Evaluate the level set function for the shadow mask
        return T - extended_T - self._shadow_mask_thresh

    @threaded_cached_property
    def _phi(self):
        return np.maximum(self._phi_shadow, self._phi_valid_angle)

    @threaded_cached_property
    def mask(self):
        return self._phi <= 0

    def get_direct_viz_level_set_func_for_cell(self, lc):
        cell = self.domain.mesh.cells[lc]
        X = self.domain.mesh.verts[cell]
        T = jmm.bb.Bb33.from_jets(X, self.eik.jet[cell])
        extended_T = jmm.bb.Bb33.from_jets(X, self.extended_eik.jet[cell])
        eps = self._shadow_mask_thresh
        return lambda b: T.f(b) - extended_T.f(b) - eps

    def get_level_set_func_for_cell(self, lc):
        phi_vis = self.get_direct_viz_level_set_func_for_cell(lc)
        phi_theta = self.get_valid_angle_level_set_func_for_cell(lc)
        return lambda b: np.maximum(phi_vis(b), phi_theta(b))

    @property
    def time(self):
        time = np.empty(self.domain.num_verts)
        time[...] = np.inf
        T = self.eik.T[self.mask]
        if not np.isfinite(T).all():
            raise RuntimeError('found bad eikonal values while computing time')
        time[self.mask] = T/Field.speed_of_sound
        return time

class PointSourceField(Field):
    def __init__(self, domain, src_index):
        jet = jmm.jet.Jet3.make_point_source()

        ftype = jmm.defs.Ftype.PointSource
        bd_inds = np.array([src_index])
        bd_T = np.array([jet.f])
        bd_grad_T = np.array([jet.fx, jet.fy, jet.fz])
        super().__init__(domain, None, ftype, bd_inds, bd_T, bd_grad_T)

        self.eik.add_pt_src_BCs(src_index, jet)
        self.extended_eik.add_pt_src_BCs(src_index, jet)

    def __reduce__(self):
        args = (self.domain, self.bd_inds[0])
        return (self.__class__, args)

    @threaded_cached_property
    def _phi(self):
        return self._phi_shadow

    def get_level_set_func_for_cell(self, lc):
        return self.get_direct_viz_level_set_func_for_cell(lc)

class ReflectedField(Field):
    def __init__(self, domain, index, bd_faces, bd_T, bd_grad_T,
                 bd_t_in, parent=None):
        num_faces = bd_faces.shape[0]
        if bd_T.shape[0] != num_faces or bd_grad_T.shape[0] != num_faces:
            raise ValueError('boundary faces and BCs must have the same shape')

        if bd_faces.shape[0] == 0:
            raise ValueError('no boundary faces were passed!')

        ftype = jmm.defs.Ftype.Reflection

        super().__init__(domain, index, ftype, bd_faces, bd_T,
                         bd_grad_T, parent=parent)
        self.bd_t_in = bd_t_in

        for lf, T, grad_T, t_in in zip(bd_faces, bd_T, bd_grad_T, bd_t_in):
            jets = [jmm.jet.Jet3(t, *dt) for t, dt in zip(T, grad_T)]
            self.eik.add_refl_BCs(*lf, *jets, t_in)
            self.extended_eik.add_refl_BCs(*lf, *jets, t_in)

    def __reduce__(self):
        args = (self.domain, self.bd_inds, self.bd_T, self.bd_grad_T, self.bd_t_in)
        kwargs = {'parent': self.parent}
        return (self.__class__, args, kwargs)

    @threaded_cached_property
    def reflector_face_normal(self):
        return self.domain.mesh.get_face_normal(*self.bd_inds[0])

    @property
    def _phi_valid_angle(self):
        # Mask away the points that have BCs, since t_out will be NAN
        # for these points
        mask = np.ones(self.domain.mesh.num_verts, dtype=np.bool_)
        mask[np.unique(self.bd_inds)] = False

        # Compute the sine of the angle the in and out tangent vectors
        # make with the line spanned by the face normal
        nu = self.reflector_face_normal
        proj = np.eye(3) - np.outer(nu, nu)

        t_in_proj = self.eik.t_in[mask]@proj
        t_in_proj /= np.sqrt(np.sum(t_in_proj**2, axis=1)).reshape(-1, 1)

        t_out_proj = self.eik.t_out[mask]@proj
        t_out_proj /= np.sqrt(np.sum(t_out_proj**2, axis=1)).reshape(-1, 1)

        dots = t_in_proj*t_out_proj
        dots = np.clip(dots.sum(1), -1, 1)

        dists = np.empty(self.domain.mesh.num_verts)
        dists[mask] = np.arccos(dots)
        dists[~mask] = 0

        # Get the base-10 exponent of the errors
        log10_dists = np.log10(dists[dists > 0])
        self._valid_angle_mask_thresh = self.h

        return dists - self._valid_angle_mask_thresh

    def get_valid_angle_level_set_func_for_cell(self, lc):
        '''Get a function defined on the cell `lc` that will allow us to
        evaluate the level set function determining whether a point is
        in the "physical zone" or not.

        '''
        nu = self.reflector_face_normal
        cell = self.domain.mesh.cells[lc]
        t_in, t_out = self.eik.t_in[cell], self.eik.t_out[cell]

        # If we're missing some `t_out` vectors with a reflected
        # field, it will be because we're at the boundary. For these
        # points, we can safely just use `grad_T`.
        mask = np.isnan(t_out).any(1)
        assert all(self.eik.has_BCs(_) for _ in cell[mask])
        t_out[mask] = self.eik.grad_T[cell[mask]]

        # Verify that all the `t_in` and `t_out` vectors are OK
        if not np.isfinite(t_in).ravel().all() or \
           not np.isfinite(t_out).ravel().all():
            raise RuntimeError(f'bad cell: {lc}')

        def phi(b):
            lhs = -jmm.slerp.slerp(t_in, b, self.eik.slerp_tol)@nu
            rhs = jmm.slerp.slerp(t_out, b, self.eik.slerp_tol)@nu
            return abs(lhs - rhs) - self._valid_angle_mask_thresh

        return phi

class DiffractedField(Field):
    def __init__(self, domain, index, bd_edges, bd_T, bd_grad_T,
                 parent=None):
        num_edges = bd_edges.shape[0]
        if bd_T.shape[0] != num_edges or bd_grad_T.shape[0] != num_edges:
            raise ValueError('boundary edges and BCs must have the same shape')

        if bd_edges.shape[0] == 0:
            raise ValueError('no boundary edges were passed!')

        ftype = jmm.defs.Ftype.EdgeDiffraction
        super().__init__(domain, index, ftype, bd_edges, bd_T,
                         bd_grad_T, parent=parent)

        for le, T, grad_T in zip(bd_edges, bd_T, bd_grad_T):
            jets = [jmm.jet.Jet3(t, *dt) for t, dt in zip(T, grad_T)]
            try:
                self.eik.add_diff_edge_BCs(*le, *jets)
                self.extended_eik.add_diff_edge_BCs(*le, *jets)
            except:
                import pdb; pdb.set_trace()

    @threaded_cached_property
    def diffractor_tangent_vector(self):

        '''The tangent vector of the diffracting edge'''
        x0, x1 = self.domain.mesh.verts[self.bd_inds[0]]
        t = x1 - x0
        return t/np.linalg.norm(t)

    @property
    def _phi_valid_angle(self):
        # Compute a mask indicating which points aren't on the part of
        # the diffracting edge with BCs
        mask = np.ones(self.domain.mesh.num_verts, dtype=np.bool_)
        mask[self.bd_inds] = False

        # Compute angle in and out vectors make with `t`
        t = self.diffractor_tangent_vector
        lhs = np.arccos(self.eik.t_in[mask]@t)
        rhs = np.arccos(self.eik.t_out[mask]@t)

        # Get the base-10 exponent of the absolute angle difference
        log_diff = np.log10(np.maximum(1e-16, abs(lhs - rhs)))

        # Compute threshold using Otsu's method---result will be an
        # exponent, so convert back
        # p = jmm.util.otsu(log_diff, bins=Field.num_mask_threshold_bins)
        # self._valid_angle_mask_thresh = 10**p
        self._valid_angle_mask_thresh = self.h

        diff = np.empty(mask.shape, dtype=np.float64)
        diff[~mask] = 0
        diff[mask] = lhs - rhs

        # Compute the valid angle mask, taking care to set the
        # boundary indices to `True`, since we don't have `t_out`
        # vectors available along the edge to do the calculation for
        # these indices.
        return abs(diff) - self._valid_angle_mask_thresh

    def get_valid_angle_level_set_func_for_cell(self, lc):
        t = self.diffractor_tangent_vector

        cell = self.domain.mesh.cells[lc]
        t_in, t_out = self.eik.t_in[cell], self.eik.t_out[cell]

        if not np.isfinite(t_in).ravel().all() or \
           not np.isfinite(t_out).ravel().all():
            raise RuntimeError(f'bad cell: {lc}')

        def phi(b):
            lhs = np.arccos(jmm.slerp.slerp(t_in, b, self.eik.slerp_tol)@t)
            rhs = np.arccos(jmm.slerp.slerp(t_out, b, self.eik.slerp_tol)@t)
            return abs(lhs - rhs) - self._valid_angle_mask_thresh

        return phi

    def _impute_diffracted_t_in(self, le, t_in):
        '''This is a specialized function that figures out how to set `t_in`
        BCs in one specific case. When an edge-diffracted field
        generates another edge-diffracted field which is incident on
        the generating field (e.g., two adjacent diffracting edges on
        a polyhedral mesh), there will be one node on the new edge
        which will be incident on the old edge. This node will be
        missing gradient data, which is correct, but we need a value
        for `t_in` at that point. Although not exactly correct, a
        reasonable value to use is the `t_in` vector at the same
        point, taken from the generating field, and rotated around the
        generating diffracting edge until it's aligned with the new
        edge.

        '''
        # Order the indices of `le` so that `l0` indexes the node
        # that's on the diffracting edge (i.e., is missing gradient
        # data), and `l1` indexes the node that *isn't* on the
        # diffracting edge (i.e., isn't).
        l0, l1 = le[np.argsort(np.isfinite(self.eik.grad_T[le]).all(1))]

        # Find the diffracting edge which `l0` is incident on. We'll
        # use this as the rotation axis.
        inc_diff_edges = self.domain.mesh.get_inc_diff_edges(l0)
        if inc_diff_edges.shape[0] < 2:
            raise RuntimeError('entered weird state while imputing t_in')
        le_diff = next(_ for _ in inc_diff_edges
                       if np.isnan(self.eik.grad_T[_]).all())

        # Compute the rotation axis from the edge endpoints.
        q = -np.subtract(*self.domain.mesh.verts[le_diff])
        q /= np.linalg.norm(q)

        # Get the `t_in` vector at `l0` and make sure it's finite
        t_in_0 = self.eik.t_in[l0]
        assert np.isfinite(t_in_0).all()

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
        t_in_imputed = rot@t_in_0

        return t_in_imputed

class MultipleArrivals(Logger):
    def __init__(self, domain, root_field, num_arrivals):
        super().__init__()

        self.domain = domain
        self.root_field = root_field
        self.num_arrivals = num_arrivals

        num_verts = self.domain.num_verts

        self._time = np.empty((num_verts, num_arrivals + 1))
        self._time[...] = np.inf

        self._amplitude = np.empty((num_verts, num_arrivals + 1))
        self._amplitude[...] = np.nan

        self._direction = np.empty((num_verts, num_arrivals + 1, 3))
        self._direction[...] = np.nan

        self.log.info('traversing scattered fields')

        for field in self.root_field.scattered_fields:
            time = field.time

            min_time = time.min()
            min_time_str = str(datetime.timedelta(seconds=min_time))

            num_physical = field.mask.sum()
            perc_physical = 100*num_physical/self.domain.mesh.num_verts

            self.log.info(
                'accepted %s (earliest arrival: %s, physical: %1.1f%%)',
                type(field).__name__, min_time_str, perc_physical)

            if not np.isfinite(min_time):
                raise RuntimeError('bad field')

            self._time[:, -1] = time

            perm = np.argsort(self._time, axis=1)
            self._time = np.take_along_axis(self._time, perm, axis=1)
            self._amplitude = np.take_along_axis(self._amplitude, perm, axis=1)
            for i in range(self._direction.shape[-1]):
                self._direction[:, :, i] = \
                    np.take_along_axis(self._direction[:, :, i], perm, axis=1)

            perc_inserted = \
                100*(perm[field.mask, -1] != num_arrivals).sum()/num_physical
            self.log.info('inserted %1.1f%% of new arrivals', perc_inserted)

            finite_time = np.isfinite(self.time)

            perc_finite = 100*finite_time.sum()/np.prod(finite_time.shape)
            self.log.info('computed %1.1f%% of all arrivals', perc_finite)

            if finite_time.all():
                break

    @property
    def time(self):
        return self._time[:, :-1]

    @property
    def amplitude(self):
        return self._amplitude[:, :-1]

    @property
    def direction(self):
        return self._direction[:, :-1, :]
