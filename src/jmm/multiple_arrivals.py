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
    minimum_amplitude = 1e-3 # == -60 dB

    def __init__(self, domain, index, ftype, bd_inds, bd_T, bd_grad_T,
                 parent=None):
        super().__init__()

        self.domain = domain
        self.index = index
        self.bd_inds = bd_inds
        self.bd_T = bd_T
        self.bd_grad_T = bd_grad_T
        self.parent = parent

        orig = None
        if self.parent is not None:
            orig = self.parent.eik

        # TODO: instantiate these lazily to save memory!
        self.eik = jmm.eik.Eik3.from_mesh_and_ftype(domain.mesh, ftype, orig)

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

        mesh, fields = self.domain.mesh, []

        for index, faces in enumerate(mesh.reflectors):
            if isinstance(self, ReflectedField) and index == self.index:
                continue
            BCs = self._get_reflection_BCs(faces)
            if not BCs:
                continue
            fields.append(ReflectedField(self.domain, index, *BCs, parent=self))

        for index, edges in enumerate(mesh.diffractors):
            if isinstance(self, DiffractedField) and index == self.index:
                continue
            BCs = self._get_diffraction_BCs(edges)
            if not BCs:
                continue
            fields.append(DiffractedField(self.domain, index, *BCs,parent=self))

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
    def direction(self):
        direction = np.empty((self.domain.num_verts, 3))
        direction[...] = self.eik.grad_T[...]
        if not np.isfinite(direction[~self._boundary_mask]).all():
            raise RuntimeError('found bad arrival directions')
        return direction

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
        tmp = np.log(Field.minimum_amplitude)
        scale = np.exp(tmp*((arc_length - self.h)/self.h)**2)
        scale[self._boundary_mask | (arc_length < self.h)] = 1

        return scale

class PointSourceField(Field):
    def __init__(self, domain, src_index):
        jet = jmm.jet.Jet3.make_point_source()

        ftype = jmm.defs.Ftype.PointSource
        bd_inds = np.array([src_index])
        bd_T = np.array([jet.f])
        bd_grad_T = np.array([jet.fx, jet.fy, jet.fz])
        super().__init__(domain, None, ftype, bd_inds, bd_T, bd_grad_T)

        self.eik.add_pt_src_BCs(src_index, jet)

    def __reduce__(self):
        args = (self.domain, self.bd_inds[0])
        return (self.__class__, args)

    @property
    def amplitude(self):
        return self._scale/self.r

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

    def __reduce__(self):
        args = (self.domain, self.bd_inds, self.bd_T, self.bd_grad_T, self.bd_t_in)
        kwargs = {'parent': self.parent}
        return (self.__class__, args, kwargs)

    @threaded_cached_property
    def reflector_face_normal(self):
        return self.domain.mesh.get_face_normal(*self.bd_inds[0])

    @property
    def _scale(self):
        nu = self.reflector_face_normal
        proj = np.eye(3) - np.outer(nu, nu)

        t_in_proj = self.eik.t_in@proj
        t_in_proj /= np.sqrt(np.sum(t_in_proj**2, axis=1)).reshape(-1, 1)

        t_out_proj = self.eik.t_out@proj
        t_out_proj /= np.sqrt(np.sum(t_out_proj**2, axis=1)).reshape(-1, 1)

        dot = np.multiply(t_in_proj, t_out_proj).sum(1)
        dot = np.clip(dot, -1, 1)

        arc_length = np.arccos(dot)

        tmp = np.log(Field.minimum_amplitude)
        scale = np.exp(tmp*((arc_length - self.h)/self.h)**2)
        scale[self._boundary_mask | (arc_length < self.h)] = 1

        return scale*super()._scale

    @property
    def amplitude(self):
        rho = self.domain.refl_coef # TODO: make this varying
        return rho*self._scale/self.r

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
            except:
                import pdb; pdb.set_trace()

    @property
    def diffractor_tangent_vector(self):
        '''The tangent vector of the diffracting edge'''
        x0, x1 = self.domain.mesh.verts[self.bd_inds[0]]
        e = x1 - x0
        return e/np.linalg.norm(e)

    @property
    def _o_face_normal(self):
        le = self.bd_inds[0]
        l0, l1 = le
        x0, x1 = self.domain.mesh.verts[le]

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

        n = self.domain.mesh.get_face_normal(lf1)
        if n@no > 0:
            return no

        l2 = next(_ for _ in lf2 if _ != l0 and _ != l1)
        x2 = self.domain.mesh.verts[l2]
        t = x2 - x0
        t /= np.linalg.norm(x2 - x0)
        no = np.cross(e, t)

        n = self.domain.mesh.get_face_normal(lf2)
        if n@no > 0:
            return no

        raise RuntimeError('entered a weird state')

    @property
    def _scale(self):
        t = self.diffractor_tangent_vector

        arc_length_in = np.arccos(np.clip(self.eik.t_in@t, -1, 1))
        arc_length_out = np.arccos(np.clip(self.eik.t_out@t, -1, 1))

        d = arc_length_in - arc_length_out

        tmp = np.log(field.minimum_amplitude)
        scale = np.exp(tmp*((d - self.h)/self.h)**2)
        scale[self._boundary_mask | (d < self.h)] = 1

        return scale*super()._scale

    @property
    def amplitude(self):
        A_in, rho_in = self.eik.A_in, self.eik.s, self.eik.rho_in, self.eik.s

        D = jmm.util.D_from_geometry(self.om, self.wedge_angle,
                                     self._o_face_normal,
                                     self.diffractor_tangent_vector,
                                     self.eik.t_out, self.eik.t_in,
                                     self.refl_coef)

        return A_in*D*np.sqrt(rho_in/(s*(s + rho_in)))

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

        self._fields = []

    def traverse(self):
        self.log.info('traversing scattered fields')

        for field in self.root_field.scattered_fields:
            self._fields.append(field)

            time = field.time

            min_time = time.min()
            min_time_str = str(datetime.timedelta(seconds=min_time))

            self.log.info('accepted %s (earliest arrival: %s)',
                          type(field).__name__, min_time_str)

            if not np.isfinite(min_time):
                raise RuntimeError('bad field')

            amp_mask = field.amplitude >= Field.minimum_amplitude

            self._time[:, -1] = time
            self._time[~amp_mask, -1] = np.inf

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

    @property
    def time(self):
        return self._time[:, :-1]

    @property
    def amplitude(self):
        return self._amplitude[:, :-1]

    @property
    def direction(self):
        return self._direction[:, :-1, :]
