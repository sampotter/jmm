import datetime
import logging
import numpy as np
import os
import pickle

from abc import ABC
from cached_property import threaded_cached_property

import jmm.eik
import jmm.jet
import jmm.mesh
import jmm.util

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
        return self.mesh.get_active_reflectors(mask)

    def get_active_diffractors(self, mask):
        return self.mesh.get_active_diffractors(mask)

class Field(ABC, Logger):
    speed_of_sound = 343
    num_mask_threshold_bins = 513

    def __init__(self, domain, bd_inds, bd_T, bd_grad_T):
        super().__init__()

        self.domain = domain
        self.bd_inds = bd_inds
        self.bd_T = bd_T
        self.bd_grad_T = bd_grad_T
        self.eik = jmm.eik.Eik3(self.domain.mesh)
        self.extended_eik = jmm.eik.Eik3(self.domain.extended_mesh)
        self.solved = False
        self._scattered_fields = []

    def __reduce__(self):
        args = (self.domain, self.bd_inds, self.bd_T, self.bd_grad_T)
        return (self.__class__, args)

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
            'solving eikonal equation on extended domain [%1.2fs]',
            jmm.util.toc())

    def _get_reflection_BCs(self, faces):
        # Get the surface normal and reflection matrix for the
        # reflector
        nu = self.domain.mesh.get_face_normal(*faces[0])
        refl = np.eye(nu.size) - 2*np.outer(nu, nu)

        # Traverse the faces and pull out the BCs for the reflection
        bd_faces, bd_T, bd_grad_T, bd_t_in = [], [], [], []
        for lf in faces:
            if np.logical_not(self.mask[lf]).all():
                continue

            # Get the eikonal gradients at the face vertices. Skip
            # this face if any of the gradients graze the surface.
            t_in = np.array([(_[1], _[2], _[3]) for _ in self.eik.jet[lf]])
            if (np.dot(t_in, nu) >= -self.domain.mesh.eps).any():
                continue

            # Reflect the ray directions over the reflector to get the
            # BCs for the eikonal gradient for the reflection
            t_out = np.dot(t_in, refl)

            bd_faces.append(lf)
            bd_T.append([_[0] for _ in self.eik.jet[lf]])
            bd_grad_T.append(t_out)
            bd_t_in.append(t_in)

        return np.array(bd_faces), np.array(bd_T), np.array(bd_grad_T), \
            np.array(bd_t_in)

    def _get_diffraction_BCs(self, edges):
        bd_edges, bd_T, bd_grad_T = [], [], []
        for le in edges:
            if np.logical_not(self.mask[le]).all():
                continue
            bd_edges.append(le)
            bd_T.append([_[0] for _ in self.eik.jet[le]])
            bd_grad_T.append([(_[1], _[2], _[3]) for _ in self.eik.jet[le]])
        return np.array(bd_edges), np.array(bd_T), np.array(bd_grad_T)

    def _init_scattered_fields(self):
        if not self.solved:
            self.solve()

        fields = []

        for _, faces in self.domain.get_active_reflectors(self.mask):
            BCs = self._get_reflection_BCs(faces)
            reflected_field = ReflectedField(self.domain, *BCs)
            fields.append(reflected_field)

        for _, edges in self.domain.get_active_diffractors(self.mask):
            BCs = self._get_diffraction_BCs(edges)
            diffracted_field = DiffractedField(self.domain, *BCs)
            fields.append(diffracted_field)

        self._scattered_fields = fields

    @property
    def scattered_fields(self):
        fields = [self]

        while fields:
            field = fields.pop(0)

            # DEBUG: dump the field before solving it for debugging
            pickled_field_path = 'field.pickle'
            if os.path.exists(pickled_field_path):
                os.remove(pickled_field_path)
            with open(pickled_field_path, 'wb') as f:
                pickle.dump(field, f)

            field.solve()
            yield field

            if not self._scattered_fields:
                self._init_scattered_fields()
            fields.extend(field._scattered_fields)

            fields = sorted(fields, key=lambda _: _.eik.T.min())

    @threaded_cached_property
    def shadow_mask(self):
        T = self.eik.T
        extended_T = self.extended_eik.T[:self.domain.num_verts]

        # Get the base-10 exponent of the absolute domain errors
        log_diff = np.log10(np.maximum(1e-16, abs(T - extended_T)))

        # Compute a threshold in order to separate the distribution of
        # exponents using Otsu's method
        p = jmm.util.otsu(log_diff, bins=Field.num_mask_threshold_bins)
        self._shadow_mask_thresh = 10**p

        return T - extended_T > self._shadow_mask_thresh

    @property
    def mask(self):
        return np.logical_and(~self.shadow_mask, self.valid_angle_mask)

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

        bd_inds = np.array([src_index])
        bd_T = np.array([jet.f])
        bd_grad_T = np.array([jet.fx, jet.fy, jet.fz])
        super().__init__(domain, bd_inds, bd_T, bd_grad_T)

        self.eik.add_trial(src_index, jet)
        self.extended_eik.add_trial(src_index, jet)

    def __reduce__(self):
        return (self.__class__, (self.domain, self.bd_inds[0]))

    @property
    def valid_angle_mask(self):
        return np.ones(self.domain.num_verts, dtype=bool)

class ReflectedField(Field):
    def __init__(self, domain, bd_faces, bd_T, bd_grad_T, bd_t_in):
        num_faces = bd_faces.shape[0]
        if bd_T.shape[0] != num_faces or bd_grad_T.shape[0] != num_faces:
            raise ValueError('boundary faces and BCs must have the same shape')

        super().__init__(domain, bd_faces, bd_T, bd_grad_T)
        self.bd_t_in = bd_t_in

        for lf, T, grad_T, t_in in zip(bd_faces, bd_T, bd_grad_T, bd_t_in):
            jets = [jmm.jet.Jet3(t, *dt) for t, dt in zip(T, grad_T)]
            self.eik.add_valid_bdf(*lf, *jets, t_in)
            self.extended_eik.add_valid_bdf(*lf, *jets, t_in)

    def __reduce__(self):
        return (self.__class__, (self.domain, self.bd_inds, self.bd_T,
                                 self.bd_grad_T, self.bd_t_in))

    @property
    def valid_angle_mask(self):
        nu = self.domain.mesh.get_face_normal(*self.bd_inds[0])

        # Compute the sine of the angle the in and out tangent vectors
        # make with the face normal
        lhs = np.sqrt(1 - (self.eik.t_in@nu)**2)
        rhs = np.sqrt(1 - (self.eik.t_out@nu)**2)

        # Get the base-10 exponent of the absolute difference between
        # the sines
        log_diff = np.log10(np.maximum(1e-16, abs(lhs - rhs)))

        # Compute threshold using Otsu's method---result will be an
        # exponent, so convert back
        p = jmm.util.otsu(log_diff, bins=Field.num_mask_threshold_bins)
        self._valid_angle_mask_thresh = 10**p

        return abs(lhs - rhs) < self._valid_angle_mask_thresh

class DiffractedField(Field):
    def __init__(self, domain, bd_edges, bd_T, bd_grad_T):
        num_edges = bd_edges.shape[0]
        if bd_T.shape[0] != num_edges or bd_grad_T.shape[0] != num_edges:
            raise ValueError('boundary faces and BCs must have the same shape')

        super().__init__(domain, bd_edges, bd_T, bd_grad_T)

        for le, T, grad_T in zip(bd_edges, bd_T, bd_grad_T):
            jets = [jmm.jet.Jet3(t, *dt) for t, dt in zip(T, grad_T)]
            self.eik.add_valid_bde(*le, *jets)
            self.extended_eik.add_valid_bde(*le, *jets)

    @property
    def valid_angle_mask(self):
        # Compute tangent vector of diffracting edge
        x0, x1 = self.domain.mesh.verts[self.bd_inds[0]]
        t = x1 - x0
        t /= np.linalg.norm(t)

        # Compute a mask indicating which points aren't on the part of
        # the diffracting edge with BCs
        mask = np.ones(self.domain.mesh.num_verts, dtype=np.bool)
        mask[self.bd_inds] = False

        # Compute angle in and out vectors make with `t`
        lhs = np.arccos(self.eik.t_in[mask]@t)
        rhs = np.arccos(self.eik.t_out[mask]@t)

        # Get the base-10 exponent of the absolute angle difference
        log_diff = np.log10(np.maximum(1e-16, abs(lhs - rhs)))

        # Compute threshold using Otsu's method---result will be an
        # exponent, so convert back
        p = jmm.util.otsu(log_diff, bins=Field.num_mask_threshold_bins)
        self._valid_angle_mask_thresh = 10**p

        # Compute the valid angle mask, taking care to set the
        # boundary indices to `True`, since we don't have `t_out`
        # vectors available along the edge to do the calculation for
        # these indices.
        valid_angle_mask = np.ones(self.domain.mesh.num_verts, dtype=np.bool)
        valid_angle_mask[mask] = abs(lhs - rhs) < self._valid_angle_mask_thresh
        valid_angle_mask[self.bd_inds] = True
        return valid_angle_mask

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
            perc_physical = 100*field.mask.sum()/self.domain.mesh.num_verts
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

            perc_inserted = 100*(perm[:, -1] != num_arrivals).sum()/num_verts
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
