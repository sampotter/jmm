import logging
import numpy as np
import os
import pickle

from abc import ABC
from cached_property import threaded_cached_property

import jmm.eik
import jmm.jet
import jmm.mesh

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

        self.extended_mesh = jmm.mesh.Mesh3(extended_verts, extended_cells)
        for le in self.mesh.get_diff_edges():
            self.extended_mesh.set_boundary_edge(*le, True)

        self.log.info('extended domain has %d vertices and %d cells',
                      extended_verts.shape[0], extended_cells.shape[0])

        self.h = self.extended_mesh.min_edge_length
        self.mask_threshold = np.log(1/self.h)*self.h**2

    @property
    def num_verts(self):
        return self.mesh.num_verts

    def get_active_reflectors(self, mask):
        return self.mesh.get_active_reflectors(mask)

    def get_active_diffractors(self, mask):
        return self.mesh.get_active_diffractors(mask)

class Field(ABC, Logger):
    speed_of_sound = 343
    num_shadow_mask_threshold_bins = 513

    def __init__(self, domain, bd_inds, bd_T, bd_grad_T):
        super().__init__()

        self.domain = domain
        self.bd_inds = bd_inds
        self.bd_T = bd_T
        self.bd_grad_t = bd_grad_T
        self.eik = jmm.eik.Eik3(self.domain.mesh)
        self.extended_eik = jmm.eik.Eik3(self.domain.extended_mesh)
        self.solved = False
        self._scattered_fields = []
        self._shadow_mask_threshold = None

    def __reduce__(self):
        args = (self.domain, self.bd_inds, self.bd_T, self.bd_grad_T)
        return (self.__class__, args)

    def solve(self):
        if self.solved:
            return
        self.solved = True

        self.log.info('solving eikonal equation on domain')
        self.eik.solve()

        self.log.info('solving eikonal equation on extended domain')
        self.extended_eik.solve()

    def _get_reflection_BCs(self, faces):
        with open('faces.pickle', 'wb') as f:
            pickle.dump(faces, f)
        nu = self.domain.mesh.get_face_normal(*faces[0])
        refl = np.eye(nu.size) - 2*np.outer(nu, nu)
        T, grad_T = [], []
        for lf in faces:
            print(lf)
            if np.logical_not(self.mask[lf]).any():
                continue
            T.append([_[0] for _ in self.eik.jet[lf]])
            print('  ', T[-1])
            grad_T.append([(_[1], _[2], _[3]) for _ in self.eik.jet[lf]])
            print('  ', grad_T[-1])
            grad_T[-1] = np.dot(grad_T[-1], refl)
            print('  ', grad_T[-1])
        return faces, np.array(T), np.array(grad_T)

    def _get_diffraction_BCs(self, edges):
        T, grad_T = [], []
        for le in edges:
            if np.logical_not(self.mask[le]).any():
                continue
            T.append([_[0] for _ in self.eik.jet[le]])
            grad_T.append([(_[1], _[2], _[3]) for _ in self.eik.jet[le]])
        return edges, np.array(T), np.array(grad_T)

    def _init_scattered_fields(self):
        if not self.solved:
            self.solve()

        for _, faces in self.domain.get_active_reflectors(self.mask):
            BCs = self._get_reflection_BCs(faces)
            reflected_field = ReflectedField(self.domain, *BCs)
            self._scattered_fields.append(reflected_field)

        for _, edges in self.domain.get_active_diffractors(self.mask):
            BCs = self._get_diffraction_BCs(edges)
            diffracted_field = DiffractedField(self.domain, *BCs)
            self._scattered_fields.append(diffracted_field)

        key = lambda field: field.T.min()
        self._scattered_fields = sorted(self._scattered_fields, key)

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

            fields = sorted(fields, key=lambda _: _.min_T)

    @threaded_cached_property
    def shadow_mask(self):
        # Get the eikonal for the regular and extended domains
        T = np.array([_[0] for _ in self.eik.jet])
        extended_T = np.array([
            _[0] for _ in self.extended_eik.jet[:self.domain.num_verts]])

        # Get the base-10 exponent of the absolute domain errors
        log_diff = np.log10(
            np.maximum(1e-16, abs(T - extended_T[:self.domain.num_verts])))

        # Compute the histogram of the exponents
        bin_counts, edges = np.histogram(
            log_diff, bins=Field.num_shadow_mask_threshold_bins)
        bin_centers = (edges[:-1] + edges[1:])/2

        def ssd(counts, centers):
            n = counts.sum()
            mu = np.sum(centers * counts) / n
            return np.sum(counts * ((centers - mu) ** 2))

        # Use Otsu's method to compute the optimum shadow mask threshold
        ssds = []
        for k in range(1, bin_counts.size):
            left_ssd = ssd(bin_counts[:k], bin_centers[:k])
            right_ssd = ssd(bin_counts[k:], bin_centers[k:])
            ssds.append(left_ssd + right_ssd)
        self._shadow_mask_threshold = 10**bin_centers[np.argmin(ssds)]

        return T - extended_T > self._shadow_mask_threshold

    @threaded_cached_property
    def mask(self):
        return np.logical_and(~self.shadow_mask, self.valid_angle_mask)

    @threaded_cached_property
    def time(self):
        time = np.empty(self.domain.num_verts)
        time[...] = np.inf
        time[self.mask] = np.array([_[0] for _ in self.eik.jet])[self.mask]
        time[self.mask] /= Field.speed_of_sound
        return time

    @threaded_cached_property
    def min_time(self):
        return self.time.min()

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

    @threaded_cached_property
    def valid_angle_mask(self):
        return np.ones(self.domain.num_verts, dtype=bool)

class ReflectedField(Field):
    def __init__(self, domain, bd_faces, bd_T, bd_grad_T):
        super().__init__(domain, bd_faces, bd_T, bd_grad_T)

        for lf, T, grad_T in zip(bd_faces, bd_T, bd_grad_T):
            jets = [jmm.jet.Jet3(t, *dt) for t, dt in zip(T, grad_T)]
            self.eik.add_valid_bdf(*lf, *jets)
            self.extended_eik.add_valid_bdf(*lf, *jets)

    @threaded_cached_property
    def valid_angle_mask(self):
        Q = self.domain.mesh.get_tangent_plane_basis(self.boundary_faces[0])
        lhs = np.sqrt(1 - (self.eik.t_in@Q)**2)
        rhs = np.sqrt(1 - (self.eik.t_out@Q)**2)
        return all(abs(lhs - rhs) < self.domain.mask_threshold)

class DiffractedField(Field):
    def __init__(self, domain, bd_edges, bd_T, bd_grad_T):
        super().__init__(domain, bd_faces, bd_T, bd_grad_T)

        for le, T, grad_T in zip(bd_edges, bd_T, bd_grad_T):
            jets = [jmm.jet.Jet3(t, *dt) for t, dt in zip(T, grad_T)]
            self.eik.add_valid_bde(*le, *jets)
            self.extended_eik.add_valid_bde(*le, *jets)

    @threaded_cached_property
    def valid_angle_mask(self):
        x0, x1 = self.domain.mesh.verts[self.boundary_edges[0]]
        q = x1 - x0
        q /= np.linalg.norm(q)
        lhs = np.sqrt(1 - (self.eik.t_in@Q)**2)
        rhs = np.sqrt(1 - (self.eik.t_out@Q)**2)
        return all(abs(lhs - rhs) < self.domain.mask_threshold)

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
            self.log.info('%s: %1.2fs', type(field).__name__, field.min_time)

            self._time[:, -1] = field.time

            perm = np.argsort(self._time, axis=1)
            self._time = np.take_along_axis(self._time, perm, axis=1)
            self._amplitude = np.take_along_axis(self._amplitude, perm, axis=1)
            for i in range(self._direction.shape[-1]):
                self._direction[:, :, i] = \
                    np.take_along_axis(self._direction[:, :, i], perm, axis=1)

            finite_time = np.isfinite(self.time)

            perc_finite = finite_time.sum()/np.prod(finite_time.shape)
            self.log.info('computed %1.1f%% of arrivals', 100*perc_finite)

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
