import logging
import numpy as np

from cached_property import threaded_cached_property

import jmm.eik
import jmm.jet

class Field(object):
    speed_of_sound = 343

    def __init__(self, tree):
        self.parent_field = None
        self.tree = tree
        self.eik = jmm.eik.Eik3(tree.mesh)
        self.extended_eik = jmm.eik.Eik3(tree.extended_mesh)
        self.solved = False
        self._scattered_fields = []
        self.log = logging.getLogger('Field')

    def solve(self):
        if self.solved:
            return
        self.solved = True

        self.log.info('solving eikonal equation on domain')
        self.eik.solve()

        self.log.info('solving eikonal equation on extended domain')
        self.extended_eik.solve()

    def _init_scattered_fields(self):
        if not self.solved:
            self.solve()

        for _, faces in self.tree.mesh.get_active_reflectors(self.shadow_mask):
            reflected_field = ReflectedField(self, faces)
            self._scattered_fields.append(reflected_field)

        for _, edges in self.tree.mesh.get_active_diffractors(self.shadow_mask):
            diffracted_field = DiffractedField(self, edges)
            self._scattered_fields.append(diffracted_field)

        key = lambda field: field.T.min()
        self._scattered_fields = sorted(self._scattered_fields, key)

    @property
    def scattered_fields(self):
        if not self._scattered_fields:
            self._init_scattered_fields()
        return self._scattered_fields

    @threaded_cached_property
    def shadow_mask(self):
        n = self.tree.mesh.num_verts
        T = np.array([_[0] for _ in self.eik.jet])
        extended_T = np.array([_[0] for _ in self.extended_eik.jet[:n]])
        return T - extended_T > self.tree.mask_threshold

    @threaded_cached_property
    def mask(self):
        return np.logical_and(~self.shadow_mask, self.valid_angle_mask)

    @threaded_cached_property
    def time(self):
        time = np.empty(self.tree.mesh.num_verts)
        time[...] = np.inf
        time[self.mask] = np.array([_[0] for _ in self.eik.jet])[self.mask]
        time[self.mask] /= Field.speed_of_sound
        return time

    @threaded_cached_property
    def min_time(self):
        return self.time.min()

class PointSourceField(Field):
    def __init__(self, tree, src_index):
        super().__init__(tree)
        self.eik.add_trial(src_index, jmm.jet.Jet3.make_point_source())
        self.extended_eik.add_trial(src_index, jmm.jet.Jet3.make_point_source())

    @threaded_cached_property
    def valid_angle_mask(self):
        return np.ones(self.tree.mesh.num_verts, dtype=bool)

class ReflectedField(Field):
    def __init__(self, parent_field, faces):
        super().__init__(parent_field.tree)

        self.parent_field = parent_field
        self.faces = faces

        bd_inds = np.unique(self.faces)
        bd_inds = bd_inds[self.parent_field.mask[bd_inds]]

        bd_jets = self.parent_field.eik.jet[bd_inds]

        T = np.array([_[0] for _ in bd_jets])

        nu = self.tree.mesh.get_face_normal(bd_inds[0])
        refl = np.eye(nu.size) - 2*np.outer(nu, nu)

        grad_T = np.array([(_[1], _[2], _[3]) for _ in bd_jets])
        grad_T = grad_T@refl

        for lf in self.faces:
            if ~self.parent_field.mask[lf].any():
                continue
            jets = [jmm.jet.Jet3(*_) for _ in self.parent_field.eik.jet[lf]]
            self.eik.add_valid_bdf(*lf, *jets)
            self.extended_eik.add_valid_bdf(*lf, *jets)

    @threaded_cached_property
    def valid_angle_mask(self):
        Q = self.tree.mesh.get_tangent_plane_basis(self.faces[0])
        lhs = np.sqrt(1 - (self.eik.t_in@Q)**2)
        rhs = np.sqrt(1 - (self.eik.t_out@Q)**2)
        return all(abs(lhs - rhs) < self.tree.mask_threshold)

class DiffractedField(Field):
    def __init__(self, parent_field, edges):
        super().__init__(parent_field.tree)

        self.parent_field = parent_field
        self.edges = edges

        bd_inds = np.unique(self.edges)
        bd_inds = bd_inds[self.parent_field.mask[bd_inds]]

        for le in self.edges:
            if ~self.parent_field.mask[le].any():
                continue
            jets = [jmm.jet.Jet3(*_) for _ in self.parent_field.eik.jet[le]]
            self.eik.add_valid_bde(*le, *jets)
            self.extended_eik.add_valid_bde(*le, *jets)

    @threaded_cached_property
    def valid_angle_mask(self):
        x0, x1 = self.tree.mesh.verts[self.edges[0]]
        q = x1 - x0
        q /= np.linalg.norm(q)
        lhs = np.sqrt(1 - (self.eik.t_in@Q)**2)
        rhs = np.sqrt(1 - (self.eik.t_out@Q)**2)
        return all(abs(lhs - rhs) < self.tree.mask_threshold)

class Tree(object):
    def __init__(self, mesh, extended_mesh):
        self.mesh = mesh
        self.extended_mesh = extended_mesh
        self.h = self.extended_mesh.min_edge_length
        self.mask_threshold = np.log(1/self.h)*self.h**2
        self.root_field = None

    def fields(self):
        fields = [self.root_field]
        while fields:
            field = fields.pop(0)
            field.solve()
            yield field
            fields.extend(field.scattered_fields)
            fields = sorted(fields, key=lambda _: _.min_T)

class PointSourceTree(Tree):
    def __init__(self, mesh, extended_mesh, src_index):
        super().__init__(mesh, extended_mesh)
        self.root_field = PointSourceField(self, src_index)

class MultipleArrivals(object):
    def __init__(self, tree, num_arrivals):
        self.tree = tree
        self.num_arrivals = num_arrivals
        self.log = logging.getLogger('MultipleArrivals')

        num_verts = self.tree.mesh.num_verts

        self._time = np.empty((num_verts, num_arrivals + 1))
        self._time[...] = np.inf

        self._amplitude = np.empty((num_verts, num_arrivals + 1))
        self._amplitude[...] = np.nan

        self._direction = np.empty((num_verts, num_arrivals + 1, 3))
        self._direction[...] = np.nan

        self.log.info('traversing %s fields', type(self.tree))

        for field in self.tree.fields():
            self.log.info(
                'getting arrivals from %s (min time: %1.2fs)',
                type(field), field.min_time)

            self._time[:, -1] = field.time

            perm = np.argsort(self._time, axis=1)
            self._time = np.take_along_axis(self._time, perm, axis=1)
            self._amplitude = np.take_along_axis(self._amplitude, perm, axis=1)
            self._direction = np.take_along_axis(self._direction, perm, axis=1)

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
        return self._direction[:, :-1]
