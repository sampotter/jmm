import numpy as np

from enum import Enum

from libc.math cimport isfinite

from jmm.bmesh cimport Bmesh33, Bmesh33Cell
from jmm.geom cimport Ray3, Rect3
from jmm.mesh cimport Mesh2, Mesh3, Mesh2Tri, Mesh3Tetra

class RobjType(Enum):
    Bmesh33Cell = 0
    Mesh2Tri = 1
    Mesh3Tetra = 2
    Tri3 = 3
    Tetra3 = 4

class RtreeSplitStrategy(Enum):
    SurfaceArea = 0

cdef class Robj:
    @staticmethod
    cdef from_ptr(const robj *obj):
        robj = Robj()
        robj._obj = obj
        return robj

    @property
    def type(self):
        return RobjType(robj_get_type(self._obj))

    @property
    def centroid(self):
        p = np.empty((3,), dtype=np.float64)
        cdef dbl[::1] p_mv = p
        robj_get_centroid(self._obj, &p_mv[0])
        return p

    def overlaps_box(self, Rect3 box):
        return robj_isects_bbox(self._obj, &box._rect)

    def intersect(self, Ray3 ray):
        cdef dbl t
        cdef bool hit = robj_intersect(self._obj, &ray._ray, &t)
        return hit, t

    def get(self):
        return self.astype({
            RobjType.Bmesh33Cell: Bmesh33Cell,
            RobjType.Mesh2Tri: Mesh2Tri,
            RobjType.Mesh3Tetra: Mesh3Tetra,
        }[self.type])

    def astype(self, Class):
        Classes = {Bmesh33Cell, Mesh2Tri, Mesh3Tetra}
        if isinstance(Class, type):
            if Class is Bmesh33Cell:
                return Bmesh33Cell.from_robj_ptr(self._obj)
            elif Class is Mesh2Tri:
                return Mesh2Tri.from_robj_ptr(self._obj)
            elif Class is Mesh3Tetra:
                return Mesh3Tetra.from_robj_ptr(self._obj)
            else:
                raise Exception(f'`Class` should be one of: {Classes}')
        else:
            raise Exception(f'`Class` should be an instance of type')

cdef class Isect:
    def __repr__(self):
        return f'Isect(t = {self.t}, obj = {self.obj})'

    @property
    def hit(self):
        return isfinite(self._isect.t)

    @property
    def t(self):
        return self._isect.t

    @property
    def obj(self):
        if np.isfinite(self.t):
            return Robj.from_ptr(self._isect.obj)

cdef class Rtree:
    '''An R-tree data structure, intended to be used for speeding up
raytracing and other basic geometric queries. Right now, it can only
be constructed from :class:`Mesh2` instances, but this limitation will
be removed in the near future.

    :param mesh: A triangle mesh stored as a :class:`Mesh2`
        instance. An R-tree will be built for this mesh.
    '''
    def __cinit__(self, leaf_thresh=32,
                  split_strategy=RtreeSplitStrategy.SurfaceArea,
                  should_alloc_and_init=True):
        if should_alloc_and_init:
            self.ptr_owner = True
            rtree_alloc(&self._rtree)
            rtree_init(self._rtree, leaf_thresh, split_strategy.value)

    def __dealloc__(self):
        if self.ptr_owner:
            rtree_deinit(self._rtree)
            rtree_dealloc(&self._rtree)

    @staticmethod
    cdef from_ptr(rtree *rtree_ptr, ptr_owner=False):
        rtree = Rtree(should_alloc_and_init=False)
        rtree.ptr_owner = ptr_owner
        rtree._rtree = rtree_ptr
        return rtree

    @staticmethod
    def from_mesh2(Mesh2 mesh, leaf_thresh=32,
                   split_strategy=RtreeSplitStrategy.SurfaceArea):
        rtree = Rtree(leaf_thresh, split_strategy)
        rtree_insert_mesh2(rtree._rtree, mesh.mesh)
        rtree_build(rtree._rtree)
        return rtree

    @property
    def bounding_box(self):
        return Rect3(rtree_get_bbox(self._rtree))

    @property
    def num_leaf_nodes(self):
        return rtree_get_num_leaf_nodes(self._rtree)

    def copy(self):
        cdef rtree *rtree = rtree_copy(self._rtree)
        return Rtree.from_ptr(rtree, ptr_owner=True)

    cdef insert_bmesh33(self, Bmesh33 bmesh):
        rtree_insert_bmesh33(self._rtree, bmesh.bmesh)

    cdef insert_mesh2(self, Mesh2 mesh):
        rtree_insert_mesh2(self._rtree, mesh.mesh)

    cdef insert_mesh3(self, Mesh3 mesh):
        rtree_insert_mesh3(self._rtree, mesh.mesh)

    def insert(self, obj):
        # TODO: this is a gross way to do this, but it works for
        # now. Would be nice to come up with a simpler way to dispatch
        # on the type of obj. One thing that would probably work is
        # implementing an `insert_into_rtree` method for each of the
        # classes listed below.
        Classes = {'Bmesh33', 'Mesh2', 'Mesh3'}
        if isinstance(obj, Bmesh33):
            self.insert_bmesh33(obj)
        elif isinstance(obj, Mesh2):
            self.insert_mesh2(obj)
        elif isinstance(obj, Mesh3):
            self.insert_mesh3(obj)
        else:
            raise Exception(f'obj should be an instance of: {Classes}')

    def build(self):
        rtree_build(self._rtree)

    def intersect(self, Ray3 ray):
        isect = Isect()
        rtree_intersect(self._rtree, &ray._ray, &isect._isect)
        return isect
