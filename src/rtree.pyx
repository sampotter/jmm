from geom import Rect3
from mesh2 cimport Mesh2

cdef class Isect:
    def __cinit__(self, isect isect):
        self._isect = isect

cdef class Rtree:
    def __cinit__(self, Mesh2 mesh):
        rtree_alloc(&self._rtree)
        rtree_init_from_tri_mesh(self._rtree, mesh.get_mesh_ptr())

    def __dealloc__(self):
        rtree_deinit(self._rtree)
        rtree_dealloc(&self._rtree)

    @property
    def bounding_box(self):
        return Rect3(rtree_get_bbox(self._rtree))

    @property
    def num_leaf_nodes(self):
        return rtree_get_num_leaf_nodes(self._rtree)

    def intersect(self, Ray3 ray):
        cdef isect isect
        if rtree_intersect(self._rtree, &ray._ray, &isect):
            return Isect()
