from defs cimport bool, dbl
from geom cimport rect3, ray3, Ray3
from mesh2 cimport mesh2

cdef extern from "rtree.h":
    cdef struct isect:
        ray3 ray
        dbl t
        void *obj
    cdef struct rtree_node:
        pass
    cdef struct rtree:
        pass
    void rtree_alloc(rtree **rtree)
    void rtree_dealloc(rtree **rtree)
    void rtree_init_from_tri_mesh(rtree *rtree, const mesh2 *mesh)
    void rtree_deinit(rtree *rtree)
    rect3 rtree_get_bbox(const rtree *rtree)
    size_t rtree_get_num_leaf_nodes(const rtree *rtree)
    bool rtree_query_bbox(const rtree *rtree, const rect3 *bbox)
    bool rtree_intersect(const rtree *rtree, const ray3 *ray, isect *isect)

cdef class Isect:
    cdef:
        isect _isect

cdef class Rtree:
    cdef:
        rtree *_rtree
