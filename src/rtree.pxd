from defs cimport bool, dbl
from geom cimport rect3, ray3
from mesh2 cimport mesh2

cdef extern from "rtree.h":
    cdef enum robj_type:
        ROBJ_MESH2_TRI
        ROBJ_TRI3
        ROBJ_TETRA3

    ctypedef struct isect:
        const ray3 *ray
        dbl t
        robj *obj

    cdef struct robj:
        pass

    robj_type robj_get_type(const robj *obj)
    void robj_insert_into_bbox(const robj *obj, rect3 *bbox)
    void robj_get_centroid(const robj *obj, dbl c[3])
    bool robj_isects_bbox(const robj *obj, const rect3 *bbox)
    bool robj_intersect(const robj *obj, const ray3 *ray, dbl *t)

    cdef struct rtree:
        pass

    void rtree_alloc(rtree **rtree)
    void rtree_dealloc(rtree **rtree)
    void rtree_init(rtree *rtree)
    void rtree_deinit(rtree *rtree)
    void rtree_insert_mesh2(rtree *rtree, const mesh2 *mesh)
    void rtree_build(rtree *rtree)
    rect3 rtree_get_bbox(const rtree *rtree)
    size_t rtree_get_num_leaf_nodes(const rtree *rtree)
    bool rtree_query_bbox(const rtree *rtree, const rect3 *bbox)
    void rtree_intersect(const rtree *rtree, const ray3 *ray, isect *isect)
    void rtree_intersectN(const rtree *rtree, const ray3 *ray, size_t n, isect *isects)
