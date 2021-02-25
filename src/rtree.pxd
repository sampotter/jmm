from defs cimport bool, dbl
from geom cimport rect3, ray3
from mesh2 cimport mesh2
from mesh3 cimport mesh3

cdef extern from "rtree.h":
    cdef enum robj_type:
        ROBJ_MESH2_TRI
        ROBJ_MESH3_TETRA
        ROBJ_TRI3
        ROBJ_TETRA3

    ctypedef struct isect:
        dbl t
        robj *obj

    cdef enum rtree_split_strategy:
        RTREE_SPLIT_STRATEGY_SURFACE_AREA

    cdef struct robj:
        pass

    robj_type robj_get_type(const robj *obj)
    const void *robj_get_data(const robj *obj)
    void robj_insert_into_bbox(const robj *obj, rect3 *bbox)
    void robj_get_centroid(const robj *obj, dbl c[3])
    bool robj_isects_bbox(const robj *obj, const rect3 *bbox)
    bool robj_intersect(const robj *obj, const ray3 *ray, dbl *t)

    cdef struct rtree:
        pass

    void rtree_alloc(rtree **rtree)
    void rtree_dealloc(rtree **rtree)
    void rtree_init(rtree *rtree, size_t leaf_thresh,
                    rtree_split_strategy split_strategy)
    void rtree_deinit(rtree *rtree)
    rtree *rtree_copy(const rtree *rtree)
    void rtree_insert_mesh2(rtree *rtree, const mesh2 *mesh)
    void rtree_insert_mesh3(rtree *rtree, const mesh3 *mesh)
    void rtree_build(rtree *rtree)
    rect3 rtree_get_bbox(const rtree *rtree)
    size_t rtree_get_num_leaf_nodes(const rtree *rtree)
    bool rtree_query_bbox(const rtree *rtree, const rect3 *bbox)
    void rtree_intersect(const rtree *rtree, const ray3 *ray, isect *isect)
    void rtree_intersectN(const rtree *rtree, const ray3 *ray, size_t n, isect *isects)
