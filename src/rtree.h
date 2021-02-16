#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "geom.h"

typedef struct mesh2 mesh2_s;

typedef enum robj_type {
  ROBJ_MESH2_TRI,
  ROBJ_TRI3,
  ROBJ_TETRA3
} robj_type_e;

typedef struct robj robj_s;

robj_type_e robj_get_type(robj_s const *obj);
void const *robj_get_data(robj_s const *obj);
void robj_insert_into_bbox(robj_s const *obj, rect3 *bbox);
void robj_get_centroid(robj_s const *obj, dbl c[3]);
bool robj_isects_bbox(robj_s const *obj, rect3 const *bbox);
bool robj_intersect(robj_s const *obj, ray3 const *ray, dbl *t);

typedef struct {
  /**
   * The intersection parameter: i.e., a scalar such that location of
   * the intersection is `ray->org + t*ray->dir`.
   */
  dbl t;

  /**
   * A pointer to the object intersecting the ray.
   */
  robj_s const *obj;
} isect;

typedef enum rtree_split_strategy {
  RTREE_SPLIT_STRATEGY_SURFACE_AREA
} rtree_split_strategy_e;

/**
 * A "static" R-tree, in the sense that all contained objects must be
 * inserted before constructing the R-tree. Once the R-tree is
 * constructed, it is fixed and cannot have objects added to it or
 * removed from it.
 */
typedef struct rtree rtree_s;

void rtree_alloc(rtree_s **rtree);
void rtree_dealloc(rtree_s **rtree);
void rtree_init(rtree_s *rtree, size_t leaf_thresh,
                rtree_split_strategy_e split_strategy);
void rtree_deinit(rtree_s *rtree);
void rtree_insert_mesh2(rtree_s *rtree, mesh2_s const *mesh);
void rtree_build(rtree_s *rtree);
rect3 rtree_get_bbox(rtree_s const *rtree);
size_t rtree_get_num_leaf_nodes(rtree_s const *rtree);
bool rtree_query_bbox(rtree_s const *rtree, rect3 const *bbox);
void rtree_intersect(rtree_s const *rtree, ray3 const *ray, isect *isect);
void rtree_intersectN(rtree_s const *rtree, ray3 const *ray, size_t n,
                      isect *isects);

#ifdef __cplusplus
}
#endif
