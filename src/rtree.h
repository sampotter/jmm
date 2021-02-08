#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "geom.h"

typedef struct mesh2 mesh2_s;

typedef struct {
  /**
   * The intersecting ray.
   */
  ray3 const *ray;

  /**
   * The intersection parameter: i.e., a scalar such that location of
   * the intersection is `ray->org + t*ray->dir`.
   */
  dbl t;

  /**
   * A pointer to the object intersecting the ray.
   */
  void *obj;
} isect;

typedef struct rtree_node rtree_node_s;

typedef struct rtree rtree_s;

void rtree_alloc(rtree_s **rtree);
void rtree_dealloc(rtree_s **rtree);
void rtree_init_from_tri_mesh(rtree_s *rtree, mesh2_s const *mesh);
void rtree_deinit(rtree_s *rtree);
rect3 rtree_get_bbox(rtree_s const *rtree);
size_t rtree_get_num_leaf_nodes(rtree_s const *rtree);
bool rtree_query_bbox(rtree_s const *rtree, rect3 const *bbox);
void rtree_intersect(rtree_s const *rtree, ray3 const *ray, isect *isect);
void rtree_intersectN(rtree_s const *rtree, ray3 const *rays, isect *isects,
                      size_t num_rays);

#ifdef __cplusplus
}
#endif
