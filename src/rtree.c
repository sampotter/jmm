#include "rtree.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "bmesh.h"
#include "def.h"
#include "log.h"
#include "macros.h"
#include "mesh2.h"
#include "mesh3.h"
#include "pool.h"
#include "stats.h"
#include "vec.h"

// Section: robj_s

struct robj {
  robj_type_e type;
  void const *data;
};

robj_type_e robj_get_type(robj_s const *obj) {
  return obj->type;
}

void const *robj_get_data(robj_s const *obj) {
  return obj->data;
}

static void bmesh33_cell_insert_into_bbox(bmesh33_cell_s const *cell, rect3 *bbox) {
  tetra3 tetra = mesh3_get_tetra(cell->mesh, cell->l);
  rect3_insert_tetra3(bbox, &tetra);
}

static void mesh2_tri_insert_into_bbox(mesh2_tri_s const *mesh_tri, rect3 *bbox) {
  rect3_insert_mesh2_tri(bbox, mesh_tri);
}

static void mesh3_tetra_insert_into_bbox(mesh3_tetra_s const *tetra, rect3 *bbox) {
  rect3_insert_mesh3_tetra(bbox, tetra);
}

static void tri3_insert_into_bbox(tri3 const *tri, rect3 *bbox) {
  rect3_insert_tri3(bbox, tri);
}

static void tetra3_insert_into_bbox(tetra3 const *tetra, rect3 *bbox) {
  rect3_insert_tetra3(bbox, tetra);
}

typedef void (*robj_insert_into_bbox_t)(void const *, rect3 *);

robj_insert_into_bbox_t _robj_insert_into_bbox[] = {
  (robj_insert_into_bbox_t)bmesh33_cell_insert_into_bbox,
  (robj_insert_into_bbox_t)mesh2_tri_insert_into_bbox,
  (robj_insert_into_bbox_t)mesh3_tetra_insert_into_bbox,
  (robj_insert_into_bbox_t)tri3_insert_into_bbox,
  (robj_insert_into_bbox_t)tetra3_insert_into_bbox
};

void robj_insert_into_bbox(robj_s const *obj, rect3 *bbox) {
  _robj_insert_into_bbox[obj->type](obj->data, bbox);
}

static void bmesh33_cell_get_centroid(bmesh33_cell_s const *cell, dbl c[3]) {
  mesh3_get_centroid(cell->mesh, cell->l, c);
}

static void mesh2_tri_get_centroid(mesh2_tri_s const *mesh_tri, dbl c[3]) {
  mesh2_get_centroid(mesh_tri->mesh, mesh_tri->l, c);
}

static void mesh3_tetra_get_centroid(mesh3_tetra_s const *tetra, dbl c[3]) {
  mesh3_get_centroid(tetra->mesh, tetra->l, c);
}

typedef void (*robj_get_centroid_t)(robj_s const *, dbl[3]);

robj_get_centroid_t _robj_get_centroid[] = {
  (robj_get_centroid_t)bmesh33_cell_get_centroid,
  (robj_get_centroid_t)mesh2_tri_get_centroid,
  (robj_get_centroid_t)mesh3_tetra_get_centroid,
  (robj_get_centroid_t)tri3_get_centroid,
  (robj_get_centroid_t)tetra3_get_centroid
};

void robj_get_centroid(robj_s const *obj, dbl c[3]) {
  _robj_get_centroid[obj->type](obj->data, c);
}

static bool bmesh33_cell_isects_bbox(bmesh33_cell_s const *cell, rect3 const *bbox) {
  (void)cell;
  (void)bbox;
  die();
  return false;
}

static bool mesh2_tri_isects_bbox(mesh2_tri_s const *mesh_tri, rect3 const *bbox) {
  dbl center[3], half[3];
  rect3_get_centroid(bbox, center);
  rect3_get_half_extent(bbox, half);
  tri3 tri = mesh2_get_tri(mesh_tri->mesh, mesh_tri->l);
  return triBoxOverlap(center, half, tri.v);
}

static bool mesh3_tetra_isects_bbox(mesh3_tetra_s const *tetra, rect3 const *bbox) {
  (void)tetra;
  (void)bbox;
  die();
  return false;
}

static bool tri3_isects_bbox(tri3 const *tri, rect3 const *bbox) {
  dbl center[3], half[3];
  rect3_get_centroid(bbox, center);
  rect3_get_half_extent(bbox, half);
  return triBoxOverlap(center, half, tri->v);
}

static bool tetra3_isects_bbox(tetra3 const *tetra, rect3 const *bbox) {
  (void)tetra;
  (void)bbox;
  die();
  return false;
}

typedef bool (*robj_isects_bbox_t)(void const *, rect3 const *bbox);

robj_isects_bbox_t _robj_isects_bbox[] = {
  (robj_isects_bbox_t)bmesh33_cell_isects_bbox,
  (robj_isects_bbox_t)mesh2_tri_isects_bbox,
  (robj_isects_bbox_t)mesh3_tetra_isects_bbox,
  (robj_isects_bbox_t)tri3_isects_bbox,
  (robj_isects_bbox_t)tetra3_isects_bbox
};

bool robj_isects_bbox(robj_s const *obj, rect3 const *bbox) {
  return _robj_isects_bbox[obj->type](obj->data, bbox);
}

/** robj_intersect: */

static bool mesh2_tri_intersect(mesh2_tri_s const *mesh_tri, ray3 const *ray,
                                dbl *t) {
  tri3 tri = mesh2_get_tri(mesh_tri->mesh, mesh_tri->l);
  return ray3_intersects_tri3(ray, &tri, t);
}

static bool mesh3_tetra_intersect(mesh3_tetra_s const *mesh_tetra,
                                  ray3 const *ray, dbl *t) {
  tetra3 tetra = mesh3_get_tetra(mesh_tetra->mesh, mesh_tetra->l);
  return ray3_intersects_tetra3(ray, &tetra, t);
}

static bool tri3_intersect(tri3 const *tri, ray3 const *ray, dbl *t) {
  return ray3_intersects_tri3(ray, tri, t);
}

static bool tetra3_intersect(tetra3 const *tetra, ray3 const *ray, dbl *t) {
  return ray3_intersects_tetra3(ray, tetra, t);
}

typedef bool (*robj_intersect_t)(void const *, ray3 const *, dbl *);

robj_intersect_t _robj_intersect[] = {
  (robj_intersect_t)bmesh33_cell_intersect,
  (robj_intersect_t)mesh2_tri_intersect,
  (robj_intersect_t)mesh3_tetra_intersect,
  (robj_intersect_t)tri3_intersect,
  (robj_intersect_t)tetra3_intersect
};

bool robj_intersect(robj_s const *obj, ray3 const *ray, dbl *t) {
  return _robj_intersect[obj->type](obj->data, ray, t);
}

/** robj_equal: */

typedef bool (*robj_equal_t)(void const *, void const *);

robj_equal_t _robj_equal[] = {
  (robj_equal_t)bmesh33_cell_equal,
  (robj_equal_t)mesh2_tri_equal,
  (robj_equal_t)mesh3_tetra_equal,
  (robj_equal_t)tri3_equal,
  (robj_equal_t)tetra3_equal
};

bool robj_equal(robj_s const *obj1, robj_s const *obj2) {
  if (obj1 == NULL || obj2 == NULL)
    return false;

  if (obj1->type != obj2->type)
    return false;

  return _robj_equal[obj1->type](obj1->data, obj2->data);
}

// Section: rnode_s

typedef enum rnode_type {
  RNODE_TYPE_INTERNAL,
  RNODE_TYPE_LEAF
} rnode_type_e;

/**
 * The default capacity of an R-tree node.
 */
#define RNODE_DEFAULT_CAPACITY 32

typedef struct rnode {
  rect3 bbox;
  rnode_type_e type;
  union {
    struct rnode *child[2];
    struct {
      robj_s *obj;
      size_t size;
      size_t capacity;
    } leaf_data;
  };
} rnode_s;

void rnode_init(rnode_s *node, rnode_type_e type) {
  node->bbox = rect3_make_empty();
  node->type = type;
  if (type == RNODE_TYPE_INTERNAL) {
    node->child[0] = NULL;
    node->child[1] = NULL;
  } else if (type == RNODE_TYPE_LEAF) {
    node->leaf_data.obj = malloc(RNODE_DEFAULT_CAPACITY*sizeof(robj_s));
    node->leaf_data.size = 0;
    node->leaf_data.capacity = RNODE_DEFAULT_CAPACITY;
  }
}

void rnode_deinit(rnode_s *node) {
  if (node == NULL)
    return;
  if (node->type == RNODE_TYPE_INTERNAL) {
    for (int ch = 0; ch < 2; ++ch) {
      rnode_deinit(node->child[ch]);
      free(node->child[ch]);
      node->child[ch] = NULL;
    }
  } else if (node->type == RNODE_TYPE_LEAF) {
    free(node->leaf_data.obj);
    node->leaf_data.obj = NULL;
    node->leaf_data.size = 0;
    node->leaf_data.capacity = 0;
  }
}

bool rnode_empty(rnode_s const *node) {
  assert(node->type == RNODE_TYPE_LEAF);
  return node->leaf_data.size == 0;
}

void rnode_copy_deep(rnode_s const *node, rnode_s *copy) {
  copy->bbox = node->bbox;
  copy->type = node->type;
  if (node->type == RNODE_TYPE_INTERNAL) {
    for (int i = 0; i < 2; ++i) {
      copy->child[i] = malloc(sizeof(rnode_s));
      rnode_copy_deep(node->child[i], copy->child[i]);
    }
  } else if (node->type == RNODE_TYPE_LEAF) {
    size_t size = node->leaf_data.size;
    size_t capacity = node->leaf_data.capacity;
    copy->leaf_data.obj = malloc(capacity*sizeof(robj_s));
    memcpy(copy->leaf_data.obj, node->leaf_data.obj, size*sizeof(robj_s));
    copy->leaf_data.size = size;
    copy->leaf_data.capacity = capacity;
  }
}

void rnode_recompute_bbox(rnode_s *node) {
  assert(node->type == RNODE_TYPE_LEAF);
  node->bbox = rect3_make_empty();
  for (size_t i = 0; i < node->leaf_data.size; ++i)
    robj_insert_into_bbox(&node->leaf_data.obj[i], &node->bbox);
}

void rnode_grow(rnode_s *node) {
  assert(node->type == RNODE_TYPE_LEAF);
  assert(0 < node->leaf_data.capacity);
  assert(node->leaf_data.size <= node->leaf_data.capacity);

  node->leaf_data.capacity *= 2;

  node->leaf_data.obj = realloc(
    node->leaf_data.obj, node->leaf_data.capacity*sizeof(robj_s));
}

void rnode_append_robj(rnode_s *node, robj_s obj) {
  assert(node->type == RNODE_TYPE_LEAF);
  if (node->leaf_data.size == node->leaf_data.capacity)
    rnode_grow(node);
  node->leaf_data.obj[node->leaf_data.size++] = obj;
  assert(node->leaf_data.size <= node->leaf_data.capacity);
}

void rnode_append_robjs(rnode_s *node, robj_s const *obj, size_t n) {
  assert(node->type == RNODE_TYPE_LEAF);
  while (node->leaf_data.size + n > node->leaf_data.capacity)
    rnode_grow(node);
  memcpy(&node->leaf_data.obj[node->leaf_data.size], obj, n*sizeof(robj_s));
  node->leaf_data.size += n;
  assert(node->leaf_data.size <= node->leaf_data.capacity);
}

size_t rnode_leaf_size(rnode_s const *node) {
  assert(node->type == RNODE_TYPE_LEAF);
  return node->leaf_data.size;
}

// TODO: better to do this using a stack instead of recursively
static bool rnode_query_bbox(rnode_s const *node, rect3 const *bbox) {
  if (!rect3_overlaps(&node->bbox, bbox))
    return false;
  if (node->type == RNODE_TYPE_INTERNAL)
    return rnode_query_bbox(node->child[0], bbox)
        || rnode_query_bbox(node->child[1], bbox);
  for (size_t i = 0; i < node->leaf_data.size; ++i)
    if (robj_isects_bbox(&node->leaf_data.obj[i], bbox))
      return true;
  return false;
}

static void
rnode_intersect_leaf(rnode_s const *node, ray3 const *ray, isect *isect,
                     robj_s const *skip_robj) {
  dbl t;
  robj_s const *obj;
  for (size_t i = 0; i < node->leaf_data.size; ++i) {
    obj = &node->leaf_data.obj[i];
    if (robj_equal(skip_robj, obj))
      continue;
    if (robj_intersect(obj, ray, &t) && 0 <= t && t < isect->t) {
      isect->t = t;
      isect->obj = obj;
    }
  }
}

static void
rnode_intersect(rnode_s const *node, ray3 const *ray, isect *isect,
                robj_s const *skip_robj) {
  /* First, check whether the ray intersects the current node's
   * bounding box. If it doesn't, we don't want to update the
   * intersection's t value or check for intersections below this
   * node, so we return early. */
  dbl t = ray3_intersect_rect3(ray, &node->bbox);
  if (isinf(t))
    return;

  /* If this is a leaf node, then we've bottomed out in the recursion:
   * time to check for intersections with each contained object. */
  if (node->type == RNODE_TYPE_LEAF) {
    rnode_intersect_leaf(node, ray, isect, skip_robj);
    return;
  }

  /* Intersect the bounding boxes of the two child nodes. */
  dbl t_child[2] = {
    ray3_intersect_rect3(ray, &node->child[0]->bbox),
    ray3_intersect_rect3(ray, &node->child[1]->bbox)
  };

  /* Sort the pair of intersection parameters, and grab pointers to
   * the child nodes in the same order. */
  rnode_s *child[2] = {node->child[0], node->child[1]};
  if (t_child[1] < t_child[0]) {
    SWAP(t_child[0], t_child[1]);
    SWAP(child[0], child[1]);
  }

  /* Intersect the child nodes in sorted order. */
  for (int i = 0; i < 2; ++i)
    if (isfinite(t_child[i]))
      rnode_intersect(child[i], ray, isect, skip_robj);
}

/**
 * The number of bins to use when using one pass of Tibshirani's
 * binmedian algorithm to approximate the median along each axis when
 * splitting a node.
 */
#define NUM_BINS 256

// TODO: factor this out as a particular "split strategy"
static
bool rnode_split_surface_area(rnode_s const *node, dbl const (*p)[3], int d,
                              rnode_s *child[2]) {
  size_t leaf_size = rnode_leaf_size(node);

  dbl *pd = malloc(leaf_size*sizeof(dbl));
  for (size_t i = 0; i < leaf_size; ++i)
    pd[i] = p[i][d];

  dbl med = dblN_median(leaf_size, pd);

  // Separate the objects in the leaf node being split into the two
  // children depending on which side of the median the centroid of
  // each objects lies.
  for (size_t i = 0; i < leaf_size; ++i) {
    int which = med < pd[i];
    rnode_append_robj(child[which], node->leaf_data.obj[i]);
  }

  // Check whether either of the children are empty at this point. If
  // they are, return false to signal that we shouldn't use this
  // split.
  bool success = true;
  if (rnode_empty(child[0]) || rnode_empty(child[1])) {
    success = false;
    goto cleanup;
  }

  // Recompute the bounding boxes of the children. This is fine to do
  // now, because when we split these nodes, their bounding boxes
  // won't change.
  for (int which = 0; which < 2; ++which)
    rnode_recompute_bbox(child[which]);

cleanup:
  free(pd);

  return success;
}

// Section: rtree_s

#define RTREE_POOL_INITIAL_CAPACITY 4096

struct rtree {
  rnode_s root;
  size_t leaf_thresh; // Maximum size of a leaf node
  rtree_split_strategy_e split_strategy;
  pool_s *pool;
  bool pool_owner;
};

void rtree_alloc(rtree_s **rtree) {
  *rtree = malloc(sizeof(rtree_s));
}

void rtree_dealloc(rtree_s **rtree) {
  free(*rtree);
  *rtree = NULL;
}

void rtree_init(rtree_s *rtree, size_t leaf_thresh,
                rtree_split_strategy_e split_strategy) {
  rnode_init(&rtree->root, RNODE_TYPE_LEAF);
  rtree->leaf_thresh = leaf_thresh;
  rtree->split_strategy = split_strategy;

  pool_alloc(&rtree->pool);
  pool_init(rtree->pool, RTREE_POOL_INITIAL_CAPACITY);
  rtree->pool_owner = true;
}

void rtree_deinit(rtree_s *rtree) {
  rnode_deinit(&rtree->root);

  if (rtree->pool_owner) {
    pool_deinit(rtree->pool);
    pool_dealloc(&rtree->pool);
  }
}

rtree_s *rtree_copy(rtree_s const *rtree) {
  rtree_s *copy = malloc(sizeof(rtree_s));
  rnode_copy_deep(&rtree->root, &copy->root);
  copy->leaf_thresh = rtree->leaf_thresh;
  copy->split_strategy = rtree->split_strategy;
  copy->pool = rtree->pool;
  copy->pool_owner = false; // The original rtree is responsible for
                            // freeing the pool
  return copy;
}

void rtree_insert_bmesh33(rtree_s *rtree, bmesh33_s const *bmesh) {
  rnode_s *node = &rtree->root;
  assert(node->type == RNODE_TYPE_LEAF);
  size_t num_cells = bmesh33_num_cells(bmesh);
  robj_s obj = {.type = ROBJ_BMESH33_CELL};
  for (size_t l = 0; l < num_cells; ++l) {
    bmesh33_cell_s *cell = pool_get(rtree->pool, sizeof(bmesh33_cell_s));
    *cell = bmesh33_get_cell(bmesh, l);
    obj.data = cell;
    rnode_append_robj(node, obj);
  }
  rnode_recompute_bbox(node);
}

void rtree_insert_mesh2(rtree_s *rtree, mesh2_s const *mesh) {
  rnode_s *node = &rtree->root;
  assert(node->type == RNODE_TYPE_LEAF);
  size_t num_faces = mesh2_nfaces(mesh);
  robj_s obj = {.type = ROBJ_MESH2_TRI};
  for (size_t l = 0; l < num_faces; ++l) {
    mesh2_tri_s *tri = pool_get(rtree->pool, sizeof(mesh2_tri_s));
    tri->mesh = mesh;
    tri->l = l;
    obj.data = tri;
    rnode_append_robj(node, obj);
  }
  rnode_recompute_bbox(node);
}

void rtree_insert_mesh3(rtree_s *rtree, mesh3_s const *mesh) {
  rnode_s *node = &rtree->root;
  assert(node->type == RNODE_TYPE_LEAF);
  size_t num_cells = mesh3_ncells(mesh);
  robj_s obj = {.type = ROBJ_MESH3_TETRA};
  for (size_t l = 0; l < num_cells; ++l) {
    mesh3_tetra_s *tetra = pool_get(rtree->pool, sizeof(mesh3_tetra_s));
    tetra->mesh = mesh;
    tetra->l = l;
    obj.data = tetra;
    rnode_append_robj(node, obj);
  }
  rnode_recompute_bbox(node);
}

static void refine_node_surface_area(rtree_s const *rtree, rnode_s *node) {
  assert(node->type == RNODE_TYPE_LEAF);

  // Refine this leaf node if it contains more than too many
  // objects. Otherwise, make this an internal node and split it
  // below.
  size_t leaf_size = rnode_leaf_size(node);
  if (leaf_size <= rtree->leaf_thresh)
    return;

  // Compute centroids to compute the splits as well as determine the
  // membership of each triangle.
  dbl (*p)[3] = malloc(leaf_size*sizeof(dbl[3]));
  for (size_t i = 0; i < leaf_size; ++i)
    robj_get_centroid(&node->leaf_data.obj[i], p[i]);

  // Create leaf nodes for each split.
  rnode_s *child[3][2];
  for (int d = 0; d < 3; ++d) {
    for (int i = 0; i < 2; ++i) {
      child[d][i] = malloc(sizeof(rnode_s));
      rnode_init(child[d][i], RNODE_TYPE_LEAF);
    }
  }

  // Do split along each axis, computing the surface area of each
  // split along the way.
  dbl surf_area[3] = {0, 0, 0};
  for (int d = 0; d < 3; ++d) {
    if (rnode_split_surface_area(node, p, d, child[d]))
      for (int i = 0; i < 2; ++i)
        surf_area[d] += rect3_surface_area(&child[d][i]->bbox);
  }

  // The best split is the one with the minimum sum of surface areas
  // (i.e., using the "surface area heuristic")
  int dmin = NO_INDEX;
  dbl min_surf_area = INFINITY;
  for (int d = 0; d < 3; ++d) {
    if (0 < surf_area[d] && surf_area[d] < min_surf_area) {
      min_surf_area = surf_area[d];
      dmin = d;
    }
  }
  assert(0 <= dmin && dmin < 3);

  // Free data for unused splits
  for (int d = 0; d < 3; ++d) {
    if (d == dmin) continue;
    for (int i = 0; i < 2; ++i) {
      rnode_deinit(child[d][i]);
      free(child[d][i]);
    }
  }
  free(p); // Free centroids

  // Convert the current node to an internal node.
  node->type = RNODE_TYPE_INTERNAL;

  // Free the data for unsplit leaf node that we're now replacing.
  free(node->leaf_data.obj);

  // Set children and recursively refine them.
  for (int i = 0; i < 2; ++i) {
    node->child[i] = child[dmin][i];
    assert(node->child[i] != NULL);
    refine_node_surface_area(rtree, node->child[i]);
  }
}

typedef void (*refine_node_t)(rtree_s const *, rnode_s *);

refine_node_t _refine_node[] = {
  refine_node_surface_area
};

void refine_node(rtree_s const *rtree, rnode_s *node) {
  _refine_node[rtree->split_strategy](rtree, node);
}

void rtree_build(rtree_s *rtree) {
  refine_node(rtree, &rtree->root);
}

rect3 rtree_get_bbox(rtree_s const *rtree) {
  return rtree->root.bbox;
}

static void count_leaf_nodes(rnode_s const *node, size_t *num_leaf_nodes) {
  if (node->type == RNODE_TYPE_LEAF) {
    ++*num_leaf_nodes;
  } else {
    count_leaf_nodes(node->child[0], num_leaf_nodes);
    count_leaf_nodes(node->child[1], num_leaf_nodes);
  }
}

size_t rtree_get_num_leaf_nodes(rtree_s const *rtree) {
  size_t num_leaf_nodes = 0;
  count_leaf_nodes(&rtree->root, &num_leaf_nodes);
  return num_leaf_nodes;
}

bool rtree_query_bbox(rtree_s const *rtree, rect3 const *bbox) {
  return rnode_query_bbox(&rtree->root, bbox);
}

void rtree_intersect(rtree_s const *rtree, ray3 const *ray, isect *isect,
                     robj_s const *skip_robj) {
  isect->t = INFINITY;
  isect->obj = NULL;
  rnode_intersect(&rtree->root, ray, isect, skip_robj);
}

void rtree_intersectN(rtree_s const *rtree, ray3 const *ray, size_t n,
                      isect *isect) {
  for (size_t i = 0; i < n; ++i)
    rtree_intersect(rtree, &ray[i], &isect[i], NULL);
}
