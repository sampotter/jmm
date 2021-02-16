#include "rtree.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gc/gc.h>

#include "def.h"
#include "log.h"
#include "macros.h"
#include "mesh2.h"
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

static void mesh2_tri_insert_into_bbox(mesh2_tri_s const *mesh_tri, rect3 *bbox) {
  rect3_insert_mesh2_tri(bbox, mesh_tri);
}

static void tri3_insert_into_bbox(tri3 const *tri, rect3 *bbox) {
  rect3_insert_tri3(bbox, tri);
}

static void tetra3_insert_into_bbox(tetra3 const *tetra, rect3 *bbox) {
  rect3_insert_tetra3(bbox, tetra);
}

typedef void (*robj_insert_into_bbox_t)(void const *, rect3 *);

robj_insert_into_bbox_t _robj_insert_into_bbox[] = {
  (robj_insert_into_bbox_t)mesh2_tri_insert_into_bbox,
  (robj_insert_into_bbox_t)tri3_insert_into_bbox,
  (robj_insert_into_bbox_t)tetra3_insert_into_bbox
};

void robj_insert_into_bbox(robj_s const *obj, rect3 *bbox) {
  _robj_insert_into_bbox[obj->type](obj->data, bbox);
}

static void mesh2_tri_get_centroid(mesh2_tri_s const *mesh_tri, dbl c[3]) {
  mesh2_get_centroid(mesh_tri->mesh, mesh_tri->l, c);
}

typedef void (*robj_get_centroid_t)(robj_s const *, dbl[3]);

robj_get_centroid_t _robj_get_centroid[] = {
  (robj_get_centroid_t)mesh2_tri_get_centroid,
  (robj_get_centroid_t)tri3_get_centroid,
  (robj_get_centroid_t)tetra3_get_centroid
};

void robj_get_centroid(robj_s const *obj, dbl c[3]) {
  _robj_get_centroid[obj->type](obj->data, c);
}

static bool mesh2_tri_isects_bbox(mesh2_tri_s const *mesh_tri, rect3 const *bbox) {
  dbl center[3], half[3];
  rect3_get_centroid(bbox, center);
  rect3_get_half_extent(bbox, half);
  tri3 tri = mesh2_get_tri(mesh_tri->mesh, mesh_tri->l);
  return triBoxOverlap(center, half, tri.v);
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
  assert(false);
  return false;
}

typedef bool (*robj_isects_bbox_t)(void const *, rect3 const *bbox);

robj_isects_bbox_t _robj_isects_bbox[] = {
  (robj_isects_bbox_t)mesh2_tri_isects_bbox,
  (robj_isects_bbox_t)tri3_isects_bbox,
  (robj_isects_bbox_t)tetra3_isects_bbox
};

bool robj_isects_bbox(robj_s const *obj, rect3 const *bbox) {
  return _robj_isects_bbox[obj->type](obj->data, bbox);
}

static bool mesh2_tri_intersect(mesh2_tri_s const *mesh_tri, ray3 const *ray,
                                dbl *t) {
  tri3 tri = mesh2_get_tri(mesh_tri->mesh, mesh_tri->l);
  return ray3_intersects_tri3(ray, &tri, t);
}

static bool tri3_intersect(tri3 const *tri, ray3 const *ray, dbl *t) {
  return ray3_intersects_tri3(ray, tri, t);
}

static bool tetra3_intersect(tetra3 const *tetra, ray3 const *ray, dbl *t) {
  return ray3_intersects_tetra3(ray, tetra, t);
}

typedef bool (*robj_intersect_t)(void const *, ray3 const *, dbl *);

robj_intersect_t _robj_intersect[] = {
  (robj_intersect_t)mesh2_tri_intersect,
  (robj_intersect_t)tri3_intersect,
  (robj_intersect_t)tetra3_intersect
};

bool robj_intersect(robj_s const *obj, ray3 const *ray, dbl *t) {
  return _robj_intersect[obj->type](obj->data, ray, t);
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
  if (node->type == RNODE_TYPE_INTERNAL) {
    for (int ch = 0; ch < 2; ++ch) {
      rnode_deinit(node->child[ch]);
      node->child[ch] = NULL;
    }
  } else if (node->type == RNODE_TYPE_LEAF) {
    free(node->leaf_data.obj);
    node->leaf_data.size = 0;
    node->leaf_data.capacity = 0;
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
  assert(node->leaf_data.capacity > 0);
  node->leaf_data.capacity *= 2;
  node->leaf_data.obj = realloc(
    node->leaf_data.obj, node->leaf_data.capacity*sizeof(robj_s));
}

void rnode_append_robj(rnode_s *node, robj_s obj) {
  assert(node->type == RNODE_TYPE_LEAF);
  if (node->leaf_data.size == node->leaf_data.capacity)
    rnode_grow(node);
  node->leaf_data.obj[node->leaf_data.size++] = obj;
}

void rnode_append_robjs(rnode_s *node, robj_s const *obj, size_t n) {
  assert(node->type == RNODE_TYPE_LEAF);
  while (node->leaf_data.size + n > node->leaf_data.capacity)
    rnode_grow(node);
  memcpy(&node->leaf_data.obj[node->leaf_data.size], obj, n*sizeof(robj_s));
  node->leaf_data.size += n;
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

// TODO: better to do this using a stack instead of recursively
static bool rnode_intersect(rnode_s const *node, ray3 const *ray, isect *isect) {
  // TODO: remove (redundant)...
  if (!rect3_occludes_ray3(&node->bbox, ray))
    return false;

  // If the node is a leaf node, intersect each of the contained
  // objects.
  if (node->type == RNODE_TYPE_LEAF) {
    dbl t;
    robj_s const *obj;
    for (size_t i = 0; i < node->leaf_data.size; ++i) {
      obj = &node->leaf_data.obj[i];
      if (robj_intersect(obj, ray, &t) && 0 <= t && t < isect->t) {
        isect->t = t;
        isect->obj = obj;
      }
    }
    return isfinite(isect->t);
  }

  // Intersect the bounding boxes of the two child nodes.
  dbl t[2] = {INFINITY, INFINITY};
  ray3_intersects_rect3(ray, &node->child[0]->bbox, &t[0]);
  ray3_intersects_rect3(ray, &node->child[1]->bbox, &t[1]);

  // Sort the pair of intersection parameters, and grab pointers to
  // the child nodes that are sorted into the same order.
  rnode_s *child[2] = {node->child[0], node->child[1]};
  if (t[1] < t[0]) {
    SWAP(t[0], t[1]);
    SWAP(child[0], child[1]);
  }

  // Intersect the child nodes in sorted order, returning early if we
  // find an intersection. Finally, return false if we don't manage to
  // intersect anything.
  for (int i = 0; i < 2; ++i)
    if (isfinite(t[i]) && rnode_intersect(child[i], ray, isect))
      return true;
  return false;
}

/**
 * The number of bins to use when using one pass of Tibshirani's
 * binmedian algorithm to approximate the median along each axis when
 * splitting a node.
 */
#define NUM_BINS 256

// TODO: factor this out as a particular "split strategy"
static
void rnode_split_surface_area(rnode_s const *node,
                                   dbl const (*p)[3], int d,
                                   rnode_s *child[2]) {
  size_t leaf_size = rnode_leaf_size(node);

  // Compute the mean and standard deviation of the centroids
  // components along the split direction.
  runstd_s runstd;
  runstd_init(&runstd);
  for (size_t i = 0; i < leaf_size; ++i)
    runstd_update(&runstd, p[i][d]);
  dbl mu = runstd_get_mean(&runstd);
  dbl sigma = runstd_get_std(&runstd);

  // Compute a rough approximation of the median by binning (based on
  // the ideas in R. Tibshirani's "Fast computation of the median by
  // successive binning").
  dbl binwidth = 2*sigma/NUM_BINS;
  size_t bins[NUM_BINS];
  size_t bincount = 0;
  memset(bins, 0x0, NUM_BINS*sizeof(size_t));
  for (size_t i = 0; i < leaf_size; ++i) {
    dbl c = p[i][d] - mu;
    if (c < -sigma || c >= sigma)
      continue;
    int k = floor((c + sigma)/binwidth);
    if (k < 0 || k >= NUM_BINS)
      continue;
    ++bins[k];
    ++bincount;
  }
  dbl binmedian;
  size_t cumsum = 0;
  for (size_t k = 0; k < NUM_BINS; ++k) {
    cumsum += bins[k];
    if (cumsum >= bincount/2) {
      binmedian = mu - sigma + binwidth*(k + 0.5);
      break;
    }
  }

  // Find the number of triangles in each part of the split so that we
  // can allocate space for them
  size_t child_size[2] = {[0] = 0};
  for (size_t i = 0; i < leaf_size; ++i)
    if (p[i][d] <= binmedian)
      ++child_size[0];
  child_size[1] = leaf_size - child_size[0];

  // Allocate space for the two child nodes
  for (int i = 0; i < 2; ++i) {
    child[i] = malloc(sizeof(rnode_s));
    child[i]->type = RNODE_TYPE_LEAF;
    child[i]->leaf_data.size = child_size[i];
    child[i]->leaf_data.obj = malloc(child_size[i]*sizeof(robj_s));
  }

  // Traverse the parent node objects and separate them based on which
  // side of the binmedian their centroid is on
  size_t j[2] = {0, 0};
  for (size_t i = 0; i < leaf_size; ++i) {
    int c = (size_t)(p[i][d] > binmedian);
    child[c]->leaf_data.obj[j[c]++] = node->leaf_data.obj[i];
  }

  // Recompute the bounding boxes of the children. This is fine to do
  // now, because when we split these nodes, their bounding boxes
  // won't change.
  rnode_recompute_bbox(child[0]);
  rnode_recompute_bbox(child[1]);
}

// Section: rtree_s

struct rtree {
  rnode_s root;
  size_t leaf_thresh; // Maximum size of a leaf node
  rtree_split_strategy_e split_strategy;
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
}

void rtree_deinit(rtree_s *rtree) {
  rnode_deinit(&rtree->root);
}

void rtree_insert_mesh2(rtree_s *rtree, mesh2_s const *mesh) {
  rnode_s *node = &rtree->root;
  assert(node->type == RNODE_TYPE_LEAF);
  size_t num_faces = mesh2_get_num_faces(mesh);
  robj_s obj = {.type = ROBJ_MESH2_TRI};
  for (size_t l = 0; l < num_faces; ++l) {
    mesh2_tri_s *tri = GC_MALLOC(sizeof(mesh2_tri_s));
    tri->mesh = mesh;
    tri->l = l;
    obj.data = tri;
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
  node->type = RNODE_TYPE_INTERNAL;

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
    rnode_split_surface_area(node, p, d, child[d]);
    for (int i = 0; i < 2; ++i)
      surf_area[d] += rect3_surface_area(&child[d][i]->bbox);
  }

  // The best split is the one with the minimum sum of surface areas
  // (i.e., using the "surface area heuristic")
  int dmin;
  dbl min_surf_area = INFINITY;
  for (int d = 0; d < 3; ++d) {
    if (surf_area[d] < min_surf_area) {
      min_surf_area = surf_area[d];
      dmin = d;
    }
  }

  // Free data for unused splits
  for (int d = 0; d < 3; ++d) {
    if (d == dmin) continue;
    for (int i = 0; i < 2; ++i) {
      rnode_deinit(child[d][i]);
      free(child[d][i]);
    }
  }
  free(p); // Free centroids

  // Set children and recursively refine them.
  for (int i = 0; i < 2; ++i) {
    node->child[i] = child[dmin][i];
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

void rtree_intersect(rtree_s const *rtree, ray3 const *ray, isect *isect) {
  isect->t = INFINITY;
  if (!rnode_intersect(&rtree->root, ray, isect))
    isect->t = NAN;
}

void rtree_intersectN(rtree_s const *rtree, ray3 const *ray, size_t n,
                      isect *isect) {
  for (size_t i = 0; i < n; ++i)
    rtree_intersect(rtree, &ray[i], &isect[i]);
}
