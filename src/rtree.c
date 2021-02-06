#include "rtree.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "def.h"
#include "log.h"
#include "mesh2.h"
#include "stats.h"

#define NUM_BINS 256
#define THRESH 128

struct rtree {
  mesh2_s const *mesh;
  rtree_node_s *root;
};

typedef enum rtree_node_type {
  RTREE_INTERNAL_NODE,
  RTREE_LEAF_NODE
} rtree_node_type_e;

struct rtree_node {
  rect3 bbox;
  rtree_node_type_e type;
  union {
    rtree_node_s *children[2];
    struct {
      size_t *face_inds;
      size_t num_faces;
    } leaf_data;
  };
};

static
void split_node_along_axis(mesh2_s const *mesh, size_t num_faces,
						   size_t const *face_inds, dbl const *centroids,
						   int isplit, rtree_node_s **children) {
  // Compute the mean and standard deviation of the centroids
  // components along the split direction.
  runstd_s runstd;
  runstd_init(&runstd);
  for (size_t i = 0; i < num_faces; ++i) {
    runstd_update(&runstd, centroids[3*i + isplit]);
  }
  dbl mu = runstd_get_mean(&runstd);
  dbl sigma = runstd_get_std(&runstd);

  // Compute a rough approximation of the median by binning (based on
  // the ideas in R. Tibshirani's "Fast computation of the median by
  // successive binning").
  dbl binwidth = 2*sigma/NUM_BINS;
  size_t bins[NUM_BINS];
  size_t bincount = 0;
  memset(bins, 0x0, NUM_BINS*sizeof(size_t));
  for (size_t i = 0; i < num_faces; ++i) {
    dbl c = centroids[3*i + isplit] - mu;
    if (c < -sigma || c >= sigma) {
      continue;
    }
    int k = floor((c + sigma)/binwidth);
    if (k < 0 || k >= NUM_BINS) {
      continue;
    }
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

  // Find the number of triangles in each part of the split.
  size_t count = 0;
  for (size_t i = 0; i < num_faces; ++i) {
    if (centroids[3*i + isplit] <= binmedian) {
      ++count;
    }
  }

  rtree_node_s *child = NULL;

  child = children[0] = malloc(sizeof(rtree_node_s));
  child->type = RTREE_LEAF_NODE;
  child->leaf_data.num_faces = count;
  child->leaf_data.face_inds = malloc(sizeof(size_t)*child->leaf_data.num_faces);

  child = children[1] = malloc(sizeof(rtree_node_s));
  child->type = RTREE_LEAF_NODE;
  child->leaf_data.num_faces = num_faces - count;
  child->leaf_data.face_inds = malloc(sizeof(size_t)*child->leaf_data.num_faces);

  size_t j[2] = {0, 0};
  for (size_t i = 0, argj; i < num_faces; ++i) {
    argj = (size_t)(centroids[3*i + isplit] > binmedian);
    child = children[argj];
    child->leaf_data.face_inds[j[argj]] = face_inds[i];
    child->leaf_data.num_faces = ++j[argj];
  }

  dbl vert[3];
  for (size_t ch = 0; ch < 2; ++ch) {
    child = children[ch];
    child->bbox = rect3_make_empty();
    for (size_t i = 0; i < child->leaf_data.num_faces; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        mesh2_get_vertex(mesh, child->leaf_data.face_inds[i], j, vert);
        rect3_insert_point(&child->bbox, vert);
      }
    }
  }
}

static void refine_node(mesh2_s const *mesh, rtree_node_s *node) {
  size_t *face_inds = node->leaf_data.face_inds;
  size_t num_faces = node->leaf_data.num_faces;

  // Only refine this node if it contains more than THRESH
  // triangles. Otherwise, make this an internal node.
  if (num_faces <= THRESH) {
    return;
  }
  node->type = RTREE_INTERNAL_NODE;

  // Compute centroids for use in finding splits and determining
  // membership of each triangle.
  dbl *centroids = malloc(3*sizeof(dbl)*num_faces);
  for (size_t i = 0; i < num_faces; ++i) {
    mesh2_get_centroid(mesh, face_inds[i], &centroids[3*i]);
  }

  // Do split along each axis and find the SAH for each split.
  dbl SAH[3];
  rtree_node_s *children[3][2];
  for (size_t isplit = 0; isplit < 3; ++isplit) {
    for (size_t ch = 0; ch < 2; ++ch) {
      children[isplit][ch] = malloc(sizeof(rtree_node_s));
    }
    split_node_along_axis(
      mesh, num_faces, face_inds, centroids, isplit, children[isplit]);
    SAH[isplit] = rect3_surface_area(&children[isplit][0]->bbox)
      + rect3_surface_area(&children[isplit][1]->bbox);
  }

  // Find the best split according to the SAH
  size_t isplit;
  dbl min_SAH = INFINITY;
  for (size_t i = 0; i < 3; ++i) {
    if (SAH[i] < min_SAH) {
      min_SAH = SAH[i];
      isplit = i;
    }
  }

  // Free data for unused splits
  for (size_t i = 0; i < 3; ++i) {
    if (i == isplit) {
      continue;
    } else {
      for (size_t ch = 0; ch < 2; ++ch) {
        free(children[i][ch]->leaf_data.face_inds);
        free(children[i][ch]);
      }
    }
  }

  free(centroids);
  free(face_inds);

  for (size_t ch = 0; ch < 2; ++ch) {
    node->children[ch] = children[isplit][ch];
    refine_node(mesh, node->children[ch]);
  }
}

void rtree_alloc(rtree_s **rtree) {
  *rtree = malloc(sizeof(rtree_s));
}

void rtree_dealloc(rtree_s **rtree) {
  free(*rtree);
  *rtree = NULL;
}

void rtree_init_from_tri_mesh(rtree_s *rtree, mesh2_s const *mesh) {
  rtree->mesh = mesh;

  // Initially set root node to be a leaf node containing the entire
  // triangle mesh.
  rtree_node_s *root = rtree->root = malloc(sizeof(rtree_node_s));
  root->bbox = mesh2_get_bounding_box(mesh);
  root->type = RTREE_LEAF_NODE;
  root->leaf_data.num_faces = mesh2_get_num_faces(mesh);
  root->leaf_data.face_inds = malloc(sizeof(size_t)*root->leaf_data.num_faces);
  for (size_t i = 0; i < root->leaf_data.num_faces; ++i) {
    root->leaf_data.face_inds[i] = i;
  }

  // Build the R-tree by refining the leaf nodes until hitting the
  // size threshold.
  refine_node(mesh, root);

  log_info("finished building rtree");
  log_info("* %llu leaves", rtree_get_num_leaf_nodes(rtree));
}

void rtree_deinit(rtree_s *rtree) {
  // TODO: ...
  (void)rtree;
  log_warn("rtree_deinit is a no-op right now! implement it!");
}

rect3 rtree_get_bbox(rtree_s const *rtree) {
  return rtree->root->bbox;
}

static void count_leaf_nodes(rtree_node_s *node, size_t *num_leaf_nodes) {
  if (node->type == RTREE_LEAF_NODE) {
    ++*num_leaf_nodes;
  } else {
    count_leaf_nodes(node->children[0], num_leaf_nodes);
    count_leaf_nodes(node->children[1], num_leaf_nodes);
  }
}

size_t rtree_get_num_leaf_nodes(rtree_s const *rtree) {
  size_t num_leaf_nodes = 0;
  count_leaf_nodes(rtree->root, &num_leaf_nodes);
  return num_leaf_nodes;
}

static bool query_bbox(rtree_node_s const *node, mesh2_s const *mesh,
                       rect3 const *bbox) {
  if (!rect3_overlaps(&node->bbox, bbox)) {
    return false;
  }
  if (node->type == RTREE_LEAF_NODE) {
    for (size_t i = 0; i < node->leaf_data.num_faces; ++i) {
      size_t f = node->leaf_data.face_inds[i];
      if (mesh2_tri_bbox_overlap(mesh, f, bbox)) {
        return true;
      }
    }
    return false;
  } else {
    return query_bbox(node->children[0], mesh, bbox)
      || query_bbox(node->children[1], mesh, bbox);
  }
}

bool rtree_query_bbox(rtree_s const *rtree, rect3 const *bbox) {
  return query_bbox(rtree->root, rtree->mesh, bbox);
}

static bool intersect(rtree_s const *rtree, rtree_node_s const *node,
                      ray3 const *ray, isect_s *isect) {
  if (!ray3_intersects_rect3(ray, &node->bbox))
    return false;
  if (node->type == RTREE_LEAF_NODE) {
    typedef struct {
      dbl t;
      size_t i;
    } hit;
    hit *hits = malloc(node->leaf_data.num_faces*sizeof(hit));
    int num_hits = 0;
    for (size_t i = 0; i < node->leaf_data.num_faces; ++i) {
      size_t f = node->leaf_data.face_inds[i];
      tri3 tri = mesh2_get_tri(rtree->mesh, f);
      if (ray3_intersects_tri3(ray, &tri, &hits[num_hits].t))
        hits[num_hits++].i = f;
    }
    if (num_hits == 0) {
      free(hits);
      return false;
    } else {
      hit closest_hit = {.t = INFINITY};
      for (int i = 0; i < num_hits; ++i)
        if (hits[i].t < closest_hit.t)
          closest_hit = hits[i];
      isect->ray = ray;
      isect->t = closest_hit.t;
      *(size_t *)isect->obj = closest_hit.i;
      free(hits);
      return true;
    }
  } else {
    return intersect(rtree, node->children[0], ray, isect)
      || intersect(rtree, node->children[1], ray, isect);
  }
}

bool rtree_intersect(rtree_s const *rtree, ray3 const *ray, isect_s *isect) {
  return intersect(rtree, rtree->root, ray, isect);
}
