#include "edgemap.h"

#include <assert.h>
#include <stdlib.h>

#include "alist.h"

// TODO: starting out with a dumb and simple implementation of this
// just to get going. We can easily optimize this by sorting these
// arrays or using hash maps later. Don't want to waste time on this
// until we can actually build the shadow and reflection cutsets
// correctly...

#define DEFAULT_CAPACITY 32

struct edgemap {
  size_t eltsize;

  // An alist whose keys are size_t and whose values are alist_s *,
  // which are alists with size_t keys and then values which are
  // whatever the element type of the edgemap is. So, basically, this
  // is a deque with two fixed levels, and where the keys are unsorted
  // (for now).
  alist_s *nodes;
};

struct edgemap_iter {
  edgemap_s const *edgemap;
  size_t i, j;
  alist_s *node;
};

void edgemap_iter_alloc(edgemap_iter_s **iter) {
  *iter = malloc(sizeof(edgemap_iter_s));
}

void edgemap_iter_dealloc(edgemap_iter_s **iter) {
  assert(*iter != NULL);
  free(*iter);
  *iter = NULL;
}

void edgemap_iter_init(edgemap_iter_s *iter, edgemap_s const *edgemap) {
  iter->edgemap = edgemap;
  iter->i = iter->j = 0;
  alist_get_by_index(iter->edgemap->nodes, iter->i, &iter->node);
}

bool edgemap_iter_next(edgemap_iter_s *iter, edge_s *edge, void *elt) {
  if (iter->i == alist_size(iter->edgemap->nodes))
    return false;

  if (iter->j == alist_size(iter->node)) {
    iter->j = 0;
    return alist_get_by_index(iter->edgemap->nodes, ++iter->i, iter->node)
      && edgemap_iter_next(iter, edge, elt);
  }

  alist_get_key(iter->edgemap->nodes, iter->i, &edge->l[0]);
  alist_get_pair(iter->node, iter->j++, &edge->l[1], elt);

  return iter->j < alist_size(iter->node)
    || (iter->j == alist_size(iter->node) &&
        iter->i < alist_size(iter->edgemap->nodes));
}

void edgemap_alloc(edgemap_s **edgemap) {
  *edgemap = malloc(sizeof(edgemap_s));
}

void edgemap_dealloc(edgemap_s **edgemap) {
  assert(*edgemap != NULL);
  free(*edgemap);
  *edgemap = NULL;
}

void edgemap_init(edgemap_s *edgemap, size_t eltsize) {
  edgemap->eltsize = eltsize;
  alist_alloc(&edgemap->nodes);
  alist_init(edgemap->nodes, sizeof(size_t), sizeof(alist_s *),
             DEFAULT_CAPACITY);
}

static void deinit_node(edgemap_s *edgemap, size_t i) {
  alist_s *node;
  alist_get_by_index(edgemap->nodes, i, &node);
  alist_deinit(node);
  alist_dealloc(&node);
}

void edgemap_deinit(edgemap_s *edgemap) {
  size_t num_nodes = alist_size(edgemap->nodes);
  for (size_t i = 0; i < num_nodes; ++i)
    deinit_node(edgemap, i);
  alist_deinit(edgemap->nodes);
  alist_dealloc(&edgemap->nodes);
}

bool edgemap_is_empty(edgemap_s const *edgemap) {
  return alist_is_empty(edgemap->nodes);
}

bool edgemap_get(edgemap_s const *edgemap, edge_s edge, void *elt) {
  alist_s *node;
  if (!alist_contains(edgemap->nodes, &edge.l[0]))
    return false;
  alist_get_by_key(edgemap->nodes, &edge.l[0], &node);
  if (!alist_contains(node, &edge.l[1]))
    return false;
  alist_get_by_key(node, &edge.l[1], elt);
  return true;
}

void insert_node(edgemap_s *edgemap, size_t l0) {
  assert(!alist_contains(edgemap->nodes, &l0));
  alist_s *node;
  alist_alloc(&node);
  alist_init(node, sizeof(size_t), edgemap->eltsize, DEFAULT_CAPACITY);
  alist_append(edgemap->nodes, &l0, &node);
}

void edgemap_set(edgemap_s *edgemap, edge_s edge, void const *elt) {
  alist_s *node;
  if (!alist_contains(edgemap->nodes, &edge.l[0]))
    insert_node(edgemap, edge.l[0]);
  alist_get_by_key(edgemap->nodes, &edge.l[0], &node);
  alist_set_by_key(node, &edge.l[1], elt);
}

bool edgemap_contains(edgemap_s const *edgemap, edge_s edge) {
  alist_s *node;
  if (!alist_contains(edgemap->nodes, &edge.l[0]))
    return false;
  alist_get_by_key(edgemap->nodes, &edge.l[0], &node);
  return alist_contains(node, &edge.l[1]);
}

size_t edgemap_size(edgemap_s const *edgemap) {
  size_t size = 0;
  alist_s *node;
  for (size_t i = 0; i < alist_size(edgemap->nodes); ++i) {
    alist_get_by_index(edgemap->nodes, i, &node);
    size += alist_size(node);
  }
  return size;
}

void edgemap_filter(edgemap_s const *edgemap, edgemap_s *edgemap_out,
                    edgemap_prop_t keep, void const *aux) {
  edgemap_iter_s *iter;
  edgemap_iter_alloc(&iter);
  edgemap_iter_init(iter, edgemap);

  edge_s edge;
  void *elt = NULL;
  while (edgemap_iter_next(iter, &edge, elt))
    if (keep(edge, elt, aux))
      edgemap_set(edgemap_out, edge, elt);

  edgemap_iter_dealloc(&iter);
}
