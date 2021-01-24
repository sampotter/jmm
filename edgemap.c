#include "edgemap.h"

#include <assert.h>
#include <stdlib.h>

#include "alist.h"
#include "macros.h"

edge_s make_edge(size_t l0, size_t l1) {
  return (edge_s) {.l = {MIN(l0, l1), MAX(l0, l1)}};
}

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

void edgemap_alloc(edgemap_s **edgemap) {
  *edgemap = malloc(sizeof(edgemap_s));
}

void edgemap_dealloc(edgemap_s **edgemap) {
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

void edgemap_get(edgemap_s const *edgemap, edge_s edge, void *elt) {
  alist_s *node;
  if (!alist_contains(edgemap->nodes, &edge.l[0])) return;
  alist_get_by_key(edgemap->nodes, &edge.l[0], &node);
  if (!alist_contains(edgemap->nodes, &edge.l[1])) return;
  alist_get_by_key(node, &edge.l[1], elt);
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
