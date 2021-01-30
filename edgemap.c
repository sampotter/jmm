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
  alist_s *lst;
};

struct edgemap_iter {
  edgemap_s const *edgemap;
  size_t i;
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
  iter->i = 0;
}

bool edgemap_iter_next(edgemap_iter_s *iter, edge_s *edge, void *elt) {
  if (iter->i == alist_size(iter->edgemap->lst))
    return false;
  alist_get_pair(iter->edgemap->lst, iter->i++, edge, elt);
  return true;
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
  alist_alloc(&edgemap->lst);
  alist_init(edgemap->lst, sizeof(edge_s), eltsize, DEFAULT_CAPACITY);
}

void edgemap_deinit(edgemap_s *edgemap) {
  alist_deinit(edgemap->lst);
  alist_dealloc(&edgemap->lst);
}

bool edgemap_is_empty(edgemap_s const *edgemap) {
  return alist_is_empty(edgemap->lst);
}

bool edgemap_get(edgemap_s const *edgemap, edge_s edge, void *elt) {
  return alist_get_by_key(edgemap->lst, &edge, elt);
}

void edgemap_set(edgemap_s *edgemap, edge_s edge, void const *elt) {
  alist_set_by_key(edgemap->lst, &edge, elt);
}

bool edgemap_contains(edgemap_s const *edgemap, edge_s edge) {
  return alist_contains(edgemap->lst, &edge);
}

size_t edgemap_size(edgemap_s const *edgemap) {
  return alist_size(edgemap->lst);
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

void edgemap_clear(edgemap_s *edgemap) {
  alist_clear(edgemap->lst);
}
