#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stddef.h>

#include "edge.h"

typedef struct edgemap edgemap_s;
typedef struct edgemap_iter edgemap_iter_s;

typedef bool (* edgemap_prop_t)(edge_s, void const *, void const *);

void edgemap_iter_alloc(edgemap_iter_s **iter);
void edgemap_iter_dealloc(edgemap_iter_s **iter);
void edgemap_iter_init(edgemap_iter_s *iter, edgemap_s const *edgemap);
bool edgemap_iter_next(edgemap_iter_s *iter, edge_s *edge, void *elt);

void edgemap_alloc(edgemap_s **edgemap);
void edgemap_dealloc(edgemap_s **edgemap);
void edgemap_init(edgemap_s *edgemap, size_t eltsize);
void edgemap_deinit(edgemap_s *edgemap);
bool edgemap_is_empty(edgemap_s const *edgemap);
bool edgemap_get(edgemap_s const *edgemap, edge_s edge, void *elt);
void edgemap_set(edgemap_s *edgemap, edge_s edge, void const *elt);
bool edgemap_contains(edgemap_s const *edgemap, edge_s edge);
size_t edgemap_size(edgemap_s const *edgemap);
void edgemap_filter(edgemap_s const *edgemap, edgemap_s *edgemap_out,
                    edgemap_prop_t keep, void const *aux);

#ifdef __cplusplus
}
#endif
