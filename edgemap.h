#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stddef.h>

#include "edge.h"

typedef struct edgemap edgemap_s;

void edgemap_alloc(edgemap_s **edgemap);
void edgemap_dealloc(edgemap_s **edgemap);
void edgemap_init(edgemap_s *edgemap, size_t eltsize);
void edgemap_deinit(edgemap_s *edgemap);
void edgemap_get(edgemap_s const *edgemap, edge_s edge, void *elt);
void edgemap_set(edgemap_s *edgemap, edge_s edge, void const *elt);
bool edgemap_contains(edgemap_s const *edgemap, edge_s edge);

#ifdef __cplusplus
}
#endif
