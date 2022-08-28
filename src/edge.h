#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "def.h"

typedef struct edge {
  size_t l[2];
} edge_s;

edge_s make_edge(size_t l0, size_t l1);
void edge_get_x0_and_dx(edge_s edge, mesh3_s const *mesh, dbl x0[3], dbl dx[3]);
void edge_get_xt(edge_s edge, mesh3_s const *mesh, dbl t, dbl xt[3]);

#ifdef __cplusplus
}
#endif
