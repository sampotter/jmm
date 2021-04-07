#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct eik3 eik3_s;
typedef struct mesh3 mesh3_s;

typedef struct edge {
  size_t l[2];
} edge_s;

edge_s make_edge(size_t l0, size_t l1);
size_t edge_get_valid_index(edge_s edge, eik3_s const *eik);
size_t edge_get_shadow_index(edge_s edge, eik3_s const *eik);
void edge_get_x0_and_dx(edge_s edge, mesh3_s const *mesh, dbl x0[3], dbl dx[3]);
void edge_get_xt(edge_s edge, mesh3_s const *mesh, dbl t, dbl xt[3]);

#ifdef __cplusplus
}
#endif
