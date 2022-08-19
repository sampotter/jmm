#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "def.h"

void mesh1_alloc(mesh1_s **mesh);
void mesh1_dealloc(mesh1_s **mesh);
void mesh1_init(mesh1_s *mesh,
                dbl3 const *verts, size_t nverts, policy_e verts_policy,
                uint2 const *edges, size_t nedges, policy_e edges_policy);
void mesh1_deinit(mesh1_s *mesh);

size_t mesh1_nedges(mesh1_s const *mesh);

void mesh1_ev(mesh1_s const *mesh, size_t le, uint2 l);
void mesh1_ve(mesh1_s const *mesh, size_t l, size_t le[2]);

#ifdef __cplusplus
}
#endif
