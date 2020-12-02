#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

// TODO: make it so that the verts/cells can be copied or not (in
// which case mesh3 is just a "view")

typedef struct mesh3 mesh3_s;

void mesh3_alloc(mesh3_s **mesh);
void mesh3_dealloc(mesh3_s **mesh);
void mesh3_init(mesh3_s *mesh,
                dbl const *verts, size_t nverts,
                size_t const *cells, size_t ncells);
void mesh3_deinit(mesh3_s *mesh);
size_t mesh3_nverts(mesh3_s const *mesh);
int mesh3_nvc(mesh3_s const *mesh, size_t i);
void mesh3_vc(mesh3_s const *mesh, size_t i, size_t *vc);
int mesh3_nvv(mesh3_s const *mesh, size_t i);
void mesh3_vv(mesh3_s const *mesh, size_t i, size_t *vv);
int mesh3_nve(mesh3_s const *mesh, size_t i);
void mesh3_ve(mesh3_s const *mesh, size_t i, size_t *ve);
int mesh3_nvf(mesh3_s const *mesh, size_t i);
void mesh3_vf(mesh3_s const *mesh, size_t i, size_t *vf);
int mesh3_ncc(mesh3_s const *mesh, size_t i);
void mesh3_cc(mesh3_s const *mesh, size_t i, size_t *cc);
bool mesh3_bdc(mesh3_s const *mesh, size_t i);

#ifdef __cplusplus
}
#endif
