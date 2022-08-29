#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "def.h"
#include "vec.h"

void mesh22_alloc(mesh22_s **mesh);
void mesh22_dealloc(mesh22_s **mesh);
void mesh22_init(mesh22_s *mesh, dbl2 const *verts, size_t nverts,
                 uint3 const *faces, size_t nfaces);
void mesh22_init_from_binary_files(mesh22_s *mesh, char const *verts_path,
                                   char const *faces_path);
void mesh22_deinit(mesh22_s *mesh);
size_t mesh22_nverts(mesh22_s const *mesh);
size_t mesh22_nfaces(mesh22_s const *mesh);
dbl2 const *mesh22_get_verts_ptr(mesh22_s const *mesh);
uint3 const *mesh22_get_faces_ptr(mesh22_s const *mesh);
void mesh22_get_vert(mesh22_s const *mesh, size_t l, dbl2 x);
void mesh22_fv(mesh22_s const *mesh, size_t l, size_t fv[3]);
size_t mesh22_nvf(mesh22_s const *mesh, size_t l);
void mesh22_vf(mesh22_s const *mesh, size_t l, size_t *vf);
size_t mesh22_nvv(mesh22_s const *mesh, size_t l);
void mesh22_vv(mesh22_s const *mesh, size_t l, size_t *vv);

#ifdef __cplusplus
}
#endif
