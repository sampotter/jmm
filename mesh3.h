#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "vec.h"

// Some ideas for improving the design of mesh3:
//
// - TODO: make it so that the verts/cells can be copied or not (in
//   which case mesh3 is just a "view")
//
// - TODO: replace e.g. pairs like mesh3_nvc and mesh3_vc with a
//   single function, mesh3_vc, which builds and returns a pointer to
//   an array_s
//
// - TODO: write generator-style functions which will elements from
//   e.g. mesh3_vc one-at-a-time

typedef struct mesh3 mesh3_s;

void mesh3_alloc(mesh3_s **mesh);
void mesh3_dealloc(mesh3_s **mesh);
void mesh3_init(mesh3_s *mesh,
                dbl const *verts, size_t nverts,
                size_t const *cells, size_t ncells);
void mesh3_deinit(mesh3_s *mesh);
dvec3 mesh3_get_vert(mesh3_s const *mesh, size_t i);
dbl const *mesh3_get_vert_ptr(mesh3_s const *mesh, size_t i);
void mesh3_copy_vert(mesh3_s const *mesh, size_t i, dbl *v);
size_t mesh3_ncells(mesh3_s const *mesh);
size_t mesh3_nverts(mesh3_s const *mesh);
int mesh3_nvc(mesh3_s const *mesh, size_t i);
void mesh3_vc(mesh3_s const *mesh, size_t i, size_t *vc);
int mesh3_nvf(mesh3_s const *mesh, size_t i);
void mesh3_vf(mesh3_s const *mesh, size_t i, size_t (*vf)[3]);
int mesh3_nvv(mesh3_s const *mesh, size_t i);
void mesh3_vv(mesh3_s const *mesh, size_t i, size_t *vv);
int mesh3_ncc(mesh3_s const *mesh, size_t i);
void mesh3_cc(mesh3_s const *mesh, size_t i, size_t *cc);
void mesh3_cv(mesh3_s const *mesh, size_t i, size_t *cv);
int mesh3_nec(mesh3_s const *mesh, size_t i, size_t j);
void mesh3_ec(mesh3_s const *mesh, size_t i, size_t j, size_t *ec);
bool *mesh3_get_bdc_ptr(mesh3_s *mesh);
bool mesh3_bdc(mesh3_s const *mesh, size_t i);
bool *mesh3_get_bdv_ptr(mesh3_s *mesh);
bool mesh3_bdv(mesh3_s const *mesh, size_t i);

#ifdef __cplusplus
}
#endif
