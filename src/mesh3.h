#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "geom.h"
#include "index.h"
#include "vec.h"

typedef struct mesh2 mesh2_s;

bool face_in_cell(size_t const f[3], size_t const c[4]);
bool point_in_face(size_t l, size_t const f[3]);
bool point_in_cell(size_t l, size_t const c[4]);

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

typedef struct mesh3_tetra {
  mesh3_s const *mesh;
  size_t l; // index of tetrahedron
} mesh3_tetra_s;

void mesh3_alloc(mesh3_s **mesh);
void mesh3_dealloc(mesh3_s **mesh);
void mesh3_init(mesh3_s *mesh,
                dbl const *verts, size_t nverts,
                size_t const *cells, size_t ncells);
void mesh3_deinit(mesh3_s *mesh);
dbl const *mesh3_get_verts_ptr(mesh3_s const *mesh);
size_t const *mesh3_get_cells_ptr(mesh3_s const *mesh);
dvec3 mesh3_get_vert(mesh3_s const *mesh, size_t i);
dbl const *mesh3_get_vert_ptr(mesh3_s const *mesh, size_t i);
void mesh3_get_vert_ptrs(mesh3_s const *mesh, size_t const *l, int n, dbl const **x);
void mesh3_copy_vert(mesh3_s const *mesh, size_t i, dbl *v);
tetra3 mesh3_get_tetra(mesh3_s const *mesh, size_t lc);
void mesh3_get_centroid(mesh3_s const *mesh, size_t lc, dbl centroid[3]);
size_t mesh3_ncells(mesh3_s const *mesh);
size_t mesh3_nverts(mesh3_s const *mesh);
void mesh3_get_bbox(mesh3_s const *mesh, rect3 *bbox);
void mesh3_get_cell_bbox(mesh3_s const *mesh, size_t i, rect3 *bbox);
bool mesh3_dbl3_in_cell(mesh3_s const *mesh, size_t i, dbl const x[3], dbl b[4]);
int mesh3_nvc(mesh3_s const *mesh, size_t i);
void mesh3_vc(mesh3_s const *mesh, size_t i, size_t *vc);
int mesh3_nve(mesh3_s const *mesh, size_t lv);
void mesh3_ve(mesh3_s const *mesh, size_t lv, size_t (*ve)[2]);
int mesh3_nvf(mesh3_s const *mesh, size_t i);
void mesh3_vf(mesh3_s const *mesh, size_t i, size_t (*vf)[3]);
int mesh3_nvv(mesh3_s const *mesh, size_t i);
void mesh3_vv(mesh3_s const *mesh, size_t i, size_t *vv);
int mesh3_ncc(mesh3_s const *mesh, size_t i);
void mesh3_cc(mesh3_s const *mesh, size_t i, size_t *cc);
void mesh3_cv(mesh3_s const *mesh, size_t i, size_t *cv);
int mesh3_nec(mesh3_s const *mesh, size_t i, size_t j);
void mesh3_ec(mesh3_s const *mesh, size_t i, size_t j, size_t *ec);
int mesh3_nee(mesh3_s const *mesh, size_t const e[2]);
void mesh3_ee(mesh3_s const *mesh, size_t const e[2], size_t (*ee)[2]);
int mesh3_nfc(mesh3_s const *mesh, size_t const f[3]);
void mesh3_fc(mesh3_s const *mesh, size_t const f[3], size_t *fc);
bool mesh3_cee(mesh3_s const *mesh, size_t c, size_t const e[2],
               size_t e_out[2]);
bool mesh3_cfv(mesh3_s const *mesh, size_t lc, size_t const lf[3], size_t *lv);
bool mesh3_ccfv(mesh3_s const *mesh, size_t lc, size_t const lf[3],
                size_t *lv_out);
bool mesh3_cvf(mesh3_s const *mesh, size_t lc, size_t lv, size_t lf[3]);
bool *mesh3_get_bdc_ptr(mesh3_s *mesh);
bool mesh3_bdc(mesh3_s const *mesh, size_t i);
bool *mesh3_get_bdv_ptr(mesh3_s *mesh);
bool mesh3_bdv(mesh3_s const *mesh, size_t i);
bool mesh3_bde(mesh3_s const *mesh, size_t const l[2]);
bool mesh3_bdf(mesh3_s const *mesh, size_t const l[3]);
bool mesh3_is_edge(mesh3_s const *mesh, size_t const l[2]);
bool mesh3_is_diff_edge(mesh3_s const *mesh, size_t const l[2]);
bool mesh3_vert_incident_on_diff_edge(mesh3_s const *mesh, size_t l);
dbl mesh3_get_min_tetra_alt(mesh3_s const *mesh);
mesh2_s *mesh3_get_surface_mesh(mesh3_s const *mesh);

#ifdef __cplusplus
}
#endif
