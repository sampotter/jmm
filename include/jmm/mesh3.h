#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "geom.h"
#include "index.h"
#include "par.h"
#include "vec.h"

bool face_in_cell(size_t const f[3], size_t const c[4]);
bool point_in_face(size_t l, size_t const f[3]);
bool point_in_cell(size_t l, size_t const c[4]);
bool edge_in_face(size_t const le[2], size_t const lf[3]);

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

struct mesh3_tetra {
  mesh3_s const *mesh;
  size_t l; // index of tetrahedron
};

tri3 mesh3_tetra_get_face(mesh3_tetra_s const *tetra, int face_inds[3]);
bool mesh3_tetra_equal(mesh3_tetra_s const *t1, mesh3_tetra_s const *t2);

typedef struct mesh3_data {
  size_t nverts;
  dbl3 *verts;
  size_t ncells;
  uint4 *cells;
} mesh3_data_s;

void mesh3_data_init_from_bin(mesh3_data_s *data, char const *verts_path, char const *cells_path);
void mesh3_data_init_from_off_file(mesh3_data_s *data, char const *path, dbl maxvol, bool verbose);
error_e mesh3_data_insert_vert(mesh3_data_s *data, dbl3 const x, dbl eps);

void mesh3_alloc(mesh3_s **mesh);
void mesh3_dealloc(mesh3_s **mesh);
void mesh3_init(mesh3_s *mesh, mesh3_data_s const *data, bool compute_bd_info, dbl const *eps);
void mesh3_deinit(mesh3_s *mesh);
dbl3 const *mesh3_get_verts_ptr(mesh3_s const *mesh);
size_t const *mesh3_get_cells_ptr(mesh3_s const *mesh);
dbl const *mesh3_get_vert_ptr(mesh3_s const *mesh, size_t i);
void mesh3_get_vert_ptrs(mesh3_s const *mesh, size_t const *l, int n, dbl const **x);
void mesh3_copy_vert(mesh3_s const *mesh, size_t i, dbl *v);
tetra3 mesh3_get_tetra(mesh3_s const *mesh, size_t lc);
tri3 mesh3_get_tri(mesh3_s const *mesh, size_t const lf[3]);
void mesh3_get_centroid(mesh3_s const *mesh, size_t lc, dbl centroid[3]);
void mesh3_get_edge_centroid(mesh3_s const *mesh, size_t e[2], dbl c[3]);
size_t mesh3_ncells(mesh3_s const *mesh);
size_t mesh3_nverts(mesh3_s const *mesh);
void mesh3_get_bbox(mesh3_s const *mesh, rect3 *bbox);
void mesh3_get_cell_bbox(mesh3_s const *mesh, size_t i, rect3 *bbox);
bool mesh3_cell_contains_point(mesh3_s const *mesh, size_t i, dbl const x[3]);
bool mesh3_contains_ball(mesh3_s const *mesh, dbl3 const x, dbl r);
size_t mesh3_find_cell_containing_point(mesh3_s const *mesh, dbl const x[3], size_t lc);
bool mesh3_contains_point(mesh3_s const *mesh, dbl3 const x);
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
void mesh3_cf(mesh3_s const *mesh, size_t lc, size_t lf[4][3]);
void mesh3_cv(mesh3_s const *mesh, size_t i, size_t *cv);
int mesh3_nec(mesh3_s const *mesh, size_t const le[2]);
void mesh3_ec(mesh3_s const *mesh, size_t const le[2], size_t *lc);
int mesh3_nee(mesh3_s const *mesh, size_t const le[2]);
void mesh3_ee(mesh3_s const *mesh, size_t const e[2], size_t (*ee)[2]);
size_t mesh3_nev(mesh3_s const *mesh, size_t const e[2]);
void mesh3_ev(mesh3_s const *mesh, size_t const e[2], size_t *v);
int mesh3_nfc(mesh3_s const *mesh, size_t const f[3]);
void mesh3_fc(mesh3_s const *mesh, size_t const f[3], size_t *fc);
bool mesh3_cee(mesh3_s const *mesh, size_t c, size_t const e[2],
               size_t e_out[2]);
bool mesh3_cfv(mesh3_s const *mesh, size_t lc, size_t const lf[3], size_t *lv);
bool mesh3_ccfv(mesh3_s const *mesh, size_t lc, size_t const lf[3],
                size_t *lv_out);
bool mesh3_cvf(mesh3_s const *mesh, size_t lc, size_t lv, size_t lf[3]);
bool mesh3_has_bd_info(mesh3_s const *mesh);
bool *mesh3_get_bdc_ptr(mesh3_s *mesh);
bool mesh3_bdc(mesh3_s const *mesh, size_t i);
bool *mesh3_get_bdv_ptr(mesh3_s *mesh);
bool mesh3_bdv(mesh3_s const *mesh, size_t i);
size_t mesh3_nbde(mesh3_s const *mesh);
bool mesh3_bde(mesh3_s const *mesh, size_t const l[2]);
size_t mesh3_nbdf(mesh3_s const *mesh);
void mesh3_get_bdf_inds(mesh3_s const *mesh, size_t l, size_t lf[3]);
void mesh3_set_bdf(mesh3_s *mesh, size_t lf[3]);
bool mesh3_is_bdf(mesh3_s const *mesh, size_t const lf[3]);
bool mesh3_is_edge(mesh3_s const *mesh, size_t const l[2]);
bool mesh3_is_diff_edge(mesh3_s const *mesh, size_t const l[2]);
bool mesh3_vert_incident_on_diff_edge(mesh3_s const *mesh, size_t l);
bool mesh3_vert_is_terminal_diff_edge_vert(mesh3_s const *mesh, size_t l);
bool mesh3_local_ray_in_vertex_cone(mesh3_s const *mesh, dbl3 const p, size_t lv);
bool mesh3_ray_prop_from_edge_is_occluded(mesh3_s const *mesh, dbl3 const t,
                                          uint2 const le);
bool mesh3_local_ray_is_occluded(mesh3_s const *mesh, size_t lhat, par3_s const *par);
bool mesh3_cell_incident_on_diff_edge(mesh3_s const *mesh, size_t lc);
dbl mesh3_get_min_tetra_alt(mesh3_s const *mesh);
dbl mesh3_get_min_edge_length(mesh3_s const *mesh);
dbl mesh3_get_mean_edge_length(mesh3_s const *mesh);
mesh2_s *mesh3_get_surface_mesh(mesh3_s const *mesh);
size_t mesh3_get_num_inc_diff_edges(mesh3_s const *mesh, size_t l);
void mesh3_get_inc_diff_edges(mesh3_s const *mesh, size_t l, size_t (*le)[2]);
size_t mesh3_get_num_inc_bdf(mesh3_s const *mesh, size_t l);
void mesh3_get_inc_bdf(mesh3_s const *mesh, size_t l, size_t (*lf)[3]);
size_t mesh3_get_num_reflectors(mesh3_s const *mesh);
size_t mesh3_get_reflector_size(mesh3_s const *mesh, size_t i);
void mesh3_get_reflector(mesh3_s const *mesh, size_t i, size_t (*lf)[3]);
mesh2_s *mesh3_get_refl_mesh(mesh3_s const *mesh, size_t i);
size_t mesh3_get_num_diffractors(mesh3_s const *mesh);
size_t mesh3_get_diffractor_size(mesh3_s const *mesh, size_t i);
void mesh3_get_diffractor(mesh3_s const *mesh, size_t i, size_t (*le)[2]);
mesh1_s *mesh3_get_diff_mesh(mesh3_s const *mesh, size_t diff_index);
size_t mesh3_get_num_diffs_inc_on_refl(mesh3_s const *mesh, size_t refl_index);
void mesh3_get_diffs_inc_on_refl(mesh3_s const *mesh, size_t refl_index, size_t *diff_index);
size_t mesh3_get_diff_index_for_bde(mesh3_s const *mesh, size_t const le[2]);
size_t mesh3_nbde(mesh3_s const *mesh);
void mesh3_get_bde_inds(mesh3_s const *mesh, size_t l, size_t le[2]);
void mesh3_set_bde(mesh3_s *mesh, size_t const le[2], bool diff);
dbl mesh3_get_eps(mesh3_s const *mesh);
void mesh3_get_face_normal(mesh3_s const *mesh, size_t const lf[3], dbl normal[3]);
void mesh3_get_R_for_face(mesh3_s const *mesh, size_t const lf[3], dbl33 R);
void mesh3_get_R_for_reflector(mesh3_s const *mesh, size_t refl_index, dbl33 R);
void mesh3_get_R_for_interior_reflector_vertex(mesh3_s const *mesh, size_t l, dbl33 R);
void mesh3_get_diff_edge_tangent(mesh3_s const *mesh, size_t const le[2],
                                 dbl t[3]);
dbl mesh3_get_edge_ext_angle(mesh3_s const *mesh, size_t const le[2]);
bool mesh3_edge_contains_point(mesh3_s const *mesh, size_t const le[2], dbl const x[3]);
size_t mesh3_get_num_bdf_inc_on_edge(mesh3_s const *mesh, size_t const le[2]);
void mesh3_get_bdf_inc_on_edge(mesh3_s const *mesh, size_t const le[2], size_t (*lf)[3]);
void mesh3_get_edge_tangent(mesh3_s const *mesh, size_t const le[2], dbl t[3]);
void mesh3_get_edge_midpoint(mesh3_s const *mesh, size_t const le[2], dbl p[3]);
void mesh3_get_face_centroid(mesh3_s const *mesh, size_t const lf[3], dbl p[3]);
void mesh3_dump_verts(mesh3_s const *mesh, char const *path);
void mesh3_dump_cells(mesh3_s const *mesh, char const *path);
bool mesh3_has_vertex(mesh3_s const *mesh, dbl3 const x);
size_t mesh3_get_vert_index(mesh3_s const *mesh, dbl3 const x);
dbl mesh3_linterp(mesh3_s const *mesh, dbl const *values, dbl3 const x);
dbl mesh3_diam_2approx(mesh3_s const *mesh, size_t l);
dbl mesh3_diam_2approx_rand(mesh3_s const *mesh, size_t trials, size_t const *seed);
dbl mesh3_get_diam(mesh3_s const *mesh);
dbl mesh3_get_edge_tol(mesh3_s const *mesh, uint2 const le);
dbl mesh3_get_face_tol(mesh3_s const *mesh, uint3 const lf);

#ifdef __cplusplus
}
#endif
