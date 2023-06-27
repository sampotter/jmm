#pragma once

#include "common.h"
#include "geom.h"

struct mesh2_tri {
  mesh2_s const *mesh;
  size_t l; // index of triangle
};

bool mesh2_tri_equal(mesh2_tri_s const *t1, mesh2_tri_s const *t2);

void mesh2_alloc(mesh2_s **mesh);
void mesh2_dealloc(mesh2_s **mesh);
void mesh2_init(mesh2_s *mesh,
                dbl3 const *verts, size_t nverts, policy_e verts_policy,
                size_t const (*faces)[3], size_t nfaces, policy_e faces_policy,
                dbl3 const *face_normals, policy_e face_normals_policy);
void mesh2_init_from_binary_files(mesh2_s *mesh, char const *verts_path,
                                  char const *faces_path);
void mesh2_deinit(mesh2_s *mesh);
void mesh2_dump_verts(mesh2_s const *mesh, char const *path);
void mesh2_dump_faces(mesh2_s const *mesh, char const *path);
size_t mesh2_nverts(mesh2_s const *mesh);
size_t mesh2_nfaces(mesh2_s const *mesh);
dbl3 const *mesh2_get_verts_ptr(mesh2_s const *mesh);
uint3 const *mesh2_get_faces_ptr(mesh2_s const *mesh);

size_t mesh2_find_face(mesh2_s const *mesh, uint3 const vf);

size_t mesh2_nvf(mesh2_s const *mesh, size_t l);
void mesh2_vf(mesh2_s const *mesh, size_t l, size_t *lf);
void mesh2_fv(mesh2_s const *mesh, size_t lf, uint3 l);
size_t mesh2_nve(mesh2_s const *mesh, size_t l);
void mesh2_ve(mesh2_s const *mesh, size_t l, uint2 *le);
void mesh2_fve(mesh2_s const *mesh, size_t lf, size_t l, uint2 le);
void mesh2_ef(mesh2_s const *mesh, uint2 const le, size_t lf[2]);
size_t mesh2_fvf(mesh2_s const *mesh, size_t lf, size_t l);

rect3 mesh2_get_bounding_box(mesh2_s const *mesh);
void mesh2_get_centroid(mesh2_s const *mesh, size_t i, dbl *centroid);
void mesh2_get_vertex(mesh2_s const *mesh, size_t i, size_t j, dbl *v);
bool mesh2_tri_bbox_overlap(mesh2_s const *mesh, size_t i, rect3 const *bbox);
tri3 mesh2_get_tri(mesh2_s const *mesh, size_t i);
void mesh2_get_unit_surface_normal(mesh2_s const *mesh, size_t lf, dbl3 n);
void mesh2_get_R_for_face(mesh2_s const *mesh, size_t lf, dbl33 R);
