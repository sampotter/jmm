#pragma once

#ifdef __cplusplus
extern "C" {
#endif

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
                dbl3 const *verts, size_t nverts, bool copy_verts,
                size_t const (*faces)[3], size_t nfaces,
                dbl3 const *face_normals);
void mesh2_init_from_binary_files(mesh2_s *mesh, char const *verts_path,
                                  char const *faces_path);
void mesh2_deinit(mesh2_s *mesh);
void mesh2_dump_verts(mesh2_s const *mesh, char const *path);
void mesh2_dump_faces(mesh2_s const *mesh, char const *path);
size_t mesh2_get_num_points(mesh2_s const *mesh);
dbl3 const *mesh2_get_points_ptr(mesh2_s const *mesh);
size_t mesh2_get_num_faces(mesh2_s const *mesh);
uint3 const *mesh2_get_faces_ptr(mesh2_s const *mesh);
rect3 mesh2_get_bounding_box(mesh2_s const *mesh);
void mesh2_get_centroid(mesh2_s const *mesh, size_t i, dbl *centroid);
void mesh2_get_vertex(mesh2_s const *mesh, size_t i, size_t j, dbl *v);
bool mesh2_tri_bbox_overlap(mesh2_s const *mesh, size_t i, rect3 const *bbox);
tri3 mesh2_get_tri(mesh2_s const *mesh, size_t i);
void mesh2_get_unit_surface_normal(mesh2_s const *mesh, size_t lf, dbl3 n);

#ifdef __cplusplus
}
#endif
