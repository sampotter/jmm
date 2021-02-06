#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "geom.h"

typedef struct mesh2 mesh2_s;

void mesh2_alloc(mesh2_s **mesh);
void mesh2_dealloc(mesh2_s **mesh);
void mesh2_init_from_binary_files(mesh2_s *mesh, char const *verts_path,
                                  char const *faces_path);
void mesh2_deinit(mesh2_s *mesh);
size_t mesh2_get_num_points(mesh2_s const *mesh);
dbl *mesh2_get_points_ptr(mesh2_s const *mesh);
size_t mesh2_get_num_faces(mesh2_s const *mesh);
size_t *mesh2_get_faces_ptr(mesh2_s const *mesh);
rect3 mesh2_get_bounding_box(mesh2_s const *mesh);
void mesh2_get_centroid(mesh2_s const *mesh, size_t i, dbl *centroid);
void mesh2_get_vertex(mesh2_s const *mesh, size_t i, size_t j, dbl *v);
bool mesh2_tri_bbox_overlap(mesh2_s const *mesh, size_t i, rect3 const *bbox);
tri3 mesh2_get_tri(mesh2_s const *mesh, size_t i);

#ifdef __cplusplus
}
#endif
