#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "bb.h"
#include "common.h"
#include "geom.h"
#include "jet.h"

typedef struct {
  bb33 const *bb;
  mesh3_s const *mesh;
  size_t l;
} bmesh33_cell;

bool bmesh33_cell_ray_intersects_level(bmesh33_cell const *cell,
                                       ray3 const *ray, dbl level, dbl b[4]);

typedef struct bmesh33 bmesh33_s;

void bmesh33_alloc(bmesh33_s **bmesh);
void bmesh33_dealloc(bmesh33_s **bmesh);
void bmesh33_init_from_mesh3_and_jets(bmesh33_s *bmesh, mesh3_s const *mesh,
                                      jet3 const *jet);
void bmesh33_deinit(bmesh33_s *bmesh);
size_t bmesh33_num_cells(bmesh33_s const *bmesh);
mesh3_s const *bmesh33_get_mesh_ptr(bmesh33_s const *bmesh);
bmesh33_s *bmesh33_get_level_bmesh(bmesh33_s const *bmesh, dbl level);
bmesh33_cell bmesh33_get_cell(bmesh33_s const *bmesh, size_t l);

#ifdef __cplusplus
}
#endif
