#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "jet.h"

typedef struct mesh3 mesh3_s;

typedef struct bmesh33 bmesh33_s;

void bmesh33_alloc(bmesh33_s **bmesh);
void bmesh33_dealloc(bmesh33_s **bmesh);
void bmesh33_init_from_mesh3_and_jets(bmesh33_s *bmesh, mesh3_s const *mesh,
                                      jet3 const *jet);
void bmesh33_deinit(bmesh33_s *bmesh);
size_t bmesh33_num_cells(bmesh33_s const *bmesh);
mesh3_s const *bmesh33_get_mesh_ptr(bmesh33_s const *bmesh);
bmesh33_s *bmesh33_get_level_bmesh(bmesh33_s const *bmesh, dbl level);

#ifdef __cplusplus
}
#endif
