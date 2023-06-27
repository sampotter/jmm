#pragma once

#include "bb.h"
#include "common.h"
#include "geom.h"
#include "jet.h"

struct bmesh33_cell {
  bb33 const *bb;
  bmesh33_s const *bmesh;
  mesh3_s const *mesh;
  size_t l;
  dbl level;
};

bool bmesh33_cell_intersect(bmesh33_cell_s const *cell, ray3 const *ray, dbl *t);
void bmesh33_cell_Df(bmesh33_cell_s const *cell, dbl3 const x, dbl3 Df);
bool bmesh33_cell_equal(bmesh33_cell_s const *c1, bmesh33_cell_s const *c2);

void bmesh33_alloc(bmesh33_s **bmesh);
void bmesh33_dealloc(bmesh33_s **bmesh);
void bmesh33_init_from_mesh3_and_jets(bmesh33_s *bmesh, mesh3_s const *mesh,
                                      jet31t const *jet);
void bmesh33_deinit(bmesh33_s *bmesh);
size_t bmesh33_num_cells(bmesh33_s const *bmesh);
dbl bmesh33_get_level(bmesh33_s const *bmesh);
mesh3_s const *bmesh33_get_mesh_ptr(bmesh33_s const *bmesh);
bmesh33_s *bmesh33_restrict_to_level(bmesh33_s const *bmesh, dbl level);
bmesh33_cell_s bmesh33_get_cell(bmesh33_s const *bmesh, size_t l);
dbl bmesh33_f(bmesh33_s const *bmesh, dbl3 const x);
bb33 *bmesh33_get_bb_ptr(bmesh33_s const *bmesh, size_t lc);
