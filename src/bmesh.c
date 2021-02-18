#include "bmesh.h"

#include <assert.h>
#include <stdlib.h>

#include "bb.h"
#include "mesh3.h"

struct bmesh33 {
  mesh3_s const *mesh;
  bool mesh_owner;
  size_t num_cells;
  bb33 *bb;
};

void bmesh33_alloc(bmesh33_s **bmesh) {
  *bmesh = malloc(sizeof(bmesh33_s));
}

void bmesh33_dealloc(bmesh33_s **bmesh) {
  free(*bmesh);
  *bmesh = NULL;
}

void bmesh33_init_from_mesh3_and_jets(bmesh33_s *bmesh, mesh3_s const *mesh,
                                      jet3 const *jet) {
  bmesh->mesh = mesh;
  bmesh->mesh_owner = false;
  bmesh->num_cells = mesh3_ncells(mesh);

  // Interpolate jets to create Bezier tetrahedra for each cell
  bmesh->bb = malloc(bmesh->num_cells*sizeof(bb33));
  for (size_t l = 0; l < bmesh->num_cells; ++l)
    bb33_init_from_cell_and_jets(&bmesh->bb[l], mesh, jet, l);
}

void bmesh33_deinit(bmesh33_s *bmesh) {
  free(bmesh->bb);
  bmesh->bb = NULL;

  if (bmesh->mesh_owner) {
    mesh3_deinit((mesh3_s *)bmesh->mesh);
    mesh3_dealloc((mesh3_s **)&bmesh->mesh);
  }
}

size_t bmesh33_num_cells(bmesh33_s const *bmesh) {
  return bmesh->num_cells;
}

mesh3_s const *bmesh33_get_mesh_ptr(bmesh33_s const *bmesh) {
  return bmesh->mesh;
}

bmesh33_s *bmesh33_get_level_bmesh(bmesh33_s const *bmesh, dbl level) {
  // First, precompute which Bezier tetra *might* bracket the
  // level. Since the graph of each Bezier tetra is contained in the
  // convex hull of the defining control points, we can quickly rule
  // out points outside of this interval (just check if the level is
  // greater than the maximum or smaller than the minimum Bezier
  // ordinate), but we may still have some false positives in what we
  // compute this way.
  size_t num_brack = 0;
  bool *brack = malloc(bmesh->num_cells*sizeof(bool));
  for (size_t l = 0; l < bmesh->num_cells; ++l) {
    brack[l] = bb33_convex_hull_brackets_value(&bmesh->bb[l], level);
    num_brack += brack[l];
  }

  bmesh33_s *level_bmesh = malloc(sizeof(bmesh33_s));
  level_bmesh->bb = malloc(num_brack*sizeof(bb33));

  dbl *verts = malloc(4*num_brack*sizeof(dbl[3]));
  size_t *cells = malloc(num_brack*sizeof(size_t[4]));

  // For each bracketed cell, we go ahead and splat all of the
  // incident vertices into `verts`. There will be duplicates. For now
  // we're trading off memory use for simplicity. Later on we can come
  // back and optimize this. While we're doing this, we also copy over
  // Bezier tetra data so that we don't have to recompute it.
  size_t lv = 0, lc = 0;
  for (size_t l = 0; l < bmesh->num_cells; ++l) {
    if (!brack[l]) continue;
    level_bmesh->bb[lc] = bmesh->bb[l]; // Copy Bezier tetra data
    mesh3_cv(bmesh->mesh, l, &cells[4*lc]);
    for (int i = 0; i < 4; ++i)
      mesh3_copy_vert(bmesh->mesh, cells[4*lc + i], &verts[4*3*lc + 3*i]);
    lv += 4;
    ++lc;
  }

  mesh3_alloc((mesh3_s **)&level_bmesh->mesh);
  mesh3_init((mesh3_s *)level_bmesh->mesh, verts, lv, cells, lc, false);

  level_bmesh->mesh_owner = true;
  level_bmesh->num_cells = mesh3_ncells(level_bmesh->mesh);
  assert(level_bmesh->num_cells == lc);

  free(cells);
  free(verts);

  return level_bmesh;
}
