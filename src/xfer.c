#include "xfer.h"

#include <assert.h>

#include "bb.h"
#include "grid3.h"
#include "index.h"
#include "mesh3.h"
#include "vec.h"

typedef struct {
  mesh3_s const *mesh;
  grid3_s const *grid;
  dbl *y;
  grid3_s subgrid;
  int offset[3];
  bb33 bb;
  size_t lc;
} xfer_tetra_wkspc_t;

static void xfer_tetra(int const subgrid_ind[3], xfer_tetra_wkspc_t *wkspc) {
  dbl const atol = 1e-14;

  // Get the Cartesian coordinates of the point in the subgrid indexed
  // by `subgrid_ind`.
  dbl point[3];
  grid3_get_point(&wkspc->subgrid, subgrid_ind, point);

  // Get the indices of this point in the original containing grid.
  ivec3 dim = ivec3_from_int3(wkspc->grid->dim);
  ivec3 grid_ind = ivec3_add(
    ivec3_from_int3(subgrid_ind), ivec3_from_int3(wkspc->offset));

  // Return early if this point is outside the original grid.
  //
  // TODO: could remove this check by ensuring that a generated
  // subgrid is contained in the parent grid.
  if (!grid3_inbounds(wkspc->grid, &grid_ind.data[0])) {
    return;
  }

  /**
   * Check if the point is in cell `lc`, computing its barycentric
   * coordinates along the way.
   */
  dbl b[4];
  tetra3 tetra = mesh3_get_tetra(wkspc->mesh, wkspc->lc);
  if (tetra3_contains_point(&tetra, point, &atol)) {
    tetra3_get_bary_coords(&tetra, point, b);
    size_t l = ind2l3(dim, grid_ind);
    dbl y = bb33_f(&wkspc->bb, b);
    // If the grid value is NaN, just set it. Otherwise, set it to the
    // average of the new value and existing grid value. This is a bit
    // arbitrary, but if we're doing everything else right, it should
    // be okay.
    if (isnan(wkspc->y[l])) {
      wkspc->y[l] = y;
    } else {
      wkspc->y[l] = (y + wkspc->y[l])/2;
    }
  }
}

void xfer(mesh3_s const *mesh, jet3 const *jet, grid3_s const *grid, dbl *y) {
  // Initialize values on grid to NaN.
  for (size_t i = 0; i < grid3_size(grid); ++i) {
    y[i] = NAN;
  }

  /**
   * Traverse the cells of the tetrahedron mesh, restrict the grid so
   * that it just covers each cell, and then evaluate the Bezier
   * tetrahedron at each grid node if it contains the grid node.
   */
  xfer_tetra_wkspc_t wkspc  = {.mesh = mesh, .grid = grid, .y = y};
  for (wkspc.lc = 0; wkspc.lc < mesh3_ncells(mesh); ++wkspc.lc) {
    rect3 bbox;
    mesh3_get_cell_bbox(mesh, wkspc.lc, &bbox);
    wkspc.subgrid = grid3_restrict_to_rect(grid, &bbox, wkspc.offset);
    bb33_init_from_cell_and_jets(&wkspc.bb, mesh, jet, wkspc.lc);
    grid3_map(&wkspc.subgrid, (grid3_map_func_t)xfer_tetra, &wkspc);
  }
}
