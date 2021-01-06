#include "grid3.h"

#include <assert.h>
#include <math.h>

grid3_s grid3_restrict_to_rect(grid3_s const *grid, rect3 const *rect, int *offset) {
  dbl h = grid->h;

  grid3_s subgrid;

  int imin[3] = {
    floor((rect->min[0] - grid->min[0])/h),
    floor((rect->min[1] - grid->min[1])/h),
    floor((rect->min[2] - grid->min[2])/h)
  };

  subgrid.min[0] = grid->min[0] + h*imin[0];
  subgrid.min[1] = grid->min[1] + h*imin[1];
  subgrid.min[2] = grid->min[2] + h*imin[2];

  assert(subgrid.min[0] <= rect->min[0]);
  assert(subgrid.min[1] <= rect->min[1]);
  assert(subgrid.min[2] <= rect->min[2]);

  int imax[3] = {
    ceil((rect->max[0] - grid->min[0])/h),
    ceil((rect->max[1] - grid->min[1])/h),
    ceil((rect->max[2] - grid->min[2])/h)
  };

  assert(rect->max[0] <= grid->min[0] + h*imax[0]);
  assert(rect->max[1] <= grid->min[1] + h*imax[1]);
  assert(rect->max[2] <= grid->min[2] + h*imax[2]);

  subgrid.dim[0] = imax[0] - imin[0] + 1;
  subgrid.dim[1] = imax[1] - imin[1] + 1;
  subgrid.dim[2] = imax[2] - imin[2] + 1;

  subgrid.h = h;

  if (offset != NULL) {
    offset[0] = imin[0];
    offset[1] = imin[1];
    offset[2] = imin[2];
  }

  return subgrid;
}

void grid3_get_point(grid3_s const *grid, int const ind[3], dbl point[3]) {
  point[0] = grid->min[0] + grid->h*ind[0];
  point[1] = grid->min[1] + grid->h*ind[1];
  point[2] = grid->min[2] + grid->h*ind[2];
}

void grid3_map(grid3_s const *grid, grid3_map_func_t func, void *ptr) {
  int ind[3];
  for (ind[0] = 0; ind[0] < grid->dim[0]; ++ind[0]) {
    for (ind[1] = 0; ind[1] < grid->dim[1]; ++ind[1]) {
      for (ind[2] = 0; ind[2] < grid->dim[2]; ++ind[2]) {
        func(ind, ptr);
      }
    }
  }
}

size_t grid3_size(grid3_s const *grid) {
  return grid->dim[0]*grid->dim[1]*grid->dim[2];
}

bool grid3_inbounds(grid3_s const *grid, int const ind[3]) {
  return 0 <= ind[0] && ind[0] < grid->dim[0]
    && 0 <= ind[1] && ind[1] < grid->dim[1]
    && 0 <= ind[2] && ind[2] < grid->dim[2];
}
