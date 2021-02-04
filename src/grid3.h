#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "geom.h"

typedef struct grid3 {
  int dim[3];
  dbl min[3];
  dbl h;
} grid3_s;

typedef void (*grid3_map_func_t)(int const [3], void *);

grid3_s grid3_restrict_to_rect(grid3_s const *grid, rect3 const *rect, int *offset);
void grid3_get_point(grid3_s const *grid, int const ind[3], dbl point[3]);
void grid3_map(grid3_s const *grid, grid3_map_func_t func, void *ptr);
size_t grid3_size(grid3_s const *grid);
bool grid3_inbounds(grid3_s const *grid, int const ind[3]);

#ifdef __cplusplus
}
#endif
