#pragma once

#include "common.h"
#include "def.h"

#define GRID2_NUM_NB 8

/* If `order == ORDER_ROW_MAJOR`, then the indexing/coordinate
 * convention that is used is:
 *
 *       +j/+y
 *   o---------->
 *   |
 *   |
 *   |
 *   |+i/+x
 *   |
 *   v
 *
 * Here, the "o" denotes the origin, where (i, j) = (0, 0) and (x, y)
 * = (0, 0). We have shown the coordinate system oriented so that it
 * matches the typical matrix indexing convention for languages that
 * are row major. E.g., a numpy array `arr`, indexed `arr[i, j]`, is
 * laid out as shown on the page. To get the usual x-y plane, we can
 * rotate this picture counterclockwise 90 degrees.
 *
 * Note: `order == ORDER_COLUMN_MAJOR` isn't implemented yet, but the
 * idea is that once it is, it will be straightforward to wrap this
 * library for languages which are column major by convention (MATLAB,
 * R, Fortran, Julia, ...). */
typedef struct grid2 {
  int2 shape;
  dbl2 xymin;
  dbl h;
  order_e order;
} grid2_s;

typedef struct grid2info grid2info_s;

void grid2_save(grid2_s const *grid, char const *path);
size_t grid2_nind(grid2_s const *grid);
size_t grid2_nindc(grid2_s const *grid);
int grid2_ind2l(grid2_s const *grid, int2 const ind);
int grid2_ind2lc(grid2_s const *grid, int2 const ind);
int grid2_indc2l(grid2_s const *grid, int2 const indc);
int grid2_indc2lc(grid2_s const *grid, int2 const indc);
void grid2_l2ind(grid2_s const *grid, int l, int2 ind);
void grid2_l2indc(grid2_s const *grid, int l, int2 indc);
void grid2_lc2ind(grid2_s const *grid, int lc, int2 ind);
void grid2_lc2indc(grid2_s const *grid, int lc, int2 indc);
int grid2_l2lc(grid2_s const *grid, int l);
int grid2_lc2l(grid2_s const *grid, int lc);
void grid2_l2xy(grid2_s const *grid, int l, dbl2 xy);
int grid2_xy2lc(grid2_s const *grid, dbl2 const xy, dbl2 cc);
bool grid2_isind(grid2_s const *grid, int2 ind);
bool grid2_isindc(grid2_s const *grid, int2 indc);
void grid2_get_inbounds(grid2_s const *grid, grid2info_s const *info,
                        int l, bool inbounds[GRID2_NUM_NB + 1]);
void grid2_get_nb(grid2_s const *grid, grid2info_s const *info,
                  int l, int l_nb[GRID2_NUM_NB], bool inbounds[GRID2_NUM_NB]);

struct grid2info {
  int2 offsets[GRID2_NUM_NB + 1];
  int nb_dl[GRID2_NUM_NB + 1];
};

void grid2info_init(grid2info_s *info, grid2_s const *grid);

typedef struct grid2_to_mesh3_mapping {
  /* The source grid for the mapping. */
  grid2_s const *grid;

  /* The target mesh for the mapping. */
  mesh3_s const *mesh;

  /* The cell index for each grid point. */
  size_t *lc;

  /* The cell vertices for each grid point. */
  size_t (*cv)[4];

  /* The barycentric coordinates of each grid point with respect to
   * the containing cell. */
  dbl4 *b;
} grid2_to_mesh3_mapping_s;

void grid2_to_mesh3_mapping_init_xy(grid2_to_mesh3_mapping_s *mapping,
                                    grid2_s const *grid, mesh3_s const *mesh,
                                    dbl z);
void grid2_to_mesh3_mapping_deinit(grid2_to_mesh3_mapping_s *mapping);
