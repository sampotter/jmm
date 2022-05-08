#include "grid2.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "vec.h"

void grid2_save(grid2_s const *grid, char const *path) {
  FILE *fp = fopen(path, "w");
  fprintf(fp, "shape: %d, %d\n", grid->shape[0], grid->shape[1]);
  fprintf(fp, "xymin: %g, %g\n", grid->xymin[0], grid->xymin[1]);
  fprintf(fp, "h: %g\n", grid->h);
  fprintf(fp, "order: %s",
          grid->order == ORDER_ROW_MAJOR ? "row major" : "column major");
  fclose(fp);
}

size_t grid2_nind(grid2_s const *grid) {
  return grid->shape[0]*grid->shape[1];
}

size_t grid2_nindc(grid2_s const *grid) {
  return (grid->shape[0] - 1)*(grid->shape[1] - 1);
}

int grid2_ind2l(grid2_s const *grid, int2 const ind) {
  return grid->order == ORDER_ROW_MAJOR ?
    ind[1] + grid->shape[1]*ind[0] :
    grid->shape[0]*ind[1] + ind[0];
}

int grid2_ind2lc(grid2_s const *grid, int2 const ind) {
  return grid->order == ORDER_ROW_MAJOR ?
    ind[1] + (grid->shape[1] - 1)*ind[0] :
    (grid->shape[0] - 1)*ind[1] + ind[0];
}

int grid2_indc2l(grid2_s const *grid, int2 const indc) {
  return grid->order == ORDER_ROW_MAJOR ?
    indc[1] + grid->shape[1]*indc[0] :
    grid->shape[0]*indc[1] + indc[0];
}

int grid2_indc2lc(grid2_s const *grid, int2 const indc) {
  return grid->order == ORDER_ROW_MAJOR ?
    indc[1] + (grid->shape[1] - 1)*indc[0] :
    (grid->shape[0] - 1)*indc[1] + indc[0];
}

void grid2_l2ind(grid2_s const *grid, int l, int2 ind) {
  if (grid->order == ORDER_ROW_MAJOR) {
    ind[0] = l/grid->shape[1];
    ind[1] = l % grid->shape[1];
  } else {
    ind[0] = l % grid->shape[0];
    ind[1] = l/grid->shape[0];
  }
}

void grid2_l2indc(grid2_s const *grid, int l, int2 indc) {
  if (grid->order == ORDER_ROW_MAJOR) {
    indc[0] = l/grid->shape[1];
    indc[1] = l % grid->shape[1];
  } else {
    indc[0] = l % grid->shape[0];
    indc[1] = l/grid->shape[0];
  }
}

void grid2_lc2ind(grid2_s const *grid, int lc, int2 ind) {
  if (grid->order == ORDER_ROW_MAJOR) {
    ind[0] = lc/(grid->shape[1] - 1);
    ind[1] = lc % (grid->shape[1] - 1);
  } else {
    ind[0] = lc % (grid->shape[0] - 1);
    ind[1] = lc/(grid->shape[0] - 1);
  }
}

void grid2_lc2indc(grid2_s const *grid, int lc, int2 indc) {
  if (grid->order == ORDER_ROW_MAJOR) {
    indc[0] = lc/(grid->shape[1] - 1);
    indc[1] = lc % (grid->shape[1] - 1);
  } else {
    indc[0] = lc % (grid->shape[0] - 1);
    indc[1] = lc/(grid->shape[0] - 1);
  }
}

int grid2_l2lc(grid2_s const *grid, int l) {
  if (grid->order == ORDER_ROW_MAJOR) {
    return l - l/grid->shape[1];
  } else {
    return l - l/grid->shape[0];
  }
}

int grid2_lc2l(grid2_s const *grid, int lc) {
  if (grid->order == ORDER_ROW_MAJOR) {
    return lc + lc/(grid->shape[1] - 1);
  } else {
    return lc + lc/(grid->shape[0] - 1);
  }
}

void grid2_l2xy(grid2_s const *grid, int l, dbl2 xy) {
  int2 ind;
  grid2_l2ind(grid, l, ind);

  xy[0] = grid->h*ind[0] + grid->xymin[0];
  xy[1] = grid->h*ind[1] + grid->xymin[1];
}

int grid2_xy2lc(grid2_s const *grid, dbl2 const xy, dbl2 cc) {
#if SJS_DEBUG
  assert(cc != NULL);
#endif

  dbl2_sub(xy, grid->xymin, cc);
  dbl2_dbl_div_inplace(cc, grid->h);
  dbl2 ind_; dbl2_floor(cc, ind_);
  dbl2_sub_inplace(cc, ind_);
  int2 ind = {ind_[0], ind_[1]};

  if (ind[0] < 0) {
    ind[0] = 0;
    cc[0] = 0.0;
  }

  if (ind[1] < 0) {
    ind[1] = 0;
    cc[1] = 0.0;
  }

  if (ind[0] >= grid->shape[0] - 1) {
    --ind[0];
    cc[0] = 1.0;
  }

  if (ind[1] >= grid->shape[1] - 1) {
    --ind[1];
    cc[1] = 1.0;
  }

  return grid2_ind2lc(grid, ind);
}

bool grid2_isind(grid2_s const *grid, int2 ind) {
  return 0 <= ind[0] && ind[0] < grid->shape[0] &&
         0 <= ind[1] && ind[1] < grid->shape[1];
}

bool grid2_isindc(grid2_s const *grid, int2 indc) {
  return 0 <= indc[0] && indc[0] < grid->shape[0] - 1 &&
         0 <= indc[1] && indc[1] < grid->shape[1] - 1;
}

void grid2_get_inbounds(grid2_s const *grid, grid2info_s const *info,
                        int l, bool inbounds[GRID2_NUM_NB + 1]) {
  int2 ind, ind0;
  grid2_l2ind(grid, l, ind);
  for (int i0 = 0; i0 < GRID2_NUM_NB + 1; ++i0) {
    int2_add(ind, info->offsets[i0], ind0);
    inbounds[i0] = grid2_isind(grid, ind0);
  }
}

void grid2_get_nb(grid2_s const *grid, grid2info_s const *info,
                  int l, int l_nb[GRID2_NUM_NB],
                  bool inbounds[GRID2_NUM_NB]) {
  int2 ind, ind_nb;
  grid2_l2ind(grid, l, ind);
  for (int i_nb = 0; i_nb < GRID2_NUM_NB; ++i_nb) {
    int2_add(ind, info->offsets[i_nb], ind_nb);
    inbounds[i_nb] = grid2_isind(grid, ind_nb);
    l_nb[i_nb] = inbounds[i_nb] ? grid2_ind2l(grid, ind_nb) : NO_INDEX;
  }
}

void grid2info_init(grid2info_s *info, grid2_s const *grid) {
  static int2 offsets[GRID2_NUM_NB + 1] = {
    {-1, -1},
    {-1,  0},
    {-1,  1},
    { 0,  1},
    { 1,  1},
    { 1,  0},
    { 1, -1},
    { 0, -1},
    {-1, -1}
  };
  memcpy(info->offsets, offsets, 9*sizeof(int2));

  for (int i = 0; i < GRID2_NUM_NB + 1; ++i)
    info->nb_dl[i] = grid2_ind2l(grid, offsets[i]);
}
