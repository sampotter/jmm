#include <jmm/grid2.h>

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <jmm/mesh3.h>
#include <jmm/vec.h>

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

/* Initialize a mapping from `grid` to `mesh` by treating `grid` as a
 * horizontal slice through the space containing `mesh` with the given
 * z value. */
void grid2_to_mesh3_mapping_init_xy(grid2_to_mesh3_mapping_s *mapping,
                                    grid2_s const *grid, mesh3_s const *mesh,
                                    dbl z) {
  mapping->grid = grid;
  mapping->mesh = mesh;

  size_t nind = grid2_nind(grid);

  /* Allocate space for cell indices, incident vertices, and
   * barycentric coordinates */
  mapping->lc = malloc(nind*sizeof(size_t));
  mapping->cv = malloc(nind*sizeof(size_t[4]));
  mapping->b = malloc(nind*sizeof(dbl4));

  /* Initialize everything to a bad value */
  for (size_t i = 0; i < nind; ++i) {
    mapping->lc[i] = (size_t)NO_INDEX;
    for (size_t j = 0; j < 4; ++j) {
      mapping->cv[i][j] = (size_t)NO_INDEX;
      mapping->b[i][j] = NAN;
    }
  }

  /* Iterate over each tetrahedron in the mesh, pull out the subgrid
   * of `grid` which is incident on the tetrahedron, and then find the
   * barycentric coordinates of each grid point which actually lies
   * inside the tetrahedral cell. */
  for (size_t lc = 0; lc < mesh3_ncells(mesh); ++lc) {
    tetra3 tetra = mesh3_get_tetra(mesh, lc);

    size_t offset[2];
    grid2_s subgrid = tetra3_get_covering_xy_subgrid(&tetra, grid, offset);

    for (size_t l = 0; l < grid2_nind(&subgrid); ++l) {
      int ind[2];
      grid2_l2ind(&subgrid, l, ind);

      int ind_orig[2] = {offset[0] + ind[0], offset[1] + ind[1]};
      size_t l_orig = grid2_ind2l(grid, ind_orig);
      if (mapping->lc[l_orig] != (size_t)NO_INDEX)
        continue;

      dbl3 x = {[2] = z};
      grid2_l2xy(grid, l_orig, x);
      if (!mesh3_cell_contains_point(mesh, lc, x))
        continue;

      mapping->lc[l_orig] = lc;
      mesh3_cv(mesh, lc, mapping->cv[l_orig]);
      tetra3_get_bary_coords(&tetra, x, mapping->b[l_orig]);
    }
  }
}

void grid2_to_mesh3_mapping_deinit(grid2_to_mesh3_mapping_s *mapping) {
  mapping->grid = NULL;
  mapping->mesh = NULL;

  free(mapping->lc);
  mapping->lc = NULL;

  free(mapping->cv);
  mapping->cv = NULL;

  free(mapping->b);
  mapping->b = NULL;
}
