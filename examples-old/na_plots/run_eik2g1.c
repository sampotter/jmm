#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <eik2g1.h>
#include <vec.h>

jet22t get_jet_gt(dbl x, dbl y) {
  dbl r = sqrt(x*x + y*y);
  dbl rx = x/r, ry = y/r;
  return (jet22t) {
    .f = r,
    .Df = {rx, ry},
    .D2f = {
      {(1 - rx*rx)/r, -ry*rx/r},
      {-rx*ry/r, (1 - ry*ry)/r}
    }
  };
}

int main(int argc, char *argv[]) {
  if (argc != 7) {
    fprintf(stderr, "usage: %s <n> <rfac> <jet_gt_path> <jet_path> <l_path> <lam_path>\n", argv[0]);
    return EXIT_FAILURE;
  }

  int n = atoi(argv[1]);
  dbl rfac = atof(argv[2]);

  grid2_s grid = {
    .shape = {n, n},
    .xymin = {-1, -1},
    .h = 2.0/(n - 1),
    .order = ORDER_ROW_MAJOR
  };

  grid2info_s info;
  grid2info_init(&info, &grid);

  eik2g1_s *eik;
  eik2g1_alloc(&eik);
  eik2g1_init(eik, &grid);

  dbl2 xy;
  int2 ind, ind_nb;
  int l_nb[GRID2_NUM_NB];
  bool inbounds[GRID2_NUM_NB];

  jet22t *jet_gt = malloc(grid2_nind(&grid)*sizeof(jet22t));
  for (size_t l = 0; l < grid2_nind(&grid); ++l) {
    grid2_l2xy(&grid, l, xy);
    jet_gt[l] = get_jet_gt(xy[0], xy[1]);
  }

  /* First, initialize the point source and its 8 neighbors with VALID
     data to make sure we don't try to update from anything singular  */
  ind[0] = n/2;
  ind[1] = n/2;
  grid2_get_nb(&grid, &info, grid2_ind2l(&grid, ind), l_nb, inbounds);
  for (size_t i = 0; i < GRID2_NUM_NB; ++i) {
    if (inbounds[i]) {
      grid2_l2ind(&grid, l_nb[i], ind_nb);
      eik2g1_add_valid(eik, ind_nb, jet_gt[l_nb[i]]);
    }
  }

  /* Next, initialize a disk of radius rfac surrounding the origin
     with VALID nodes */
  for (size_t l = 0; l < grid2_nind(&grid); ++l) {
    grid2_l2ind(&grid, l, ind);
    if (!eik2g1_is_valid(eik, ind)) {
      grid2_l2xy(&grid, l, xy);
      if (dbl2_norm(xy) < rfac) {
        eik2g1_add_valid(eik, ind, jet_gt[l]);
      }
    }
  }

  /* Next, initialize any nodes neighboring the VALID nodes with TRIAL
     nodes */
  for (size_t l = 0; l < grid2_nind(&grid); ++l) {
    grid2_l2ind(&grid, l, ind);
    if (eik2g1_is_valid(eik, ind)) {
      grid2_get_nb(&grid, &info, l, l_nb, inbounds);
      for (size_t i = 0; i < GRID2_NUM_NB; ++i) {
        grid2_l2ind(&grid, l_nb[i], ind_nb);
        if (inbounds[i] && eik2g1_is_far(eik, ind_nb)) {
          eik2g1_add_trial(eik, ind_nb, jet_gt[l_nb[i]]);
        }
      }
    }
  }

  eik2g1_solve(eik);

  FILE *fp = fopen(argv[3], "w");
  fwrite(jet_gt, sizeof(jet22t), grid2_nind(&grid), fp);
  fclose(fp);

  fp = fopen(argv[4], "w");
  fwrite(eik2g1_get_jet_ptr(eik), sizeof(jet22t), grid2_nind(&grid), fp);
  fclose(fp);

  fp = fopen(argv[5], "w");
  for (size_t l = 0; l < grid2_nind(&grid); ++l) {
    grid2_l2ind(&grid, l, ind);
    par2_s par = eik2g1_get_par(eik, ind);
    fwrite(par.l, sizeof(size_t[2]), 1, fp);
  }
  fclose(fp);

  fp = fopen(argv[6], "w");
  for (size_t l = 0; l < grid2_nind(&grid); ++l) {
    grid2_l2ind(&grid, l, ind);
    par2_s par = eik2g1_get_par(eik, ind);
    fwrite(&par.b[1], sizeof(dbl), 1, fp);
  }
  fclose(fp);

  eik2g1_deinit(eik);
  eik2g1_dealloc(&eik);

  free(jet_gt);

  return EXIT_SUCCESS;
}
