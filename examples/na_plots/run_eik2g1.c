#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <jmm/array.h>
#include <jmm/eik2g1.h>
#include <jmm/vec.h>

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

void dump_sol_info_rec(eik2g1_s const *eik, grid2_s const *grid,
                       int2 const ind, int depth, FILE *fp,
                       array_s *visited) {
  int l = grid2_ind2l(grid, ind);
  if (array_contains(visited, &l))
    return;

  array_append(visited, &l);

  fprintf(fp, "%d %d %d ", depth, ind[0], ind[1]);

  eik2g1_sol_info_s sol_info = eik2g1_get_sol_info(eik, ind);
  fprintf(fp, "%0.15g ", sol_info.lam_T);
  fprintf(fp, "%0.15g ", sol_info.lam_tau);
  fprintf(fp, "%0.15g ", sol_info.lam_star);
  fprintf(fp, "%0.15g ", sol_info.E0);
  fprintf(fp, "%0.15g ", sol_info.E0_term1);
  fprintf(fp, "%0.15g ", sol_info.E0_term2);
  fprintf(fp, "%0.15g ", sol_info.E0_term3);

  bool has_parent = eik2g1_has_par(eik, ind);
  fprintf(fp, "%d\n", has_parent ? 1 : 0);

  if (!has_parent)
    return;

  par2_s par = eik2g1_get_par(eik, ind);

  int2 ind0, ind1;
  grid2_l2ind(grid, par.l[0], ind0);
  grid2_l2ind(grid, par.l[1], ind1);

  dump_sol_info_rec(eik, grid, ind0, depth + 1, fp, visited);
  dump_sol_info_rec(eik, grid, ind1, depth + 1, fp, visited);
}

void dump_sol_infos(eik2g1_s const *eik, grid2_s const *grid,
                    dbl2 xy_start, char const *path) {
  int2 ind = {
    (int)floor((xy_start[0] - grid->xymin[0])/grid->h),
    (int)floor((xy_start[1] - grid->xymin[1])/grid->h)
  };

  FILE *fp = fopen(path, "w");

  array_s *visited;
  array_alloc(&visited);
  array_init(visited, sizeof(int), ARRAY_DEFAULT_CAPACITY);

  dump_sol_info_rec(eik, grid, ind, 0, fp, visited);

  array_deinit(visited);
  array_dealloc(&visited);

  fclose(fp);
}

int main(int argc, char *argv[]) {
  if (argc != 8) {
    fprintf(stderr,
            "usage: %s <n> <rfac> <jet_gt_path> <jet_path> <l_path> "
            "<lam_path> <sol_infos_path>\n", argv[0]);
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

  dump_sol_infos(eik, &grid, (dbl2) {-0.2, 0.95}, argv[7]);

  // quick test...

  {
    ind[0] = (int)floor(0.45*n);
    ind[1] = (int)floor(0.25*n);
    eik2g1_sol_info_s sol_info = eik2g1_get_sol_info(eik, ind);
    printf("h = %0.15g\n", grid.h);
    printf("h^2 = %0.15g\n", pow(grid.h, 2));
    printf("h^3 = %0.15g\n", pow(grid.h, 3));
    printf("h^3.5 = %0.15g\n", pow(grid.h, 3.5));
    printf("h^4 = %0.15g\n", pow(grid.h, 4));
    printf("sol_info:\n");
    printf("  lam_T = %0.15g\n", sol_info.lam_T);
    printf("  lam_T_check = %0.15g\n", sol_info.lam_T_check);
    printf("  lam_tau = %0.15g\n", sol_info.lam_tau);
    printf("  lam_star = %0.15g\n", sol_info.lam_star);
    printf("  |lam_T - lam_tau| = %0.15g\n", fabs(sol_info.lam_T - sol_info.lam_tau));
    printf("  |lam_tau - lam_star| = %0.15g\n", fabs(sol_info.lam_tau - sol_info.lam_star));
    printf("  |lam_star - lam_T| = %0.15g\n", fabs(sol_info.lam_star - sol_info.lam_T));
    printf("  FT(lamT) = %0.15g\n", sol_info.FT_lamT);
    printf("  Ftau(lamT) = %0.15g\n", sol_info.Ftau_lamT);
    printf("  Ftau(lamtau) = %0.15g\n", sol_info.Ftau_lamtau);
    printf("  E0 = %0.15g\n", sol_info.E0);
    printf("  That = %0.15g\n", sol_info.That);
    printf("  E0_check0 = %0.15g\n", sol_info.E0_check0);
    printf("  E0_check1 = %0.15g\n", sol_info.E0_check1);
    printf("  E0_check2 = %0.15g\n", sol_info.E0_check2);
    printf("  E0_check3 = %0.15g\n", sol_info.E0_check3);
    printf("  E0_term1 = %0.15g\n", sol_info.E0_term1);
    printf("  E0_term2 = %0.15g\n", sol_info.E0_term2);
    printf("  E0_term3 = %0.15g\n", sol_info.E0_term3);
    printf("  E0_term1_v2 = %0.15g\n", sol_info.E0_term1_v2);
    printf("  E0_term1 - E0_term1_v2 = %0.15g\n", sol_info.E0_term1 - sol_info.E0_term1_v2);
    printf("  T0_error = %0.15g\n", sol_info.T0_error);
    printf("  T1_error = %0.15g\n", sol_info.T1_error);
    printf("  DT0_error = %0.15g\n", sol_info.DT0_error);
    printf("  DT1_error = %0.15g\n", sol_info.DT1_error);
  }

  eik2g1_deinit(eik);
  eik2g1_dealloc(&eik);

  free(jet_gt);

  return EXIT_SUCCESS;
}
