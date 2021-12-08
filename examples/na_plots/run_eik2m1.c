#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "eik2m1.h"
#include "log.h"
#include "vec.h"

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
  if (argc != 8) {
    fprintf(stderr,
            "usage: %s <verts_path> <faces_path> <rfac> <jet_gt_path> "
            "<jet_path> <l_path> <lam_path>\n",
            argv[0]);
    return EXIT_FAILURE;
  }

  mesh22_s *mesh;
  mesh22_alloc(&mesh);
  mesh22_init_from_binary_files(mesh, argv[1], argv[2]);

  dbl rfac = atof(argv[3]);

  eik2m1_s *eik;
  eik2m1_alloc(&eik);
  eik2m1_init(eik, mesh);

  dbl2 xy;

  size_t nverts = mesh22_nverts(mesh);
  log_info("%lu vertices in mesh", nverts);

  jet22t *jet_gt = malloc(nverts*sizeof(jet22t));
  for (size_t l = 0; l < nverts; ++l) {
    mesh22_get_vert(mesh, l, xy);
    jet_gt[l] = get_jet_gt(xy[0], xy[1]);
  }

  /* Next, initialize a disk of radius rfac surrounding the origin
     with VALID nodes */
  size_t nfac = 0;
  for (size_t l = 0; l < nverts; ++l) {
    if (!eik2m1_is_valid(eik, l)) {
      mesh22_get_vert(mesh, l, xy);
      if (dbl2_norm(xy) < rfac) {
        ++nfac;
        eik2m1_add_valid(eik, l, jet_gt[l]);
      }
    }
  }
  log_info("factored %lu nodes", nfac);

  /* Next, initialize any nodes neighboring the VALID nodes with TRIAL
     nodes */

  size_t ntrial = 0;
  for (size_t l = 0; l < nverts; ++l) {
    if (!eik2m1_is_valid(eik, l))
      continue;

    size_t nvv = mesh22_nvv(mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh22_vv(mesh, l, vv);

    for (size_t i = 0; i < nvv; ++i) {
      if (eik2m1_is_far(eik, vv[i])) {
        ++ntrial;
        eik2m1_add_trial(eik, vv[i], jet_gt[vv[i]]);
      }
    }

    free(vv);
  }
  log_info("added %lu TRIAL nodes", ntrial);

  log_set_level(LOG_DEBUG);
  eik2m1_solve(eik);
  log_set_level(LOG_FATAL);

  FILE *fp = fopen(argv[4], "w");
  fwrite(jet_gt, sizeof(jet22t), nverts, fp);
  fclose(fp);

  fp = fopen(argv[5], "w");
  fwrite(eik2m1_get_jet_ptr(eik), sizeof(jet22t), nverts, fp);
  fclose(fp);

  fp = fopen(argv[6], "w");
  for (size_t l = 0; l < nverts; ++l) {
    par2_s par = eik2m1_get_par(eik, l);
    fwrite(par.l, sizeof(uint2), 1, fp);
  }
  fclose(fp);

  fp = fopen(argv[7], "w");
  for (size_t l = 0; l < nverts; ++l) {
    par2_s par = eik2m1_get_par(eik, l);
    fwrite(&par.b[1], sizeof(dbl), 1, fp);
  }
  fclose(fp);

  eik2m1_deinit(eik);
  eik2m1_dealloc(&eik);

  free(jet_gt);

  mesh22_deinit(mesh);
  mesh22_dealloc(&mesh);

  return EXIT_SUCCESS;
}
