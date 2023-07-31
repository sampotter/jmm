#include <omp.h>
#include <stdio.h>
#include <string.h>

#include <jmm/bmesh.h>
#include <jmm/eik3.h>
#include <jmm/log.h>
#include <jmm/mesh3.h>
#include <jmm/util.h>

static void
get_evaluation_grid(size_t num_az, size_t num_el, dbl **phi, dbl **theta) {
  *phi = malloc(num_az*sizeof(dbl));
  for (size_t i = 0; i < num_az; ++i) {
    dbl t = ((dbl)i)/((dbl)(num_az - 1));
    (*phi)[i] = 2*JMM_PI*t;
  }

  *theta = malloc(num_el*sizeof(dbl));
  for (size_t i = 0; i < num_el; ++i) {
    dbl t = ((dbl)i)/((dbl)(num_el - 1));
    (*theta)[i] = JMM_PI/2 + asin(-1 + 2*t);
  }
}

int main(int argc, char const *argv[]) {
  if (argc < 4) {
    printf("usage: %s <n> <off_path> [nthreads]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  dbl eps = 1e-5;
  dbl maxvol = 8e5;
  bool verbose = true;
  dbl c = 340.3e3; // mm/s

  size_t num_el = atoi(argv[1]);
  size_t num_az = 2*num_el;
  dbl r_grid = 500; // mm

  char const *off_path = argv[2];

  size_t nthreads = omp_get_num_threads();
  if (argc >= 4)
    nthreads = atoi(argv[3]);
  omp_set_num_threads(nthreads);

  printf("running ITD simulation w/ %lu threads\n", nthreads);

  /* NOTE: units in mm! */

  /* For meshes in `./old_meshes`: */
  // dbl rfac = 4; // mm
  // dbl3 xsrc_R = {0.0, -69.910568, 0.0}; // right ear
  // dbl3 xsrc_L = {0.0, 69.821861, 0.0}; // left ear

  /* For meshes in `./meshes`: */
  dbl rfac = 3; // mm
  dbl3 xsrc_R = {0.0, -74.5, 0.0}; // right ear
  dbl3 xsrc_L = {0.0, 63.5, 0.0}; // left ear

  /* Set up evaluation grid */
  dbl *phi_grid, *theta_grid;
  get_evaluation_grid(num_az, num_el, &phi_grid, &theta_grid);

  toc();

  mesh3_data_s data;
  mesh3_data_init_from_off_file(&data, off_path, maxvol, verbose);

  /* Add the point sources to the mesh */
  mesh3_data_insert_vert(&data, xsrc_R, eps);
  mesh3_data_insert_vert(&data, xsrc_L, eps);

  /* Set up tetrahedron mesh */
  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, &eps);

  /* Write vertices and cells to disk in row-major order */
  mesh3_dump_verts(mesh, "verts.bin");
  mesh3_dump_cells(mesh, "cells.bin");

  printf("set up tetrahedron mesh [%.2fs]\n", toc());

  /* Solve the eikonal equation for the left ear */
  eik3_s *eik_L;
  eik3_alloc(&eik_L);
  eik3_init(eik_L, mesh, &SFUNC_CONSTANT);
  eik3_add_pt_src_bcs(eik_L, xsrc_L, rfac);
  eik3_solve(eik_L);
  eik3_dump_jet(eik_L, "jet_L.bin");

  printf("computed eikonal for left ear [%.2fs]\n", toc());

  /* Solve the eikonal equation for the right ear */
  eik3_s *eik_R;
  eik3_alloc(&eik_R);
  eik3_init(eik_R, mesh, &SFUNC_CONSTANT);
  eik3_add_pt_src_bcs(eik_R, xsrc_R, rfac);
  eik3_solve(eik_R);
  eik3_dump_jet(eik_R, "jet_R.bin");

  printf("computed eikonal for right ear [%.2fs]\n", toc());

  jet31t const *jet_L = eik3_get_jet_ptr(eik_L);
  jet31t const *jet_R = eik3_get_jet_ptr(eik_R);

  /* Set up tetrahedral spline interpolating jet data for left ear */
  bmesh33_s *tau_L;
  bmesh33_alloc(&tau_L);
  bmesh33_init_from_mesh3_and_jets(tau_L, mesh, jet_L);

  /* Set up tetrahedral spline interpolating jet data for right ear */
  bmesh33_s *tau_R;
  bmesh33_alloc(&tau_R);
  bmesh33_init_from_mesh3_and_jets(tau_R, mesh, jet_R);

  printf("set up tetrahedral splines for interpolating L/R eikonals [%.2fs]\n", toc());

  FILE *fp = NULL;

  dbl *itd = malloc(num_el*num_az*sizeof(dbl));

  { size_t k = 0;
#pragma omp parallel for
    for (size_t i = 0; i < num_el; ++i) {
      dbl el = theta_grid[i];
      for (size_t j = 0; j < num_az; ++j) {
        dbl az = phi_grid[j];
        dbl3 point = {cos(az)*sin(el), sin(az)*sin(el), cos(el)};
        dbl3_dbl_mul_inplace(point, r_grid);
        dbl T_L = bmesh33_f(tau_L, point)/c;
        dbl T_R = bmesh33_f(tau_R, point)/c;
        itd[k++] = T_R - T_L;
      }
    } }

  printf("interpolated ITD to grid [%.2fs]\n", toc());

  fp = fopen("itd_grid.bin", "wb");
  fwrite(itd, sizeof(dbl), num_el*num_az, fp);
  fclose(fp);

  printf("wrote ITD to itd_grid.bin [%.2fs]\n", toc());

  free(itd);

  bmesh33_deinit(tau_R);
  bmesh33_dealloc(&tau_R);

  bmesh33_deinit(tau_L);
  bmesh33_dealloc(&tau_L);

  eik3_deinit(eik_L);
  eik3_dealloc(&eik_L);

  eik3_deinit(eik_R);
  eik3_dealloc(&eik_R);

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  mesh3_data_deinit(&data);
}
