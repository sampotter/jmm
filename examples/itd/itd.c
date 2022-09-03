#include <stdio.h>
#include <string.h>

#include <jmm/bmesh.h>
#include <jmm/eik3.h>
#include <jmm/mesh3.h>

int main(void) {
  dbl eps = 1e-5;
  dbl maxvol = 8e5;
  bool verbose = true;
  dbl c = 340.3e3; // mm/s

  size_t num_el = 256;
  size_t num_az = 2*num_el;

  /* NOTE: units in mm! */

  dbl rfac = 4; // mm
  dbl3 xsrc_R = {0.0, -69.910568, 0.0}; // right ear
  dbl3 xsrc_L = {0.0, 69.821861, 0.0}; // left ear

  dbl fliege_grid_r = 500; // mm

  char const *itd_bin_path = "itd.bin";

  /* Load Fliege grid points */
  char const *fliege_grid_txt_path = "fliege64.txt";
  char const *fliege_grid_bin_path = "fliege64.bin";

  array_s *fliege_grid_arr;
  array_alloc(&fliege_grid_arr);
  array_init(fliege_grid_arr, sizeof(dbl3), ARRAY_DEFAULT_CAPACITY);

  /* Read CSV text file containing Fliege grid node elevation (first
   * column) and azimuth (second column) */
  FILE *fp = fopen(fliege_grid_txt_path, "r");
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  while ((read = getline(&line, &len, fp)) != -1) {
    dbl el, az;
    char *saveptr, *pc;
    pc = strtok_r(line, " ", &saveptr);
    el = atof(pc);
    pc = strtok_r(NULL, " ", &saveptr);
    az = atof(pc);
    dbl3 point = {cos(az)*sin(el), sin(az)*sin(el), cos(el)};
    dbl3_dbl_mul_inplace(point, fliege_grid_r);
    array_append(fliege_grid_arr, &point);
  }
  fclose(fp);

  /* Dump Cartesian coordinates of Fliege grid to a binary file */
  fp = fopen(fliege_grid_bin_path, "wb");
  for (size_t i = 0; i < array_size(fliege_grid_arr); ++i) {
    dbl3 point;
    array_get(fliege_grid_arr, i, &point);
    fwrite(point, sizeof(dbl3), 1, fp);
  }
  fclose(fp);

  /* Units are in mm */
  char const *off_path = "HUTUB_pp2_in_cube_50k_origin_fixed.off";

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

  /* Solve the eikonal equation for the left ear */
  eik3_s *eik_L;
  eik3_alloc(&eik_L);
  eik3_init(eik_L, mesh);
  eik3_add_pt_src_bcs(eik_L, xsrc_L, rfac);
  eik3_solve(eik_L);
  eik3_dump_jet(eik_L, "jet_L.bin");

  /* Solve the eikonal equation for the right ear */
  eik3_s *eik_R;
  eik3_alloc(&eik_R);
  eik3_init(eik_R, mesh);
  eik3_add_pt_src_bcs(eik_R, xsrc_R, rfac);
  eik3_solve(eik_R);
  eik3_dump_jet(eik_R, "jet_R.bin");

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

  /* Compute ITD at Fliege grid nodes and save to disk */
  fp = fopen(itd_bin_path, "wb");
  for (size_t i = 0; i < array_size(fliege_grid_arr); ++i) {
    dbl3 point;
    array_get(fliege_grid_arr, i, &point);
    dbl T_L = bmesh33_f(tau_L, point)/c;
    dbl T_R = bmesh33_f(tau_R, point)/c;
    dbl itd = T_R - T_L;
    fwrite(&itd, sizeof(dbl), 1, fp);
  }
  fclose(fp);

  fp = fopen("itd_uniform.bin", "wb");
  for (size_t i = 0; i < num_el; ++i) {
    dbl el = JMM_PI*i/(dbl)num_el;
    for (size_t j = 0; j < num_az; ++j) {
      dbl az = 2*JMM_PI*j/(dbl)num_az;
      dbl3 point = {cos(az)*sin(el), sin(az)*sin(el), cos(el)};
      dbl3_dbl_mul_inplace(point, fliege_grid_r);
      dbl T_L = bmesh33_f(tau_L, point)/c;
      dbl T_R = bmesh33_f(tau_R, point)/c;
      dbl itd = T_R - T_L;
      fwrite(&itd, sizeof(dbl), 1, fp);
    }
  }
  fclose(fp);

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

  array_deinit(fliege_grid_arr);
  array_dealloc(&fliege_grid_arr);
}
