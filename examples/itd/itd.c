#include <jmm/eik3.h>
#include <jmm/mesh3.h>

int main(void) {
  dbl eps = 1e-5;
  dbl maxvol = 8e5;
  bool verbose = true;
  dbl3 xsrc = {-0.9e3, -0.9e3, -0.9e3};

  /* Units are in mm */
  char const *off_path = "HUTUB_pp2_in_cube_50k.off";

  mesh3_data_s data;
  mesh3_data_from_off_file(&data, off_path, maxvol, verbose);
  mesh3_data_insert_vert(&data, xsrc, eps);

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, &eps);

  mesh3_dump_verts(mesh, "verts.bin");
  mesh3_dump_cells(mesh, "cells.bin");

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);
}
