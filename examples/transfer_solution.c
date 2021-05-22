#include <stdio.h>
#include <stdlib.h>

#include "def.h"
#include "eik3.h"
#include "grid3.h"
#include "mesh3.h"
#include "xfer.h"

int main(int argc, char *argv[]) {
  if (argc != 5) {
    printf("usage: %s <verts.bin> <cells.bin> <indsrc> <T.bin>\n",
           argv[0]);
    exit(EXIT_FAILURE);
  }

  size_t vertsize = 3*sizeof(dbl), nverts;
  size_t cellsize = 4*sizeof(size_t), ncells;
  dbl *verts;
  size_t *cells;

  FILE *fp = NULL;

  /**
   * Set up mesh
   */

  // Read verts from binary file at argv[1]
  fp = fopen(argv[1], "rb");
  fseek(fp, 0, SEEK_END);
  nverts = ftell(fp)/vertsize;
  rewind(fp);
  verts = malloc(vertsize*nverts);
  fread(verts, vertsize, nverts, fp);
  fclose(fp);

  // Read cells from binary file at argv[2]
  fp = fopen(argv[2], "rb");
  fseek(fp, 0, SEEK_END);
  ncells = ftell(fp)/cellsize;
  rewind(fp);
  cells = malloc(cellsize*ncells);
  fread(cells, cellsize, ncells, fp);
  fclose(fp);

  // Create tetrahedron mesh for solver
  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, verts, nverts, cells, ncells, true);

  /**
   * Solve point source problem
   */

  // Get index of point source
  size_t indsrc = atoi(argv[3]);

  // Set up and run solver
  eik3_s *eik;
  eik3_alloc(&eik);
  eik3_init(eik, mesh, FTYPE_POINT_SOURCE);
  eik3_add_pt_src_BCs(eik, indsrc, jet3_make_point_source(0));
  eik3_solve(eik);

  // Set up grid for transfer
  rect3 bbox;
  mesh3_get_bbox(mesh, &bbox);

  grid3_s grid;
  grid.h = 0.25;
  for (int i = 0; i < 3; ++i) {
    grid.min[i] = bbox.min[i];
    grid.dim[i] = ceil(bbox.max[i] - bbox.min[i])/grid.h;
  }

  dbl *T = malloc(grid3_size(&grid)*sizeof(dbl));
  xfer(mesh, eik3_get_jet_ptr(eik), &grid, T);

  fp = fopen(argv[4], "wb");
  fwrite(T, sizeof(dbl), grid3_size(&grid), fp);
  fclose(fp);

  free(T);
  free(cells);
  free(verts);
}
