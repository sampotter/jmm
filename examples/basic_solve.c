#include <stdio.h>
#include <stdlib.h>

#include <eik3.h>
#include <mesh3.h>
#include <vec.h>

#define R0 0.1

int main(int argc, char *argv[]) {
  if (argc != 5) {
    printf("usage: %s <verts.bin> <cells.bin> <indsrc> <jets.bin>\n",
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
  mesh3_init(mesh, verts, nverts, cells, ncells);

  /**
   * Solve point source problem
   */

  // Get index of point source
  size_t indsrc = atoi(argv[3]);

  // Set up and run solver
  eik3_s *eik;
  eik3_alloc(&eik);
  eik3_init(eik, mesh);
  eik3_add_trial(eik, indsrc, jet3_make_point_source(0));
  eik3_solve(eik);

  // Write jets to disk
  fp = fopen(argv[4], "wb");
  fwrite(eik3_get_jet_ptr(eik), sizeof(jet3), nverts, fp);
  fclose(fp);

  /**
   * Teardown
   */

  eik3_deinit(eik);
  eik3_dealloc(&eik);

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  free(verts);
  free(cells);
}
