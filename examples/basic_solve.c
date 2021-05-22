#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <eik3.h>
#include <mesh3.h>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("usage: %s <verts.bin> <cells.bin>\n",
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
  nverts = 7613; // OVERRIDE!!!
  verts = malloc(vertsize*nverts);
  fread(verts, vertsize, nverts, fp);
  fclose(fp);

  // Read cells from binary file at argv[2]
  fp = fopen(argv[2], "rb");
  fseek(fp, 0, SEEK_END);
  ncells = ftell(fp)/cellsize;
  rewind(fp);
  ncells = 31419; // OVERRIDE!!!
  cells = malloc(cellsize*ncells);
  fread(cells, cellsize, ncells, fp);
  fclose(fp);

  // Create tetrahedron mesh for solver
  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, verts, nverts, cells, ncells, true);

  // Set up solver
  eik3_s *eik;
  eik3_alloc(&eik);
  eik3_init(eik, mesh, FTYPE_POINT_SOURCE);
  eik3_add_pt_src_BCs(eik, 0, jet3_make_point_source(0));
  eik3_solve(eik);

  /* Clean everything up */

  eik3_deinit(eik);
  eik3_dealloc(&eik);

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  free(verts);
  free(cells);
}
