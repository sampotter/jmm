#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <eik3.h>
#include <mesh3.h>

#define R0 0.1

int main(int argc, char *argv[]) {
  if (argc != 4) {
    printf("usage: %s <verts.bin> <cells.bin> <indsrc or bcpath>\n",
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

  // Set up solver
  eik3_s *eik;
  eik3_alloc(&eik);
  eik3_init(eik, mesh);

  // Set boundary conditions from passed index or file

  errno = 0;
  char *endptr;
  size_t indsrc = strtol(argv[3], &endptr, 10);
  if (endptr != argv[3]) {
    eik3_add_trial(eik, indsrc, jet3_make_point_source(0));
  } else {
    fp = fopen(argv[3], "r");
    jet3 jet;
    while (fscanf(fp, "%lu %lf %lf %lf %lf\n",
                  &indsrc, &jet.f, &jet.fx, &jet.fy, &jet.fz)
           != EOF)
      eik3_add_trial(eik, indsrc, jet);
  }

  eik3_solve(eik);

  /* Clean everything up */

  eik3_deinit(eik);
  eik3_dealloc(&eik);

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  free(verts);
  free(cells);
}
