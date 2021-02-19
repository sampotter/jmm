#include <stdio.h>
#include <stdlib.h>

#include <bmesh.h>
#include <eik3.h>
#include <mesh2.h>
#include <mesh3.h>
#include <rtree.h>

#define LEVEL 0.5

int main(int argc, char *argv[]) {
  if (argc != 4) {
    printf("usage: %s <verts.bin> <cells.bin> <indsrc>\n", argv[0]);
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
  eik3_init(eik, mesh);
  eik3_add_trial(eik, indsrc, jet3_make_point_source(0));
  eik3_solve(eik);

  // Get Bezier tetra mesh interpolation solution
  bmesh33_s *bmesh;
  bmesh33_alloc(&bmesh);
  bmesh33_init_from_mesh3_and_jets(bmesh, mesh, eik3_get_jet_ptr(eik));

  // Get surface mesh
  mesh2_s *surf_mesh = mesh3_get_surface_mesh(mesh);

  // Create R-tree and insert surface mesh
  rtree_s *rtree;
  rtree_alloc(&rtree);
  rtree_init(rtree, 32, RTREE_SPLIT_STRATEGY_SURFACE_AREA);
  rtree_insert_mesh2(rtree, surf_mesh);

  bmesh33_s *level_bmesh = bmesh33_get_level_bmesh(bmesh, LEVEL);

  rtree_s *level_rtree = rtree_copy(rtree);
  mesh3_s const *level_mesh = bmesh33_get_mesh_ptr(level_bmesh);
  rtree_insert_mesh3(level_rtree, level_mesh);
  rtree_build(level_rtree);

  rtree_deinit(level_rtree);
  rtree_dealloc(&level_rtree);

  bmesh33_deinit(level_bmesh);
  bmesh33_dealloc(&level_bmesh);

  /**
   * Teardown
   */

  bmesh33_deinit(bmesh);
  bmesh33_dealloc(&bmesh);

  rtree_deinit(rtree);
  rtree_dealloc(&rtree);

  mesh2_deinit(surf_mesh);
  mesh2_dealloc(&surf_mesh);

  eik3_deinit(eik);
  eik3_dealloc(&eik);

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  free(verts);
  free(cells);
}
