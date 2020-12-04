#include <stdio.h>
#include <stdlib.h>

#include <eik3.h>
#include <mesh3.h>
#include <vec.h>

#define R0 0.1

int main(int argc, char *argv[]) {
  if (argc != 4) {
    printf("usage: <%s> <verts.bin> <cells.bin> <indsrc> <jets.txt>\n",
           argv[0]);
    exit(EXIT_FAILURE);
  }

  size_t vertsize = 3*sizeof(dbl), nverts;
  size_t cellsize = 4*sizeof(size_t), ncells;
  dbl *verts;
  size_t *cells;

  FILE *fp = NULL;

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

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, verts, nverts, cells, ncells);

  size_t indsrc = atoi(argv[3]);
  dvec3 xsrc;
  mesh3_get_vert(mesh, indsrc, xsrc.data);

  eik3_s *eik;
  eik3_alloc(&eik);
  eik3_init(eik, mesh);

  jet3 jet;
  dvec3 x, x_minus_xsrc;
  for (size_t l = 0; l < nverts; ++l) {
    mesh3_get_vert(mesh, l, x.data);
    x_minus_xsrc = dvec3_sub(x, xsrc);
    jet.f = dvec3_norm(x_minus_xsrc);
    if (jet.f <= R0) {
      x_minus_xsrc = dvec3_dbl_div(x_minus_xsrc, jet.f);
      jet.fx = x_minus_xsrc.data[0];
      jet.fy = x_minus_xsrc.data[1];
      jet.fz = x_minus_xsrc.data[2];
      eik3_add_valid(eik, l, jet);
    }
  }

  int nvv;
  size_t *vv;
  for (size_t l = 0; l < nverts; ++l) {
    if (eik3_is_valid(eik, l)) {
      nvv = mesh3_nvv(mesh, l);
      vv = malloc(sizeof(size_t)*nvv);
      mesh3_vv(mesh, l, vv);
      for (int i = 0; i < nvv; ++i) {
        if (!eik3_is_valid(eik, l)) {
          mesh3_get_vert(mesh, l, x.data);
          x_minus_xsrc = dvec3_sub(x, xsrc);
          jet.f = dvec3_norm(x_minus_xsrc);
          x_minus_xsrc = dvec3_dbl_div(x_minus_xsrc, jet.f);
          jet.fx = x_minus_xsrc.data[0];
          jet.fy = x_minus_xsrc.data[1];
          jet.fz = x_minus_xsrc.data[2];
          eik3_add_trial(eik, l, jet);
        }
      }
      free(vv);
    }
  }

  eik3_solve(eik);

  fp = fopen(argv[4], "wb");
  fwrite(eik3_get_jet_ptr(eik), sizeof(jet3), nverts, fp);
  fclose(fp);

  eik3_deinit(eik);
  eik3_dealloc(&eik);

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  free(verts);
  free(cells);
}
