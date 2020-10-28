#include "mesh3.h"

#include <stdio.h>
#include <stdlib.h>

#include "vec.h"

struct mesh3 {
  dvec3 *verts;
  size_t nverts;
  ivec3 *cells;
  size_t ncells;
};

void mesh3_alloc(mesh3_s **mesh) {
  printf("mesh3_alloc(%p)\n", mesh);

  *mesh = malloc(sizeof(mesh3_s));

  printf("  %p\n", *mesh);
}

void mesh3_dealloc(mesh3_s **mesh) {
  printf("mesh3_dealloc(%p)\n", mesh);

  // free(*mesh);
  // *mesh = NULL;
}

void mesh3_init(mesh3_s * mesh,
                dbl const *verts, size_t nverts,
                int const *cells, size_t ncells) {
  printf("mesh3_init(%p)\n", mesh);

  mesh->verts = malloc(sizeof(dvec3)*nverts);
  for (size_t i = 0; i < nverts; ++i) {
    for (int j = 0; j < 3; ++j) {
      mesh->verts[i].data[j] = verts[3*i + j];
    }
  }
  mesh->nverts = nverts;

  printf("  %f, %f, %f\n", mesh->verts[0].data[0], mesh->verts[0].data[1], mesh->verts[0].data[2]);

  mesh->cells = malloc(sizeof(ivec3)*nverts);
  for (size_t i = 0; i < ncells; ++i) {
    for (int j = 0; j < 3; ++j) {
      mesh->cells[i].data[j] = cells[3*i + j];
    }
  }
  mesh->ncells = ncells;

  printf("  %p\n  %p\n", mesh->verts, mesh->cells);

  printf("  returning from mesh3_init\n");
}

void mesh3_deinit(mesh3_s *mesh) {
  printf("mesh3_deinit(%p)\n  %p\n  %p\n", mesh, mesh->verts, mesh->cells);

  // free(mesh->verts);
  // free(mesh->cells);
  // mesh->verts = NULL;
  // mesh->cells = NULL;
}
