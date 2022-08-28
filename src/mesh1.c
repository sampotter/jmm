#include "mesh1.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

struct mesh1 {
  dbl3 const *verts;
  size_t nverts;
  policy_e verts_policy;

  uint2 const *edges;
  size_t nedges;
  policy_e edges_policy;

  size_t (*ve)[2];
};

void mesh1_alloc(mesh1_s **mesh) {
  *mesh = malloc(sizeof(mesh1_s));
}

void mesh1_dealloc(mesh1_s **mesh) {
  free(*mesh);
  *mesh = NULL;
}

static void init_ve(mesh1_s *mesh) {
  /* Allocate space for `ve` and initalize each entry to
   * `NO_INDEX`. */
  mesh->ve = malloc(mesh->nverts*sizeof(size_t[2]));
  for (size_t l = 0; l < mesh->nverts; ++l)
    for (size_t i = 0; i < 2; ++i)
      mesh->ve[l][i] = (size_t)NO_INDEX;

  /* Initialize `ve` */
  for (size_t le = 0; le < mesh->nedges; ++le) {
    for (size_t i = 0, l; i < 2; ++i) {
      l = mesh->edges[le][i];
      if (mesh->ve[l][0] == (size_t)NO_INDEX) {
        mesh->ve[l][0] = le;
      } else {
        assert(mesh->ve[l][1] == (size_t)NO_INDEX);
        mesh->ve[l][1] = le;
      }
    }
  }
}

void mesh1_init(mesh1_s *mesh,
                dbl3 const *verts, size_t nverts, policy_e verts_policy,
                uint2 const *edges, size_t nedges, policy_e edges_policy) {
  /* Set up `verts`: */
  switch (verts_policy) {
  case POLICY_COPY: {
    mesh->verts = malloc(nverts*sizeof(dbl3));
    memcpy((dbl3 *)mesh->verts, verts, nverts*sizeof(dbl3));
    break;
  }
  case POLICY_XFER:
  case POLICY_VIEW: {
    mesh->verts = (dbl3 *)verts;
    break;
  }
  default:
    assert(false);
  }
  mesh->nverts = nverts;
  mesh->verts_policy = verts_policy;

  /* Set up `edges`: */
  switch (edges_policy) {
  case POLICY_COPY: {
    mesh->edges = malloc(nedges*sizeof(uint2));
    memcpy((uint2 *)mesh->edges, edges, nedges*sizeof(uint2));
    break;
  }
  case POLICY_XFER:
  case POLICY_VIEW: {
    mesh->edges = edges;
    break;
  }
  default:
    assert(false);
  }
  mesh->nedges = nedges;
  mesh->edges_policy = edges_policy;

  /* Initialize the LUT used to find the edges incident on each
   * vertex (`mesh->ve`). */
  init_ve(mesh);
}

void mesh1_deinit(mesh1_s *mesh) {
  /* Deinit `verts`: */
  switch (mesh->verts_policy) {
  case POLICY_COPY:
  case POLICY_XFER: {
    free((dbl3 *)mesh->verts);
    break;
  }
  case POLICY_VIEW:
    break;
  default:
    assert(false);
  }
  mesh->nverts = 0;
  mesh->verts = NULL;
  mesh->verts_policy = POLICY_INVALID;

  /* Deinit `edges`: */
  switch (mesh->edges_policy) {
  case POLICY_COPY:
  case POLICY_XFER: {
    free((uint2 *)mesh->edges);
    break;
  }
  case POLICY_VIEW:
    break;
  default:
    assert(false);
  }
  mesh->edges = NULL;
  mesh->nedges = 0;
  mesh->edges_policy = POLICY_INVALID;
}

size_t mesh1_nedges(mesh1_s const *mesh) {
  return mesh->nedges;
}

void mesh1_ev(mesh1_s const *mesh, size_t le, uint2 l) {
  l[0] = mesh->edges[le][0];
  l[1] = mesh->edges[le][1];
}

void mesh1_ve(mesh1_s const *mesh, size_t l, size_t le[2]) {
  le[0] = mesh->ve[l][0];
  le[1] = mesh->ve[l][1];
}
