#include <jmm/edge.h>

#include <assert.h>

#include <jmm/eik3.h>
#include <jmm/mesh3.h>

#include "macros.h"

edge_s make_edge(size_t l0, size_t l1) {
  return (edge_s) {.l = {MIN(l0, l1), MAX(l0, l1)}};
}

void edge_get_x0_and_dx(edge_s edge, mesh3_s const *mesh, dbl x0[3], dbl dx[3]) {
  mesh3_copy_vert(mesh, edge.l[0], x0);
  dbl3_sub(mesh3_get_vert_ptr(mesh, edge.l[1]), x0, dx);
}

void edge_get_xt(edge_s edge, mesh3_s const *mesh, dbl t, dbl xt[3]) {
  dbl const *x0 = mesh3_get_vert_ptr(mesh, edge.l[0]);
  dbl dx[3];
  dbl3_sub(mesh3_get_vert_ptr(mesh, edge.l[1]), x0, dx);
  dbl3_saxpy(t, dx, x0, xt);
}
