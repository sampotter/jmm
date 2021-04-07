#include "edge.h"

#include <assert.h>

#include "eik3.h"
#include "macros.h"
#include "mesh3.h"

edge_s make_edge(size_t l0, size_t l1) {
  return (edge_s) {.l = {MIN(l0, l1), MAX(l0, l1)}};
}

size_t edge_get_valid_index(edge_s edge, eik3_s const *eik) {
  size_t const *l = edge.l;
  size_t lv = eik3_is_valid(eik, l[0]) ? l[0] : l[1];
  assert(eik3_is_valid(eik, lv));
  return lv;
}

size_t edge_get_shadow_index(edge_s edge, eik3_s const *eik) {
  size_t const *l = edge.l;
  size_t ls = eik3_is_shadow(eik, l[0]) ? l[0] : l[1];
  assert(eik3_is_shadow(eik, ls));
  return ls;
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
