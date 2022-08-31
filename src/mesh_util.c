#include "mesh_util.h"

#include <assert.h>

#include <jmm/mat.h>
#include <jmm/util.h>
#include <jmm/vec.h>

bool face_in_cell(size_t const f[3], size_t const c[4]) {
  return point_in_cell(f[0], c) && point_in_cell(f[1], c) &&
    point_in_cell(f[2], c);
}

bool point_in_face(size_t l, size_t const f[3]) {
  return f[0] == l || f[1] == l || f[2] == l;
}

bool point_in_cell(size_t l, size_t const c[4]) {
  return c[0] == l || c[1] == l || c[2] == l || c[3] == l;
}

bool edge_in_face(size_t const le[2], size_t const lf[3]) {
  assert(le[0] != le[1]);
  assert(lf[0] != lf[1]);
  assert(lf[1] != lf[2]);
  assert(lf[2] != lf[0]);

  return (le[0] == lf[0] || le[0] == lf[1] || le[0] == lf[2])
      && (le[1] == lf[0] || le[1] == lf[1] || le[1] == lf[2]);
}

int edge_cmp(size_t const e1[2], size_t const e2[2]) {
  int cmp = compar_size_t(&e1[0], &e2[0]);
  return cmp != 0 ? cmp : compar_size_t(&e1[1], &e2[1]);
}

void R_from_n(dbl3 const n, dbl33 R) {
  dbl33 nnT;
  dbl3_outer(n, n, nnT);

  dbl33 eye;
  dbl33_eye(eye);

  dbl33_saxpy(-2, nnT, eye, R);
}
