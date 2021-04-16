#include "par.h"

#include <assert.h>
#include <math.h>

#include "eik3.h"
#include "macros.h"
#include "mesh3.h"

par3_s make_par3(size_t l[3], dbl b[3]) {
  par3_s par = {.l = {l[0], l[1], l[2]}, .b = {b[0], b[1], b[2]}};
  for (int i = 0; i < 2; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      if (par.b[i] < par.b[j]) {
        SWAP(par.l[i], par.l[j]);
        SWAP(par.b[i], par.b[j]);
      }
    }
  }
  dbl const atol = 1e-15;
  assert(par.b[0] > atol);
  if (par.b[1] <= atol) {
    par.l[1] = NO_PARENT;
    par.b[1] = NAN;
  }
  if (par.b[2] <= atol) {
    par.l[2] = NO_PARENT;
    par.b[2] = NAN;
  }
  return par;
}

void par3_init_empty(par3_s *par) {
  par->l[0] = par->l[1] = par->l[2] = NO_PARENT;
  par->b[0] = par->b[1] = par->b[2] = NAN;
}

void par3_set(par3_s *par, size_t const *l, dbl const *b, int n) {
  for (int i = 0; i < n; ++i) {
    par->l[i] = l[i];
    par->b[i] = b[i];
  }
}

size_t par3_size(par3_s const *par) {
  size_t size = (int)(par->l[0] != NO_PARENT) + (int)(par->l[1] != NO_PARENT)
    + (int)(par->l[2] != NO_PARENT);
  assert(0 <= size && size <= 3);
  return size;
}

void par3_get_xb(par3_s const *par, eik3_s const *eik, dbl xb[3]) {
  dbl const *x[3];
  mesh3_get_vert_ptrs(eik3_get_mesh(eik), par->l, 3, x);

  for (int i = 0; i < 3; ++i) {
    xb[i] = 0;
    for (int j = 0; j < 3; ++j) {
      xb[i] += par->b[j]*x[j][i];
    }
  }
}
