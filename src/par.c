#include <jmm/par.h>

#include <assert.h>
#include <math.h>

#include <jmm/array.h>
#include <jmm/eik3.h>
#include <jmm/mesh3.h>

#include "macros.h"

void par2_init_empty(par2_s *par) {
  par->l[0] = par->l[1] = NO_PARENT;
  par->b[0] = par->b[1] = NAN;
}

bool par2_is_empty(par2_s const *par) {
  return par->l[0] == NO_PARENT;
}

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
  assert(size <= 3);
  return size;
}

void par3_get_xb(par3_s const *par, mesh3_s const *mesh, dbl xb[3]) {
  size_t l[3];
  dbl b[3];
  size_t num_active = par3_get_active(par, l, b);

  dbl const *x[3];
  mesh3_get_vert_ptrs(mesh, l, num_active, x);

  dbl3_zero(xb);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < num_active; ++j)
      xb[i] += b[j]*x[j][i];
}

bool par3_is_empty(par3_s const *par) {
  return par->l[0] == NO_PARENT;
}

size_t par3_num_active(par3_s const *par) {
  dbl const atol = 1e-14;
  size_t num_active = 0;
  for (size_t i = 0; i < 3; ++i)
    num_active += par->l[i] != NO_PARENT && par->b[i] > atol;
  return num_active;
}

size_t par3_get_active_inds(par3_s const *par, size_t l[3]) {
  dbl const atol = 1e-14;
  size_t num_active = 0;
  for (size_t i = 0; i < 3; ++i)
    if (par->l[i] != NO_PARENT && par->b[i] > atol)
      l[num_active++] = par->l[i];
  for (size_t i = num_active; i < 3; ++i)
    l[i] = (size_t)NO_PARENT;
  return num_active;
}

size_t par3_get_active_and_inactive_inds(par3_s const *par, uint3 la, uint3 li) {
  dbl const atol = 1e-14;

  bool active[3] = {false, false, false};

  size_t num_active = 0;
  for (size_t i = 0; i < 3; ++i) {
    if (par->l[i] != NO_PARENT && par->b[i] > atol) {
      active[i] = true;
      la[num_active++] = par->l[i];
    }
  }
  for (size_t i = num_active; i < 3; ++i)
    la[i] = (size_t)NO_PARENT;

  size_t num_inactive = 0;
  for (size_t i = 0; i < 3; ++i)
    if (!active[i])
      li[num_inactive++] = par->l[i];
  for (size_t i = num_inactive; i < 3; ++i)
    li[i] = (size_t)NO_PARENT;

  return num_active;
}

size_t par3_get_active(par3_s const *par, size_t *l, dbl *b) {
  dbl const atol = 1e-14;
  size_t num_active = 0;
  for (size_t i = 0; i < 3; ++i)
    if (par->l[i] != NO_PARENT && par->b[i] > atol) {
      l[num_active] = par->l[i];
      b[num_active++] = par->b[i];
    }

  /* Normalize the active coefficients */
  dbl b_norm1 = 0;
  for (size_t i = 0; i < num_active; ++i)
    b_norm1 += fabs(b[i]);
  for (size_t i = 0; i < num_active; ++i)
    b[i] /= b_norm1;

  return num_active;
}

bool par3_has_active_parent(par3_s const *par, size_t l) {
  if (l == (size_t)NO_PARENT)
    return false;
  size_t i;
  for (i = 0; i < 3; ++i)
    if (par->l[i] == l && fabs(par->b[i]) > EPS)
      break;
  return i < 3;
}
