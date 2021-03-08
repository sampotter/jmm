#include "cubic.h"

#include <assert.h>
#include <string.h>

#include <gsl/gsl_poly.h>

#include "mat.h"

static dmat44 V_inv = {
  .data = {
    { 1,  0,  0,  0},
    { 0,  0,  1,  0},
    {-3,  3, -2, -1},
    { 2, -2,  1,  1}
  }
};

/**
 * This function is a bit of a stop gap until we move the dvec* and
 * dmat* types.
 */
cubic_s cubic_from_data(dbl f[2], dbl Df[2]) {
  cubic_s cubic;
  cubic_set_data(&cubic, (dvec4) {.data = {f[0], f[1], Df[0], Df[1]}});
  return cubic;
}

void cubic_set_data(cubic_s *cubic, dvec4 data) {
  cubic->a = dmat44_dvec4_mul(V_inv, data);
}

void cubic_set_data_from_ptr(cubic_s *cubic, dbl const *data_ptr) {
  dvec4 data;
  memcpy((void *)data.data, (void *)data_ptr, 4*sizeof(dbl));
  cubic_set_data(cubic, data);
}

void cubic_reverse_on_unit_interval(cubic_s *cubic) {
  static dmat44 M = {
    .data = {
      {1,  1,  1,  1},
      {0, -1, -2, -3},
      {0,  0,  1,  3},
      {0,  0,  0, -1}
    }
  };
  cubic->a = dmat44_dvec4_mul(M, cubic->a);
}

dbl cubic_f(cubic_s const *cubic, dbl lam) {
  return dvec4_dot(cubic->a, dvec4_m(lam));
}

dbl cubic_df(cubic_s const *cubic, dbl lam) {
  return dvec4_dot(cubic->a, dvec4_dm(lam));
}

dbl cubic_d2f(cubic_s const *cubic, dbl lam) {
  return dvec4_dot(cubic->a, dvec4_d2m(lam));
}

void cubic_make_monic(cubic_s *cubic) {
  dbl a3 = cubic->a.data[3];
  for (int i = 0; i < 4; ++i)
    cubic->a.data[i] /= a3;
}

/**
 * Compute the real roots of `cubic`, storing them in `lam`, and
 * returning the number of real roots.
 */
int cubic_get_real_roots(cubic_s const *cubic, dbl lam[3]) {
  cubic_s monic = {.a = cubic->a};

  cubic_make_monic(&monic);

  dbl const *a = &monic.a.data[0];

  /* Note that we reverse the order of `a` here when transferring the
   * coefficients to GSL's `gsl_poly_solve_cubic`. */
  return gsl_poly_solve_cubic(a[2], a[1], a[0], &lam[0], &lam[1], &lam[2]);
}

void cubic_add_constant(cubic_s *cubic, dbl f) {
  cubic->a.data[0] += f;
}
