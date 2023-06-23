#include <jmm/cubic.h>

#include <assert.h>
#include <string.h>

#include <gsl/gsl_poly.h>

#include <jmm/mat.h>

static dbl44 V_inv = {
  { 1,  0,  0,  0},
  { 0,  0,  1,  0},
  {-3,  3, -2, -1},
  { 2, -2,  1,  1}
};

cubic_s cubic_from_lagrange_data(dbl4 f) {
  static dbl44 Vinv = {
    { 1. ,   0. ,   0. ,   0. },
    {-5.5,   9. ,  -4.5,   1. },
    { 9. , -22.5,  18. ,  -4.5},
    {-4.5,  13.5, -13.5,   4.5}
  };

  cubic_s cubic;
  dbl44_dbl4_mul(Vinv, f, cubic.a);
  return cubic;
}

/**
 * This function is a bit of a stop gap until we move the dvec* and
 * dbl* types.
 */
cubic_s cubic_from_data(dbl2 f, dbl2 Df) {
  cubic_s cubic;
  cubic_set_data(&cubic, DBL4(f[0], f[1], Df[0], Df[1]));
  return cubic;
}

void cubic_set_data(cubic_s *cubic, dbl4 data) {
  dbl44_dbl4_mul(V_inv, data, cubic->a);
}

void cubic_set_data_from_ptr(cubic_s *cubic, dbl const *data_ptr) {
  cubic_set_data(cubic, DBL4_FROM_PTR(data_ptr));
}

void cubic_reverse_on_unit_interval(cubic_s *cubic) {
  static dbl44 M = {
    {1,  1,  1,  1},
    {0, -1, -2, -3},
    {0,  0,  1,  3},
    {0,  0,  0, -1}
  };
  dbl4 tmp;
  dbl44_dbl4_mul(M, cubic->a, tmp);
  dbl4_copy(tmp, cubic->a);
}

dbl cubic_f(cubic_s const *cubic, dbl lam) {
  dbl4 m; dbl4_m(lam, m);
  return dbl4_dot(cubic->a, m);
}

dbl cubic_df(cubic_s const *cubic, dbl lam) {
  dbl4 dm; dbl4_dm(lam, dm);
  return dbl4_dot(cubic->a, dm);
}

dbl cubic_d2f(cubic_s const *cubic, dbl lam) {
  dbl4 d2m; dbl4_d2m(lam, d2m);
  return dbl4_dot(cubic->a, d2m);
}

void cubic_make_monic(cubic_s *cubic) {
  dbl a3 = cubic->a[3];
  for (int i = 0; i < 4; ++i)
    cubic->a[i] /= a3;
}

int
/* gsl_poly_ */ solve_cubic (double a, double b, double c,
							 double *x0, double *x1, double *x2);

/**
 * Compute the real roots of `cubic`, storing them in `lam`, and
 * returning the number of real roots.
 */
int cubic_get_real_roots(cubic_s const *cubic, dbl lam[3]) {
  cubic_s monic = *cubic;
  cubic_make_monic(&monic);

  /* Note that we reverse the order of `a` here when transferring the
   * coefficients to GSL's `gsl_poly_solve_cubic`. */
  return solve_cubic(monic.a[2], monic.a[1], monic.a[0],
                     &lam[0], &lam[1], &lam[2]);
}

int cubic_real_roots_in_interval(cubic_s const *cubic, dbl lam[3], dbl a, dbl b) {
  dbl lam_[3];
  size_t num_real_roots = cubic_get_real_roots(cubic, lam_);
  size_t num_real_roots_in_interval = 0;
  for (size_t i = 0; i < num_real_roots; ++i)
    if (a <= lam_[i] && lam_[i] <= b)
      lam[num_real_roots_in_interval++] = lam_[i];
  return num_real_roots_in_interval;
}

void cubic_add_constant(cubic_s *cubic, dbl f) {
  cubic->a[0] += f;
}
