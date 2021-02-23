#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "vec.h"

#define NUM_RANDOM_TRIALS 10

Describe(vec);

BeforeEach(vec) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(vec) {}

Ensure(vec, dblN_nsum_extreme_example_works) {
  dbl x[4] = {1, 1e100, 1, -1e100};
  dbl sum = dblN_nsum(x, 4);
  assert_that(sum, is_equal_to(2));
}

Ensure(vec, dblN_nsum_works_on_bad_bb32_df_example) {
  dbl x[3] = {
    -0.055555555577250698,
    -1.0847586719884042E-11,
     0.055555555523012758
  };
  int perm[6][3] = {
    {0, 1, 2}, {0, 2, 1},
    {1, 0, 2}, {2, 0, 1},
    {1, 2, 0}, {2, 1, 0}
  };
  dbl xperm[3], nsum[6];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j)
      xperm[j] = x[perm[i][j]];
    nsum[i] = dblN_nsum(xperm, 3);
  }
  // Assert that all six of these summations are *exactly* equal to
  // one another. This should hold for this choice of x. Not sure
  // about other choices (probably a small (~1ULP) error can be
  // introduced).
  assert_that(nsum[0], is_equal_to(nsum[1]));
  assert_that(nsum[0], is_equal_to(nsum[2]));
  assert_that(nsum[0], is_equal_to(nsum[3]));
  assert_that(nsum[0], is_equal_to(nsum[4]));
  assert_that(nsum[0], is_equal_to(nsum[5]));
}

Ensure(vec, dblN_ndot_works_on_bad_bb32_df_example) {
  dbl x[3] = {
    -0.16666666669920935,
     0.16666666663412383,
    -3.2542760153297934E-11
  };
  dbl y[3] = {
    0.33333333339841881,
    0.33333333320316222,
    0.33333333339841892
  };
  int perm[6][3] = {
    {0, 1, 2}, {0, 2, 1},
    {1, 0, 2}, {2, 0, 1},
    {1, 2, 0}, {2, 1, 0}
  };
  dbl xperm[3], yperm[3], ndot[6];
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      xperm[j] = x[perm[i][j]];
      yperm[j] = y[perm[i][j]];
    }
    ndot[i] = dblN_ndot(xperm, yperm, 3);
  }
  // See comment in `dblN_nsum_works_on_bad_bb32_df_example`.
  assert_that(ndot[0], is_equal_to(ndot[1]));
  assert_that(ndot[0], is_equal_to(ndot[2]));
  assert_that(ndot[0], is_equal_to(ndot[3]));
  assert_that(ndot[0], is_equal_to(ndot[4]));
  assert_that(ndot[0], is_equal_to(ndot[5]));
}
