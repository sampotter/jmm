#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "bb.h"
#include "vec.h"

#define NUM_RANDOM_TRIALS 10

Describe(bb3);

BeforeEach(bb3) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(bb3) {}

Ensure (bb3, has_quadratic_precision) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  /**
   * Endpoints of interval in original domain.
   */
  dbl x[2] = {
    gsl_ran_gaussian(rng, 1.0),
    gsl_ran_gaussian(rng, 1.0)
  };

  /**
   * Coefficients of quadratic.
   */
  dbl A = gsl_ran_gaussian(rng, 1.0);
  dbl B = gsl_ran_gaussian(rng, 1.0);
  dbl C = gsl_ran_gaussian(rng, 1.0);

  /**
   * Compute Bezier ordinates for interpolant.
   */
  dbl f[2] = {A*x[0]*x[0] + B*x[0] + C, A*x[1]*x[1] + B*x[1] + C};
  dbl Df[2] = {2*A*x[0] + B, 2*A*x[1] + B};
  bb3 bb;
  bb3_init_from_1d_data(&bb, f, Df, x);

  /**
   * Do some random tests.
   */
  dbl a[2] = {-1, 1};
  dbl y, f_gt, f_bb, Df_gt, Df_bb, b[2];
  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    b[0] = gsl_ran_flat(rng, 0, 1);
    b[1] = 1 - b[0];
    y = dbl2_dot(b, x);

    f_gt = A*y*y + B*y + C;
    f_bb = bb3_f(&bb, b);
    assert_that_double(f_gt, is_nearly_double(f_bb));

    Df_gt = 2*A*y + B;
    Df_bb = bb3_df(&bb, b, a)/(x[1] - x[0]);
    assert_that_double(Df_gt, is_nearly_double(Df_bb));
  }

  gsl_rng_free(rng);
}
