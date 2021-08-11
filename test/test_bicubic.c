#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define NUM_RANDOM_TRIALS 10

#include "bicubic.h"

Describe(bicubic);

BeforeEach(bicubic) {
  double_absolute_tolerance_is(1e-15);
  double_absolute_tolerance_is(1e-15);
}

AfterEach(bicubic) {}

Ensure(bicubic, get_fx_on_edge_works) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  bicubic_s bicubic;
  cubic_s cubic;

  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        bicubic.A[i][j] = gsl_ran_gaussian(rng, 1.0);

    cubic = bicubic_get_fx_on_edge(&bicubic, LAMBDA, 0);

    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      dbl lam = gsl_ran_flat(rng, 0, 1);
      dbl fx_cubic = cubic_f(&cubic, lam);
      dbl fx_bicubic = bicubic_fx(&bicubic, DBL2(lam, 0));
      assert_that_double(fx_cubic, is_nearly_double(fx_bicubic));
    }
  }

  gsl_rng_free(rng);
}
