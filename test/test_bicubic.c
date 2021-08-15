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

Ensure(bicubic, set_data_works) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  dbl44 data;
  bicubic_s bicubic;

  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        data[i][j] = gsl_ran_gaussian(rng, 1.0);

    bicubic_set_data(&bicubic, data);

    assert_that_double(bicubic_f(&bicubic, DBL2(0, 0)), is_nearly_double(data[0][0]));
    assert_that_double(bicubic_f(&bicubic, DBL2(0, 1)), is_nearly_double(data[0][1]));
    assert_that_double(bicubic_f(&bicubic, DBL2(1, 0)), is_nearly_double(data[1][0]));
    assert_that_double(bicubic_f(&bicubic, DBL2(1, 1)), is_nearly_double(data[1][1]));

    assert_that_double(bicubic_fx(&bicubic, DBL2(0, 0)), is_nearly_double(data[2][0]));
    assert_that_double(bicubic_fx(&bicubic, DBL2(0, 1)), is_nearly_double(data[2][1]));
    assert_that_double(bicubic_fx(&bicubic, DBL2(1, 0)), is_nearly_double(data[3][0]));
    assert_that_double(bicubic_fx(&bicubic, DBL2(1, 1)), is_nearly_double(data[3][1]));

    assert_that_double(bicubic_fy(&bicubic, DBL2(0, 0)), is_nearly_double(data[0][2]));
    assert_that_double(bicubic_fy(&bicubic, DBL2(0, 1)), is_nearly_double(data[0][3]));
    assert_that_double(bicubic_fy(&bicubic, DBL2(1, 0)), is_nearly_double(data[1][2]));
    assert_that_double(bicubic_fy(&bicubic, DBL2(1, 1)), is_nearly_double(data[1][3]));

    assert_that_double(bicubic_fxy(&bicubic, DBL2(0, 0)), is_nearly_double(data[2][2]));
    assert_that_double(bicubic_fxy(&bicubic, DBL2(0, 1)), is_nearly_double(data[2][3]));
    assert_that_double(bicubic_fxy(&bicubic, DBL2(1, 0)), is_nearly_double(data[3][2]));
    assert_that_double(bicubic_fxy(&bicubic, DBL2(1, 1)), is_nearly_double(data[3][3]));
  }

  gsl_rng_free(rng);
}

Ensure(bicubic, get_f_on_edge_works) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  bicubic_s bicubic;
  cubic_s cubic;
  dbl lam, mu, f_cubic, f_bicubic;

  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        bicubic.A[i][j] = gsl_ran_gaussian(rng, 1.0);

    cubic = bicubic_get_f_on_edge(&bicubic, LAMBDA, 0);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      lam = gsl_ran_flat(rng, 0, 1);
      f_cubic = cubic_f(&cubic, lam);
      f_bicubic = bicubic_f(&bicubic, DBL2(lam, 0));
      assert_that_double(f_cubic, is_nearly_double(f_bicubic));
    }

    cubic = bicubic_get_f_on_edge(&bicubic, LAMBDA, 1);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      lam = gsl_ran_flat(rng, 0, 1);
      f_cubic = cubic_f(&cubic, lam);
      f_bicubic = bicubic_f(&bicubic, DBL2(lam, 1));
      assert_that_double(f_cubic, is_nearly_double(f_bicubic));
    }

    cubic = bicubic_get_f_on_edge(&bicubic, MU, 0);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      mu = gsl_ran_flat(rng, 0, 1);
      f_cubic = cubic_f(&cubic, mu);
      f_bicubic = bicubic_f(&bicubic, DBL2(0, mu));
      assert_that_double(f_cubic, is_nearly_double(f_bicubic));
    }

    cubic = bicubic_get_f_on_edge(&bicubic, MU, 1);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      mu = gsl_ran_flat(rng, 0, 1);
      f_cubic = cubic_f(&cubic, mu);
      f_bicubic = bicubic_f(&bicubic, DBL2(1, mu));
      assert_that_double(f_cubic, is_nearly_double(f_bicubic));
    }
  }

  gsl_rng_free(rng);
}

Ensure(bicubic, get_fx_on_edge_works) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  bicubic_s bicubic;
  cubic_s cubic;
  dbl lam, mu, fx_cubic, fx_bicubic;

  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        bicubic.A[i][j] = gsl_ran_gaussian(rng, 1.0);

    cubic = bicubic_get_fx_on_edge(&bicubic, LAMBDA, 0);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      lam = gsl_ran_flat(rng, 0, 1);
      fx_cubic = cubic_f(&cubic, lam);
      fx_bicubic = bicubic_fx(&bicubic, DBL2(lam, 0));
      assert_that_double(fx_cubic, is_nearly_double(fx_bicubic));
    }

    cubic = bicubic_get_fx_on_edge(&bicubic, LAMBDA, 1);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      lam = gsl_ran_flat(rng, 0, 1);
      fx_cubic = cubic_f(&cubic, lam);
      fx_bicubic = bicubic_fx(&bicubic, DBL2(lam, 1));
      assert_that_double(fx_cubic, is_nearly_double(fx_bicubic));
    }

    cubic = bicubic_get_fx_on_edge(&bicubic, MU, 0);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      mu = gsl_ran_flat(rng, 0, 1);
      fx_cubic = cubic_f(&cubic, mu);
      fx_bicubic = bicubic_fx(&bicubic, DBL2(0, mu));
      assert_that_double(fx_cubic, is_nearly_double(fx_bicubic));
    }

    cubic = bicubic_get_fx_on_edge(&bicubic, MU, 1);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      mu = gsl_ran_flat(rng, 0, 1);
      fx_cubic = cubic_f(&cubic, mu);
      fx_bicubic = bicubic_fx(&bicubic, DBL2(1, mu));
      assert_that_double(fx_cubic, is_nearly_double(fx_bicubic));
    }
  }

  gsl_rng_free(rng);
}

Ensure(bicubic, get_fy_on_edge_works) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  bicubic_s bicubic;
  cubic_s cubic;
  dbl lam, mu, fy_cubic, fy_bicubic;

  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        bicubic.A[i][j] = gsl_ran_gaussian(rng, 1.0);

    cubic = bicubic_get_fy_on_edge(&bicubic, LAMBDA, 0);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      lam = gsl_ran_flat(rng, 0, 1);
      fy_cubic = cubic_f(&cubic, lam);
      fy_bicubic = bicubic_fy(&bicubic, DBL2(lam, 0));
      assert_that_double(fy_cubic, is_nearly_double(fy_bicubic));
    }

    cubic = bicubic_get_fy_on_edge(&bicubic, LAMBDA, 1);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      lam = gsl_ran_flat(rng, 0, 1);
      fy_cubic = cubic_f(&cubic, lam);
      fy_bicubic = bicubic_fy(&bicubic, DBL2(lam, 1));
      assert_that_double(fy_cubic, is_nearly_double(fy_bicubic));
    }

    cubic = bicubic_get_fy_on_edge(&bicubic, MU, 0);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      mu = gsl_ran_flat(rng, 0, 1);
      fy_cubic = cubic_f(&cubic, mu);
      fy_bicubic = bicubic_fy(&bicubic, DBL2(0, mu));
      assert_that_double(fy_cubic, is_nearly_double(fy_bicubic));
    }

    cubic = bicubic_get_fy_on_edge(&bicubic, MU, 1);
    for (int __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      mu = gsl_ran_flat(rng, 0, 1);
      fy_cubic = cubic_f(&cubic, mu);
      fy_bicubic = bicubic_fy(&bicubic, DBL2(1, mu));
      assert_that_double(fy_cubic, is_nearly_double(fy_bicubic));
    }
  }

  gsl_rng_free(rng);
}
