#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "mat.h"
#include "opt.h"

#define NUM_RANDOM_TRIALS 10

Describe(triqp2);

BeforeEach(triqp2) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(triqp2) {}

Ensure (triqp2, solve_works) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  triqp2_s qp;

  /**
   * Use this Hessian for all of the test problems in this set of
   * tests.
   */
  qp.A[0][0] = 3;
  qp.A[1][0] = qp.A[0][1] = 1;
  qp.A[1][1] = 2;

  dbl x[2];

  // Interior point minimizers

  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {

    x[0] = gsl_ran_flat(rng, 0, 1);
    x[1] = gsl_ran_flat(rng, 0, 1);
    if (x[0] + x[1] > 1) {
      x[0] = 1 - x[0];
      x[1] = 1 - x[1];
    }

    dbl22_dbl2_mul(qp.A, x, qp.b);
    dbl2_negate(qp.b);

    triqp2_solve(&qp);

    assert_that_double(qp.x[0], is_nearly_double(x[0]));
    assert_that_double(qp.x[1], is_nearly_double(x[1]));
  }

  // Vertex minimizers

  qp.b[0] = qp.b[1] = 0;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(0));

  qp.b[0] = -3;
  qp.b[1] = -1;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(1));
  assert_that_double(qp.x[1], is_nearly_double(0));

  qp.b[0] = -1;
  qp.b[1] = -2;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(1));

  // Minimizers on the boundary

  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    x[0] = gsl_ran_flat(rng, 0, 1);
    x[1] = 0;

    dbl22_dbl2_mul(qp.A, x, qp.b);
    dbl2_negate(qp.b);

    triqp2_solve(&qp);

    assert_that_double(qp.x[0], is_nearly_double(x[0]));
    assert_that_double(qp.x[1], is_nearly_double(x[1]));
  }

  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    x[0] = 0;
    x[1] = gsl_ran_flat(rng, 0, 1);

    dbl22_dbl2_mul(qp.A, x, qp.b);
    dbl2_negate(qp.b);

    triqp2_solve(&qp);

    assert_that_double(qp.x[0], is_nearly_double(x[0]));
    assert_that_double(qp.x[1], is_nearly_double(x[1]));
  }

  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    x[0] = gsl_ran_flat(rng, 0, 1);
    x[1] = 1 - x[0];

    dbl22_dbl2_mul(qp.A, x, qp.b);
    dbl2_negate(qp.b);

    triqp2_solve(&qp);

    assert_that_double(qp.x[0], is_nearly_double(x[0]));
    assert_that_double(qp.x[1], is_nearly_double(x[1]));
  }

  // Minimizers with nonzero Lagrange multipliers

  qp.b[0] = 2;
  qp.b[1] = 1.5;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(0));

  qp.b[0] = -1;
  qp.b[1] = 0.5;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(1./3));
  assert_that_double(qp.x[1], is_nearly_double(0));

  qp.b[0] = 1;
  qp.b[1] = -0.5;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(0.25));

  qp.b[0] = -3.5;
  qp.b[1] = -2;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(5./6));
  assert_that_double(qp.x[1], is_nearly_double(1./6));

  gsl_rng_free(rng);
}

Ensure(triqp2, collected_problems) {
  triqp2_s qp;

  qp.A[0][0] = 3.7455872107577894;
  qp.A[0][1] = qp.A[1][0] = 1.8797755666030906;
  qp.A[1][1] = 1.4232583973466495;

  qp.b[0] = -1.8806971755346154;
  qp.b[1] = -0.71064019331246775;

  triqp2_solve(&qp);

  assert_that_double(qp.x[0], is_nearly_double(0.5021101017573482));
  assert_that_double(qp.x[1], is_nearly_double(0));

  qp.A[0][0] = 0.27593389435955507;
  qp.A[0][1] = qp.A[1][0] = -0.00022250249910501976;
  qp.A[1][1] = 0.30352699961859408;

  qp.b[0] = -0.27593389435955518;
  qp.b[1] = 0.00022250249910489028;

  triqp2_solve(&qp);

  assert_that_double(qp.x[0], is_nearly_double(1));
  assert_that_double(qp.x[1], is_nearly_double(0));
}
