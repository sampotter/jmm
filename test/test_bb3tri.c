#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "bb.h"
#include "mat.h"
#include "util.h"
#include "vec.h"

#define NUM_RANDOM_TRIALS 10

void get_random_lambda(gsl_rng *rng, dbl lam[3]) {
  lam[1] = gsl_ran_flat(rng, 0, 1);
  lam[2] = gsl_ran_flat(rng, 0, 1);
  if (lam[1] + lam[2] > 1) {
    lam[1] = 1 - lam[1];
    lam[2] = 1 - lam[2];
  }
  lam[0] = 1 - lam[1] - lam[2];
}

void get_conv_comb(dbl const x[3][3], dbl const lam[3], dbl y[3]) {
  for (int i = 0; i < 3; ++i) {
    y[i] = x[0][i]*lam[0] + x[1][i]*lam[1] + x[2][i]*lam[2];
  }
}

Describe(bb32);

BeforeEach(bb32) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(bb32) {}

Ensure (bb32, has_linear_precision) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  /**
   * Sample a random linear function (domain is R3).
   */

  dbl A[3];
  for (int i = 0; i < 3; ++i) {
    A[i] = gsl_ran_gaussian(rng, 1.0);
  }

  dbl B = gsl_ran_gaussian(rng, 1.0);

  dbl x[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      x[i][j] = gsl_ran_gaussian(rng, 1.0);
    }
  }

  dbl f[3] = {
    dbl3_dot(A, x[0]) + B,
    dbl3_dot(A, x[1]) + B,
    dbl3_dot(A, x[2]) + B
  };

  dbl Df[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Df[i][j] = A[j];
    }
  }

  bb32 bb;
  bb32_init_from_3d_data(&bb, f, Df, x);

  dbl lam[3], y[3];
  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    get_random_lambda(rng, lam);
    get_conv_comb(x, lam, y);

    dbl f_gt = dbl3_dot(A, y) + B;
    dbl f_bb = bb32_f(&bb, lam);
    assert_that_double(f_gt, is_nearly_double(f_bb));
  }
}

Ensure (bb32, has_quadratic_precision) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  /**
   * Sample a random quadratic function (domain is R3).
   */

  dbl A[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      A[i][j] = gsl_ran_gaussian(rng, 1.0);
    }
  }

  dbl b[3];
  for (int i = 0; i < 3; ++i) {
    b[i] = gsl_ran_gaussian(rng, 1.0);
  }

  dbl d = gsl_ran_gaussian(rng, 1.0);

  dbl x[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      x[i][j] = gsl_ran_gaussian(rng, 1.0);
    }
  }

  dbl dx[2][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);

  dbl f[3] = {
    q(A, b, d, x[0]),
    q(A, b, d, x[1]),
    q(A, b, d, x[2])
  };

  dbl Df[3][3];
  for (int i = 0; i < 3; ++i) {
    Dq(A, b, x[i], Df[i]);
  }

  bb32 bb;
  bb32_init_from_3d_data(&bb, f, Df, x);

  dbl a[2][3] = {
    {-1, 1, 0},
    {-1, 0, 1}
  };

  dbl lam[3];
  dbl q_gt, q_bb, tmp[3], Dq_gt[2], Dq_bb[2];

  dbl At[3][3], A_plus_At[3][3];
  dbl33_transposed(A, At);
  dbl33_add(A, At, A_plus_At);

  dbl D2q_gt[2][2];
  dbl33_dbl3_mul(A_plus_At, dx[0], tmp);
  D2q_gt[0][0] = dbl3_dot(tmp, dx[0]);
  D2q_gt[0][1] = D2q_gt[1][0] = dbl3_dot(tmp, dx[1]);
  dbl33_dbl3_mul(A_plus_At, dx[1], tmp);
  D2q_gt[1][1] = dbl3_dot(tmp, dx[1]);

  dbl D2q_bb[2][2];

  /**
   * Check nodal values first
   */
  for (int i = 0; i < 3; ++i) {
    lam[0] = lam[1] = lam[2] = 0;
    lam[i] = 1;

    q_gt = q(A, b, d, x[i]);
    q_bb = bb32_f(&bb, lam);

    assert_that_double(q_gt, is_nearly_double(q_bb));

    /**
     * Note: Df(xi)'(xj - xi) = Df(X*ei)'(ej - ei)
     */

    Dq(A, b, x[i], tmp);
    Dq_gt[0] = dbl3_dot(tmp, dx[0]);
    Dq_gt[1] = dbl3_dot(tmp, dx[1]);

    Dq_bb[0] = bb32_df(&bb, lam, a[0]);
    Dq_bb[1] = bb32_df(&bb, lam, a[1]);

    assert_that_double(Dq_gt[0], is_nearly_double(Dq_bb[0]));
    assert_that_double(Dq_gt[1], is_nearly_double(Dq_bb[1]));

    D2q_bb[0][0] = bb32_d2f(&bb, lam, a[0], a[0]);
    D2q_bb[1][0] = bb32_d2f(&bb, lam, a[1], a[0]);
    D2q_bb[0][1] = bb32_d2f(&bb, lam, a[0], a[1]);
    D2q_bb[1][1] = bb32_d2f(&bb, lam, a[1], a[1]);

    assert_that_double(D2q_gt[0][0], is_nearly_double(D2q_bb[0][0]));
    assert_that_double(D2q_gt[1][0], is_nearly_double(D2q_bb[1][0]));
    assert_that_double(D2q_gt[0][1], is_nearly_double(D2q_bb[0][1]));
    assert_that_double(D2q_gt[1][1], is_nearly_double(D2q_bb[1][1]));
  }

  /**
   * Check random locations elsewhere in the triangle
   */
  dbl y[3];
  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    get_random_lambda(rng, lam);
    get_conv_comb(x, lam, y);

    dbl q_gt = q(A, b, d, y);
    dbl q_bb = bb32_f(&bb, lam);

    assert_that_double(q_gt, is_nearly_double(q_bb));
  }

  gsl_rng_free(rng);
}

Ensure(bb32, works_for_simple_olim6_update) {
  dbl f[3] = {1, 1, 1};

  dbl Df[3][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
  };

  dbl x[3][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
  };

  dbl dx[2][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);

  dbl a[2][3] = {
    {-1, 1, 0},
    {-1, 0, 1}
  };

  dbl Dfa[3][2] = {
    {dbl3_dot(dx[0], Df[0]), dbl3_dot(dx[1], Df[0])},
    {dbl3_dot(dx[0], Df[1]), dbl3_dot(dx[1], Df[1])},
    {dbl3_dot(dx[0], Df[2]), dbl3_dot(dx[1], Df[2])}
  };

  bb32 bb;
  bb32_init_from_3d_data(&bb, f, Df, x);

  {
    dbl b[3] = {1, 0, 0};
    assert_that_double(bb32_f(&bb, b), is_nearly_double(1.0));
    assert_that_double(bb32_df(&bb, b, a[0]), is_nearly_double(Dfa[0][0]));
    assert_that_double(bb32_df(&bb, b, a[1]), is_nearly_double(Dfa[0][1]));
  }

  {
    dbl b[3] = {0, 1, 0};
    assert_that_double(bb32_f(&bb, b), is_nearly_double(1.0));
    assert_that_double(bb32_df(&bb, b, a[0]), is_nearly_double(Dfa[1][0]));
    assert_that_double(bb32_df(&bb, b, a[1]), is_nearly_double(Dfa[1][1]));
  }

  {
    dbl b[3] = {0, 0, 1};
    assert_that_double(bb32_f(&bb, b), is_nearly_double(1.0));
    assert_that_double(bb32_df(&bb, b, a[0]), is_nearly_double(Dfa[2][0]));
    assert_that_double(bb32_df(&bb, b, a[1]), is_nearly_double(Dfa[2][1]));
  }
}
