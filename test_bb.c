#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "bb.h"
#include "mat.h"
#include "vec.h"

#define NUM_RANDOM_TRIALS 10

Describe(bb);

BeforeEach(bb) {
  significant_figures_for_assert_double_are(13);
}

AfterEach(bb) {}

Ensure (bb, bb3_has_quadratic_precision) {
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
  dbl c[4];
  bb3_interp(f, Df, x, c);

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
    f_bb = bb3(c, b);
    assert_that_double(f_gt, is_equal_to_double(f_bb));

    Df_gt = 2*A*y + B;
    Df_bb = dbb3(c, b, a)/(x[1] - x[0]);
    assert_that_double(Df_gt, is_equal_to_double(Df_bb));
  }

  gsl_rng_free(rng);
}

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

Ensure (bb, bb3tri_has_linear_precision) {
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

  dbl c[10];
  bb3tri_interp3(f, Df, x, c);

  dbl lam[3], y[3];
  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    get_random_lambda(rng, lam);
    get_conv_comb(x, lam, y);

    dbl f_gt = dbl3_dot(A, y) + B;
    dbl f_bb = bb3tri(c, lam);
    assert_that_double(f_gt, is_equal_to_double(f_bb));
  }
}

dbl q(dbl A[3][3], dbl b[3], dbl d, dbl const x[3]) {
  dbl Ax[3];
  dbl33_dbl3_mul(A, x, Ax);
  return dbl3_dot(Ax, x) + dbl3_dot(b, x) + d;
}

void Dq(dbl A[3][3], dbl b[3], dbl const x[3], dbl Dq[3]) {
  dbl At[3][3];
  dbl33_transposed(A, At);

  dbl Ax[3], Atx[3];
  dbl33_dbl3_mul(A, x, Ax);
  dbl33_dbl3_mul(At, x, Atx);

  dbl tmp[3];
  dbl3_add(Ax, Atx, tmp);

  dbl3_add(tmp, b, Dq);
}

Ensure (bb, bb3tri_has_quadratic_precision) {
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

  dbl c[10];
  bb3tri_interp3(f, Df, x, c);

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
    q_bb = bb3tri(c, lam);

    assert_that_double(q_gt, is_equal_to_double(q_bb));

    /**
     * Note: Df(xi)'(xj - xi) = Df(X*ei)'(ej - ei)
     */

    Dq(A, b, x[i], tmp);
    Dq_gt[0] = dbl3_dot(tmp, dx[0]);
    Dq_gt[1] = dbl3_dot(tmp, dx[1]);

    Dq_bb[0] = dbb3tri(c, lam, a[0]);
    Dq_bb[1] = dbb3tri(c, lam, a[1]);

    assert_that_double(Dq_gt[0], is_equal_to_double(Dq_bb[0]));
    assert_that_double(Dq_gt[1], is_equal_to_double(Dq_bb[1]));

    D2q_bb[0][0] = d2bb3tri(c, lam, a[0], a[0]);
    D2q_bb[1][0] = d2bb3tri(c, lam, a[1], a[0]);
    D2q_bb[0][1] = d2bb3tri(c, lam, a[0], a[1]);
    D2q_bb[1][1] = d2bb3tri(c, lam, a[1], a[1]);

    assert_that_double(D2q_gt[0][0], is_equal_to_double(D2q_bb[0][0]));
    assert_that_double(D2q_gt[1][0], is_equal_to_double(D2q_bb[1][0]));
    assert_that_double(D2q_gt[0][1], is_equal_to_double(D2q_bb[0][1]));
    assert_that_double(D2q_gt[1][1], is_equal_to_double(D2q_bb[1][1]));
  }

  /**
   * Check random locations elsewhere in the triangle
   */
  dbl y[3];
  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    get_random_lambda(rng, lam);
    get_conv_comb(x, lam, y);

    dbl q_gt = q(A, b, d, y);
    dbl q_bb = bb3tri(c, lam);

    assert_that_double(q_gt, is_equal_to_double(q_bb));
  }

  gsl_rng_free(rng);
}

Ensure(bb, bb3tri_works_for_simple_olim6_update) {
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

  dbl c[10];

  dbl a[2][3] = {
    {-1, 1, 0},
    {-1, 0, 1}
  };

  dbl Dfa[3][2] = {
    {dbl3_dot(dx[0], Df[0]), dbl3_dot(dx[1], Df[0])},
    {dbl3_dot(dx[0], Df[1]), dbl3_dot(dx[1], Df[1])},
    {dbl3_dot(dx[0], Df[2]), dbl3_dot(dx[1], Df[2])}
  };

  bb3tri_interp3(f, Df, x, c);

  {
    dbl b[3] = {1, 0, 0};
    assert_that_double(bb3tri(c, b), is_equal_to_double(1.0));
    assert_that_double(dbb3tri(c, b, a[0]), is_equal_to_double(Dfa[0][0]));
    assert_that_double(dbb3tri(c, b, a[1]), is_equal_to_double(Dfa[0][1]));
  }

  {
    dbl b[3] = {0, 1, 0};
    assert_that_double(bb3tri(c, b), is_equal_to_double(1.0));
    assert_that_double(dbb3tri(c, b, a[0]), is_equal_to_double(Dfa[1][0]));
    assert_that_double(dbb3tri(c, b, a[1]), is_equal_to_double(Dfa[1][1]));
  }

  {
    dbl b[3] = {0, 0, 1};
    assert_that_double(bb3tri(c, b), is_equal_to_double(1.0));
    assert_that_double(dbb3tri(c, b, a[0]), is_equal_to_double(Dfa[2][0]));
    assert_that_double(dbb3tri(c, b, a[1]), is_equal_to_double(Dfa[2][1]));
  }
}
