#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "bb.h"
#include "macros.h"
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

void get_random_affine_coefs(gsl_rng *rng, dbl a[3]) {
  /* Sample a standard normal Gaussian vector */
  for (size_t i = 0; i < 3; ++i)
    a[i] = gsl_ran_gaussian(rng, 1.0);

  /* Project it into the orthogonal complement of the span of the
   * vector (1, 1, 1). */
  dbl shift = dbl3_sum(a)/3;
  for (size_t i = 0; i < 3; ++i)
    a[i] -= shift;
}

void get_conv_comb(dbl const x[3][3], dbl const lam[3], dbl y[3]) {
  for (int i = 0; i < 3; ++i) {
    y[i] = x[0][i]*lam[0] + x[1][i]*lam[1] + x[2][i]*lam[2];
  }
}

Describe(bb32);

BeforeEach(bb32) {
  double_absolute_tolerance_is(2.5e-15);
  double_relative_tolerance_is(2.5e-15);
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

Ensure(bb32, df_is_symmetric) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  bb32 bb = {
    .c = {
      1, 0.66666666666666674, 0.66666666666666674, 1, 0.66666666666666674,
      0.50000000000000022, 0.66666666666666674, 0.66666666666666674,
      0.66666666666666674, 1}
  };

  dbl b[3], df[2], a[2][3] = {{-1, 1, 0}, {-1, 0, 1}};

  get_random_lambda(rng, b);

  b[0] = 0.33333333339841881;
  b[1] = 0.33333333320316222;
  b[2] = 0.33333333339841892;

  df[0] = bb32_df(&bb, b, a[0]);
  SWAP(b[1], b[2]);
  SWAP(a[0][1], a[0][2]);
  df[1] = bb32_df(&bb, b, a[0]);
  SWAP(b[1], b[2]);
  SWAP(a[0][1], a[0][2]);
  assert_that_double(df[0], is_nearly_double(df[1]));

  df[0] = bb32_df(&bb, b, a[1]);
  SWAP(b[1], b[2]);
  SWAP(a[1][1], a[1][2]);
  df[1] = bb32_df(&bb, b, a[1]);
  SWAP(b[1], b[2]);
  SWAP(a[1][1], a[1][2]);
  assert_that_double(df[0], is_nearly_double(df[1]));
}

Ensure(bb32, adjacent_bb32_are_C0) {
  jet3 jet1[3] = {
    {
      .f = 6.0515782990567839,
      .fx = 0.10494745590458605,
      .fy = -0.12623235128221608,
      .fz = -0.98643369011247684
    },
    {
      .f = 5.9505211776379818, .fx = NAN, .fy = NAN, .fz = NAN
    },
    {
      .f = 5.7756471295123735,
      .fx = 0.30240346235068588,
      .fy = 0.014147977152708482,
      .fz = -0.95307501315521004
    }
  };

  dbl x1[3][3] = {
    {1.3604636896622806, 2.8496775671458718, -5.1111864057935801},
    {1, 3.5, -5.1999988555908203},
    {1.4641119106968401, 3.4797503769549065, -4.8855490197739462}
  };

  bb32 bb1;
  bb32_init_from_jets(&bb1, jet1, x1);

  jet3 jet2[3] = {
    {
      .f = 6.0515782990567839,
      .fx = 0.10494745590458605,
      .fy = -0.12623235128221608,
      .fz = -0.98643369011247684
    },
    {
      .f = 5.8776740144675363,
      .fx = 0.42650997855742201,
      .fy = -0.23420586793585149,
      .fz = -0.87363427680887262
    },
    {
      .f = 5.7756471295123735,
      .fx = 0.30240346235068588,
      .fy = 0.014147977152708482,
      .fz = -0.95307501315521004
    }
  };

  dbl x2[3][3] = {
    {1.3604636896622806, 2.8496775671458718, -5.1111864057935801},
    {2.0820123257785261, 2.9290136105001272, -4.6698515448704896},
    {1.4641119106968401, 3.4797503769549065, -4.8855490197739462}
  };

  bb32 bb2;
  bb32_init_from_jets(&bb2, jet2, x2);

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  dbl b[3];
  b[1] = 0;

  for (size_t _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    b[0] = gsl_ran_flat(rng, 0, 1);
    b[2] = 1 - b[0];

    dbl f1 = bb32_f(&bb1, b);
    dbl f2 = bb32_f(&bb2, b);

    assert_that_double(f1, is_nearly_double(f2));
  }

  gsl_rng_free(rng);
}

Ensure(bb32, init_from_3d_data_and_init_from_jets_are_equivalent) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  dbl T[3], DT[3][3], x[3][3];
  jet3 jet[3];

  bb32 bb_from_3d_data, bb_from_jets;

  dbl b[3], a1[3], a2[3];
  dbl f_from_3d_data, df_from_3d_data, d2f_from_3d_data;
  dbl f_from_jets, df_from_jets, d2f_from_jets;

  for (size_t _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (size_t i = 0; i < 3; ++i)
      T[i] = gsl_ran_gaussian(rng, 1.0);

    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        DT[i][j] = gsl_ran_gaussian(rng, 1.0);

    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        x[i][j] = gsl_ran_gaussian(rng, 1.0);

    for (size_t i = 0; i < 3; ++i) {
      jet[i].f = T[i];
      jet[i].fx = DT[i][0];
      jet[i].fy = DT[i][1];
      jet[i].fz = DT[i][2];
    }

    bb32_init_from_3d_data(&bb_from_3d_data, T, DT, x);
    bb32_init_from_jets(&bb_from_jets, jet, x);

    dbl const *c_from_3d_data = bb_from_3d_data.c;
    dbl const *c_from_jets = bb_from_jets.c;

    for (size_t i = 0; i < 10; ++i)
      assert_that_double(c_from_3d_data[i], is_nearly_double(c_from_jets[i]));

    for (size_t __ = 0; __ < NUM_RANDOM_TRIALS; ++__) {
      get_random_lambda(rng, b);
      get_random_affine_coefs(rng, a1);
      get_random_affine_coefs(rng, a1);

      f_from_3d_data = bb32_f(&bb_from_3d_data, b);
      df_from_3d_data = bb32_df(&bb_from_3d_data, b, a1);
      d2f_from_3d_data = bb32_d2f(&bb_from_3d_data, b, a1, a2);

      f_from_jets = bb32_f(&bb_from_jets, b);
      df_from_jets = bb32_df(&bb_from_jets, b, a1);
      d2f_from_jets = bb32_d2f(&bb_from_jets, b, a1, a2);

      assert_that_double(f_from_3d_data, is_nearly_double(f_from_jets));
      assert_that_double(df_from_3d_data, is_nearly_double(df_from_jets));
      assert_that_double(d2f_from_3d_data, is_nearly_double(d2f_from_jets));
    }
  }

  gsl_rng_free(rng);
}
