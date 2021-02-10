#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "bb.h"
#include "mat.h"
#include "util.h"
#include "vec.h"

#define NUM_RANDOM_TRIALS 10

#define SET_UP_RANDOM_CONVEX_COMBINATION()                              \
  dbl b[4], xb[3];                                                      \
  do {                                                                  \
    for (int j = 1; j < 4; ++j) {                                       \
      b[j] = gsl_ran_flat(rng, 0, 1);                                   \
    }                                                                   \
    b[0] = 1 - dbl3_sum(&b[1]);                                         \
  } while (!dbl4_nonneg(b));                                            \
  xb[0] = b[0]*x[0][0] + b[1]*x[1][0] + b[2]*x[2][0] + b[3]*x[3][0];    \
  xb[1] = b[0]*x[0][1] + b[1]*x[1][1] + b[2]*x[2][1] + b[3]*x[3][1];    \
  xb[2] = b[0]*x[0][2] + b[1]*x[1][2] + b[2]*x[2][2] + b[3]*x[3][2];

Describe(bb33);

BeforeEach(bb33) {
  double_absolute_tolerance_is(5e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(bb33) {}

Ensure (bb33, has_linear_precision) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  /**
   * Sample a random linear function (domain is R3).
   */

  dbl A[3];
  for (int i = 0; i < 3; ++i) {
    A[i] = gsl_ran_gaussian(rng, 1.0);
  }

  dbl B = gsl_ran_gaussian(rng, 1.0);

  dbl x[4][3];
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 3; ++j) {
      x[i][j] = gsl_ran_gaussian(rng, 1.0);
    }
  }

  dbl f[4] = {
    dbl3_dot(A, x[0]) + B,
    dbl3_dot(A, x[1]) + B,
    dbl3_dot(A, x[2]) + B,
    dbl3_dot(A, x[3]) + B
  };

  dbl Df[4][3];
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 3; ++j) {
      Df[i][j] = A[j];
    }
  }

  bb33 bb;
  bb33_init_from_3d_data(&bb, f, Df, x);

  dbl f_gt, f_bb;
  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    SET_UP_RANDOM_CONVEX_COMBINATION();

    f_gt = dbl3_dot(A, xb) + B;
    f_bb = bb33_f(&bb, b);
    assert_that_double(f_gt, is_nearly_double(f_bb));
  }
}

Ensure (bb33, has_quadratic_precision) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  dbl qA[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      qA[i][j] = gsl_ran_gaussian(rng, 1.0);
    }
  }

  dbl qB[3];
  for (int i = 0; i < 3; ++i) {
    qB[i] = gsl_ran_gaussian(rng, 1.0);
  }

  dbl qC = gsl_ran_gaussian(rng, 1.0);

  dbl qAt[3][3], D2q[3][3];
  dbl33_transposed(qA, qAt);
  dbl33_add(qA, qAt, D2q);

  dbl x[4][3];
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 3; ++j) {
      x[i][j] = gsl_ran_gaussian(rng, 1.0);
    }
  }

  dbl f[4];
  for (int i = 0; i < 4; ++i) {
    f[i] = q(qA, qB, qC, x[i]);
  }

  dbl Df[4][3];
  for (int i = 0; i < 4; ++i) {
    Dq(qA, qB, x[i], Df[i]);
  }

  bb33 bb;
  bb33_init_from_3d_data(&bb, f, Df, x);

  dbl f_gt, f_bb, Df_gt, Df_bb, D2f_gt, D2f_bb;

  /**
   * Check values at tetrahedron vertices
   */
  for (int i = 0; i < 4; ++i) {
    dbl b[4] = {0, 0, 0, 0};
    b[i] = 1;

    // Check function values
    f_gt = q(qA, qB, qC, x[i]);
    f_bb = bb33_f(&bb, b);
    assert_that_double(f_gt, is_nearly_double(f_bb));

    // Check directional derivatives
    dbl Dfi[3], dx[2][3], a[2][4], tmp_[3];
    Dq(qA, qB, x[i], Dfi);
    for (int j = 0; j < 4; ++j) {
      if (i == j) continue;

      a[0][0] = a[0][1] = a[0][2] = a[0][3] = 0;
      a[0][i] = -1;
      a[0][j] = 1;

      dbl3_sub(x[j], x[i], dx[0]);

      Df_gt = dbl3_dot(Dfi, dx[0]);
      Df_bb = bb33_df(&bb, b, a[0]);
      assert_that_double(Df_gt, is_nearly_double(Df_bb));

      // Check mixed directional derivatives
      dbl33_dbl3_mul(D2q, dx[0], tmp_);
      for (int k = 0; k < 4; ++k) {
        if (i == k) continue;

        a[1][0] = a[1][1] = a[1][2] = a[1][3] = 0;
        a[1][i] = -1;
        a[1][k] = 1;

        dbl3_sub(x[k], x[i], dx[1]);

        D2f_gt = dbl3_dot(tmp_, dx[1]);
        D2f_bb = bb33_d2f(&bb, b, a);
        assert_that_double(D2f_gt, is_nearly_double(D2f_bb));
      }
    }
  }

  /**
   * Check random locations in the tetrahedron
   */
  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    SET_UP_RANDOM_CONVEX_COMBINATION();

    dbl Dfb[3], dx[2][3], a[2][4], tmp_[3];

    // Check function values
    f_gt = q(qA, qB, qC, xb);
    f_bb = bb33_f(&bb, b);
    assert_that_double(f_gt, is_nearly_double(f_bb));

    Dq(qA, qB, xb, Dfb);
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (i == j) continue;

        // Check directional derivatives
        a[0][0] = a[0][1] = a[0][2] = a[0][3] = 0;
        a[0][i] = -1;
        a[0][j] = 1;
        dbl3_sub(x[j], x[i], dx[0]);
        Df_gt = dbl3_dot(Dfb, dx[0]);
        Df_bb = bb33_df(&bb, b, a[0]);
        assert_that_double(Df_gt, is_nearly_double(Df_bb));

        dbl33_dbl3_mul(D2q, dx[0], tmp_);

        for (int k = 0; k < 4; ++k) {
          if (i == k) continue;

          // Check mixed directional derivatives
          a[1][0] = a[1][1] = a[1][2] = a[1][3] = 0;
          a[1][i] = -1;
          a[1][k] = 1;
          dbl3_sub(x[k], x[i], dx[1]);
          D2f_gt = dbl3_dot(tmp_, dx[1]);
          D2f_bb = bb33_d2f(&bb, b, a);
          assert_that_double(D2f_gt, is_nearly_double(D2f_bb));
        }
      }
    }
  }

  gsl_rng_free(rng);
}
