#include <cgreen/cgreen.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "mat.h"

#define NUM_RANDOM_TRIALS 10

Describe(dbl44);

BeforeEach(dbl44) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(dbl44) {}

int solve_gt(dbl const A[4][4], dbl const b[4], dbl x[4]) {
  int error;

  gsl_matrix *A_gsl = gsl_matrix_alloc(4, 4);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      gsl_matrix_set(A_gsl, i, j, A[i][j]);
    }
  }

  gsl_vector *b_gsl = gsl_vector_alloc(4);
  for (int i = 0; i < 4; ++i) {
    gsl_vector_set(b_gsl, i, b[i]);
  }

  gsl_permutation *p_gsl = gsl_permutation_alloc(4);
  int signum;
  error = gsl_linalg_LU_decomp(A_gsl, p_gsl, &signum);

  gsl_vector *x_gsl = gsl_vector_alloc(4);
  error = gsl_linalg_LU_solve(A_gsl, p_gsl, b_gsl, x_gsl);

  for (int i = 0; i < 4; ++i) {
    x[i] = gsl_vector_get(x_gsl, i);
  }

  gsl_vector_free(x_gsl);
  gsl_permutation_free(p_gsl);
  gsl_vector_free(b_gsl);
  gsl_matrix_free(A_gsl);

  return error;
}

Ensure (dbl44, dbl4_solve_works) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  dbl A[4][4], b[4], x[4], x_gt[4];

  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        A[i][j] = gsl_ran_gaussian(rng, 1.0);
      }
      b[i] = gsl_ran_gaussian(rng, 1.0);
    }

    int error = solve_gt(A, b, x_gt);
    assert_that(!error);

    dbl44_dbl4_solve(A, b, x);

    for (int i = 0; i < 4; ++i) {
      assert_that_double(x[i], is_nearly_double(x_gt[i]));
    }
  }

  gsl_rng_free(rng);
}

int get_det_gt(dbl const A[4][4], dbl *det) {
  int error;

  gsl_matrix *A_gsl = gsl_matrix_alloc(4, 4);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      gsl_matrix_set(A_gsl, i, j, A[i][j]);
    }
  }

  gsl_permutation *p_gsl = gsl_permutation_alloc(4);
  int signum;
  error = gsl_linalg_LU_decomp(A_gsl, p_gsl, &signum);

  *det = gsl_linalg_LU_det(A_gsl, signum);

  gsl_permutation_free(p_gsl);
  gsl_matrix_free(A_gsl);

  return error;
}

Ensure (dbl44, det_works) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  dbl A[4][4];

  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        A[i][j] = gsl_ran_gaussian(rng, 1.0);
      }
    }

    dbl det_gt;
    int error = get_det_gt(A, &det_gt);
    assert_that(!error);

    dbl det = dbl44_det(A);

    assert_that_double(det, is_nearly_double(det_gt));
  }

  gsl_rng_free(rng);
}
