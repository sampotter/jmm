#include <cgreen/cgreen.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "mat.h"

#define NUM_RANDOM_TRIALS 10

Describe(dbl44);

BeforeEach(dbl44) {
  double_absolute_tolerance_is(1e-13);
  double_relative_tolerance_is(1e-13);
}

AfterEach(dbl44) {}

int dbl44_det_gsl(dbl const A[4][4], dbl *det) {
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

    dbl det_gsl;
    int error = dbl44_det_gsl(A, &det_gsl);
    assert_that(!error);

    dbl det = dbl44_det(A);

    assert_that_double(det, is_nearly_double(det_gsl));
  }

  gsl_rng_free(rng);
}

int dbl44_dbl4_solve_gsl(dbl const A[4][4], dbl const b[4], dbl x[4]) {
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

  dbl A[4][4], b[4], x[4], x_gsl[4];

  for (int _ = 0; _ < NUM_RANDOM_TRIALS; ++_) {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        A[i][j] = gsl_ran_gaussian(rng, 1.0);
      }
      b[i] = gsl_ran_gaussian(rng, 1.0);
    }

    int error = dbl44_dbl4_solve_gsl(A, b, x_gsl);
    assert_that(!error);

    dbl44_dbl4_solve(A, b, x);

    for (int i = 0; i < 4; ++i) {
      assert_that_double(x[i], is_nearly_double(x_gsl[i]));
    }
  }

  gsl_rng_free(rng);
}

Ensure (dbl44, dbl4_mul_works) {
  dbl44 A = {
    {0.04440129622672145, -1.791725852886146, -0.5465298382048241, -0.8416259191899776},
    {-1.287411706678259, -1.774436226554693, 0.9067743971259874, -0.4538715589706994},
    {-2.658977126008932, -0.4894769507574523, 0.2398244941599099, 0.5530479398676308},
    {-1.664763207458297, 0.496870525703022, -1.210463986309432, -1.215804454782426}
  };

  dbl4 b = {-1.428737780331591, 0.05818709668635221, 0.2934308272841706, 0.03704906574396089};
  dbl4 x_gt = {-0.359243291423951, 1.985384496175262, 3.861361643515079, 2.007189675043842};
  dbl4 y_gt = {-0.9802523382021779, 2.331438087912691, 0.8591357373548978, 1.293290174714141};

  dbl4 x;
  dbl44_dbl4_mul(A, b, x);

  assert_that_double(x[0], is_nearly_double(x_gt[0]));
  assert_that_double(x[1], is_nearly_double(x_gt[1]));
  assert_that_double(x[2], is_nearly_double(x_gt[2]));
  assert_that_double(x[3], is_nearly_double(x_gt[3]));

  dbl4 y;
  dbl4_dbl44_mul(b, A, y);

  assert_that_double(y[0], is_nearly_double(y_gt[0]));
  assert_that_double(y[1], is_nearly_double(y_gt[1]));
  assert_that_double(y[2], is_nearly_double(y_gt[2]));
  assert_that_double(y[3], is_nearly_double(y_gt[3]));
}
