#include "mat.h"

#include <assert.h>
#include <string.h>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>

#include "macros.h"

bool dbl22_isfinite(dbl22 const A) {
  return isfinite(A[0][0]) &&
         isfinite(A[1][0]) &&
         isfinite(A[0][1]) &&
         isfinite(A[1][1]);
}

bool dbl22_all_nan(dbl22 const A) {
  return isnan(A[0][0]) && isnan(A[0][1]) && isnan(A[1][0]) && isnan(A[1][1]);
}

dbl dbl22_det(dbl22 const X) {
  return X[0][0]*X[1][1] - X[1][0]*X[0][1];
}

dbl dbl22_trace(dbl22 const A) {
  return A[0][0] + A[1][1];
}

void dbl22_add(dbl22 const A, dbl22 const B, dbl22 C) {
  C[0][0] = A[0][0] + B[0][0];
  C[0][1] = A[0][1] + B[0][1];
  C[1][0] = A[1][0] + B[1][0];
  C[1][1] = A[1][1] + B[1][1];
}

void dbl22_add_inplace(dbl22 A, dbl22 const B) {
  A[0][0] += B[0][0];
  A[0][1] += B[0][1];
  A[1][0] += B[1][0];
  A[1][1] += B[1][1];
}

void dbl22_copy(dbl22 const A, dbl22 B) {
  memcpy(B, A, sizeof(dbl22));
}

void dbl22_dbl2_mul(dbl22 const A, dbl2 const x, dbl2 b) {
  b[0] = A[0][0]*x[0] + A[0][1]*x[1];
  b[1] = A[1][0]*x[0] + A[1][1]*x[1];
}

void dbl22_dbl2_solve(dbl22 const A, dbl2 const b, dbl2 x) {
  dbl det = dbl22_det(A);
  x[0] = (A[1][1]*b[0] - A[1][0]*b[1])/det;
  x[1] = (A[0][0]*b[1] - A[0][1]*b[0])/det;
}

void dbl22_dbl_div(dbl22 const A, dbl a, dbl22 B) {
  B[0][0] = A[0][0]/a;
  B[0][1] = A[0][1]/a;
  B[1][0] = A[1][0]/a;
  B[1][1] = A[1][1]/a;
}

void dbl22_dbl_div_inplace(dbl22 A, dbl a) {
  A[0][0] /= a;
  A[0][1] /= a;
  A[1][0] /= a;
  A[1][1] /= a;
}

void dbl22_dbl_mul(dbl22 const A, dbl a, dbl22 B) {
  B[0][0] = a*A[0][0];
  B[0][1] = a*A[0][1];
  B[1][0] = a*A[1][0];
  B[1][1] = a*A[1][1];
}

void dbl22_dbl_mul_inplace(dbl22 A, dbl a) {
  A[0][0] *= a;
  A[0][1] *= a;
  A[1][0] *= a;
  A[1][1] *= a;
}

void dbl22_eigvals(dbl22 const A, dbl2 lam) {
  dbl tr = dbl22_trace(A);
  dbl disc = tr*tr - 4*dbl22_det(A);
  assert(disc >= 0);
  dbl tmp = sqrt(disc);
  lam[0] = (tr + tmp)/2;
  lam[1] = (tr - tmp)/2;
}

void dbl22_eye(dbl22 eye) {
  eye[0][0] = eye[1][1] = 1;
  eye[1][0] = eye[0][1] = 0;
}

void dbl22_inv(dbl22 const X, dbl22 Y) {
  dbl det = dbl22_det(X);
  Y[0][0] = X[1][1]/det;
  Y[1][1] = X[0][0]/det;
  Y[1][0] = -X[1][0]/det;
  Y[0][1] = -X[0][1]/det;
}

void dbl22_invert(dbl22 X) {
  dbl det = dbl22_det(X);
  SWAP(X[0][0], X[1][1]);
  X[1][0] *= -1;
  X[0][1] *= -1;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      X[i][j] /= det;
}

void dbl22_mul(dbl22 const A, dbl22 const B, dbl22 C) {
  dbl22_zero(C);
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k)
        C[i][j] += A[i][k]*B[k][j];
}

void dbl22_perturb(dbl22 A, dbl eps) {
  A[0][0] += eps;
  A[1][1] += eps;
}

void dbl22_saxpy(dbl a, dbl22 const X, dbl22 const Y, dbl22 Z) {
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      Z[i][j] = a*X[i][j] + Y[i][j];
}

void dbl22_sub(dbl22 const A, dbl22 const B, dbl22 C) {
  C[0][0] = A[0][0] - B[0][0];
  C[0][1] = A[0][1] - B[0][1];
  C[1][0] = A[1][0] - B[1][0];
  C[1][1] = A[1][1] - B[1][1];
}

void dbl22_transpose(dbl22 A) {
  SWAP(A[0][1], A[1][0]);
}

void dbl22_zero(dbl22 A) {
  A[0][0] = A[1][0] = A[0][1] = A[1][1] = 0;
}

void dbl2_outer(dbl2 const u, dbl2 const v, dbl22 uv) {
  uv[0][0] = u[0]*v[0];
  uv[0][1] = u[0]*v[1];
  uv[1][0] = u[1]*v[0];
  uv[1][1] = u[1]*v[1];
}

bool dbl33_isfinite(dbl33 const A) {
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (!isfinite(A[i][j]))
        return false;
  return true;
}

bool dbl33_isnan(dbl33 const A) {
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (!isnan(A[i][j]))
        return false;
  return true;
}

dbl dbl33_det(dbl33 const A) {
  return A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
       - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
       + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

dbl dbl3_dbl33_dbl3_dot(dbl3 const x, dbl33 const A, dbl3 const y) {
  dbl tmp = 0;
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      tmp += x[i]*A[i][j]*y[j];
  return tmp;
}

void dbl33_add(dbl33 const A, dbl33 const B, dbl33 C) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }
}

void dbl33_add_inplace(dbl33 A, dbl33 const B) {
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      A[i][j] += B[i][j];
}

void dbl33_copy(dbl33 const in, dbl33 out) {
  memcpy(out, in, sizeof(dbl33));
}

void dbl33_dbl3_mul(dbl33 const A, dbl3 const x, dbl3 b) {
  dbl3_zero(b);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      b[i] += A[i][j]*x[j];
}

void dbl33_dbl3_mul_inplace(dbl33 const A, dbl3 x) {
  dbl3 b;
  dbl33_dbl3_mul(A, x, b);
  dbl3_copy(b, x);
}

void dbl33_dbl3_nmul(dbl33 const A, dbl3 const x, dbl3 b) {
  for (int i = 0; i < 3; ++i)
    b[i] = dbl3_ndot(A[i], x);
}

void dbl33_dbl3_solve(dbl33 const A, dbl3 const b, dbl3 x) {
  dbl33 LU;
  memcpy(LU, A, 3*3*sizeof(dbl));

  dbl const n = 3;

  int perm[3] = {0, 1, 2};

  // LU decomposition with partial pivoting
  for (int k = 0; k < n - 1; ++k) {
    // Find max |Ui:| for i >= k
    int argi = k;
    dbl absmax = fabs(LU[k][k]);
    for (int i = k + 1; i < n; ++i) {
      dbl tmp = fabs(LU[i][k]);
      if (tmp > absmax) {
        absmax = tmp;
        argi = i;
      }
    }

    if (argi != k) {
      // Swap rows i and k of L and U
      for (int j = 0; j < n; ++j)
        SWAP(LU[k][j], LU[argi][j]);

      // Swap entries i and k of permutation
      SWAP(perm[k], perm[argi]);
    }

    if (fabs(LU[k][k]) < 1e-15)
      continue;

    // Normalize column k of L
    for (int i = k + 1; i < n; ++i)
      LU[i][k] /= LU[k][k];

    // Subtract outer product
    for (int i = k + 1; i < n; ++i)
      for (int j = k + 1; j < n; ++j)
        LU[i][j] -= LU[i][k]*LU[k][j];
  }

  // Solve Ly = Pb
  for (int i = 0; i < n; ++i) {
    x[i] = b[perm[i]];
    for (int j = 0; j < i; ++j) {
      x[i] -= LU[i][j]*x[j];
    }
  }

  // Solve Ux = y
  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j) {
      x[i] -= LU[i][j]*x[j];
    }
    x[i] /= LU[i][i];
  }
}

void dbl33_dbl_div(dbl33 const A, dbl a, dbl33 B) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      B[i][j] = A[i][j]/a;
}

void dbl33_dbl_div_inplace(dbl33 A, dbl a) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      A[i][j] /= a;
}

void dbl33_eigvals_sym(dbl33 const A, dbl3 lam) {
  gsl_matrix *A_gsl = gsl_matrix_alloc(3, 3);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      gsl_matrix_set(A_gsl, i, j, (A[i][j] + A[j][i])/2);

  gsl_vector *lam_gsl = gsl_vector_alloc(3);

  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(3);
  gsl_eigen_symm(A_gsl, lam_gsl, w);

  for (size_t i = 0; i < 3; ++i)
    lam[i] = gsl_vector_get(lam_gsl, i);

  gsl_eigen_symm_free(w);
  gsl_vector_free(lam_gsl);
  gsl_matrix_free(A_gsl);
}

void dbl33_eye(dbl33 A) {
  A[0][0] = 1; A[0][1] = 0; A[0][2] = 0;
  A[1][0] = 0; A[1][1] = 1; A[1][2] = 0;
  A[2][0] = 0; A[2][1] = 0; A[2][2] = 1;
}

void dbl33_get_column( dbl33 const A, int i, dbl3 x) {
  x[0] = A[0][i];
  x[1] = A[1][i];
  x[2] = A[2][i];
}

void dbl33_invert(dbl33 A) {
  dbl33 eye;
  dbl33_eye(eye);

  // TODO: very inefficient: should do LU with partial pivoting with
  // multiple RHSs
  dbl33 tmp;
  for (size_t i = 0; i < 3; ++i)
    dbl33_dbl3_solve(A, eye[i], tmp[i]);

  dbl33_transposed(tmp, A);
}

void dbl33_make_axis_angle_rotation_matrix(dbl3 axis, dbl angle, dbl33 rot) {
  assert(fabs(1 - dbl3_norm(axis)) < 1e-14);

  dbl x = axis[0], y = axis[1], z = axis[2];
  dbl c = cos(angle);
  dbl s = sin(angle);

  dbl ci = 1 - c;

  rot[0][0] = ci*x*x + c;
  rot[0][1] = ci*x*y - z*s;
  rot[0][2] = ci*x*z + y*s;

  rot[1][0] = ci*y*x + z*s;
  rot[1][1] = ci*y*y + c;
  rot[1][2] = ci*y*z - x*s;

  rot[2][0] = ci*z*x - y*s;
  rot[2][1] = ci*z*y + x*s;
  rot[2][2] = ci*z*z + c;
}

void dbl33_mul(dbl33 const A, dbl33 const B, dbl33 C) {
  dbl33_zero(C);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        C[i][j] += A[i][k]*B[k][j];
}

void dbl33_nan(dbl33 A) {
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      A[i][j] = NAN;
}

void dbl33_saxpy(dbl a, dbl33 const X, dbl33 const Y, dbl33 Z) {
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      Z[i][j] = a*X[i][j] + Y[i][j];
}

void dbl33_set_column(dbl33 A, int i, dbl3 const x) {
  A[0][i] = x[0];
  A[1][i] = x[1];
  A[2][i] = x[2];
}

void dbl33_sub(dbl33 const A, dbl33 const B, dbl33 C) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      C[i][j] = A[i][j] - B[i][j];
}

void dbl33_sub_inplace(dbl33 A, dbl33 const B) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      A[i][j] -= B[i][j];
}

void dbl33_symmetrize(dbl33 A) {
  A[1][0] = A[0][1] = (A[1][0] + A[0][1])/2;
  A[2][0] = A[0][2] = (A[2][0] + A[0][2])/2;
  A[2][1] = A[1][2] = (A[2][1] + A[1][2])/2;
}

void dbl33_transpose(dbl33 A) {
  SWAP(A[1][0], A[0][1]);
  SWAP(A[2][0], A[0][2]);
  SWAP(A[2][1], A[1][2]);
}

void dbl33_transposed(dbl33 const A, dbl33 At) {
  memcpy((void *)At, (void *)A, sizeof(dbl)*3*3);
  SWAP(At[1][0], At[0][1]);
  SWAP(At[2][0], At[0][2]);
  SWAP(At[2][1], At[1][2]);
}

void dbl33_zero(dbl33 A) {
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      A[i][j] = 0;
}

void dbl3_dbl33_mul(dbl3 const x, dbl33 const A, dbl3 b) {
  dbl3_zero(b);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      b[j] += x[i]*A[i][j];
}

void dbl3_outer(dbl3 const u, dbl3 const v, dbl33 uv) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      uv[i][j] = u[i]*v[j];
}

void dbl4_dbl43_mul(dbl4 const b, dbl43 const A, dbl3 x) {
  dbl3_zero(x);
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 4; ++i)
      x[j] += b[i]*A[i][j];
}

void dbl43_dbl3_mul(dbl43 const A, dbl3 const b, dbl4 x) {
  dbl4_zero(x);
  for (size_t i = 0; i < 4; ++i)
    for (size_t j = 0; j < 3; ++j)
      x[i] = A[i][j]*b[j];
}

dbl dbl44_det(dbl44 const A) {
  dbl tmp1[6], tmp2[4];

  tmp1[0] =  A[1][1]*A[2][2]*A[3][3];
  tmp1[1] = -A[1][1]*A[3][2]*A[2][3];
  tmp1[2] = -A[1][2]*A[2][1]*A[3][3];
  tmp1[3] =  A[1][2]*A[3][1]*A[2][3];
  tmp1[4] =  A[1][3]*A[2][1]*A[3][2];
  tmp1[5] = -A[1][3]*A[3][1]*A[2][2];

  tmp2[0] =  A[0][0]*dblN_nsum(tmp1, 6);

  tmp1[0] =  A[1][0]*A[2][2]*A[3][3];
  tmp1[1] = -A[1][0]*A[3][2]*A[2][3];
  tmp1[2] = -A[1][2]*A[2][0]*A[3][3];
  tmp1[3] =  A[1][2]*A[3][0]*A[2][3];
  tmp1[4] =  A[1][3]*A[2][0]*A[3][2];
  tmp1[5] = -A[1][3]*A[3][0]*A[2][2];

  tmp2[1] = -A[0][1]*dblN_nsum(tmp1, 6);

  tmp1[0] =  A[1][0]*A[2][1]*A[3][3];
  tmp1[1] = -A[1][0]*A[3][1]*A[2][3];
  tmp1[2] = -A[1][1]*A[2][0]*A[3][3];
  tmp1[3] =  A[1][1]*A[3][0]*A[2][3];
  tmp1[4] =  A[1][3]*A[2][0]*A[3][1];
  tmp1[5] = -A[1][3]*A[3][0]*A[2][1];

  tmp2[2] =  A[0][2]*dblN_nsum(tmp1, 6);

  tmp1[0] =  A[1][0]*A[2][1]*A[3][2];
  tmp1[1] = -A[1][0]*A[3][1]*A[2][2];
  tmp1[2] = -A[1][1]*A[2][0]*A[3][2];
  tmp1[3] =  A[1][1]*A[3][0]*A[2][2];
  tmp1[4] =  A[1][2]*A[2][0]*A[3][1];
  tmp1[5] = -A[1][2]*A[3][0]*A[2][1];

  tmp2[3] = -A[0][3]*dblN_nsum(tmp1, 6);

  return dblN_nsum(tmp2, 4);
}

dbl dbl4_dbl44_dbl4_dot(dbl4 const x, dbl44 const A, dbl4 const y) {
  dbl z = 0;
  for (size_t i = 0; i < 4; ++i)
    for (size_t j = 0; j < 4; ++j)
      z += x[i]*A[i][j]*y[j];
  return z;
}

void dbl44_dbl4_mul(dbl44 const A, dbl4 const x, dbl4 b) {
  dbl4_zero(b);
  for (size_t i = 0; i < 4; ++i)
    for (size_t j = 0; j < 4; ++j)
      b[i] += A[i][j]*x[j];
}

void dbl44_dbl4_solve(dbl44 const A, dbl4 const b, dbl4 x) {
  dbl44 LU;
  memcpy(LU, A, 4*4*sizeof(dbl));

  dbl const n = 4;

  int perm[4] = {0, 1, 2, 3};

  // LU decomposition with partial pivoting
  for (int k = 0; k < n - 1; ++k) {
    // Find max |Ui:| for i >= k
    int argi = k;
    dbl absmax = fabs(LU[k][k]);
    for (int i = k + 1; i < n; ++i) {
      dbl tmp = fabs(LU[i][k]);
      if (tmp > absmax) {
        absmax = tmp;
        argi = i;
      }
    }

    if (argi != k) {
      // Swap rows i and k of L and U
      for (int j = 0; j < n; ++j)
        SWAP(LU[k][j], LU[argi][j]);

      // Swap entries i and k of permutation
      SWAP(perm[k], perm[argi]);
    }

    if (fabs(LU[k][k]) < 1e-15)
      continue;

    // Normalize column k of L
    for (int i = k + 1; i < n; ++i)
      LU[i][k] /= LU[k][k];

    // Subtract outer product
    for (int i = k + 1; i < n; ++i)
      for (int j = k + 1; j < n; ++j)
        LU[i][j] -= LU[i][k]*LU[k][j];
  }

  // Solve Ly = Pb
  for (int i = 0; i < n; ++i) {
    x[i] = b[perm[i]];
    for (int j = 0; j < i; ++j) {
      x[i] -= LU[i][j]*x[j];
    }
  }

  // Solve Ux = y
  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j) {
      x[i] -= LU[i][j]*x[j];
    }
    x[i] /= LU[i][i];
  }
}

void dbl44_get_column(dbl44 const A, int j, dbl4 a) {
  for (int i = 0; i < 4; ++i) {
    a[i] = A[i][j];
  }
}

void dbl44_mul(dbl44 const A, dbl44 const B, dbl44 C) {
  dbl44_zero(C);
  for (size_t i = 0; i < 4; ++i)
    for (size_t j = 0; j < 4; ++j)
      for (size_t k = 0; k < 4; ++k)
        C[i][j] += A[i][k]*B[k][j];
}

void dbl44_set_column(dbl44 A, int j, dbl4 const a) {
  for (int i = 0; i < 4; ++i) {
    A[i][j] = a[i];
  }
}

void dbl4_dbl44_mul(dbl4 const x, dbl44 const A, dbl4 b) {
  dbl4_zero(b);
  for (size_t i = 0; i < 4; ++i)
    for (size_t j = 0; j < 4; ++j)
      b[j] += x[i]*A[i][j];
}

void dbl44_zero(dbl44 A) {
  memset(A, 0x0, sizeof(dbl44));
}
