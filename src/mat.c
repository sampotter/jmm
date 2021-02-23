#include "mat.h"

#include <assert.h>
#include <string.h>

#include "macros.h"

void dbl22_add(dbl A[2][2], dbl B[2][2], dbl C[2][2]) {
  C[0][0] = A[0][0] + B[0][0];
  C[0][1] = A[0][1] + B[0][1];
  C[1][0] = A[1][0] + B[1][0];
  C[1][1] = A[1][1] + B[1][1];
}

void dbl22_dbl2_solve(dbl A[2][2], dbl b[2], dbl x[2]) {
  dbl det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
  x[0] = (A[1][1]*b[0] - A[1][0]*b[1])/det;
  x[1] = (A[0][0]*b[1] - A[0][1]*b[0])/det;
}

dbl dbl22_trace(dbl const A[2][2]) {
  return A[0][0] + A[1][1];
}

void dbl22_dbl2_mul(dbl const A[2][2], dbl const x[2], dbl b[2]) {
  b[0] = A[0][0]*x[0] + A[0][1]*x[1];
  b[1] = A[1][0]*x[0] + A[1][1]*x[1];
}

bool dbl22_isfinite(dbl const A[2][2]) {
  return isfinite(A[0][0]) && isfinite(A[1][0]) && isfinite(A[0][1])
    && isfinite(A[1][1]);
}

void dbl3_outer(dbl u[3], dbl v[3], dbl uv[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      uv[i][j] = u[i]*v[j];
    }
  }
}

void dbl33_add(dbl const A[3][3], dbl const B[3][3], dbl C[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }
}

void dbl33_mul(dbl A[3][3], dbl B[3][3], dbl C[3][3]) {
  memset((void *)C, 0x0, sizeof(dbl)*3*3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
}

void dbl33_sub(dbl A[3][3], dbl B[3][3], dbl C[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      C[i][j] = A[i][j] - B[i][j];
    }
  }
}

void dbl33_dbl3_mul(dbl const A[3][3], dbl const x[3], dbl b[3]) {
  memset((void *)b, 0x0, sizeof(dbl)*3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      b[i] += A[i][j]*x[j];
    }
  }
}

void dbl33_dbl3_mul_inplace(dbl const A[3][3], dbl x[3]) {
  dbl x_[3] = {x[0], x[1], x[2]};
  dbl33_dbl3_mul(A, x_, x);
}

void dbl33_dbl3_nmul(dbl const A[3][3], dbl const x[3], dbl b[3]) {
  for (int i = 0; i < 3; ++i)
    b[i] = dbl3_ndot(A[i], x);
}

void dbl33_transpose(dbl A[3][3]) {
  SWAP(A[1][0], A[0][1]);
  SWAP(A[2][0], A[0][2]);
  SWAP(A[2][1], A[1][2]);
}

void dbl33_transposed(dbl const A[3][3], dbl At[3][3]) {
  memcpy((void *)At, (void *)A, sizeof(dbl)*3*3);
  SWAP(At[1][0], At[0][1]);
  SWAP(At[2][0], At[0][2]);
  SWAP(At[2][1], At[1][2]);
}

void dbl33_dbl_div(dbl A[3][3], dbl a, dbl B[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      B[i][j] = A[i][j]/a;
    }
  }
}

dbl dbl33_det(dbl const A[3][3]) {
  return A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
       - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
       + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

void dbl33_dbl3_solve(dbl const A[3][3], dbl const b[3], dbl x[3]) {
  // TODO: bad---solving using Cramer's rule... need to implement
  // pivoted LU instead
  dbl det = dbl33_det(A);
  dbl tmp[3], A_[3][3];
  memcpy(A_, A, sizeof(dbl[3][3]));
  for (int j = 0; j < 3; ++j) {
    dbl33_get_column(A_, j, tmp);
    dbl33_set_column(A_, j, b);
    x[j] = dbl33_det(A_)/det;
    dbl33_set_column(A_, j, tmp);
  }
}

void dbl33_get_column(dbl const A[3][3], int i, dbl x[3]) {
  x[0] = A[0][i];
  x[1] = A[1][i];
  x[2] = A[2][i];
}

void dbl33_set_column(dbl A[3][3], int i, dbl const x[3]) {
  A[0][i] = x[0];
  A[1][i] = x[1];
  A[2][i] = x[2];
}

void dbl44_dbl4_solve(dbl const A[4][4], dbl const b[4], dbl x[4]) {
  dbl LU[4][4];
  memcpy(LU, A, 4*4*sizeof(dbl));

  dbl const n = 4;

  // LU decomposition
  for (int k = 0; k < n - 1; ++k) {
    for (int i = k + 1; i < n; ++i) {
      LU[i][k] /= LU[k][k];
    }
    for (int i = k + 1; i < n; ++i) {
      for (int j = k + 1; j < n; ++j) {
        LU[i][j] -= LU[i][k]*LU[k][j];
      }
    }
  }

  // forward sub
  for (int i = 0; i < n; ++i) {
    x[i] = b[i];
    for (int j = 0; j < i; ++j) {
      x[i] -= LU[i][j]*x[j];
    }
  }

  // backward sub
  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j) {
      x[i] -= LU[i][j]*x[j];
    }
    x[i] /= LU[i][i];
  }
}

dbl dbl44_det(dbl const A[4][4]) {
  return A[0][0]*(
    A[1][1]*(A[2][2]*A[3][3] - A[3][2]*A[2][3])
  - A[1][2]*(A[2][1]*A[3][3] - A[3][1]*A[2][3])
  + A[1][3]*(A[2][1]*A[3][2] - A[3][1]*A[2][2])
  ) - A[0][1]*(
    A[1][0]*(A[2][2]*A[3][3] - A[3][2]*A[2][3])
  - A[1][2]*(A[2][0]*A[3][3] - A[3][0]*A[2][3])
  + A[1][3]*(A[2][0]*A[3][2] - A[3][0]*A[2][2])
  ) + A[0][2]*(
    A[1][0]*(A[2][1]*A[3][3] - A[3][1]*A[2][3])
  - A[1][1]*(A[2][0]*A[3][3] - A[3][0]*A[2][3])
  + A[1][3]*(A[2][0]*A[3][1] - A[3][0]*A[2][1])
  ) - A[0][3]*(
    A[1][0]*(A[2][1]*A[3][2] - A[3][1]*A[2][2])
  - A[1][1]*(A[2][0]*A[3][2] - A[3][0]*A[2][2])
  + A[1][2]*(A[2][0]*A[3][1] - A[3][0]*A[2][1])
  );
}

void dbl44_get_col(dbl const A[4][4], int j, dbl a[4]) {
  for (int i = 0; i < 4; ++i) {
    a[i] = A[i][j];
  }
}

void dbl44_set_col(dbl A[4][4], int j, dbl const a[4]) {
  for (int i = 0; i < 4; ++i) {
    A[i][j] = a[i];
  }
}

dmat22 dmat22_add(dmat22 A, dmat22 B) {
  return (dmat22) {
    .rows = {
      {A.rows[0].x + B.rows[0].x, A.rows[0].y + B.rows[0].y},
      {A.rows[1].x + B.rows[1].x, A.rows[1].y + B.rows[1].y},
    }
  };
}

dmat22 dmat22_sub(dmat22 A, dmat22 B) {
  return (dmat22) {
    .rows = {
      {A.rows[0].x - B.rows[0].x, A.rows[0].y - B.rows[0].y},
      {A.rows[1].x - B.rows[1].x, A.rows[1].y - B.rows[1].y},
    }
  };
}

dmat22 dmat22_mul(dmat22 A, dmat22 B) {
  dmat22 C;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      C.data[i][j] = 0;
      for (int k = 0; k < 2; ++k) {
        C.data[i][j] += A.data[i][k]*B.data[k][j];
      }
    }
  }
  return C;
}

dmat22 dmat22_dbl_mul(dmat22 A, dbl a) {
  return (dmat22) {
    .rows = {
      {a*A.rows[0].x, a*A.rows[0].y},
      {a*A.rows[1].x, a*A.rows[1].y}
    }
  };
}

dmat22 dmat22_dbl_div(dmat22 A, dbl a) {
  return (dmat22) {
    .rows = {
      {A.rows[0].x/a, A.rows[0].y/a},
      {A.rows[1].x/a, A.rows[1].y/a}
    }
  };
}

dvec2 dmat22_dvec2_mul(dmat22 A, dvec2 x) {
  return (dvec2) {
    A.rows[0].x*x.x + A.rows[0].y*x.y,
    A.rows[1].x*x.x + A.rows[1].y*x.y
  };
}

dvec2 dmat22_dvec2_solve(dmat22 A, dvec2 b) {
  dbl det = A.data[0][0]*A.data[1][1] - A.data[0][1]*A.data[1][0];
  return (dvec2) {
    .x = (A.data[1][1]*b.x - A.data[0][1]*b.y)/det,
    .y = (A.data[0][0]*b.y - A.data[1][0]*b.x)/det
  };
}

dmat22 dvec2_outer(dvec2 u, dvec2 v) {
  return (dmat22) {
    .rows = {
      {u.x*v.x, u.x*v.y},
      {u.y*v.x, u.y*v.y}
    }
  };
}

void dmat22_invert(dmat22 *A) {
  dbl det = A->data[0][0]*A->data[1][1] - A->data[1][0]*A->data[0][1];
  *A = (dmat22) {
    .rows = {
      { A->rows[1].y/det, -A->rows[0].y/det},
      {-A->rows[1].x/det,  A->rows[0].x/det}
    }
  };
}

dbl dmat22_trace(dmat22 const *A) {
  return A->data[0][0] + A->data[1][1];
}

dbl dmat22_det(dmat22 const *A) {
  return A->data[0][0]*A->data[1][1] - A->data[0][1]*A->data[1][0];
}

void dmat22_eigvals(dmat22 const *A, dbl *lam1, dbl *lam2) {
  dbl tr = dmat22_trace(A);
  dbl disc = tr*tr - 4*dmat22_det(A);
  assert(disc >= 0);
  dbl tmp = sqrt(disc);
  *lam1 = (tr + tmp)/2;
  *lam2 = (tr - tmp)/2;
}

void dmat22_transpose(dmat22 *A) {
  dbl tmp = A->data[1][0];
  A->data[1][0] = A->data[0][1];
  A->data[0][1] = tmp;
}

#define a(i, j) A->rows[i].data[j]

dbl dmat33_det(dmat33 const *A) {
  return a(0, 0)*(a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)) -
    a(0, 1)*(a(1, 0)*a(2, 2) - a(1, 2)*a(2, 0)) +
    a(0, 2)*(a(1, 0)*a(2, 1) - a(1, 1)*a(2, 0));
}

#undef a

dmat33 dmat33_eye(void) {
  dmat33 I;
  memset(&I, 0x0, sizeof(dmat33));
  I.rows[0].data[0] = 1;
  I.rows[1].data[1] = 1;
  I.rows[2].data[2] = 1;
  return I;
}

dvec3 dmat33_getcol(dmat33 const *A, int j) {
  dvec3 a;
  for (int i = 0; i < 3; ++i) {
    a.data[i] = A->rows[j].data[i];
  }
  return a;
}

void dmat33_setcol(dmat33 *A, dvec3 a, int j) {
  for (int i = 0; i < 3; ++i) {
    A->rows[j].data[i] = a.data[i];
  }
}

dvec3 dmat33_dvec3_mul(dmat33 A, dvec3 x) {
  return (dvec3) {
    .data = {
      dvec3_dot(A.rows[0], x),
      dvec3_dot(A.rows[1], x),
      dvec3_dot(A.rows[2], x)
    }
  };
}

dmat33 dmat33_dbl_div(dmat33 A, dbl a) {
  return (dmat33) {
    .rows = {
      dvec3_dbl_div(A.rows[0], a),
      dvec3_dbl_div(A.rows[1], a),
      dvec3_dbl_div(A.rows[2], a)
    }
  };
}

dvec3 dmat33_dvec3_solve(dmat33 A, dvec3 b) {
  dbl det = dmat33_det(&A);
  dvec3 x, tmp;
  for (int j = 0; j < 3; ++j) {
    tmp = dmat33_getcol(&A, j);
    dmat33_setcol(&A, b, j);
    x.data[j] = dmat33_det(&A)/det;
    dmat33_setcol(&A, tmp, j);
  }
  return x;
}

dmat33 dmat33_mul(dmat33 A, dmat33 B) {
  dmat33 C;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      C.rows[i].data[j] = dvec3_dot(A.rows[i], dmat33_getcol(&B, j));
    }
  }
  return C;
}

dmat33 dmat33_sub(dmat33 A, dmat33 B) {
  return (dmat33) {
    .rows = {
      dvec3_sub(A.rows[0], B.rows[0]),
      dvec3_sub(A.rows[1], B.rows[1]),
      dvec3_sub(A.rows[2], B.rows[2])
    }
  };
}

dmat33 dvec3_outer(dvec3 u, dvec3 v) {
  dmat33 uv;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      uv.rows[i].data[j] = u.data[i]*v.data[j];
    }
  }
  return uv;
}

void dmat33_transpose(dmat33 *A) {
  SWAP(A->rows[1].data[0], A->rows[0].data[1]);
  SWAP(A->rows[2].data[0], A->rows[0].data[2]);
  SWAP(A->rows[2].data[1], A->rows[1].data[2]);
}

dvec4 dmat44_dvec4_mul(dmat44 const A, dvec4 const x) {
  dvec4 y;
  for (int i = 0; i < 4; ++i) {
    y.data[i] = dvec4_dot(A.rows[i], x);
  }
  return y;
}

dvec4 dvec4_dmat44_mul(dvec4 const x, dmat44 const A) {
  dvec4 y;
  for (int j = 0; j < 4; ++j) {
    y.data[j] = dvec4_dot(x, dmat44_col(A, j));
  }
  return y;
}

dmat44 dmat44_dmat44_mul(dmat44 const A, dmat44 const B) {
  dmat44 C;
  for (int i = 0; i < 4; ++i) {
    C.rows[i] = dvec4_dmat44_mul(A.rows[i], B);
  }
  return C;
}

dvec4 dmat44_col(dmat44 const A, int j) {
  dvec4 a;
  for (int i = 0; i < 4; ++i) {
    a.data[i] = A.rows[i].data[j];
  }
  return a;
}
