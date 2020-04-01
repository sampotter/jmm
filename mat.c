#include "mat.h"

#include <assert.h>

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
