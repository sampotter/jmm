#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

void dbl22_add(dbl A[2][2], dbl B[2][2], dbl C[2][2]);
void dbl22_dbl2_solve(dbl A[2][2], dbl b[2], dbl x[2]);
dbl dbl22_trace(dbl const A[2][2]);
void dbl22_dbl2_mul(dbl const A[2][2], dbl const x[2], dbl b[2]);
bool dbl22_isfinite(dbl const A[2][2]);

void dbl3_outer(dbl u[3], dbl v[3], dbl uv[3][3]);
void dbl33_add(dbl const A[3][3], dbl const B[3][3], dbl C[3][3]);
void dbl33_mul(dbl A[3][3], dbl B[3][3], dbl C[3][3]);
void dbl33_sub(dbl A[3][3], dbl B[3][3], dbl C[3][3]);
void dbl33_dbl3_mul(dbl const A[3][3], dbl const x[3], dbl b[3]);
void dbl33_transpose(dbl A[3][3]);
void dbl33_transposed(dbl const A[3][3], dbl At[3][3]);
void dbl33_dbl_div(dbl A[3][3], dbl a, dbl B[3][3]);
dbl dbl33_det(dbl const A[3][3]);

void dbl44_dbl4_solve(dbl const A[4][4], dbl const b[4], dbl x[4]);
dbl dbl44_det(dbl const A[4][4]);
void dbl44_get_col(dbl const A[4][4], int j, dbl a[4]);
void dbl44_set_col(dbl A[4][4], int j, dbl const a[4]);

typedef struct {
  union {
    dbl data[2][2];
    dvec2 rows[2];
  };
} dmat22;

dmat22 dmat22_add(dmat22 A, dmat22 B);
dmat22 dmat22_sub(dmat22 A, dmat22 B);
dmat22 dmat22_mul(dmat22 A, dmat22 B);
dmat22 dmat22_dbl_mul(dmat22 A, dbl a);
dmat22 dmat22_dbl_div(dmat22 A, dbl a);
dvec2 dmat22_dvec2_mul(dmat22 A, dvec2 x);
dvec2 dmat22_dvec2_solve(dmat22 A, dvec2 b);
dmat22 dvec2_outer(dvec2 u, dvec2 v);
void dmat22_invert(dmat22 *A);
dbl dmat22_trace(dmat22 const *A);
dbl dmat22_det(dmat22 const *A);
void dmat22_eigvals(dmat22 const *A, dbl *lam1, dbl *lam2);
void dmat22_transpose(dmat22 *A);

typedef struct {
  dvec3 rows[3];
} dmat33;

dbl dmat33_det(dmat33 const *A);
dmat33 dmat33_eye();
dvec3 dmat33_getcol(dmat33 const *A, int j);
void dmat33_setcol(dmat33 *A, dvec3 a, int j);
dvec3 dmat33_dvec3_mul(dmat33 A, dvec3 x);
dmat33 dmat33_dbl_div(dmat33 A, dbl a);
dvec3 dmat33_dvec3_solve(dmat33 A, dvec3 b);
dmat33 dmat33_mul(dmat33 A, dmat33 B);
dmat33 dmat33_sub(dmat33 A, dmat33 B);
dmat33 dvec3_outer(dvec3 u, dvec3 v);
void dmat33_transpose(dmat33 *A);

typedef struct {
  union {
    dbl data[4][4];
    dvec4 rows[4];
  };
} dmat44;

dvec4 dmat44_dvec4_mul(dmat44 const A, dvec4 const x);
dvec4 dvec4_dmat44_mul(dvec4 const x, dmat44 const A);
dmat44 dmat44_dmat44_mul(dmat44 const A, dmat44 const B);
dvec4 dmat44_col(dmat44 const A, int j);

#ifdef __cplusplus
}
#endif
