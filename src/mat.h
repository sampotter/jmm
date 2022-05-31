#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

#define DBL22(a00, a01, a10, a11) (dbl22) {{a00, a01}, {a10, a11}}

bool dbl22_isfinite(dbl22 const A);
bool dbl22_all_nan(dbl22 const A);
dbl dbl22_det(dbl22 const X);
dbl dbl22_trace(dbl22 const A);
void dbl22_add(dbl22 const A, dbl22 const B, dbl22 C);
void dbl22_add_inplace(dbl22 A, dbl22 const B);
void dbl22_copy(dbl22 const A, dbl22 B);
void dbl22_dbl2_mul(dbl22 const A, dbl2 const x, dbl2 b);
void dbl22_dbl2_solve(dbl22 const A, dbl2 const b, dbl2 x);
void dbl22_dbl_div(dbl22 const A, dbl a, dbl22 B);
void dbl22_dbl_div_inplace(dbl22 A, dbl a);
void dbl22_dbl_mul(dbl22 const A, dbl a, dbl22 B);
void dbl22_dbl_mul_inplace(dbl22 A, dbl a);
void dbl22_eigvals(dbl22 const A, dbl2 lam);
void dbl22_eye(dbl22 eye);
void dbl22_inv(dbl22 const X, dbl22 Y);
void dbl22_invert(dbl22 X);
void dbl22_mul(dbl22 const A, dbl22 const B, dbl22 C);
void dbl22_perturb(dbl22 A, dbl eps);
void dbl22_saxpy(dbl a, dbl22 const X, dbl22 const Y, dbl22 Z);
void dbl22_sub(dbl22 const A, dbl22 const B, dbl22 C);
void dbl22_transpose(dbl22 A);
void dbl22_zero(dbl22 A);
void dbl2_outer(dbl2 const u, dbl2 const v, dbl22 uv);

bool dbl33_isfinite(dbl33 const A);
bool dbl33_isnan(dbl33 const A);
dbl dbl33_det(dbl33 const A);
dbl dbl3_dbl33_dbl3_dot(dbl3 const x, dbl33 const A, dbl3 const y);
void dbl33_add(dbl33 const A, dbl33 const B, dbl33 C);
void dbl33_add_inplace(dbl33 A, dbl33 const B);
void dbl33_copy(dbl33 const in, dbl33 out);
void dbl33_dbl3_mul(dbl33 const A, dbl3 const x, dbl3 b);
void dbl33_dbl3_mul_inplace(dbl33 const A, dbl3 x);
void dbl33_dbl3_nmul(dbl33 const A, dbl3 const x, dbl3 b);
void dbl33_dbl3_solve(dbl33 const A, dbl3 const b, dbl3 x);
void dbl33_dbl_div(dbl33 const A, dbl a, dbl33 B);
void dbl33_dbl_div_inplace(dbl33 A, dbl a);
void dbl33_eigvals_sym(dbl33 const A, dbl3 lam);
void dbl33_eye(dbl33 A);
void dbl33_get_column(dbl33 const A, int i, dbl3 a);
void dbl33_invert(dbl33 A);
void dbl33_make_axis_angle_rotation_matrix(dbl3 axis, dbl angle, dbl33 rot);
void dbl33_mul(dbl33 const A, dbl33 const B, dbl33 C);
void dbl33_nan(dbl33 A);
void dbl33_saxpy(dbl a, dbl33 const X, dbl33 const Y, dbl33 Z);
void dbl33_set_column(dbl33 A, int i, dbl3 const a);
void dbl33_sub(dbl33 const A, dbl33 const B, dbl33 C);
void dbl33_sub_inplace(dbl33 A, dbl33 const B);
void dbl33_symmetrize(dbl33 A);
void dbl33_transpose(dbl33 A);
void dbl33_transposed(dbl33 const A, dbl33 At);
void dbl33_zero(dbl33 A);
void dbl3_dbl33_mul(dbl3 const x, dbl33 const A, dbl3 b);
void dbl3_outer(dbl3 const u, dbl3 const v, dbl33 uv);

void dbl4_dbl43_mul(dbl4 const b, dbl43 const A, dbl3 x);
void dbl43_dbl3_mul(dbl43 const A, dbl3 const b, dbl4 x);

dbl dbl44_det(dbl44 const A);
dbl dbl4_dbl44_dbl4_dot(dbl4 const x, dbl44 const A, dbl4 const y);
void dbl44_dbl4_mul(dbl44 const A, dbl4 const x, dbl4 b);
void dbl44_dbl4_solve(dbl44 const A, dbl4 const b, dbl4 x);
void dbl44_copy(dbl44 const A, dbl44 B);
void dbl44_get_column(dbl44 const A, int j, dbl4 a);
void dbl44_mul(dbl44 const A, dbl44 const B, dbl44 C);
void dbl44_set_column(dbl44 A, int j, dbl4 const a);
void dbl4_dbl44_mul(dbl4 const x, dbl44 const A, dbl4 b);
void dbl44_zero(dbl44 A);
void dbl44_invert(dbl44 A);
void dbl44_transpose(dbl44 A);

#ifdef __cplusplus
}
#endif
