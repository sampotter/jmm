#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

#include <immintrin.h>
#include <math.h>

#define DBL2(x, y) (dbl2) {x, y}
#define DBL3(x, y, z) (dbl3) {x, y, z}
#define DBL4(x, y, z, w) (dbl4) {x, y, z, w}

#define DBL4_FROM_PTR(p) (dbl4) {p[0], p[1], p[2], p[3]}

bool dbl2_bary(dbl2 const u);
bool dbl2_isfinite(dbl2 const u);
bool dbl2_all_nan(dbl2 const u);
dbl dbl2_dist(dbl2 const u, dbl2 const v);
dbl dbl2_dot(dbl2 const u, dbl2 const v);
dbl dbl2_maxdist(dbl2 const u, dbl2 const v);
dbl dbl2_maxnorm(dbl2 const u);
dbl dbl2_norm(dbl2 const u);
dbl dbl2_normsq(dbl2 const u);
dbl dbl2_prod(dbl2 const u);
dbl dbl2_sum(dbl2 const u);
void dbl2_add(dbl2 const u, dbl2 const v, dbl2 w);
void dbl2_avg(dbl2 const u, dbl2 const v, dbl2 w);
void dbl2_copy(dbl2 const u, dbl2 v);
void dbl2_cproj(dbl2 const u, dbl2 const v, dbl2 w);
void dbl2_dbl_div(dbl2 const u, dbl a, dbl2 v);
void dbl2_dbl_div_inplace(dbl2 u, dbl a);
void dbl2_dbl_mul(dbl2 const u, dbl a, dbl2 v);
void dbl2_floor(dbl2 const u, dbl2 v);
void dbl2_negate(dbl2 u);
void dbl2_normalize(dbl2 u);
void dbl2_normalized(dbl2 const u, dbl2 v);
void dbl2_saxpy(dbl a, dbl2 const x, dbl2 const y, dbl2 z);
void dbl2_sub(dbl2 const u, dbl2 const v, dbl2 w);
void dbl2_sub_inplace(dbl2 u, dbl2 const v);
void dbl2_zero(dbl2 u);
void dbl2_lincomb(dbl a, dbl2 const u, dbl b, dbl2 const v, dbl2 w);

bool dbl3_equal(dbl3 const x, dbl3 const y);
bool dbl3_is_normalized(dbl3 const u);
bool dbl3_is_zero(dbl3 const u);
bool dbl3_isfinite(dbl3 const x);
bool dbl3_all_nan(dbl3 const x);
bool dbl3_nonneg(dbl3 const x);
bool dbl3_valid_bary_coord(dbl3 const b);
dbl dbl3_dist(dbl3 const u, dbl3 const v);
dbl dbl3_dist_sq(dbl3 const u, dbl3 const v);
dbl dbl3_dot(dbl3 const u, dbl3 const v);
dbl dbl3_maxdist(dbl3 const u, dbl3 const v);
dbl dbl3_maxnorm(dbl3 const u);
dbl dbl3_minimum(dbl3 const u);
dbl dbl3_ndot(dbl3 const u, dbl3 const v);
dbl dbl3_norm(dbl3 const u);
dbl dbl3_normalize(dbl3 u);
dbl dbl3_normsq(dbl3 const u);
dbl dbl3_nsum(dbl3 const u);
dbl dbl3_sum(dbl3 const u);
dbl dbl3_wnormsq(dbl33 const A, dbl3 const x);
size_t dbl3_argmax(dbl3 const u);
void dbl3_abs(dbl3 const u, dbl3 v);
void dbl3_argsort(dbl3 const u, size_t perm[3]);
void dbl3_add(dbl3 const u, dbl3 const v, dbl3 w);
void dbl3_add_inplace(dbl3 u, dbl3 const v);
void dbl3_avg(dbl3 const u, dbl3 const v, dbl3 w);
void dbl3_cc(dbl3 const u0, dbl3 const u1, dbl t0, dbl3 ut);
void dbl3_copy(dbl3 const u, dbl3 v);
void dbl3_cross(dbl3 const u, dbl3 const v, dbl3 w);
void dbl3_dbl_div(dbl3 const u, dbl a, dbl3 v);
void dbl3_dbl_div_inplace(dbl3 u, dbl a);
void dbl3_dbl_mul(dbl3 const u, dbl a, dbl3 v);
void dbl3_dbl_mul_inplace(dbl3 u, dbl a);
void dbl3_get_rand_ortho(dbl3 const x, dbl3 y);
void dbl3_inf(dbl3 u);
void dbl3_max(dbl3 const u, dbl3 const v, dbl3 w);
void dbl3_min(dbl3 const u, dbl3 const v, dbl3 w);
void dbl3_nan(dbl3 u);
void dbl3_negate(dbl3 u);
void dbl3_neginf(dbl3 u);
void dbl3_normalize1(dbl3 x);
void dbl3_normalized(dbl3 const u, dbl3 v);
void dbl3_one(dbl3 u);
void dbl3_saxpy(dbl a, dbl3 const x, dbl3 const y, dbl3 z);
void dbl3_saxpy_inplace(dbl a, dbl3 const x, dbl3 y);
void dbl3_sort(dbl3 u);
void dbl3_sub(dbl3 const u, dbl3 const v, dbl3 w);
void dbl3_sub_inplace(dbl3 u, dbl3 const v);
void dbl3_zero(dbl3 u);

bool dbl4_is_rgba(dbl4 const u);
bool dbl4_nonneg(dbl4 const u);
bool dbl4_valid_bary_coord(dbl4 const b);
dbl dbl4_dist(dbl4 const u, dbl4 const v);
dbl dbl4_dot(dbl4 const u, dbl4 const v);
dbl dbl4_norm(dbl4 const u);
dbl dbl4_nsum(dbl4 const u);
dbl dbl4_sum(dbl4 const u);
void dbl4_add(dbl4 const u, dbl4 const v, dbl4 w);
void dbl4_copy(dbl4 const u, dbl4 v);
void dbl4_d2m(dbl x, dbl4 d2m);
void dbl4_dbl_div(dbl4 const u, dbl a, dbl4 v);
void dbl4_dbl_div_inplace(dbl4 u, dbl a);
void dbl4_dm(dbl x, dbl4 dm);
void dbl4_e1(dbl4 e1);
void dbl4_e(dbl4 e, size_t i);
void dbl4_iota(dbl4 u);
void dbl4_m(dbl x, dbl4 m);
void dbl4_normalize1(dbl4 u);
void dbl4_one(dbl4 u);
void dbl4_saxpy(dbl a, dbl4 const x, dbl4 const y, dbl4 z);
void dbl4_sub(dbl4 const u, dbl4 const v, dbl4 w);
void dbl4_zero(dbl4 u);

dbl dblN_mean(dbl const *x, size_t n);
dbl dblN_ndot(dbl const *x, dbl const *y, size_t n);
dbl dblN_nsum(dbl const *x, size_t n);
void dblN_minmax(dbl const *x, size_t n, dbl *min, dbl *max);
dbl dblN_binmedian(size_t n, dbl const *x, size_t num_bins);
dbl dblN_median(size_t n, dbl const *x);

void int2_add(int2 const p, int2 const q, int2 r);
void int2_copy(int2 const p, int2 q);

bool int3_equal(int3 const p, int3 const q);
int int3_prod(int3 const p);
void int3_add(int3 const p, int3 const q, int3 r);
void int3_dbl_mul(int3 const p, dbl a, dbl3 x);
void int3_int_div(int3 const p, int q, int3 r);

bool uint3_equal(uint3 const i, uint3 const j);

#ifdef __cplusplus
}
#endif
