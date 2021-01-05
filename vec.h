#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

#include <immintrin.h>
#include <math.h>

void dbl2_add(dbl u[2], dbl v[2], dbl w[2]);
void dbl2_sub(dbl const *u, dbl const *v, dbl *w);
dbl dbl2_dot(dbl const *u, dbl const *v);
void dbl2_negate(dbl u[2]);
dbl dbl2_sum(dbl const u[2]);
void dbl2_saxpy(dbl a, dbl const x[2], dbl const y[2], dbl z[2]);
dbl dbl2_dist(dbl const u[2], dbl const v[2]);
dbl dbl2_norm(dbl const u[2]);
dbl dbl2_norm_sq(dbl const u[2]);
dbl dbl2_maxdist(dbl const u[2], dbl const v[2]);
dbl dbl2_maxnorm(dbl const u[2]);
bool dbl2_isfinite(dbl const u[2]);

void dbl3_add(dbl const u[3], dbl const v[3], dbl w[3]);
void dbl3_sub(dbl const *u, dbl const *v, dbl *w);
dbl dbl3_dot(dbl const *u, dbl const *v);
dbl dbl3_dist(dbl const u[3], dbl const v[3]);
dbl dbl3_norm(dbl const u[3]);
void dbl3_dbl_div(dbl u[3], dbl a, dbl v[3]);
dbl dbl3_normalize(dbl u[3]);
void dbl3_saxpy(dbl a, dbl const x[3], dbl const y[3], dbl z[3]);
bool dbl3_is_zero(dbl const u[3]);
dbl dbl3_maxnorm(dbl const u[3]);
void dbl3_cross(dbl const u[3], dbl const v[3], dbl w[3]);
dbl dbl3_sum(dbl const u[3]);

bool dbl4_nonneg(dbl const u[4]);
dbl dbl4_sum(dbl const u[4]);

typedef struct {
  dbl x;
  dbl y;
} dvec2;

void int3_add(int const p[3], int const q[3], int r[3]);

dvec2 dvec2_zero();
dvec2 dvec2_ccomb(dvec2 v0, dvec2 v1, dbl t);
dvec2 dvec2_from_ptr(dbl const *u);
dbl dvec2_maxdist(dvec2 u, dvec2 v);
dbl dvec2_dist(dvec2 v0, dvec2 v1);
dbl dvec2_maxnorm(dvec2 v);
dbl dvec2_norm(dvec2 v);
dbl dvec2_norm_sq(dvec2 v);
dbl dvec2_dot(dvec2 u, dvec2 v);
dvec2 dvec2_add(dvec2 u, dvec2 v);
dvec2 dvec2_sub(dvec2 u, dvec2 v);
dvec2 dvec2_saxpy(dbl a, dvec2 x, dvec2 y);
dvec2 dvec2_dbl_div(dvec2 v, dbl a);
dvec2 dvec2_dbl_mul(dvec2 v, dbl a);
dvec2 dvec2_floor(dvec2 v);
void dvec2_negate(dvec2 *v);
void dvec2_normalize(dvec2 *v);
dvec2 dvec2_cproj(dvec2 u, dvec2 v);
dvec2 dvec2_avg(dvec2 u, dvec2 v);
dbl dvec2_sum(dvec2 u);

typedef struct {
  union __attribute__((aligned(16))) {
    dbl data[4];
    __m256d packed;
  };
} dvec3;

dvec3 dvec3_add(dvec3 u, dvec3 v);
dvec3 dvec3_dbl_div(dvec3 u, dbl a);
dvec3 dvec3_dbl_mul(dvec3 u, dbl a);
dbl dvec3_dist(dvec3 u, dvec3 v);
dbl dvec3_dist_sq(dvec3 u, dvec3 v);
dbl dvec3_dot(dvec3 u, dvec3 v);
dvec3 dvec3_infinity();
dbl dvec3_maxdist(dvec3 u, dvec3 v);
dbl dvec3_maxnorm(dvec3 u);
dvec3 dvec3_nan();
dbl dvec3_norm(dvec3 u);
dbl dvec3_norm_sq(dvec3 u);
dvec3 dvec3_normalized(dvec3 u);
dvec3 dvec3_one();
dvec3 dvec3_saxpy(dbl a, dvec3 x, dvec3 y);
dvec3 dvec3_sub(dvec3 u, dvec3 v);
dvec3 dvec3_zero();
int dvec3_argmax(dvec3 u);
void dvec3_normalize(dvec3 *u);

typedef struct {
  union {
    dbl data[4];
    struct {
      dbl x;
      dbl y;
      dbl z;
      dbl w;
    } xyzw;
  };
} dvec4;

dbl dvec4_dot(dvec4 v0, dvec4 v1);
dbl dvec4_sum(dvec4 v);
dvec4 dvec4_add(dvec4 u, dvec4 v);
dvec4 dvec4_dbl_div(dvec4 u, dbl a);
dvec4 dvec4_m(dbl x);
dvec4 dvec4_dm(dbl x);
dvec4 dvec4_d2m(dbl x);
dvec4 dvec4_e1();
dvec4 dvec4_one();
dvec4 dvec4_iota();

typedef struct {
  int i;
  int j;
} ivec2;

ivec2 ivec2_add(ivec2 p, ivec2 q);

typedef struct {
  union __attribute__((aligned(16))) {
    int data[4];
    __m128i packed;
  };
} ivec3;

ivec3 ivec3_from_int3(int const p[3]);
ivec3 ivec3_add(ivec3 p, ivec3 q);
int ivec3_prod(ivec3 p);
dvec3 ivec3_dbl_mul(ivec3 p, dbl a);
bool ivec3_equal(ivec3 p, ivec3 q);
ivec3 ivec3_int_div(ivec3 p, int q);

ivec2 dvec2_to_ivec2(dvec2 v);
ivec3 dvec3_to_ivec3(dvec3 x);
dvec3 ivec3_to_dvec3(ivec3 p);

#ifdef __cplusplus
}
#endif
