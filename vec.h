#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "vec.h"

#include <math.h>

typedef struct {
  dbl x;
  dbl y;
} dvec2;

dvec2 dvec2_zero();
dvec2 dvec2_ccomb(dvec2 v0, dvec2 v1, dbl t);
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

typedef struct {
  dbl x;
  dbl y;
  dbl z;
} dvec3;

dvec3 dvec3_dbl_div(dvec3 u, dbl a);
dbl dvec3_dist(dvec3 u, dvec3 v);
dbl dvec3_dot(dvec3 u, dvec3 v);
dbl dvec3_maxdist(dvec3 u, dvec3 v);
dbl dvec3_maxnorm(dvec3 u);
dvec3 dvec3_nan();
dbl dvec3_norm(dvec3 u);
dvec3 dvec3_normalized(dvec3 u);
dvec3 dvec3_saxpy(dbl a, dvec3 x, dvec3 y);
dvec3 dvec3_sub(dvec3 u, dvec3 v);

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
  int i;
  int j;
  int k;
} ivec3;

int ivec3_prod(ivec3 p);
dvec3 ivec3_dbl_mul(ivec3 p, dbl a);
bool ivec3_equal(ivec3 p, ivec3 q);
ivec3 ivec3_int_div(ivec3 p, int q);

ivec2 dvec2_to_ivec2(dvec2 v);

#ifdef __cplusplus
}
#endif
