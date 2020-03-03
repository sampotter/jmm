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

dvec2 dvec2_ccomb(dvec2 v0, dvec2 v1, dbl t);
dbl dvec2_dist(dvec2 v0, dvec2 v1);
dbl dvec2_norm(dvec2 v);
dbl dvec2_dot(dvec2 u, dvec2 v);
dvec2 dvec2_sub(dvec2 u, dvec2 v);
dvec2 dvec2_dbl_div(dvec2 v, dbl a);
dvec2 dvec2_floor(dvec2 v);

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
dvec4 dvec4_m(dbl x);
dvec4 dvec4_dm(dbl x);
dvec4 dvec4_e1();
dvec4 dvec4_one();
dvec4 dvec4_iota();

typedef struct {
  int i;
  int j;
} ivec2;

ivec2 dvec2_to_ivec2(dvec2 v);

#ifdef __cplusplus
}
#endif
