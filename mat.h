#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

typedef struct {
  union {
    dbl data[2][2];
    dvec2 rows[2];
  };
} dmat22;

dmat22 dmat22_add(dmat22 A, dmat22 B);
dmat22 dmat22_sub(dmat22 A, dmat22 B);
dmat22 dmat22_dbl_div(dmat22 A, dbl a);
dvec2 dmat22_dvec2_mul(dmat22 A, dvec2 x);
dvec2 dmat22_dvec2_solve(dmat22 A, dvec2 b);
dmat22 dvec2_outer(dvec2 u, dvec2 v);

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
