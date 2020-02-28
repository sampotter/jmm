#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

typedef struct {
  union {
    dbl data[16];
    dvec4 rows[4];
  };
} dmat44;

dvec4 dmat44_dvec4_mul(dmat44 const A, dvec4 const x);
dmat44 dmat44_dmat44_mul(dmat44 const A, dmat44 const B);
dvec4 dmat44_col(dmat44 const A, int j);

#ifdef __cplusplus
}
#endif
