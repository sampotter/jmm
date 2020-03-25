#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "cubic.h"
#include "mat.h"
#include "vec.h"

dvec2 get_lp(dvec2 xy, dvec2 xy0, dvec2 xy1, dbl eta);

typedef struct {
  cubic_s cubic;
  dvec2 xy, xy0, xy1;
} F3_context;

dbl F3(dbl eta, void *context);
dbl dF3_deta(dbl eta, void *context);

typedef struct {
  cubic_s T, Tx, Ty;
  dvec2 xy, xy0, xy1;
} F4_context;

dbl F4(dbl eta, dbl th, void *context);
dbl dF4_deta(dbl eta, dbl th, void *context);
dbl dF4_dth(dbl eta, dbl th, void *context);
dvec2 grad_F4(dbl eta, dbl th, void *context);
dmat22 hess_F4(dbl eta, dbl th, void *context);

#ifdef __cplusplus
}
#endif
