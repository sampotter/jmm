#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "field.h"

typedef struct {
  field2_s const *slow;
  dbl s, s0;
  dvec2 lp;
  dbl L;
  dvec2 xy_xy0_avg;
  dvec2 t0;
  dbl S4_th;
  dbl S4;
} S4_context;

void S4_compute(dbl th, S4_context *context);

#ifdef __cplusplus
}
#endif
