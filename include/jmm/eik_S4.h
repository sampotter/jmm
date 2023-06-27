#pragma once

#include "field.h"

typedef struct {
  field2_s const *slow;
  dbl s, s0;
  dbl2 lp;
  dbl L;
  dbl2 xy_xy0_avg;
  dbl2 t0;
  dbl S4_th;
  dbl S4;
} S4_context;

void S4_compute(dbl th, S4_context *context);
