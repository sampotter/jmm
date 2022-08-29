#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "cubic.h"
#include "field.h"
#include "vec.h"

typedef struct {
  // Inputs:
  cubic_s T_cubic;
  dbl2 xy, xy0, xy1;
  field2_s const *slow;

  // Outputs:
  dbl F3;
  dbl F3_eta;
} F3_context;

void F3_compute(dbl eta, F3_context *context);

#ifdef __cplusplus
}
#endif
