#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "cubic.h"
#include "vec.h"

typedef struct {
  cubic_s cubic;
  dvec2 xy, xy0, xy1;
  dbl h;
} F_data;

dbl F(dbl lam, void *context);
dbl dF_dt(dbl lam, void *context);

#ifdef __cplusplus
}
#endif
