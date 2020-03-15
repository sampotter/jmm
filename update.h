#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "cubic.h"
#include "vec.h"

typedef struct {
  cubic_s cubic;
  dvec2 xy, xy0, xy1;
} F3_context;

dbl F3(dbl eta, void *context);
dbl dF3_deta(dbl eta, void *context);

#ifdef __cplusplus
}
#endif
