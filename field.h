#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

typedef struct field2 {
  dbl(*f)(dbl, dbl, void*);
  dvec2(*grad_f)(dbl, dbl, void*);
  void *context;
} field2_s;

dbl field2_f(field2_s const *field, dvec2 xy);
dvec2 field2_grad_f(field2_s const *field, dvec2 xy);

#ifdef __cplusplus
}
#endif
