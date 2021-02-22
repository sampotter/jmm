#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "vec.h"

struct field2 {
  dbl(*f)(dbl, dbl, void*);
  dvec2(*grad_f)(dbl, dbl, void*);
  void *context;
};

dbl field2_f(field2_s const *field, dvec2 xy);
dvec2 field2_grad_f(field2_s const *field, dvec2 xy);

typedef struct field3 {
  dbl(*f)(dvec3, void*);
  dvec3(*grad_f)(dvec3, void*);
  void *context;
} field3_s;

dbl field3_f(field3_s const *field, dvec3 p);
dvec3 field3_grad_f(field3_s const *field, dvec3 p);

#ifdef __cplusplus
}
#endif
