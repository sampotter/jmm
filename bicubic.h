#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "cubic.h"
#include "mat.h"

typedef enum {LAMBDA, MU} bicubic_variable;

typedef struct bicubic {
  dmat44 A;
} bicubic_s;

void bicubic_set_data(bicubic_s *bicubic, dmat44 data);
void bicubic_set_data_from_ptr(bicubic_s *bicubic, dbl const *data_ptr);
cubic_s bicubic_restrict(bicubic_s const *bicubic, bicubic_variable var, int edge);
dbl bicubic_f(bicubic_s const *bicubic, dvec2 cc);
dbl bicubic_fx(bicubic_s const *bicubic, dvec2 cc);
dbl bicubic_fy(bicubic_s const *bicubic, dvec2 cc);
dbl bicubic_fxy(bicubic_s const *bicubic, dvec2 cc);

#ifdef __cplusplus
}
#endif
