#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "cubic.h"
#include "mat.h"

#if SJS_DEBUG
#include <stdbool.h>
#endif

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
dvec4 interpolate_fxy_at_verts(dvec4 fx, dvec4 fy, dbl h);
#if SJS_DEBUG
bool bicubic_valid(bicubic_s const *bicubic);
void bicubic_invalidate(bicubic_s *bicubic);
#endif

#ifdef __cplusplus
}
#endif
