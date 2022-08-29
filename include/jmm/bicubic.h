#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "cubic.h"
#include "mat.h"

typedef enum {LAMBDA, MU} bicubic_variable;

typedef struct bicubic {
  dbl44 A;
} bicubic_s;

void bicubic_set_data(bicubic_s *bicubic, dbl44 data);
void bicubic_set_data_from_ptr(bicubic_s *bicubic, dbl const *data_ptr);
cubic_s bicubic_get_f_on_edge(bicubic_s const *bicubic, bicubic_variable var, int edge);
cubic_s bicubic_get_fx_on_edge(bicubic_s const *bicubic, bicubic_variable var, int edge);
cubic_s bicubic_get_fy_on_edge(bicubic_s const *bicubic, bicubic_variable var, int edge);
cubic_s bicubic_get_fxx_on_edge(bicubic_s const *bicubic, bicubic_variable var, int edge);
cubic_s bicubic_get_fyy_on_edge(bicubic_s const *bicubic, bicubic_variable var, int edge);
dbl bicubic_f(bicubic_s const *bicubic, dbl2 cc);
dbl bicubic_fx(bicubic_s const *bicubic, dbl2 cc);
dbl bicubic_fy(bicubic_s const *bicubic, dbl2 cc);
dbl bicubic_fxx(bicubic_s const *bicubic, dbl2 cc);
dbl bicubic_fxy(bicubic_s const *bicubic, dbl2 cc);
dbl bicubic_fyy(bicubic_s const *bicubic, dbl2 cc);
void interpolate_fxy_at_verts(dbl4 const fx, dbl4 const fy, dbl h, dbl4 fxy);
bool bicubic_valid(bicubic_s const *bicubic);
void bicubic_invalidate(bicubic_s *bicubic);

#ifdef __cplusplus
}
#endif
