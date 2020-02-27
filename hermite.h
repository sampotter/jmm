#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "vec.h"

typedef enum {LAMBDA, MU} bicubic_variable;

typedef struct {
  dbl a[4];
} cubic;

dbl cubic_f(cubic const *cubic, dbl lam);
dbl cubic_df(cubic const *cubic, dbl lam);

typedef struct {
  dbl A[4][4];
} bicubic;

void bicubic_set_A(bicubic *bicubic, dbl data[4][4]);
cubic bicubic_restrict(bicubic const *bicubic, bicubic_variable var, int edge);
dbl bicubic_f(bicubic const *bicubic, dvec2 cc);
dbl bicubic_fx(bicubic const *bicubic, dvec2 cc);
dbl bicubic_fy(bicubic const *bicubic, dvec2 cc);
dbl bicubic_fxy(bicubic const *bicubic, dvec2 cc);

#ifdef __cplusplus
}
#endif
