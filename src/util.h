#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

dbl clamp(dbl x, dbl a, dbl b);
int sgn(dbl x);
dbl q(dbl const A[3][3], dbl const b[3], dbl c, dbl const x[3]);
void Dq(dbl const A[3][3], dbl const b[3], dbl const x[3], dbl Dq[3]);
int compar_size_t(size_t const *i, size_t const *j);
int signum(dbl x);
dbl shrink(dbl x, dbl eps);

#ifdef __cplusplus
}
#endif
