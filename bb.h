#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

dbl bb3(dbl const *c, dbl const *b);
dbl dbb3(dbl const *c, dbl const *b, dbl const *a);
dbl d2bb3(dbl const *c, dbl const *b, dbl const *a);

void bb3tri_interp3(dbl f[3], dbl Df[3][3], dbl x[3][3], dbl *c);
dbl bb3tri(dbl const *c, dbl const *b);
dbl dbb3tri(dbl const *c, dbl const *b, dbl const *a);
dbl d2bb3tri(dbl const *c, dbl const *b, dbl const *a1, dbl const *a2);

dbl bb3tet(dbl const *c, dbl const *b);

#ifdef __cplusplus
}
#endif
