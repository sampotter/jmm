#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

void bb3_interp(dbl const f[2], dbl const Df[2], dbl const x[2], dbl c[4]);
void bb3_interp3(dbl const f[2], dbl const Df[2][3], dbl const x[2][3], dbl c[4]);
dbl bb3(dbl const *c, dbl const *b);
dbl dbb3(dbl const *c, dbl const *b, dbl const *a);
dbl d2bb3(dbl const *c, dbl const *b, dbl const *a);

void bb3tri_interp3(dbl const f[3], dbl const Df[3][3], dbl const x[3][3], dbl c[10]);
dbl bb3tri(dbl const *c, dbl const *b);
dbl dbb3tri(dbl const *c, dbl const *b, dbl const *a);
dbl d2bb3tri(dbl const *c, dbl const *b, dbl const *a1, dbl const *a2);

dbl bb3tet(dbl const *c, dbl const *b);

#ifdef __cplusplus
}
#endif
