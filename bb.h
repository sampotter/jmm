#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "jet.h"

typedef struct mesh3 mesh3_s;

void bb3_interp(dbl const f[2], dbl const Df[2], dbl const x[2], dbl c[4]);
void bb3_interp3(dbl const f[2], dbl const Df[2][3], dbl const x[2][3], dbl c[4]);
dbl bb3(dbl const *c, dbl const *b);
dbl dbb3(dbl const *c, dbl const *b, dbl const *a);
dbl d2bb3(dbl const *c, dbl const *b, dbl const *a);

void bb3tri_interp3(dbl const f[3], dbl const Df[3][3], dbl const x[3][3], dbl c[10]);
dbl bb3tri(dbl const *c, dbl const *b);
dbl dbb3tri(dbl const *c, dbl const *b, dbl const *a);
dbl d2bb3tri(dbl const *c, dbl const *b, dbl const *a1, dbl const *a2);

void bb3tet_interp3(dbl const f[4], dbl const Df[4][3], dbl const x[4][3],
                    dbl c[20]);
void bb3tet_for_cell(mesh3_s const *mesh, jet3 const *jet, size_t lc,
                     dbl c[20]);
dbl bb3tet(dbl const c[20], dbl const b[4]);
dbl dbb3tet(dbl const c[20], dbl const b[4], dbl const a[4]);
dbl d2bb3tet(dbl const c[20], dbl const b[4], dbl const a[2][4]);

#ifdef __cplusplus
}
#endif
