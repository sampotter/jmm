#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "jet.h"

typedef struct mesh3 mesh3_s;

typedef struct {
  dbl c[4];
} bb3;

void bb3_init_from_1d_data(bb3 *bb, dbl const f[2], dbl const Df[2], dbl const x[2]);
void bb3_init_from_3d_data(bb3 *bb, dbl const f[2], dbl const Df[2][3], dbl const x[2][3]);
dbl bb3_f(bb3 const *bb, dbl const *b);
dbl bb3_df(bb3 const *bb, dbl const *b, dbl const *a);
dbl bb3_d2f(bb3 const *bb, dbl const *b, dbl const *a);

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
