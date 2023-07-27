#pragma once

#include "common.h"
#include "cubic.h"
#include "jet.h"

typedef struct {
  dbl c[4];
} bb31;

void bb31_init_from_1d_data(bb31 *bb, dbl const f[2], dbl const Df[2], dbl const x[2]);
void bb31_init_from_3d_data(bb31 *bb, dbl const f[2], dbl const Df[2][3], dbl const x[2][3]);
void bb31_init_from_jets(bb31 *bb, jet31t const jet[2], dbl const x[2][3]);
void bb31_init_from_jet21t(bb31 *bb, jet21t const jet[2], dbl2 const x[2]);
void bb31_init_from_cubic(bb31 *bb, cubic_s const *cubic, dbl const x[2]);
dbl bb31_f(bb31 const *bb, dbl const *b);
dbl bb31_df(bb31 const *bb, dbl const *b, dbl const *a);
dbl bb31_d2f(bb31 const *bb, dbl const *b, dbl const *a);
void bb31_reverse(bb31 *bb);

typedef struct {
  dbl c[10];
} bb32;

void bb32_init_from_3d_data(bb32 *bb, dbl3 const f, dbl3 const Df[3], dbl3 const x[3]);
void bb32_init_from_jets(bb32 *bb, jet31t const jet[3], dbl const x[3][3]);
dbl bb32_f(bb32 const *bb, dbl const *b);
dbl bb32_df(bb32 const *bb, dbl const *b, dbl const *a);
dbl bb32_d2f(bb32 const *bb, dbl const *b, dbl const *a1, dbl const *a2);

typedef struct {
  dbl c[20];
} bb33;

void bb33_init_from_3d_data(bb33 *bb, dbl const f[4], dbl const Df[4][3], dbl const x[4][3]);
void bb33_init_from_cell_and_jets(bb33 *bb, mesh3_s const *mesh, jet31t const *jet, size_t lc);
void bb33_init_from_jets(bb33 *bb, jet31t const jet[4], dbl const x[4][3]);
dbl bb33_f(bb33 const *bb, dbl const b[4]);
dbl bb33_df(bb33 const *bb, dbl const b[4], dbl const a[4]);
dbl bb33_d2f(bb33 const *bb, dbl const b[4], dbl4 const a[2]);
bool bb33_convex_hull_brackets_value(bb33 const *bb, dbl value);
cubic_s bb33_restrict_along_interval(bb33 const *bb, dbl b0[4], dbl b1[4]);
