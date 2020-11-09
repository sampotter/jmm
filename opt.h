#pragma once

#include "def.h"

typedef struct costfunc2 {
  void (*func)(dbl const *, dbl *, void *);
  void (*grad)(dbl const *, dbl *, void *);
  void (*hess)(dbl const *, dbl *, void *);
  void *wkspc;
} costfunc2_s;

void eqp_bary_2_0(dbl const *G, dbl const *c, dbl *x);
void eqp_bary_2_1(dbl const *G, dbl const *c, dbl *x);
void eqp_bary_2_2(dbl const *G, dbl const *c, dbl *x);
void iqp_bary_2(dbl const *G, dbl const *c, dbl const *x0, dbl *x, bool *error,
                dbl tol, int niters);
void sqp_bary_3_2(costfunc2_s const *costfunc, dbl const *xinit, dbl *x,
                  dbl *f, bool *error, dbl tol, int niters);
