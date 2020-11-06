#pragma once

#include "mat.h"
#include "vec.h"

typedef struct costfunc3 {
  void (*func)(dvec3 const *, dbl *, void *);
  void (*grad)(dvec3 const *, dvec3 *, void *);
  void (*hess)(dvec3 const *, dmat33 *, void *);
  void *wkspc;
} costfunc3_s;

typedef struct baryopt_wkspc {
  costfunc3_s *costfunc;
  dbl beta;
  dbl sigma;
  dbl eps;
  dbl ftol;
  dbl xtol;
  dvec3 *x;
  dbl f;
} baryopt_wkspc_s;

bool baryopt_step(baryopt_wkspc_s *wkspc);
void baryopt_solve(baryopt_wkspc_s *wkspc);
