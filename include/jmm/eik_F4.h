#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#include "cubic.h"
#include "field.h"
#include "mat.h"
#include "vec.h"

// TODO: we want to make the following changes:
//
// instead of having these contexts, we just pass this data to
// functions below which directly compute "complete 2-jets" (we can
// think of the Newton iteration as being an algorithm that
// iteratively refines a 2-jet of a scalar field)
//
// the benefit: everything below becomes "pure", and there's no
// mutating state

typedef struct {
  // Inputs:
  cubic_s T_cubic, Tx_cubic, Ty_cubic;
  dbl2 xy, xy0, xy1;
  field2_s const *slow;

  // Outputs:
  dbl F4;
  dbl F4_eta;
  dbl F4_th;
} F4_context;

void F4_compute(dbl eta, dbl th, F4_context *context);
void F4_get_grad(F4_context const *context, dbl2 g);
void F4_hess_fd(dbl eta, dbl th, dbl eps, F4_context *context, dbl22 H);
void F4_bfgs_init(dbl eta, dbl th, dbl2 x0, dbl2 g0, dbl22 H0,
                  F4_context *context);
bool F4_bfgs_step(dbl2 const xk, dbl2 const gk, dbl22 const Hk,
                  dbl2 xk1, dbl2 gk1, dbl22 Hk1, F4_context *context);

#ifdef __cplusplus
}
#endif
