#include "eik_F3.h"

// TODO: change naming conventions to make this take up less space:
//
// eta -> x
// th -> y
// d2S_deta_dth -> Sxy
//
// should make reading and debugging this much simplerj
//
// TODO: also make the change:
//
// xy -> p or r or z or something
//
// TODO: accumulate temporary variables in the context so that we
// don't recompute tons and tons of things
//
// TODO: make sure we're doing things as simply as possibly in terms
// of evaluating derivatives recursively and with minimal work

void F3_compute(dbl eta, F3_context *context) {
  dbl T = cubic_f(&context->T_cubic, eta);
  dbl T_eta = cubic_df(&context->T_cubic, eta);

  dvec2 dxy = dvec2_sub(context->xy1, context->xy0);
  dvec2 xyeta = dvec2_saxpy(eta, dxy, context->xy0);

  dvec2 lp = dvec2_sub(context->xy, xyeta);
  dbl L = dvec2_norm(lp);
  lp = dvec2_dbl_div(lp, L);
  dbl L_eta = -dvec2_dot(lp, dxy);

  dbl s0 = field2_f(context->slow, xyeta);
  dbl s1 = field2_f(context->slow, context->xy);
  dbl s0_eta = dvec2_dot(field2_grad_f(context->slow, xyeta), dxy);

  context->F3 = T + (s0 + s1)*L/2;
  context->F3_eta = T_eta + (s0_eta*L + (s0 + s1)*L_eta)/2;
}
