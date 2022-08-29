#include <jmm/eik_F3.h>

// TODO: change naming conventions to make this take up less space:
//
// eta -> x
// th -> y
// d2S_deta_dth -> Sxy
//
// should make reading and debugging this much simpler
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

  dbl2 dxy; dbl2_sub(context->xy1, context->xy0, dxy);
  dbl2 xyeta; dbl2_saxpy(eta, dxy, context->xy0, xyeta);

  dbl2 lp; dbl2_sub(context->xy, xyeta, lp);
  dbl L = dbl2_norm(lp);
  dbl2_dbl_div_inplace(lp, L);
  dbl L_eta = -dbl2_dot(lp, dxy);

  dbl s0 = field2_f(context->slow, xyeta);
  dbl s1 = field2_f(context->slow, context->xy);
  dbl2 grad_s_eta; field2_grad_f(context->slow, xyeta, grad_s_eta);
  dbl s0_eta = dbl2_dot(grad_s_eta, dxy);

  context->F3 = T + (s0 + s1)*L/2;
  context->F3_eta = T_eta + (s0_eta*L + (s0 + s1)*L_eta)/2;
}
