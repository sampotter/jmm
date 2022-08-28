#include "eik_S4.h"

void S4_compute(dbl th, S4_context *context) {
  dbl2 t = {cos(th), sin(th)};
  dbl2 t_th = {-t[1], t[0]};

  dbl2 t_minus_t0; dbl2_sub(t, context->t0, t_minus_t0);

  dbl2 xym;
  dbl2_saxpy(-context->L/8, t_minus_t0, context->xy_xy0_avg, xym);

  dbl sm = field2_f(context->slow, xym);

  dbl2 t_plus_t0; dbl2_add(t, context->t0, t_minus_t0);

  dbl2 tm;
  dbl2_lincomb(1.5, context->lp, -0.25, t_plus_t0, tm);

  dbl tmnorm = dbl2_norm(tm);

  dbl2 gsm; field2_grad_f(context->slow, xym, gsm);
  dbl2 xym_th; dbl2_dbl_mul(t_th, -context->L/8, xym_th);
  dbl sm_th = dbl2_dot(gsm, xym_th);

  dbl2 tm_th; dbl2_dbl_mul(t_th, -0.25, tm_th);
  dbl tmnorm_th = dbl2_dot(tm, tm_th)/tmnorm;

  context->S4 = (context->s + 4*sm*tmnorm + context->s0)/6;
  context->S4_th = 2*(sm_th*tmnorm + sm*tmnorm_th)/3;
}
