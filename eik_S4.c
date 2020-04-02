#include "eik_S4.h"

void S4_compute(dbl th, S4_context *context) {
  dvec2 t = {.x = cos(th), .y = sin(th)};
  dvec2 t_th = {.x = -t.y, .y = t.x};

  dvec2 t_minus_t0 = dvec2_sub(t, context->t0);

  dvec2 xym = dvec2_sub(
    context->xy_xy0_avg,
    dvec2_dbl_mul(t_minus_t0, context->L/8)
  );

  dbl sm = field2_f(context->slow, xym);
  dvec2 tm = dvec2_sub(
    dvec2_dbl_mul(context->lp, 1.5),
    dvec2_dbl_mul(dvec2_add(context->t0, t), 0.25)
  );
  dbl tmnorm = dvec2_norm(tm);

  dvec2 gsm = field2_grad_f(context->slow, xym);
  dvec2 xym_th = dvec2_dbl_mul(t_th, -context->L/8);
  dbl sm_th = dvec2_dot(gsm, xym_th);

  dvec2 tm_th = dvec2_dbl_mul(t_th, -0.25);
  dbl tmnorm_th = dvec2_dot(tm, tm_th)/tmnorm;

  context->S4 = (context->s + 4*sm*tmnorm + context->s0)/6;
  context->S4_th = 2*(sm_th*tmnorm + sm*tmnorm_th)/3;
}
