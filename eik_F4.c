#include "eik_F4.h"

#include <assert.h>
#include <stdio.h>

#include "math.h"

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

void F4_compute(dbl eta, dbl th, F4_context *context) {
  dbl T = cubic_f(&context->T_cubic, eta);
  dbl T_eta = cubic_df(&context->T_cubic, eta);
  // dbl T_eta_eta = cubic_d2f(&context->T_cubic, eta);

  // t0 is normalized by definition
  dvec2 t0 = {
    .x = cubic_f(&context->Tx_cubic, eta),
    .y = cubic_f(&context->Ty_cubic, eta)
  };
  dbl gradTnorm = dvec2_norm(t0);
  t0 = dvec2_dbl_div(t0, gradTnorm);

  dvec2 t0_eta = {
    .x = cubic_df(&context->Tx_cubic, eta),
    .y = cubic_df(&context->Ty_cubic, eta)
  };
  t0_eta = dvec2_cproj(t0, dvec2_dbl_div(t0_eta, gradTnorm));

  // t1 is normalized by definition
  dvec2 t1 = {.x = cos(th), .y = sin(th)};

  // avoid recomputing sin and cos of th
  dvec2 t1_th = {.x = -t1.y, .y = t1.x};

  dvec2 dxy = dvec2_sub(context->xy1, context->xy0);
  dvec2 xyeta = dvec2_saxpy(eta, dxy, context->xy0);

  dvec2 lp = dvec2_sub(context->xy, xyeta);
  dbl L = dvec2_norm(lp);
  lp = dvec2_dbl_div(lp, L);
  dbl L_eta = -dvec2_dot(lp, dxy);
  // dbl L_eta_eta = (dvec2_norm_sq(dxy) - L_eta*L_eta)/L;

  dvec2 t1_minus_t0 = dvec2_sub(t1, t0);

  dvec2 xym = dvec2_sub(
    dvec2_dbl_div(dvec2_add(context->xy, xyeta), 2),
    dvec2_dbl_mul(t1_minus_t0, L/8)
  );

  dbl s0 = field2_f(context->slow, xyeta);
  dbl s1 = field2_f(context->slow, context->xy);
  dbl sm = field2_f(context->slow, xym);

  // tm is unnormalized by definition, but its norm will be close to 1
  // because of the quasiuniform parametrization
  dvec2 tm = dvec2_sub(
    dvec2_dbl_mul(lp, 1.5),
    dvec2_dbl_mul(dvec2_add(t0, t1), 0.25)
  );
  dbl tmnorm = dvec2_norm(tm);

  dvec2 lp_eta = dvec2_cproj(lp, dvec2_dbl_div(dxy, -L));

  dvec2 tm_eta = dvec2_sub(
    dvec2_dbl_mul(lp_eta, 1.5),
    dvec2_dbl_mul(t0_eta, 0.25)
  );
  dvec2 tm_th = dvec2_dbl_mul(t1_th, -0.25);

  dbl tmnorm_eta = dvec2_dot(tm, tm_eta)/tmnorm;
  dbl tmnorm_th = dvec2_dot(tm, tm_th)/tmnorm;

  dvec2 xym_eta = dvec2_add(
    dvec2_dbl_div(dxy, 2),
    dvec2_dbl_div(
      dvec2_sub(
        dvec2_dbl_mul(t0_eta, L),
        dvec2_dbl_mul(t1_minus_t0, L_eta)
      ),
      8
    )
  );
  dvec2 xym_th = dvec2_dbl_mul(t1_th, -L/8);

  dvec2 gseta = field2_grad_f(context->slow, xyeta);
  dvec2 gsm = field2_grad_f(context->slow, xym);

  // dmat22 Hsm = field2_hess_f(context->slow, xym);

  dbl s0_eta = dvec2_dot(gseta, dxy);

  dbl sm_eta = dvec2_dot(gsm, xym_eta);
  dbl sm_th = dvec2_dot(gsm, xym_th);
  // dbl sm_eta_eta = dvec2_dot(dmat22_dvec2_mul(Hsm, xym_eta), xym_eta)
  //   + dvec2_dot(gsm, xm_eta_eta);

  dbl S = (s0 + s1 + 4*sm*tmnorm)/6;
  dbl S_eta = (s0_eta + 4*(sm_eta*tmnorm + sm*tmnorm_eta))/6;
  dbl S_th = 2*(sm_th*tmnorm + sm*tmnorm_th)/3;
  // dbl S_eta_eta = (
  //   s0_eta_eta
  //   + 4*(sm_eta_eta*tmnorm + 2*sm_eta*tmnorm_eta + sm*tmnorm_eta_eta))/6;
  // dbl S_eta_th = 2*(
  //   sm_eta_th*tmnorm + sm_th*tmnorm_eta +
  //   sm_eta*tmnorm_th + sm*tmnorm_eta_th)/3;
  // dbl S_th_th = 2*(sm_th_th*tmnorm + 2*sm_th*tmnorm_th + sm*tmnorm_th)/3;

  context->F4 = T + L*S;
  context->F4_eta = T_eta + L*S_eta + S*L_eta;
  context->F4_th = L*S_th;
  // context->F4_eta_eta = T_eta_eta + L_eta_eta*S + 2*L_eta*S_eta + L*S_eta_eta;
  // context->F4_eta_th = L_eta*S_th + L*S_eta_th;
  // context->F4_th_th = L*S_th_th;
}

dvec2 F4_get_grad(F4_context const *context) {
  return (dvec2) {context->F4_eta, context->F4_th};
}

dmat22 F4_hess_fd(dbl eta, dbl th, dbl eps, F4_context *context) {
  dmat22 hess;

  F4_compute(eta + eps, th, context);
  hess.rows[0] = F4_get_grad(context);
  F4_compute(eta - eps, th, context);
  hess.rows[0] = dvec2_sub(hess.rows[0], F4_get_grad(context));
  hess.rows[0] = dvec2_dbl_div(hess.rows[0], 2*eps);

  F4_compute(eta, th + eps, context);
  hess.rows[1] = F4_get_grad(context);
  F4_compute(eta, th - eps, context);
  hess.rows[1] = dvec2_sub(hess.rows[1], F4_get_grad(context));
  hess.rows[1] = dvec2_dbl_div(hess.rows[1], 2*eps);

  return hess;
}

void F4_bfgs_init(dbl eta, dbl th, dvec2 *x0, dvec2 *g0, dmat22 *H0,
                  F4_context *context) {
  *x0 = (dvec2) {.x = eta, .y = th};
  F4_compute(x0->x, x0->y, context);
  *g0 = F4_get_grad(context);
  *H0 = F4_hess_fd(eta, th, 1e-7, context);
  /**
   * If the hessian is indefinite, we want to perturb it by a constant
   * multiple times the identity so that it it's positive definite. We
   * assume that the hessian can't be negative definite (hopefully
   * this can't happen... something really weird would have had to
   * have happened.
   */
  {
    dbl lam1, lam2;
    dmat22_eigvals(H0, &lam1, &lam2);
    assert(lam1 > 0);
    if (lam2 < 0) {
      H0->data[0][0] -= 2*lam2;
      H0->data[1][1] -= 2*lam2;
    }
  }
  dmat22_invert(H0);
}

static void update_dfp(dvec2 xk1, dvec2 xk, dvec2 gk1, dvec2 gk, dmat22 Hk,
                       dmat22 *Hk1) {
  if (xk1.x == xk.x && xk1.y == xk.y) {
    return;
  }
  dvec2 sk = dvec2_sub(xk1, xk);
  dvec2 yk = dvec2_sub(gk1, gk);
  dvec2 tmp1 = dmat22_dvec2_mul(Hk, yk);
  dbl tmp2 = dvec2_dot(yk, tmp1);
  dmat22 tmp3 = dvec2_outer(tmp1, tmp1);
  *Hk1 = dmat22_sub(Hk, dmat22_dbl_div(tmp3, tmp2));
  tmp3 = dvec2_outer(sk, sk);
  tmp2 = dvec2_dot(yk, sk);
  *Hk1 = dmat22_add(*Hk1, dmat22_dbl_div(tmp3, tmp2));
}

// static void update_bfgs(dvec2 xk1, dvec2 xk, dvec2 gk1, dvec2 gk, dmat22 Hk,
//                         dmat22 *Hk1) {
//   dvec2 sk = dvec2_sub(xk1, xk);
//   dvec2 yk = dvec2_sub(gk1, gk);
//   dbl rhok = 1/dvec2_dot(yk, sk);

//   dmat22 tmp = dvec2_outer(sk, yk);
//   tmp = dmat22_dbl_mul(tmp, -rhok);
//   tmp.data[0][0] += 1;
//   tmp.data[1][1] += 1;

//   *Hk1 = dmat22_mul(tmp, Hk);
//   dmat22_transpose(&tmp);
//   *Hk1 = dmat22_mul(*Hk1, tmp);

//   tmp = dvec2_outer(sk, sk);
//   tmp = dmat22_dbl_mul(tmp, rhok);

//   *Hk1 = dmat22_add(*Hk1, tmp);
// }

bool F4_bfgs_step(dvec2 xk, dvec2 gk, dmat22 Hk,
                  dvec2 *xk1, dvec2 *gk1, dmat22 *Hk1,
                  F4_context *context) {
  dvec2 pk = dmat22_dvec2_mul(Hk, gk);
  dvec2_negate(&pk);

  // Verify that pk is a descent direction.
  dbl pk_dot_gk = dvec2_dot(pk, gk);
  assert(pk_dot_gk < 0);

  // Scale the step so that 0 <= eta <= 1.
  dbl t = pk.x > 0. ? 1. : 0.;
  t -= xk.x;
  t /= pk.x;
  t = clamp(t, 0, 1);

  /**
   * Do an inexact backtracking line search to find `t` such that the
   * sufficient decrease conditions are satisfied.
   */
  if (t > EPS && pk_dot_gk < -1e-13) {
    dbl const c1 = 1e-4;
    dbl const rho = 0.9;

    dbl fk = context->F4, fk1;
    while (true) {
      *xk1 = dvec2_add(xk, dvec2_dbl_mul(pk, t));
      F4_compute(xk1->x, xk1->y, context);
      fk1 = context->F4;
      if (fk1 <= fk + c1*t*pk_dot_gk) {
        break;
      } else {
        t *= rho;
      }
    }
  } else {
    *xk1 = dvec2_add(xk, dvec2_dbl_mul(pk, t));
    F4_compute(xk1->x, xk1->y, context);
  }

  // Now, compute a new gradient and do the DFP update to update our
  // approximation of the inverse Hessian.
  *gk1 = F4_get_grad(context);
  update_dfp(*xk1, xk, *gk1, gk, Hk, Hk1);

  if (t < 1 && (fabs(xk1->x) < EPS || fabs(1 - xk1->x) < EPS)) {
    dbl F4_th_th = Hk1->data[0][0]/dmat22_det(Hk1);
    dbl F4_th = gk1->y;
    xk1->y -= F4_th/F4_th_th;

    F4_compute(xk1->x, xk1->y, context);
    *gk1 = F4_get_grad(context);
    update_dfp(*xk1, xk, *gk1, gk, Hk, Hk1);
  }

  return true;
}
