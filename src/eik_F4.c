#include "eik_F4.h"

#include <assert.h>
#include <stdio.h>

#include "macros.h"
#include "util.h"

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
  dbl2 tmp;

  dbl T = cubic_f(&context->T_cubic, eta);
  dbl T_eta = cubic_df(&context->T_cubic, eta);
  // dbl T_eta_eta = cubic_d2f(&context->T_cubic, eta);

  // t0 is normalized by definition
  dbl2 t0 = {
    cubic_f(&context->Tx_cubic, eta),
    cubic_f(&context->Ty_cubic, eta)
  };
  dbl gradTnorm = dbl2_norm(t0);
  dbl2_dbl_div_inplace(t0, gradTnorm);

  dbl2 t0_eta = {
    cubic_df(&context->Tx_cubic, eta),
    cubic_df(&context->Ty_cubic, eta)
  };
  dbl2_dbl_div(t0_eta, gradTnorm, tmp);
  dbl2_cproj(t0, tmp, t0_eta);

  // t1 is normalized by definition
  dbl2 t1 = {cos(th), sin(th)};

  // avoid recomputing sin and cos of th
  dbl2 t1_th = {-t1[1], t1[0]};

  dbl2 dxy; dbl2_sub(context->xy1, context->xy0, dxy);
  dbl2 xyeta; dbl2_saxpy(eta, dxy, context->xy0, xyeta);

  dbl2 lp; dbl2_sub(context->xy, xyeta, lp);
  dbl L = dbl2_norm(lp);
  dbl2_dbl_div_inplace(lp, L);
  dbl L_eta = -dbl2_dot(lp, dxy);
  // dbl L_eta_eta = (dvec2_norm_sq(dxy) - L_eta*L_eta)/L;

  dbl2 t1_minus_t0; dbl2_sub(t1, t0, t1_minus_t0);

  dbl2_add(context->xy, xyeta, tmp);
  dbl2 xym; dbl2_lincomb(0.5, tmp, -L/8, t1_minus_t0, xym);

  dbl s0 = field2_f(context->slow, xyeta);
  dbl s1 = field2_f(context->slow, context->xy);
  dbl sm = field2_f(context->slow, xym);

  // tm is unnormalized by definition, but its norm will be close to 1
  // because of the quasiuniform parametrization
  dbl2_add(t0, t1, tmp);
  dbl2 tm; dbl2_lincomb(1.5, lp, -0.25, tmp, tm);
  dbl tmnorm = dbl2_norm(tm);

  dbl2_dbl_div(dxy, -L, tmp);
  dbl2 lp_eta; dbl2_cproj(lp, tmp, lp_eta);

  dbl2 tm_eta; dbl2_lincomb(1.5, lp_eta, -0.25, t0_eta, tm_eta);
  dbl2 tm_th; dbl2_dbl_mul(t1_th, -0.25, tm_th);

  dbl tmnorm_eta = dbl2_dot(tm, tm_eta)/tmnorm;
  dbl tmnorm_th = dbl2_dot(tm, tm_th)/tmnorm;

  dbl2_lincomb(L, t0_eta, -L_eta, t1_minus_t0, tmp);
  dbl2 xym_eta; dbl2_lincomb(0.5, dxy, 0.125, tmp, xym_eta);
  dbl2 xym_th; dbl2_dbl_mul(t1_th, -L/8, xym_th);

  dbl2 gseta; field2_grad_f(context->slow, xyeta, gseta);
  dbl2 gsm; field2_grad_f(context->slow, xym, gsm);

  // dmat22 Hsm = field2_hess_f(context->slow, xym);

  dbl s0_eta = dbl2_dot(gseta, dxy);

  dbl sm_eta = dbl2_dot(gsm, xym_eta);
  dbl sm_th = dbl2_dot(gsm, xym_th);
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

void F4_get_grad(F4_context const *context, dbl2 g) {
  g[0] = context->F4_eta;
  g[1] = context->F4_th;
}

void F4_hess_fd(dbl eta, dbl th, dbl eps, F4_context *context, dbl22 H) {
  dbl2 gp, gm;

  F4_compute(eta + eps, th, context);
  F4_get_grad(context, gp);

  F4_compute(eta - eps, th, context);
  F4_get_grad(context, gm);

  dbl2_sub(gp, gm, H[0]);
  dbl2_dbl_div_inplace(H[0], 2*eps);

  F4_compute(eta, th + eps, context);
  F4_get_grad(context, gp);

  F4_compute(eta, th - eps, context);
  F4_get_grad(context, gm);

  dbl2_sub(gp, gm, H[1]);
  dbl2_dbl_div_inplace(H[1], 2*eps);

  // symmetrize the hessian before returning
  H[0][1] = H[1][0] = (H[0][1] + H[1][0])/2;
}

static void make_hessian_positive_definite(dbl22 H) {
  dbl2 lam;
  dbl22_eigvals(H, lam);
  assert(lam[0] > 0);
  if (lam[1] < 0)
    dbl22_perturb(H, -2*lam[1]);
}

void F4_bfgs_init(dbl eta, dbl th, dbl2 x0, dbl2 g0, dbl22 H0,
                  F4_context *context) {
  x0[0] = eta;
  x0[1] = th;

  F4_compute(x0[0], x0[1], context);

  F4_get_grad(context, g0);

  F4_hess_fd(eta, th, 1e-7, context, H0);
  make_hessian_positive_definite(H0);
  dbl22_invert(H0);
}

static void get_hess_bfgs(dbl2 const xk, dbl2 const xk1,
                          dbl2 const gk, dbl2 const gk1,
                          dbl22 const Hk, dbl22 Hk1) {
  if (dbl2_dist(xk1, xk) < EPS) {
    dbl22_copy(Hk, Hk1);
    return;
  }

  dbl2 sk; dbl2_sub(xk1, xk, sk);
  dbl2 yk; dbl2_sub(gk1, gk, yk);
  dbl rhok = 1/dbl2_dot(yk, sk);

  dbl22 tmp1, tmp2;

  dbl2_outer(sk, yk, tmp1);
  dbl22_dbl_mul_inplace(tmp1, -rhok);
  dbl22_perturb(tmp1, 1);

  dbl22_mul(tmp1, Hk, tmp2);
  dbl22_transpose(tmp1);
  dbl22_mul(tmp2, tmp1, Hk1);

  dbl2_outer(sk, sk, tmp1);
  dbl22_dbl_mul_inplace(tmp1, rhok);
  dbl22_add_inplace(Hk1, tmp1);

  make_hessian_positive_definite(Hk1);
}

bool F4_bfgs_step(dbl2 const xk, dbl2 const gk, dbl22 const Hk,
                  dbl2 xk1, dbl2 gk1, dbl22 Hk1,
                  F4_context *context) {
  assert(dbl2_maxnorm(gk) > EPS);

  dbl2 pk; dbl22_dbl2_mul(Hk, gk, pk);
  dbl2_negate(pk);

  // Verify that pk is a descent direction.
  dbl pk_dot_gk = dbl2_dot(pk, gk);
  assert(pk_dot_gk < 0);

  // Scale the step so that 0 <= eta <= 1.
  dbl t = pk[0] > 0. ? 1. : 0.;
  t -= xk[0];
  t /= pk[0];
  t = clamp(t, 0, 1);

//  dbl const c1 = 1e-4;
  dbl const rho = 0.5;

  /**
   * Do an inexact backtracking line search to find `t` such that the
   * sufficient decrease conditions are satisfied.
   */
  if (t > EPS && pk_dot_gk < 0) {
    dbl fk = context->F4, fk1;
    while (true) {
      dbl2_saxpy(t, pk, xk, xk1);
      F4_compute(SPLAT2(xk1), context);
      fk1 = context->F4;
      if (fk1 <= fk + EPS) {
        break;
      } else {
        t *= rho;
      }
    }
  } else {
    // TODO: don't do this if t = 0! no need, just set xk1 = xk
    dbl2_saxpy(t, pk, xk, xk1);
    F4_compute(SPLAT2(xk1), context);
  }

  // Now, compute a new gradient and do the DFP update to update our
  // approximation of the inverse Hessian.
  //
  // TODO: this is also wasteful when t = 0! see above
  F4_get_grad(context, gk1);
  update_bfgs(xk1, xk, gk1, gk, Hk, Hk1);

  if (t < 1 && (fabs(xk1[0]) < EPS || fabs(1 - xk1[0]) < EPS)) {
    dbl2_copy(gk, gk1);
    dbl2_copy(xk, xk1);

    // Reset t to be a line search parameter along `th`
    t = 1;

    // Reset pk and pk_dot_gk for descent along `th`
    dbl F4_th_th = Hk1[0][0]/dbl22_det(Hk1);
    dbl F4_th = gk[1];
    pk[1] = -F4_th/F4_th_th;
    pk_dot_gk = pk[1]*gk[1];

    // Backtracking line search
    dbl fk = context->F4, fk1;
    while (true) {
      xk1[1] = xk[1] + t*pk[1];
      F4_compute(xk1[0], xk1[1], context);
      fk1 = context->F4;
      if (fk1 <= fk + EPS) {
        break;
      } else {
        t *= rho;
      }
    }

    F4_get_grad(context, gk1);
    update_bfgs(xk1, xk, gk1, gk, Hk, Hk1);
  }

  return true;
}
