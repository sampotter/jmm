#include "update.h"

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

dvec2 get_lp(dvec2 xy, dvec2 xy0, dvec2 xy1, dbl eta) {
  return dvec2_sub(xy, dvec2_ccomb(xy0, xy1, eta));
}

dbl F3(dbl eta, void *context) {
  F3_context *ctx = (F3_context *)context;

  dbl T = cubic_f(&ctx->cubic, eta);

  dvec2 xyeta = dvec2_ccomb(ctx->xy0, ctx->xy1, eta);
  dbl L = dvec2_dist(ctx->xy, xyeta);

  return T + L;
}

dbl dF3_deta(dbl eta, void *context) {
  F3_context *ctx = (F3_context *)context;

  dbl dT_deta = cubic_df(&ctx->cubic, eta);

  dvec2 xyeta = dvec2_ccomb(ctx->xy0, ctx->xy1, eta);
  dvec2 xyeta_minus_xy = dvec2_sub(xyeta, ctx->xy);
  dbl L = dvec2_norm(xyeta_minus_xy);

  dvec2 dxy = dvec2_sub(ctx->xy1, ctx->xy0);
  dbl dL_deta = dvec2_dot(dxy, xyeta_minus_xy)/L;

  return dT_deta + dL_deta;
}

dbl F4(dbl eta, dbl th, void *context) {
  F4_context *ctx = (F4_context *)context;

  dbl T = cubic_f(&ctx->T, eta);

  dvec2 lp = get_lp(ctx->xy, ctx->xy0, ctx->xy1, eta);
  dbl L = dvec2_norm(lp);
  lp = dvec2_dbl_div(lp, L);

  dvec2 t0 = {
    .x = cubic_f(&ctx->Tx, eta),
    .y = cubic_f(&ctx->Ty, eta)
  };
  dvec2_normalize(&t0);

  dvec2 t1 = {.x = cos(th), .y = sin(th)};

  dvec2 t = dvec2_sub(
    dvec2_dbl_mul(lp, 1.5),
    dvec2_dbl_mul(dvec2_add(t0, t1), 0.25)
  );

  dbl S = (1.0 + 2.0*dvec2_norm(t))/3.0;

  return T + L*S;
}

dbl dF4_deta(dbl eta, dbl th, void *context) {
  F4_context *ctx = (F4_context *)context;

  dbl dT_deta = cubic_df(&ctx->T, eta);

  dvec2 lp = get_lp(ctx->xy, ctx->xy0, ctx->xy1, eta);
  dbl L = dvec2_norm(lp);
  lp = dvec2_dbl_div(lp, L);

  dvec2 dxy = dvec2_sub(ctx->xy1, ctx->xy0);
  dvec2 dlp_deta = dvec2_cproj(lp, dvec2_dbl_div(dxy, -L));

  dbl dL_deta = -dvec2_dot(dxy, lp);

  dvec2 t0 = {
    .x = cubic_f(&ctx->Tx, eta),
    .y = cubic_f(&ctx->Ty, eta)
  };
  dbl t0_norm = dvec2_norm(t0);
  t0 = dvec2_dbl_div(t0, t0_norm);

  dvec2 t1 = {.x = cos(th), .y = sin(th)};

  dvec2 t = dvec2_sub(
    dvec2_dbl_mul(lp, 1.5),
    dvec2_dbl_mul(dvec2_add(t0, t1), 0.25)
  );
  dbl tnorm = dvec2_norm(t);
  t = dvec2_dbl_div(t, tnorm);

  dbl S = (1.0 + 2.0*tnorm)/3.0;

  dvec2 dt0_deta = {
    .x = cubic_df(&ctx->Tx, eta),
    .y = cubic_df(&ctx->Ty, eta)
  };
  dt0_deta = dvec2_cproj(t0, dvec2_dbl_div(dt0_deta, t0_norm));

  dbl dS_deta = dvec2_dot(
    t,
    dvec2_sub(
      dlp_deta,
      dvec2_dbl_div(dt0_deta, 6.0)
    )
  );

  return dT_deta + dL_deta*S + L*dS_deta;
}

dbl dF4_dth(dbl eta, dbl th, void *context) {
  F4_context *ctx = (F4_context *)context;

  dvec2 dt1_dth = {.x = -sin(th), .y = cos(th)};

  dvec2 lp = get_lp(ctx->xy, ctx->xy0, ctx->xy1, eta);
  dbl L = dvec2_norm(lp);
  lp = dvec2_dbl_div(lp, L);

  dvec2 t0 = {
    .x = cubic_f(&ctx->Tx, eta),
    .y = cubic_f(&ctx->Ty, eta)
  };
  dvec2_normalize(&t0);

  dvec2 t1 = {.x = cos(th), .y = sin(th)};

  dvec2 t = dvec2_sub(
    dvec2_dbl_mul(lp, 1.5),
    dvec2_dbl_mul(dvec2_add(t0, t1), 0.25));
  dvec2_normalize(&t);

  dbl dS_dth = (2.0/3.0)*dvec2_dot(
    dvec2_dbl_mul(dt1_dth, -0.25),
    t
  );

  return L*dS_dth;
}

dvec2 grad_F4(dbl eta, dbl th, void *context) {
  F4_context *ctx = (F4_context *)context;

  dbl dT_deta = cubic_df(&ctx->T, eta);

  dvec2 lp = get_lp(ctx->xy, ctx->xy0, ctx->xy1, eta);
  dbl L = dvec2_norm(lp);
  lp = dvec2_dbl_div(lp, L);

  dvec2 dxy = dvec2_sub(ctx->xy1, ctx->xy0);
  dvec2 dlp_deta = dvec2_cproj(lp, dvec2_dbl_div(dxy, -L));

  dbl dL_deta = -dvec2_dot(dxy, lp);

  dvec2 t0 = {
    .x = cubic_f(&ctx->Tx, eta),
    .y = cubic_f(&ctx->Ty, eta)
  };
  dbl t0_norm = dvec2_norm(t0);
  t0 = dvec2_dbl_div(t0, t0_norm);

  dvec2 t1 = {.x = cos(th), .y = sin(th)};

  dvec2 t = dvec2_sub(
    dvec2_dbl_mul(lp, 1.5),
    dvec2_dbl_mul(dvec2_add(t0, t1), 0.25)
  );
  dbl tnorm = dvec2_norm(t);
  t = dvec2_dbl_div(t, tnorm);

  dbl S = (1.0 + 2.0*tnorm)/3.0;

  dvec2 dt0_deta = {
    .x = cubic_df(&ctx->Tx, eta),
    .y = cubic_df(&ctx->Ty, eta)
  };
  dt0_deta = dvec2_cproj(t0, dvec2_dbl_div(dt0_deta, t0_norm));

  dbl dS_deta = dvec2_dot(
    t,
    dvec2_sub(
      dlp_deta,
      dvec2_dbl_div(dt0_deta, 6.0)
    )
  );

  dvec2 dt1_dth = {.x = -sin(th), .y = cos(th)};

  dbl dS_dth = (2.0/3.0)*dvec2_dot(
    dvec2_dbl_mul(dt1_dth, -0.25),
    t
  );

  return (dvec2) {
    .x = dT_deta + dL_deta*S + L*dS_deta,
    .y = L*dS_dth
  };
}

dmat22 hess_F4(dbl eta, dbl th, void *context) {
  F4_context *ctx = (F4_context *)context;

  dvec2 dxy = dvec2_sub(ctx->xy1, ctx->xy0);
  dvec2 lp = get_lp(ctx->xy, ctx->xy0, ctx->xy1, eta);

  dbl L = dvec2_norm(lp);
  dbl dL_deta = -dvec2_dot(dxy, lp);
  dbl d2L_deta2 = (dvec2_norm_sq(dxy) - dL_deta)/L;

  dvec2 dlp_deta = dvec2_cproj(lp, dvec2_dbl_div(dxy, -L));
  dvec2 d2lp_deta2 = dvec2_add(
    dvec2_dbl_mul(lp, -d2L_deta2),
    dvec2_dbl_mul(dlp_deta, -2*dL_deta)
  );

  dbl d2T_deta2 = cubic_d2f(&ctx->T, eta);

  dvec2 t0 = {
    .x = cubic_f(&ctx->Tx, eta),
    .y = cubic_f(&ctx->Ty, eta)
  };

  dvec2 g = {
    .x = cubic_f(&ctx->Tx, eta),
    .y = cubic_f(&ctx->Ty, eta)
  };
  dvec2 dg_deta = {
    .x = cubic_df(&ctx->Tx, eta),
    .y = cubic_df(&ctx->Ty, eta)
  };
  dvec2 d2g_deta2 = {
    .x = cubic_d2f(&ctx->Tx, eta),
    .y = cubic_d2f(&ctx->Ty, eta)
  };

  dbl gnorm = dvec2_norm(g);
  dbl dgnorm_deta = dvec2_dot(g, dg_deta)/gnorm;
  dbl d2gnorm_deta2 = (
    dvec2_norm_sq(dg_deta) + dvec2_dot(g, d2g_deta2) -
    dgnorm_deta*dgnorm_deta)/gnorm;

  dvec2 dt0_deta = dvec2_dbl_div(
    dvec2_sub(
      dg_deta,
      dvec2_dbl_mul(t0, dgnorm_deta)
    ),
    gnorm
  );
  dvec2 d2t0_deta2 = dvec2_dbl_div(
    dvec2_sub(
      d2g_deta2,
      dvec2_add(
        dvec2_dbl_mul(dt0_deta, 2*dgnorm_deta),
        dvec2_dbl_mul(t0, d2gnorm_deta2)
      )
    ),
    gnorm
  );

  dvec2 dt_deta = dvec2_sub(
    dvec2_dbl_mul(dlp_deta, 1.5),
    dvec2_dbl_mul(dt0_deta, 0.25)
  );
  dvec2 d2t_deta2 = dvec2_sub(
    dvec2_dbl_mul(d2lp_deta2, 1.5),
    dvec2_dbl_mul(d2t0_deta2, 0.25)
  );

  dvec2 t1 = {.x = cos(th), .y = sin(th)};

  dvec2 t = dvec2_sub(
    dvec2_dbl_mul(lp, 1.5),
    dvec2_dbl_mul(dvec2_add(t0, t1), 0.25)
  );
  dbl tnorm = dvec2_norm(t);
  dbl dtnorm_deta = dvec2_dot(t, dt_deta)/tnorm;

  dbl S = (1.0 + 2.0*tnorm)/3.0;
  dbl dS_deta = 2*dtnorm_deta/3;

  dvec2 dt1_dth = {.x = -sin(th), .y = cos(th)};
  dvec2 d2t1_dth2 = {.x = -cos(th), .y = -sin(th)};

  dvec2 dt_dth = dvec2_dbl_mul(dt1_dth, -0.25);
  dvec2 dt2_dth2 = dvec2_dbl_mul(d2t1_dth2, -0.25);

  dbl dS_dth = (2.0/3.0)*dvec2_dot(dt_dth, t);

  dbl d2S_deta2 = (
    2*(dvec2_norm_sq(dt_deta) + dvec2_dot(t, d2t_deta2))/3 -
    dtnorm_deta*dS_deta)/tnorm;

  dbl d2S_deta_dth = (2.0*dvec2_dot(dt_deta, dt_dth))/(3*dvec2_dot(dt_deta, t));

  dbl d2S_dth2 =
    2*(dvec2_norm_sq(dt_dth) + dvec2_dot(t, dt2_dth2))/3 -
    3*dS_dth*dS_dth/2;

  dmat22 hess;
  hess.data[0][0] = d2T_deta2 + d2L_deta2*S + 2*dL_deta*dS_deta + L*d2S_deta2;
  hess.data[1][0] = hess.data[0][1] = dL_deta*dS_dth + L*d2S_deta_dth;
  hess.data[1][1] = L*d2S_dth2;

  return hess;
}
