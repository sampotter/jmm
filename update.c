#include "update.h"

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

  dvec2 xyeta = dvec2_ccomb(ctx->xy0, ctx->xy1, eta);
  dvec2 lp = dvec2_sub(ctx->xy, xyeta);
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

  dvec2 xyeta = dvec2_ccomb(ctx->xy0, ctx->xy1, eta);
  dvec2 lp = dvec2_sub(ctx->xy, xyeta);
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

  dvec2 xyeta = dvec2_ccomb(ctx->xy0, ctx->xy1, eta);
  dvec2 lp = dvec2_sub(ctx->xy, xyeta);
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
