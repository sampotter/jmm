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
