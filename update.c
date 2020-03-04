#include "update.h"

dbl F(dbl lam, void *context) {
  update_data *data = (update_data *)context;
  dbl T = cubic_f(&data->cubic, lam);
  dvec2 xylam = dvec2_ccomb(data->xy0, data->xy1, lam);
  dbl L = dvec2_dist(data->xy, xylam);
  return T + L;
}

dbl dF_dt(dbl lam, void *context) {
  update_data *data = (update_data *)context;
  dbl dT_dt = cubic_df(&data->cubic, lam);
  dvec2 xylam = dvec2_ccomb(data->xy0, data->xy1, lam);
  dvec2 xylam_minus_xy = dvec2_sub(xylam, data->xy);
  dbl L = dvec2_norm(xylam_minus_xy);
  dvec2 dxy = dvec2_sub(data->xy1, data->xy0);
  dbl dL_dt = dvec2_dot(dxy, xylam_minus_xy)/L;
  return dT_dt + dL_dt;
}
