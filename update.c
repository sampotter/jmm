#include "update.h"

dbl F(dbl lam, void *context) {
  update_data *data = (update_data *)context;
  dbl T = cubic_f(&data->cubic, lam);
  dvec2 xylam = dvec2_ccomb(data->xy0, data->xy1, lam);
  dbl L = dvec2_dist(data->xy, xylam);
  return T + data->h*L;
}

dbl dF_dt(dbl lam, void *context) {
  update_data *data = (update_data *)context;
  dbl dT_dt = cubic_df(&data->cubic, lam);
  dvec2 xylam = dvec2_ccomb(data->xy0, data->xy1, lam);
  dvec2 dxy = dvec2_sub(data->xy, xylam);
  dbl L = dvec2_norm(dxy);
  dbl dL_dt = dvec2_dot(dvec2_sub(data->xy0, data->xy1), dxy)/L;
  return dT_dt + data->h*dL_dt;
}
