#include "geom.h"

#include <math.h>

#include "mat.h"

void lin_comb_unit_vec_3(dbl const t_in[3][3], dbl const b[3], dbl t_out[3]) {
  for (int i = 0; i < 3; ++i) {
    t_out[i] = 0;
    for (int j = 0; j < 3; ++j) {
      t_out[i] += b[i]*t_in[i][j];
    }
  }
  dbl3_normalize(t_out);

bool points_are_coplanar(dbl const **x) {
  dbl dx[3][3];
  dbl3_sub(x[1], x[0], dx[0]);
  dbl3_sub(x[2], x[0], dx[1]);
  dbl3_sub(x[3], x[2], dx[2]);
  dbl det = dbl33_det(dx);
  return fabs(det) < 1e-15;
}
}
