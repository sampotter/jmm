#include "geom.h"

#include "vec.h"

void lin_comb_unit_vec_3(dbl const t_in[3][3], dbl const b[3], dbl t_out[3]) {
  for (int i = 0; i < 3; ++i) {
    t_out[i] = 0;
    for (int j = 0; j < 3; ++j) {
      t_out[i] += b[i]*t_in[i][j];
    }
  }
  dbl3_normalize(t_out);
}
