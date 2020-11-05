#include "bary.h"

#include "mat.h"
#include "vec.h"

#include <stdio.h>

void bary3_coef(dbl const **x, dbl const *y, dbl *c) {
  // TODO: this doesn't seem that efficient. It's probably faster to
  // set up the linear system and solve it directly.

  dmat33 A;
  dbl3_sub(x[0], x[3], A.rows[0].data);
  dbl3_sub(x[1], x[3], A.rows[1].data);
  dbl3_sub(x[2], x[3], A.rows[2].data);

  dbl vol = dmat33_det(&A);

  dvec3 tmp;
  for (int i = 0; i < 3; ++i) {
    tmp = A.rows[i];
    dbl3_sub(y, x[3], A.rows[i].data);
    c[i] = dmat33_det(&A)/vol;
    A.rows[i] = tmp;
  }

  dbl3_sub(x[0], y, A.rows[0].data);
  dbl3_sub(x[1], y, A.rows[1].data);
  dbl3_sub(x[2], y, A.rows[2].data);
  c[3] = dmat33_det(&A)/vol;
}
