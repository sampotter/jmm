#include "util.h"

#include <math.h>

#include "mat.h"
#include "vec.h"

dbl clamp(dbl x, dbl a, dbl b) {
  return fmax(a, fmin(x, b));
}

int sgn(dbl x) {
  if (x > 0) {
    return 1;
  } else if (x < 0) {
    return -1;
  } else {
    return 0;
  }
}

dbl q(dbl const A[3][3], dbl const b[3], dbl c, dbl const x[3]) {
  dbl Ax[3];
  dbl33_dbl3_mul(A, x, Ax);
  return dbl3_dot(Ax, x) + dbl3_dot(b, x) + c;
}

void Dq(dbl const A[3][3], dbl const b[3], dbl const x[3], dbl Dq[3]) {
  dbl At[3][3];
  dbl33_transposed(A, At);

  dbl Ax[3], Atx[3];
  dbl33_dbl3_mul(A, x, Ax);
  dbl33_dbl3_mul(At, x, Atx);

  dbl tmp[3];
  dbl3_add(Ax, Atx, tmp);

  dbl3_add(tmp, b, Dq);
}
