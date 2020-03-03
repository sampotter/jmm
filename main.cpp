#include "sjs_eik.h"

#include <stdlib.h>

dbl T(dbl x, dbl y) {
  return sqrt(x*x + y*y);
}

dbl Tx(dbl x, dbl y) {
  return x/T(x, y);
}

dbl Ty(dbl x, dbl y) {
  return y/T(x, y);
}

dbl Txy(dbl x, dbl y) {
  return -Tx(x, y)*Ty(x, y)/T(x, y);
}

int main() {
  sjs * scheme;
  sjs_alloc(&scheme);

  int N = 101;
  int i0 = N/2;
  ivec2 shape = {N, N};
  dvec2 xymin = {-1, -1};
  dbl h = 2/(N-1);
  sjs_init(scheme, shape, xymin, h);

  int R = 5;
  for (int i = 0; i < N; ++i) {
    int abs_di = abs(i - i0);
    dbl x = h*i + xymin.x;
    for (int j = 0; j < N; ++j) {
      int abs_dj = abs(j - i0);
      dbl y = h*j + xymin.y;
      int r = abs_di + abs_dj;
      if (r <= R) {
        ivec2 ind = {i, j};
        jet J = {T(x, y), Tx(x, y), Ty(x, y), Txy(x, y)};
        if (r < R) {
          sjs_add_valid(scheme, ind, J);
        } else {
          sjs_add_trial(scheme, ind, J);
        }
      }
    }
  }
  sjs_build_cells(scheme);

  sjs_solve(scheme);

  sjs_deinit(scheme);
  sjs_dealloc(&scheme);
}
