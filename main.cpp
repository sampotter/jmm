#include "eik.h"

#include <stdlib.h>
#include <stdio.h>

#define MAX(x, y) x > y ? x : y

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

dbl s(dbl x, dbl y, void *context) {
  (void) x;
  (void) y;
  (void) context;
  return 1.0;
}

dvec2 grad_s(dbl x, dbl y, void *context) {
  (void) x;
  (void) y;
  (void) context;
  return dvec2_zero();
}

int main() {
  eik * scheme;
  eik_alloc(&scheme);

  field2_s slow = {
    .f = s,
    .grad_f = grad_s,
    .context = NULL
  };

  int N = 101;
  int i0 = N/2;
  ivec2 shape = {N, N};
  dvec2 xymin = {-1, -1};
  dbl h = 2.0/(N-1);
  eik_init(scheme, &slow, shape, xymin, h);

  int R = 5;
  for (int i = 0; i < N; ++i) {
    int abs_di = abs(i - i0);
    dbl x = h*i + xymin.x;
    for (int j = 0; j < N; ++j) {
      int abs_dj = abs(j - i0);
      dbl y = h*j + xymin.y;
      int r = MAX(abs_di, abs_dj);
      if (r <= R) {
        ivec2 ind = {i, j};
        jet J = {T(x, y), Tx(x, y), Ty(x, y), Txy(x, y)};
        if (r < R) {
          eik_add_valid(scheme, ind, J);
        } else {
          eik_add_trial(scheme, ind, J);
        }
      }
    }
  }
  eik_build_cells(scheme);

  eik_solve(scheme);

  eik_deinit(scheme);
  eik_dealloc(&scheme);
}
