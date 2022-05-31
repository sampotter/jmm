#include "util.h"

#include <math.h>
#include <string.h>
#include <time.h>

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

int compar_size_t(size_t const *i, size_t const *j) {
  if (*i < *j) {
    return -1;
  } else if (*i > *j) {
    return 1;
  } else {
    return 0;
  }
}

int signum(dbl x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}

dbl shrink(dbl x, dbl eps) {
  return fmax(0, x - eps) + fmin(0, x + eps);
}

bool contains(void const *arr, size_t len, void const *elt, size_t size) {
  char *ptr = (char *)arr;
  for (size_t i = 0; i < len; ++i) {
    if (!memcmp((void *)(ptr + size*i), elt, size)) {
      return true;
    }
  }
  return false;
}

dbl toc() {
  static clock_t t1 = 0;
  clock_t t0 = t1;
  t1 = clock();
  return ((double)t1 - (double)t0)/CLOCKS_PER_SEC;
}
