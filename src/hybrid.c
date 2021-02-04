#include "hybrid.h"

#include <math.h>

#include "util.h"

dbl hybrid(dbl (*f)(dbl, void *), dbl a, dbl b, void *context) {
  dbl c, d, fa, fb, fc, fd, dm, df, ds, dd, tmp;

  fa = f(a, context);
  if (fabs(fa) <= EPS) {
    return a;
  }

  fb = f(b, context);
  if (fabs(fb) <= EPS) {
    return b;
  }

  if (sgn(fa) == sgn(fb)) {
    return sgn(fa) == 1 ? a : b;
  }

  c = a;
  fc = fa;
  for (;;) {
    if (fabs(fc) < fabs(fb)) {
      tmp = b; b = c; c = tmp;
      tmp = fb; fb = fc; fc = tmp;
      a = c;
      fa = fc;
    }
    if (fabs(b - c) <= EPS) {
      break;
    }
    dm = (c - b)/2;
    df = fa - fb;
    ds = df == 0 ? dm : -fb*(a - b)/df;
    dd = sgn(ds) != sgn(dm) || fabs(ds) > fabs(dm) ? dm : ds;
    if (fabs(dd) < EPS) {
      dd = EPS*sgn(dm)/2;
    }
    d = b + dd;
    fd = f(d, context);
    if (fd == 0) {
      c = d;
      b = c;
      fc = fd;
      fb = fc;
      break;
    }
    a = b;
    b = d;
    fa = fb;
    fb = fd;
    if (sgn(fb) == sgn(fc)) {
      c = a;
      fc = fa;
    }
  }
  return (b + c)/2;
}
