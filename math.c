#include "math.h"

#include <math.h>

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
