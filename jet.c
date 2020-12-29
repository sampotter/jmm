#include "jet.h"

#include <math.h>
#include <string.h>

jet3 jet3_make_point_source(dbl tau) {
  return (jet3) {.f = tau, .fx = NAN, .fy = NAN, .fz = NAN};
}

bool jet3_eq(jet3 const *jet1, jet3 const *jet2) {
  return !memcmp((void *)jet1, (void *)jet2, sizeof(jet3));
}

bool jet3_is_finite(jet3 const *jet) {
  return isfinite(jet->f)
    && isfinite(jet->fx) && isfinite(jet->fy) && isfinite(jet->fz);
}
