#include "jet.h"

#include <assert.h>
#include <math.h>
#include <string.h>

jet3 jet3_make_empty() {
  return (jet3) {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN};
}

jet3 jet3_make_point_source(dbl tau) {
  return (jet3) {.f = tau, .fx = NAN, .fy = NAN, .fz = NAN};
}

bool jet3_approx_eq(jet3 const *jet1, jet3 const *jet2, dbl atol) {
  return fabs(jet1->f - jet2->f) <= atol &&
    fabs(jet1->fx - jet2->fx) <= atol &&
    fabs(jet1->fy - jet2->fy) <= atol &&
    fabs(jet1->fz - jet2->fz) <= atol;
}

bool jet3_eq(jet3 const *jet1, jet3 const *jet2) {
  return !memcmp((void *)jet1, (void *)jet2, sizeof(jet3));
}

bool jet3_is_finite(jet3 const *jet) {
  return isfinite(jet->f)
    && isfinite(jet->fx) && isfinite(jet->fy) && isfinite(jet->fz);
}

bool jet3_is_nan(jet3 const *jet) {
  bool is_nan = isnan(jet->fx);
  if (is_nan) {
    assert(isnan(jet->fy)); // Consistency check... If one of the
    assert(isnan(jet->fz)); // partials is NaN, all of them should be.
  }
  return is_nan;
}

bool jet3_is_point_source(jet3 const *jet) {
  return isfinite(jet->f) && jet3_is_nan(jet);
}
