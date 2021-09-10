#include "jet.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "mat.h"

bool jet2_is_finite(jet2 const *jet) {
  return isfinite(jet->f) && dbl2_isfinite(jet->Df) && isfinite(jet->fxy);
}

bool jet2_is_point_source(jet2 const *jet) {
  return isfinite(jet->f) && dbl2_all_nan(jet->Df) && isnan(jet->fxy);
}

jet22t jet22t_make_empty() {
  return (jet22t) {
    .f = INFINITY,
    .Df = {NAN, NAN},
    .D2f = {{NAN, NAN}, {NAN, NAN}}
  };
}

bool jet22t_is_point_source(jet22t const *jet) {
  return isfinite(jet->f) && dbl2_all_nan(jet->Df) && dbl22_all_nan(jet->D2f);
}

jet3 jet3_make_empty() {
  return (jet3) {.f = INFINITY, .Df = {NAN, NAN, NAN}};
}

jet3 jet3_make_point_source(dbl tau) {
  return (jet3) {.f = tau, .Df = {NAN, NAN, NAN}};
}

bool jet3_approx_eq(jet3 const *jet1, jet3 const *jet2, dbl atol) {
  return fabs(jet1->f - jet2->f) <= atol &&
    dbl3_maxdist(jet1->Df, jet2->Df) <= atol;
}

bool jet3_eq(jet3 const *jet1, jet3 const *jet2) {
  return !memcmp((void *)jet1, (void *)jet2, sizeof(jet3));
}

bool jet3_is_finite(jet3 const *jet) {
  return isfinite(jet->f) && dbl3_isfinite(jet->Df);
}

bool jet3_is_point_source(jet3 const *jet) {
  return isfinite(jet->f) && dbl3_all_nan(jet->Df);
}
