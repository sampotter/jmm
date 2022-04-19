#include "jet.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "mat.h"

bool jet21p_is_finite(jet21p const *jet) {
  return isfinite(jet->f) && dbl2_isfinite(jet->Df) && isfinite(jet->fxy);
}

bool jet21p_is_point_source(jet21p const *jet) {
  return isfinite(jet->f) && dbl2_all_nan(jet->Df) && isnan(jet->fxy);
}

jet21t jet21t_make_empty() {
  return (jet21t) {
    .f = INFINITY,
    .Df = {NAN, NAN},
    .D2f = {{NAN, NAN}, {NAN, NAN}}
  };
}

bool jet21t_is_point_source(jet21t const *jet) {
  return isfinite(jet->f) && dbl2_all_nan(jet->Df) && dbl22_all_nan(jet->D2f);
}

jet31t jet31t_make_empty() {
  return (jet31t) {.f = INFINITY, .Df = {NAN, NAN, NAN}};
}

jet31t jet31t_make_point_source(dbl tau) {
  return (jet31t) {.f = tau, .Df = {NAN, NAN, NAN}};
}

bool jet31t_approx_eq(jet31t const *jet1, jet31t const *jet2, dbl atol) {
  return fabs(jet1->f - jet2->f) <= atol &&
    dbl3_maxdist(jet1->Df, jet2->Df) <= atol;
}

bool jet31t_eq(jet31t const *jet1, jet31t const *jet2) {
  return !memcmp((void *)jet1, (void *)jet2, sizeof(jet31t));
}

bool jet31t_is_finite(jet31t const *jet) {
  return isfinite(jet->f) && dbl3_isfinite(jet->Df);
}

bool jet31t_is_point_source(jet31t const *jet) {
  return isfinite(jet->f) && dbl3_all_nan(jet->Df);
}

jet32t jet32t_make_empty() {
  return (jet32t) {
    .f = INFINITY,
    .Df = {NAN, NAN, NAN},
    .D2f = {{NAN, NAN, NAN},
            {NAN, NAN, NAN},
            {NAN, NAN, NAN}}
  };
}

jet32t jet32t_make_point_source(dbl tau) {
  return (jet32t) {
    .f = tau,
    .Df = {NAN, NAN, NAN},
    .D2f = {{NAN, NAN, NAN},
            {NAN, NAN, NAN},
            {NAN, NAN, NAN}}
  };
}

bool jet32t_is_finite(jet32t const *jet) {
  return isfinite(jet->f)
    && dbl3_isfinite(jet->Df)
    && dbl33_isfinite(jet->D2f);
}

bool jet32t_is_point_source(jet32t const *jet) {
  return isfinite(jet->f)
    && dbl3_all_nan(jet->Df)
    && dbl33_isnan(jet->D2f);
}

bool jet32t_is_singular(jet32t const *jet) {
  return !isfinite(jet->f) || !dbl3_isfinite(jet->Df)
    || !dbl33_isfinite(jet->D2f);
}

jet31t jet31t_from_jet32t(jet32t jet) {
  return (jet31t) {.f = jet.f, .Df = {jet.Df[0], jet.Df[1], jet.Df[2]}};
}
