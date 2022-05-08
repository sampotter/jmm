#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct jet {
  dbl f;
  dbl2 Df;
  dbl fxy;
} jet21p;

bool jet21p_is_finite(jet21p const *jet);
bool jet21p_is_point_source(jet21p const *jet);

typedef struct {
  dbl f;
  dbl2 Df;
  dbl22 D2f;
} jet21t;

jet21t jet21t_make_empty();
bool jet21t_is_point_source(jet21t const *jet);

typedef struct {
  dbl f;
  dbl3 Df;
} jet31t;

jet31t jet31t_make_empty();
jet31t jet31t_make_point_source(dbl tau);
bool jet31t_approx_eq(jet31t const *jet1, jet31t const *jet2, dbl atol);
bool jet31t_eq(jet31t const *jet1, jet31t const *jet2);
bool jet31t_is_finite(jet31t const *jet);
bool jet31t_is_point_source(jet31t const *jet);
void jet31t_sub(jet31t const *in1, jet31t const *in2, jet31t *out);

typedef struct {
  dbl f;
  dbl3 Df;
  dbl33 D2f;
} jet32t;

jet32t jet32t_make_empty();
jet32t jet32t_make_point_source(dbl tau);
bool jet32t_is_finite(jet32t const *jet);
bool jet32t_is_point_source(jet32t const *jet);
bool jet32t_is_singular(jet32t const *jet);

// Conversions:

jet31t jet31t_from_jet32t(jet32t jet);

#ifdef __cplusplus
}
#endif
