#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct jet {
  dbl f;
  dbl2 Df;
  dbl fxy;
} jet2;

bool jet2_is_finite(jet2 const *jet);
bool jet2_is_point_source(jet2 const *jet);

typedef struct {
  dbl f;
  dbl2 Df;
  dbl22 D2f;
} jet22t;

jet22t jet22t_make_empty();
bool jet22t_is_finite(jet22t const *jet);
bool jet22t_is_point_source(jet22t const *jet);
jet22t jet22t_sub(jet22t const *J, jet22t const *K);

typedef struct {
  dbl f;
  dbl3 Df;
} jet3;

jet3 jet3_make_empty();
jet3 jet3_make_point_source(dbl tau);
bool jet3_approx_eq(jet3 const *jet1, jet3 const *jet2, dbl atol);
bool jet3_eq(jet3 const *jet1, jet3 const *jet2);
bool jet3_is_finite(jet3 const *jet);
bool jet3_is_point_source(jet3 const *jet);

#ifdef __cplusplus
}
#endif
