#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

typedef struct cubic {
  /**
   * The coefficients `a` are stored so that `a.data[i]` is the
   * coefficient for the `i`th monomial.
   */
  dvec4 a;
} cubic_s;

cubic_s cubic_from_data(dbl f[2], dbl Df[2]);
void cubic_set_data(cubic_s *cubic, dvec4 data);
void cubic_set_data_from_ptr(cubic_s *cubic, dbl const *data_ptr);
void cubic_reverse_on_unit_interval(cubic_s *cubic);
dbl cubic_f(cubic_s const *cubic, dbl lam);
dbl cubic_df(cubic_s const *cubic, dbl lam);
dbl cubic_d2f(cubic_s const *cubic, dbl lam);
void cubic_make_monic(cubic_s *cubic);
int cubic_get_real_roots(cubic_s const *cubic, dbl lam[3]);
void cubic_add_constant(cubic_s *cubic, dbl f);

#ifdef __cplusplus
}
#endif
