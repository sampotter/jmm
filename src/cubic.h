#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

typedef struct cubic {
  dvec4 a;
} cubic_s;

cubic_s cubic_from_data(dbl f[2], dbl Df[2]);
void cubic_set_data(cubic_s *cubic, dvec4 data);
void cubic_set_data_from_ptr(cubic_s *cubic, dbl const *data_ptr);
void cubic_reverse_on_unit_interval(cubic_s *cubic);
dbl cubic_f(cubic_s const *cubic, dbl lam);
dbl cubic_df(cubic_s const *cubic, dbl lam);
dbl cubic_d2f(cubic_s const *cubic, dbl lam);

#ifdef __cplusplus
}
#endif
