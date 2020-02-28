#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "vec.h"

typedef struct cubic {
  dvec4 a;
} cubic_s;

void cubic_set_data(cubic_s *cubic, dvec4 data);
void cubic_set_data_from_ptr(cubic_s *cubic, dbl const *data_ptr);
dbl cubic_f(cubic_s const *cubic, dbl lam);
dbl cubic_df(cubic_s const *cubic, dbl lam);

#ifdef __cplusplus
}
#endif
