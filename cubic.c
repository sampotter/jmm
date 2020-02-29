#include "cubic.h"

#include <string.h>

#include "mat.h"

static dmat44 V_inv = {
  .data = {
    { 1,  0,  0,  0},
    { 0,  0,  1,  0},
    {-3,  3, -2, -1},
    { 2, -2,  1,  1}
  }
};

void cubic_set_data(cubic_s *cubic, dvec4 data) {
  cubic->a = dmat44_dvec4_mul(V_inv, data);
}

void cubic_set_data_from_ptr(cubic_s *cubic, dbl const *data_ptr) {
  dvec4 data;
  memcpy((void *)data.data, (void *)data_ptr, 4*sizeof(dbl));
  cubic_set_data(cubic, data);
}

dbl cubic_f(cubic_s const *cubic, dbl lam) {
  dbl const *a = cubic->a.data;
  return a[0] + lam*(a[1] + lam*(a[2] + lam*a[3]));
}

dbl cubic_df(cubic_s const *cubic, dbl lam) {
  dbl const *a = cubic->a.data;
  return a[1] + lam*(2*a[2] + 3*lam*a[3]);
}
