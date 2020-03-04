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

void cubic_reverse_on_unit_interval(cubic_s *cubic) {
  static dmat44 M = {
    .data = {
      {1,  1,  1,  1},
      {0, -1, -2, -3},
      {0,  0,  1,  3},
      {0,  0,  0, -1}
    }
  };
  cubic->a = dmat44_dvec4_mul(M, cubic->a);
}

dbl cubic_f(cubic_s const *cubic, dbl lam) {
  return dvec4_dot(cubic->a, dvec4_m(lam));
}

dbl cubic_df(cubic_s const *cubic, dbl lam) {
  return dvec4_dot(cubic->a, dvec4_dm(lam));
}
