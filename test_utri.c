#include <cgreen/cgreen.h>

#include "utri.h"

Describe(utri);

BeforeEach(utri) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(utri) {}

Ensure (utri, tri11_works) {
  dbl x[3] = {1, 1, 0};
  dbl Xt[2][3] = {{1, 0, 0}, {0, 1, 0}};
  dbl jet_data[2][4] = {{1, 1, 0, 0}, {1, 0, 1, 0}};
  jet3 jet[2];

  utri_s *utri;
  utri_alloc(&utri);

  dbl lam;

  memcpy(&jet[0], jet_data[0], 4*sizeof(dbl));
  memcpy(&jet[1], jet_data[1], 4*sizeof(dbl));

  utri_init(utri, x, Xt, jet);
  assert_that(utri_is_causal(utri));

  utri_solve(utri);
  lam = utri_get_lambda(utri);

  assert_that_double(lam, is_nearly_double(0.5));

  utri_dealloc(&utri);
}
