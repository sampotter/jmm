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

  int perm[6][3] = {
    {0, 1, 2},
    {0, 2, 1},
    {1, 0, 2},
    {2, 0, 1},
    {1, 2, 0},
    {2, 1, 0}
  };

  dbl f, lam;
  dbl const f_gt = 1.4571067811865475;
  dbl const lam_gt = 0.5;

  dbl3 x_perm, Xt_perm[2];
  dbl jet_data_perm[2][4];

  jet_data_perm[0][0] = jet_data[0][0];
  jet_data_perm[1][0] = jet_data[1][0];

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      x_perm[j] = x[perm[i][j]];
      Xt_perm[0][j] = Xt[0][perm[i][j]];
      Xt_perm[1][j] = Xt[1][perm[i][j]];
      jet_data_perm[0][1 + j] = jet_data[0][1 + perm[i][j]];
      jet_data_perm[1][1 + j] = jet_data[1][1 + perm[i][j]];
    }

    memcpy(&jet[0], jet_data_perm[0], 4*sizeof(dbl));
    memcpy(&jet[1], jet_data_perm[1], 4*sizeof(dbl));

    // TODO: can't for the life of me figure out what the problem here
    // is... seems to be some weird issue with GCC 11.2.1 on Fedora...
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
    utri_spec_s spec = utri_spec_from_raw_data(x_perm, Xt_perm, jet);
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
    utri_init(utri, &spec);
    assert_that(utri_is_causal(utri));

    utri_solve(utri);

    f = utri_get_value(utri);
    assert_that_double(f, is_nearly_double(f_gt));

    lam = utri_get_lambda(utri);
    assert_that_double(lam, is_nearly_double(lam_gt));
  }

  utri_dealloc(&utri);
}
