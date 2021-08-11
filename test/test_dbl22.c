#include <cgreen/cgreen.h>

#include "mat.h"

Describe(dbl22);

BeforeEach(dbl22) {
  double_absolute_tolerance_is(1e-13);
  double_relative_tolerance_is(1e-13);
}

AfterEach(dbl22) {}

Ensure(dbl22, isfinite_works) {
  dbl22 A = {{1, 1}, {1, 1}};
  assert_that(dbl22_isfinite(A));
  A[0][0] = INFINITY;
  assert_false(dbl22_isfinite(A));
  A[0][0] = NAN;
  assert_false(dbl22_isfinite(A));
}

Ensure(dbl22, det_works) {
  dbl22 A = {{4, 3}, {2, -3}};
  dbl det_gt = -18;
  assert_that_double(dbl22_det(A), is_nearly_double(det_gt));
}

Ensure(dbl22, trace_works) {
  dbl22 A = {{9, -1}, {0, 2}};
  dbl trace_gt = 11;
  assert_that_double(dbl22_trace(A), is_nearly_double(trace_gt));
}
