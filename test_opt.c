#include <cgreen/cgreen.h>

#include "opt.h"

Describe(triqp2);

BeforeEach(triqp2) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(triqp2) {}

Ensure (triqp2, solve_works) {
  triqp2_s qp;

  /**
   * Use this Hessian for all of the test problems in this set of
   * tests.
   */
  qp.A[0][0] = qp.A[1][1] = 2;
  qp.A[1][0] = qp.A[0][1] = 1;

  /**
   * A test case with an interior point solution...
   */
  qp.b[0] = qp.b[1] = -1;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(1./3));
  assert_that_double(qp.x[1], is_nearly_double(1./3));

  /**
   * ... with a solution satisfying 1 - x[0] - x[1] = 0...
   */
  qp.b[0] = qp.b[1] = -1.5;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(1./2));
  assert_that_double(qp.x[1], is_nearly_double(1./2));

  /**
   * ... with a solution satisfying 1 - x[0] - x[1] <= 0...
   */
  qp.b[0] = qp.b[1] = -3;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(1./2));
  assert_that_double(qp.x[1], is_nearly_double(1./2));

  /**
   * ... some QPs with the solution at vertex (0, 0)...
   */
  qp.b[0] = 2;
  qp.b[1] = 1;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(0));
  qp.b[0] = 1;
  qp.b[1] = 2;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(0));
  qp.b[0] = 2;
  qp.b[1] = 2;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(0));

  /**
   * ... with 0 < x[0] < 1 and x[1] == 1...
   */
  qp.b[0] = -1;
  qp.b[1] = -1./2;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0.5));
  assert_that_double(qp.x[1], is_nearly_double(0));
  qp.b[0] = -1.5;
  qp.b[1] = -0.75;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0.75));
  assert_that_double(qp.x[1], is_nearly_double(0));

  /**
   * ... with x[0] == 0 and 0 < x[1] < 1...
   */
  qp.b[0] = -0.5;
  qp.b[1] = -0.75;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(0.5));
  qp.b[0] = -0.75;
  qp.b[1] = -1.5;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(0.75));

  /**
   * ... at the vertex (1, 0)...
   */
  qp.b[0] = -3;
  qp.b[1] = -1;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(1));
  assert_that_double(qp.x[1], is_nearly_double(0));

  /**
   * ... at the vertex (0, 1)...
   */
  qp.b[0] = -1;
  qp.b[1] = -3;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0));
  assert_that_double(qp.x[1], is_nearly_double(1));

  /**
   * ... and with x[0] > 0, x[1] > 0, and x[0] + x[1] == 1.
   */
  qp.b[0] = -2;
  qp.b[1] = -2.2;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0.4));
  assert_that_double(qp.x[1], is_nearly_double(0.6));
  qp.b[0] = -2.2;
  qp.b[1] = -2;
  triqp2_solve(&qp);
  assert_that_double(qp.x[0], is_nearly_double(0.6));
  assert_that_double(qp.x[1], is_nearly_double(0.4));
}
