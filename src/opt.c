#include <jmm/opt.h>

#include <assert.h>

#include <jmm/mat.h>
#include <jmm/util.h>

#include "macros.h"

void triqp2_solve(triqp2_s *qp, dbl tol) {
  dbl2 const *A = qp->A;
  dbl const *b = &qp->b[0];
  dbl *x = &qp->x[0];

  assert(dbl2_isfinite(b));
  assert(dbl2_isfinite(A[0]));
  assert(dbl2_isfinite(A[1]));

  bool quadratic_part_is_zero =
    fabs(A[0][0]) < tol && fabs(A[0][1]) < tol &&
    fabs(A[1][0]) < tol && fabs(A[1][1]) < tol;

  bool linear_part_is_zero = fabs(b[0]) < tol && fabs(b[1]) < tol;

  /* If both the quadratic and linear parts are zero, we just set the
   * optimum to be the centroid and return... this may not be the most
   * sensible thing to do. */
  if (quadratic_part_is_zero && linear_part_is_zero) {
    x[0] = 1./3;
    x[1] = 1./3;
    return;
  }

  /* If the quadratic part is zero, return the optimum value at the
   * vertex for the linear part */
  if (quadratic_part_is_zero) {
    dbl p00 = 0, p10 = b[0], p01 = b[1];
    if (p00 < p10 && p00 < p01) {
      x[0] = 0;
      x[1] = 0;
    } else if (p10 < p00 && p10 < p01) {
      x[0] = 1;
      x[1] = 0;
    } else if (p01 < p00 && p01 < p10) {
      x[0] = 0;
      x[1] = 1;
    } else if (p00 == p10 && p10 == p01) {
      x[0] = x[1] = 1./3;
    } else if (p00 == p01) {
      x[0] = 0;
      x[1] = 1./2;
    } else if (p01 == p10) {
      x[0] = 1./2;
      x[1] = 1./2;
    } else if (p10 == p00) {
      x[0] = 1./2;
      x[1] = 1;
    } else {
      assert(false);
    }
    return;
  }

  /* If the linear part is zero, then the optimum is the origin */
  if (linear_part_is_zero) {
    x[0] = 0;
    x[1] = 0;
    return;
  }

  /** Next, try to compute the global minimizer: */

  dbl22_dbl2_solve(A, b, x);
  dbl2_negate(x);
  if (x[0] >= 0 && x[1] >= 0 && x[0] + x[1] <= 1) {
    /* If it's in the interior, we're done. */
    return;
  }

  /** If that failed, minimize the quadratic over each edge: */

  /* We parametrize the edge [1, 0] x {0} using x -> (x, 0) with 0 <=
   * x <= 1. The constrained optimum is x10_opt.
   *
   * The edge {0} x [0, 1] is parametrized y -> (0, y), likewise with
   * 0 <= y <= 1. Optimum is y01_opt.
   *
   * The diagonal edge is parametrized x -> (x, 1 - x), with the
   * optimum x value being x11_opt.
   *
   * We compute these below. We're just minimizing quadratics, so
   * minima are available analytically over RR. Afterwards we project
   * each back to the unit interval [0, 1].
   *
   * Note: when we compute these, we need to check and see if these
   * are linear functions we're trying to minimize. If they are, we
   * need to manually find the endpoint. */

  dbl x10_opt = NAN, y01_opt = NAN, x11_opt = NAN;

  dbl p_x10_opt = NAN, p_y01_opt = NAN, p_x11_opt = NAN;

  if (A[0][0] == 0) {
    x10_opt = b[0] >= 0 ? 1 : 0;
  } else {
    x10_opt = clamp(-b[0]/A[0][0], 0, 1);
  }
  p_x10_opt = (b[0] + A[0][0]*x10_opt/2)*x10_opt;

  if (A[1][1] == 0) {
    y01_opt = b[1] >= 0 ? 1 : 0;
  } else {
    y01_opt = clamp(-b[1]/A[1][1], 0, 1);
  }
  p_y01_opt = (b[1] + A[1][1]*y01_opt/2)*y01_opt;

  /* coefficients in x for polynomial on diagonal
   * (i.e., p(x) = q(x, 1 - x)) */
  dbl p_x11_c[3] = {
    A[1][1]/2 + b[1],                // x^0
    A[0][1] - A[1][1] + b[0] - b[1], // x^1
    A[0][0]/2 - A[0][1] + A[1][1]/2  // x^2
  };

  if (p_x11_c[2] == 0) {
    x11_opt = p_x11_c[1] >= 0 ? 1 : 0;
  } else {
    x11_opt = clamp(-p_x11_c[1]/(2*p_x11_c[2]), 0, 1);
  }
  p_x11_opt = p_x11_c[0] + p_x11_c[1]*x11_opt + p_x11_c[2]*x11_opt*x11_opt;

  if (p_x10_opt < p_y01_opt && p_x10_opt < p_x11_opt) {
    x[0] = x10_opt;
    x[1] = 0;
  }

  else if (p_y01_opt < p_x10_opt && p_y01_opt < p_x11_opt) {
    x[0] = 0;
    x[1] = y01_opt;
  }

  else if (p_x11_opt < p_x10_opt && p_x11_opt < p_y01_opt) {
    x[0] = x11_opt;
    x[1] = 1 - x11_opt;
  }

  else if (p_x10_opt == p_y01_opt) {
    assert(x10_opt == 0 && y01_opt == 0);
    x[0] = 0;
    x[1] = 0;
  }

  else if (p_x10_opt == p_x11_opt) {
    assert(fabs(1 - x10_opt) < tol && x11_opt == 1);
    x[0] = 1;
    x[1] = 0;
  }

  else if (p_y01_opt == p_x11_opt) {
    const double cost1 = fabs(x11_opt) + fabs(y01_opt - 1);
    const double cost2 = fabs(x11_opt - 1) + fabs(y01_opt);
    if (cost1 < cost2) {
      x[0] = 0;
      x[1] = 1;
    } else {
      assert(x10_opt == 0.5);
      x[0] = 0.5;
      x[1] = 0;
    }
  }

  else assert(false);

  return;
}
