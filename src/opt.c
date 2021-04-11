#include "opt.h"

#include <assert.h>

#include "mat.h"

void triqp2_solve(triqp2_s *qp) {
  assert(dbl2_isfinite(qp->b));
  assert(dbl2_isfinite(qp->A[0]));
  assert(dbl2_isfinite(qp->A[1]));
  assert(dbl2_isfinite(qp->x));

  dbl const atol = 1e-15;

  dbl x[2], lam[2]; // Temporary variables used below

  /**
   * First, try to compute global minimizer. If it satifies the
   * constraints, we're done.
   */
  dbl22_dbl2_solve(qp->A, qp->b, qp->x);
  dbl2_negate(qp->x);
  if (qp->x[0] >= 0 && qp->x[1] >= 0 && qp->x[0] + qp->x[1] <= 1) {
    return;
  }

  // TODO: if this implementation works, I can simplify what follows
  // quite a bit by looking at b *first*.

  /**
   * Try to minimize over x[0].
   */
  if (qp->x[1] <= 0) {
    x[0] = -qp->b[0]/qp->A[0][0];
    if (0 <= x[0] && x[0] <= 1) {
      // Check Lagrange multiplier
      if (qp->A[1][0]*x[0] + qp->b[1] >= 0) {
        qp->x[0] = x[0];
        qp->x[1] = 0;
        return;
      }
    }
  }

  /**
   * Try to minimize over x[1].
   */
  if (qp->x[0] <= 0) {
    x[1] = -qp->b[1]/qp->A[1][1];
    if (0 <= x[1] && x[1] <= 1) {
      // Check Lagrange multiplier
      if (qp->A[0][1]*x[1] + qp->b[0] >= 0) {
        qp->x[0] = 0;
        qp->x[1] = x[1];
        return;
      }
    }
  }

  /**
   * Try to minimize over 1 - x[0] - x[1].
   */
  if (dbl2_sum(qp->x) >= 1) {
    dbl s = qp->A[0][0] - qp->A[1][0] + qp->b[0] - qp->b[1];
    s /= qp->A[0][0] - qp->A[1][0] - qp->A[0][1] + qp->A[1][1];
    if (-atol <= s && s <= 1 + atol) {
      s = fmax(0, fmin(1, s));
      x[0] = 1 - s;
      x[1] = s;
      // Check Lagrange multiplier
      if (dbl2_dot(qp->A[0], x) + dbl2_dot(qp->A[1], x)
          + dbl2_sum(qp->b) <= atol) {
        qp->x[0] = x[0];
        qp->x[1] = x[1];
        return;
      }
    }
  }

  /**
   * Check if (0, 0) is optimal.
   */
  if (qp->b[0] >= 0 && qp->b[1] >= 0) {
    qp->x[0] = qp->x[1] = 0;
    return;
  }

  // TODO: we should be able to simplify these both significantly and
  // remove the lambda variable entirely... but let's do that later

  /**
   * Check if (1, 0) is optimal.
   */
  {
    lam[0] = qp->A[0][0] + qp->b[0];
    if (lam[0] <= 0) {
      lam[1] = lam[0] - qp->A[1][0] - qp->b[1];
      if (lam[1] <= 0) {
        qp->x[0] = 1;
        qp->x[1] = 0;
        return;
      }
    }
  }

  /**
   * Check if (0, 1) is optimal.
   */
  {
    lam[0] = qp->A[1][1] + qp->b[1];
    if (lam[0] <= 0) {
      lam[1] = lam[0] - qp->A[0][1] - qp->b[0];
      if (lam[1] <= 0) {
        qp->x[0] = 0;
        qp->x[1] = 1;
        return;
      }
    }
  }

  // We shouldn't reach this point! We want to different conditionals
  // above to exhaust every case and positively identify an argmin for
  // the quadratic.
  assert(false);
}
