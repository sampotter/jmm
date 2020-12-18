#include "opt.h"

#include "mat.h"

void triqp2_solve(triqp2_s *qp) {
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
  qp->x[0] = -qp->A[0][0]/qp->b[0];
  if (0 <= qp->x[0] && qp->x[0] <= 1) {
    // Check Lagrange multiplier
    if (qp->b[1] >= 0) {
      qp->x[1] = 0;
      return;
    }
  }

  /**
   * Try to minimize over x[1].
   */
  qp->x[1] = -qp->A[1][1]/qp->b[1];
  if (0 <= qp->x[1] && qp->x[1] <= 1) {
    // Check Lagrange multiplier
    if (qp->b[0] >= 0) {
      qp->x[0] = 0;
      return;
    }
  }

  /**
   * Try to minimize over 1 - x[0] - x[1].
   */
  qp->x[1] = (qp->A[0][0] - qp->A[1][0] + qp->b[0] - qp->b[1])/2;
  if (0 <= qp->x[1] && qp->x[1] <= 1) {
    // Check Lagrange multiplier
    qp->x[0] = 1 - qp->x[1];
    if (dbl2_dot(qp->A[0], qp->x) + dbl2_dot(qp->A[1], qp->x)
        + dbl2_sum(qp->b) <= 0) {
      return;
    }
  }

  /**
   * Check if (0, 0) is optimal.
   */
  qp->x[0] = qp->x[1] = 0;
  if (qp->b[0] + qp->b[1] >= 0) {
    return;
  }

  /**
   * Check if (1, 0) is optimal.
   */
  qp->x[0] = 1;
  dbl lam0 = qp->A[1][0] + qp->b[1];
  if (lam0 <= 0 && lam0 >= qp->A[0][0] + qp->b[0]) {
    return;
  }

  /**
   * Check if (0, 1) is optimal.
   */
  qp->x[0] = 0;
  qp->x[1] = 1;
  lam0 = qp->A[0][0] + qp->b[0];
  if (lam0 <= 0 && lam0 >= qp->A[1][0] + qp->b[1]) {
    return;
  }
}
