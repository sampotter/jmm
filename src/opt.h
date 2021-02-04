#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

/**
 * A structure defining the following quadratic program:
 *
 *     minimize  q(x) = x'*A*x/2 + b'*x + c
 *   subject to  x[0] >= 0
 *               x[1] >= 0
 *               1 - x[0] - x[1] >= 0
 *               A = A'
 *
 * This is a simple inequality constrained quadratic program on a
 * triangle. For this problem, A is assumed to be symmetric so that
 * grad(q) = A*x + b. Note also that since we only find the minimizing
 * argument, we ignore the parameter `c`; hence, it is not provided as
 * a field in the `triqp2` struct.
 */
typedef struct triqp2 {
  dbl b[2], A[2][2];
  dbl x[2];
} triqp2_s;

/**
 * Solve the quadratic program specified by `qp`. The variables `c`,
 * `b`, and `A` should be specified. Afterwards, the minimizing
 * argument will be stored in `x`. The minimum value is not computed.
 */
void triqp2_solve(triqp2_s *qp);

#ifdef __cplusplus
}
#endif
