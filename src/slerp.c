#include "slerp.h"

#include <assert.h>

#include "mat.h"

/* Spherical weighted average of two unit vectors, such that the ratio
 * of the arc length from `p0` to `q` to the arc length from `p0` t
 * *`p1` is `w[0]`, and vice versa. */
void slerp2(dbl const p0[3], dbl const p1[3], dbl const w[2], dbl q[3]) {
  dbl const atol = 1e-14;

  assert(fabs(1 - dbl3_norm(p0)) < atol);
  assert(fabs(1 - dbl3_norm(p1)) < atol);
  assert(dbl2_bary(w));

  dbl p0_dot_p1 = dbl3_dot(p0, p1);

  /* If `p0` and `p1` are nearly the same, return their average
   * now. Note that if `p0_dot_p1 > 1`, the following calculations
   * will yield `nan`s, so this is also a safeguard. */
  if (p0_dot_p1 > 1 - atol) {
    dbl3_cc(p0, p1, 0.5, q);
    return;
  }

  dbl om = acos(p0_dot_p1);
  dbl sin_om = sin(om);
  dbl sin_w0_om = sin(w[0]*om);
  dbl sin_w1_om = sin(w[1]*om);

  for (size_t i = 0; i < 3; ++i)
    q[i] = (sin_w0_om*p0[i] + sin_w1_om*p1[i])/sin_om;

  assert(fabs(1 - dbl3_norm(q)) < atol);
}

/* Compute the arc length between between unit vectors `p` and `q`. */
static dbl arclength(dbl const p[3], dbl const q[3]) {
  dbl dot = dbl3_dot(p, q);

  /* If `q` and `q0` are very close together, `acos` is singular, so
   * we just compute the Euclidean distance instead. */
  if (dot > 1 - 1e-7)
    return dbl3_dist(p, q);
  else if (dot < -1 + 1e-7)
    assert(false);
  else
    return acos(dot);
}

static dbl slerp3_f(dbl const q[3], dbl const p[3][3], dbl const w[3]) {
  dbl f = 0, d;
  for (size_t i = 0; i < 3; ++i) {
    d = arclength(q, p[i]);
    f += w[i]*d*d;
  }
  return f/2;
}

/* Compute the spherical weighted average of three unit vectors, up to
 * tolerance `tol`, where `tol` measures the arc length distance
 * between `q` and the true average. The rows of `p` are the unit
 * vectors to average, `w` is a vector of barycentric coordinates, and
 * the result is written to `q`. */
void slerp3(dbl const p[3][3], dbl const w[3], dbl q[3], dbl tol) {
  dbl const atol = 1e-14;

  assert(dbl3_valid_bary_coord(w));

  /* If one of the `w[i]` is nearly 1, then we just copy `p[i]` to `q`
   * and return early. */
  for (size_t i = 0; i < 3; ++i) {
    if (w[i] > 1 - atol) {
      dbl3_copy(p[i], q);
      return;
    }
  }

  for (size_t i = 0; i < 3; ++i)
    assert(fabs(1 - dbl3_norm(p[i])) < atol);

  /* Initialize iterate is the normalized weighted combination of the
   * `p[i]`s. */
  dbl q0[3];
  dbl3_dbl33_mul(w, p, q0);
  dbl3_normalize(q0);

  dbl f0 = slerp3_f(q0, p, w), f;

  size_t niter = 0, max_niter = 20;

  /* Fixed point iteration to find `q`: */
  dbl alpha[3], qt[3], P[3][3], dist, t;
  do {
    assert(niter++ < max_niter);

    /* Compute the normalize factor used to project each `p[i]`. */
    for (size_t i = 0; i < 3; ++i)
      alpha[i] = dbl3_dot(p[i], q0);

    /* Check that `p[i]` lies in the hemisphere determined by `q0`. */
    for (size_t i = 0; i < 3; ++i)
      assert(alpha[i] > 0);

    /* Map each `p[i]` to the tangent plane of S^2 at `q0`. */
    for (size_t i = 0; i < 3; ++i)
      dbl3_dbl_div(p[i], alpha[i], P[i]);

    /* Set `q` to be the weighted average of `P` computed in the tangent
     * plane, then retract to S^2 (i.e., just normalize `q`). */
    dbl3_dbl33_mul(w, P, q);
    dbl3_normalize(q);

    /* Do backtracking line search. */
    t = 1;
    f = slerp3_f(q, p, w);
    while (f > f0) {
      t /= 2;
      slerp2(q0, q, (dbl[2]) {1 - t, t}, qt);
      if (t < 1e-4)
        break;
      f = slerp3_f(qt, p, w);
    }
    if (t < 1)
      dbl3_copy(qt, q);

    /* Use the arc length from `q0` to `q` to check
     * for convergence. */
    dist = arclength(q0, q);

    /* Update `q0` for next iteration */
    dbl3_copy(q, q0);
  } while (dist > tol);
}
