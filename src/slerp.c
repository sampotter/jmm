#include "slerp.h"

#include <assert.h>

#include "mat.h"
#include "util.h"

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

  assert(dot >= -1 + 1e-7);

  /* If `q` and `q0` are very close together, `acos` is singular, so
   * we just compute the Euclidean distance instead. */
  return dot > 1 - 1e-7 ? dbl3_dist(p, q) : acos(dot);
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

  /* Count and check which coefficients are approximately nonzero */
  bool active[3];
  size_t num_active = 0, active_inds[3];
  for (size_t i = 0; i < 3; ++i)
    if ((active[i] = w[i] > tol))
      active_inds[num_active++] = i;

  if (num_active == 1) {
    dbl3_copy(p[active_inds[0]], q);
    return;
  }

  if (num_active == 2) {
    dbl v[2];
    for (size_t i = 0; i < 2; ++i)
      v[i] = w[active_inds[i]];

    dbl normalize = fabs(v[0]) + fabs(v[1]);
    v[0] /= normalize;
    v[1] /= normalize;

    slerp2(p[active_inds[0]], p[active_inds[1]], v, q);
    return;
  }

  for (size_t i = 0; i < 3; ++i)
    assert(fabs(1 - dbl3_norm(p[i])) < atol);

  /* Initialize iterate is the normalized weighted combination of the
   * `p[i]`s. */
  dbl q0[3];
  dbl3_dbl33_mul(w, p, q0);
  dbl3_normalize(q0);

  dbl f0 = slerp3_f(q0, p, w), f;

#if JMM_DEBUG
  size_t niter = 0, max_niter = 20;
#endif

  /* Fixed point iteration to find `q`: */
  dbl alpha[3], qt[3], P[3][3], dist, t;
  do {
#if JMM_DEBUG
    assert(niter++ < max_niter);
#endif

    /* Compute the normalize factor used to project each `p[i]`. */
    dbl33_dbl3_mul(p, q0, alpha);

    /* Check that none of `p[i]`s are `q0`'s antipodal point. */
    for (size_t i = 0; i < 3; ++i)
      assert(alpha[i] > -1 + atol);

    /* Map each `p[i]` to the tangent plane of S^2 at `q0` (or rather,
     * the tangent plane offset by `q0`). We do this by evaluating the
     * exponential map of S^2 at `q0`. */
    for (size_t i = 0; i < 3; ++i) {
      dbl dot = dbl3_dot(q0, p[i]); // dot = q0'*p[i]
      dbl tmp[3]; dbl3_saxpy(-dot, q0, p[i], tmp);
      dbl norm = dbl3_norm(tmp);
      if (norm > atol) {
        dbl theta = acos(clamp(dot, -1, 1));
        dbl3_saxpy(theta/norm, tmp, q0, P[i]);
      } else {
        /* If `tmp` is approximately 0, just set `P[i]` to `q0` */
        dbl3_copy(q0, P[i]);
      }
    }

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

static dbl slerp4_f(dbl const q[3], dbl const p[4][3], dbl const w[4]) {
  dbl f = 0, d;
  for (size_t i = 0; i < 4; ++i) {
    d = arclength(q, p[i]);
    f += w[i]*d*d;
  }
  return f/2;
}

void slerp4(dbl const p[4][3], dbl const w[4], dbl q[3], dbl tol) {
  dbl const atol = 1e-14;

  assert(dbl4_valid_bary_coord(w));

  /* Count and check which coefficients are approximately nonzero */
  bool active[4];
  size_t num_active = 0, active_inds[4];
  for (size_t i = 0; i < 4; ++i)
    if ((active[i] = w[i] > tol))
      active_inds[num_active++] = i;

  if (num_active == 1) {
    dbl3_copy(p[active_inds[0]], q);
    return;
  }

  if (num_active == 2) {
    dbl v[2];
    for (size_t i = 0; i < 2; ++i)
      v[i] = w[active_inds[i]];

    dbl normalize = fabs(v[0]) + fabs(v[1]);
    v[0] /= normalize;
    v[1] /= normalize;

    slerp2(p[active_inds[0]], p[active_inds[1]], v, q);
    return;
  }

  if (num_active == 3) {
    dbl p_active[3][3];
    for (size_t i = 0; i < 3; ++i)
      dbl3_copy(p[active_inds[i]], p_active[i]);

    dbl w_active[3];
    for (size_t i = 0; i < 3; ++i)
      w_active[i] = w[active_inds[i]];
    dbl3_normalize1(w_active);

    slerp3(p_active, w_active, q, tol);
    return;
  }

  for (size_t i = 0; i < 4; ++i)
    assert(fabs(1 - dbl3_norm(p[i])) < atol);

  /* Initialize iterate is the normalized weighted combination of the
   * `p[i]`s. */
  dbl q0[3];
  dbl4_dbl43_mul(w, p, q0);
  dbl3_normalize(q0);

  dbl f0 = slerp4_f(q0, p, w), f;

#if JMM_DEBUG
  size_t niter = 0, max_niter = 20;
#endif

  /* Fixed point iteration to find `q`: */
  dbl alpha[4], qt[3], P[4][3], dist, t;
  do {
#if JMM_DEBUG
    assert(niter++ < max_niter);
#endif

    /* Compute the normalize factor used to project each `p[i]`. */
    dbl43_dbl3_mul(p, q0, alpha);

    /* Check that none of `p[i]`s are `q0`'s antipodal point. */
    for (size_t i = 0; i < 4; ++i)
      assert(alpha[i] > -1 + atol);

    /* Map each `p[i]` to the tangent plane of S^2 at `q0` (or rather,
     * the tangent plane offset by `q0`). We do this by evaluating the
     * exponential map of S^2 at `q0`. */
    for (size_t i = 0; i < 4; ++i) {
      dbl dot = dbl3_dot(q0, p[i]); // dot = q0'*p[i]
      dbl tmp[3]; dbl3_saxpy(-dot, q0, p[i], tmp);
      dbl norm = dbl3_norm(tmp);
      if (norm > atol) {
        dbl theta = acos(clamp(dot, -1, 1));
        dbl3_saxpy(theta/norm, tmp, q0, P[i]);
      } else {
        /* If `tmp` is approximately 0, just set `P[i]` to `q0` */
        dbl3_copy(q0, P[i]);
      }
    }

    /* Set `q` to be the weighted average of `P` computed in the tangent
     * plane, then retract to S^2 (i.e., just normalize `q`). */
    dbl4_dbl43_mul(w, P, q);
    dbl3_normalize(q);

    /* Do backtracking line search. */
    t = 1;
    f = slerp4_f(q, p, w);
    while (f > f0) {
      t /= 2;
      slerp2(q0, q, (dbl[2]) {1 - t, t}, qt);
      if (t < 1e-4)
        break;
      f = slerp4_f(qt, p, w);
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
