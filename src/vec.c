#include "vec.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "stats.h"

bool dbl2_bary(dbl2 const u) {
  dbl const atol = 1e-14;
  return u[0] > -atol && u[1] > -atol && fabs(1 - u[0] - u[1]) < atol;
}

bool dbl2_isfinite(dbl2 const u) {
  return isfinite(u[0]) && isfinite(u[1]);
}

bool dbl2_all_nan(dbl2 const u) {
  return isnan(u[0]) && isnan(u[1]);
}

dbl dbl2_dist(dbl2 const u, dbl2 const v) {
  dbl tmp[2] = {v[0] - u[0], v[1] - u[1]};
  return sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
}

dbl dbl2_dot(dbl2 const u, dbl2 const v) {
  return u[0]*v[0] + u[1]*v[1];
}

dbl dbl2_maxdist(dbl2 const u, dbl2 const v) {
  return fmax(fabs(v[0] - u[0]), fabs(v[1] - u[1]));
}

dbl dbl2_maxnorm(dbl2 const u) {
  return fmax(fabs(u[0]), fabs(u[1]));
}

dbl dbl2_norm(dbl2 const u) {
  return sqrt(u[0]*u[0] + u[1]*u[1]);
}

dbl dbl2_normsq(dbl2 const u) {
  return u[0]*u[0] + u[1]*u[1];
}

dbl dbl2_prod(dbl2 const u) {
  return u[0]*u[1];
}

dbl dbl2_sum(dbl2 const u) {
  return u[0] + u[1];
}

void dbl2_add(dbl2 const u, dbl2 const v, dbl2 w) {
  w[0] = u[0] + v[0];
  w[1] = u[1] + v[1];
}

void dbl2_avg(dbl2 const u, dbl2 const v, dbl2 w) {
  w[0] = (u[0] + v[0])/2;
  w[1] = (u[1] + v[1])/2;
}

void dbl2_copy(dbl2 const u, dbl2 v) {
  v[0] = u[0];
  v[1] = u[1];
}

void dbl2_cproj(dbl2 const u, dbl2 const v, dbl2 w) {
  w[0] = (1 - u[0]*u[0])*v[0] - u[0]*u[1]*v[1];
  w[1] = -u[0]*u[1]*v[0] + (1 - u[1]*u[1])*v[1];
}

void dbl2_dbl_div(dbl2 const u, dbl a, dbl2 v) {
  v[0] = u[0]/a;
  v[1] = u[1]/a;
}

void dbl2_dbl_div_inplace(dbl2 u, dbl a) {
  u[0] /= a;
  u[1] /= a;
}

void dbl2_dbl_mul(dbl2 const u, dbl a, dbl2 v) {
  v[0] = a*u[0];
  v[1] = a*u[1];
}

void dbl2_floor(dbl2 const u, dbl2 v) {
  v[0] = floor(u[0]);
  v[1] = floor(u[1]);
}

void dbl2_negate(dbl2 u) {
  u[0] = -u[0];
  u[1] = -u[1];
}

void dbl2_normalize(dbl2 u) {
  dbl u_norm = dbl2_norm(u);
  u[0] /= u_norm;
  u[1] /= u_norm;
}

void dbl2_normalized(dbl2 const u, dbl2 v) {
  dbl2_copy(u, v);
  dbl2_normalize(v);
}

void dbl2_saxpy(dbl a, dbl2 const x, dbl2 const y, dbl2 z) {
  z[0] = a*x[0] + y[0];
  z[1] = a*x[1] + y[1];
}

void dbl2_sub(dbl2 const u, dbl2 const v, dbl2 w) {
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
}

void dbl2_sub_inplace(dbl2 u, dbl2 const v) {
  u[0] -= v[0];
  u[1] -= v[1];
}

void dbl2_zero(dbl2 u) {
  u[0] = u[1] = 0;
}

void dbl2_lincomb(dbl a, dbl2 const u, dbl b, dbl2 const v, dbl2 w) {
  w[0] = a*u[0] + b*v[0];
  w[1] = a*u[1] + b*v[1];
}

bool dbl3_equal(dbl3 const x, dbl3 const y) {
  return x[0] == y[0] && x[1] == y[1] && x[2] == y[2];
}

bool dbl3_is_normalized(dbl3 const u) {
  return fabs(1 - dbl3_norm(u)) < 1e-13;
}

bool dbl3_is_zero(dbl3 const u) {
  return u[0] == 0 && u[1] == 0 && u[2] == 0;
}

bool dbl3_isfinite(dbl3 const x) {
  return isfinite(x[0]) && isfinite(x[1]) && isfinite(x[2]);
}

bool dbl3_nonneg(dbl3 const x) {
  return x[0] >= 0 && x[1] >= 0 && x[2] >= 0;
}

bool dbl3_all_nan(dbl3 const x) {
  return isnan(x[0]) && isnan(x[1]) && isnan(x[2]);
}

bool dbl3_valid_bary_coord(dbl3 const b) {
  dbl const atol = 1e-14;
  return b[0] > -atol && b[1] > -atol && b[2] > -atol
    && fabs(1 - dbl3_sum(b)) < atol;
}

dbl dbl3_dist(dbl3 const u, dbl3 const v) {
  dbl tmp[3] = {v[0] - u[0], v[1] - u[1], v[2] - u[2]};
  return sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
}

dbl dbl3_dist_sq(dbl3 const u, dbl3 const v) {
  dbl3 w; dbl3_sub(u, v, w);
  return dbl3_normsq(w);
}

dbl dbl3_dot(dbl3 const u, dbl3 const v) {
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

dbl dbl3_maxdist(dbl3 const u, dbl3 const v) {
  return fmax(fabs(u[0] - v[0]),
              fmax(fabs(u[1] - v[1]),
                   fabs(u[2] - v[2])));
}

dbl dbl3_maxnorm(dbl3 const u) {
  return fmax(fabs(u[0]),
              fmax(fabs(u[1]),
                   fabs(u[2])));
}

dbl dbl3_minimum(dbl3 const u) {
  return fmin(u[0], fmin(u[1], u[2]));
}

dbl dbl3_nanmin(dbl3 const u) {
  dbl nanmin = INFINITY;
  for (size_t i = 0; i < 3; ++i)
    if (!isnan(u[i]))
      nanmin = fmin(nanmin, u[i]);
  return isfinite(nanmin) ? nanmin : NAN;
}

dbl dbl3_nanmax(dbl3 const u) {
  dbl nanmax = -INFINITY;
  for (size_t i = 0; i < 3; ++i)
    if (!isnan(u[i]))
      nanmax = fmax(nanmax, u[i]);
  return isfinite(nanmax) ? nanmax : NAN;
}

dbl dbl3_nanmean(dbl3 const u) {
  size_t n = 0;
  dbl sum;
  for (size_t i = 0; i < 3; ++i) {
    if (!isnan(u[i])) {
      sum += u[i];
      ++n;
    }
  }
  return sum/n;
}

dbl dbl3_ndot(dbl3 const u, dbl3 const v) {
  dbl uv[3] = {u[0]*v[0], u[1]*v[1], u[2]*v[2]};
  return dbl3_nsum(uv);
}

dbl dbl3_norm(dbl3 const u) {
  return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

dbl dbl3_normalize(dbl3 u) {
  dbl unorm = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  assert(unorm != 0);
  u[0] /= unorm;
  u[1] /= unorm;
  u[2] /= unorm;
  return unorm;
}

dbl dbl3_normsq(dbl3 const u) {
  return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

dbl dbl3_nsum(dbl3 const u) {
  volatile dbl sum = 0, c = 0, t, z;

  z = u[0];
  t = sum + z;
  c += fabs(sum) < fabs(z) ? (z - t) + sum : (sum - t) + z;
  sum = t;

  z = u[1];
  t = sum + z;
  c += fabs(sum) < fabs(z) ? (z - t) + sum : (sum - t) + z;
  sum = t;

  z = u[2];
  t = sum + z;
  c += fabs(sum) < fabs(z) ? (z - t) + sum : (sum - t) + z;
  sum = t;

  return sum + c;
}

dbl dbl3_sum(dbl3 const u) {
  return u[0] + u[1] + u[2];
}

dbl dbl3_wnormsq(dbl33 const A, dbl3 const x) {
  dbl y[3];
  for (size_t i = 0; i < 3; ++i)
    y[i] = dbl3_dot(A[i], x);
  return dbl3_dot(y, x);
}

size_t dbl3_argmax(dbl3 const u) {
  dbl umax = -INFINITY;
  size_t argmax;
  for (size_t i = 0; i < 3; ++i) {
    if (u[i] > umax) {
      umax = u[i];
      argmax = i;
    }
  }
  return argmax;
}

void dbl3_abs(dbl3 const u, dbl3 v) {
  v[0] = fabs(u[0]);
  v[1] = fabs(u[1]);
  v[2] = fabs(u[2]);
}

void dbl3_argsort(dbl3 const u, size_t perm[3]) {
  if (u[0] < u[1] && u[0] < u[2]) {
    perm[0] = 0;
    perm[1] = 1;
    perm[2] = 2;
    if (u[1] > u[2]) SWAP(perm[1], perm[2]);
  } else if (u[1] < u[0] && u[1] < u[2]) {
    perm[0] = 1;
    perm[1] = 0;
    perm[2] = 2;
    if (u[0] > u[2]) SWAP(perm[1], perm[2]);
  } else {
    perm[0] = 2;
    perm[1] = 0;
    perm[2] = 1;
    if (u[0] > u[1]) SWAP(perm[1], perm[2]);
  }
}

void dbl3_add(dbl3 const u, dbl3 const v, dbl3 w) {
  w[0] = u[0] + v[0];
  w[1] = u[1] + v[1];
  w[2] = u[2] + v[2];
}

void dbl3_add_inplace(dbl3 u, dbl3 const v) {
  u[0] += v[0];
  u[1] += v[1];
  u[2] += v[2];
}

void dbl3_avg(dbl3 const u, dbl3 const v, dbl3 w) {
  w[0] = (u[0] + v[0])/2;
  w[1] = (u[1] + v[1])/2;
  w[2] = (u[2] + v[2])/2;
}

void dbl3_cc(dbl3 const u0, dbl3 const u1, dbl t0, dbl3 ut) {
  dbl t1 = 1 - t0;
  ut[0] = t0*u0[0] + t1*u1[0];
  ut[1] = t0*u0[1] + t1*u1[1];
  ut[2] = t0*u0[2] + t1*u1[2];
}

void dbl3_copy(dbl3 const u, dbl v[3]) {
  v[0] = u[0];
  v[1] = u[1];
  v[2] = u[2];
}

void dbl3_cross(dbl3 const u, dbl3 const v, dbl3 w) {
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

void dbl3_dbl_div(dbl3 const u, dbl a, dbl3 v) {
  v[0] = u[0]/a;
  v[1] = u[1]/a;
  v[2] = u[2]/a;
}

void dbl3_dbl_div_inplace(dbl3 u, dbl a) {
  u[0] /= a;
  u[1] /= a;
  u[2] /= a;
}

void dbl3_dbl_mul(dbl3 const u, dbl a, dbl3 v) {
  v[0] = a*u[0];
  v[1] = a*u[1];
  v[2] = a*u[2];
}

void dbl3_dbl_mul_inplace(dbl3 u, dbl a) {
  u[0] *= a;
  u[1] *= a;
  u[2] *= a;
}

/* Set `y` to be a random vector which is orthogonal to `x`. */
void dbl3_get_rand_ortho(dbl3 const x, dbl3 y) {
  do {
    for (size_t i = 0; i < 3; ++i)
      y[i] = 2*drand48() - 1;
    dbl3_normalize(y);
  } while (dbl3_dot(x, y) < 1e-1);

  /* Project `y` into the orthogonal complement of `x` and
   * normalize */
  dbl tmp[3]; dbl3_saxpy(-dbl3_dot(x, y), x, y, tmp);
  dbl3_normalized(tmp, y);
}

void dbl3_inf(dbl3 u) {
  u[0] = u[1] = u[2] = INFINITY;
}

void dbl3_max(dbl3 const u, dbl3 const v, dbl3 w) {
  w[0] = fmax(u[0], v[0]);
  w[1] = fmax(u[1], v[1]);
  w[2] = fmax(u[2], v[2]);
}

void dbl3_min(dbl3 const u, dbl3 const v, dbl3 w) {
  w[0] = fmin(u[0], v[0]);
  w[1] = fmin(u[1], v[1]);
  w[2] = fmin(u[2], v[2]);
}

void dbl3_nan(dbl3 u) {
  u[0] = u[1] = u[2] = NAN;
}

void dbl3_negate(dbl3 u) {
  u[0] = -u[0];
  u[1] = -u[1];
  u[2] = -u[2];
}

void dbl3_neginf(dbl3 u) {
  u[0] = u[1] = u[2] = -INFINITY;
}

void dbl3_normalize1(dbl3 x) {
  dbl xnorm1 = fabs(x[0]) + fabs(x[1]) + fabs(x[2]);
  dbl3_dbl_div_inplace(x, xnorm1);
}

void dbl3_normalized(dbl3 const u, dbl3 v) {
  dbl unorm = dbl3_norm(u);
  v[0] = u[0]/unorm;
  v[1] = u[1]/unorm;
  v[2] = u[2]/unorm;
}

void dbl3_one(dbl3 u) {
  u[0] = u[1] = u[2] = 1;
}

void dbl3_gather(dbl const *x, uint3 J, dbl3 xJ) {
  for (size_t i = 0; i < 3; ++i)
    if (J[i] != (size_t)NO_INDEX)
      xJ[i] = x[J[i]];
}

void dbl3_saxpy(dbl a, dbl3 const x, dbl3 const y, dbl3 z) {
  z[0] = a*x[0] + y[0];
  z[1] = a*x[1] + y[1];
  z[2] = a*x[2] + y[2];
}

void dbl3_saxpy_inplace(dbl a, dbl3 const x, dbl3 y) {
  y[0] += a*x[0];
  y[1] += a*x[1];
  y[2] += a*x[2];
}

void dbl3_sort(dbl3 u) {
  dbl tmp;

  // Sort `u` using a sorting network of size 3.

  tmp = fmin(u[0], u[1]);
  u[1] = fmax(u[0], u[1]);
  u[0] = tmp;

  tmp = fmin(u[0], u[2]);
  u[2] = fmax(u[0], u[2]);
  u[0] = tmp;

  tmp = fmin(u[1], u[2]);
  u[2] = fmax(u[1], u[2]);
  u[1] = tmp;
}

void dbl3_sub(dbl3 const u, dbl3 const v, dbl3 w) {
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
  w[2] = u[2] - v[2];
}

void dbl3_sub_inplace(dbl3 u, dbl3 const v) {
  u[0] -= v[0];
  u[1] -= v[1];
  u[2] -= v[2];
}

void dbl3_zero(dbl3 u) {
  u[0] = u[1] = u[2] = 0;
}

bool dbl4_is_rgba(dbl4 const u) {
  return u[0] >= 0 && u[1] >= 0 && u[2] >= 0 && u[3] >= 0
    && u[0] <= 1 && u[1] <= 1 && u[2] <= 1 && u[3] <= 1;
}

bool dbl4_nonneg(dbl4 const u) {
  return u[0] >= 0 && u[1] >= 0 && u[2] >= 0 && u[3] >= 0;
}

bool dbl4_valid_bary_coord(dbl4 const b) {
  dbl const atol = 1e-14;
  return b[0] > -atol && b[1] > -atol && b[2] > -atol && b[3] > -atol
    && fabs(1 - dbl4_sum(b)) < atol;
}

dbl dbl4_dist(dbl4 const u, dbl4 const v) {
  dbl tmp[4] = {u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]};
  return sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2] + tmp[3]*tmp[3]);
}

dbl dbl4_dot(dbl4 const u, dbl4 const v) {
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3];
}

dbl dbl4_norm(dbl4 const u) {
  return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
}

dbl dbl4_nsum(dbl4 const u) {
  // TODO: implement fast inline version of this
  return dblN_nsum(u, 4);
}

dbl dbl4_sum(dbl4 const u) {
  return u[0] + u[1] + u[2] + u[3];
}

void dbl4_add(dbl4 const u, dbl4 const v, dbl4 w) {
  w[0] = u[0] + v[0];
  w[1] = u[1] + v[1];
  w[2] = u[2] + v[2];
  w[3] = u[3] + v[3];
}

void dbl4_copy(dbl4 const u, dbl4 v) {
  v[0] = u[0];
  v[1] = u[1];
  v[2] = u[2];
  v[3] = u[3];
}

void dbl4_d2m(dbl x, dbl4 d2m) {
  d2m[0] = 0.0;
  d2m[1] = 0.0;
  d2m[2] = 2.0;
  d2m[3] = 6.0*x;
}

void dbl4_dbl_div(dbl4 const u, dbl a, dbl4 v) {
  v[0] = u[0]/a;
  v[1] = u[1]/a;
  v[2] = u[2]/a;
  v[3] = u[3]/a;
}

void dbl4_dbl_div_inplace(dbl4 u, dbl a) {
  u[0] /= a;
  u[1] /= a;
  u[2] /= a;
  u[3] /= a;
}

void dbl4_dm(dbl x, dbl4 dm) {
  dm[0] = 0;
  dm[1] = 1;
  dm[2] = 2*x;
  dm[3] = 3*x*x;
}

void dbl4_e1(dbl4 e1) {
  e1[0] = 1;
  e1[1] = 0;
  e1[2] = 0;
  e1[3] = 0;
}

void dbl4_e(dbl4 e, size_t i) {
  for (size_t j = 0; j < 4; ++j)
    e[j] = i == j ? 1 : 0;
}

void dbl4_iota(dbl4 iota) {
  iota[0] = 0;
  iota[1] = 1;
  iota[2] = 2;
  iota[3] = 3;
}

void dbl4_m(dbl x, dbl4 m) {
  m[0] = 1;
  m[1] = x;
  m[2] = x*x;
  m[3] = x*x*x;
}

void dbl4_normalize1(dbl4 u) {
  dbl4_dbl_div_inplace(u, fabs(dbl4_nsum(u)));
}

void dbl4_one(dbl4 one) {
  one[0] = 1;
  one[1] = 1;
  one[2] = 1;
  one[3] = 1;
}

void dbl4_saxpy(dbl a, dbl4 const x, dbl4 const y, dbl4 z) {
  z[0] = a*x[0] + y[0];
  z[1] = a*x[1] + y[1];
  z[2] = a*x[2] + y[2];
  z[3] = a*x[3] + y[3];
}

void dbl4_sub(dbl4 const u, dbl4 const v, dbl4 w) {
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
  w[2] = u[2] - v[2];
  w[3] = u[3] - v[3];
}

void dbl4_zero(dbl4 u) {
  u[0] = u[1] = u[2] = u[3] = 0;
}

dbl dblN_mean(dbl const *x, size_t n) {
  dbl mean = 0;
  for (size_t i = 0; i < n; ++i)
    mean += x[i];
  return mean/n;
}

/**
 * Dot product implemented using a Neumaier sum (see `dblN_nsum`).
 */
dbl dblN_ndot(dbl const * restrict x, dbl const * restrict y, size_t n) {
  volatile dbl sum = 0, c = 0, t, z;
  for (size_t i = 0; i < n; ++i) {
    z = x[i]*y[i];
    t = sum + z;
    if (fabs(sum) < fabs(z))
      c += (z - t) + sum;
    else
      c += (sum - t) + z;
    sum = t;
  }
  return sum + c;
}

/**
 * Neumaier sum of N doubles (see
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm).
 */
dbl dblN_nsum(dbl const *x, size_t n) {
  volatile dbl sum = 0, c = 0, t;
  for (size_t i = 0; i < n; ++i) {
    t = sum + x[i];
    if (fabs(sum) < fabs(x[i]))
      c += (x[i] - t) + sum;
    else
      c += (sum - t) + x[i];
    sum = t;
  }
  return sum + c;
}

void dblN_minmax(dbl const *x, size_t n, dbl *min, dbl *max) {
  *min = INFINITY;
  *max = -INFINITY;
  for (size_t i = 0; i < n; ++i) {
    *min = fmin(*min, x[i]);
    *max = fmax(*max, x[i]);
  }
}

int dbl_compar(dbl const *a, dbl const *b) {
  if (*a < *b) {
    return -1;
  } else if (*a > *b) {
    return 1;
  } else {
    return 0;
  }
}

dbl dblN_median(size_t n, dbl const *x) {
  dbl *x_sorted = malloc(n*sizeof(dbl));
  memcpy(x_sorted, x, n*sizeof(dbl));
  qsort(x_sorted, n, sizeof(dbl), (compar_t)dbl_compar);
  dbl median = n % 2 == 0 ?
    (x_sorted[n/2 - 1] + x_sorted[n/2])/2 :
    x_sorted[n/2];
  free(x_sorted);
  return median;
}

dbl dblN_binmedian(size_t n, dbl const *x, size_t num_bins) {
  // Compute the mean and standard deviation of the centroids
  // components along the split direction.
  runstd_s runstd;
  runstd_init(&runstd);
  for (size_t i = 0; i < n; ++i)
    runstd_update(&runstd, x[i]);
  dbl mu = runstd_get_mean(&runstd);
  dbl sigma = runstd_get_std(&runstd);

  // Compute a rough approximation of the median by binning (based on
  // the ideas in R. Tibshirani's "Fast computation of the median by
  // successive binning").
  dbl binwidth = 2*sigma/num_bins;
  size_t bincount = 0;
  size_t *bins = malloc(num_bins*sizeof(size_t));
  memset(bins, 0x0, num_bins*sizeof(size_t));
  for (size_t i = 0; i < n; ++i) {
    dbl c = x[i] - mu;
    if (c < -sigma || c >= sigma)
      continue;
    int k = floor((c + sigma)/binwidth);
    if (k < 0 || k >= (int)num_bins)
      continue;
    ++bins[k];
    ++bincount;
  }
  dbl binmedian;
  size_t cumsum = 0;
  for (size_t k = 0; k < num_bins; ++k) {
    cumsum += bins[k];
    if (cumsum >= bincount/2) {
      binmedian = mu - sigma + binwidth*(k + 0.5);
      break;
    }
  }
  free(bins);

  return binmedian;
}

void int2_add(int2 const p, int2 const q, int2 r) {
  r[0] = p[0] + q[0];
  r[1] = p[1] + q[1];
}

void int2_copy(int2 const p, int2 q) {
  q[0] = p[0];
  q[1] = p[1];
}

bool int3_equal(int3 const p, int3 const q) {
  return p[0] == q[0] && p[1] == q[1] && p[2] == q[2];
}

int int3_prod(int3 const p) {
  return p[0]*p[1]*p[2];
}

void int3_add(int3 const p, int3 const q, int3 r) {
  r[0] = p[0] + q[0];
  r[1] = p[1] + q[1];
  r[2] = p[2] + q[2];
}

void int3_dbl_mul(int3 const p, dbl a, dbl3 x) {
  x[0] = a*p[0];
  x[1] = a*p[1];
  x[2] = a*p[2];
}

void int3_int_div(int3 const p, int q, int3 r) {
  r[0] = p[0]/q;
  r[1] = p[1]/q;
  r[2] = p[2]/q;
}

bool uint3_equal(uint3 const i, uint3 const j) {
  return i[0] == j[0] && i[1] == j[1] && i[2] == j[2];
}

bool uint3_contains_uint2(uint3 const i, uint2 const j) {
  return (j[0] == i[0] || j[0] == i[1] || j[0] == i[2])
    && (j[1] == i[0] || j[1] == i[1] || j[1] == i[2]);
}

size_t uint3_diff_uint2(uint3 const i, uint2 const j, uint3 k) {
  size_t n = 0;
  for (size_t a = 0; a < 3; ++a)
    if (i[a] != j[0] && i[a] != j[1])
      k[n++] = i[a];
  return n;
}

size_t uint3_find(uint3 const i, size_t j) {
  for (size_t a = 0; a < 3; ++a)
    if (i[a] == j)
      return a;
  return NO_INDEX;
}

bool uint3_is_sorted(uint3 const i) {
  return i[0] <= i[1] && i[1] <= i[2];
}

void uint3_set(uint3 i, size_t value) {
  i[0] = i[1] = i[2] = value;
}
