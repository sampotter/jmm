#include "vec.h"

#include <assert.h>

#include "macros.h"

void dbl2_add(dbl u[2], dbl v[2], dbl w[2]) {
  w[0] = u[0] + v[0];
  w[1] = u[1] + v[1];
}

void dbl2_sub(dbl const *u, dbl const *v, dbl *w) {
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
}

dbl dbl2_dot(dbl const *u, dbl const *v) {
  return u[0]*v[0] + u[1]*v[1];
}

void dbl2_negate(dbl u[2]) {
  u[0] = -u[0];
  u[1] = -u[1];
}

dbl dbl2_sum(dbl const u[2]) {
  return u[0] + u[1];
}

void dbl2_saxpy(dbl a, dbl const x[2], dbl const y[2], dbl z[2]) {
  z[0] = a*x[0] + y[0];
  z[1] = a*x[1] + y[1];
}

dbl dbl2_dist(dbl const u[2], dbl const v[2]) {
  dbl tmp[2] = {v[0] - u[0], v[1] - u[1]};
  return sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
}

dbl dbl2_norm(dbl const u[2]) {
  return sqrt(u[0]*u[0] + u[1]*u[1]);
}

dbl dbl2_norm_sq(dbl const u[2]) {
  return u[0]*u[0] + u[1]*u[1];
}

dbl dbl2_maxdist(dbl const u[2], dbl const v[2]) {
  return fmax(fabs(v[0] - u[0]), fabs(v[1] - u[1]));
}

dbl dbl2_maxnorm(dbl const u[2]) {
  return fmax(fabs(u[0]), fabs(u[1]));
}

bool dbl2_isfinite(dbl const u[2]) {
  return isfinite(u[0]) && isfinite(u[1]);
}

void dbl3_add(dbl const u[3], dbl const v[3], dbl w[3]) {
  w[0] = u[0] + v[0];
  w[1] = u[1] + v[1];
  w[2] = u[2] + v[2];
}

void dbl3_add_inplace(dbl u[3], dbl const v[3]) {
  u[0] += v[0];
  u[1] += v[1];
  u[2] += v[2];
}

void dbl3_sub(dbl const *u, dbl const *v, dbl *w) {
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
  w[2] = u[2] - v[2];
}

void dbl3_sub_inplace(dbl u[3], dbl const v[3]) {
  u[0] -= v[0];
  u[1] -= v[1];
  u[2] -= v[2];
}

dbl dbl3_dot(dbl const *u, dbl const *v) {
  return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

dbl dbl3_dist(dbl const u[3], dbl const v[3]) {
  dbl tmp[3] = {v[0] - u[0], v[1] - u[1], v[2] - u[2]};
  return sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
}

dbl dbl3_norm(dbl const u[3]) {
  return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

dbl dbl3_normsq(dbl const u[3]) {
  return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

void dbl3_dbl_div(dbl u[3], dbl a, dbl v[3]) {
  v[0] = u[0]/a;
  v[1] = u[1]/a;
  v[2] = u[2]/a;
}

void dbl3_dbl_div_inplace(dbl u[3], dbl a) {
  u[0] /= a;
  u[1] /= a;
  u[2] /= a;
}

dbl dbl3_normalize(dbl u[3]) {
  dbl unorm = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  assert(unorm != 0);
  u[0] /= unorm;
  u[1] /= unorm;
  u[2] /= unorm;
  return unorm;
}

void dbl3_normalized(dbl const u[3], dbl v[3]) {
  dbl unorm = dbl3_norm(u);
  v[0] = u[0]/unorm;
  v[1] = u[1]/unorm;
  v[2] = u[2]/unorm;
}

void dbl3_saxpy(dbl a, dbl const x[3], dbl const y[3], dbl z[3]) {
  z[0] = a*x[0] + y[0];
  z[1] = a*x[1] + y[1];
  z[2] = a*x[2] + y[2];
}

bool dbl3_is_zero(dbl const u[3]) {
  return u[0] == 0 && u[1] == 0 && u[2] == 0;
}

dbl dbl3_maxnorm(dbl const u[3]) {
  return fmax(fabs(u[0]), fmax(fabs(u[1]), fabs(u[2])));
}

void dbl3_cross(dbl const u[3], dbl const v[3], dbl w[3]) {
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
}

dbl dbl3_sum(dbl const u[3]) {
  return u[0] + u[1] + u[2];
}

void dbl3_negate(dbl u[3]) {
  u[0] = -u[0];
  u[1] = -u[1];
  u[2] = -u[2];
}

void dbl3_inf(dbl u[3]) {
  u[0] = u[1] = u[2] = INFINITY;
}

void dbl3_neginf(dbl u[3]) {
  u[0] = u[1] = u[2] = -INFINITY;
}

void dbl3_min(dbl const u[3], dbl const v[3], dbl w[3]) {
  w[0] = fmin(u[0], v[0]);
  w[1] = fmin(u[1], v[1]);
  w[2] = fmin(u[2], v[2]);
}

void dbl3_max(dbl const u[3], dbl const v[3], dbl w[3]) {
  w[0] = fmax(u[0], v[0]);
  w[1] = fmax(u[1], v[1]);
  w[2] = fmax(u[2], v[2]);
}

void dbl3_copy(dbl const u[3], dbl v[3]) {
  v[0] = u[0];
  v[1] = u[1];
  v[2] = u[2];
}

dbl dbl3_ndot(dbl const u[3], dbl const v[3]) {
  return dblN_ndot(u, v, 3); // TODO: inline and optimize me!
}

bool dbl3_valid_bary_coord(dbl const b[3]) {
  dbl const atol = 1e-14;
  return b[0] > -atol && b[1] > -atol && b[2] > -atol
    && fabs(1 - dbl3_sum(b)) < atol;
}

void dbl3_sort(dbl u[3]) {
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

void dbl3_zero(dbl u[3]) {
  u[0] = u[1] = u[2] = 0;
}

void dbl3_cc(dbl const u0[3], dbl const u1[3], dbl t0, dbl ut[3]) {
  dbl t1 = 1 - t0;
  ut[0] = t0*u0[0] + t1*u1[0];
  ut[1] = t0*u0[1] + t1*u1[1];
  ut[2] = t0*u0[2] + t1*u1[2];
}

void dbl3_normalize1(dbl x[3]) {
  dbl xnorm1 = fabs(x[0]) + fabs(x[1]) + fabs(x[2]);
  dbl3_dbl_div_inplace(x, xnorm1);
}

bool dbl3_isfinite(dbl const x[3]) {
  return isfinite(x[0]) && isfinite(x[1]) && isfinite(x[2]);
}

/* Direct, bitwise comparison of `x` and `y`. Returns `true` if all
 * components are equal. */
bool dbl3_equal(dbl const x[3], dbl const y[3]) {
  return x[0] == y[0] && x[1] == y[1] && x[2] == y[2];
}

void dbl4_dbl_div_inplace(dbl u[4], dbl a) {
  u[0] /= a;
  u[1] /= a;
  u[2] /= a;
  u[3] /= a;
}

dbl dbl4_dist(dbl const u[4], dbl const v[4]) {
  dbl tmp[4] = {u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]};
  return sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2] + tmp[3]*tmp[3]);
}

bool dbl4_nonneg(dbl const u[4]) {
  return u[0] >= 0 && u[1] >= 0 && u[2] >= 0 && u[3] >= 0;
}

dbl dbl4_norm(dbl const u[4]) {
  return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
}

void dbl4_saxpy(dbl a, dbl const x[4], dbl const y[4], dbl z[4]) {
  z[0] = a*x[0] + y[0];
  z[1] = a*x[1] + y[1];
  z[2] = a*x[2] + y[2];
  z[3] = a*x[3] + y[3];
}

void dbl4_sub(dbl const u[4], dbl const v[4], dbl w[4]) {
  w[0] = u[0] - v[0];
  w[1] = u[1] - v[1];
  w[2] = u[2] - v[2];
  w[3] = u[3] - v[3];
}

dbl dbl4_sum(dbl const u[4]) {
  return u[0] + u[1] + u[2] + u[3];
}

bool dbl4_valid_bary_coord(dbl const b[4]) {
  dbl const atol = 1e-15;
  return b[0] > -atol && b[1] > -atol && b[2] > -atol && b[3] > -atol
    && fabs(1 - dbl4_sum(b)) < atol;
}

void dbl4_normalize1(dbl u[4]) {
  dbl4_dbl_div_inplace(u, fabs(dbl4_nsum(u)));
}

dbl dbl4_nsum(dbl const u[4]) {
  // TODO: implement fast inline version of this
  return dblN_nsum(u, 4);
}

dbl dblN_mean(dbl const *x, size_t n) {
  dbl mean = 0;
  for (size_t i = 0; i < n; ++i)
    mean += x[i];
  return mean/n;
}

void dblN_minmax(dbl const *x, size_t n, dbl *min, dbl *max) {
  *min = INFINITY;
  *max = -INFINITY;
  for (size_t i = 0; i < n; ++i) {
    *min = fmin(*min, x[i]);
    *max = fmax(*max, x[i]);
  }
}

/**
 * Dot product implemented using a Neumaier sum (see `dblN_nsum`).
 */
dbl dblN_ndot(dbl const *x, dbl const *y, size_t n) {
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

void int3_add(int const p[3], int const q[3], int r[3]) {
  r[0] = p[0] + q[0];
  r[1] = p[1] + q[1];
  r[2] = p[2] + q[2];
}

dvec2 dvec2_zero(void) {
  return (dvec2) {0.0, 0.0};
}

dvec2 dvec2_ccomb(dvec2 v0, dvec2 v1, dbl t) {
  dvec2 vt = {(1 - t)*v0.x + t*v1.x, (1 - t)*v0.y + t*v1.y};
  return vt;
}

dvec2 dvec2_from_ptr(dbl const *u) {
  return (dvec2) {.x = u[0], .y = u[1]};
}

dbl dvec2_maxdist(dvec2 u, dvec2 v) {
  return fmax(fabs(v.x - u.x), fabs(v.y - u.y));
}

dbl dvec2_dist(dvec2 v0, dvec2 v1) {
  dbl dx = v1.x - v0.x, dy = v1.y - v0.y;
  return sqrt(dx*dx + dy*dy);
}

dbl dvec2_maxnorm(dvec2 v) {
  return fmax(fabs(v.x), fabs(v.y));
}

dbl dvec2_norm(dvec2 v) {
  return sqrt(v.x*v.x + v.y*v.y);
}

dbl dvec2_norm_sq(dvec2 v) {
  return v.x*v.x + v.y*v.y;
}

dbl dvec2_dot(dvec2 u, dvec2 v) {
  return u.x*v.x + u.y*v.y;
}

dvec2 dvec2_add(dvec2 u, dvec2 v) {
  dvec2 w = {.x = u.x + v.x, .y = u.y + v.y};
  return w;
}

dvec2 dvec2_sub(dvec2 u, dvec2 v) {
  dvec2 w = {.x = u.x - v.x, .y = u.y - v.y};
  return w;
}

dvec2 dvec2_saxpy(dbl a, dvec2 x, dvec2 y) {
  return (dvec2) {.x = a*x.x + y.x, .y = a*x.y + y.y};
}

dvec2 dvec2_dbl_div(dvec2 v, dbl a) {
  dvec2 w = {.x = v.x/a, .y = v.y/a};
  return w;
}

dvec2 dvec2_dbl_mul(dvec2 v, dbl a) {
  dvec2 w = {.x = a*v.x, .y = a*v.y};
  return w;
}

dvec2 dvec2_floor(dvec2 v) {
  dvec2 w = {.x = floor(v.x), .y = floor(v.y)};
  return w;
}

void dvec2_negate(dvec2 *v) {
  v->x = -v->x;
  v->y = -v->y;
}

void dvec2_normalize(dvec2 *v) {
  dbl vnorm = sqrt(v->x*v->x + v->y*v->y);
  v->x /= vnorm;
  v->y /= vnorm;
}

dvec2 dvec2_cproj(dvec2 u, dvec2 v) {
  dvec2 w = {
    .x = (1 - u.x*u.x)*v.x - u.x*u.y*v.y,
    .y = -u.x*u.y*v.x + (1 - u.y*u.y)*v.y
  };
  return w;
}

dvec2 dvec2_avg(dvec2 u, dvec2 v) {
  return (dvec2) {(u.x + v.x)/2, (u.y + v.y)/2};
}

dbl dvec2_sum(dvec2 u) {
  return u.x + u.y;
}

dvec3 dvec3_add(dvec3 u, dvec3 v) {
  return (dvec3) {.packed = _mm256_add_pd(u.packed, v.packed)};
}

dvec3 dvec3_dbl_div(dvec3 u, dbl a) {
  dvec3 v;
  v.packed = _mm256_broadcast_sd(&a);
  v.packed = _mm256_div_pd(u.packed, v.packed);
  return v;
}

dvec3 dvec3_dbl_mul(dvec3 u, dbl a) {
  dvec3 v;
  v.packed = _mm256_broadcast_sd(&a);
  v.packed = _mm256_mul_pd(u.packed, v.packed);
  return v;
}

dbl dvec3_dist(dvec3 u, dvec3 v) {
  return dvec3_norm(dvec3_sub(u, v));
}

dbl dvec3_dist_sq(dvec3 u, dvec3 v) {
  return dvec3_norm_sq(dvec3_sub(u, v));
}

dbl dvec3_dot(dvec3 u, dvec3 v) {
  u.packed = _mm256_mul_pd(u.packed, v.packed);
  return u.data[0] + u.data[1] + u.data[2];
}

dvec3 dvec3_infinity(void) {
  return (dvec3) {
    .data = {
      INFINITY,
      INFINITY,
      INFINITY
    }
  };
}

dbl dvec3_maxdist(dvec3 u, dvec3 v) {
  return dvec3_maxnorm(dvec3_sub(u, v));
}

dbl dvec3_maxnorm(dvec3 u) {
  return fmax(fabs(u.data[0]), fmax(fabs(u.data[1]), fabs(u.data[2])));
}

dvec3 dvec3_nan(void) {
  return (dvec3) {
    .data = {NAN, NAN, NAN}
  };
}

dbl dvec3_norm(dvec3 u) {
  u.packed = _mm256_mul_pd(u.packed, u.packed);
  return sqrt(u.data[0] + u.data[1] + u.data[2]);
}

dbl dvec3_norm_sq(dvec3 u) {
  u.packed = _mm256_mul_pd(u.packed, u.packed);
  return u.data[0] + u.data[1] + u.data[2];
}

dvec3 dvec3_normalized(dvec3 u) {
  return dvec3_dbl_div(u, dvec3_norm(u));
}

dvec3 dvec3_one(void) {
  return (dvec3) {.data = {1, 1, 1}};
}

dvec3 dvec3_saxpy(dbl a, dvec3 u, dvec3 v) {
  dvec3 w;
  w.packed = _mm256_broadcast_sd(&a);
  w.packed = _mm256_mul_pd(w.packed, u.packed);
  w.packed = _mm256_add_pd(w.packed, v.packed);
  return w;
}

dvec3 dvec3_sub(dvec3 u, dvec3 v) {
  return (dvec3) {.packed = _mm256_sub_pd(u.packed, v.packed)};
}

dvec3 dvec3_zero(void) {
  return (dvec3) {.data = {0, 0, 0}};
}

int dvec3_argmax(dvec3 u) {
  dbl umax = -INFINITY;
  int argmax;
  for (int i = 0; i < 3; ++i) {
    if (u.data[i] > umax) {
      umax = u.data[i];
      argmax = i;
    }
  }
  return argmax;
}

void dvec3_normalize(dvec3 *u) {
  dbl unorm = u->data[0]*u->data[0] + u->data[1]*u->data[1]
    + u->data[2]*u->data[2];
  unorm /= sqrt(unorm);
  u->data[0] /= unorm;
  u->data[1] /= unorm;
  u->data[2] /= unorm;
}

dbl dvec4_dot(dvec4 v0, dvec4 v1) {
  dbl tmp = 0;
  for (int i = 0; i < 4; ++i) {
    tmp += v0.data[i]*v1.data[i];
  }
  return tmp;
}

dbl dvec4_sum(dvec4 v) {
  dbl tmp = 0;
  for (int i = 0; i < 4; ++i) {
    tmp += v.data[i];
  }
  return tmp;
}

dvec4 dvec4_add(dvec4 u, dvec4 v) {
  dvec4 w;
  w.data[0] = u.data[0] + v.data[0];
  w.data[1] = u.data[1] + v.data[1];
  w.data[2] = u.data[2] + v.data[2];
  w.data[3] = u.data[3] + v.data[3];
  return w;
}

dvec4 dvec4_dbl_div(dvec4 u, dbl a) {
  dvec4 v;
  v.data[0] = u.data[0]/a;
  v.data[1] = u.data[1]/a;
  v.data[2] = u.data[2]/a;
  v.data[3] = u.data[3]/a;
  return v;
}

dvec4 dvec4_m(dbl x) {
  dvec4 m;
  m.data[0] = 1.0;
  m.data[1] = x;
  m.data[2] = x*x;
  m.data[3] = x*x*x;
  return m;
}

dvec4 dvec4_dm(dbl x) {
  dvec4 dm;
  dm.data[0] = 0.0;
  dm.data[1] = 1.0;
  dm.data[2] = 2.0*x;
  dm.data[3] = 3.0*x*x;
  return dm;
}

dvec4 dvec4_d2m(dbl x) {
  return (dvec4) {
    .data = {
      0.0,
      0.0,
      2.0,
      6.0*x
    }
  };
}

dvec4 dvec4_e1(void) {
  dvec4 e1;
  e1.data[0] = 1.0;
  e1.data[1] = 0.0;
  e1.data[2] = 0.0;
  e1.data[3] = 0.0;
  return e1;
}

dvec4 dvec4_one(void) {
  dvec4 one;
  one.data[0] = 1.0;
  one.data[1] = 1.0;
  one.data[2] = 1.0;
  one.data[3] = 1.0;
  return one;
}

dvec4 dvec4_iota(void) {
  dvec4 iota;
  iota.data[0] = 0.0;
  iota.data[1] = 1.0;
  iota.data[2] = 2.0;
  iota.data[3] = 3.0;
  return iota;
}

ivec2 ivec2_add(ivec2 p, ivec2 q) {
  return (ivec2) {p.i + q.i, p.j + q.j};
}

ivec3 ivec3_from_int3(int const p[3]) {
  return (ivec3) {.data = {p[0], p[1], p[2]}};
}

ivec3 ivec3_add(ivec3 p, ivec3 q) {
  return (ivec3) {
    .data = {
      p.data[0] + q.data[0],
      p.data[1] + q.data[1],
      p.data[2] + q.data[2]
    }
  };
}

int ivec3_prod(ivec3 p) {
  return p.data[0]*p.data[1]*p.data[2];
}

dvec3 ivec3_dbl_mul(ivec3 p, dbl a) {
  return (dvec3) {
    .packed = _mm256_mul_pd(ivec3_to_dvec3(p).packed, _mm256_broadcast_sd(&a))
  };
}

bool ivec3_equal(ivec3 p, ivec3 q) {
  return p.data[0] == q.data[0] && p.data[1] == q.data[1] && p.data[2] == q.data[2];
}

ivec3 ivec3_int_div(ivec3 p, int q) {
  return (ivec3) {.data = {p.data[0]/q, p.data[1]/q, p.data[2]/q}};
}

ivec2 dvec2_to_ivec2(dvec2 v) {
  ivec2 ij = {.i = (int)v.x, .j = (int)v.y};
  return ij;
}

ivec3 dvec3_to_ivec3(dvec3 x) {
  return (ivec3) {
    .data = {(int)x.data[0], (int)x.data[1], (int)x.data[2]}
  };
}

dvec3 ivec3_to_dvec3(ivec3 p) {
  return (dvec3) {
    .data = {
      (dbl) p.data[0],
      (dbl) p.data[1],
      (dbl) p.data[2],
      (dbl) p.data[3]
    }
  };
}
