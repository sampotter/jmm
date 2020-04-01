#include "vec.h"

dvec2 dvec2_zero() {
  return (dvec2) {0.0, 0.0};
}

dvec2 dvec2_ccomb(dvec2 v0, dvec2 v1, dbl t) {
  dvec2 vt = {(1 - t)*v0.x + t*v1.x, (1 - t)*v0.y + t*v1.y};
  return vt;
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

dvec4 dvec4_e1() {
  dvec4 e1;
  e1.data[0] = 1.0;
  e1.data[1] = 0.0;
  e1.data[2] = 0.0;
  e1.data[3] = 0.0;
  return e1;
}

dvec4 dvec4_one() {
  dvec4 one;
  one.data[0] = 1.0;
  one.data[1] = 1.0;
  one.data[2] = 1.0;
  one.data[3] = 1.0;
  return one;
}

dvec4 dvec4_iota() {
  dvec4 iota;
  iota.data[0] = 0.0;
  iota.data[1] = 1.0;
  iota.data[2] = 2.0;
  iota.data[3] = 3.0;
  return iota;
}

ivec2 dvec2_to_ivec2(dvec2 v) {
  ivec2 ij = {.i = (int)v.x, .j = (int)v.y};
  return ij;
}
