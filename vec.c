#include "vec.h"

dvec2 dvec2_ccomb(dvec2 v0, dvec2 v1, dbl t) {
  dvec2 vt = {(1 - t)*v0.x + t*v1.x, (1 - t)*v0.y + t*v1.y};
  return vt;
}

dbl dvec2_dist(dvec2 v0, dvec2 v1) {
  dbl dx = v1.x - v0.x, dy = v1.y - v0.y;
  return sqrt(dx*dx + dy*dy);
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
