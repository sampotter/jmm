#include "bicubic.h"

#include <string.h>

static dmat44 V_inv = {
  .data = {
   { 1,  0,  0,  0},
   { 0,  0,  1,  0},
   {-3,  3, -2, -1},
   { 2, -2,  1,  1}
  }
};

static dmat44 V_inv_tr = {
  .data = {
    {1, 0, -3,  2},
    {0, 0,  3, -2},
    {0, 1, -2,  1},
    {0, 0, -1,  1}
  }
};

void bicubic_set_data(bicubic_s *bicubic, dmat44 data) {
  bicubic->A = dmat44_dmat44_mul(V_inv, data);
  bicubic->A = dmat44_dmat44_mul(bicubic->A, V_inv_tr);
}

void bicubic_set_data_from_ptr(bicubic_s *bicubic, dbl const *data_ptr) {
  memcpy((void *)bicubic->A.data, (void *)data_ptr, 16*sizeof(dbl));
}

cubic_s
bicubic_restrict(bicubic_s const *bicubic, bicubic_variable var, int edge) {
  // TODO: a bunch of this is wrong---fix
  cubic_s cubic;
  if (var == LAMBDA) {
    if (edge == 0) {
      cubic.a = bicubic->A.rows[0];
    } else {
      for (int i = 0; i < 4; ++i) {
        cubic.a.data[i] = dvec4_sum(bicubic->A.rows[i]);
      }
    }
  } else {
    if (edge == 0) {
      cubic.a = dmat44_col(bicubic->A, 0);
    } else {
      for (int i = 0; i < 4; ++i) {
        cubic.a.data[i] = dvec4_sum(bicubic->A.rows[i]);
      }
    }
  }
  return cubic;
}

dbl bicubic_f(bicubic_s const *bicubic, dvec2 cc) {
  (void) bicubic;
  (void) cc;
  dbl f = 0;
  return f;
}

dbl bicubic_fx(bicubic_s const *bicubic, dvec2 cc) {
  (void) bicubic;
  (void) cc;
  dbl fx = 0;
  return fx;
}

dbl bicubic_fy(bicubic_s const *bicubic, dvec2 cc) {
  (void) bicubic;
  (void) cc;
  dbl fy = 0;
  return fy;
}

dbl bicubic_fxy(bicubic_s const *bicubic, dvec2 cc) {
  (void) bicubic;
  (void) cc;
  dbl fxy = 0;
  return fxy;
}
