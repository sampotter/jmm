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
  cubic_s cubic;
  cubic.a = var == LAMBDA ?
      dmat44_dvec4_mul(bicubic->A, edge == 0 ? dvec4_e1() : dvec4_one()) :
      dvec4_dmat44_mul(edge == 0 ? dvec4_e1() : dvec4_one(), bicubic->A);
  return cubic;
}

dbl bicubic_f(bicubic_s const *bicubic, dvec2 cc) {
  return dvec4_dot(
    dvec4_m(cc.x),
    dmat44_dvec4_mul(bicubic->A, dvec4_m(cc.y))
  );
}

dbl bicubic_fx(bicubic_s const *bicubic, dvec2 cc) {
  return dvec4_dot(
    dvec4_dm(cc.x),
    dmat44_dvec4_mul(bicubic->A, dvec4_m(cc.y))
  );
}

dbl bicubic_fy(bicubic_s const *bicubic, dvec2 cc) {
  return dvec4_dot(
    dvec4_m(cc.x),
    dmat44_dvec4_mul(bicubic->A, dvec4_dm(cc.y))
  );
}

dbl bicubic_fxy(bicubic_s const *bicubic, dvec2 cc) {
  return dvec4_dot(
    dvec4_dm(cc.x),
    dmat44_dvec4_mul(bicubic->A, dvec4_dm(cc.y))
  );
}

/**
 * TODO: move this into the `bicubic` module, add some unit tests,
 * etc.
 */
dvec4 interpolate_fxy_at_verts(dvec4 fx, dvec4 fy, dbl h) {
  static dmat44 Ax = {
    .data = {
      {-3.0/4.0,  3.0/4.0,  1.0/4.0, -1.0/4.0},
      { 1.0/4.0, -1.0/4.0, -3.0/4.0,  3.0/4.0},
      {-3.0/4.0,  3.0/4.0,  1.0/4.0, -1.0/4.0},
      { 1.0/4.0, -1.0/4.0, -3.0/4.0,  3.0/4.0}
    }
  };
  static dmat44 Ay = {
    .data = {
      {-3.0/4.0,  3.0/4.0,  1.0/4.0, -1.0/4.0},
      {-3.0/4.0,  3.0/4.0,  1.0/4.0, -1.0/4.0},
      { 1.0/4.0, -1.0/4.0, -3.0/4.0,  3.0/4.0},
      { 1.0/4.0, -1.0/4.0, -3.0/4.0,  3.0/4.0}
    }
  };
  return dvec4_dbl_div(
    dvec4_add(
      dmat44_dvec4_mul(Ax, fx),
      dmat44_dvec4_mul(Ay, fy)
    ),
    h
  );
}

#if SJS_DEBUG
bool bicubic_valid(bicubic_s const *bicubic) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (isnan(bicubic->A.data[i][j])) {
        return true;
      }
    }
  }
  return true;
}

void bicubic_invalidate(bicubic_s *bicubic) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      bicubic->A.data[i][j] = NAN;
    }
  }
}
#endif
