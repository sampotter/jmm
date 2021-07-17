#include "bicubic.h"

#include <string.h>

static dbl44 V_inv = {
  { 1,  0,  0,  0},
  { 0,  0,  1,  0},
  {-3,  3, -2, -1},
  { 2, -2,  1,  1}
};

static dbl44 V_inv_tr = {
  {1, 0, -3,  2},
  {0, 0,  3, -2},
  {0, 1, -2,  1},
  {0, 0, -1,  1}
};

static dbl44 D = {
  {0, 0, 0, 0},
  {1, 0, 0, 0},
  {0, 2, 0, 0},
  {0, 0, 3, 0}
};

static dbl44 D_tr = {
  {0, 1, 0, 0},
  {0, 0, 2, 0},
  {0, 0, 0, 3},
  {0, 0, 0, 0}
};

void bicubic_set_data(bicubic_s *bicubic, dbl44 data) {
  dbl44_mul(V_inv, data, bicubic->A);
  dbl44_mul(bicubic->A, V_inv_tr, bicubic->A);
}

void bicubic_set_data_from_ptr(bicubic_s *bicubic, dbl const *data_ptr) {
  memcpy((void *)bicubic->A, (void *)data_ptr, sizeof(dbl44));
}

static void restrict_A(dbl44 const A, bicubic_variable var, int edge, dbl4 a) {
  dbl4 x;
  if (edge == 0)
    dbl4_e1(x);
  else
    dbl4_one(x);
  if (var == LAMBDA)
    dbl44_dbl4_mul(A, x, a);
  else
    dbl4_dbl44_mul(x, A, a);
}

cubic_s
bicubic_get_f_on_edge(bicubic_s const *bicubic, bicubic_variable var, int edge) {
  cubic_s cubic;
  restrict_A(bicubic->A, var, edge, cubic.a);
  return cubic;
}

cubic_s
bicubic_get_fx_on_edge(bicubic_s const *bicubic, bicubic_variable var, int edge) {
  dbl44 Ax; dbl44_mul(D_tr, bicubic->A, Ax);
  cubic_s cubic;
  restrict_A(Ax, var, edge, cubic.a);
  return cubic;
}

cubic_s
bicubic_get_fy_on_edge(bicubic_s const *bicubic, bicubic_variable var, int edge) {
  dbl44 Ay; dbl44_mul(bicubic->A, D, Ay);
  cubic_s cubic;
  restrict_A(Ay, var, edge, cubic.a);
  return cubic;
}

dbl bicubic_f(bicubic_s const *bicubic, dbl2 cc) {
  dbl4 mx; dbl4_m(cc[0], mx);
  dbl4 my; dbl4_m(cc[1], my);
  return dbl4_dbl44_dbl4_dot(mx, bicubic->A, my);
}

dbl bicubic_fx(bicubic_s const *bicubic, dbl2 cc) {
  dbl4 dmx; dbl4_dm(cc[0], dmx);
  dbl4 my; dbl4_m(cc[1], my);
  return dbl4_dbl44_dbl4_dot(dmx, bicubic->A, my);
}

dbl bicubic_fy(bicubic_s const *bicubic, dbl2 cc) {
  dbl4 mx; dbl4_m(cc[0], mx);
  dbl4 dmy; dbl4_dm(cc[1], dmy);
  return dbl4_dbl44_dbl4_dot(mx, bicubic->A, dmy);
}

dbl bicubic_fxy(bicubic_s const *bicubic, dbl2 cc) {
  dbl4 dmx; dbl4_dm(cc[0], dmx);
  dbl4 dmy; dbl4_dm(cc[1], dmy);
  return dbl4_dbl44_dbl4_dot(dmx, bicubic->A, dmy);
}

void interpolate_fxy_at_verts(dbl4 const fx, dbl4 const fy, dbl h, dbl4 fxy) {
  static dbl44 Ax = {
    {-3./4,  3./4,  1./4, -1./4},
    { 1./4, -1./4, -3./4,  3./4},
    {-3./4,  3./4,  1./4, -1./4},
    { 1./4, -1./4, -3./4,  3./4}
  };
  static dbl44 Ay = {
    {-3./4,  3./4,  1./4, -1./4},
    {-3./4,  3./4,  1./4, -1./4},
    { 1./4, -1./4, -3./4,  3./4},
    { 1./4, -1./4, -3./4,  3./4}
  };
  dbl4 tmp1; dbl44_dbl4_mul(Ax, fx, tmp1);
  dbl4 tmp2; dbl44_dbl4_mul(Ay, fy, tmp2);
  dbl4_add(tmp1, tmp2, fxy);
  dbl4_dbl_div_inplace(fxy, h);
}

bool bicubic_valid(bicubic_s const *bicubic) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      if (!isfinite(bicubic->A[i][j]))
        return false;
  return true;
}

void bicubic_invalidate(bicubic_s *bicubic) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      bicubic->A[i][j] = NAN;
}
