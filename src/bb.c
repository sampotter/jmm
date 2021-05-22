#include "bb.h"

#include <assert.h>

#include "mesh3.h"
#include "vec.h"

#define TRI000 0

/**
 * [0] [1]
 * [2]
 *
 * 100 010
 * 001
 */
#define TRI100 0
#define TRI010 1
#define TRI001 2

/**
 * [0] [1] [2]
 * [3] [4]
 * [5]
 *
 * 200 110 020
 * 101 011
 * 002
 */
#define TRI200 0
#define TRI110 1
#define TRI020 2
#define TRI101 3
#define TRI011 4
#define TRI002 5

/**
 * [0] [1] [2] [3]
 * [4] [5] [6]
 * [7] [8]
 * [9]
 *
 * 300 210 120 030
 * 201 111 021
 * 102 012
 * 003
 */
#define TRI300 0
#define TRI210 1
#define TRI120 2
#define TRI030 3
#define TRI201 4
#define TRI111 5
#define TRI021 6
#define TRI102 7
#define TRI012 8
#define TRI003 9

/**
 * [0]  [1]
 * [2]
 *
 * [3]
 *
 * 1000 0100
 * 0010
 *
 * 0001
 */
#define TET1000 0
#define TET0100 1
#define TET0010 2
#define TET0001 3

/**
 * [0]  [1]  [2]
 * [3]  [4]
 * [5]
 *
 * [6]  [7]
 * [8]
 *
 * [9]
 *
 * 2000 1100 0200
 * 1010 0110
 * 0020
 *
 * 1001 0101
 * 0011
 *
 * 0002
 */
#define TET2000 0
#define TET1100 1
#define TET0200 2
#define TET1010 3
#define TET0110 4
#define TET0020 5
#define TET1001 6
#define TET0101 7
#define TET0011 8
#define TET0002 9

/**
 * [0]  [1]  [2]  [3]
 * [4]  [5]  [6]
 * [7]  [8]
 * [9]
 *
 * [10] [11] [12]
 * [13] [14]
 * [15]
 *
 * [16] [17]
 * [18]
 *
 * [19]
 *
 * 3000 2100 1200 0300
 * 2010 1110 0210
 * 1020 0120
 * 0030
 *
 * 2001 1101 0201
 * 1011 0111
 * 0021
 *
 * 1002 0102
 * 0012
 *
 * 0003
 */
#define TET3000 0
#define TET2100 1
#define TET1200 2
#define TET0300 3
#define TET2010 4
#define TET1110 5
#define TET0210 6
#define TET1020 7
#define TET0120 8
#define TET0030 9
#define TET2001 10
#define TET1101 11
#define TET0201 12
#define TET1011 13
#define TET0111 14
#define TET0021 15
#define TET1002 16
#define TET0102 17
#define TET0012 18
#define TET0003 19

void bb31_init_from_1d_data(bb31 *bb, dbl const f[2], dbl const Df[2], dbl const x[2]) {
  bb->c[1] = (bb->c[0] = f[0]) + Df[0]*(x[1] - x[0])/3;
  bb->c[2] = (bb->c[3] = f[1]) + Df[1]*(x[0] - x[1])/3;
}

void bb31_init_from_3d_data(bb31 *bb, dbl const f[2], dbl const Df[2][3], dbl const x[2][3]) {
  dbl dx[3];
  dbl3_sub(x[1], x[0], dx);
  bb->c[1] = (bb->c[0] = f[0]) + dbl3_dot(dx, Df[0])/3;
  bb->c[2] = (bb->c[3] = f[1]) - dbl3_dot(dx, Df[1])/3;
}

dbl bb31_f(bb31 const *bb, dbl const *b) {
  dbl tmp[3];

  tmp[0] = b[0]*bb->c[0] + b[1]*bb->c[1];
  tmp[1] = b[0]*bb->c[1] + b[1]*bb->c[2];
  tmp[2] = b[0]*bb->c[2] + b[1]*bb->c[3];

  tmp[0] = b[0]*tmp[0] + b[1]*tmp[1];
  tmp[1] = b[0]*tmp[1] + b[1]*tmp[2];

  return b[0]*tmp[0] + b[1]*tmp[1];
}

dbl bb31_df(bb31 const *bb, dbl const *b, dbl const *a) {
  dbl tmp[3];

  tmp[0] = a[0]*bb->c[0] + a[1]*bb->c[1];
  tmp[1] = a[0]*bb->c[1] + a[1]*bb->c[2];
  tmp[2] = a[0]*bb->c[2] + a[1]*bb->c[3];

  tmp[0] = b[0]*tmp[0] + b[1]*tmp[1];
  tmp[1] = b[0]*tmp[1] + b[1]*tmp[2];

  return 3*(b[0]*tmp[0] + b[1]*tmp[1]);
}

dbl bb31_d2f(bb31 const *bb, dbl const *b, dbl const *a) {
  dbl tmp[3];

  tmp[0] = a[0]*bb->c[0] + a[1]*bb->c[1];
  tmp[1] = a[0]*bb->c[1] + a[1]*bb->c[2];
  tmp[2] = a[0]*bb->c[2] + a[1]*bb->c[3];

  tmp[0] = a[0]*tmp[0] + a[1]*tmp[1];
  tmp[1] = a[0]*tmp[1] + a[1]*tmp[2];

  return 3*(b[0]*tmp[0] + b[1]*tmp[1]);
}

void bb31_reverse(bb31 *bb) {
  dbl tmp[4];
  for (size_t i = 0; i < 4; ++i) tmp[i] = bb->c[3 - i];
  for (size_t i = 0; i < 4; ++i) bb->c[i] = tmp[i];
}

void bb32_init_from_3d_data(bb32 *bb, dbl const f[3], dbl const Df[3][3], dbl const x[3][3]) {
  dbl *c = bb->c;

  c[TRI300] = f[0];
  c[TRI030] = f[1];
  c[TRI003] = f[2];

  dbl tmp[6];

  dbl3_sub(x[TRI010], x[TRI100], tmp);
  c[TRI210] = c[TRI300] + dbl3_ndot(Df[TRI100], tmp)/3;
  c[TRI120] = c[TRI030] - dbl3_ndot(Df[TRI010], tmp)/3;

  dbl3_sub(x[TRI001], x[TRI100], tmp);
  c[TRI201] = c[TRI300] + dbl3_ndot(Df[TRI100], tmp)/3;
  c[TRI102] = c[TRI003] - dbl3_ndot(Df[TRI001], tmp)/3;

  dbl3_sub(x[TRI001], x[TRI010], tmp);
  c[TRI021] = c[TRI030] + dbl3_ndot(Df[TRI010], tmp)/3;
  c[TRI012] = c[TRI003] - dbl3_ndot(Df[TRI001], tmp)/3;

  tmp[0] = c[TRI210];
  tmp[1] = c[TRI201];
  tmp[2] = c[TRI120];
  tmp[3] = c[TRI021];
  tmp[4] = c[TRI012];
  tmp[5] = c[TRI102];
  c[TRI111] = dblN_nsum(tmp, 6)/4;

  tmp[0] = c[TRI300];
  tmp[1] = c[TRI030];
  tmp[2] = c[TRI003];
  c[TRI111] -= dbl3_nsum(tmp)/6;
}

void bb32_init_from_jets(bb32 *bb, jet3 const jet[3], dbl const x[3][3],
                         bool const is_point_source[3]) {
  /* Assume at least one of the jets has gradient information */
  assert(!(is_point_source[0] && is_point_source[1] && is_point_source[2]));

  dbl *c = bb->c;

  c[TRI300] = jet[0].f;
  c[TRI030] = jet[1].f;
  c[TRI003] = jet[2].f;

  dbl tmp[6];

  /* First, compute the data for the ordinates neighboring each of the
   * jets with gradient data */

  if (!is_point_source[0]) {
    dbl3_sub(x[TRI010], x[TRI100], tmp);
    c[TRI210] = c[TRI300] + dbl3_ndot(&jet[TRI100].fx, tmp)/3;

    dbl3_sub(x[TRI001], x[TRI100], tmp);
    c[TRI201] = c[TRI300] + dbl3_ndot(&jet[TRI100].fx, tmp)/3;
  }

  if (!is_point_source[1]) {
    dbl3_sub(x[TRI001], x[TRI010], tmp);
    c[TRI021] = c[TRI030] + dbl3_ndot(&jet[TRI010].fx, tmp)/3;

    dbl3_sub(x[TRI100], x[TRI010], tmp);
    c[TRI120] = c[TRI030] + dbl3_ndot(&jet[TRI010].fx, tmp)/3;
  }

  if (!is_point_source[2]) {
    dbl3_sub(x[TRI100], x[TRI001], tmp);
    c[TRI102] = c[TRI003] - dbl3_ndot(&jet[TRI001].fx, tmp)/3;

    dbl3_sub(x[TRI010], x[TRI001], tmp);
    c[TRI012] = c[TRI003] + dbl3_ndot(&jet[TRI001].fx, tmp)/3;
  }

  /* Next, use condensation of parameters to compute the ordinates
   * neighboring the jets *without* gradient data */

  if (is_point_source[0]) {
    c[TRI210] = is_point_source[1] ?
      (2*c[TRI300] + c[TRI030])/3 : (c[TRI300] + c[TRI120])/2;
    c[TRI201] = is_point_source[2] ?
      (2*c[TRI300] + c[TRI003])/3 : (c[TRI300] + c[TRI102])/2;
  }

  if (is_point_source[1]) {
    c[TRI021] = is_point_source[2] ?
      (2*c[TRI030] + c[TRI003])/3 : (c[TRI030] + c[TRI012])/2;
    c[TRI120] = is_point_source[0] ?
      (2*c[TRI030] + c[TRI300])/3 : (c[TRI030] + c[TRI210])/2;
  }

  if (is_point_source[2]) {
    c[TRI102] = is_point_source[0] ?
      (2*c[TRI003] + c[TRI300])/3 : (c[TRI003] + c[TRI201])/2;
    c[TRI012] = is_point_source[1] ?
      (2*c[TRI003] + c[TRI030])/3 : (c[TRI003] + c[TRI021])/2;
  }

  /* Use condensation of parameters to compute the 10th parameter from
   * the rest */

  tmp[0] = c[TRI210];
  tmp[1] = c[TRI201];
  tmp[2] = c[TRI120];
  tmp[3] = c[TRI021];
  tmp[4] = c[TRI012];
  tmp[5] = c[TRI102];
  c[TRI111] = dblN_nsum(tmp, 6)/4;

  tmp[0] = c[TRI300];
  tmp[1] = c[TRI030];
  tmp[2] = c[TRI003];
  c[TRI111] -= dbl3_nsum(tmp)/6;
}

static void bb32_reduce1(dbl const *in1, dbl const *in2, dbl *out) {
  dbl tmp[3];

  tmp[0] = in1[TRI100]*in2[TRI300];
  tmp[1] = in1[TRI010]*in2[TRI210];
  tmp[2] = in1[TRI001]*in2[TRI201];
  out[TRI200] = dbl3_nsum(tmp);

  tmp[0] = in1[TRI100]*in2[TRI210];
  tmp[1] = in1[TRI010]*in2[TRI120];
  tmp[2] = in1[TRI001]*in2[TRI111];
  out[TRI110] = dbl3_nsum(tmp);

  tmp[0] = in1[TRI100]*in2[TRI120];
  tmp[1] = in1[TRI010]*in2[TRI030];
  tmp[2] = in1[TRI001]*in2[TRI021];
  out[TRI020] = dbl3_nsum(tmp);

  tmp[0] = in1[TRI100]*in2[TRI201];
  tmp[1] = in1[TRI010]*in2[TRI111];
  tmp[2] = in1[TRI001]*in2[TRI102];
  out[TRI101] = dbl3_nsum(tmp);

  tmp[0] = in1[TRI100]*in2[TRI111];
  tmp[1] = in1[TRI010]*in2[TRI021];
  tmp[2] = in1[TRI001]*in2[TRI012];
  out[TRI011] = dbl3_nsum(tmp);

  tmp[0] = in1[TRI100]*in2[TRI102];
  tmp[1] = in1[TRI010]*in2[TRI012];
  tmp[2] = in1[TRI001]*in2[TRI003];
  out[TRI002] = dbl3_nsum(tmp);
}

static void bb32_reduce2(dbl const *in1, dbl const *in2, dbl *out) {
  dbl tmp[3];

  tmp[0] = in1[TRI100]*in2[TRI200];
  tmp[1] = in1[TRI010]*in2[TRI110];
  tmp[2] = in1[TRI001]*in2[TRI101];
  out[TRI100] = dbl3_nsum(tmp);

  tmp[0] = in1[TRI100]*in2[TRI110];
  tmp[1] = in1[TRI010]*in2[TRI020];
  tmp[2] = in1[TRI001]*in2[TRI011];
  out[TRI010] = dbl3_nsum(tmp);

  tmp[0] = in1[TRI100]*in2[TRI101];
  tmp[1] = in1[TRI010]*in2[TRI011];
  tmp[2] = in1[TRI001]*in2[TRI002];
  out[TRI001] = dbl3_nsum(tmp);
}

static void bb32_reduce3(dbl const *in1, dbl const *in2, dbl *out) {
  dbl tmp[3];

  tmp[0] = in1[TRI100]*in2[TRI100];
  tmp[1] = in1[TRI010]*in2[TRI010];
  tmp[2] = in1[TRI001]*in2[TRI001];
  out[TRI000] = dbl3_nsum(tmp);
}

dbl bb32_f(bb32 const *bb, dbl const *b) {
  dbl tmp[6];
  bb32_reduce1(b, bb->c, tmp);
  bb32_reduce2(b, tmp, tmp);
  bb32_reduce3(b, tmp, tmp);
  return tmp[TRI000];
}

dbl bb32_df(bb32 const *bb, dbl const *b, dbl const *a) {
  dbl tmp[6];
  bb32_reduce1(a, bb->c, tmp);
  bb32_reduce2(b, tmp, tmp);
  bb32_reduce3(b, tmp, tmp);
  return 3*tmp[TRI000];
}

dbl bb32_d2f(bb32 const *bb, dbl const *b, dbl const *a1, dbl const *a2) {
  dbl tmp[6];
  bb32_reduce1(a1, bb->c, tmp);
  bb32_reduce2(a2, tmp, tmp);
  bb32_reduce3(b, tmp, tmp);
  return 6*tmp[TRI000];
}

void bb33_init_from_3d_data(bb33 *bb, dbl const f[4], dbl const Df[4][3], dbl const x[4][3]) {
  dbl dx[3];

  dbl *c = bb->c;

  /**
   * Start by computing the Bezier ordinates for each vertex of the
   * tetrahedron. (These are just the function values.)
   */

  c[TET3000] = f[0];
  c[TET0300] = f[1];
  c[TET0030] = f[2];
  c[TET0003] = f[3];

  /**
   * Next, compute Bezier ordinates along tetrahedron edges first.
   */

  // 1 <-> 2
  dbl3_sub(x[TET0100], x[TET1000], dx);
  c[TET2100] = c[TET3000] + dbl3_ndot(Df[TET1000], dx)/3;
  c[TET1200] = c[TET0300] - dbl3_ndot(Df[TET0100], dx)/3;

  // 1 <-> 3
  dbl3_sub(x[TET0010], x[TET1000], dx);
  c[TET2010] = c[TET3000] + dbl3_ndot(Df[TET1000], dx)/3;
  c[TET1020] = c[TET0030] - dbl3_ndot(Df[TET0010], dx)/3;

  // 1 <-> 4
  dbl3_sub(x[TET0001], x[TET1000], dx);
  c[TET2001] = c[TET3000] + dbl3_ndot(Df[TET1000], dx)/3;
  c[TET1002] = c[TET0003] - dbl3_ndot(Df[TET0001], dx)/3;

  // 2 <-> 3
  dbl3_sub(x[TET0010], x[TET0100], dx);
  c[TET0210] = c[TET0300] + dbl3_ndot(Df[TET0100], dx)/3;
  c[TET0120] = c[TET0030] - dbl3_ndot(Df[TET0010], dx)/3;

  // 2 <-> 4
  dbl3_sub(x[TET0001], x[TET0100], dx);
  c[TET0201] = c[TET0300] + dbl3_ndot(Df[TET0100], dx)/3;
  c[TET0102] = c[TET0003] - dbl3_ndot(Df[TET0001], dx)/3;

  // 3 <-> 4
  dbl3_sub(x[TET0001], x[TET0010], dx);
  c[TET0021] = c[TET0030] + dbl3_ndot(Df[TET0010], dx)/3;
  c[TET0012] = c[TET0003] - dbl3_ndot(Df[TET0001], dx)/3;

  /**
   * Finally, compute Bezier ordinates in the center of each
   * tetrahedron face.
   */

  c[TET1110] = (c[TET2100] + c[TET2010] + c[TET1200] + c[TET0210] + c[TET0120] + c[TET1020])/4
    - (c[TET3000] + c[TET0300] + c[TET0030])/6;

  c[TET1101] = (c[TET2100] + c[TET2001] + c[TET1200] + c[TET0201] + c[TET0102] + c[TET1002])/4
    - (c[TET3000] + c[TET0300] + c[TET0003])/6;

  c[TET1011] = (c[TET2010] + c[TET2001] + c[TET1020] + c[TET0021] + c[TET0012] + c[TET1002])/4
    - (c[TET3000] + c[TET0030] + c[TET0003])/6;

  c[TET0111] = (c[TET0210] + c[TET0201] + c[TET0120] + c[TET0021] + c[TET0012] + c[TET0102])/4
    - (c[TET0300] + c[TET0030] + c[TET0003])/6;

  // (We're done---there are no interior Bezier ordinates for the 3D
  // generalization of the 9-parameter interpolant.)
}

void bb33_init_from_cell_and_jets(bb33 *bb, mesh3_s const *mesh, jet3 const *jet, size_t lc) {
  dbl y[4], Dy[4][3], x[4][3];
  size_t lv[4];
  mesh3_cv(mesh, lc, lv);
  for (int i = 0; i < 4; ++i) {
    y[i] = jet[lv[i]].f;
    Dy[i][0] = jet[lv[i]].fx;
    Dy[i][1] = jet[lv[i]].fy;
    Dy[i][2] = jet[lv[i]].fz;
    mesh3_copy_vert(mesh, lv[i], x[i]);
  }
  bb33_init_from_3d_data(bb, y, Dy, x);
}

dbl bb33_f(bb33 const *bb, dbl const b[4]) {
  dbl tmp[10];

  dbl const *c = bb->c;

  tmp[TET2000] = b[TET1000]*c[TET3000] + b[TET0100]*c[TET2100] + b[TET0010]*c[TET2010] + b[TET0001]*c[TET2001];
  tmp[TET1100] = b[TET1000]*c[TET2100] + b[TET0100]*c[TET1200] + b[TET0010]*c[TET1110] + b[TET0001]*c[TET1101];
  tmp[TET0200] = b[TET1000]*c[TET1200] + b[TET0100]*c[TET0300] + b[TET0010]*c[TET0210] + b[TET0001]*c[TET0201];
  tmp[TET1010] = b[TET1000]*c[TET2010] + b[TET0100]*c[TET1110] + b[TET0010]*c[TET1020] + b[TET0001]*c[TET1011];
  tmp[TET0110] = b[TET1000]*c[TET1110] + b[TET0100]*c[TET0210] + b[TET0010]*c[TET0120] + b[TET0001]*c[TET0111];
  tmp[TET0020] = b[TET1000]*c[TET1020] + b[TET0100]*c[TET0120] + b[TET0010]*c[TET0030] + b[TET0001]*c[TET0021];
  tmp[TET1001] = b[TET1000]*c[TET2001] + b[TET0100]*c[TET1101] + b[TET0010]*c[TET1011] + b[TET0001]*c[TET1002];
  tmp[TET0101] = b[TET1000]*c[TET1101] + b[TET0100]*c[TET0201] + b[TET0010]*c[TET0111] + b[TET0001]*c[TET0102];
  tmp[TET0011] = b[TET1000]*c[TET1011] + b[TET0100]*c[TET0111] + b[TET0010]*c[TET0021] + b[TET0001]*c[TET0012];
  tmp[TET0002] = b[TET1000]*c[TET1002] + b[TET0100]*c[TET0102] + b[TET0010]*c[TET0012] + b[TET0001]*c[TET0003];

  tmp[TET1000] = b[TET1000]*tmp[TET2000] + b[TET0100]*tmp[TET1100] + b[TET0010]*tmp[TET1010] + b[TET0001]*tmp[TET1001];
  tmp[TET0100] = b[TET1000]*tmp[TET1100] + b[TET0100]*tmp[TET0200] + b[TET0010]*tmp[TET0110] + b[TET0001]*tmp[TET0101];
  tmp[TET0010] = b[TET1000]*tmp[TET1010] + b[TET0100]*tmp[TET0110] + b[TET0010]*tmp[TET0020] + b[TET0001]*tmp[TET0011];
  tmp[TET0001] = b[TET1000]*tmp[TET1001] + b[TET0100]*tmp[TET0101] + b[TET0010]*tmp[TET0011] + b[TET0001]*tmp[TET0002];

  return b[TET1000]*tmp[TET1000] + b[TET0100]*tmp[TET0100] + b[TET0010]*tmp[TET0010] + b[TET0001]*tmp[TET0001];
}

dbl bb33_df(bb33 const *bb, dbl const b[4], dbl const a[4]) {
  dbl tmp[10];

  dbl const *c = bb->c;

  tmp[TET2000] = a[TET1000]*c[TET3000] + a[TET0100]*c[TET2100] + a[TET0010]*c[TET2010] + a[TET0001]*c[TET2001];
  tmp[TET1100] = a[TET1000]*c[TET2100] + a[TET0100]*c[TET1200] + a[TET0010]*c[TET1110] + a[TET0001]*c[TET1101];
  tmp[TET0200] = a[TET1000]*c[TET1200] + a[TET0100]*c[TET0300] + a[TET0010]*c[TET0210] + a[TET0001]*c[TET0201];
  tmp[TET1010] = a[TET1000]*c[TET2010] + a[TET0100]*c[TET1110] + a[TET0010]*c[TET1020] + a[TET0001]*c[TET1011];
  tmp[TET0110] = a[TET1000]*c[TET1110] + a[TET0100]*c[TET0210] + a[TET0010]*c[TET0120] + a[TET0001]*c[TET0111];
  tmp[TET0020] = a[TET1000]*c[TET1020] + a[TET0100]*c[TET0120] + a[TET0010]*c[TET0030] + a[TET0001]*c[TET0021];
  tmp[TET1001] = a[TET1000]*c[TET2001] + a[TET0100]*c[TET1101] + a[TET0010]*c[TET1011] + a[TET0001]*c[TET1002];
  tmp[TET0101] = a[TET1000]*c[TET1101] + a[TET0100]*c[TET0201] + a[TET0010]*c[TET0111] + a[TET0001]*c[TET0102];
  tmp[TET0011] = a[TET1000]*c[TET1011] + a[TET0100]*c[TET0111] + a[TET0010]*c[TET0021] + a[TET0001]*c[TET0012];
  tmp[TET0002] = a[TET1000]*c[TET1002] + a[TET0100]*c[TET0102] + a[TET0010]*c[TET0012] + a[TET0001]*c[TET0003];

  tmp[TET1000] = b[TET1000]*tmp[TET2000] + b[TET0100]*tmp[TET1100] + b[TET0010]*tmp[TET1010] + b[TET0001]*tmp[TET1001];
  tmp[TET0100] = b[TET1000]*tmp[TET1100] + b[TET0100]*tmp[TET0200] + b[TET0010]*tmp[TET0110] + b[TET0001]*tmp[TET0101];
  tmp[TET0010] = b[TET1000]*tmp[TET1010] + b[TET0100]*tmp[TET0110] + b[TET0010]*tmp[TET0020] + b[TET0001]*tmp[TET0011];
  tmp[TET0001] = b[TET1000]*tmp[TET1001] + b[TET0100]*tmp[TET0101] + b[TET0010]*tmp[TET0011] + b[TET0001]*tmp[TET0002];

  return 3*(b[TET1000]*tmp[TET1000] + b[TET0100]*tmp[TET0100] + b[TET0010]*tmp[TET0010] + b[TET0001]*tmp[TET0001]);
}

dbl bb33_d2f(bb33 const *bb, dbl const b[4], dbl const a[2][4]) {
  dbl tmp[10];

  dbl const *c = bb->c;

  tmp[TET2000] = a[0][TET1000]*c[TET3000] + a[0][TET0100]*c[TET2100] + a[0][TET0010]*c[TET2010] + a[0][TET0001]*c[TET2001];
  tmp[TET1100] = a[0][TET1000]*c[TET2100] + a[0][TET0100]*c[TET1200] + a[0][TET0010]*c[TET1110] + a[0][TET0001]*c[TET1101];
  tmp[TET0200] = a[0][TET1000]*c[TET1200] + a[0][TET0100]*c[TET0300] + a[0][TET0010]*c[TET0210] + a[0][TET0001]*c[TET0201];
  tmp[TET1010] = a[0][TET1000]*c[TET2010] + a[0][TET0100]*c[TET1110] + a[0][TET0010]*c[TET1020] + a[0][TET0001]*c[TET1011];
  tmp[TET0110] = a[0][TET1000]*c[TET1110] + a[0][TET0100]*c[TET0210] + a[0][TET0010]*c[TET0120] + a[0][TET0001]*c[TET0111];
  tmp[TET0020] = a[0][TET1000]*c[TET1020] + a[0][TET0100]*c[TET0120] + a[0][TET0010]*c[TET0030] + a[0][TET0001]*c[TET0021];
  tmp[TET1001] = a[0][TET1000]*c[TET2001] + a[0][TET0100]*c[TET1101] + a[0][TET0010]*c[TET1011] + a[0][TET0001]*c[TET1002];
  tmp[TET0101] = a[0][TET1000]*c[TET1101] + a[0][TET0100]*c[TET0201] + a[0][TET0010]*c[TET0111] + a[0][TET0001]*c[TET0102];
  tmp[TET0011] = a[0][TET1000]*c[TET1011] + a[0][TET0100]*c[TET0111] + a[0][TET0010]*c[TET0021] + a[0][TET0001]*c[TET0012];
  tmp[TET0002] = a[0][TET1000]*c[TET1002] + a[0][TET0100]*c[TET0102] + a[0][TET0010]*c[TET0012] + a[0][TET0001]*c[TET0003];

  tmp[TET1000] = a[1][TET1000]*tmp[TET2000] + a[1][TET0100]*tmp[TET1100] + a[1][TET0010]*tmp[TET1010] + a[1][TET0001]*tmp[TET1001];
  tmp[TET0100] = a[1][TET1000]*tmp[TET1100] + a[1][TET0100]*tmp[TET0200] + a[1][TET0010]*tmp[TET0110] + a[1][TET0001]*tmp[TET0101];
  tmp[TET0010] = a[1][TET1000]*tmp[TET1010] + a[1][TET0100]*tmp[TET0110] + a[1][TET0010]*tmp[TET0020] + a[1][TET0001]*tmp[TET0011];
  tmp[TET0001] = a[1][TET1000]*tmp[TET1001] + a[1][TET0100]*tmp[TET0101] + a[1][TET0010]*tmp[TET0011] + a[1][TET0001]*tmp[TET0002];

  return 6*(b[TET1000]*tmp[TET1000] + b[TET0100]*tmp[TET0100] + b[TET0010]*tmp[TET0010] + b[TET0001]*tmp[TET0001]);
}

bool bb33_convex_hull_brackets_value(bb33 const *bb, dbl value) {
  dbl min, max;
  dblN_minmax(bb->c, 20, &min, &max);
  return min <= value && value <= max;
}

/**
 * Restrict `bb` to the interval [`b0`, `b1`] returning the cubic
 * polynomial `p(t) = q((1 - t)*b0 + t*b1)`, where `q(b)` is the
 * trivariate polynomial represented by `bb`.
 */
cubic_s bb33_restrict_along_interval(bb33 const *bb, dbl b0[4], dbl b1[4]) {
  dbl db[4];
  dbl4_sub(b1, b0, db);
  dbl f[2] = {bb33_f(bb, b0), bb33_f(bb, b1)};
  dbl Df[2] = {bb33_df(bb, b0, db), bb33_df(bb, b1, db)};
  return cubic_from_data(f, Df);
}

// Local Variables:
// column-enforce-column: 160
// End:
