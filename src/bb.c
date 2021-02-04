#include "bb.h"

#include "mesh3.h"
#include "vec.h"

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

void bb3_interp(dbl const f[2], dbl const Df[2], dbl const x[2], dbl c[4]) {
  c[1] = (c[0] = f[0]) + Df[0]*(x[1] - x[0])/3;
  c[2] = (c[3] = f[1]) + Df[1]*(x[0] - x[1])/3;
}

void bb3_interp3(dbl const f[2], dbl const Df[2][3], dbl const x[2][3], dbl c[4]) {
  dbl dx[3];
  dbl3_sub(x[1], x[0], dx);
  c[1] = (c[0] = f[0]) + dbl3_dot(dx, Df[0])/3;
  c[2] = (c[3] = f[1]) - dbl3_dot(dx, Df[1])/3;
}

dbl bb3(dbl const *c, dbl const *b) {
  dbl tmp[3];

  tmp[0] = b[0]*c[0] + b[1]*c[1];
  tmp[1] = b[0]*c[1] + b[1]*c[2];
  tmp[2] = b[0]*c[2] + b[1]*c[3];

  tmp[0] = b[0]*tmp[0] + b[1]*tmp[1];
  tmp[1] = b[0]*tmp[1] + b[1]*tmp[2];

  return b[0]*tmp[0] + b[1]*tmp[1];
}

dbl dbb3(dbl const *c, dbl const *b, dbl const *a) {
  dbl tmp[3];

  tmp[0] = a[0]*c[0] + a[1]*c[1];
  tmp[1] = a[0]*c[1] + a[1]*c[2];
  tmp[2] = a[0]*c[2] + a[1]*c[3];

  tmp[0] = b[0]*tmp[0] + b[1]*tmp[1];
  tmp[1] = b[0]*tmp[1] + b[1]*tmp[2];

  return 3*(b[0]*tmp[0] + b[1]*tmp[1]);
}

dbl d2bb3(dbl const *c, dbl const *b, dbl const *a) {
  dbl tmp[3];

  tmp[0] = a[0]*c[0] + a[1]*c[1];
  tmp[1] = a[0]*c[1] + a[1]*c[2];
  tmp[2] = a[0]*c[2] + a[1]*c[3];

  tmp[0] = a[0]*tmp[0] + a[1]*tmp[1];
  tmp[1] = a[0]*tmp[1] + a[1]*tmp[2];

  return 3*(b[0]*tmp[0] + b[1]*tmp[1]);
}

void bb3tri_interp3(dbl const f[3], dbl const Df[3][3], dbl const x[3][3], dbl c[10]) {
  c[TRI300] = f[0];
  c[TRI030] = f[1];
  c[TRI003] = f[2];

  dbl tmp[3];

  dbl3_sub(x[TRI010], x[TRI100], tmp);
  c[TRI210] = c[TRI300] + dbl3_dot(Df[TRI100], tmp)/3;
  c[TRI120] = c[TRI030] - dbl3_dot(Df[TRI010], tmp)/3;

  dbl3_sub(x[TRI001], x[TRI100], tmp);
  c[TRI201] = c[TRI300] + dbl3_dot(Df[TRI100], tmp)/3;
  c[TRI102] = c[TRI003] - dbl3_dot(Df[TRI001], tmp)/3;

  dbl3_sub(x[TRI001], x[TRI010], tmp);
  c[TRI021] = c[TRI030] + dbl3_dot(Df[TRI010], tmp)/3;
  c[TRI012] = c[TRI003] - dbl3_dot(Df[TRI001], tmp)/3;

  c[TRI111] = (c[TRI210] + c[TRI201] + c[TRI120] + c[TRI021] + c[TRI012] + c[TRI102])/4
    - (c[TRI300] + c[TRI030] + c[TRI003])/6;
}

dbl bb3tri(dbl const *c, dbl const *b) {
  dbl tmp[6];

  tmp[TRI200] = b[TRI100]*c[TRI300] + b[TRI010]*c[TRI210] + b[TRI001]*c[TRI201];
  tmp[TRI110] = b[TRI100]*c[TRI210] + b[TRI010]*c[TRI120] + b[TRI001]*c[TRI111];
  tmp[TRI020] = b[TRI100]*c[TRI120] + b[TRI010]*c[TRI030] + b[TRI001]*c[TRI021];
  tmp[TRI101] = b[TRI100]*c[TRI201] + b[TRI010]*c[TRI111] + b[TRI001]*c[TRI102];
  tmp[TRI011] = b[TRI100]*c[TRI111] + b[TRI010]*c[TRI021] + b[TRI001]*c[TRI012];
  tmp[TRI002] = b[TRI100]*c[TRI102] + b[TRI010]*c[TRI012] + b[TRI001]*c[TRI003];

  tmp[TRI100] = b[TRI100]*tmp[TRI200] + b[TRI010]*tmp[TRI110] + b[TRI001]*tmp[TRI101];
  tmp[TRI010] = b[TRI100]*tmp[TRI110] + b[TRI010]*tmp[TRI020] + b[TRI001]*tmp[TRI011];
  tmp[TRI001] = b[TRI100]*tmp[TRI101] + b[TRI010]*tmp[TRI011] + b[TRI001]*tmp[TRI002];

  return b[TRI100]*tmp[TRI100] + b[TRI010]*tmp[TRI010] + b[TRI001]*tmp[TRI001];
}

dbl dbb3tri(dbl const *c, dbl const *b, dbl const *a) {
  dbl tmp[6];

  tmp[TRI200] = a[TRI100]*c[TRI300] + a[TRI010]*c[TRI210] + a[TRI001]*c[TRI201];
  tmp[TRI110] = a[TRI100]*c[TRI210] + a[TRI010]*c[TRI120] + a[TRI001]*c[TRI111];
  tmp[TRI020] = a[TRI100]*c[TRI120] + a[TRI010]*c[TRI030] + a[TRI001]*c[TRI021];
  tmp[TRI101] = a[TRI100]*c[TRI201] + a[TRI010]*c[TRI111] + a[TRI001]*c[TRI102];
  tmp[TRI011] = a[TRI100]*c[TRI111] + a[TRI010]*c[TRI021] + a[TRI001]*c[TRI012];
  tmp[TRI002] = a[TRI100]*c[TRI102] + a[TRI010]*c[TRI012] + a[TRI001]*c[TRI003];

  tmp[TRI100] = b[TRI100]*tmp[TRI200] + b[TRI010]*tmp[TRI110] + b[TRI001]*tmp[TRI101];
  tmp[TRI010] = b[TRI100]*tmp[TRI110] + b[TRI010]*tmp[TRI020] + b[TRI001]*tmp[TRI011];
  tmp[TRI001] = b[TRI100]*tmp[TRI101] + b[TRI010]*tmp[TRI011] + b[TRI001]*tmp[TRI002];

  return 3*(b[TRI100]*tmp[TRI100] + b[TRI010]*tmp[TRI010] + b[TRI001]*tmp[TRI001]);
}

dbl d2bb3tri(dbl const *c, dbl const *b, dbl const *a1, dbl const *a2) {
  dbl tmp[6];

  tmp[TRI200] = a1[TRI100]*c[TRI300] + a1[TRI010]*c[TRI210] + a1[TRI001]*c[TRI201];
  tmp[TRI110] = a1[TRI100]*c[TRI210] + a1[TRI010]*c[TRI120] + a1[TRI001]*c[TRI111];
  tmp[TRI020] = a1[TRI100]*c[TRI120] + a1[TRI010]*c[TRI030] + a1[TRI001]*c[TRI021];
  tmp[TRI101] = a1[TRI100]*c[TRI201] + a1[TRI010]*c[TRI111] + a1[TRI001]*c[TRI102];
  tmp[TRI011] = a1[TRI100]*c[TRI111] + a1[TRI010]*c[TRI021] + a1[TRI001]*c[TRI012];
  tmp[TRI002] = a1[TRI100]*c[TRI102] + a1[TRI010]*c[TRI012] + a1[TRI001]*c[TRI003];

  tmp[TRI100] = a2[TRI100]*tmp[TRI200] + a2[TRI010]*tmp[TRI110] + a2[TRI001]*tmp[TRI101];
  tmp[TRI010] = a2[TRI100]*tmp[TRI110] + a2[TRI010]*tmp[TRI020] + a2[TRI001]*tmp[TRI011];
  tmp[TRI001] = a2[TRI100]*tmp[TRI101] + a2[TRI010]*tmp[TRI011] + a2[TRI001]*tmp[TRI002];

  return 6*(b[TRI100]*tmp[TRI100] + b[TRI010]*tmp[TRI010] + b[TRI001]*tmp[TRI001]);
}

void bb3tet_interp3(dbl const f[4], dbl const Df[4][3], dbl const x[4][3], dbl c[20]) {
  dbl dx[3];

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
  c[TET2100] = c[TET3000] + dbl3_dot(Df[TET1000], dx)/3;
  c[TET1200] = c[TET0300] - dbl3_dot(Df[TET0100], dx)/3;

  // 1 <-> 3
  dbl3_sub(x[TET0010], x[TET1000], dx);
  c[TET2010] = c[TET3000] + dbl3_dot(Df[TET1000], dx)/3;
  c[TET1020] = c[TET0030] - dbl3_dot(Df[TET0010], dx)/3;

  // 1 <-> 4
  dbl3_sub(x[TET0001], x[TET1000], dx);
  c[TET2001] = c[TET3000] + dbl3_dot(Df[TET1000], dx)/3;
  c[TET1002] = c[TET0003] - dbl3_dot(Df[TET0001], dx)/3;

  // 2 <-> 3
  dbl3_sub(x[TET0010], x[TET0100], dx);
  c[TET0210] = c[TET0300] + dbl3_dot(Df[TET0100], dx)/3;
  c[TET0120] = c[TET0030] - dbl3_dot(Df[TET0010], dx)/3;

  // 2 <-> 4
  dbl3_sub(x[TET0001], x[TET0100], dx);
  c[TET0201] = c[TET0300] + dbl3_dot(Df[TET0100], dx)/3;
  c[TET0102] = c[TET0003] - dbl3_dot(Df[TET0001], dx)/3;

  // 3 <-> 4
  dbl3_sub(x[TET0001], x[TET0010], dx);
  c[TET0021] = c[TET0030] + dbl3_dot(Df[TET0010], dx)/3;
  c[TET0012] = c[TET0003] - dbl3_dot(Df[TET0001], dx)/3;

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
  // generalzation of the 9-parameter interpolant.)
}

void bb3tet_for_cell(mesh3_s const *mesh, jet3 const *jet, size_t lc,
                     dbl c[20]) {
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
  bb3tet_interp3(y, Dy, x, c);
}

dbl bb3tet(dbl const c[20], dbl const b[4]) {
  dbl tmp[10];

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

dbl dbb3tet(dbl const c[20], dbl const b[4], dbl const a[4]) {
  dbl tmp[10];

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

dbl d2bb3tet(dbl const c[20], dbl const b[4], dbl const a[2][4]) {
  dbl tmp[10];

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

// Local Variables:
// column-enforce-column: 160
// End:
