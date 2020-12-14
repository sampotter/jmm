#include "bb.h"

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

  tmp[0] = b[0]*c[0] + b[1]*c[1];
  tmp[1] = b[0]*c[1] + b[1]*c[2];
  tmp[2] = b[0]*c[2] + b[1]*c[3];

  tmp[0] = b[0]*tmp[0] + b[1]*tmp[1];
  tmp[1] = b[0]*tmp[1] + b[1]*tmp[2];

  return a[0]*tmp[0] + a[1]*tmp[1];
}

dbl d2bb3(dbl const *c, dbl const *b, dbl const *a) {
  dbl tmp[3];

  tmp[0] = b[0]*c[0] + b[1]*c[1];
  tmp[1] = b[0]*c[1] + b[1]*c[2];
  tmp[2] = b[0]*c[2] + b[1]*c[3];

  tmp[0] = a[0]*tmp[0] + a[1]*tmp[1];
  tmp[1] = a[0]*tmp[1] + a[1]*tmp[2];

  return a[0]*tmp[0] + a[1]*tmp[1];
}

void bb3tri_interp3(dbl f[3], dbl Df[3][3], dbl x[3][3], dbl *c) {
  c[TRI300] = f[0];
  c[TRI030] = f[1];
  c[TRI003] = f[2];

  dbl tmp[3];

  dbl3_sub(x[TRI010], x[TRI100], tmp);
  c[TRI210] = c[TRI300] + dbl3_dot(Df[TRI100], tmp)/3;

  dbl3_sub(x[TRI001], x[TRI100], tmp);
  c[TRI201] = c[TRI300] + dbl3_dot(Df[TRI100], tmp)/3;

  dbl3_sub(x[TRI100], x[TRI010], tmp);
  c[TRI120] = c[TRI030] + dbl3_dot(Df[TRI010], tmp)/3;

  dbl3_sub(x[TRI001], x[TRI010], tmp);
  c[TRI021] = c[TRI030] + dbl3_dot(Df[TRI010], tmp)/3;

  dbl3_sub(x[TRI010], x[TRI001], tmp);
  c[TRI012] = c[TRI003] + dbl3_dot(Df[TRI001], tmp)/3;

  dbl3_sub(x[TRI100], x[TRI001], tmp);
  c[TRI102] = c[TRI003] + dbl3_dot(Df[TRI001], tmp)/3;

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

  tmp[TRI200] = b[TRI100]*c[TRI300] + b[TRI010]*c[TRI210] + b[TRI001]*c[TRI201];
  tmp[TRI110] = b[TRI100]*c[TRI210] + b[TRI010]*c[TRI120] + b[TRI001]*c[TRI111];
  tmp[TRI020] = b[TRI100]*c[TRI120] + b[TRI010]*c[TRI030] + b[TRI001]*c[TRI021];
  tmp[TRI101] = b[TRI100]*c[TRI201] + b[TRI010]*c[TRI111] + b[TRI001]*c[TRI102];
  tmp[TRI011] = b[TRI100]*c[TRI111] + b[TRI010]*c[TRI021] + b[TRI001]*c[TRI012];
  tmp[TRI002] = b[TRI100]*c[TRI102] + b[TRI010]*c[TRI012] + b[TRI001]*c[TRI003];

  tmp[TRI100] = b[TRI100]*tmp[TRI200] + b[TRI010]*tmp[TRI110] + b[TRI001]*tmp[TRI101];
  tmp[TRI010] = b[TRI100]*tmp[TRI110] + b[TRI010]*tmp[TRI020] + b[TRI001]*tmp[TRI011];
  tmp[TRI001] = b[TRI100]*tmp[TRI101] + b[TRI010]*tmp[TRI011] + b[TRI001]*tmp[TRI002];

  return a[TRI100]*tmp[TRI100] + a[TRI010]*tmp[TRI010] + a[TRI001]*tmp[TRI001];
}

dbl d2bb3tri(dbl const *c, dbl const *b, dbl const *a1, dbl const *a2) {
  dbl tmp[6];

  tmp[TRI200] = b[TRI100]*c[TRI300] + b[TRI010]*c[TRI210] + b[TRI001]*c[TRI201];
  tmp[TRI110] = b[TRI100]*c[TRI210] + b[TRI010]*c[TRI120] + b[TRI001]*c[TRI111];
  tmp[TRI020] = b[TRI100]*c[TRI120] + b[TRI010]*c[TRI030] + b[TRI001]*c[TRI021];
  tmp[TRI101] = b[TRI100]*c[TRI201] + b[TRI010]*c[TRI111] + b[TRI001]*c[TRI102];
  tmp[TRI011] = b[TRI100]*c[TRI111] + b[TRI010]*c[TRI021] + b[TRI001]*c[TRI012];
  tmp[TRI002] = b[TRI100]*c[TRI102] + b[TRI010]*c[TRI012] + b[TRI001]*c[TRI003];

  tmp[TRI100] = a1[TRI100]*tmp[TRI200] + a1[TRI010]*tmp[TRI110] + a1[TRI001]*tmp[TRI101];
  tmp[TRI010] = a1[TRI100]*tmp[TRI110] + a1[TRI010]*tmp[TRI020] + a1[TRI001]*tmp[TRI011];
  tmp[TRI001] = a1[TRI100]*tmp[TRI101] + a1[TRI010]*tmp[TRI011] + a1[TRI001]*tmp[TRI002];

  return a2[TRI100]*tmp[TRI100] + a2[TRI010]*tmp[TRI010] + a2[TRI001]*tmp[TRI001];
}

dbl bb3tet(dbl const *c, dbl const *b) {
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

// Local Variables:
// column-enforce-column: 160
// End:
