#include "index.h"

#include "def.h"

#include <assert.h>
#include <stddef.h>

int ind2l(ivec2 shape, ivec2 ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  return ind.j + shape.j*ind.i;
#else
  return shape.i*ind.j + ind.i;
#endif
}

int ind2lc(ivec2 shape, ivec2 ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  return ind.j + (shape.j - 1)*ind.i;
#else
  return (shape.i - 1)*ind.j + ind.i;
#endif
}

int indc2l(ivec2 shape, ivec2 indc) {
#if ORDERING == ROW_MAJOR_ORDERING
  return indc.j + shape.j*indc.i;
#else
  return shape.i*indc.j + indc.i;
#endif
}

int indc2lc(ivec2 shape, ivec2 indc) {
#if ORDERING == ROW_MAJOR_ORDERING
  return indc.j + (shape.j - 1)*indc.i;
#else
  return (shape.i - 1)*indc.j + indc.i;
#endif
}

ivec2 l2ind(ivec2 shape, int l) {
#if ORDERING == ROW_MAJOR_ORDERING
  ivec2 ind = {.i = l/shape.j, .j = l % shape.j};
#else
  ivec2 ind = {.i = l % shape.i, .j = l/shape.i};
#endif
  return ind;
}

ivec2 l2indc(ivec2 shape, int l) {
#if ORDERING == ROW_MAJOR_ORDERING
  ivec2 indc = {.i = l/shape.j, .j = l % shape.j};
#else
  ivec2 indc = {.i = l % shape.i, .j = l/shape.i};
#endif
  return indc;
}

ivec2 lc2ind(ivec2 shape, int lc) {
#if ORDERING == ROW_MAJOR_ORDERING
  ivec2 ind = {.i = lc/(shape.j - 1), .j = lc % (shape.j - 1)};
#else
  ivec2 ind = {.i = lc % (shape.i - 1), .j = lc/(shape.i - 1)};
#endif
  return ind;
}

ivec2 lc2indc(ivec2 shape, int lc) {
#if ORDERING == ROW_MAJOR_ORDERING
  ivec2 indc = {.i = lc/(shape.j - 1), .j = lc % (shape.j - 1)};
#else
  ivec2 indc = {.i = lc % (shape.i - 1), .j = lc/(shape.i - 1)};
#endif
  return indc;
}

int l2lc(ivec2 shape, int l) {
#if ORDERING == ROW_MAJOR_ORDERING
  return l - l/shape.j;
#else
  return l - l/shape.i;
#endif
}

int lc2l(ivec2 shape, int lc) {
#if ORDERING == ROW_MAJOR_ORDERING
  return lc + lc/(shape.j - 1);
#else
  return lc + lc/(shape.i - 1);
#endif
}

int xy_to_lc_and_cc(ivec2 shape, dvec2 xymin, dbl h, dvec2 xy, dvec2 *cc) {
#if SJS_DEBUG
  assert(cc != NULL);
#endif

  *cc = dvec2_sub(xy, xymin);
  *cc = dvec2_dbl_div(*cc, h);
  dvec2 ind_ = dvec2_floor(*cc);
  *cc = dvec2_sub(*cc, ind_);
  ivec2 ind = dvec2_to_ivec2(ind_);

  if (ind.i < 0) {
    ind.i = 0;
    cc->x = 0.0;
  }

  if (ind.j < 0) {
    ind.j = 0;
    cc->y = 0.0;
  }

  if (ind.i >= shape.i - 1) {
    --ind.i;
    cc->x = 1.0;
  }

  if (ind.j >= shape.j - 1) {
    --ind.j;
    cc->y = 1.0;
  }

  return ind2lc(shape, ind);
}

int ind2l3(ivec3 shape, ivec3 ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  return ind.data[2] + shape.data[2]*(ind.data[1] + shape.data[1]*ind.data[0]);
#else
#  error not implemented yet
#endif
}

ivec3 l2ind3(ivec3 shape, int l) {
#if ORDERING == ROW_MAJOR_ORDERING
  return (ivec3) {
    .data = {
      l/(shape.data[2]*shape.data[1]),
      l/shape.data[2] % shape.data[1],
      l % shape.data[2]
    }
  };
#else
#  error not implemented yet
#endif
}
