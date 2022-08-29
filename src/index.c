#include <jmm/index.h>

#include <assert.h>
#include <stddef.h>

#include <jmm/def.h>

int ind2l(int2 const shape, int2 const ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  return ind[1] + shape[1]*ind[0];
#else
  return shape[0]*ind[1] + ind[0];
#endif
}

int ind2lc(int2 const shape, int2 const ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  return ind[1] + (shape[1] - 1)*ind[0];
#else
  return (shape[0] - 1)*ind[1] + ind[0];
#endif
}

int indc2l(int2 const shape, int2 const indc) {
#if ORDERING == ROW_MAJOR_ORDERING
  return indc[1] + shape[1]*indc[0];
#else
  return shape[0]*indc[1] + indc[0];
#endif
}

int indc2lc(int2 const shape, int2 const indc) {
#if ORDERING == ROW_MAJOR_ORDERING
  return indc[1] + (shape[1] - 1)*indc[0];
#else
  return (shape[0] - 1)*indc[1] + indc[0];
#endif
}

void l2ind(int2 const shape, int l, int2 ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  ind[0] = l/shape[1];
  ind[1] = l % shape[1];
#else
  ind[0] = l % shape[0];
  int[1] = l/shape[0];
#endif
}

void l2indc(int2 const shape, int l, int2 indc) {
#if ORDERING == ROW_MAJOR_ORDERING
  indc[0] = l/shape[1];
  indc[1] = l % shape[1];
#else
  indc[0] = l % shape[0];
  indc[1] = l/shape[0];
#endif
}

void lc2ind(int2 const shape, int lc, int2 ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  ind[0] = lc/(shape[1] - 1);
  ind[1] = lc % (shape[1] - 1);
#else
  ind[0] = lc % (shape[0] - 1);
  ind[1] = lc/(shape[0] - 1);
#endif
}

void lc2indc(int2 const shape, int lc, int2 indc) {
#if ORDERING == ROW_MAJOR_ORDERING
  indc[0] = lc/(shape[1] - 1);
  indc[1] = lc % (shape[1] - 1);
#else
  indc[0] = lc % (shape[0] - 1);
  indc[1] = lc/(shape[0] - 1);
#endif
}

int l2lc(int2 const shape, int l) {
#if ORDERING == ROW_MAJOR_ORDERING
  return l - l/shape[1];
#else
  return l - l/shape[0];
#endif
}

int lc2l(int2 const shape, int lc) {
#if ORDERING == ROW_MAJOR_ORDERING
  return lc + lc/(shape[1] - 1);
#else
  return lc + lc/(shape[0] - 1);
#endif
}

int xy_to_lc_and_cc(int2 const shape, dbl2 const xymin, dbl h, dbl2 const xy,
                    dbl2 cc) {
#if SJS_DEBUG
  assert(cc != NULL);
#endif

  dbl2_sub(xy, xymin, cc);
  dbl2_dbl_div_inplace(cc, h);
  dbl2 ind_; dbl2_floor(cc, ind_);
  dbl2_sub_inplace(cc, ind_);
  int2 ind = {ind_[0], ind_[1]};

  if (ind[0] < 0) {
    ind[0] = 0;
    cc[0] = 0.0;
  }

  if (ind[1] < 0) {
    ind[1] = 0;
    cc[1] = 0.0;
  }

  if (ind[0] >= shape[0] - 1) {
    --ind[0];
    cc[0] = 1.0;
  }

  if (ind[1] >= shape[1] - 1) {
    --ind[1];
    cc[1] = 1.0;
  }

  return ind2lc(shape, ind);
}

int ind2l3(int3 const shape, int3 const ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  return ind[2] + shape[2]*(ind[1] + shape[1]*ind[0]);
#else
#  error not implemented yet
#endif
}

void l2ind3(int3 const shape, int l, int3 ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  ind[0] = l/(shape[2]*shape[1]);
  ind[1] = l/shape[2] % shape[1];
  ind[2] = l % shape[2];
#else
#  error not implemented yet
#endif
}
