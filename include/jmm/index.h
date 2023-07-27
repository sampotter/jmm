#pragma once

#include "vec.h"

int ind2l(int2 const shape, int2 const ind);
int ind2lc(int2 const shape, int2 const ind);
int indc2l(int2 const shape, int2 const indc);
int indc2lc(int2 const shape, int2 const indc);
void l2ind(int2 const shape, int l, int2 ind);
void l2indc(int2 const shape, int l, int2 indc);
void lc2ind(int2 const shape, int lc, int2 ind);
void lc2indc(int2 const shape, int lc, int2 indc);
int l2lc(int2 const shape, int l);
int lc2l(int2 const shape, int lc);
int xy_to_lc_and_cc(int2 const shape, dbl2 const xymin, dbl h, dbl2 const xy,
                    dbl2 cc);

int ind2l3(int3 const shape, int3 const ind);
void l2ind3(int3 const shape, int l, int3 ind);
