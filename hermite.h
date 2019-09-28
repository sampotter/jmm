#pragma once

#include "def.h"

typedef enum {LAMBDA, MU} bicubic_variable;

typedef struct {
  dbl a[4];
} cubic;

dbl cubic_f(cubic *cubic, dbl lam);
dbl cubic_df(cubic *cubic, dbl lam);

typedef struct {
  dbl A[4][4];
} bicubic;

void bicubic_set_A(bicubic *bicubic, dbl data[4][4]);
cubic bicubic_restrict(bicubic *bicubic, bicubic_variable var, int edge);
