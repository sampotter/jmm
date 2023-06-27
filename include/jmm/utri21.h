#pragma once

#include "bb.h"

typedef struct utri21 {
  dbl2 xhat, x[2], dx;
  bb31 T;
  jet21t jet;
} utri21_s;

void utri21_init(utri21_s *utri, dbl2 const xhat, dbl2 const x[2],
                 jet21t const jet[2]);
dbl utri21_F(utri21_s const *wkspc, dbl lambda);
dbl utri21_dF(utri21_s const *wkspc, dbl lam);
bool utri21_solve(utri21_s *utri, dbl *lam);
