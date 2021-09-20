#include "utri21.h"

#include <assert.h>

#include "hybrid.h"
#include "mat.h"

void utri21_init(utri21_s *utri, dbl2 const xhat, dbl2 const x[2],
                 jet22t const jet[2]) {
  dbl2_copy(xhat, utri->xhat);
  dbl2_copy(x[0], utri->x[0]);
  dbl2_copy(x[1], utri->x[1]);

  dbl2_sub(utri->x[1], utri->x[0], utri->dx);

  bb31_init_from_jet22t(&utri->T, jet, utri->x);
}

dbl utri21_F(utri21_s const *wkspc, dbl lam) {
  dbl T = bb31_f(&wkspc->T, (dbl2) {1 - lam, lam});

  dbl2 xlam;
  dbl2_saxpy(lam, wkspc->dx, wkspc->x[0], xlam);

  dbl2 xhat_minus_xlam;
  dbl2_sub(wkspc->xhat, xlam, xhat_minus_xlam);

  dbl L = dbl2_norm(xhat_minus_xlam);

  return T + L;
}

dbl utri21_dF(utri21_s const *wkspc, dbl lam) {
  dbl T_lam = bb31_df(&wkspc->T, (dbl2) {1 - lam, lam}, (dbl2) {-1, 1});

  dbl2 xlam;
  dbl2_saxpy(lam, wkspc->dx, wkspc->x[0], xlam);

  dbl2 xhat_minus_xlam;
  dbl2_sub(wkspc->xhat, xlam, xhat_minus_xlam);

  dbl L = dbl2_norm(xhat_minus_xlam);
  dbl L_lam = -dbl2_dot(wkspc->dx, xhat_minus_xlam)/L;

  return T_lam + L_lam;
}

bool utri21_solve(utri21_s *utri, dbl *lam) {
  bool found = hybrid((hybrid_cost_func_t)utri21_dF, 0, 1, utri, lam);

  dbl T0 = utri21_F(utri, 0);
  dbl T1 = utri21_F(utri, 1);

  dbl T;
  if (found) {
    T = utri21_F(utri, *lam);
  } else {
    *lam = T0 < T1 ? 0 : 1;

    // Compute the Lagrange multiplier for the active constraint and
    // return if it's positive
    dbl mu = pow(-1, *lam)*utri21_dF(utri, *lam);
    if (mu > 0)
      return false;

    T = fmin(T0, T1);
  }

  assert(T > T0);
  assert(T > T1);

  jet22t *jet = &utri->jet;

  jet->f = T;

  dbl2 xlam;
  dbl2_saxpy(*lam, utri->dx, utri->x[0], xlam);
  dbl2_sub(utri->xhat, xlam, jet->Df);
  dbl L = dbl2_norm(jet->Df);
  dbl2_dbl_div_inplace(jet->Df, L);

  dbl22 eye, tt;
  dbl22_eye(eye);
  dbl2_outer(jet->Df, jet->Df, tt);
  dbl22_sub(eye, tt, jet->D2f);
  dbl22_dbl_div_inplace(jet->D2f, jet->f);

  return true;
}
