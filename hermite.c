#include "hermite.h"

dbl cubic_f(cubic const *cubic, dbl lam) {
  dbl const *a = cubic->a;
  return a[0] + lam*(a[1] + lam*(a[2] + lam*a[3]));
}

dbl cubic_df(cubic const *cubic, dbl lam) {
  dbl const *a = cubic->a;
  return a[1] + lam*(2*a[2] + 3*lam*a[3]);
}

dbl V_inv[4][4] = {
  { 1,  0,  0,  0},
  { 0,  0,  1,  0},
  {-3,  3, -2, -1},
  { 2, -2,  1,  1}
};

void bicubic_set_A(bicubic *bicubic, dbl data[4][4]) {
  dbl tmp[4][4];
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      tmp[i][j] = 0;
      for (int k = 0; k < 4; ++k) {
        tmp[i][j] += V_inv[i][k]*data[k][j];
      }
    }
  }

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      bicubic->A[i][j] = 0;
      for (int k = 0; k < 4; ++k) {
        bicubic->A[i][j] += tmp[i][k]*V_inv[i][k];
      }
    }
  }
}

cubic bicubic_restrict(bicubic const *bicubic, bicubic_variable var, int edge) {
  cubic cubic;
  if (var == LAMBDA) {
    if (edge == 0) {
      for (int alpha = 0; alpha < 4; ++alpha) {
        cubic.a[alpha] = bicubic->A[alpha][0];
      }
    } else {
      for (int alpha = 0; alpha < 4; ++alpha) {
        cubic.a[alpha] = 0;
        for (int beta = 0; beta < 4; ++beta) {
          cubic.a[alpha] += bicubic->A[alpha][beta];
        }
      }
    }
  } else {
    if (edge == 0) {
      for (int beta = 0; beta < 4; ++beta) {
        cubic.a[beta] = bicubic->A[0][beta];
      }
    } else {
      for (int beta = 0; beta < 4; ++beta) {
        cubic.a[beta] = 0;
        for (int alpha = 0; alpha < 4; ++alpha) {
          cubic.a[beta] += bicubic->A[alpha][beta];
        }
      }
    }
  }
  return cubic;
}

dbl bicubic_f(bicubic const *bicubic, dvec2 cc) {
  dbl const (*A)[4] = bicubic->A;
  dbl f = 0, lam = cc.x, mu = cc.y, lam_pow = 1, mu_pow;
  for (int alpha = 0; alpha < 4; ++alpha) {
    mu_pow = 1;
    for (int beta = 0; beta < 4; ++beta) {
      f += A[alpha][beta]*lam_pow*mu_pow;
      mu_pow *= mu;
    }
    lam_pow *= lam;
  }
  return f;
}

dbl bicubic_fx(bicubic const *bicubic, dvec2 cc) {
  dbl const (*A)[4] = bicubic->A;
  dbl fx = 0, lam = cc.x, mu = cc.y, lam_pow = 1, mu_pow;
  for (int alpha = 1; alpha < 4; ++alpha) {
    mu_pow = 1;
    for (int beta = 0; beta < 4; ++beta) {
      fx += alpha*A[alpha][beta]*lam_pow*mu_pow;
      mu_pow *= mu;
    }
    lam_pow *= lam;
  }
  return fx;
}

dbl bicubic_fy(bicubic const *bicubic, dvec2 cc) {
  dbl const (*A)[4] = bicubic->A;
  dbl fy = 0, lam = cc.x, mu = cc.y, lam_pow = 1, mu_pow;
  for (int alpha = 0; alpha < 4; ++alpha) {
    mu_pow = 1;
    for (int beta = 1; beta < 4; ++beta) {
      fy += beta*A[alpha][beta]*lam_pow*mu_pow;
      mu_pow *= mu;
    }
    lam_pow *= lam;
  }
  return fy;
}

dbl bicubic_fxy(bicubic const *bicubic, dvec2 cc) {
  dbl const (*A)[4] = bicubic->A;
  dbl fxy = 0, lam = cc.x, mu = cc.y, lam_pow = 1, mu_pow;
  for (int alpha = 1; alpha < 4; ++alpha) {
    mu_pow = 1;
    for (int beta = 1; beta < 4; ++beta) {
      fxy += alpha*beta*A[alpha][beta]*lam_pow*mu_pow;
      mu_pow *= mu;
    }
    lam_pow *= lam;
  }
  return fxy;
}
