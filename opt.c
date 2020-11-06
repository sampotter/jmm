#include "opt.h"

#include "vec.h"

bool baryopt_step(baryopt_wkspc_s *wkspc) {
    // ik = np.argmax(xk)
  int i = dvec3_argmax(*wkspc->x);

    // Tk = np.eye(xk.size)
    // Tk[ik, :] = 1
    // Tk = np.linalg.inv(Tk)
  dmat33 T = dmat33_eye();
  for (int j = 0; j < 3; ++j) {
    if (i != j) T.rows[i].data[j] = -1;
  }

    // yk = np.linalg.solve(Tk, xk)
  dvec3 y = *wkspc->x;
  y.data[i] = 1;

    // fk = F(xk)
    // gk = Tk.T@DF(xk)
    // Hk = Tk.T@D2F(xk)@Tk
  dvec3 g;
  dmat33 H;
  wkspc->costfunc->func(wkspc->x, &wkspc->f, wkspc->costfunc->wkspc);
  wkspc->costfunc->grad(wkspc->x, &g, wkspc->costfunc->wkspc);
  wkspc->costfunc->hess(wkspc->x, &H, wkspc->costfunc->wkspc);

    // pk = np.maximum(0, yk - gk/np.diag(Hk))
  dvec3 p = y;
  for (int j = 0; j < 3; ++j) {
    p.data[j] = fmax(0, p.data[j] - g.data[j]/H.rows[j].data[j]);
  }

    // I = np.setdiff1d(np.arange(xk.size), [ik])
    // wk = np.linalg.norm((yk - pk)[I])
  dbl w = 0;
  for (int j = 0; j < 3; ++j) {
    if (i == j) continue;
    w += (y.data[j] - p.data[j])*(y.data[j] - p.data[j]);
  }
  w = sqrt(w);

    // epsk = min(eps, wk)
  dbl eps = fmin(wkspc->eps, w);

    // Ikp = np.where((0 <= yk) & (yk <= epsk) & (gk > 0))[0]
  bool I[3];
  for (int j = 0; j < 3; ++j) {
    I[j] = 0 <= y.data[j] && y.data[j] <= eps && g.data[j] > 0;
  }

    // for i in Ikp:
    //     for j in Ikp:
    //         if i == j:
    //             continue
    //         Hk[i, j] = 0
  for (int j = 0; j < 3; ++j) {
    for (int k = 0; k < 3; ++k) {
      if (j == k) continue;
      H.rows[j].data[k] = 0;
    }
  }

    // pk = np.linalg.solve(Hk, gk)
  p = dmat33_dvec3_solve(H, g);

    // Ikpc = np.setdiff1d(np.arange(xk.size), Ikp)

  dbl f1;
  dvec3 y1, x1;

  bool converged = false;

    // m = 0

  dbl alpha = 1, lhs, rhs;
  for (int m = 0; ; ++m) {

    // def get_yk1(alpha):
    //     yk1 = np.maximum(0, yk - alpha*pk)
    //     yk1[ik] = 1
    //     return yk1
    // yk1 = get_yk1(beta**m)
    y1 = dvec3_sub(y, dvec3_dbl_mul(p, alpha));
    for (int j = 0; j < 3; ++j) {
      y1.data[j] = fmax(0, y1.data[j]);
    }

    x1 = y1;
    for (int j = 0; j < 3; ++j) {
      if (i != j) x1.data[i] -= y1.data[j];
    }

    for (int j = 0; j < 3; ++j) {
      if (x1.data[j] < 0) {
        goto step;
      }
    }

    // fk1 = F(Tk@yk1)
    wkspc->costfunc->func(&x1, &f1, wkspc->costfunc->wkspc);

    if (dvec3_maxdist(x1, *wkspc->x) <= wkspc->xtol) {
      converged = true;
      break;
    }

    // lhs = fk - fk1
    lhs = wkspc->f - f1;

    if (fabs(lhs)/fabs(wkspc->f) <= wkspc->ftol) {
      converged = true;
      break;
    }

    // rhs = sigma*(gk[Ikpc]@pk[Ikpc]*beta**m + gk[Ikp]@(yk - yk1)[Ikp])
    rhs = 0;
    for (int j = 0; j < 3; ++j) {
      rhs += I[j] ?
        g.data[j]*(y.data[j] - y1.data[j]) :
        alpha*g.data[j]*p.data[j];
    }
    rhs *= wkspc->sigma;

    // while abs(fk - fk1)/abs(fk) > np.finfo(np.dtype(fk)).eps and \
    //       ((Tk@yk1 < 0).any() or lhs < rhs):
    if (lhs >= rhs) {
      break;
    }

    //     m += 1
    //     yk1 = get_yk1(beta**m)
    //     fk1 = F(Tk@yk1)
    //     lhs = fk - fk1
    //     rhs = sigma*(gk[Ikpc]@pk[Ikpc]*beta**m + gk[Ikp]@(yk - yk1)[Ikp])
  step:
    alpha *= wkspc->beta;
  }

  *wkspc->x = x1;
  wkspc->f = f1;

    // xk1 = Tk@yk1

    // return xk1, m

  return converged;
}

void baryopt_solve(baryopt_wkspc_s *wkspc) {
  while (!baryopt_step(wkspc)) {}
}
