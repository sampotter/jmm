#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "jet.h"

typedef struct mesh3 mesh3_s;

typedef struct utetra {
  dbl f;
  dbl g[2];
  dbl H[2][2];

  dbl p[2]; // Newton step

  dbl x[3]; // x[l]
  dbl X[3][3]; // X = [x[l0] x[l1] x[l2]]
  dbl Xt[3][3]; // X'
  dbl XtX[3][3]; // X'*X

  // B-coefs for 9-point triangle interpolation T on base of update
  dbl Tc[10];

  dbl x_minus_xb[3];

  dbl g_dot_p_active; // Dot product between g and p "restricted to
                      // the active set". If a constraint is active,
                      // we only want to "restrict this dot product"
                      // to the active subspace.

  int niter;
} utetra_s;

void utetra_init(utetra_s *cf, mesh3_s const *mesh, jet3 const *jet,
                 size_t l, size_t l0, size_t l1, size_t l2);
void utetra_set_lambda(utetra_s *cf, dbl const lam[2]);
void utetra_solve(utetra_s *cf, dbl lam[2], jet3 *jet);

#ifdef __cplusplus
}
#endif
