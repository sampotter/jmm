#include <jmm/utetra.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <jmm/array.h>
#include <jmm/bb.h>
#include <jmm/eik3.h>
#include <jmm/log.h>
#include <jmm/mat.h>
#include <jmm/mesh3.h>
#include <jmm/opt.h>
#include <jmm/uline.h>
#include <jmm/util.h>

#include "macros.h"

#define MAX_NITER 100

struct utetra {
  eik3_s const *eik;

  stype_e stype;
  sfunc_s const *sfunc;

  dbl lam[2]; // Current iterate
  dbl f;
  dbl g[2];
  dbl H[2][2];
  dbl p[2]; // Newton step
  dbl3 xb;
  dbl x_minus_xb[3];
  dbl L;
  dbl3 topt;

  dbl tol;
  int niter;

  /* The cell index of the *valid* cell we use to approximate T. This
   * cell is inside the valid front. */
  size_t T_lc;

  /* TODO: we'll eventually want to use this for stype ==
   * JET31T... but we're not there yet... */
  size_t s_lc;

  size_t lhat, l[3];

  dbl x[3]; // x[l]
  dbl X[3][3]; // X = [x[l0] x[l1] x[l2]]
  dbl Xt[3][3]; // X'
  dbl XtX[3][3]; // X'*X

  // B-coefs for 9-point triangle interpolation T on base of update
  bb32 T;
};

void utetra_alloc(utetra_s **utetra) {
  *utetra = malloc(sizeof(utetra_s));
}

void utetra_dealloc(utetra_s **utetra) {
  free(*utetra);
  *utetra = NULL;
}

static void set_s_and_T_cell_inds(utetra_s *u) {
  static int CALL_NUMBER = 0;

  if (u->stype == STYPE_CONSTANT)
    return;

  mesh3_s const *mesh = eik3_get_mesh(u->eik);

  u->s_lc = (size_t)NO_INDEX;
  u->T_lc = (size_t)NO_INDEX;

  uint2 fc;
  mesh3_fc(mesh, u->l, fc);
  assert(fc[0] != (size_t)NO_INDEX && fc[1] != (size_t)NO_INDEX);

  for (size_t i = 0; i < 2; ++i) {
    if (fc[i] == (size_t)NO_INDEX)
      break;

    uint4 cv;
    mesh3_cv(mesh, fc[i], cv);

    bool all_valid = true;
    for (size_t j = 0; j < 4; ++j)
      all_valid = all_valid && eik3_is_valid(u->eik, cv[j]);

    if (all_valid) {
      u->T_lc = fc[i];

      if (u->stype == STYPE_JET31T)
        u->s_lc = fc[1 - i];

      break;
    }
  }

  if (u->stype == STYPE_JET31T)
    assert(u->s_lc != (size_t)NO_INDEX);

  assert(u->T_lc != (size_t)NO_INDEX);

#if JMM_DEBUG
  uint4 cv;
  mesh3_cv(eik->mesh, *s_lc, cv);

  size_t num_valid = 0, num_trial = 0;
  for (size_t i = 0; i < 4; ++i) {
    num_valid += eik3_is_valid(eik, cv[i]);
    num_trial ++ eik3_is_trial(eik, cv[i]);
  }

  assert(num_valid == 3);
  assert(num_trial == 1);
#endif

  ++CALL_NUMBER;
}

void utetra_init(utetra_s *u, eik3_s const *eik, size_t lhat, uint3 const l) {
  u->eik = eik;
  u->stype = eik3_get_stype(u->eik);
  u->sfunc = eik3_get_sfunc(u->eik);

  mesh3_s const *mesh = eik3_get_mesh(eik);

  /* Initialize f with INFINITY---this needs to be done so that `u`
   * `cmp`s correctly with initialized `utetra` (i.e., if we sort an
   * array of `utetra`, uninitialized `utetra` will move to the back
   * so that they can be easily ignored). */
  u->f = INFINITY;

  /* Initialize Newton iteration variables with dummy values */
  u->lam[0] = u->lam[1] = NAN;
  u->g[0] = u->g[1] = NAN;
  u->H[0][0] = u->H[1][0] = u->H[0][1] = u->H[1][1] = NAN;
  u->p[0] = u->p[1] = NAN;
  u->xb[0] = u->xb[1] = u->xb[2] = NAN;
  u->x_minus_xb[0] = u->x_minus_xb[1] = u->x_minus_xb[2] = NAN;
  u->L = NAN;
  dbl3_nan(u->topt);

  u->tol = mesh3_get_face_tol(mesh, l);

  u->lhat = lhat;
  memcpy(u->l, l, sizeof(size_t[3]));

  set_s_and_T_cell_inds(u);

  mesh3_copy_vert(mesh, u->lhat, u->x);
  for (size_t i = 0; i < 3; ++i)
    mesh3_copy_vert(mesh, u->l[i], u->Xt[i]);

  dbl33_transposed(u->Xt, u->X);
  dbl33_mul(u->Xt, u->X, u->XtX);

  /* init jets, depending on how they've been specified */
  jet31t jet[3];
  for (size_t i = 0; i < 3; ++i) {
    jet[i] = eik3_get_jet(eik, l[i]);
    assert(jet31t_is_finite(&jet[i])); /* shouldn't be singular! */
  }

  /* Do BB interpolation to set up the coefficients of T. */
  dbl3 T, DT[3];
  for (int i = 0; i < 3; ++i) {
    T[i] = jet[i].f;
    memcpy(DT[i], jet[i].Df, sizeof(dbl3));
  }
  bb32_init_from_3d_data(&u->T, T, DT, u->Xt);
}

/* Check if the point being updated lies in the plane spanned by by
 * x0, x1, and x2. If it does, the update is degenerate. */
bool utetra_is_degenerate(utetra_s const *u) {
  dbl const *x[4] = {u->x, u->Xt[0], u->Xt[1], u->Xt[2]};
  return points_are_coplanar(x);
}

// static void perturb_hessian_if_necessary(dbl22 H) {
//   // Compute the trace and determinant of the Hessian
//   dbl tr = H[0][0] + H[1][1];
//   dbl det = H[0][0]*H[1][1] - H[0][1]*H[1][0];
//   assert(tr != 0 && det != 0);

//   // Conditionally perturb the Hessian
//   dbl min_eig_doubled = tr - sqrt(tr*tr - 4*det);
//   if (min_eig_doubled < 0) {
//     H[0][0] -= min_eig_doubled;
//     H[1][1] -= min_eig_doubled;
//   }
// }

// static void get_p(dbl22 const H, dbl2 const g, dbl2 const lam, dbl2 p) {
//   dbl2 tmp;

//   // Set up QP for next iterate
//   triqp2_s qp;
//   dbl22_dbl2_mul(H, lam, tmp);
//   dbl2_sub(g, tmp, qp.b);
//   memcpy(qp.A, H, sizeof(dbl22));

//   // ... solve it
//   triqp2_solve(&qp, EPS);

//   // Compute the projected Newton step
//   dbl2_sub(qp.x, lam, p);
// }

// // TODO: question... would it make more sense to use different
// // vectors for a1 and a2? This choice seems to result in a lot of
// // numerical instability. For now I'm fixing this by replacing sums
// // and dot products involving a1 or a2 with the Neumaier equivalent.
// static dbl const a1[3] = {-1, 1, 0};
// static dbl const a2[3] = {-1, 0, 1};

// static void set_L_and_x_minus_xb(utetra_s *u, dbl3 const b) {
//   dbl33_dbl3_mul(u->X, b, u->xb);
//   dbl3_sub(u->x, u->xb, u->x_minus_xb);
//   u->L = dbl3_norm(u->x_minus_xb);
//   assert(u->L > 0);
// }

// static void set_f_g_and_H_stype_constant(utetra_s *u, dbl3 const b) {
//   dbl3 tmp1;
//   dbl33 tmp2;
//   dbl2 DL, DT;
//   dbl22 D2L, D2T;

//   dbl33_dbl3_mul(u->Xt, u->x_minus_xb, tmp1);
//   dbl3_dbl_div(tmp1, -u->L, tmp1);

//   DL[0] = dbl3_ndot(a1, tmp1);
//   DL[1] = dbl3_ndot(a2, tmp1);
//   assert(dbl2_isfinite(DL));

//   dbl3_outer(tmp1, tmp1, tmp2);
//   dbl33_sub(u->XtX, tmp2, tmp2);
//   dbl33_dbl_div(tmp2, u->L, tmp2);

//   dbl33_dbl3_nmul(tmp2, a1, tmp1);
//   D2L[0][0] = dbl3_ndot(tmp1, a1);
//   D2L[1][0] = D2L[0][1] = dbl3_ndot(tmp1, a2);
//   dbl33_dbl3_nmul(tmp2, a2, tmp1);
//   D2L[1][1] = dbl3_ndot(tmp1, a2);
//   assert(dbl22_isfinite(D2L));

//   DT[0] = bb32_df(&u->T, b, a1);
//   DT[1] = bb32_df(&u->T, b, a2);

//   D2T[0][0] = bb32_d2f(&u->T, b, a1, a1);
//   D2T[1][0] = D2T[0][1] = bb32_d2f(&u->T, b, a1, a2);
//   D2T[1][1] = bb32_d2f(&u->T, b, a2, a2);

//   u->f = u->L + bb32_f(&u->T, b);
//   assert(isfinite(u->f));

//   dbl2_add(DL, DT, u->g);
//   assert(dbl2_isfinite(u->g));

//   dbl22_add(D2L, D2T, u->H);
//   assert(dbl22_isfinite(u->H));
// }

// static void lift_b_from_face_to_cell(uint3 const lf, dbl3 const bf, uint4 lc, dbl4 bc) {
//   dbl4_zero(bc);
//   for (uint i = 0; i < 4; ++i)
//     for (uint j = 0; j < 3; ++j)
//       if (lc[i] == lf[j])
//         bc[i] = bf[j];
// }

// static void get_tb(utetra_s const *u, dbl3 const bf, dbl3 tb, dbl33 hessT, dbl *gradTnorm) {
//   mesh3_s const *mesh = eik3_get_mesh(u->eik);

//   /** Set up BB poly approximating T: */

//   size_t cv[4];
//   mesh3_cv(mesh, u->T_lc, cv);

//   dbl4 bc;
//   lift_b_from_face_to_cell(u->l, bf, cv, bc);

//   jet31t jet[4];
//   for (uint i = 0; i < 4; ++i)
//     jet[i] = eik3_get_jet(u->eik, cv[i]);

//   dbl3 x[4];
//   for (uint i = 0; i < 4; ++i)
//     mesh3_copy_vert(mesh, cv[i], x[i]);

//   bb33 T;
//   bb33_init_from_jets(&T, jet, x);

//   /** Evaluate in Cartesian coordinates: */

//   dbl4 A[3];
//   for (size_t i = 0; i < 3; ++i) {
//     dbl4_zero(A[i]);
//     A[i][i] = 1;
//     A[i][3] = -1;
//   }

//   /* set up dX */
//   dbl33 dX;
//   for (size_t i = 0; i < 3; ++i)
//     dbl3_sub(x[i], x[3], dX[i]);

//   dbl33 dXinv, dXinvT;
//   dbl33_copy(dX, dXinv);
//   dbl33_invert(dXinv);
//   dbl33_transposed(dXinv, dXinvT);

//   /* compute gradient in affine space */
//   dbl3 gradTaff;
//   for (size_t i = 0; i < 3; ++i)
//     gradTaff[i] = bb33_df(&T, bc, A[i]);

//   /* transform gradient back to Cartesian */
//   dbl3 gradT;
//   dbl33_dbl3_mul(dXinv, gradTaff, gradT);

//   if (hessT != NULL) {
//     /* compute Hessian in affine space */
//     dbl33 hessTaff;
//     for (size_t i = 0; i < 3; ++i) {
//       for (size_t j = 0; j < 3; ++j) {
//         dbl4 a[2];
//         dbl4_copy(A[i], a[0]);
//         dbl4_copy(A[j], a[1]);
//         hessTaff[i][j] = bb33_d2f(&T, bc, a);
//       }
//     }

//     /* transform Hessian back to Cartesian */
//     dbl33 tmp;
//     dbl33_mul(dXinv, hessTaff, tmp);
//     dbl33_mul(tmp, dXinvT, hessT);
//   }

//   /** Normalize to get t */

//   if (gradTnorm == NULL) {
//     dbl3_normalized(gradT, tb);
//   } else {
//     *gradTnorm = dbl3_norm(gradT);
//     dbl3_dbl_div(gradT, *gradTnorm, tb);
//   }
// }

// static void get_D2T(utetra_s const *u, dbl3 const bf, dbl33 D2T) {
//   mesh3_s const *mesh = eik3_get_mesh(u->eik);

//   size_t cv[4];
//   mesh3_cv(mesh, u->T_lc, cv);

//   dbl4 bc;
//   lift_b_from_face_to_cell(u->l, bf, cv, bc);

//   jet31t jet[4];
//   for (uint i = 0; i < 4; ++i)
//     jet[i] = eik3_get_jet(u->eik, cv[i]);

//   dbl3 x[4];
//   for (uint i = 0; i < 4; ++i)
//     mesh3_copy_vert(mesh, cv[i], x[i]);

//   bb33 T;
//   bb33_init_from_jets(&T, jet, x);

//   static dbl44 A;
//   for (size_t i = 0; i < 4; ++i) {
//     for (size_t j = 0; j < 3; ++j)
//       A[j][i] = x[i][j];
//     A[3][i] = 1;
//   }
//   dbl44_invert(A);
//   dbl44_transpose(A);

//   dbl4 a[2];
//   for (uint i = 0; i < 3; ++i) {
//     dbl4_copy(A[i], a[0]);
//     for (uint j = 0; j < 3; ++j) {
//       dbl4_copy(A[j], a[1]);
//       D2T[i][j] = bb33_d2f(&T, bc, a);
//     }
//   }
// }

// static void get_phim_and_phipm(utetra_s const *u, dbl3 const tb, dbl3 phim, dbl3 phipm) {
//   for (size_t i = 0; i < 3; ++i)
//     phim[i] = (u->xb[i] + u->x[i])/2 - u->L*(u->topt[i] - tb[i])/8;

//   for (size_t i = 0; i < 3; ++i)
//     phipm[i] = 1.5*u->x_minus_xb[i]/u->L - (tb[i] + u->topt[i])/4;
// }

// static void set_f_g_and_H_stype_func_ptr(utetra_s *u, dbl3 const bf) {
//   /** Compute intermediate quantities: */

//   dbl3 tb;
//   dbl33 hessTb;
//   dbl gradTnorm;
//   get_tb(u, bf, tb, hessTb, &gradTnorm);

//   dbl3 phim, phipm;
//   get_phim_and_phipm(u, tb, phim, phipm);

//   dbl phipmnorm = dbl3_norm(phipm);

//   dbl3 gradL; {
//     dbl3_sub(u->xb, u->x, gradL);
//     dbl3_dbl_div_inplace(gradL, u->L); }

//   dbl3 ellp; /* ellp = -gradL... no real real to compute this... */
//   dbl3_copy(gradL, ellp);
//   dbl3_negate(ellp);

//   dbl33 gradellp; {
//     for (size_t i = 0; i < 3; ++i) {
//       for (size_t j = 0; j < 3; ++j) {
//         dbl delta = i == j ? 1 : 0;
//         gradellp[i][j] = -(delta - ellp[i]*ellp[j])/u->L; } } }

//   dbl33 gradtb; {
//     for (size_t i = 0; i < 3; ++i) {
//       for (size_t j = 0; j < 3; ++j) {
//         dbl delta = i == j ? 1 : 0;
//         gradtb[i][j] = (delta - tb[i]*tb[j])*hessTb[i][j]/gradTnorm; } } }

//   dbl33 gradphim; {
//     dbl33_eye(gradphim);
//     dbl33_dbl_div_inplace(gradphim, 2);
//     for (size_t i = 0; i < 3; ++i)
//       for (size_t j = 0; j < 3; ++j)
//         gradphim[i][j] -= ((u->topt[i] - tb[i])*gradL[j] - u->L*gradtb[i][j])/8; }

//   dbl33 gradphipm; {
//     dbl33_zero(gradphipm);
//     dbl33_saxpy_inplace(3./2, gradellp, gradphipm);
//     dbl33_saxpy_inplace(-1./4, gradtb, gradphipm); }

//   /** Compute s and grads at each point: */

//   dbl sb, sm, shat;
//   dbl3 gradsb, gradsm;

//   sb = u->sfunc->funcs.s(u->xb);
//   sm = u->sfunc->funcs.s(phim);
//   shat = u->sfunc->funcs.s(u->x);

//   u->sfunc->funcs.Ds(u->xb, gradsb);
//   u->sfunc->funcs.Ds(phim, gradsm);

//   /** Quick sanity check: */

//   dbl3 tmp;
//   dbl33_dbl3_mul(hessTb, tb, tmp);

//   /** Compute Q: */

//   dbl Q = u->L*(sb + 4*sm*phipmnorm + shat)/6;

//   /** Compute gradQ: */

//   dbl3 gradQ = {0, 0, 0};

//   // + Q*gradL/L
//   dbl3_saxpy_inplace(Q/u->L, gradL, gradQ);

//   // + (L/6)*gradsb
//   dbl3_saxpy_inplace(u->L/6, gradsb, gradQ);

//   // + (2*L/3)*phipmnorm*dot(gradphipm, gradsm)
//   { dbl3 tmp;
//     dbl33 gradphimT;
//     dbl33_transposed(gradphim, gradphimT);
//     dbl33_dbl3_mul(gradphimT, gradsm, tmp);
//     dbl3_saxpy_inplace(2*u->L*phipmnorm/3, tmp, gradQ); }

//   // + (2*L/3)*s(phim)*dot(gradphipm, phipm)/phipmnorm
//   { dbl3 tmp;
//     dbl33 gradphipmT;
//     dbl33_transposed(gradphipm, gradphipmT);
//     dbl33_dbl3_mul(gradphipmT, phipm, tmp);
//     dbl3_saxpy_inplace(2*u->L*sm/(3*phipmnorm), tmp, gradQ); }

//   /** Compute T and DT: */

//   dbl T = bb32_f(&u->T, bf);
//   dbl2 DT = {bb32_df(&u->T, bf, a1), bb32_df(&u->T, bf, a2)};

//   /** Compute f: */

//   u->f = T + Q;

//   /** Compute Df (a.k.a. "g"): */

//   dbl2 DQ; {
//     dbl3 tmp;
//     dbl33_dbl3_mul(u->Xt, gradQ, tmp);
//     DQ[0] = dbl3_ndot(a1, tmp);
//     DQ[1] = dbl3_ndot(a2, tmp); }

//   dbl2_add(DT, DQ, u->g);

//   /** Compute D2f (a.k.a. "H"): */

//   // dbl22_add(D2T, D2Q, u->H);

//   // TODO: for now, we just set the Hessian to I, reducing this to
//   // projected gradient descent... just want to test this for now
//   dbl22_eye(u->H);
// }

// static void get_s_and_Ds(utetra_s const *u, uint4 const cv, dbl4 const bc, dbl *s, dbl3 Ds) {
//   // TODO: can't remember why I'm requiring this here...
//   assert(u->stype == STYPE_JET31T);

//   mesh3_s const *mesh = eik3_get_mesh(u->eik);

//   jet31t jet[4];
//   for (uint i = 0; i < 4; ++i)
//     jet[i] = u->sfunc->data_jet31t[cv[i]];

//   dbl3 x[4];
//   for (uint i = 0; i < 4; ++i)
//     mesh3_copy_vert(mesh, cv[i], x[i]);

//   bb33 s_bb;
//   bb33_init_from_jets(&s_bb, jet, x);

//   static dbl44 A;
//   for (size_t i = 0; i < 4; ++i) {
//     for (size_t j = 0; j < 3; ++j)
//       A[j][i] = x[i][j];
//     A[3][i] = 1;
//   }
//   dbl44_invert(A);
//   dbl44_transpose(A);

//   *s = bb33_f(&s_bb, bc);

//   for (uint i = 0; i < 3; ++i)
//     Ds[i] = bb33_df(&s_bb, bc, A[i]);
// }

// static void set_f_g_and_H_stype_jet31t(utetra_s *u, dbl3 const bf) {
//   assert(false); // TODO: this is just a rough draft and hasn't been
//                   // tested! probably doesn't work!

//   mesh3_s const *mesh = eik3_get_mesh(u->eik);

//   dbl T = bb32_f(&u->T, bf);

//   dbl3 tb;
//   dbl gradTnorm;
//   get_tb(u, bf, tb, NULL, &gradTnorm);

//   dbl3 phim, phipm;
//   get_phim_and_phipm(u, tb, phim, phipm);

//   dbl phipmnorm = dbl3_norm(phipm);

//   dbl2 Dphipmnorm;
//   { dbl3 ellp;
//     for (size_t i = 0; i < 3; ++i)
//       ellp[i] = u->x_minus_xb[i]/u->L;

//     dbl33 Dx_ellp;
//     for (size_t i = 0; i < 3; ++i)
//       for (size_t j = 0; j < 3; ++j)
//         Dx_ellp[i][j] = -(i == j ? 1 : 0 - ellp[i]*ellp[j])/u->L;

//     dbl33 Dx_tb;
//     { dbl33 D2T;
//       get_D2T(u, bf, D2T);

//       dbl33 tmp;
//       for (size_t i = 0; i < 3; ++i)
//         for (size_t j = 0; j < 3; ++j)
//           Dx_tb[i][j] = (i == j ? 1 : 0 - tb[i]*tb[j])/gradTnorm;

//       dbl33_mul(tmp, D2T, Dx_tb); }

//     dbl33 Dx_phipm;
//     for (size_t i = 0; i < 3; ++i)
//       for (size_t j = 0; j < 3; ++j)
//         Dx_phipm[i][j] = 3*Dx_ellp[i][j]/2 - Dx_tb[i][j]/4;

//     dbl3 Dx_phipmnorm;
//     dbl33_dbl3_mul(Dx_phipm, phipm, Dx_phipmnorm);
//     dbl3_dbl_div_inplace(Dx_phipmnorm, phipmnorm);

//     dbl3 tmp;
//     dbl33_dbl3_mul(u->Xt, Dx_phipmnorm, tmp);

//     Dphipmnorm[0] = dbl3_ndot(a1, tmp);
//     Dphipmnorm[1] = dbl3_ndot(a2, tmp); }

//   uint4 s_cv;
//   mesh3_cv(mesh, u->s_lc, s_cv);

//   dbl4 bc;
//   lift_b_from_face_to_cell(u->l, bf, s_cv, bc);

//   dbl4 bcm;
//   { tetra3 s_tetra = mesh3_get_tetra(mesh, u->s_lc);
//     tetra3_get_bary_coords(&s_tetra, phim, bcm); }

//   dbl s, sm;
//   dbl3 Ds, Dsm;
//   get_s_and_Ds(u, s_cv, bc, &s, Ds);
//   get_s_and_Ds(u, s_cv, bcm, &sm, Dsm);

//   dbl shat = u->sfunc->data_jet31t[u->lhat].f;

//   dbl Q = u->L*(s + 4*sm*phipmnorm + shat)/6;

//   dbl3 tmp1;
//   dbl33_dbl3_mul(u->Xt, u->x_minus_xb, tmp1);
//   dbl3_dbl_div(tmp1, -u->L, tmp1);

//   dbl2 DL;
//   DL[0] = dbl3_ndot(a1, tmp1);
//   DL[1] = dbl3_ndot(a2, tmp1);
//   assert(dbl2_isfinite(DL));

//   dbl2 DT = {bb32_df(&u->T, bf, a1), bb32_df(&u->T, bf, a2)};

//   dbl L_over_6 = u->L/6;

//   dbl2 DQ = {0, 0};
//   dbl2_saxpy_inplace(Q/u->L, DL, DQ);
//   dbl2_saxpy_inplace(L_over_6, Ds, DQ);
//   dbl2_saxpy_inplace(L_over_6*(4*phipmnorm), Dsm, DQ);
//   dbl2_saxpy_inplace(L_over_6*(4*sm), Dphipmnorm, DQ);

//   // dbl22 D2Q;

//   // dbl22 D2T;
//   // D2T[0][0] = bb32_d2f(&u->T, b, a1, a1);
//   // D2T[1][0] = D2T[0][1] = bb32_d2f(&u->T, b, a1, a2);
//   // D2T[1][1] = bb32_d2f(&u->T, b, a2, a2);

//   u->f = T + Q;
//   dbl2_add(DT, DQ, u->g);
//   // dbl22_add(D2T, D2Q, u->H);

//   // TODO: for now, we just set the Hessian to I, reducing this to
//   // projected gradient descent... just want to test this for now
//   dbl22_eye(u->H);
// }

// static void set_f_g_and_H(utetra_s *u, dbl3 const bf) {
//   switch(u->stype) {
//   case STYPE_CONSTANT:
//     set_f_g_and_H_stype_constant(u, bf);
//     break;
//   case STYPE_FUNC_PTR:
//     set_f_g_and_H_stype_func_ptr(u, bf);
//     break;
//   case STYPE_JET31T:
//     set_f_g_and_H_stype_jet31t(u, bf);
//     break;
//   default:
//     assert(false);
//   }
// }

// static void get_dGdt_stype_func_ptr(utetra_s const *u, dbl3 const bf, dbl3 dGdt) {
//   dbl3 tb;
//   get_tb(u, bf, tb, NULL, NULL);

//   dbl3 phim, phipm;
//   get_phim_and_phipm(u, tb, phim, phipm);

//   dbl sm = u->sfunc->funcs.s(phim);
//   dbl3 Dsm; u->sfunc->funcs.Ds(phim, Dsm);

//   dbl phipmnorm = dbl3_norm(phipm);

//   dbl3_zero(dGdt);
//   dbl3_saxpy_inplace(u->L*phipmnorm/2, Dsm, dGdt);
//   dbl3_saxpy_inplace(sm/phipmnorm, phipm, dGdt);
// }

// static void get_dGdt_stype_jet31t(utetra_s const *u, dbl3 const bf, dbl3 dGdt) {
//   mesh3_s const *mesh = eik3_get_mesh(u->eik);

//   dbl3 tb;
//   get_tb(u, bf, tb, NULL, NULL);

//   dbl3 phim, phipm;
//   get_phim_and_phipm(u, tb, phim, phipm);

//   dbl4 bcm;
//   { tetra3 s_tetra = mesh3_get_tetra(mesh, u->s_lc);
//     tetra3_get_bary_coords(&s_tetra, phim, bcm); }

//   dbl sm;
//   dbl3 Dsm;
//   { uint4 s_cv;
//     mesh3_cv(mesh, u->s_lc, s_cv);
//     get_s_and_Ds(u, s_cv, bcm, &sm, Dsm); }

//   dbl phipmnorm = dbl3_norm(phipm);

//   dbl3_zero(dGdt);
//   dbl3_saxpy_inplace(u->L*phipmnorm/2, Dsm, dGdt);
//   dbl3_saxpy_inplace(sm/phipmnorm, phipm, dGdt);
// }

// static void get_dGdt(utetra_s const *u, dbl3 const bf, dbl3 dGdt) {
//   assert(u->stype != STYPE_CONSTANT);
//   switch (u->stype) {
//   case STYPE_FUNC_PTR:
//     get_dGdt_stype_func_ptr(u, bf, dGdt);
//     break;
//   case STYPE_JET31T:
//     get_dGdt_stype_jet31t(u, bf, dGdt);
//     break;
//   default:
//     assert(false);
//   }
// }

// static void set_topt(utetra_s *u, dbl3 const b) {
//   assert(u->stype == STYPE_CONSTANT);

//   /* Initialize topt to tangent vector of straight line connecting xb
//    * and xhat the first time this is called */
//   if (!dbl3_isfinite(u->topt))
//     dbl3_normalized(u->x_minus_xb, u->topt);

//   /* Just use a fixed point iteration here to find optimal topt */
//   dbl angle;
//   do {
//     dbl3 dGdt;
//     get_dGdt(u, b, dGdt);
//     dbl3_normalize(dGdt);
//     angle = fabs(dbl3_angle(u->topt, dGdt));
//     dbl3_copy(dGdt, u->topt);
//   } while (angle > u->tol);

//   /* TODO: Not sure what a good tolerance here is, but they hopefully
//    * shouldn't deviate too much! */
//   assert(dbl3_angle(u->x_minus_xb, u->topt) < JMM_PI/2);
// }

// static void set_lambda(utetra_s *u, dbl2 const lam) {
//   assert(u->stype == STYPE_CONSTANT);

//   u->lam[0] = lam[0];
//   u->lam[1] = lam[1];

//   dbl3 b = {1 - lam[0] - lam[1], lam[0], lam[1]};
//   assert(dbl3_valid_bary_coord(b));
//   dbl33_dbl3_mul(u->X, b, u->xb);

//   set_L_and_x_minus_xb(u, b);

//   /* Set topt before updating f, g, and H */
//   set_topt(u, b);

//   /* Set the cost function value, its gradient, and its Hessian */
//   set_f_g_and_H(u, b);

//   /* Now, compute Newton step solving the minimization problem:
//    *
//    *     minimize  y’*H*y/2 + [g - H*x]’*y + [x’*H*x/2 - g’*x + f(x)]
//    *   subject to  x >= 0
//    *               sum(x) <= 1
//    *
//    * perturbing the Hessian below should ensure a descent
//    * direction. (It would be interesting to see if we can remove the
//    * perturbation entirely.) */

//   // Possibly perturb the Hessian to make it positive definite
//   perturb_hessian_if_necessary(u->H);

//   // Compute the projected Newton step from the current iterate to the
//   // next iterate
//   get_p(u->H, u->g, lam, u->p);
// }

// static void step(utetra_s *u) {
//   assert(u->stype == STYPE_CONSTANT);

//   dbl const atol = 1e-15, c1 = 1e-4;

//   dbl lam1[2], f, c1_times_g_dot_p, beta;

//   // Get values for current iterate
//   dbl lam[2] = {u->lam[0], u->lam[1]};
//   dbl p[2] = {u->p[0], u->p[1]};
//   f = u->f;

//   // Do backtracking line search
//   beta = 1;
//   c1_times_g_dot_p = c1*dbl2_dot(p, u->g);
//   dbl2_saxpy(beta, p, lam, lam1);
//   set_lambda(u, lam1);
//   while (u->f > f + beta*c1_times_g_dot_p + atol) {
//     beta *= 0.9;
//     assert(beta > 1e-16);
//     dbl2_saxpy(beta, p, lam, lam1);
//     set_lambda(u, lam1);
//   }

//   ++u->niter;
// }

dbl eval_poly(dbl const *a, dbl const *lam) {
  dbl x = lam[0], y = lam[1];
  return a[0] + a[1]*x + a[2]*y + a[3]*x*x + a[4]*x*y + a[5]*y*y;
}

void contract(dbl2 const lam_center, dbl factor, dbl2 lam) {
  for (size_t i = 0; i < 2; ++i)
    lam[i] = (lam[i] - lam_center[i])/factor + lam_center[i];
}

/**
 * Do a tetrahedron update starting at `lam`, writing the result to
 * `jet`. If `lam` is `NULL`, then the first iterate will be selected
 * automatically.
 */
void utetra_solve(utetra_s *u, dbl const *lam) {
  static int CALL_NUMBER = 0;


  // DEBUGGING


  // FILE *fp = fopen("f.bin", "wb");
  // for (size_t i = 0; i <= 100; ++i) {
  //   for (size_t j = 0; j <= 100; ++j) {
  //     dbl2 lam = {i/100., j/100.};
  //     dbl f = NAN;
  //     if (lam[0] + lam[1] <= 1) {
  //       set_lambda(u, lam);
  //       f = u->f;
  //     }
  //     fwrite(&f, sizeof(dbl), 1, fp);
  //   }
  // }
  // fclose(fp);


  // if (u->stype == STYPE_FUNC_PTR) {

    dbl2 lam_prev = {NAN, NAN}, lam_opt = {NAN, NAN};

    dbl2 lam_node[6] = {
      {0, 0},   {0.5, 0},   {1, 0},
      {0, 0.5}, {0.5, 0.5},
      {0, 1}
    };

    dbl beta = 10.0;
    dbl factor = (beta + 1)/beta;
    dbl prev_error = NAN;

    size_t num_iter = 0;

    // printf("DOING S != 1 UTETRA #%d\n", CALL_NUMBER);

    while (true) {

      // printf("* it = %lu\n", num_iter);

      dbl f[6] = {NAN, NAN, NAN, NAN, NAN, NAN};
      for (size_t i = 0; i < 6; ++i) {
        dbl const *lam_ = lam_node[i];
        dbl3 x_node;
        dbl3 b = {1 - lam_[0] - lam_[1], lam_[0], lam_[1]};
        dbl33_dbl3_mul(u->X, b, x_node);

        dbl T = bb32_f(&u->T, b);

        uline_s *uline;
        uline_alloc(&uline);
        uline_init_from_points(uline, u->eik, u->x, x_node, u->tol, T);
        uline_solve(uline);

        f[i] = uline_get_value(uline);

        uline_dealloc(&uline);
      }

      dbl const invV[6][6] = {
        { 1,  0,  0,  0,  0,  0},
        {-3,  4, -1,  0,  0,  0},
        {-3,  0,  0,  4,  0, -1},
        { 2, -4,  2,  0,  0,  0},
        { 4, -4,  0, -4,  4,  0},
        { 2,  0,  0, -4,  0,  2}
      };

      dbl a[6];
      for (size_t i = 0; i < 6; ++i) {
        a[i] = 0;
        for (size_t j = 0; j < 6; ++j) {
          a[i] += invV[i][j]*f[j];
        }
      }

      /* Check that everything is correct at the nodal values... */
      dbl2 const lam_node_orig[6] = {
        {0, 0},   {0.5, 0},   {1, 0},
        {0, 0.5}, {0.5, 0.5},
        {0, 1}
      };
      for (size_t i = 0; i < 6; ++i)
        assert(fabs(eval_poly(a, lam_node_orig[i]) - f[i]) < 1e-12);

      // /* Write to disk... */
      // FILE *fp = fopen("f_poly.bin", "w");
      // for (size_t i = 0; i <= 100; ++i) {
      //   for (size_t j = 0; j <= 100; ++j) {
      //     dbl2 lam = {i/100., j/100.};
      //     dbl f = NAN;
      //     if (lam[0] + lam[1] <= 1) {
      //       f = eval_poly(a, lam);
      //     }
      //     fwrite(&f, sizeof(dbl), 1, fp);
      //   }
      // }
      // fclose(fp);

      triqp2_s qp = {
        .b = {a[1], a[2]},
        .A = {{2*a[3], a[4]}, {a[4], 2*a[5]}},
        .x = {NAN, NAN}
      };

      triqp2_solve(&qp, pow(u->tol, 2));

      // printf("  - lam_node = np.array([[%g, %g]", lam_node[0][0], lam_node[0][1]);
      // for (size_t i = 1; i < 6; ++i)
      //   printf(", [%g, %g]", lam_node[i][0], lam_node[i][1]);
      // printf("])\n");
      // printf("  - lam = np.array([%g, %g]), error = %g\n", qp.x[0], qp.x[1], dbl2_dist(qp.x, lam_prev));

      lam = &qp.x[0];
      dbl error = dbl2_dist(lam, lam_prev);
      if (error <= u->tol) {
        dbl2_copy(lam, lam_opt);
        break;
      } else {
        dbl2_copy(lam, lam_prev);
      }

      if (error > 2*prev_error) {
        beta += 1;
        factor = (beta + 1)/beta;
        // printf("  ! reduced factor to %g\n", factor);
      }

      for (size_t i = 0; i < 6; ++i) {
        contract(lam_prev, factor, lam_node[i]);
        assert(lam_node[i][0] >= -EPS);
        assert(lam_node[i][1] >= -EPS);
        assert(lam_node[i][0] + lam_node[i][1] <= 1 + EPS);
      }

      prev_error = error;
      ++num_iter;

      if (num_iter == MAX_NITER) {
        log_warn("utetra_solve: reached max no. iters");
        dbl2_copy(lam, lam_opt);
        break;
      }
    }

    /* make sure to set u->lam now */
    dbl2_copy(lam_opt, u->lam);

    /** Set f and topt now */

    assert(isinf(u->f));
    assert(dbl3_all_nan(u->topt));

    dbl3 bopt = {1 - lam_opt[0] - lam_opt[1], lam_opt[0], lam_opt[1]};
    dbl Topt = bb32_f(&u->T, bopt);

    dbl3 xopt;
    dbl33_dbl3_mul(u->X, bopt, xopt);

    uline_s *uline;
    uline_alloc(&uline);
    uline_init_from_points(uline, u->eik, u->x, xopt, u->tol, Topt);
    uline_solve(uline);

    u->f = uline_get_value(uline);
    uline_get_topt(uline, u->topt);

    uline_dealloc(&uline);
  // }

  // /////

  // else if (u->stype == STYPE_CONSTANT) {
  //   u->niter = 0;

  //   if (lam == NULL)
  //     set_lambda(u, (dbl[2]) {1./3, 1./3});
  //   else
  //     set_lambda(u, lam);

  //   int const max_niter = 100;
  //   for (int _ = 0; _ < max_niter; ++_) {
  //     if (dbl2_norm(u->p) <= u->tol)
  //       break;
  //     step(u);
  //   }
  //   if (dbl2_norm(u->p) > u->tol)
  //     log_warn("utetra_solve: exceeded max no. of iters.");
  // }

  // else { assert(false); } // TODO: stype not implemented

  ++CALL_NUMBER;
}

static void get_b(utetra_s const *u, dbl b[3]) {
  assert(!isnan(u->lam[0]) && !isnan(u->lam[1]));
  b[0] = 1 - u->lam[0] - u->lam[1];
  b[1] = u->lam[0];
  b[2] = u->lam[1];
}

dbl utetra_get_value(utetra_s const *u) {
  return u->f;
}

// void get_that_stype_constant(utetra_s const *u, dbl that[3]) {
//   dbl3_normalized(u->x_minus_xb, that);
// }

// void get_that_stype_func_ptr(utetra_s const *u, dbl that[3]) {
//   dbl shat = u->sfunc->funcs.s((dbl *)u->x);
//   dbl3_dbl_mul(u->topt, shat, that);
// }

// static void get_that(utetra_s const *u, dbl that[3]) {
//   switch(u->stype) {
//   case STYPE_CONSTANT:
//     get_that_stype_constant(u, that);
//     break;
//   case STYPE_FUNC_PTR:
//     get_that_stype_func_ptr(u, that);
//     break;
//   default:
//     assert(false); // TODO: implement remaining stypes
//   }
// }

void utetra_get_jet31t(utetra_s const *u, jet31t *jet) {
  jet->f = u->f;
  dbl3_copy(u->topt, jet->Df);
}

/**
 * Compute the Lagrange multipliers for the constraint optimization
 * problem corresponding to this type of update
 */
static void get_lag_mults(utetra_s const *u, dbl alpha[3]) {
  dbl const atol = 5e-15;
  dbl b[3] = {1 - dbl2_sum(u->lam), u->lam[0], u->lam[1]};
  alpha[0] = alpha[1] = alpha[2] = 0;
  // TODO: optimize this
  if (fabs(b[0] - 1) < atol) {
    alpha[0] = 0;
    alpha[1] = -u->g[0];
    alpha[2] = -u->g[1];
  } else if (fabs(b[1] - 1) < atol) {
    alpha[0] = u->g[0];
    alpha[1] = 0;
    alpha[2] = u->g[0] - u->g[1];
  } else if (fabs(b[2] - 1) < atol) {
    alpha[0] = u->g[0];
    alpha[1] = u->g[0] - u->g[1];
    alpha[2] = 0;
  } else if (fabs(b[0]) < atol) { // b[1] != 0 && b[2] != 0
    alpha[0] = dbl2_sum(u->g)/2;
    alpha[1] = 0;
    alpha[2] = 0;
  } else if (fabs(b[1]) < atol) { // b[0] != 0 && b[2] != 0
    alpha[0] = 0;
    alpha[1] = -u->g[0];
    alpha[2] = 0;
  } else if (fabs(b[2]) < atol) { // b[0] != 0 && b[1] != 0
    alpha[0] = 0;
    alpha[1] = 0;
    alpha[2] = -u->g[1];
  } else {
    assert(dbl3_valid_bary_coord(b));
  }
}

bool utetra_has_interior_point_solution(utetra_s const *u) {
  /* check whether minimizer is an interior point minimizer */
  dbl alpha[3];
  get_lag_mults(u, alpha);
  return dbl3_maxnorm(alpha) <= 1e-15;
}

bool utetra_is_backwards(utetra_s const *utetra, eik3_s const *eik) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  dbl3 x0, dx[2];
  mesh3_copy_vert(mesh, utetra->l[0], x0);
  dbl3_sub(mesh3_get_vert_ptr(mesh, utetra->l[1]), x0, dx[0]);
  dbl3_sub(mesh3_get_vert_ptr(mesh, utetra->l[2]), x0, dx[1]);

  dbl3 t;
  dbl3_cross(dx[0], dx[1], t);

  dbl3 dxhat0;
  dbl3_sub(mesh3_get_vert_ptr(mesh, utetra->lhat), x0, dxhat0);

  if (dbl3_dot(dxhat0, t) < 0)
    dbl3_negate(t);

  jet31t const *jet = eik3_get_jet_ptr(eik);

  for (size_t i = 0; i < 3; ++i)
    if (dbl3_dot(t, jet[utetra->l[i]].Df) <= 0)
      return true;

  return false;
}

#if JMM_DEBUG
static bool update_inds_are_set(utetra_s const *utetra) {
  return !(
    utetra->l[0] == (size_t)NO_INDEX ||
    utetra->l[1] == (size_t)NO_INDEX ||
    utetra->l[2] == (size_t)NO_INDEX);
}
#endif

#if JMM_DEBUG
static bool all_inds_are_set(utetra_s const *utetra) {
  return utetra->lhat != (size_t)NO_INDEX && update_inds_are_set(utetra);
}
#endif

bool utetra_ray_is_occluded(utetra_s const *utetra, eik3_s const *eik) {
#if JMM_DEBUG
  assert(all_inds_are_set(utetra));
#endif
  mesh3_s const *mesh = eik3_get_mesh(eik);
  par3_s par = utetra_get_parent(utetra);
  return mesh3_local_ray_is_occluded(mesh, utetra->lhat, &par);
}

int utetra_get_num_interior_coefs(utetra_s const *utetra) {
  dbl const atol = 1e-14;
  dbl b[3];
  get_b(utetra, b);
  return (atol < b[0] && b[0] < 1 - atol)
       + (atol < b[1] && b[1] < 1 - atol)
       + (atol < b[2] && b[2] < 1 - atol);
}

size_t utetra_get_l(utetra_s const *utetra) {
  assert(utetra->lhat != (size_t)NO_INDEX);
  return utetra->lhat;
}

size_t utetra_get_active_inds(utetra_s const *utetra, uint3 la) {
#if JMM_DEBUG
  assert(update_inds_are_set(utetra));
#endif

  la[0] = la[1] = la[2] = (size_t)NO_INDEX;

  dbl3 b;
  get_b(utetra, b);

  size_t na = 0;
  for (size_t i = 0; i < 3; ++i)
    if (b[i] > EPS)
      la[na++] = utetra->l[i];
  SORT_UINT3(la);

  return na;
}

par3_s utetra_get_parent(utetra_s const *utetra) {
  par3_s par; get_b(utetra, par.b);
  assert(dbl3_valid_bary_coord(par.b));
  memcpy(par.l, utetra->l, sizeof(size_t[3]));
  return par;
}

void utetra_get_other_inds(utetra_s const *utetra, size_t li, size_t l[2]) {
  assert(utetra->l[0] == li || utetra->l[1] == li || utetra->l[2]);

  for (size_t i = 0, j = 0; i < 3; ++i)
    if (utetra->l[i] != li)
      l[j++] = utetra->l[i];
}

bool utetra_index_is_active(utetra_s const *utetra, size_t l) {
  assert(l == utetra->l[0] || l == utetra->l[1] || l == utetra->l[2]);

  dbl lam0 = 1 - utetra->lam[0] - utetra->lam[1];

  if (l == utetra->l[0])
    return lam0 != 0;
  else if (l == utetra->l[1])
    return utetra->lam[0] != 0;
  else if (l == utetra->l[2])
    return utetra->lam[1] != 0;
  else
    assert(false);
}

bool utetra_has_inds(utetra_s const *u, size_t lhat, uint3 const l) {
  return u->lhat == lhat && uint3_equal(u->l, l);
}

static void get_other_inds(array_s const *utetras, size_t la, size_t i, uint2 l_other) {
  utetra_s const *utetra_other;
  array_get(utetras, i, &utetra_other);
  utetra_get_other_inds(utetra_other, la, l_other);
}

static bool
utetras_bracket_ray_1(utetra_s const *utetra, size_t la, array_s const *utetras) {
  size_t num_utetra_other = array_size(utetras);

  size_t i_current = 0;

  uint2 l_other = {NO_INDEX, NO_INDEX};;
  get_other_inds(utetras, la, i_current, l_other);

  size_t l_start = l_other[0];
  size_t l_current = l_other[1];

  /* Next, we walk around the ring of utetra until we've visited all
   * of them... */
  size_t num_visited = 1;
  while (num_visited++ < num_utetra_other) {
    /* Search through the utetra for the one we should step to next */
    for (size_t i = 0; i < num_utetra_other; ++i) {
      if (i == i_current)
        continue;

      get_other_inds(utetras, la, i, l_other);
      if (l_current != l_other[0] && l_current != l_other[1])
        continue;

      /* Update l_current and i_current and break */
      l_current = l_current == l_other[0] ? l_other[1] : l_other[0];
      i_current = i;
      break;

      /* TODO: actually check whether these vertices form a ring
       * around the active index after projecting onto the normal
       * plane */
      (void)utetra;
    }
  }

  /* We're done walking... check if we're back where we started */
  return l_start == l_current;
}

static void get_unit_normal_for_verts_by_index(mesh3_s const *mesh,
                                               uint3 const l, dbl3 n) {
  dbl3 x0;
  mesh3_copy_vert(mesh, l[0], x0);

  dbl3 dx[2];
  for (size_t i = 0; i < 2; ++i)
    dbl3_sub(mesh3_get_vert_ptr(mesh, l[i + 1]), x0, dx[i]);

  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
}

static bool
utetras_bracket_ray_2(utetra_s const *utetra, array_s const *utetras) {
  mesh3_s const *mesh = eik3_get_mesh(utetra->eik);

  size_t l = utetra_get_l(utetra);

  par3_s par = utetra_get_parent(utetra);
  uint3 la, li;
  size_t na = par3_get_active_and_inactive_inds(&par, la, li);
  assert(na == 2);
  SORT_UINT3(la);
  SORT_UINT3(li);

  dbl3 xa, xi;
  mesh3_copy_vert(mesh, la[0], xa);
  mesh3_copy_vert(mesh, li[0], xi);

  for (size_t i = 0; i < array_size(utetras); ++i) {
    utetra_s const *utetra_other;
    array_get(utetras, i, &utetra_other);

    assert(l == utetra_get_l(utetra_other));

    par3_s par_ = utetra_get_parent(utetra_other);

    uint3 la_, li_;
    size_t na_ = par3_get_active_and_inactive_inds(&par_, la_, li_);
    SORT_UINT3(la_);
    SORT_UINT3(li_);

    /* Check if both `utetra` have the same active indices */
    if (na_ != 2 && la_[0] != la[0] && la_[1] != la[1])
      continue;

    /* Get the unit normal for the plane determined by the update
     * index and the active edge */
    dbl3 n;
    get_unit_normal_for_verts_by_index(mesh, (uint3) {l, la[0], la[1]}, n);

    /* Get the inactive vertex for the current other utetra */
    dbl3 xi_;
    mesh3_copy_vert(mesh, li_[0], xi_);

    /* Check whether the two inactive vertices lie on either side of
     * the plane */
    dbl n_dot_xa = dbl3_dot(n, xa);
    dbl dp = dbl3_dot(n, xi) - n_dot_xa;
    dbl dp_ = dbl3_dot(n, xi_) - n_dot_xa;
    if (sgn(dp) != sgn(dp_))
      return true;
  }

  return false;
}

bool utetra_is_bracketed_by_utetras(utetra_s const *utetra, array_s const *utetras) {
  if (array_is_empty(utetras))
    return false;

  par3_s par = utetra_get_parent(utetra);
  uint3 la;
  size_t na = par3_get_active_inds(&par, la);
  assert(na == 1 || na == 2);

  if (na == 1)
    return utetras_bracket_ray_1(utetra, la[0], utetras);
  else if (na == 2)
    return utetras_bracket_ray_2(utetra, utetras);
  else
    assert(false);
}

static dbl get_edge_lam(utetra_s const *u, uint2 const le) {
  assert(uint3_contains_uint2(u->l, le));

  /* Get the barycentric coordinates of the optimum */
  dbl3 b; get_b(u, b);

#if JMM_DEBUG
  { uint3 l_diff;
    size_t n = uint3_diff_uint2(u->l, le, l_diff);
    assert(n == 1);
    assert(fabs(b[uint3_find(u->l, l_diff[0])]) < EPS); }
#endif

  /* Get the index of `le[0]` in `u->l` */
  size_t i0 = uint3_find(u->l, le[0]);

  /* Compute the convex coefficent `lam = 1 - b[i0]` such that `x_opt`
   * is `(1 - lam)*x[le[0]] + lam*x[le[1]]`. */
  dbl lam = 1 - b[i0];

  return lam;
}

size_t get_common_inds(uint3 const la1, uint3 const la2, uint3 le) {
  assert(uint3_is_sorted(la1));
  assert(uint3_is_sorted(la2));

  le[0] = le[1] = le[2] = (size_t)NO_INDEX;

  size_t ne = 0;
  while (la1[ne] == la2[ne] && la1[ne] != (size_t)NO_INDEX) {
    le[ne] = la1[ne];
    ++ne;
  }
  for (size_t i = ne; i < 3; ++i) {
    if ((la1[ne] == (size_t)NO_INDEX) ^ (la2[ne] == (size_t)NO_INDEX)) {
      le[ne] = la1[ne] == (size_t)NO_INDEX ? la2[ne] : la1[ne];
      ++ne;
    }
  }

  return ne;
}

bool utetras_have_same_minimizer(utetra_s const *u1, utetra_s const *u2) {
  assert(u1->eik == u2->eik);

  uint3 la1, la2;
  size_t na1 = utetra_get_active_inds(u1, la1);
  size_t na2 = utetra_get_active_inds(u2, la2);
  assert(na1 < 3 && na2 < 3);

  if (na1 == 1 && na2 == 1)
    return la1[0] == la2[0];

  if (na1 == 2 && na2 == 2 && !uint3_equal(la1, la2))
    return false;

  uint3 le;
  size_t ne = get_common_inds(la1, la2, le);
  assert(ne <= 2);
  if (ne == 0)
    return false;

  /* The effective tolerance used to check the closeness of the two
   * minimizers. Its value depends on which indices are active for
   * each tetrahedron update. */
  dbl edge_tol = mesh3_get_edge_tol(eik3_get_mesh(u1->eik), le);
  edge_tol = fmax(edge_tol, fmax(u1->tol, u2->tol));

  dbl lam1;
  if (na1 == 2)
    lam1 = get_edge_lam(u1, le);
  else
    lam1 = la1[0] == le[0] ? 0 : 1;

  dbl lam2;
  if (na2 == 2)
    lam2 = get_edge_lam(u2, le);
  else
    lam2 = la2[0] == le[0] ? 0 : 1;

  bool same = fabs(lam1 - lam2) <= edge_tol;

  return same;
}

bool utetras_have_same_inds(utetra_s const *u1, utetra_s const *u2) {
  if (u1->lhat != u2->lhat)
    return false;

  size_t l1[3];
  memcpy(l1, u1->l, sizeof(size_t[3]));
  SORT3(l1[0], l1[1], l1[2]);

  size_t l2[3];
  memcpy(l2, u2->l, sizeof(size_t[3]));
  SORT3(l2[0], l2[1], l1[2]);

  return l1[0] == l2[0] && l1[1] == l2[1] && l1[2] == l2[2];
}

#if JMM_TEST
void utetra_step(utetra_s *u) {
  step(u);
}

void utetra_get_lambda(utetra_s const *u, dbl lam[2]) {
  lam[0] = u->lam[0];
  lam[1] = u->lam[1];
}

void utetra_set_lambda(utetra_s *u, dbl const lam[2]) {
  set_lambda(u, lam);
}

size_t utetra_get_num_iter(utetra_s const *u) {
  return u->niter;
}
#endif
