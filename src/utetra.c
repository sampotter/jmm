#include "utetra.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "array.h"
#include "bb.h"
#include "eik3.h"
#include "macros.h"
#include "mat.h"
#include "mesh3.h"
#include "opt.h"
#include "util.h"

#define MAX_NITER 100

struct utetra {
  dbl lam[2]; // Current iterate

  dbl f;
  dbl g[2];
  dbl H[2][2];

  dbl p[2]; // Newton step

  dbl angles[3];

  dbl x[3]; // x[l]
  dbl X[3][3]; // X = [x[l0] x[l1] x[l2]]
  dbl Xt[3][3]; // X'
  dbl XtX[3][3]; // X'*X

  dbl L, T[3];

  // B-coefs for 9-point triangle interpolation T on base of update
  bb32 bb_T;

  dbl x_minus_xb[3];

  size_t lhat, l[3];

  int niter;
};

void utetra_alloc(utetra_s **utetra) {
  *utetra = malloc(sizeof(utetra_s));
}

void utetra_dealloc(utetra_s **utetra) {
  free(*utetra);
  *utetra = NULL;
}

void utetra_reset(utetra_s *cf) {
  cf->niter = 0;
  cf->f = INFINITY;
}

bool utetra_init_from_eik3(utetra_s *cf, eik3_s const *eik,
                           size_t l, size_t l0, size_t l1, size_t l2) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  jet3 const *jet = eik3_get_jet_ptr(eik);
  return utetra_init_from_ptrs(cf, mesh, jet, l, l0, l1, l2);
}

bool utetra_init_from_ptrs(utetra_s *cf, mesh3_s const *mesh, jet3 const *jet,
                           size_t l, size_t l0, size_t l1, size_t l2) {
  dbl x[3];
  mesh3_copy_vert(mesh, l, x);

  dbl Xt[3][3];
  mesh3_copy_vert(mesh, l0, Xt[0]);
  mesh3_copy_vert(mesh, l1, Xt[1]);
  mesh3_copy_vert(mesh, l2, Xt[2]);

  jet3 jet_[3] = {jet[l0], jet[l1], jet[l2]};

  assert(jet3_is_finite(&jet_[0]));
  assert(jet3_is_finite(&jet_[1]));
  assert(jet3_is_finite(&jet_[2]));

  cf->lhat = l;

  cf->l[0] = l0;
  cf->l[1] = l1;
  cf->l[2] = l2;

  return utetra_init(cf, x, Xt, jet_);
}

bool utetra_init_no_inds(utetra_s *cf, dbl const x[3], dbl const Xt[3][3],
                         jet3 const jet[3]) {
  // Set the indices to `NO_INDEX` here to ensure that we don't
  // accidentally trip over weird things later that assume access to
  // indices that make sense.
  cf->lhat = cf->l[0] = cf->l[1] = cf->l[2] = NO_INDEX;

  return utetra_init(cf, x, Xt, jet);
}

bool utetra_init(utetra_s *cf, dbl const x[3], dbl const Xt[3][3],
                 jet3 const jet[3]) {
  // Initialize lambda with a dummy value so that we can check for bad
  // accesses later
  cf->lam[0] = cf->lam[1] = NAN;

  memcpy(cf->x, x, 3*sizeof(dbl));
  memcpy(cf->Xt, Xt, 3*3*sizeof(dbl));

  dbl33_transposed(cf->Xt, cf->X);
  dbl33_mul(cf->Xt, cf->X, cf->XtX);

  dbl Xt_minus_x[3][3];
  for (int i = 0; i < 3; ++i) {
    dbl3_sub(cf->Xt[i], cf->x, Xt_minus_x[i]);
    dbl3_normalize(Xt_minus_x[i]);
  }
  for (int i = 0; i < 3; ++i) {
    cf->angles[i] = dbl3_dot(Xt_minus_x[i], Xt_minus_x[(i + 1) % 3]);
  }

  dbl DT[3][3];
  for (int i = 0; i < 3; ++i) {
    cf->T[i] = jet[i].f;
    DT[i][0] = jet[i].fx;
    DT[i][1] = jet[i].fy;
    DT[i][2] = jet[i].fz;
  }
  bb32_init_from_3d_data(&cf->bb_T, cf->T, &DT[0], cf->Xt);

  cf->niter = 0;

  // Compute the surface normal for the plane spanned by (x1 - x0, x2
  // - x0), using DT[i] to determine its orientation. Return whether x
  // is on the right side of this plane.
  dbl n[3];
  dbl dx[2][3];
  dbl3_sub(cf->Xt[1], cf->Xt[0], dx[0]);
  dbl3_sub(cf->Xt[2], cf->Xt[0], dx[1]);
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
  int sgn[3] = {
    signum(dbl3_dot(DT[0], n)),
    signum(dbl3_dot(DT[1], n)),
    signum(dbl3_dot(DT[2], n))
  };
  if (sgn[0] == -1)
    dbl3_negate(n);
  dbl dot = -dbl3_dot(Xt_minus_x[0], n) > 0;
  return dot > 0;
}

bool utetra_is_degenerate(utetra_s const *cf) {
  // Check if the point being updated lies in the plane spanned by by
  // x0, x1, and x2. If it does, the update is degenerate.
  dbl const *x[4] = {cf->x, cf->Xt[0], cf->Xt[1], cf->Xt[2]};
  return points_are_coplanar(x);
}

bool utetra_is_causal(utetra_s const *cf) {
  return cf->angles[0] >= 0 && cf->angles[1] >= 0 && cf->angles[2] >= 0;
}

void utetra_step(utetra_s *cf) {
  dbl const atol = 1e-15, c1 = 1e-4;

  dbl lam1[2], f, c1_times_g_dot_p, beta;

  // Get values for current iterate
  dbl lam[2] = {cf->lam[0], cf->lam[1]};
  dbl p[2] = {cf->p[0], cf->p[1]};
  f = cf->f;

  // Do backtracking line search
  beta = 1;
  c1_times_g_dot_p = c1*dbl2_dot(p, cf->g);
  dbl2_saxpy(beta, p, lam, lam1);
  utetra_set_lambda(cf, lam1);
  while (cf->f > f + beta*c1_times_g_dot_p + atol) {
    beta *= 0.9;
    dbl2_saxpy(beta, p, lam, lam1);
    utetra_set_lambda(cf, lam1);
  }

  ++cf->niter;
}

/**
 * Do a tetrahedron update starting at `lam`, writing the result to
 * `jet`. This assumes that `utetra_set_lambda` has already been
 * called, so that `cf` is currently at `lam`.
 */
void utetra_solve(utetra_s *cf) {
  // I think we're actually supposed to use cf->g here instead of
  // cf->p, but maybe cf->p is more appropriate for constraint
  // Newton's method, since it takes into account the constraints
  // while cf->p doesn't.
  int const max_niter = 100;
  dbl const atol = 1e-15, rtol = 1e-15, tol = rtol*dbl2_norm(cf->p) + atol;
  for (int _ = 0; _ < max_niter; ++_) {
    if (dbl2_norm(cf->p) <= tol)
      break;
    utetra_step(cf);
  }
  // printf("solved: cf->niter = %d, |cf->p| = %0.16g\n", cf->niter, dbl2_norm(cf->p));
}

void utetra_get_lambda(utetra_s const *cf, dbl lam[2]) {
  assert(!isnan(cf->lam[0]) && !isnan(cf->lam[1]));
  lam[0] = cf->lam[0];
  lam[1] = cf->lam[1];
}

void utetra_set_lambda(utetra_s *cf, dbl const lam[2]) {
  // TODO: question... would it make more sense to use different
  // vectors for a1 and a2? This choice seems to result in a lot of
  // numerical instability. For now I'm fixing this by replacing sums
  // and dot products involving a1 or a2 with the Neumaier equivalent.
  static dbl a1[3] = {-1, 1, 0};
  static dbl a2[3] = {-1, 0, 1};

  dbl b[3], xb[3], tmp1[3], tmp2[3][3], DL[2], D2L[2][2], DT[2], D2T[2][2];
  dbl tmp3[2];

  cf->lam[0] = lam[0];
  cf->lam[1] = lam[1];

  b[1] = lam[0];
  b[2] = lam[1];
  b[0] = 1 - b[1] - b[2];

  assert(dbl3_valid_bary_coord(b));

  dbl33_dbl3_mul(cf->X, b, xb);
  dbl3_sub(cf->x, xb, cf->x_minus_xb);
  cf->L = dbl3_norm(cf->x_minus_xb);
  assert(cf->L > 0);

  dbl33_dbl3_mul(cf->Xt, cf->x_minus_xb, tmp1);
  dbl3_dbl_div(tmp1, -cf->L, tmp1);

  DL[0] = dbl3_ndot(a1, tmp1);
  DL[1] = dbl3_ndot(a2, tmp1);
  assert(dbl2_isfinite(DL));

  dbl3_outer(tmp1, tmp1, tmp2);
  dbl33_sub(cf->XtX, tmp2, tmp2);
  dbl33_dbl_div(tmp2, cf->L, tmp2);

  dbl33_dbl3_nmul(tmp2, a1, tmp1);
  D2L[0][0] = dbl3_ndot(tmp1, a1);
  D2L[1][0] = D2L[0][1] = dbl3_ndot(tmp1, a2);
  dbl33_dbl3_nmul(tmp2, a2, tmp1);
  D2L[1][1] = dbl3_ndot(tmp1, a2);
  assert(dbl22_isfinite(D2L));

  DT[0] = bb32_df(&cf->bb_T, b, a1);
  DT[1] = bb32_df(&cf->bb_T, b, a2);

  D2T[0][0] = bb32_d2f(&cf->bb_T, b, a1, a1);
  D2T[1][0] = D2T[0][1] = bb32_d2f(&cf->bb_T, b, a1, a2);
  D2T[1][1] = bb32_d2f(&cf->bb_T, b, a2, a2);

  cf->f = cf->L + bb32_f(&cf->bb_T, b);
  assert(isfinite(cf->f));

  dbl2_add(DL, DT, cf->g);
  assert(dbl2_isfinite(cf->g));

  dbl22_add(D2L, D2T, cf->H);
  assert(dbl22_isfinite(cf->H));

  /**
   * Finally, compute Newton step solving the minimization problem:
   *
   *     minimize  y’*H*y/2 + [g - H*x]’*y + [x’*H*x/2 - g’*x + f(x)]
   *   subject to  x >= 0
   *               sum(x) <= 1
   *
   * perturbing the Hessian below should ensure a descent
   * direction. (It would be interesting to see if we can remove the
   * perturbation entirely.)
   */

  // Conditionally perturb the Hessian
  dbl tr = cf->H[0][0] + cf->H[1][1];
  dbl det = cf->H[0][0]*cf->H[1][1] - cf->H[0][1]*cf->H[1][0];
  dbl min_eig_doubled = tr - sqrt(tr*tr - 4*det);
  if (min_eig_doubled < 0) {
    cf->H[0][0] -= min_eig_doubled;
    cf->H[1][1] -= min_eig_doubled;
  }

  // Solve a quadratic program to find the next iterate.
  triqp2_s qp;
  dbl22_dbl2_mul(cf->H, lam, tmp3);
  dbl2_sub(cf->g, tmp3, qp.b);
  memcpy((void *)qp.A, (void *)cf->H, sizeof(dbl)*2*2);
  triqp2_solve(&qp);

  // Compute the projected Newton step from the current iterate and
  // next iterate.
  dbl2_sub(qp.x, lam, cf->p);
}

void utetra_get_bary_coords(utetra_s const *cf, dbl b[3]) {
  assert(!isnan(cf->lam[0]) && !isnan(cf->lam[1]));
  b[0] = 1 - cf->lam[0] - cf->lam[1];
  b[1] = cf->lam[0];
  b[2] = cf->lam[1];
}

dbl utetra_get_value(utetra_s const *cf) {
  return cf->f;
}

void utetra_get_gradient(utetra_s const *cf, dbl g[2]) {
  g[0] = cf->g[0];
  g[1] = cf->g[1];
}

void utetra_get_jet(utetra_s const *cf, jet3 *jet) {
  jet->f = cf->f;
  jet->fx = cf->x_minus_xb[0]/cf->L;
  jet->fy = cf->x_minus_xb[1]/cf->L;
  jet->fz = cf->x_minus_xb[2]/cf->L;
}

/**
 * Compute the Lagrange multipliers for the constraint optimization
 * problem corresponding to this type of update
 */
void utetra_get_lag_mults(utetra_s const *cf, dbl alpha[3]) {
  dbl const atol = 5e-15;
  dbl b[3] = {1 - dbl2_sum(cf->lam), cf->lam[0], cf->lam[1]};
  alpha[0] = alpha[1] = alpha[2] = 0;
  // TODO: optimize this
  if (fabs(b[0] - 1) < atol) {
    alpha[0] = 0;
    alpha[1] = -cf->g[0];
    alpha[2] = -cf->g[1];
  } else if (fabs(b[1] - 1) < atol) {
    alpha[0] = cf->g[0];
    alpha[1] = 0;
    alpha[2] = cf->g[0] - cf->g[1];
  } else if (fabs(b[2] - 1) < atol) {
    alpha[0] = cf->g[0];
    alpha[1] = cf->g[0] - cf->g[1];
    alpha[2] = 0;
  } else if (fabs(b[0]) < atol) { // b[1] != 0 && b[2] != 0
    alpha[0] = dbl2_sum(cf->g)/2;
    alpha[1] = 0;
    alpha[2] = 0;
  } else if (fabs(b[1]) < atol) { // b[0] != 0 && b[2] != 0
    alpha[0] = 0;
    alpha[1] = -cf->g[0];
    alpha[2] = 0;
  } else if (fabs(b[2]) < atol) { // b[0] != 0 && b[1] != 0
    alpha[0] = 0;
    alpha[1] = 0;
    alpha[2] = -cf->g[1];
  } else {
    assert(dbl3_valid_bary_coord(b));
  }
}

int utetra_get_num_iter(utetra_s const *cf) {
  return cf->niter;
}

bool utetra_has_interior_point_solution(utetra_s const *cf) {
  dbl alpha[3];
  utetra_get_lag_mults(cf, alpha);
  return dbl3_maxnorm(alpha) <= 1e-15;
}

int utetra_cmp(utetra_s const **h1, utetra_s const **h2) {
  utetra_s const *u1 = *h1;
  utetra_s const *u2 = *h2;
  if (u1 == NULL && u2 == NULL) {
    return 0;
  } else if (u2 == NULL) {
    return -1;
  } else if (u1 == NULL) {
    return 1;
  } else {
    dbl T1 = u1->f, T2 = u2->f;
    if (T1 < T2) {
      return -1;
    } else if (T1 > T2) {
      return 1;
    } else {
      return 0;
    }
  }
}

bool utetra_adj_are_optimal(utetra_s const *u1, utetra_s const *u2) {
  dbl const atol = 1e-15;
  dbl lam1[2], lam2[2];
  utetra_get_lambda(u1, lam1);
  utetra_get_lambda(u2, lam2);
  return fabs(lam1[0] - lam2[0]) <= atol
    && fabs(lam1[1]) <= atol
    && fabs(lam2[1]) <= atol
    && fabs(utetra_get_value(u1) - utetra_get_value(u2)) <= atol;
}

ray3 utetra_get_ray(utetra_s const *utetra) {
  ray3 ray;
  dbl b[3];
  utetra_get_bary_coords(utetra, b);
  dbl33_dbl3_mul(utetra->X, b, ray.org);
  dbl3_normalized(utetra->x_minus_xb, ray.dir);
  return ray;
}

void utetra_get_point_on_ray(utetra_s const *utetra, dbl t, dbl xt[3]) {
  // TODO: optimize this by using utetra->x instead of computing xb
  dbl b[3], xb[3], L;
  utetra_get_bary_coords(utetra, b);
  dbl33_dbl3_mul(utetra->X, b, xb);
  L = dbl3_norm(utetra->x_minus_xb);
  dbl3_saxpy(t/L, utetra->x_minus_xb, xb, xt);
}

int utetra_get_interior_coefs_mask(utetra_s const *utetra, bool I[3]) {
  dbl const atol = 1e-14;
  I[0] = utetra->lam[0] + utetra->lam[1] < 1 - atol;
  I[1] = utetra->lam[0] > atol;
  I[2] = utetra->lam[1] > atol;
  return I[0] + I[1] + I[2];
}

bool utetra_update_inds_are_set(utetra_s const *utetra) {
  return !(
    utetra->l[0] == (size_t)NO_INDEX ||
    utetra->l[1] == (size_t)NO_INDEX ||
    utetra->l[2] == (size_t)NO_INDEX);
}

bool utetra_inds_are_set(utetra_s const *utetra) {
  return utetra->lhat != (size_t)NO_INDEX && utetra_update_inds_are_set(utetra);
}

bool utetra_update_ray_is_physical(utetra_s const *utetra, eik3_s const *eik) {
  assert(utetra_inds_are_set(utetra));

  size_t const *l = utetra->l;

  mesh3_s const *mesh = eik3_get_mesh(eik);

  // In the following section, we want to quickly look at the boundary
  // near the start of the ray and see if we can rule out the update
  // based on this information alone.
  //
  // TODO: some of the following checks could be fine if the tangent
  // vector at the ray origin lies in the plane spanned by the update
  // base, but we don't have to worry about that until we deal with
  // nonconstant speeds

  int num_int = utetra_get_num_interior_coefs(utetra);

  size_t l_int[3] = {NO_INDEX, NO_INDEX, NO_INDEX};
  utetra_get_interior_coefs(utetra, l_int);

  if (num_int == 2 && mesh3_is_nondiff_boundary_edge(mesh, l_int))
    return false;

  if (num_int == 3 && mesh3_bdf(mesh, l))
    return false;

  // TODO: the following section where we check to see if the stuff
  // below gives "an interior ray" can be wrapped up and reused for
  // both this and the corresponding section in utri.c...

  ray3 ray = utetra_get_ray(utetra);

  /**
   * First, check if the start of the ray is "in free space". To do
   * this, we just check if we can a cell containing a point just
   * ahead of and just behind the origin of the ray.
   */

  // Get points just before and just after the start of the ray. We
  // perturb forward and backward by one half of the minimum triangle
  // altitude (taken of the entire mesh). This is to ensure that the
  // perturbations are small enough to stay inside neighboring
  // tetrahedra, but also large enough to take into consideration the
  // characteristic length scale of the mesh.
  dbl xm[3], xp[3], t = mesh3_get_min_tetra_alt(mesh)/2;
  ray3_get_point(&ray, -t, xm);
  ray3_get_point(&ray, t, xp);

  // Find the number and location of interior coefficients.
  bool I[3];
  utetra_get_interior_coefs_mask(utetra, I);

  bool xm_in_cell = false, xp_in_cell = false;

  array_s *cells;
  array_alloc(&cells);
  array_init(cells, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  int nvc;
  size_t *vc;

  for (int i = 0; i < 3; ++i) {
    if (!I[i]) continue;

    nvc = mesh3_nvc(mesh, l[i]);
    vc = malloc(nvc*sizeof(size_t));
    mesh3_vc(mesh, l[i], vc);

    for (int j = 0; j < nvc; ++j)
      if (!array_contains(cells, &vc[j]))
        array_append(cells, &vc[j]);

    free(vc);
  }

  dbl b[4];
  size_t lc;
  for (size_t i = 0; i < array_size(cells); ++i) {
    array_get(cells, i, &lc);
    xm_in_cell |= mesh3_dbl3_in_cell(mesh, lc, xm, b);
    xp_in_cell |= mesh3_dbl3_in_cell(mesh, lc, xp, b);
    if (xm_in_cell && xp_in_cell)
      break;
  }

  /* Next, we'll pull out the boundary faces incident on these cells
     and check if the ray intersects any of them. If it does, the ray
     is unphysical. */

  array_s *bdf;
  array_alloc(&bdf);
  array_init(bdf, sizeof(size_t[3]), ARRAY_DEFAULT_CAPACITY);

  // Get the incident boundary faces
  size_t cf[4][3];
  for (size_t i = 0; i < array_size(cells); ++i) {
    array_get(cells, i, &lc);
    mesh3_cf(mesh, lc, cf);
    for (size_t j = 0; j < 4; ++j) {
      if (array_contains(bdf, cf[j]) ||
          point_in_face(utetra->lhat, cf[j])) // skip faces containing
                                              // the update point to
                                              // make the intersection
                                              // test simpler
        continue;
      if (mesh3_bdf(mesh, cf[j]))
        array_append(bdf, cf[j]);
    }
  }

  // Check for intersections with the nearby boundary faces
  bool found_bdf_isect = false;
  size_t lf[3];
  for (size_t i = 0; i < array_size(bdf); ++i) {
    array_get(bdf, i, lf);
    tri3 tri = mesh3_get_tri(mesh, lf);
    if (ray3_and_tri3_are_parallel(&ray, &tri))
      continue;
    dbl t;
    bool isect = ray3_intersects_tri3(&ray, &tri, &t);
    if (isect && t <= utetra->L) {
      found_bdf_isect = true;
      break;
    }
  }

  array_deinit(bdf);
  array_dealloc(&bdf);

  array_deinit(cells);
  array_dealloc(&cells);

  if (found_bdf_isect)
    return false;

  // If we didn't find a containing cell, we can conclude that the ray
  // is unphysical!
  if (!xm_in_cell || !xp_in_cell)
    return false;

  /**
   * Next, check and see if the point just before the end of the ray
   * lies in a cell.
   */

  size_t lhat = utetra->lhat;

  dbl xhatm[3];
  dbl3_saxpy(-t/utetra->L, utetra->x_minus_xb, utetra->x, xhatm);

  nvc = mesh3_nvc(mesh, lhat);
  vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(mesh, lhat, vc);

  bool xhatm_in_cell = false;
  for (int i = 0; i < nvc; ++i) {
    xhatm_in_cell = mesh3_dbl3_in_cell(mesh, vc[i], xhatm, b);
    if (xhatm_in_cell)
      break;
  }

  free(vc);

  return xhatm_in_cell;
}

int utetra_get_num_interior_coefs(utetra_s const *utetra) {
  dbl const atol = 1e-14;
  return (utetra->lam[0] > atol) + (utetra->lam[1] > atol)
    + (utetra->lam[0] + utetra->lam[1] < 1 - atol);
}

void utetra_get_update_inds(utetra_s const *utetra, size_t l[3]) {
  assert(utetra_update_inds_are_set(utetra));
  memcpy(l, utetra->l, sizeof(size_t[3]));
}

bool utetra_has_shadow_solution(utetra_s const *utetra, eik3_s const *eik) {
  assert(utetra_inds_are_set(utetra));

  dbl const atol = 1e-14;

  dbl b[3];
  utetra_get_bary_coords(utetra, b);

  bool has_shadow_solution = true;

  // Check if active => shadow
  for (int i = 0; i < 3; ++i)
    if (fabs(b[i]) > atol)
      has_shadow_solution &= eik3_is_shadow(eik, utetra->l[i]);

  return has_shadow_solution;
}

int utetra_get_num_shared_inds(utetra_s const *u1, utetra_s const *u2) {
  assert(utetra_inds_are_set(u1));
  assert(utetra_inds_are_set(u2));

  size_t const *l1 = u1->l, *l2 = u2->l;
  assert(l1[0] != l1[1] && l1[1] != l1[2]);
  assert(l2[0] != l2[1] && l2[1] != l2[2]);
  int num_shared_inds = (
    (l1[0] == l2[0]) + (l1[0] == l2[1]) + (l1[0] == l2[2]) +
    (l1[1] == l2[0]) + (l1[1] == l2[1]) + (l1[1] == l2[2]) +
    (l1[2] == l2[0]) + (l1[2] == l2[1]) + (l1[2] == l2[2]));
  assert(num_shared_inds <= 3);
  return num_shared_inds;
}

bool utetras_yield_same_update(utetra_s const **utetra, int n) {
  dbl const atol = 1e-14;

  dbl x[2][3];
  jet3 jet[2];

  // Prefetch the first coords and jet
  utetra_get_x(utetra[0], x[0]);
  utetra_get_jet(utetra[0], &jet[0]);

  for (int i = 1; i < n; ++i) {
    // Get the next jet and check that it's finite
    utetra_get_jet(utetra[i], &jet[1]);
    if (!jet3_is_finite(&jet[1]))
      return false;

    // Get the next coords
    utetra_get_x(utetra[i], x[1]);

    // Check if the base of the update rays coincide...
    if (dbl3_dist(x[0], x[1]) > atol)
      return false;

    // Check if the computed jets are the same...
    if (!jet3_approx_eq(&jet[0], &jet[1], atol))
      return false;

    // Swap the coords and jet that we just fetched to make way for
    // the next ones
    for (int j = 0; j < 3; ++j) SWAP(x[0][j], x[1][j]);
    SWAP(jet[0], jet[1]);
  }

  return true;
}

size_t utetra_get_l(utetra_s const *utetra) {
  assert(utetra_inds_are_set(utetra));
  return utetra->lhat;
}

void utetra_set_update_inds(utetra_s *utetra, size_t l[3]) {
  memcpy(utetra->l, l, sizeof(size_t[3]));
}

void utetra_set_l0(utetra_s *utetra, size_t l0) {
  utetra->l[0] = l0;
}

void utetra_set_l1(utetra_s *utetra, size_t l1) {
  utetra->l[1] = l1;
}

void utetra_set_l2(utetra_s *utetra, size_t l2) {
  utetra->l[2] = l2;
}

/**
 * Check whether the optimum of `u1` is incident on the base of `u2`.
 */
bool utetra_opt_inc_on_other_utetra(utetra_s const *u1, utetra_s const *u2) {
  assert(utetra_inds_are_set(u1));
  assert(utetra_inds_are_set(u2));

  dbl b1[3];
  utetra_get_bary_coords(u1, b1);

  bool zero[3] = {b1[0] == 0, b1[1] == 0, b1[2] == 0};

  int num_bd = 3 - !zero[0] - !zero[1] - !zero[2];
  assert(num_bd < 3);
  if (num_bd == 0)
    return false;

  size_t const *l1 = u1->l, *l2 = u2->l;

  bool l1_inc[3] = {
    l1[0] == l2[0] || l1[0] == l2[1] || l1[0] == l2[2],
    l1[1] == l2[0] || l1[1] == l2[1] || l1[1] == l2[2],
    l1[2] == l2[0] || l1[2] == l2[1] || l1[2] == l2[2]
  };
  if (l1_inc[0] + l1_inc[1] + l1_inc[2] == 0)
    return false;

  if (num_bd == 1) { // TODO: are the following tests just equivalent to "zero XOR l1_inc"?
    if (zero[0])
      return l1_inc[1] && l1_inc[2];
    if (zero[1])
      return l1_inc[0] && l1_inc[2];
    if (zero[2])
      return l1_inc[0] && l1_inc[1];
    assert(false);
  }

  if (num_bd == 2) {
    if (zero[1] && zero[2])
      return l1_inc[0];
    if (zero[0] && zero[2])
      return l1_inc[1];
    if (zero[0] && zero[1])
      return l1_inc[2];
    assert(false);
  }

  assert(false);
  return false;
}

void utetra_get_x(utetra_s const *u, dbl x[3]) {
  dbl b[3];
  utetra_get_bary_coords(u, b);
  dbl33_dbl3_mul(u->X, b, x);
}

void utetra_get_interior_coefs(utetra_s const *utetra, size_t *l) {
  assert(utetra_update_inds_are_set(utetra));

  dbl const atol = 1e-14;

  size_t i = 0;

  if (utetra->lam[0] > atol)
    l[i++] = utetra->l[1];

  if (utetra->lam[1] > atol)
    l[i++] = utetra->l[2];

  if (utetra->lam[0] + utetra->lam[1] < 1 - atol)
    l[i++] = utetra->l[0];
}

size_t utetra_get_active_inds(utetra_s const *utetra, size_t l[3]) {
  assert(utetra_update_inds_are_set(utetra));

  dbl const atol = 1e-14;

  dbl alpha[3];
  utetra_get_lag_mults(utetra, alpha);

  size_t num_active = 0;

  l[0] = l[1] = l[2] = NO_INDEX;
  for (int i = 0; i < 3; ++i)
    if (fabs(alpha[i]) < atol)
      l[num_active++] = utetra->l[i];

  return num_active;
}

par3_s utetra_get_parent(utetra_s const *utetra) {
  par3_s par;
  utetra_get_bary_coords(utetra, par.b);
  utetra_get_update_inds(utetra, par.l);
  return par;
}
