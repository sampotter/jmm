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

  /* The state of the vertices of the update base (i.e., `state[i]`
   * goes with index `l[i]`). */
  state_e state[3];

  int niter;

  /* The number of update indices that are `SHADOW`. Used to determine
   * what to do with the `split` updates.  */
  size_t num_shadow;

  /* "Split" updates. If any of the update indices are `SHADOW` nodes,
   * then we split the base of the update along the shadow boundary
   * and only perform the minimization over the `VALID` region. The
   * shadow boundary can either be a triangle or a quadrilateral, so
   * `split` can consist of one or two `utetra`, to keep things
   * simple. (We could implement a `upyramid` type, to handle
   * quadrilateral bases as a special case, but this isn't a high
   * priority since we can accomplish the same thing by doing two
   * tetrahedron updates and patching them together.)
   *
   * Note that when this update is split, the information stored
   * directly in this `utetra` becomes a bit redundant. All of the
   * calls to the different functions on this type will be redirected
   * to the split updates. We could try to clean this up a bit later
   * using some kind of polymorphism, but this isn't a high
   * priority. */
  utetra_s **split;
};

void utetra_alloc(utetra_s **utetra) {
  *utetra = malloc(sizeof(utetra_s));
}

void utetra_dealloc(utetra_s **utetra) {
  free(*utetra);
  *utetra = NULL;
}

/* Initialize the split updates when there's one `SHADOW` update
 * vertex. This results in a pair of split `utetra` forming a
 * quadrilateral update base. */
static void init_split_utetra_s1(utetra_s *u, eik3_s const *eik) {
  dbl const atol = 1e-15;

  /* Get indices and move the `SHADOW` index to `l[2]` */
  size_t l[3] = {u->l[0], u->l[1], u->l[2]};
  if (u->state[0] == SHADOW) SWAP(l[0], l[2]);
  if (u->state[1] == SHADOW) SWAP(l[1], l[2]);

  /* Get the cut point parameters for each edge. */
  dbl t[2];
  assert(eik3_get_cutedge_t(eik, l[0], l[2], &t[0]));
  assert(eik3_get_cutedge_t(eik, l[1], l[2], &t[1]));

  /* First, check if both cut points nearly coincide with the `VALID`
   * nodes. When this happens, reset all of the states to `SHADOW` and
   * return early. */
  // TODO: this is actually slightly wrong. Technically, when this
  // happens, the edge [x0, x1] is `VALID`, and we could conceivably
  // have a `VALID` update from this edge. One thing we could do would
  // be to replace this update with a triangle update in this
  // case. Things would start to get pretty complicated in that case,
  // though!
  if (t[0] < atol && t[1] < atol) {
    u->state[0] = u->state[1] = u->state[2] = SHADOW;
    u->num_shadow = 3;
    return;
  }

  /* Next, check if the cut points both nearly coincide with the
   * `SHADOW` vertex. When this happens, we should instead avoid
   * splitting the update, and just reset the state of the `SHADOW`
   * node to `VALID`, since this is basically the case we're dealing
   * with. */
  if (t[0] > 1 - atol && t[1] > 1 - atol) {
    u->state[0] = u->state[1] = u->state[2] = VALID;
    u->num_shadow = 0;
    return;
  }

  /* Allocate two utetra */
  u->split = malloc(2*sizeof(utetra_s*));
  utetra_alloc(&u->split[0]);
  utetra_alloc(&u->split[1]);

  mesh3_s const *mesh = eik3_get_mesh(eik);

  /* We'll initialize the first `utetra` so that its base is the
   * triangle consisting of nodes `l0`, `l1`, and the cut point
   * between `l0` and `l2`. First, we just grab nodes `l0` and
   * `l1`. We also want to arrange the nodes so that the first
   * coordinate of this and the second `utetra` agree. That is, (t, 0)
   * indexes the same point on the edge [x1, xt], shared by each split
   * `utetra`. */
  dbl x[3], Xt[3][3], dx[3];
  mesh3_copy_vert(mesh, u->lhat, x);
  mesh3_copy_vert(mesh, l[1], Xt[0]); // `l[1]` goes first! see the
  mesh3_copy_vert(mesh, l[0], Xt[2]); // note above about coordinates

  dbl const *x0 = Xt[2], *x1 = Xt[0], *x2 = mesh3_get_vert_ptr(mesh, l[2]);

  /* Get the cut point between nodes `l0` and `l2`. */
  dbl3_sub(x2, x0, dx);
  dbl3_saxpy(t[0], dx, x0, Xt[1]);

  /* Get the jets. */
  jet3 jet[3];
  jet[0] = eik3_get_jet(eik, l[1]);
  jet[2] = eik3_get_jet(eik, l[0]);
  assert(eik3_get_cutedge_jet(eik, l[0], l[2], &jet[1]));

  /* Initialize the first `utetra`... */
  utetra_init_no_inds(u->split[0], x, Xt, jet);

  /* ... and make sure the states of each base node are
   * `VALID`. (Technically, the cut points don't have a state, since
   * they aren't grid nodes, but they are in the `VALID` zone, so it
   * makes sense to label them as such.) */
  state_e *state = u->split[0]->state;
  state[0] = state[1] = state[2] = VALID;
  u->split[0]->num_shadow = 0;

  /* Initialize the second `utetra` so that its base is the triangle
   * consisting of node `l1`, the cut point between `l0` and `l2`, and
   * the cut point between `l1` and `l2`.  */
  dbl3_sub(x2, x1, dx);
  dbl3_saxpy(t[1], dx, x1, Xt[2]);

  /* Get the jet at that cut point. */
  assert(eik3_get_cutedge_jet(eik, l[1], l[2], &jet[2]));

  /* Initialize the second `utetra`... */
  utetra_init_no_inds(u->split[1], x, Xt, jet);

  /* ... and again, set the inds to `VALID`. */
  state = u->split[1]->state;
  state[0] = state[1] = state[2] = VALID;
  u->split[1]->num_shadow = 0;
}

/* Initialize the split update when there are two `SHADOW` update
 * points. This results in a single `utetra`. */
static void init_split_utetra_s2(utetra_s *u, eik3_s const *eik) {
  dbl const atol = 1e-15;

  /* Get indices and move the `VALID` index to `l[0]`. */
  size_t l[3] = {u->l[0], u->l[1], u->l[2]};
  if (u->state[1] == VALID) SWAP(l[0], l[1]);
  if (u->state[2] == VALID) SWAP(l[0], l[2]);

  /* Get the cut point parameters for each edge. */
  dbl t[2];
  assert(eik3_get_cutedge_t(eik, l[0], l[1], &t[0]));
  assert(eik3_get_cutedge_t(eik, l[0], l[2], &t[1]));

  /* Check if both of the cut points are nearly coincident with the
   * `SHADOW` nodes. When this happens, set all nodes' states to
   * `VALID` and return. */
  if (t[0] > 1 - atol && t[1] > 1 - atol) {
    u->state[0] = u->state[1] = u->state[2] = VALID;
    u->num_shadow = 0;
    return;
  }

  /* ... and if both cut points are nearly coincident wit the `VALID`
   * node, set all nodes' states to `SHADOW` and return. */
  if (t[0] < atol && t[1] < atol) {
    u->state[0] = u->state[1] = u->state[2] = SHADOW;
    u->num_shadow = 3;
    return;
  }

  mesh3_s const *mesh = eik3_get_mesh(eik);

  /* Allocate one utetra */
  u->split = malloc(sizeof(utetra_s *));
  utetra_alloc(&u->split[0]);

  dbl x[3], Xt[3][3], dx[3];

  /* Get the vertices for `l[0]` and the update point. */
  mesh3_copy_vert(mesh, u->lhat, x);
  mesh3_copy_vert(mesh, l[0], Xt[0]);

  /* Get the cut point for the `(l[0], l[1])` cut edge. */
  dbl3_sub(mesh3_get_vert_ptr(mesh, l[1]), Xt[0], dx);
  dbl3_saxpy(t[0], dx, Xt[0], Xt[1]);

  /* Get the cut point for the `(l[0], l[2])` cut edge. */
  dbl3_sub(mesh3_get_vert_ptr(mesh, l[2]), Xt[0], dx);
  dbl3_saxpy(t[1], dx, Xt[0], Xt[2]);

  /* Get the jets for each point. */
  jet3 jet[3];
  jet[0] = eik3_get_jet(eik, l[0]);
  assert(eik3_get_cutedge_jet(eik, l[0], l[1], &jet[1]));
  assert(eik3_get_cutedge_jet(eik, l[0], l[2], &jet[2]));

  /* Initialize the split `utetra`... */
  utetra_init_no_inds(u->split[0], x, Xt, jet);

  /* ... and set the state of each update vertex to `VALID`. */
  state_e *state = u->split[0]->state;
  state[0] = state[1] = state[2] = VALID;
}

/* Conditionally initialize the "split" updates. If there's more than
 * one shadow node, then we split the base of the update, restricting
 * the domain using information from the cutset. */
static void init_split_utetra(utetra_s *u, eik3_s const *eik) {
  if (u->num_shadow == 1)
    init_split_utetra_s1(u, eik);
  else if (u->num_shadow == 2)
    init_split_utetra_s2(u, eik);
  else
    u->split = NULL;
}

bool utetra_init_from_eik3(utetra_s *cf, eik3_s const *eik,
                           size_t l, size_t l0, size_t l1, size_t l2) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  jet3 const *jet = eik3_get_jet_ptr(eik);

  if (!utetra_init_from_ptrs(cf, mesh, jet, l, l0, l1, l2))
    return false;

  state_e const *state = eik3_get_state_ptr(eik);
  cf->state[0] = state[l0];
  cf->state[1] = state[l1];
  cf->state[2] = state[l2];

  cf->num_shadow = 0;
  for (size_t i = 0; i < 3; ++i)
    cf->num_shadow += cf->state[i] == SHADOW;

  init_split_utetra(cf, eik);

  // TODO: would be good to initialize the containing utetra with
  // dummy values when the update is split so it's a bit more obvious
  // while debugging... (set things to NAN etc...)

  return true;
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
  cf->f = INFINITY;

  /* Initialize some variables with dummy values. */
  cf->lam[0] = cf->lam[1] = NAN;
  cf->L = NAN;

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

  cf->state[0] = cf->state[1] = cf->state[2] = UNKNOWN;

  cf->num_shadow = 0;
  cf->split = NULL;

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

static bool is_split(utetra_s const *cf) {
  return cf->num_shadow == 1 || cf->num_shadow == 2;
}

static size_t num_split(utetra_s const *cf) {
  return cf->split == NULL ? 0 : 3 - cf->num_shadow;
}

void utetra_deinit(utetra_s *u) {
  if (is_split(u)) {
    assert(u->split != NULL);
    for (size_t i = 0; i < num_split(u); ++i) {
      utetra_deinit(u->split[i]);
      utetra_dealloc(&u->split[i]);
    }
  } else {
    assert(u->split == NULL);
  }
}

bool utetra_is_degenerate(utetra_s const *cf) {
  // Check if the point being updated lies in the plane spanned by by
  // x0, x1, and x2. If it does, the update is degenerate.
  dbl const *x[4] = {cf->x, cf->Xt[0], cf->Xt[1], cf->Xt[2]};
  return points_are_coplanar(x);
}

static void set_lambda(utetra_s *cf, dbl const lam[2]) {
  assert(!is_split(cf));

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

  // Compute the trace and determinant of the Hessian
  dbl tr = cf->H[0][0] + cf->H[1][1];
  dbl det = cf->H[0][0]*cf->H[1][1] - cf->H[0][1]*cf->H[1][0];
  assert(tr != 0 && det != 0);

  // Conditionally perturb the Hessian
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

static void step(utetra_s *cf) {
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
  set_lambda(cf, lam1);
  while (cf->f > f + beta*c1_times_g_dot_p + atol) {
    beta *= 0.9;
    dbl2_saxpy(beta, p, lam, lam1);
    set_lambda(cf, lam1);
  }

  ++cf->niter;
}

/**
 * Do a tetrahedron update starting at `lam`, writing the result to
 * `jet`. If `lam` is `NULL`, then the first iterate will be selected
 * automatically.
 */
void utetra_solve(utetra_s *cf, dbl const *lam) {
  /* When `cf` is a split update, we solve its sub-updates instead of
   * it, so `cf`'s data goes more-or-less unused. */
  if (is_split(cf)) {
    for (size_t i = 0; i < num_split(cf); ++i)
      utetra_solve(cf->split[i], NULL);
    return;
  }

  cf->niter = 0;

  if (lam == NULL)
    set_lambda(cf, (dbl[2]) {0, 0}); // see, automatically!
  else
    set_lambda(cf, lam);

  // I think we're actually supposed to use cf->g here instead of
  // cf->p, but maybe cf->p is more appropriate for constraint
  // Newton's method, since it takes into account the constraints
  // while cf->p doesn't.
  int const max_niter = 100;
  dbl const atol = 1e-15, rtol = 1e-15, tol = rtol*dbl2_norm(cf->p) + atol;
  for (int _ = 0; _ < max_niter; ++_) {
    if (dbl2_norm(cf->p) <= tol)
      break;
    step(cf);
  }
}

static void get_b(utetra_s const *cf, dbl b[3]) {
  assert(!is_split(cf));
  assert(!isnan(cf->lam[0]) && !isnan(cf->lam[1]));
  b[0] = 1 - cf->lam[0] - cf->lam[1];
  b[1] = cf->lam[0];
  b[2] = cf->lam[1];
}

/* If `u` is split, return the split `utetra` with the smaller
 * value. Otherwise, return `NULL`. This is a convenience function
 * used to speed up dispatching function calls below. */
static utetra_s *split_utetra_select(utetra_s const *u) {
  if (is_split(u))
    return num_split(u) == 1 || u->split[0]->f < u->split[1]->f ?
      u->split[0] : u->split[1];
  else
    return NULL;
}

dbl utetra_get_value(utetra_s const *cf) {
  if (is_split(cf))
    return num_split(cf) == 1 ?
      cf->split[0]->f :
      fmin(cf->split[0]->f, cf->split[1]->f);
  else
    return cf->f;
}

void utetra_get_jet(utetra_s const *cf, jet3 *jet) {
  utetra_s const *u = split_utetra_select(cf);
  if (u) {
    utetra_get_jet(u, jet);
    return;
  }

  jet->f = cf->f;
  jet->fx = cf->x_minus_xb[0]/cf->L;
  jet->fy = cf->x_minus_xb[1]/cf->L;
  jet->fz = cf->x_minus_xb[2]/cf->L;
}

/**
 * Compute the Lagrange multipliers for the constraint optimization
 * problem corresponding to this type of update
 */
static void get_lag_mults(utetra_s const *cf, dbl alpha[3]) {
  assert(!is_split(cf));
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

static bool split_utetra_has_interior_point_solution(utetra_s const *cf) {
  dbl const atol = 1e-14;

  utetra_s **u = cf->split;

  if (num_split(cf) == 1)
    return utetra_has_interior_point_solution(u[0]);

  dbl alpha[2][3];
  for (size_t i = 0; i < 2; ++i)
    get_lag_mults(u[i], alpha[i]);

  jet3 jet[2];
  for (size_t i = 0; i < 2; ++i)
    utetra_get_jet(u[i], &jet[i]);

  if (fabs(jet[0].f - jet[1].f) < atol)
    return atol <= u[0]->lam[0] && u[0]->lam[0] <= 1 - atol &&
      fabs(u[0]->lam[0] - u[1]->lam[0]) < atol &&
      u[0]->lam[1] < atol && u[1]->lam[1] < atol;

  if (jet[0].f < jet[1].f)
    return fabs(alpha[1][0]) < atol && fabs(alpha[1][1]) < atol &&
      dbl3_maxnorm(alpha[0]) < atol;

  if (jet[0].f > jet[1].f)
    return fabs(alpha[0][0]) < atol && fabs(alpha[0][1]) < atol &&
      dbl3_maxnorm(alpha[1]) < atol;

  assert(false);
}

bool utetra_has_interior_point_solution(utetra_s const *cf) {
  if (is_split(cf))
    return split_utetra_has_interior_point_solution(cf);

  dbl alpha[3];
  get_lag_mults(cf, alpha);
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
    dbl f1 = utetra_get_value(u1), f2 = utetra_get_value(u2);
    if (f1 < f2) {
      return -1;
    } else if (f1 > f2) {
      return 1;
    } else {
      return 0;
    }
  }
}

bool utetra_adj_are_optimal(utetra_s const *u1, utetra_s const *u2) {
  assert(!is_split(u1));
  assert(!is_split(u2));

  dbl const atol = 1e-15;

  return fabs(u1->lam[0] - u2->lam[0]) <= atol
    && fabs(u1->lam[1]) <= atol
    && fabs(u2->lam[1]) <= atol
    && fabs(utetra_get_value(u1) - utetra_get_value(u2)) <= atol;
}

static ray3 get_ray(utetra_s const *utetra) {
  utetra_s const *u = split_utetra_select(utetra);
  if (u)
    return get_ray(u);

  ray3 ray;
  dbl b[3];
  get_b(utetra, b);
  dbl33_dbl3_mul(utetra->X, b, ray.org);
  dbl3_normalized(utetra->x_minus_xb, ray.dir);
  return ray;
}

static bool update_inds_are_set(utetra_s const *utetra) {
  return !(
    utetra->l[0] == (size_t)NO_INDEX ||
    utetra->l[1] == (size_t)NO_INDEX ||
    utetra->l[2] == (size_t)NO_INDEX);
}

static bool all_inds_are_set(utetra_s const *utetra) {
  return utetra->lhat != (size_t)NO_INDEX && update_inds_are_set(utetra);
}

static void get_interior_coefs(utetra_s const *utetra, size_t *l) {
  assert(update_inds_are_set(utetra));

  dbl const atol = 1e-14;

  size_t i = 0;

  if (utetra->lam[0] > atol)
    l[i++] = utetra->l[1];

  if (utetra->lam[1] > atol)
    l[i++] = utetra->l[2];

  if (utetra->lam[0] + utetra->lam[1] < 1 - atol)
    l[i++] = utetra->l[0];
}

static dbl get_L(utetra_s const *u) {
  utetra_s *u_split = split_utetra_select(u);
  return u_split ? get_L(u_split) : u->L;
}

bool utetra_update_ray_is_physical(utetra_s const *utetra, eik3_s const *eik) {
  assert(all_inds_are_set(utetra));

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

  if (num_int == 2) {
    size_t l_int[3] = {NO_INDEX, NO_INDEX, NO_INDEX};
    get_interior_coefs(utetra, l_int);
    if (mesh3_is_nondiff_boundary_edge(mesh, l_int))
      return false;
  }

  if (num_int == 3 && mesh3_bdf(mesh, l))
    return false;

  // TODO: the following section where we check to see if the stuff
  // below gives "an interior ray" can be wrapped up and reused for
  // both this and the corresponding section in utri.c...

  ray3 ray = get_ray(utetra);

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
  dbl xm[3], xp[3], t = mesh3_get_min_tetra_alt(mesh)/2; // TODO: rename t...
  ray3_get_point(&ray, -t, xm);
  ray3_get_point(&ray, t, xp);

  array_s *cells;
  array_alloc(&cells);
  array_init(cells, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* Fill `cells` with all of the cells incident on any of the
   * neighboring update indices. */
  for (int i = 0; i < 3; ++i) {
    int nvc = mesh3_nvc(mesh, l[i]);
    size_t *vc = malloc(nvc*sizeof(size_t));
    mesh3_vc(mesh, l[i], vc);
    for (int j = 0; j < nvc; ++j)
      if (!array_contains(cells, &vc[j]))
        array_append(cells, &vc[j]);
    free(vc);
  }

  bool xm_in_cell = false, xp_in_cell = false;
  for (size_t i = 0, lc; i < array_size(cells); ++i) {
    array_get(cells, i, &lc);
    xm_in_cell |= mesh3_dbl3_in_cell(mesh, lc, xm, NULL);
    xp_in_cell |= mesh3_dbl3_in_cell(mesh, lc, xp, NULL);
    if (xm_in_cell && xp_in_cell)
      break;
  }
  if (!xm_in_cell || !xp_in_cell)
    return false;

  /* Next, we'll pull out the boundary faces incident on these cells
   * and check if the ray intersects any of them. If it does, the ray
   * is unphysical. */

  array_s *bdf;
  array_alloc(&bdf);
  array_init(bdf, sizeof(size_t[3]), ARRAY_DEFAULT_CAPACITY);

  // Get the incident boundary faces
  size_t cf[4][3];
  for (size_t i = 0, lc; i < array_size(cells); ++i) {
    array_get(cells, i, &lc);
    mesh3_cf(mesh, lc, cf);
    for (size_t j = 0; j < 4; ++j) {
      if (array_contains(bdf, cf[j]))
        continue;
      /* skip faces containing the update point to make the
       * intersection test simpler */
      if (point_in_face(utetra->lhat, cf[j]))
        continue;
      if (mesh3_bdf(mesh, cf[j]))
        array_append(bdf, cf[j]);
    }
  }

  dbl L = get_L(utetra);

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
    if (isect && t <= L) {
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

  /* Next, check and see if the point just before the end of the ray
   * lies in a cell. */

  dbl xhatm[3];
  ray3_get_point(&ray, L - t, xhatm);

  int nvc = mesh3_nvc(mesh, utetra->lhat);
  size_t *vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(mesh, utetra->lhat, vc);

  bool xhatm_in_cell = false;
  for (int i = 0; i < nvc; ++i) {
    xhatm_in_cell = mesh3_dbl3_in_cell(mesh, vc[i], xhatm, NULL);
    if (xhatm_in_cell)
      break;
  }

  free(vc);

  return xhatm_in_cell;
}

static int split_utetra_get_num_interior_coefs(utetra_s const *utetra) {
  dbl const atol = 1e-15;

  size_t n = num_split(utetra);

  dbl f[2];
  f[0] = utetra->split[0]->f;
  f[1] = n == 2 ? utetra->split[1]->f : NAN;

  /* Select utetra */
  utetra_s const *u;
  if (n == 1 || f[0] < f[1])
    u = utetra->split[0];
  else
    u = utetra->split[1];

  dbl b[3];
  get_b(u, b);

  /* If we're dealing with a quad split (`n == 2`) and the update
   * point is in the interior of the quad, then we should return 4. */
  if (n == 2 &&
      b[0] > atol && b[1] > atol && b[2] <= atol)
    return 4;

  /* Compute the number of interior coefficients the normal way and
   * return them. */
  return (b[0] > atol) + (b[1] > atol) + (b[2] > atol);
}

int utetra_get_num_interior_coefs(utetra_s const *utetra) {
  if (is_split(utetra))
    return split_utetra_get_num_interior_coefs(utetra);

  dbl const atol = 1e-15;
  dbl b[3];
  get_b(utetra, b);
  return (b[0] > atol) + (b[1] > atol) + (b[2] > atol);
}

bool utetra_has_shadow_solution(utetra_s const *utetra) {
  /* First, handle the easy cases... all `VALID` or all `SHADOW`
   * nodes. */

  if (utetra->num_shadow == 0)
    return false;

  if (utetra->num_shadow == 3)
    return true;

  /* Next, handle the hard case (split updates). */

  assert(is_split(utetra));

  dbl const atol = 1e-14;

  utetra_s **u = utetra->split;

  size_t n = num_split(utetra);
  if (n == 1 || u[1]->f > u[0]->f) {
    dbl b0 = 1 - u[0]->lam[0] - u[0]->lam[1];
    return b0 < atol;
  } else {
    return u[0]->lam[0] >= 1 - atol;
  }
}

static size_t get_shared_inds(utetra_s const *u1, utetra_s const *u2, size_t *l) {
  assert(update_inds_are_set(u1));
  assert(update_inds_are_set(u2));

  size_t const *l1 = u1->l, *l2 = u2->l;

  size_t i = 0;

  for (size_t j = 0; j < 3; ++j)
    if (l1[j] == l2[0] || l1[j] == l2[1] || l1[j] == l2[2])
      l[i++] = l1[j];

  return i; // == num_shared_inds
}

static bool get_point_for_index(utetra_s const *utetra, size_t l, dbl x[3]) {
  assert(update_inds_are_set(utetra));

  for (int i = 0; i < 3; ++i)
    if (utetra->l[i] == l) {
      memcpy(x, utetra->Xt[i], sizeof(dbl[3]));
      return true;
    }

  return false;
}

static bool get_op_ind(utetra_s const *utetra, size_t const le[2], size_t *l) {
  bool found[3] = {false, false, false};

  for (int i = 0; i < 3; ++i)
    if (utetra->l[i] == le[0] || utetra->l[i] == le[1])
      found[i] = true;

  if (found[0] + found[1] + found[2] == 2) {
    for (int i = 0; i < 3; ++i)
      if (!found[i]) {
        *l = utetra->l[i];
        return true;
      }
    assert(false);
  }

  return false;
}

static size_t get_num_equal(utetra_s const **u, size_t n) {
  dbl const atol = 1e-14;

  dbl x[3], y[3];

  /* The first thing to check is whether the jet computed by each
   * update and the start of the update ray parametrized by each
   * update is the same. This is necessary but not sufficient. */

  jet3 jet[2];

  // Prefetch the first coords and jet
  utetra_get_x(u[0], x);
  utetra_get_jet(u[0], &jet[0]);

  size_t neq = 1;

  // First, check that the update data for each utetra is the
  // same. This is cheaper than the topological check that follows.
  for (neq = 1; neq < n; ++neq) {
    // Get the next jet and check that it's finite
    utetra_get_jet(u[neq], &jet[1]);
    if (!jet3_is_finite(&jet[1]))
      break;

    // Get the next coords
    utetra_get_x(u[neq], y);

    // Check if the base of the update rays coincide...
    if (dbl3_dist(x, y) > atol)
      break;

    // Check if the computed jets are the same...
    if (!jet3_approx_eq(&jet[0], &jet[1], atol))
      break;

    // Swap the coords and jet that we just fetched to make way for
    // the next ones
    for (int j = 0; j < 3; ++j)
      SWAP(x[j], y[j]);
    SWAP(jet[0], jet[1]);
  }

  return neq;
}

static size_t count_unique_inds(utetra_s const **u, size_t n) {
  array_s *l;

  array_alloc(&l);
  array_init(l, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (!array_contains(l, &u[i]->l[j]))
        array_append(l, &u[i]->l[j]);

  size_t num_unique = array_size(l);

  array_deinit(l);
  array_dealloc(&l);

  return num_unique;
}

static void get_shared_and_op_inds_3(utetra_s const **u, size_t *l_shared,
                                     size_t l_op[3]) {
  assert(count_unique_inds(u, 3) == 4);

  size_t tmp1[2], tmp2[2];
  get_shared_inds(u[0], u[1], tmp1);
  get_shared_inds(u[0], u[2], tmp2);

  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      if (tmp1[i] == tmp2[j]) {
        *l_shared = tmp1[i];
        break;
      }

  size_t j = 0;
  for (size_t i = 0; i < 3; ++i)
    if (u[0]->l[i] != *l_shared)
      l_op[j++] = u[0]->l[i];
  assert(j == 2);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j)
      if (u[1]->l[i] != l_op[j]) {
        l_op[2] = u[1]->l[i];
        break;
      }
}

/* Checks whether `n` tetrahedron updates stored in `utetra[i]` yield
 * the same update. This is a pretty complicated check... Need to
 * explain the cases this handles a bit better. */
bool utetras_yield_same_update(utetra_s const **u, size_t n) {
  assert(n > 1);

  size_t neq = get_num_equal(u, n);

  if (neq == 1)
    return false;

  /* We want to check whether the mesh comprised of the `n` update
   * bases (triangles) are "thick" from the perspective of the update
   * ray. Basically, what we're trying to check is whether, after
   * projecting the union of the update bases into the plane normal to
   * the update ray, the ray hits an interior point or not. E.g., we
   * want to rule out rays that only *graze* the update triangles, or
   * weird sets of updates that aren't manifold at the intersection
   * point. As a reminder, this function should only be called when we
   * were unable to find a single update with an interior point
   * solution. */

  jet3 jet;
  utetra_get_jet(u[0], &jet);

  if (neq == 2) {
    // Get the indices of the shared edge, and validate that this edge
    // is incident on the base of each tetrahedron update
    size_t l_shared[2];
    size_t num_shared = get_shared_inds(u[0], u[1], l_shared);
    assert(num_shared < 3);

    if (num_shared == 1)
      return false;

    dbl x[3], dx[3], y[3], dy[3], normal[3];

    // Get the vectors defining the line spanned by the shared edge
    get_point_for_index(u[0], l_shared[0], y);
    get_point_for_index(u[0], l_shared[1], dy);
    dbl3_sub_inplace(dy, y);

    // Get the normal for the plane spanned by the edge and update ray,
    // taking advantage of the fact that all the update rays can be
    // assumed to be equal at this point
    dbl3_cross(dy, &jet.fx, normal);
    dbl3_normalize(normal); // (probably overkill)

    // For each tetrahedron update, pull out the points corresponding to
    // index that isn't incident on the shared edge, and translate them
    // so that the edge passes through the origin. Then, compute the dot
    // product between this point and the normal vector `n`, keeping
    // track of the minimum and maximum dot products.
    dbl amin = INFINITY, amax = -INFINITY;
    for (size_t i = 0, l_op; i < neq; ++i) {
      get_op_ind(u[i], l_shared, &l_op);
      get_point_for_index(u[i], l_op, x);
      dbl3_sub(x, y, dx);
      dbl a = dbl3_dot(normal, dx);
      amin = fmin(amin, a);
      amax = fmax(amax, a);
    }

    // We want to check that at least two of the "opposite update
    // indices" lie on either side of the plane spanned by the edge and
    // the update ray. The dot products `amin` and `amax` are the
    // extreme dot products values, so we just check this here.
    return amin < 0 && 0 < amax;
  }

  if (neq == 3) {
    size_t num_unique = count_unique_inds(u, neq);
    assert(num_unique >= 4);
    if (num_unique > 4)
      return false;

    size_t l_shared, l_op[3];
    get_shared_and_op_inds_3(u, &l_shared, l_op);

    /* Get the `tri3` indexed by `l_op`. This is a little circuitous
     * since we don't know the provenance of each index in `l_op` at
     * this point, so we need to do a quick search...  */
    tri3 tri;
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        if (get_point_for_index(u[j], l_op[i], tri.v[i]))
          continue;

    /* Get the update ray with which we'll try to intersect `tri`. */
    ray3 ray = get_ray(u[0]);

    dbl unused;

    /* First check the default orientation of the update ray... */
    if (ray3_intersects_tri3(&ray, &tri, &unused))
      return true;

    /* ... and if that doesn't intersect `tri`, flip it around and try
     * the other orientation, as well. We do this because what we
     * *really* want to do is intersect the line spanned by `ray` with
     * `tri`. */
    // TODO: implement a `line_intersects_tri3` function to simplify
    // this code...
    dbl3_negate(ray.dir);
    return ray3_intersects_tri3(&ray, &tri, &unused);
  }

  assert(false);
}

size_t utetra_get_l(utetra_s const *utetra) {
  assert(utetra->lhat != (size_t)NO_INDEX);
  return utetra->lhat;
}

void utetra_set_update_inds(utetra_s *utetra, size_t l[3]) {
  memcpy(utetra->l, l, sizeof(size_t[3]));
}

static tri3 get_tri(utetra_s const *u) {
  utetra_s const *u_split = split_utetra_select(u);
  if (u_split)
    return get_tri(u_split);

  tri3 tri;
  memcpy(tri.v, u->Xt, sizeof(dbl[3][3]));
  return tri;
}

/**
 * Check whether the optimum of `u` is incident on the base of
 * `u_other`.
 */
bool utetra_opt_inc_on_other_utetra(utetra_s const *u, utetra_s const *u_other) {
  dbl const atol = 1e-14;
  dbl xb[3];
  utetra_get_x(u, xb);
  tri3 tri = get_tri(u_other);
  return tri3_dist(&tri, xb) < atol;
}

void utetra_get_x(utetra_s const *u, dbl x[3]) {
  utetra_s const *u_split = split_utetra_select(u);
  if (u_split)
    return utetra_get_x(u_split, x);

  dbl b[3];
  get_b(u, b);
  dbl33_dbl3_mul(u->X, b, x);
}

size_t utetra_get_active_inds(utetra_s const *utetra, size_t l[3]) {
  assert(!is_split(utetra));
  assert(update_inds_are_set(utetra));

  dbl const atol = 1e-14;

  dbl alpha[3];
  get_lag_mults(utetra, alpha);

  size_t num_active = 0;

  l[0] = l[1] = l[2] = NO_INDEX;
  for (int i = 0; i < 3; ++i)
    if (fabs(alpha[i]) < atol)
      l[num_active++] = utetra->l[i];

  return num_active;
}

par3_s utetra_get_parent(utetra_s const *utetra) {
  par3_s par;

  /* If this is a split update, we pull out the optimum from that
   * update and then compute its barycentric coordinates with respect
   * to the base of the containing update (`utetra`). */
  if (is_split(utetra)) {
    utetra_s const *u = split_utetra_select(utetra);
    dbl x[3]; utetra_get_x(u, x);
    tri3 tri = get_tri(utetra); // tri from *containing* update
    tri3_get_bary_coords(&tri, x, par.b);
    dbl3_normalize1(par.b);
  } else {
    get_b(utetra, par.b);
  }

  assert(dbl3_valid_bary_coord(par.b));

  memcpy(par.l, utetra->l, sizeof(size_t[3]));

  return par;
}

#if JMM_TEST
void utetra_step(utetra_s *u) {
  step(u);
}

void utetra_get_lambda(utetra_s const *u, dbl lam[2]) {
  assert(!is_split(u));
  lam[0] = u->lam[0];
  lam[1] = u->lam[1];
}

void utetra_set_lambda(utetra_s *u, dbl const lam[2]) {
  set_lambda(u, lam);
}

size_t utetra_get_num_iter(utetra_s const *u) {
  assert(!is_split(u));
  return u->niter;
}
#endif
