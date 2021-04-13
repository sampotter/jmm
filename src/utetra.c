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

utetra_spec_s utetra_spec_empty() {
  return (utetra_spec_s) {
    .eik = NULL,
    .lhat = NO_INDEX,
    .l = {NO_INDEX, NO_INDEX, NO_INDEX},
    .state = {UNKNOWN, UNKNOWN, UNKNOWN},
    .xhat = {NAN, NAN, NAN},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
    }
  };
}

utetra_spec_s utetra_spec_from_eik_and_inds(eik3_s const *eik, size_t l,
                                            size_t l0, size_t l1, size_t l2) {
  state_e const *state = eik3_get_state_ptr(eik);
  return (utetra_spec_s) {
    .eik = eik,
    .lhat = l,
    .l = {l0, l1, l2},
    .state = {state[l0], state[l1], state[l2]},
    .xhat = {NAN, NAN, NAN},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
    }
  };
}

utetra_spec_s utetra_spec_from_eik_without_l(eik3_s const *eik, dbl const x[3],
                                             size_t l0, size_t l1, size_t l2) {
  state_e const *state = eik3_get_state_ptr(eik);
  return (utetra_spec_s) {
    .eik = eik,
    .lhat = NO_INDEX,
    .l = {l0, l1, l2},
    .state = {state[l0], state[l1], state[l2]},
    .xhat = {x[0], x[1], x[2]},
    .x = {
      {NAN, NAN, NAN},
      {NAN, NAN, NAN},
      {NAN, NAN, NAN}
    },
    .jet = {
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN},
      {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN}
    }
  };
}

utetra_spec_s utetra_spec_from_ptrs(mesh3_s const *mesh, jet3 const *jet,
                                    size_t l, size_t l0, size_t l1, size_t l2) {
  utetra_spec_s spec = {
    .eik = NULL,
    .lhat = NO_INDEX,
    .l = {NO_INDEX, NO_INDEX, NO_INDEX},
    .state = {UNKNOWN, UNKNOWN, UNKNOWN},
    .jet = {jet[l0], jet[l1], jet[l2]}
  };

  mesh3_copy_vert(mesh, l, spec.xhat);
  mesh3_copy_vert(mesh, l0, spec.x[0]);
  mesh3_copy_vert(mesh, l1, spec.x[1]);
  mesh3_copy_vert(mesh, l2, spec.x[2]);

  return spec;
}

struct utetra {
  dbl lam[2]; // Current iterate
  dbl f;
  dbl g[2];
  dbl H[2][2];
  dbl p[2]; // Newton step
  dbl x_minus_xb[3];
  dbl L;
  int niter;

  size_t lhat, l[3];

  /* The state of the vertices of the update base (i.e., `state[i]`
   * goes with index `l[i]`). */
  state_e state[3];

  dbl x[3]; // x[l]
  dbl X[3][3]; // X = [x[l0] x[l1] x[l2]]
  dbl Xt[3][3]; // X'
  dbl XtX[3][3]; // X'*X

  // B-coefs for 9-point triangle interpolation T on base of update
  bb32 T;

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

/* Initialize the split update when there are two `SHADOW` update
 * points. This results in a single `utetra`. */
static void init_split_s1(utetra_s *u_par, eik3_s const *eik) {
  dbl const atol = 1e-15;

  size_t l[3];
  memcpy(l, u_par->l, sizeof(size_t[3]));

  /* Move the `VALID` index to `l[0]`. */
  if (u_par->state[1] == VALID) SWAP(l[0], l[1]);
  if (u_par->state[2] == VALID) SWAP(l[0], l[2]);

  /* Get the cut point parameters for each edge. */
  dbl t[2];
  assert(eik3_get_cutedge_t(eik, l[0], l[1], &t[0]));
  assert(eik3_get_cutedge_t(eik, l[0], l[2], &t[1]));

  /* Check if both of the cut points are nearly coincident with the
   * `SHADOW` nodes. When this happens, set all nodes' states to
   * `VALID` and return. */
  if (t[0] > 1 - atol && t[1] > 1 - atol) {
    u_par->state[0] = u_par->state[1] = u_par->state[2] = VALID;
    u_par->num_shadow = 0;
    u_par->split = NULL;
    return;
  }

  /* ... and if both cut points are nearly coincident wit the `VALID`
   * node, set all nodes' states to `SHADOW` and return. */
  if (t[0] < atol && t[1] < atol) {
    u_par->state[0] = u_par->state[1] = u_par->state[2] = SHADOW;
    u_par->num_shadow = 3;
    u_par->split = NULL;
    return;
  }

  utetra_spec_s spec = utetra_spec_empty();

  spec.eik = eik;
  spec.state[0] = spec.state[1] = spec.state[2] = VALID;

  if (u_par->lhat == (size_t)NO_INDEX) {
    assert(dbl3_isfinite(u_par->x));
    dbl3_copy(u_par->x, spec.xhat);
  } else {
    spec.lhat = u_par->lhat;
  }

  mesh3_s const *mesh = eik3_get_mesh(spec.eik);

  mesh3_copy_vert(mesh, l[0], spec.x[0]);

  dbl dx[3];

  /* Get the cut point for the `(l[0], l[1])` cut edge. */
  dbl3_sub(mesh3_get_vert_ptr(mesh, l[1]), spec.x[0], dx);
  dbl3_saxpy(t[0], dx, spec.x[0], spec.x[1]);

  /* Get the cut point for the `(l[0], l[2])` cut edge. */
  dbl3_sub(mesh3_get_vert_ptr(mesh, l[2]), spec.x[0], dx);
  dbl3_saxpy(t[1], dx, spec.x[0], spec.x[2]);

  /* Get the jets for each point. */
  spec.jet[0] = eik3_get_jet(eik, l[0]);
  assert(eik3_get_cutedge_jet(eik, l[0], l[1], &spec.jet[1]));
  assert(eik3_get_cutedge_jet(eik, l[0], l[2], &spec.jet[2]));

  /* Set up the single split `utetra` */
  u_par->split = malloc(sizeof(utetra_s *));
  utetra_alloc(&u_par->split[0]);
  utetra_init(u_par->split[0], &spec);
}

/* Initialize the split updates when there's one `SHADOW` update
 * vertex. This results in a pair of split `utetra` forming a
 * quadrilateral update base. */
static void init_split_s2(utetra_s *u_par, eik3_s const *eik) {
  dbl const atol = 1e-15;

  /* Get indices and move the `SHADOW` index to `l[2]` */
  size_t l[3]; memcpy(l, u_par->l, sizeof(size_t[3]));
  if (u_par->state[0] == SHADOW) SWAP(l[0], l[2]);
  if (u_par->state[1] == SHADOW) SWAP(l[1], l[2]);

  /* Get the cut point parameters for each edge. */
  dbl t[2];
  assert(eik3_get_cutedge_t(eik, l[0], l[2], &t[0]));
  assert(eik3_get_cutedge_t(eik, l[1], l[2], &t[1]));

  // TODO: the next two sections are slightly wrong. Technically, when
  // this happens, the edge [x0, x1] is `VALID`, and we could
  // conceivably have a `VALID` update from this edge. One thing we
  // could do would be to replace this update with a triangle update
  // in this case. Things would start to get pretty complicated in
  // that case, though!

  /* First, check if both cut points nearly coincide with the `VALID`
   * nodes. When this happens, reset all of the states to `SHADOW` and
   * return early. */
  if (t[0] < atol && t[1] < atol) {
    u_par->state[0] = u_par->state[1] = u_par->state[2] = SHADOW;
    u_par->num_shadow = 3;
    u_par->split = NULL;
    return;
  }

  /* Next, check if the cut points both nearly coincide with the
   * `SHADOW` vertex. When this happens, we should instead avoid
   * splitting the update, and just reset the state of the `SHADOW`
   * node to `VALID`, since this is basically the case we're dealing
   * with. */
  if (t[0] > 1 - atol && t[1] > 1 - atol) {
    u_par->state[0] = u_par->state[1] = u_par->state[2] = VALID;
    u_par->num_shadow = 0;
    u_par->split = NULL;
    return;
  }

  utetra_spec_s spec[2] = {utetra_spec_empty(), utetra_spec_empty()};

  for (size_t i = 0; i < 2; ++i) {
    spec[i].eik = eik;
    spec[i].state[0] = spec[i].state[1] = spec[i].state[2] = VALID;
    spec[i].lhat = u_par->lhat;
  }

  /* If we weren't able to set each `spec[i].lhat` correctly above,
   * set `spec[i].xhat` now. */
  if (u_par->lhat == (size_t)NO_INDEX) {
    assert(dbl3_isfinite(u_par->x));
    for (size_t i = 0; i < 2; ++i)
      dbl3_copy(u_par->x, spec[i].xhat);
  }

  mesh3_s const *mesh = eik3_get_mesh(eik);

  dbl const *x[3];
  for (size_t i = 0; i < 3; ++i)
    x[i] = mesh3_get_vert_ptr(mesh, l[i]);

  dbl dx[3];

  /* We'll initialize the first `utetra` so that its base is the
   * triangle consisting of nodes `l0`, `l1`, and the cut point
   * between `l0` and `l2`. First, we just grab nodes `l0` and
   * `l1`. We also want to arrange the nodes so that the first
   * coordinate of this and the second `utetra` agree. That is, (t, 0)
   * indexes the same point on the edge [x1, xt], shared by each split
   * `utetra`. */
  dbl3_copy(x[1], spec[0].x[0]); // `l[1]` goes first! see the note
  dbl3_copy(x[0], spec[0].x[2]); // above about indexing

  /* Get the cut point between nodes `l0` and `l2`. */
  dbl3_sub(x[2], x[0], dx);
  dbl3_saxpy(t[0], dx, x[0], spec[0].x[1]);

  /* Get the jets. */
  spec[0].jet[0] = eik3_get_jet(eik, l[1]);
  spec[0].jet[2] = eik3_get_jet(eik, l[0]);
  assert(eik3_get_cutedge_jet(eik, l[0], l[2], &spec[0].jet[1]));

  /* Initialize the second `utetra` so that its base is the triangle
   * consisting of node `l1`, the cut point between `l0` and `l2`, and
   * the cut point between `l1` and `l2`.  */
  dbl3_copy(spec[0].x[0], spec[1].x[0]);
  dbl3_copy(spec[0].x[1], spec[1].x[1]);
  dbl3_sub(x[2], x[1], dx);
  dbl3_saxpy(t[1], dx, x[1], spec[1].x[2]);

  /* Get the jet at that cut point. */
  spec[1].jet[0] = spec[0].jet[0];
  spec[1].jet[1] = spec[0].jet[1];
  assert(eik3_get_cutedge_jet(eik, l[1], l[2], &spec[1].jet[2]));

  /* Set up the pair of split `utetra` before returning */
  u_par->split = malloc(2*sizeof(utetra_s*));
  for (size_t i = 0; i < 2; ++i) {
    utetra_alloc(&u_par->split[i]);
    utetra_init(u_par->split[i], &spec[i]);
  }
}

bool utetra_init(utetra_s *u, utetra_spec_s const *spec) {
  /* First, validate the spec */

  bool passed_lhat = spec->lhat != (size_t)NO_INDEX;

  bool passed_xhat = dbl3_isfinite(spec->xhat);
  assert(passed_lhat ^ passed_xhat); // pass exactly one of these

  bool passed_l0 = spec->l[0] != (size_t)NO_INDEX;
  bool passed_l1 = spec->l[1] != (size_t)NO_INDEX;
  bool passed_l2 = spec->l[2] != (size_t)NO_INDEX;
  bool passed_l = passed_l0 && passed_l1 && passed_l2;
  if (passed_l0 || passed_l1 || passed_l2)
    assert(passed_l);

  bool passed_x0 = dbl3_isfinite(spec->x[0]);
  bool passed_x1 = dbl3_isfinite(spec->x[1]);
  bool passed_x2 = dbl3_isfinite(spec->x[2]);
  bool passed_x = passed_x0 && passed_x1 && passed_x2;
  if (passed_x0 || passed_x1 || passed_x2)
    assert(passed_x);

  assert(passed_l ^ passed_x); // pass exactly one of these

  bool passed_state0 = spec->state[0] != UNKNOWN;
  bool passed_state1 = spec->state[1] != UNKNOWN;
  bool passed_state2 = spec->state[2] != UNKNOWN;
  bool passed_state = passed_state0 && passed_state1 && passed_state2;
  if (passed_state0 || passed_state1 || passed_state2)
    assert(passed_state);

  bool passed_jet0 = jet3_is_finite(&spec->jet[0]);
  bool passed_jet1 = jet3_is_finite(&spec->jet[1]);
  bool passed_jet2 = jet3_is_finite(&spec->jet[2]);
  bool passed_jet = passed_jet0 && passed_jet1 && passed_jet2;
  if (passed_jet0 || passed_jet1 || passed_jet2)
    assert(passed_jet);

  assert(passed_jet ^ passed_l); // exactly one of these

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
  u->x_minus_xb[0] = u->x_minus_xb[1] = u->x_minus_xb[2] = NAN;
  u->L = NAN;

  mesh3_s const *mesh = spec->eik ? eik3_get_mesh(spec->eik) : NULL;

  u->lhat = spec->lhat;

  /* Initialize x(hat) */
  if (passed_lhat) {
    assert(mesh);
    mesh3_copy_vert(mesh, u->lhat, u->x);
  } else { // passed_xhat
    dbl3_copy(spec->xhat, u->x);
  }

  memcpy(u->l, spec->l, sizeof(size_t[3]));

  if (passed_l)
    for (size_t i = 0; i < 3; ++i)
      mesh3_copy_vert(mesh, u->l[i], u->Xt[i]);
  else  // passed_x
    for (size_t i = 0; i < 3; ++i)
      dbl3_copy(spec->x[i], u->Xt[i]);

  dbl33_transposed(u->Xt, u->X);
  dbl33_mul(u->Xt, u->X, u->XtX);

  memcpy(u->state, spec->state, sizeof(state_e[3]));

  jet3 jet[3];
  for (size_t i = 0; i < 3; ++i)
    jet[i] = passed_jet ? spec->jet[i] : eik3_get_jet(spec->eik, spec->l[i]);

  dbl T[3], DT[3][3];
  for (int i = 0; i < 3; ++i) {
    T[i] = jet[i].f;
    memcpy(DT[i], &jet[i].fx, sizeof(dbl[3]));
  }
  bb32_init_from_3d_data(&u->T, T, &DT[0], u->Xt);

  u->num_shadow = 0;
  for (size_t i = 0; i < 3; ++i)
    u->num_shadow += u->state[i] == SHADOW;

  if (u->num_shadow == 1)
    init_split_s2(u, spec->eik);
  else if (u->num_shadow == 2)
    init_split_s1(u, spec->eik);
  else
    u->split = NULL;

  /* Check if we now have an update with all `SHADOW` indices. This
   * can happen when there's initially one or two `VALID` indices that
   * lie exactly on the shadow boundary, which get set to `SHADOW`
   * when `init_split_s1` or `init_split_s2` are called. */
  if (u->num_shadow == 3)
    return false;

  // Compute the surface normal for the plane spanned by (x1 - x0, x2
  // - x0), using DT[i] to determine its orientation. Return whether x
  // is on the right side of this plane.
  dbl n[3];
  dbl dx[2][3];
  dbl3_sub(u->Xt[1], u->Xt[0], dx[0]);
  dbl3_sub(u->Xt[2], u->Xt[0], dx[1]);
  dbl3_cross(dx[0], dx[1], n);
  dbl3_normalize(n);
  int sgn[3] = {
    signum(dbl3_dot(DT[0], n)),
    signum(dbl3_dot(DT[1], n)),
    signum(dbl3_dot(DT[2], n))
  };
  if (sgn[0] == -1)
    dbl3_negate(n);
  dbl x0_minus_x[3];
  dbl3_sub(u->Xt[0], u->x, x0_minus_x);
  dbl dot = -dbl3_dot(x0_minus_x, n) > 0;
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

/* Check if the point being updated lies in the plane spanned by by
 * x0, x1, and x2. If it does, the update is degenerate. */
bool utetra_is_degenerate(utetra_s const *u) {
  dbl const *x[4] = {u->x, u->Xt[0], u->Xt[1], u->Xt[2]};
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

  DT[0] = bb32_df(&cf->T, b, a1);
  DT[1] = bb32_df(&cf->T, b, a2);

  D2T[0][0] = bb32_d2f(&cf->T, b, a1, a1);
  D2T[1][0] = D2T[0][1] = bb32_d2f(&cf->T, b, a1, a2);
  D2T[1][1] = bb32_d2f(&cf->T, b, a2, a2);

  cf->f = cf->L + bb32_f(&cf->T, b);
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
    return dbl3_maxnorm(alpha[0]) < atol;

  if (jet[0].f > jet[1].f)
    return dbl3_maxnorm(alpha[1]) < atol;

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
  assert(dbl3_isfinite(xm));
  assert(dbl3_isfinite(xp));

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

  utetra_s const *u = split_utetra_select(utetra);

  dbl b[3];
  get_b(u, b);

  /* If we're dealing with a quad split (`n == 2`) and the update
   * point is in the interior of the quad, then we should return 4. */
  if (num_split(u) == 2 && b[0] > atol && b[1] > atol && b[2] <= atol)
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

/* Get the triangle of `u`. If `u` is split, this does *not* return
 * the triangle of the split updates (since that would not be
 * well-defined!). */
static tri3 get_tri(utetra_s const *u) {
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

static void get_alpha_for_active_inds_s1(utetra_s const *utetra, dbl alpha[3]) {
  assert(false); // TODO: take a careful look at this the first time it hits!

  utetra_s const *u = split_utetra_select(utetra);

  get_lag_mults(u, alpha);

  /* If the first Lagrange multiplier for `u` is nonzero, it indicates
   * that optimum is incident on the shadow boundary, which in the
   * interior of `utetra`. We set it to zero here to properly
   * communicate this to the caller of `utetra_get_active_inds`, which
   * properly refers to `utetra`. */
  alpha[0] = 0;

  /* When we set up `u` originally, we either swapped the position of
   * `l0` and `l1` or `l0` and `l2`, depending on the state of each
   * node, and so that ultimately `l0` is `VALID`, and `l1` and `l2`
   * are `SHADOW`. We can undo this permutation by checking whether
   * `utetra`'s `x0` equals `u`'s `x1` or `x2`.
   *
   * Note that when we compare below, we do exact, bitwise comparison,
   * since the vertices in `u` are direct copies of those in
   * `utetra`. */
  if (dbl3_equal(utetra->Xt[0], u->Xt[1]))
    SWAP(alpha[0], alpha[1]);
  else if (dbl3_equal(utetra->Xt[0], u->Xt[2]))
    SWAP(alpha[0], alpha[2]);
}

static void get_alpha_for_active_inds_s2(utetra_s const *utetra, dbl alpha[3]) {
  bool shadow_adj_split = utetra->split[0]->f > utetra->split[1]->f;
  utetra_s const *u = shadow_adj_split ? utetra->split[1] : utetra->split[0];

  get_lag_mults(u, alpha);

  /* Set Lagrange multiplier corresponding to internal edge to zero,
   * so that an interior point minimizer is treated as such. */
  alpha[2] = 0;

  /* Rearrange the Lagrange multipliers so that they match `utetra`'s
   * indices. We have to do this differently for each update. */
  if (shadow_adj_split)
    /* Only need to swap `l0` and `l1` for the split update adjacent
     * to the shadow zone... */
    SWAP(alpha[0], alpha[1]);
  else {
    assert(false); // TODO: take a careful look at this the first time it hits!

    /* ... but need to rotate for the other split. This is a bit
     * inscrutable, but to be clear: we want to set (0, 1, 2) <- (2,
     * 0, 1). (Drawing a picture will make it clear.) */
    ROTATE3(alpha[0], alpha[2], alpha[1]);
  }
}

size_t utetra_get_active_inds(utetra_s const *utetra, size_t l[3]) {
  assert(update_inds_are_set(utetra));

  dbl const atol = 1e-14;

  dbl alpha[3];
  switch (num_split(utetra)) {
  case 0: get_lag_mults(utetra, alpha); break;
  case 1: get_alpha_for_active_inds_s1(utetra, alpha); break;
  case 2: get_alpha_for_active_inds_s2(utetra, alpha); break;
  default: assert(false);
  }

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
    dbl x[3]; utetra_get_x(utetra, x);
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
