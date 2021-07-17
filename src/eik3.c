#include "eik3.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include "alist.h"
#include "array.h"
#include "bb.h"
#include "edge.h"
#include "heap.h"
#include "log.h"
#include "macros.h"
#include "mat.h"
#include "mesh3.h"
#include "slerp.h"
#include "utetra.h"
#include "util.h"
#include "utri.h"
#include "vec.h"

typedef struct bde_bc {
  size_t le[2];
  bb31 bb;
} bde_bc_s;

bde_bc_s make_bde_bc(size_t const le[2], bb31 const *bb) {
  assert(le[0] != le[1]);
  bde_bc_s bc = {.le = {le[0], le[1]}, .bb = *bb};
  if (le[0] > le[1]) {
    SWAP(bc.le[0], bc.le[1]);
    bb31_reverse(&bc.bb);
  }
  return bc;
}

// TODO:
//
// - There's a common pattern that gets reused here a lot:
//
//   1. Use a mesh iterator to enumerate the bases of a set of
//      `utetra` or `utri`
//   2. Do all these updates
//   3. Sort them
//   4. Commit an update based on what's going on with the first
//      couple of updates
//
//   I rewrote this from scratch each time I did it because I wasn't
//   sure how often I'd be doing it! But there you are. Ideally, we
//   should refactor these "updates sequences" into new module (e.g.,
//   utetras.h, utris.h), and call them from here. It will make the
//   code in this file much easier to understand and likely suss out a
//   few bugs.

/* A structure managing a jet marching method solving the eikonal
 * equation in 3D on an unstructured tetrahedron mesh.
 *
 * NOTE: this is just for s = 1 at the moment. Will extend this to
 * handle s != later. */
struct eik3 {
  mesh3_s *mesh;
  eik3_s const *orig; /* This field's originator */
  jet3 *jet;
  dbl33 *hess;
  state_e *state;
  int *pos;
  par3_s *par;
  heap_s *heap;

  /* In some cases, we'll skip old updates that might be useful at a
   * later stage. We keep track of them here. */
  array_s *old_updates;
  array_s *old_bd_utri; // old two-point boundary `utri`

  /* The number of times each vertex has had a boundary condition
   * set. This is used to quickly check where what kind of point we're
   * updating from. Each entry of `num_BCs` starts at `0` and is
   * incremented each time a boundary condition is added. For a point
   * source, there will be exactly one node with `num_BCs[l] ==
   * 1`. For a diffracting edge, each interior node of the BCs will
   * have `num_BCs[l] == 2`, the terminal nodes will have `num_BCs[l]
   * == 1`, and similarly for reflectors, where `num_BCs[l] <= 3` will
   * hold. */
  int8_t *num_BCs;

  /* Boundary conditions for a (diffracting) boundary edge, which are
   * just cubic polynomials defined over the edge. These are used to
   * perform two-point updates from diffracting edges initially when
   * solving edge diffraction problems. */
  array_s *bde_bc;

  /* Field of unit vectors where `t_in[l]` indicates the ray direction
   * of the ray leading into the point of reflection or
   * diffraction. This will be `NAN` for a point source problem. */
  dbl (*t_in)[3];

  /* Field of unit vectors where `t_out[l]` indicates the ray
   * direction at the beginning of the ray leading to node `l`.*/
  dbl (*t_out)[3];

  alist_s *rho1;

  /* Convergence tolerances. The parameter `h` is an estimate of the
   * fineness of the mesh. The variables `h2` and `h3` are convenience
   * variables containing `h^2` and `h^3`. */
  dbl h, h2, h3;

  dbl slerp_tol;

  /* Useful statistics for debugging */
  size_t num_accepted; /* number of nodes fixed by `eik3_step` */

  /* An array containing the order in which the individual nodes were
   * accepted. That is, `accepted[i] == l` means that `eik3_step()`
   * returned `l` when it was called for the `i`th time. */
  size_t *accepted;

  ftype_e ftype;
};

void eik3_alloc(eik3_s **eik) {
  *eik = malloc(sizeof(eik3_s));
}

void eik3_dealloc(eik3_s **eik) {
  free(*eik);
  *eik = NULL;
}

static dbl value(void *ptr, int l) {
  eik3_s *eik = (eik3_s *)ptr;
  assert(l >= 0);
  assert(l < (int)mesh3_nverts(eik->mesh));
  dbl T = eik->jet[l].f;
  return T;
}

static void setpos(void *ptr, int l, int pos) {
  eik3_s *eik = (eik3_s *)ptr;
  eik->pos[l] = pos;
}

void eik3_init(eik3_s *eik, mesh3_s *mesh, ftype_e ftype, eik3_s const *orig) {
  eik->mesh = mesh;
  eik->orig = orig;
  eik->ftype = ftype;

  size_t nverts = mesh3_nverts(mesh);

  eik->jet = malloc(nverts*sizeof(jet3));
  for (size_t l = 0; l < nverts; ++l) {
    eik->jet[l] = (jet3) {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN};
  }

  eik->hess = malloc(nverts*sizeof(dbl33));
  for (size_t l = 0; l < nverts; ++l)
    dbl33_nan(eik->hess[l]);

  eik->state = malloc(nverts*sizeof(state_e));
  for (size_t l = 0; l < nverts; ++l) {
    eik->state[l] = FAR;
  }

  eik->pos = malloc(nverts*sizeof(int));
  for (size_t l = 0; l < nverts; ++l) {
    eik->pos[l] = NO_INDEX;
  }

  eik->par = malloc(nverts*sizeof(par3_s));
  for (size_t l = 0; l < nverts; ++l) {
    par3_init_empty(&eik->par[l]);
  }

  /**
   * When we compute the initial heap capacity, we want to estimate
   * the number of nodes that could comprise the expanding numerical
   * front at any one time. We can't know this ahead of time, so we
   * set it to a constant multiple times (# nodes)^(1/d). In this
   * case, d=3. Even if this is an underestimate, well still reduce
   * the number of times the heap needs to be expanded at solve time.
   */
  int capacity = (int)3*cbrt(nverts);

  heap_alloc(&eik->heap);
  heap_init(eik->heap, capacity, value, setpos, (void *)eik);

  eik->num_accepted = 0;

  eik->accepted = malloc(nverts*sizeof(size_t));
  for (size_t i = 0; i < nverts; ++i)
    eik->accepted[i] = (size_t)NO_INDEX;

  array_alloc(&eik->old_updates);
  array_init(eik->old_updates, sizeof(utetra_s *), 16);

  array_alloc(&eik->old_bd_utri);
  array_init(eik->old_bd_utri, sizeof(utri_s *), 16);

  eik->num_BCs = calloc(nverts, sizeof(int8_t));

  array_alloc(&eik->bde_bc);
  array_init(eik->bde_bc, sizeof(bde_bc_s), ARRAY_DEFAULT_CAPACITY);

  if (eik->ftype != FTYPE_POINT_SOURCE) {
    eik->t_in = malloc(nverts*sizeof(dbl[3]));
    for (size_t l = 0; l < nverts; ++l)
      for (size_t i = 0; i < 3; ++i)
        eik->t_in[l][i] = NAN;
  }

  eik->t_out = malloc(nverts*sizeof(dbl[3]));
  for (size_t l = 0; l < nverts; ++l)
    for (size_t i = 0; i < 3; ++i)
      eik->t_out[l][i] = NAN;

  if (eik->ftype == FTYPE_EDGE_DIFFRACTION) {
    alist_alloc(&eik->rho1);
    alist_init(eik->rho1, sizeof(size_t), sizeof(dbl), 16);
  }

  eik->h = mesh3_get_min_edge_length(mesh);
  eik->h2 = eik->h*eik->h;
  eik->h3 = eik->h*eik->h*eik->h;

  eik->slerp_tol = eik->h3;
}

void eik3_deinit(eik3_s *eik) {
  free(eik->jet);
  eik->jet = NULL;

  free(eik->hess);
  eik->hess = NULL;

  free(eik->state);
  eik->state = NULL;

  free(eik->pos);
  eik->pos = NULL;

  free(eik->par);
  eik->par = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);

  utetra_s *utetra;
  for (size_t i = 0; i < array_size(eik->old_updates); ++i) {
    array_get(eik->old_updates, i, &utetra);
    utetra_deinit(utetra);
    utetra_dealloc(&utetra);
  }
  array_deinit(eik->old_updates);
  array_dealloc(&eik->old_updates);

  utri_s *utri;
  for (size_t i = 0; i < array_size(eik->old_bd_utri); ++i) {
    array_get(eik->old_bd_utri, i, &utri);
    utri_dealloc(&utri);
  }
  array_deinit(eik->old_bd_utri);
  array_dealloc(&eik->old_bd_utri);

  free(eik->num_BCs);
  eik->num_BCs = NULL;

  array_deinit(eik->bde_bc);
  array_dealloc(&eik->bde_bc);

  if (eik->ftype != FTYPE_POINT_SOURCE) {
    free(eik->t_in);
    eik->t_in = NULL;
  }

  free(eik->t_out);
  eik->t_out = NULL;

  if (eik->ftype == FTYPE_EDGE_DIFFRACTION) {
    alist_deinit(eik->rho1);
    alist_dealloc(&eik->rho1);
  }
}

static bool is_singular(eik3_s const *eik, size_t l) {
  bool singular_gradient = !dbl3_isfinite(&eik->jet[l].fx);
#if JMM_DEBUG
  assert(singular_gradient == !dbl33_isfinite(eik->hess[l]));
#endif
  return singular_gradient;
}

/* Check whether a point is a terminal point on a diffracting
 * edge. This is true if there is only one diffracting edge with a
 * full complement of BCs incident on `l`. */
static bool is_diff_edge_terminal_point(eik3_s const *eik, size_t l) {
  mesh3_s const *mesh = eik->mesh;

  size_t nbde = mesh3_get_num_inc_diff_edges(mesh, l);
  size_t (*le)[2] = malloc(nbde*sizeof(size_t[2]));
  mesh3_get_inc_diff_edges(mesh, l, le);

  size_t num_bde_bc = 0;

  for (size_t i = 0; i < nbde; ++i)
    num_bde_bc += eik3_has_bde_bc(eik, le[i]);

  free(le);

  return num_bde_bc == 1;
}

static bool can_update_from_point(eik3_s const *eik, size_t l) {
  return eik->state[l] == VALID && !eik3_is_point_source(eik, l);
}

static void prop_hess_from_pt_src(eik3_s *eik, size_t l, size_t l0) {
  dbl L = dbl3_dist(mesh3_get_vert_ptr(eik->mesh, l),
                    mesh3_get_vert_ptr(eik->mesh, l0));

  dbl eye[3][3]; dbl33_eye(eye);
  dbl op[3][3]; dbl3_outer(&eik->jet[l].fx, &eik->jet[l].fx, op);
  dbl33_sub(eye, op, eik->hess[l]);
  dbl33_dbl_div_inplace(eik->hess[l], L);

  assert(dbl33_isfinite(eik->hess[l]));
}

static void do_1pt_update(eik3_s *eik, size_t l, size_t l0,
                          bool require_commit) {
  // TODO: check if the updating ray is physical... should probably
  // just design a new ADT, `uline`... UGH...

  /* Compute the new jet */
  jet3 jet;
  dbl const *x = mesh3_get_vert_ptr(eik->mesh, l);
  dbl const *x0 = mesh3_get_vert_ptr(eik->mesh, l0);
  dbl3_sub(x, x0, &jet.fx);
  dbl L = dbl3_normalize(&jet.fx);
  jet.f = eik->jet[l0].f + L;

  bool should_commit = jet.f <= eik->jet[l].f;

  if (!require_commit && !should_commit)
    return;

  assert(should_commit);
  eik3_set_jet(eik, l, jet);

  /* Compute the Hessian of the eikonal at `l` */
  prop_hess_from_pt_src(eik, l, l0);

  par3_s par = {.l = {l0, NO_PARENT, NO_PARENT}, .b = {1, NAN, NAN}};
  eik3_set_par(eik, l, par);
}

/* Compute the Hessian (`hess`) of the eikonal at a point which was
 * updated by the triangle update `utri`, assuming that the base of
 * `utri` is a diffracting edge. The fact that the diffracting edge is
 * a caustic forces us to compute the Hessian directly from the
 * curvatures of the wavefront at the update point. */
static bool prop_hess_from_diff_edge(dbl const t[3], dbl L, dbl const e[3],
                                     dbl rho1, dbl hess[3][3]) {
  assert(dbl3_isfinite(t));
  assert(isfinite(L));
  assert(dbl3_isfinite(e));
  assert(isfinite(rho1));

  /* Compute the dot product between the ray vector and the edge
   * direction. We need it when we compute each of the curvatures
   * below */
  dbl t_dot_e = dbl3_dot(t, e);

  /* If `t` and `e` are colinear, we're dealing with a grazing ray. In
   * this case, we'll be unable to compute the Hessian correctly, so
   * we should bail now. */
  if (fabs(1 - fabs(t_dot_e)) < 1e-14)
    return false;

  /* Principal curvature and direction for the normal section
   * corresponding to the plane of diffraction (i.e., the plane
   * spanned by the edge and the ray updating \hat{x}. */
  dbl k1 = 1/rho1, q1[3];
  dbl3_saxpy(-t_dot_e, t, e, q1);

  /* Principal curvature and direction for the normal section
   * orthogonal to ray and `q1`. The wavefront in this section is a
   * circle about the diffracting edge that goes through \hat{x}. */
  dbl k2, q2[3];
  k2 = 1/(L*sqrt(1 - t_dot_e*t_dot_e));
  dbl3_cross(t, e, q2);

  /* Set the Hessian to be the sum of the outer product of each
   * principal direction, weighted by the principal curvatures */
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      hess[i][j] = k1*q1[i]*q1[j] + k2*q2[i]*q2[j];

  assert(dbl33_isfinite(hess));

  return true;
}

static void prop_hess_along_ray(dbl const t[3], dbl L,
                                dbl const hess0[3][3], dbl hess[3][3]) {
  /* The first row of `Q` will be the ray direction, and the second
   * and third will be unit vectors that span the orthogonal
   * complement of the ray direction, filling out the basis */
  dbl Q[3][3];

  /* The ray direction vector (i.e., `t_out`) */
  dbl3_copy(t, Q[0]);

  /* Get a random vector that's orthogonal to `t` */
  dbl3_get_rand_ortho(Q[0], Q[1]);

  /* Complete the basis */
  dbl3_cross(Q[0], Q[1], Q[2]);
  dbl3_normalize(Q[2]); // TODO: probably unnecessary

  /* Compute `P0` by restricting `hess0` to the orthogonal complement
   * of `t` (i.e., the normal plane of the ray). */
  dbl P0[2][2];
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      P0[i][j] = dbl3_dbl33_dbl3_dot(Q[i + 1], hess0, Q[j + 1]);

  /* Propagate `P0` along the ray to compute `P` */
  dbl eye[2][2]; dbl22_eye(eye);
  dbl tmp1[2][2]; dbl22_saxpy(L, P0, eye, tmp1); dbl22_invert(tmp1);
  dbl P[2][2]; dbl22_mul(P0, tmp1, P);

  /* Recover the Hessian at `xhat` from `P` */
  dbl tmp2[3][3]; dbl33_zero(tmp2); /* 3x3 matrix with lower right 2x2 subblock set to P */
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      tmp2[i + 1][j + 1] = P[i][j];
  dbl tmp3[3][3]; dbl33_mul(tmp2, Q, tmp3);
  dbl33_transposed(Q, tmp2);
  dbl33_mul(tmp2, tmp3, hess);
}

static void get_hess_cc(eik3_s const *eik, size_t const *l, dbl const *b,
                        size_t n, dbl hess[3][3]) {
  for (size_t i = 0; i < n; ++i)
    assert(dbl33_isfinite(eik->hess[l[i]]));

  dbl33_zero(hess);
  for (size_t k = 0; k < n; ++k)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        hess[i][j] += b[k]*eik->hess[l[k]][i][j];
}

static
bool get_rho1_cc(eik3_s const *eik, dbl const e[3], size_t const *l,
                 dbl const *b, size_t n, dbl *rho1_cc) {
  assert(n == 1 || n == 2);

  dbl rho1[n];
  for (size_t i = 0; i < n; ++i)
    rho1[i] = NAN;

  size_t num_set = 0;

  /* First, fill in the `rho1` values that are set as BCs
   * initially. This will only ever happen if we're solving an edge
   * diffraction problem. */
  if (eik->ftype == FTYPE_EDGE_DIFFRACTION) {
    for (size_t i = 0; i < n; ++i) {
      if (!alist_contains(eik->rho1, &l[i]))
        continue;
      assert(eik3_has_BCs(eik, l[i]));
      alist_get_by_key(eik->rho1, &l[i], &rho1[i]);
      num_set += isfinite(rho1[i]);
      assert(isfinite(rho1[i]));
    }
  }

  /* Fill in the remaining `rho1` values by computing the radius of
   * curvature using the Hessian at each edge endpoint. */
  if (num_set < n) {
    for (size_t i = 0; i < n; ++i) {
      if (isfinite(rho1[i]))
        continue;

      dbl const *t = &eik->jet[l[i]].fx;
      if (!dbl3_isfinite(t))
        continue;

      dbl q[3]; dbl3_saxpy(-dbl3_dot(t, e), t, e, q);

      if (dbl33_isfinite(eik->hess[l[i]])) {
        rho1[i] = 1/dbl3_wnormsq(eik->hess[l[i]], q);
      } else {
        assert(dbl33_isnan(eik->hess[l[i]]));
        rho1[i] = 0;
      }
    }
  }

  for (size_t i = 0; i < n; ++i)
    if (!isfinite(rho1[i]))
      return false;

  *rho1_cc = n == 1 ? rho1[0] : dbl2_dot(b, rho1);

  return *rho1_cc != 0 && isfinite(*rho1_cc);
}

static
void prop_hess_degenerate_1(eik3_s *eik, dbl const t[3], dbl L,
                            size_t l0, size_t lhat) {
  /* If there's only one active vertex and there's no Hessian at
   * that point, then we should have done a diffracting update,
   * which we need to figure out now. */

  size_t num_inc_diff_edges = mesh3_get_num_inc_diff_edges(eik->mesh, l0);
  assert(num_inc_diff_edges > 0);

  size_t (*le)[2] = malloc(num_inc_diff_edges*sizeof(size_t[2]));
  mesh3_get_inc_diff_edges(eik->mesh, l0, le);

  // TODO: this is hacky and needs to be fixed: we find the first
  // diffracting edge that has `rho1` set for each endpoint and
  // just use that. This isn't correct in the general case but
  // should be OK for now.
  size_t l[2] = {l0, (size_t)NO_INDEX};
  dbl e[3];
  for (size_t i = 0; i < num_inc_diff_edges; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      if (le[i][j] != l[0] && alist_contains(eik->rho1, &le[i][j])) {
        l[1] = le[i][j];
        mesh3_get_diff_edge_tangent(eik->mesh, l, e);
        break;
      }
    }
  }

  dbl rho1;
  assert(get_rho1_cc(eik, e, l, (dbl[2]) {1, 0}, 2, &rho1));
  prop_hess_from_diff_edge(t, L, e, rho1, eik->hess[lhat]);
}

static bool get_hess_from_utri(utri_s const *u, eik3_s const *eik, dbl33 hess) {
  dbl t[3];
  utri_get_t(u, t);

  dbl L = utri_get_L(u);

  mesh3_s const *mesh = eik->mesh;

  size_t l[2];
  utri_get_update_inds(u, l);
  assert(mesh3_is_diff_edge(mesh, l));

  dbl e[3];
  mesh3_get_diff_edge_tangent(mesh, l, e);

  dbl b = utri_get_b(u);

  dbl rho1;
  if (!get_rho1_cc(eik, e, l, (dbl[2]) {1 - b, b}, 2, &rho1))
    return false;

  prop_hess_from_diff_edge(t, L, e, rho1, hess);

  return true;
}

static
bool prop_hess_degenerate_N(eik3_s *eik, dbl const t[3], dbl L,
                            size_t const *l, size_t n, dbl const xb[3],
                            size_t lhat) {

  array_s *utri_arr;
  array_alloc(&utri_arr);
  array_init(utri_arr, sizeof(utri_s *), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < n; ++i) {
    if (!mesh3_bdv(eik->mesh, l[i]))
      continue;

    if (dbl33_isfinite(eik->hess[l[i]]))
      continue;

    size_t nde = mesh3_get_num_inc_diff_edges(eik->mesh, l[i]);
    if (nde == 0)
      continue;

    size_t (*le)[2] = malloc(nde*sizeof(size_t[2]));
    mesh3_get_inc_diff_edges(eik->mesh, l[i], le);

    for (size_t j = 0; j < nde; ++j) {
      utri_s *utri;
      utri_alloc(&utri);

      if (!eik3_has_BCs(eik, le[j][0]) ||
          !eik3_has_BCs(eik, le[j][1]))
        continue;

      if (mesh3_edge_contains_point(eik->mesh, le[j], xb))
        continue;

      utri_spec_s spec = utri_spec_from_eik_without_l(eik,xb,le[j][0],le[j][1]);
      if (utri_init(utri, &spec)) {
        utri_solve(utri);
        assert(utri_is_finite(utri));
        array_append(utri_arr, &utri);
      } else {
        utri_dealloc(&utri);
      }
    }

    free(le);
  }

  assert(!array_is_empty(utri_arr));

  array_sort(utri_arr, (compar_t)utri_cmp);

  dbl33 hess0;
  dbl33_nan(hess0);

  for (size_t i = 0; i < array_size(utri_arr); ++i) {
    utri_s *utri[2];
    array_get(utri_arr, i, &utri[0]);

    bool should_get_hess_from_utri = false;
    if (utri_has_interior_point_solution(utri[0]) &&
        utri_update_ray_is_physical(utri[0], eik)) {
      should_get_hess_from_utri = true;
    } else if (i + 1 < array_size(utri_arr)) {
      array_get(utri_arr, i + 1, &utri[1]);
      if (utris_yield_same_update(utri[0], utri[1]) &&
          utri_update_ray_is_physical(utri[0], eik)) {
        should_get_hess_from_utri = true;
      }
    }

    if (should_get_hess_from_utri &&
        get_hess_from_utri(utri[0], eik, hess0))
      break;
  }

  bool success = dbl33_isfinite(hess0);

  if (success)
    prop_hess_along_ray(t, L, hess0, eik->hess[lhat]);

  array_deinit(utri_arr);
  array_dealloc(&utri_arr);

  return success;
}

static
bool prop_hess(eik3_s *eik, size_t num_active, size_t const *l, dbl const *b,
               size_t lhat, dbl const t[3], dbl L, dbl const xb[3]) {
  /* Figure out which Hessians are finite and count them */
  size_t num_finite_hess = 0;
  bool finite_hess[num_active];
  for (size_t i = 0; i < num_active; ++i)
    num_finite_hess += finite_hess[i] = dbl33_isfinite(eik->hess[l[i]]);

  /* Check if the Hessians are consistent: */
  for (size_t i = 0; i < num_active; ++i) {
    if (finite_hess[i])
      continue;
#if JMM_DEBUG
    assert(eik3_has_BCs(eik, l[i]));
    /* for edge-diffracted waves, we should have set a value for rho1
     * as a BC if the Hessian isn't finite */
    if (eik->ftype == FTYPE_EDGE_DIFFRACTION)
      assert(alist_contains(eik->rho1, &l[i]));
#endif
    /* for a reflected wave, we aren't able to propagate the Hessian
     * using the available data, so we return now */
    if (eik->ftype == FTYPE_REFLECTION)
      return false;
  }

  bool success = false;

  if (num_finite_hess < num_active) {
    if (num_active == 1) {
      prop_hess_degenerate_1(eik, t, L, l[0], lhat);
      success = true;
    } else {
      success = prop_hess_degenerate_N(eik, t, L, l, num_active, xb, lhat);
    }
  } else {
    dbl33 hess0;
    get_hess_cc(eik, l, b, num_active, hess0);
    prop_hess_along_ray(t, L, hess0, eik->hess[lhat]);
    success = true;
  }

  return success;
}

static bool commit_tri_update(eik3_s *eik, size_t lhat, utri_s const *utri) {
  if (utri_get_value(utri) >= eik->jet[lhat].f)
    return false;

  jet3 jet;
  utri_get_jet(utri, &jet);
  eik3_set_jet(eik, lhat, jet);

  par3_s par = utri_get_par(utri);
  eik3_set_par(eik, lhat, par);

  dbl t[3];
  utri_get_t(utri, t);

  dbl L = utri_get_L(utri);

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  /* First, we want to try to compute the Hessian when the base of the
   * update is a diffracting edge.
   *
   * Before doing this, we make sure to check that the `utri` isn't
   * geometrically degenerate (i.e., x, x0, and x1 aren't
   * collinear). If this happens, the result doesn't will be
   * singular. This situation will arise when we try to march away
   * from a diffracting edge emanating from the diffracting edge from
   * which this field emanated originally. In these cases, we can just
   * use `prop_hess` below, since the Hessian will be defined and
   * usable.  */
  if (utri_inc_on_diff_edge(utri, eik->mesh) &&
      !utri_is_degenerate(utri)) {
    /* If the `utri` traced a "terminal ray" from the end of a
     * diffracting edge, it behaves the same as a point source. */
    if (num_active == 1 &&
        is_singular(eik, l[0]) &&
        utri_emits_terminal_ray(utri, eik)) {
      prop_hess_from_pt_src(eik, lhat, l[0]);
      return true;
    }

    size_t l_base[2]; utri_get_update_inds(utri, l_base);

    /* Get the diffracting edge's unit tangent vector */
    dbl e[3]; mesh3_get_diff_edge_tangent(eik->mesh, l_base, e);

    /* Get the radius of curvature in the plane of diffraction of the
     * incident wavefront at the diffracting edge. Then try to
     * propagate the Hessian from the diffracting edge using rho1. */
    dbl rho1;
    if (get_rho1_cc(eik, e, l, b, num_active, &rho1) &&
        prop_hess_from_diff_edge(t, L, e, rho1, eik->hess[lhat]))
      return true;

    /* If that fails, try to use finite differences to approximate the
     * Hessian (this involves solving several more `utri`s */
    if (utri_approx_hess(utri, eik->h2, eik->hess[lhat]))
      return true;

    /* We failed to update the Hessian---don't accept the update */
    // TODO: we should always be able to do this...
    return false;
  }

  dbl xb[3];
  utri_get_xb(utri, xb);

  prop_hess(eik, num_active, l, b, lhat, t, L, xb);

  return true;
}

static void do_2pt_bd_updates(eik3_s *eik, size_t l, size_t l0) {
  size_t nve = mesh3_nve(eik->mesh, l);
  size_t (*ve)[2] = malloc(nve*sizeof(size_t[2]));
  mesh3_ve(eik->mesh, l, ve);

  utri_s **utri = calloc(nve, sizeof(size_t *));

  size_t l1;
  size_t nup = 0;
  for (size_t i = 0; i < nve; ++i) {
    utri[i] = NULL;
    if (ve[i][0] == l0 || ve[i][1] == l0) {
      l1 = ve[i][0] == l0 ? ve[i][1] : ve[i][0];
      size_t lf[3] = {l, l0, l1};
      if (can_update_from_point(eik, l1) &&
          !eik3_is_point_source(eik, l1) &&
          mesh3_is_bdf(eik->mesh, lf)) {
        utri_alloc(&utri[i]);
        utri_spec_s spec = utri_spec_from_eik(eik, l, l0, l1);
        ++nup;
        if (utri_init(utri[i], &spec))
          utri_solve(utri[i]);
      }
    }
  }

  qsort(utri, nve, sizeof(utri_s *), (compar_t)utri_cmp);

  for (size_t i = 0; i < nup; ++i)
    assert(utri[i] != NULL);
  for (size_t i = nup; i < nve; ++i)
    assert(utri[i] == NULL);

  /* Go through old boundary triangle updates and append any that have
   * the same update index and share an edge with the current updates
   * in `utri`. We will have solved these already, so don't need to
   * resolve them. */
  size_t copied_utri = 0;
  for (size_t i = 0; i < nve; ++i) {
    if (utri[i] == NULL)
      continue;
    utri_s *u_old;
    for (size_t j = array_size(eik->old_bd_utri); j > 0; --j) {
      array_get(eik->old_bd_utri, j - 1, &u_old);
      if (utri_get_l(u_old) != l ||
          !utri_opt_inc_on_other_utri(u_old, utri[i]))
        continue;
      array_delete(eik->old_bd_utri, j - 1);
      if (nup + copied_utri + 1 > nve)
        utri = realloc(utri, (nup + copied_utri + 1)*sizeof(utri_s *));
      utri[nup + copied_utri++] = u_old;
    }
  }
  nup += copied_utri;

  qsort(utri, nup, sizeof(utri_s *), (compar_t)utri_cmp);

  // Keep track of which utri to free. If we copy any over to
  // `old_bd_utri`, we want to make sure we don't accidentally free
  // them.
  bool *free_utri = malloc(nup*sizeof(bool));
  for (size_t i = 0; i < nup; ++i)
    free_utri[i] = true;

  // Try to commit a triangle update
  //
  // TODO: as always, this is a bit complicated. Would be nice to
  // simplify this or at least factor it out somewhere else.
  for (size_t i = 0; i < nup; ++i) {
    if (!utri_is_finite(utri[i]))
      break;
    if (utri_has_interior_point_solution(utri[i])) {
      if (utri_update_ray_is_physical(utri[i], eik) &&
          commit_tri_update(eik, l, utri[i])) {
        break;
      }
    } else if (i + 1 < nup &&
               utris_yield_same_update(utri[i], utri[i + 1]) &&
               utri_update_ray_is_physical(utri[i], eik) &&
               commit_tri_update(eik, l, utri[i])) {
      break;
    } else {
      array_append(eik->old_bd_utri, &utri[i]);
      free_utri[i] = false;
    }
  }

  for (size_t i = 0; i < nup; ++i) {
    if (utri[i] == NULL || !free_utri[i]) continue;
    utri_dealloc(&utri[i]);
  }

  free(free_utri);
  free(utri);
  free(ve);
}

static bool can_update_from_face(eik3_s const *eik,
                                 size_t l0, size_t l1, size_t l2) {
  return can_update_from_point(eik, l0) &&
    can_update_from_point(eik, l1) &&
    can_update_from_point(eik, l2);
}

static
size_t get_update_fan(eik3_s const *eik, size_t l0, size_t **l1, size_t **l2) {
  size_t nvc = mesh3_nvc(eik->mesh, l0);
  size_t *vc = malloc(nvc*sizeof(size_t));
  mesh3_vc(eik->mesh, l0, vc);

  array_s *le_arr;
  array_alloc(&le_arr);
  array_init(le_arr, sizeof(size_t[2]), ARRAY_DEFAULT_CAPACITY);

  bool is_valid[4];
  size_t cv[4], le[2], num_valid;
  for (size_t i = 0; i < nvc; ++i) {
    mesh3_cv(eik->mesh, vc[i], cv);

    num_valid = 0;
    for (size_t j = 0; j < 4; ++j)
      num_valid += is_valid[j] = eik->state[cv[j]] == VALID;

    if (num_valid != 3)
      continue;

    size_t k = 0;
    for (size_t j = 0; j < 4; ++j)
      if (cv[j] != l0 && is_valid[j])
        le[k++] = cv[j];
    assert(k == 2);

    SORT2(le[0], le[1]);

    if (!can_update_from_face(eik, l0, le[0], le[1]))
      continue;

    if (array_contains(le_arr, &le))
      continue;

    array_append(le_arr, &le);
  }

  size_t num_updates = array_size(le_arr);

  *l1 = malloc(num_updates*sizeof(size_t));
  *l2 = malloc(num_updates*sizeof(size_t));

  for (size_t i = 0; i < array_size(le_arr); ++i) {
    array_get(le_arr, i, &le);
    (*l1)[i] = le[0];
    (*l2)[i] = le[1];
  }

  array_deinit(le_arr);
  array_dealloc(&le_arr);

  free(vc);

  return num_updates;
}

static
void get_valid_incident_diff_edges(eik3_s const *eik, size_t l0,
                                   array_s *l1) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  int nvv = mesh3_nvv(mesh, l0);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(mesh, l0, vv);

  size_t le[2] = {[0] = l0};
  for (int i = 0; i < nvv; ++i) {
    le[1] = vv[i];

    if (!mesh3_is_diff_edge(mesh, le))
      continue;

    if (!eik3_is_valid(eik, le[1]))
      continue;

    if (can_update_from_point(eik, le[1]) || eik3_has_bde_bc(eik, le))
      array_append(l1, &le[1]);
  }

  free(vv);
}

static bool commit_tetra_update(eik3_s *eik, size_t lhat, utetra_s const *utetra) {
  if (utetra_get_value(utetra) >= eik->jet[lhat].f)
    return false;

  jet3 jet;
  utetra_get_jet(utetra, &jet);
  eik3_set_jet(eik, lhat, jet);

  par3_s par = utetra_get_parent(utetra);
  eik3_set_par(eik, lhat, par);

  dbl t[3];
  utetra_get_t(utetra, t);

  dbl L = utetra_get_L(utetra);

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  /* If we do a tetrahedron update from the end of a diffracting edge,
   * we might need to treat that vertex as a point source and
   * propagate the Hessian accordingly. */
  if (eik->ftype == FTYPE_EDGE_DIFFRACTION &&
      num_active == 1 &&
      mesh3_vert_incident_on_diff_edge(eik->mesh, l[0]) &&
      utetra_emits_terminal_ray(utetra, eik)) {
    prop_hess_from_pt_src(eik, lhat, l[0]);
    return true;
  }

  dbl xb[3];
  utetra_get_x(utetra, xb);

  if (prop_hess(eik, num_active, l, b, lhat, t, L, xb))
    return true;

  if (utetra_approx_hess(utetra, eik->h2, eik->hess[lhat]))
    return true;

  return false;
}

static void update(eik3_s *eik, size_t l, size_t l0) {
  /**
   * The first thing we do is check if we're trying to update from a
   * point source. In this case, we just solve the two-point BVP to
   * high accuracy and return immediately.
   */
  if (eik->ftype == FTYPE_POINT_SOURCE &&
      eik3_is_point_source(eik, l0)) {
    do_1pt_update(eik, l, l0, true);
    return;
  }

  // TODO: comment me
  if (eik->ftype == FTYPE_EDGE_DIFFRACTION &&
      eik3_has_BCs(eik, l0) &&
      is_singular(eik, l0) &&
      is_diff_edge_terminal_point(eik, l0)) {
    do_1pt_update(eik, l, l0, false);
    return;
  }

  // Next, if `l` is a boundary point, we want to do any two-point
  // updates that are immersed in the boundary. (These are "creeping
  // rays", which can be physical.)
  if (mesh3_bdv(eik->mesh, l) && mesh3_bdv(eik->mesh, l0))
    do_2pt_bd_updates(eik, l, l0);

  /**
   * First, find the "update triangle fan"
   */
  size_t *l1, *l2;
  size_t num_utetra = get_update_fan(eik, l0, &l1, &l2);
  if (num_utetra == 0)
    return;

  /**
   * Before doing tetrahedron updates, we want to check if there are
   * any diffracting edges updates that aren't adjacent to `l0`. These
   * won't be covered by `do_all_diff_edge_updates_and_adjust` in
   * `update_neighbors`.
   */

  // First, check which of the l1's and l2's are adjacent to l0

  // TODO: the way we're checking for adjacent updates here is pretty
  // inefficient, but not sure if we can do better...

  bool *l_l1_adj = malloc(num_utetra*sizeof(bool));
  for (size_t i = 0; i < num_utetra; ++i)
    l_l1_adj[i] = mesh3_is_edge(eik->mesh, (size_t[2]) {l, l1[i]});

  bool *l_l2_adj = malloc(num_utetra*sizeof(bool));
  for (size_t i = 0; i < num_utetra; ++i)
    l_l2_adj[i] = mesh3_is_edge(eik->mesh, (size_t[2]) {l, l2[i]});

  // Count and mark the non-adjacent edges are diffracting edges

  size_t num_diff_edges = 0;
  bool *is_diff_edge = calloc(num_utetra, sizeof(bool));
  for (size_t i = 0; i < num_utetra; ++i) {
    if (l_l1_adj[i] || l_l2_adj[i]) {
      is_diff_edge[i] = false;
      continue;
    } else {
      size_t e[2] = {l1[i], l2[i]};
      is_diff_edge[i] = mesh3_is_diff_edge(eik->mesh, e);
      num_diff_edges += is_diff_edge[i];
    }
  }

  array_s *diff_utri;
  array_alloc(&diff_utri);
  array_init(diff_utri, sizeof(utri_s *), 2*num_diff_edges);

  utri_s *u;
  utri_spec_s spec;
  for (size_t i = 0; i < num_utetra; ++i) {
    if (!is_diff_edge[i]) continue;

    utri_alloc(&u);

    array_append(diff_utri, &u);

    spec = utri_spec_from_eik(eik, l, l1[i], l2[i]);
    spec.orig_index = i;
    if (utri_init(u, &spec))
      utri_solve(u);

    if (utri_has_interior_point_solution(u))
      continue;

    // TODO: the following section cries out for refactoring...

    size_t l_active = utri_get_active_ind(u);
    assert(l_active != (size_t)NO_INDEX);

    size_t l_inactive = utri_get_inactive_ind(u);
    assert(l_inactive != (size_t)NO_INDEX);

    size_t num_inc = mesh3_get_num_inc_diff_edges(eik->mesh, l_active);
    size_t (*le)[2] = malloc(num_inc*sizeof(size_t[2]));
    mesh3_get_inc_diff_edges(eik->mesh, l_active, le);

    for (size_t k = 0; k < num_inc; ++k) {
      /* Skip the current edge */
      if ((le[k][0] == l_active ^ le[k][0] == l_inactive) &&
          (le[k][1] == l_active ^ le[k][1] == l_inactive))
        continue;

      /* Check if this edge is already in `diff_utri`. */
      bool already_found = false;
      for (size_t m = 0; m < array_size(diff_utri); ++m) {
        utri_s *u_;
        array_get(diff_utri, m, &u_);
        if (utri_contains_update_ind(u_, le[k][0]) &&
            utri_contains_update_ind(u_, le[k][1])) {
          already_found = true;
          break;
        }
      }
      if (already_found)
        continue;

      /* Add the edge! */
      utri_alloc(&u);
      array_append(diff_utri, &u);
      spec = utri_spec_from_eik(eik, l, le[k][0], le[k][1]);
      if (utri_init(u, &spec))
        utri_solve(u);
    }

    free(le);
  }

  array_sort(diff_utri, (compar_t)utri_cmp);

  bool *updated_from_diff_edge = calloc(num_utetra, sizeof(bool));

  for (size_t i = 0, j; i < array_size(diff_utri); ++i) {
    utri_s *u, *u_;
    array_get(diff_utri, i, &u);
    if (!utri_is_finite(u))
      continue;
    if (utri_has_interior_point_solution(u)) {
      if (utri_update_ray_is_physical(u, eik) &&
          commit_tri_update(eik, l, u)) {
        if (utri_has_orig_index(u)) {
          j = utri_get_orig_index(u);
          updated_from_diff_edge[j] = true;
        }
        break;
      }
    } else if (i + 1 < array_size(diff_utri)) {
      array_get(diff_utri, i + 1, &u_);
      if (utris_yield_same_update(u, u_) &&
          utri_update_ray_is_physical(u, eik) &&
          commit_tri_update(eik, l, u)) {
        if (utri_has_orig_index(u)) {
          j = utri_get_orig_index(u);
          updated_from_diff_edge[j] = true;
        }
        if (utri_has_orig_index(u_)) {
          j = utri_get_orig_index(u_);
          updated_from_diff_edge[j] = true;
        }
      }
    }
  }

  for (size_t i = 0; i < array_size(diff_utri); ++i) {
    utri_s *u;
    array_get(diff_utri, i, &u);
    utri_dealloc(&u);
  }

  array_deinit(diff_utri);
  array_dealloc(&diff_utri);

  free(is_diff_edge);

  free(l_l1_adj);
  free(l_l2_adj);

  /**
   * Now we move on to doing tetrahedron updates
   */

  // Allocate tetrahedron updates
  utetra_s **utetra = malloc(num_utetra*sizeof(utetra_s *));

  // Do each tetrahedron update and sort
  for (size_t i = 0; i < num_utetra; ++i) {
    utetra_alloc(&utetra[i]);

    utetra_spec_s spec = utetra_spec_from_eik_and_inds(eik, l, l0, l1[i], l2[i]);
    if (!utetra_init(utetra[i], &spec))
      continue;

    // This is a gross hack. What we do here is prioritize a
    // diffracting edge that's incident on this tetrahedron update. It
    // might yield a somewhat higher value, but when we're close to a
    // diffracting edge, it's important to correct the ray to ensure
    // that it emits from the diffracting edge. So, we skip the
    // tetrahedron update here.
    if (updated_from_diff_edge[i])
      continue;

    // TODO: move into utetra_init?
    if (utetra_is_degenerate(utetra[i]))
      continue;

    utetra_solve(utetra[i], NULL);
  }

  free(updated_from_diff_edge);

  // Go through old updates and append any that have the same update
  // index (`l`) and share an edge with the updates currently in
  // `utetra`. Note: we will have already solved these updates! No
  // need to redo them.
  //
  // TODO: this is a bit of a mess :-(
  size_t copied_utetra = 0;
  for (size_t i = 0; i < num_utetra; ++i) {
    utetra_s *old_utetra;
    for (size_t j = array_size(eik->old_updates); j > 0; --j) {
      array_get(eik->old_updates, j - 1, &old_utetra);
      // TODO: should filter out the entries of `old_updates` that
      // have `lhat == l` beforehand so that we don't have to do this
      // check for each element of `utetra[]`... wasteful
      if (utetra_get_l(old_utetra) != l)
        continue;
      if (!utetra_opt_inc_on_other_utetra(old_utetra, utetra[i]))
        continue;
      // TODO: ACTUALLY, need to check if `old_utetra`'s optimum is
      // incident on `utetra[i]`'s base!
      array_delete(eik->old_updates, j - 1);
      utetra = realloc(
        utetra, (num_utetra + copied_utetra + 1)*sizeof(utetra_s *));
      utetra[num_utetra + copied_utetra++] = old_utetra;
    }
  }
  num_utetra += copied_utetra;

  // Sort the updates by their eikonal value
  qsort(utetra, num_utetra, sizeof(utetra_s *), (compar_t)utetra_cmp);

  // Keep track of which updates to free. If we copy any updates over
  // to `old_updates`, we want to make sure we don't accidentally free
  // them.
  bool *should_free_utetra = malloc(num_utetra*sizeof(bool));
  for (size_t i = 0; i < num_utetra; ++i)
    should_free_utetra[i] = true;

  // See if we can commit a tetrahedron update
  for (size_t i = 0; i < num_utetra; ++i) {
    if (!isfinite(utetra_get_value(utetra[i])))
      break;
    if (utetra_has_interior_point_solution(utetra[i]) ||
        utetra_emits_terminal_ray(utetra[i], eik)) {
      if (utetra_update_ray_is_physical(utetra[i], eik) &&
          commit_tetra_update(eik, l, utetra[i]))
        break;
    } else {
      size_t num_int = utetra_get_num_interior_coefs(utetra[i]);
      assert(num_int == 1 || num_int == 2);
      size_t num_adj = 4 - num_int;
      if (i + num_adj <= num_utetra &&
          utetras_yield_same_update((utetra_s const **)&utetra[i], num_adj) &&
          utetra_update_ray_is_physical(utetra[i], eik) &&
          commit_tetra_update(eik, l, utetra[i])) {
        break;
      } else {
        array_append(eik->old_updates, &utetra[i]);
        should_free_utetra[i] = false;
      }
    }
  }

  for (size_t i = 0; i < num_utetra; ++i)
    if (should_free_utetra[i]) {
      utetra_deinit(utetra[i]);
      utetra_dealloc(&utetra[i]);
    }
  free(utetra);

  free(should_free_utetra);

  free(l1);
  free(l2);
}

static void adjust(eik3_s *eik, size_t l) {
  assert(eik->state[l] == TRIAL);
  assert(l < mesh3_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l]);
}

size_t eik3_peek(eik3_s const *eik) {
  return heap_front(eik->heap);
}

void do_diff_edge_updates_and_adjust(eik3_s *eik, size_t l0, size_t l1,
                                     size_t *l0_nb, int l0_nnb) {
#if JMM_DEBUG
  if (!eik3_has_bde_bc(eik, (size_t[2]) {l0, l1})) {
    assert(!eik3_is_point_source(eik, l0));
    assert(!eik3_is_point_source(eik, l1));
  }
#endif

  int l1_nnb = mesh3_nvv(eik->mesh, l1);
  size_t *l1_nb = malloc(l1_nnb*sizeof(size_t));
  mesh3_vv(eik->mesh, l1, l1_nb);

  array_s *nb;
  array_alloc(&nb);
  array_init(nb, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  // TODO: should really do this work in the calling function, passing
  // in a copy of the array each time...
  for (int i = 0; i < l0_nnb; ++i)
    if (eik3_is_trial(eik, l0_nb[i]))
      array_append(nb, &l0_nb[i]);

  for (int i = 0; i < l1_nnb; ++i)
    if (eik3_is_trial(eik, l1_nb[i]) && !array_contains(nb, &l1_nb[i]))
      array_append(nb, &l1_nb[i]);

  free(l1_nb);

  int nnb = array_size(nb);

  utri_s *utri;
  utri_alloc(&utri);

  size_t l;

  // TODO: a problem with what I'm doing here: may need to do adjacent
  // diffracting edge updates. This could be a problem if I have
  // curved obstacle, or if a ray goes around the corner of an
  // obstacle. We'll see how far we can get with this for now...

  /* Do a triangle update for each neighbor of the diff. edge */
  for (int i = 0; i < nnb; ++i) {
    array_get(nb, i, &l);

    if (l == l0 || l == l1)
      continue;

    utri_spec_s spec = utri_spec_from_eik(eik, l, l0, l1);

    if (!utri_init(utri, &spec))
      continue;

    if (utri_is_degenerate(utri))
      continue;

    utri_solve(utri);

    /* If we're computing a reflection, we need to skip these kinds of
     * updates. We need to give the reflected wave a chance to
     * propagate to nearby points. The edge-diffracted triangle
     * updates can lead to artificially small arrival times. */
    if (eik->ftype == FTYPE_REFLECTION &&
        utri_inc_on_refl_BCs(utri, eik) &&
        !utri_accept_refl_BCs_update(utri, eik))
      continue;

    if ((utri_has_interior_point_solution(utri) ||
         utri_emits_terminal_ray(utri, eik)) &&
        utri_update_ray_is_physical(utri, eik) &&
        commit_tri_update(eik, l, utri))
      adjust(eik, l);
  }

  utri_dealloc(&utri);

  array_deinit(nb);
  array_dealloc(&nb);
}

void do_all_diff_edge_updates_and_adjust(eik3_s *eik, size_t l0, size_t *l0_nb,
                                         int l0_nnb) {
  array_s *l1;
  array_alloc(&l1);
  array_init(l1, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  get_valid_incident_diff_edges(eik, l0, l1);

  for (size_t i = 0; i < array_size(l1); ++i)
    do_diff_edge_updates_and_adjust(
      eik, l0, *(size_t *)array_get_ptr(l1, i), l0_nb, l0_nnb);

  array_deinit(l1);
  array_dealloc(&l1);
}

void update_neighbors(eik3_s *eik, size_t l0) {
  state_e l0_state = eik->state[l0];
  assert(l0_state == VALID);

  size_t l; // Node l is a neighbor of node l0

  // Get i0's neighboring nodes.
  int nnb = mesh3_nvv(eik->mesh, l0);
  size_t *nb = malloc(nnb*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, nb);

  if (l0_state == VALID) {
    // Set FAR nodes to TRIAL and insert them into the heap.
    for (int i = 0; i < nnb; ++i) {
      if (eik->state[l = nb[i]] == FAR) {
        eik->state[l] = TRIAL;
        heap_insert(eik->heap, l);
      }
    }
  }

  // Find newly VALID diffracting edges (l0, l1) and update all
  // neighboring nodes---these are nodes neighboring *both* l0 and
  // l1. If we don't look at all of these neighbors, we could miss
  // important edge updates.
  do_all_diff_edge_updates_and_adjust(eik, l0, nb, nnb);

  // Update neighboring nodes.
  for (int i = 0; i < nnb; ++i) {
    if (eik->state[l = nb[i]] == TRIAL) {
      update(eik, l, l0);
      adjust(eik, l); // TODO: we should avoid calling adjust
                      // repeatedly here and above. Instead we should
                      // use a flag to track which nodes actually had
                      // their values change and then adjust before
                      // finally returning from this function.
    }
  }

  free(nb);
}

static
void get_t_in_from_orig_a1(eik3_s const *eik, size_t l0, dbl3 t_in) {
  assert(!dbl3_isfinite(eik->t_in[l0]));

  mesh3_s const *mesh = eik->mesh;

  /* Find all of the diffracting edges that are incident on `l0` */
  size_t nbde = mesh3_get_num_inc_diff_edges(mesh, l0);
  size_t (*le)[2] = malloc(nbde*sizeof(size_t[2]));
  mesh3_get_inc_diff_edges(mesh, l0, le);

  /* Find the specific diff. edge for which the originating field had
   * BCs specified (and hence, is causing the missing `t_in` data that
   * we're imputing now)... */
  size_t i = 0;
  for (i = 0; i < nbde; ++i)
    if (eik3_has_bde_bc(eik->orig, le[i]))
      break;

  size_t le_orig[2];
  memcpy(le_orig, le, sizeof(size_t[2]));

#if JMM_DEBUG
  /* ... and verify that there's only one such diffracting edge */
  for (size_t j = i + 1; j < nbde; ++j)
    assert(!eik3_has_bde_bc(eik->orig, le[j]));
#endif

  for (i = 0; i < nbde; ++i)
    if (eik3_has_bde_bc(eik, le[i]))
      break;

  assert(le[i][0] == l0 ^ le[i][1] == l0);
  size_t l1 = le[i][0] == l0 ? le[i][1] : le[i][0];

#if JMM_DEBUG
  /* ... and verify that there's only one such diffracting edge */
  for (size_t j = i + 1; j < nbde; ++j)
    assert(!eik3_has_bde_bc(eik, le[j]));
#endif

  dbl const *x0 = mesh3_get_vert_ptr(mesh, l0);
  dbl const *x1 = mesh3_get_vert_ptr(mesh, l1);
  dbl const *xe = mesh3_get_vert_ptr(
    mesh, le_orig[0] == l0 ? le_orig[1] : le_orig[0]);

  dbl3 e_orig;
  dbl3_sub(xe, x0, e_orig);
  dbl3_normalize(e_orig);

  dbl3 e;
  dbl3_sub(x1, x0, e);
  dbl3_normalize(e);

  dbl const *t_in_orig = eik->orig->t_in[l0];
  assert(dbl3_isfinite(t_in_orig));

  dbl3 t_in_orig_proj;
  dbl3_saxpy(-dbl3_dot(t_in_orig, e_orig), e_orig, t_in_orig, t_in_orig_proj);
  dbl3_normalize(t_in_orig_proj);

  dbl theta = acos(clamp(dbl3_dot(t_in_orig_proj, e), -1, 1));
  dbl33 rot;
  dbl33_make_axis_angle_rotation_matrix(e_orig, -theta, rot);

  dbl33_dbl3_mul(rot, t_in_orig, t_in);
}

static
void get_t_in_from_orig_a2(eik3_s const *eik, size_t const *l, dbl const *b,
                           bool const *finite, dbl3 t_in) {
  size_t const num_active = 2;

  size_t l_nan = (size_t)NO_INDEX;
  for (size_t i = 0; i < num_active; ++i)
    if (!finite[i])
      l_nan = l[i];

  mesh3_s const *mesh = eik->mesh;

  dbl xb[3];
  dbl3_zero(xb);
  for (size_t i = 0; i < num_active; ++i) {
    dbl const *x = mesh3_get_vert_ptr(mesh, l[i]);
    for (size_t j = 0; j < 3; ++j)
      xb[j] += b[i]*x[j];
  }

  size_t nbde = mesh3_get_num_inc_diff_edges(mesh, l_nan);
  size_t (*le)[2] = malloc(nbde*sizeof(size_t[2]));
  mesh3_get_inc_diff_edges(mesh, l_nan, le);

  size_t i;
  for (i = 0; i < nbde; ++i)
    if (eik3_has_bde_bc(eik->orig, le[i]))
      break;

#if JMM_DEBUG
  for (size_t j = i + 1; j < nbde; ++j)
    assert(!eik3_has_bde_bc(eik->orig, le[j]));
#endif

  size_t l0 = le[i][0], l1 = le[i][1];
  utri_spec_s spec = utri_spec_from_eik_without_l(eik->orig, xb, l0, l1);

  utri_s *u;
  utri_alloc(&u);
  utri_init(u, &spec);
  utri_solve(u);
  utri_get_t(u, t_in);
  utri_dealloc(&u);

  free(le);
}

static
void get_t_in_from_orig(eik3_s const *eik, size_t const *l, dbl const *b,
                        bool const *finite, size_t num_active, dbl3 t_in) {
#if JMM_DEBUG
  assert(eik->orig != NULL);
  assert(1 <= num_active && num_active <= 3);
  size_t num_finite = 0;
  for (size_t i = 0; i < num_active; ++i)
    num_finite += finite[i];
  assert(num_finite == num_active - 1);
#endif
  if (num_active == 1)
    get_t_in_from_orig_a1(eik, l[0], t_in);
  if (num_active == 2)
    get_t_in_from_orig_a2(eik, l, b, finite, t_in);
  if (num_active == 3)
    assert(false); // TODO: not implemented yet
}

static void compute_t_in(eik3_s *eik, size_t l0) {
  assert(eik->ftype != FTYPE_POINT_SOURCE);
  assert(eik->state[l0] == VALID);

  par3_s par = eik3_get_par(eik, l0);

  /* If the node has no parents, it should have had some BCs
   * supplied. Verify this before returning. */
  if (par3_is_empty(&par)) {
    assert(eik->num_BCs[l0] > 0);
    return;
  }

  size_t num_active = par3_num_active(&par);
  assert(num_active > 0);

  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  dbl t_in[num_active][3];
  for (size_t i = 0; i < num_active; ++i)
    dbl3_copy(eik->t_in[l[i]], t_in[i]);

  bool finite[num_active];
  size_t num_finite = 0;
  for (size_t i = 0; i < num_active; ++i)
    num_finite += finite[i] = dbl3_isfinite(t_in[i]);

  for (size_t i = 0; i < num_active; ++i) {
    if (finite[i]) continue;
    assert(eik->num_BCs[l[i]] > 0);
    dbl3_copy(&eik->jet[l[i]].fx, t_in[i]);

    /* Note: this may fail! After running this loop, we may still have
     * num_finite < num_active. */
    num_finite += finite[i] = dbl3_isfinite(t_in[i]);
  }

  if (num_finite < num_active) {
    get_t_in_from_orig(eik, l, b, finite, num_active, eik->t_in[l0]);
    goto coda;
  }

  assert(num_finite == num_active);
  assert(1 <= num_active && num_active <= 3);

  if (num_active == 1)
    dbl3_copy(t_in[0], eik->t_in[l0]);
  if (num_active == 2)
    slerp2(t_in[0], t_in[1], b, eik->t_in[l0]);
  if (num_active == 3)
    // slerp3(t_in, b, eik->t_in[l0], eik->slerp_tol);
    nlerp3(t_in, b, eik->t_in[l0]);

coda:
  assert(dbl3_isfinite(eik->t_in[l0]));
}

static void compute_t_out(eik3_s *eik, size_t l0) {
  assert(eik->state[l0] == VALID);

  par3_s par = eik3_get_par(eik, l0);

  /* If the node has no parents, then it must have boundary
   * conditions. We stipulate that `t_out` is undefined (NAN) at
   * points with boundary conditions (although we could specify a
   * well-defined field of `t_out` vectors in the case of a
   * reflection). We do this for two reasons:
   *
   * 1) to ensure that the nodes neighboring the nodes with BCs have
   * their `t_out` vectorset correctly, ensuring that something
   * sensible is marched through the domain, and
   *
   * 2) to use `t_out[l] == NAN` as a mask signaling that the BC nodes
   * should be included in the valid angle mask */
  if (par3_is_empty(&par)) {
    assert(eik->num_BCs[l0] > 0);
    return;
  }

  size_t num_active = par3_num_active(&par);
  assert(num_active > 0);

  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  dbl t_out[num_active][3];
  for (size_t i = 0; i < num_active; ++i)
    dbl3_copy(eik->t_out[l[i]], t_out[i]);

  bool finite[num_active];
  size_t num_finite = 0;
  for (size_t i = 0; i < num_active; ++i)
    num_finite += finite[i] = dbl3_isfinite(t_out[i]);

  /* If there aren't enough finite `t_out` vectors, then (in the
   * general case) what we need to do is solve an optimization problem
   * to compute `eik->t_out[l0]`. Since we assume c == 1, then the
   * result of this 2pt BVP would just be the gradient at `l0`, so we
   * just copy that over now. */
  if (num_finite < num_active) {
    dbl3_copy(&eik->jet[l0].fx, eik->t_out[l0]);
    return;
  }

  if (num_active == 1)
    dbl3_copy(t_out[0], eik->t_out[l0]);
  else if (num_active == 2)
    slerp2(t_out[0], t_out[1], b, eik->t_out[l0]);
  else if (num_active == 3)
    // slerp3(t_out, b, eik->t_out[l0], eik->slerp_tol);
    nlerp3(t_out, b, eik->t_out[l0]);
  else
    assert(false);

  assert(eik->num_BCs[l0] > 0 || dbl3_isfinite(eik->t_out[l0]));
}

/* Remove old tetrahedron and two-point boundary updates targeting
 * `l0` from `eik`. */
static void purge_old_updates(eik3_s *eik, size_t l0) {
  utetra_s *old_utetra;
  for (size_t i = array_size(eik->old_updates); i > 0; --i) {
    array_get(eik->old_updates, i - 1, &old_utetra);
    if (utetra_get_l(old_utetra) == l0) {
      utetra_deinit(old_utetra);
      utetra_dealloc(&old_utetra);
      array_delete(eik->old_updates, i - 1);
    }
  }

  utri_s *old_utri;
  for (size_t i = array_size(eik->old_bd_utri); i > 0; --i) {
    array_get(eik->old_bd_utri, i - 1, &old_utri);
    if (utri_get_l(old_utri) == l0) {
      utri_dealloc(&old_utri);
      array_delete(eik->old_bd_utri, i - 1);
    }
  }
}

size_t eik3_step(eik3_s *eik) {
  size_t l0 = heap_front(eik->heap);

  assert(eik->state[l0] == TRIAL);

  heap_pop(eik->heap);

  eik->state[l0] = VALID;

  if (!isfinite(eik->jet[l0].f))
    goto coda;

  /* Only compute `t_in` for scattered fields... */
  if (eik->ftype != FTYPE_POINT_SOURCE)
    compute_t_in(eik, l0);

  compute_t_out(eik, l0);

  purge_old_updates(eik, l0);
  update_neighbors(eik, l0);

coda:

  /* Increment the number of nodes that have been accepted, and mark
   * that the `eik->num_accepted`th node was `l0`. */
  eik->accepted[eik->num_accepted++] = l0;

  return l0;
}

void eik3_solve(eik3_s *eik) {
  while (heap_size(eik->heap) > 0) {
    (void)eik3_step(eik);
  }
}

bool eik3_is_solved(eik3_s const *eik) {
  return eik->num_accepted == mesh3_nverts(eik->mesh);
}

void eik3_add_trial(eik3_s *eik, size_t l, jet3 jet) {
  if (eik->state[l] == VALID) {
    log_warn("failed to add TRIAL node %lu (already VALID)", l);
    return;
  }

  if (isfinite(eik->jet[l].f)) {
    log_warn("failed to add TRIAL node %lu (finite jet)", l);
    return;
  }

  assert(eik->pos[l] == NO_INDEX);
  assert(eik->state[l] == FAR);

  eik->jet[l] = jet;
  eik->state[l] = TRIAL;
  heap_insert(eik->heap, l);

  ++eik->num_BCs[l];
}

void eik3_add_trial_w_data(eik3_s *eik, size_t l, jet3 jet,
                           dbl const hess[3][3],
                           dbl const t_in[3], dbl const t_out[3]) {
  eik3_add_trial(eik, l, jet);

  dbl33_copy(hess, eik->hess[l]);
  if (eik->ftype != FTYPE_POINT_SOURCE)
    dbl3_copy(t_in, eik->t_in[l]);
  dbl3_copy(t_out, eik->t_out[l]);
}

bool eik3_is_point_source(eik3_s const *eik, size_t l) {
  /* First, check whether the jet satisfies the conditions for being a
   * point source. If it doesn't, return early. */
  if (!jet3_is_point_source(&eik->jet[l]))
    return false;

  /* Next, check if this point source is adjacent to any other jets
   * that satisfy the conditions for being a point source. If it is,
   * then this is actually a point with boundary data for an edge
   * diffraction problem. */

  size_t nvv = mesh3_nvv(eik->mesh, l);
  size_t *vv = malloc(nvv*sizeof(size_t));
  mesh3_vv(eik->mesh, l, vv);

  bool is_point_source = true;

  size_t le[2] = {[0] = l};

  for (size_t i = 0; i < nvv; ++i) {
    /* Skip `vv[i]` if its jet doesn't look like a point source. */
    if (!jet3_is_point_source(&eik->jet[vv[i]]))
      continue;

    /* If `l` and `vv[i]` have a diffracting edge BC, then this
     * definitely isn't a point source. Bail. */
    le[1] = vv[i];
    if (eik3_has_bde_bc(eik, le)) {
      assert(mesh3_is_diff_edge(eik->mesh, le)); /* sanity check */
      is_point_source = false;
      break;
    }
  }

  free(vv);

  return is_point_source;
}

bool eik3_is_far(eik3_s const *eik, size_t l) {
  return eik->state[l] == FAR;
}

bool eik3_is_trial(eik3_s const *eik, size_t l) {
  return eik->state[l] == TRIAL;
}

bool eik3_is_valid(eik3_s const *eik, size_t l) {
  return eik->state[l] == VALID;
}

mesh3_s *eik3_get_mesh(eik3_s const *eik) {
  return eik->mesh;
}

jet3 eik3_get_jet(eik3_s const *eik, size_t l) {
  return eik->jet[l];
}

void eik3_set_jet(eik3_s *eik, size_t l, jet3 jet) {
  eik->jet[l] = jet;
}

jet3 *eik3_get_jet_ptr(eik3_s const *eik) {
  return eik->jet;
}

dbl33 *eik3_get_hess_ptr(eik3_s const *eik) {
  return eik->hess;
}

state_e *eik3_get_state_ptr(eik3_s const *eik) {
  return eik->state;
}

par3_s eik3_get_par(eik3_s const *eik, size_t l) {
  return eik->par[l];
}

void eik3_set_par(eik3_s *eik, size_t l, par3_s par) {
  eik->par[l] = par;
}

bool eik3_has_par(eik3_s const *eik, size_t l) {
  return !par3_is_empty(&eik->par[l]);
}

void eik3_get_DT(eik3_s const *eik, size_t l, dbl DT[3]) {
  memcpy(DT, &eik->jet[l].fx, 3*sizeof(dbl));
}

dbl const *eik3_get_DT_ptr(eik3_s const *eik, size_t l) {
  return &eik->jet[l].fx;
}

dbl *eik3_get_t_in_ptr(eik3_s const *eik) {
  assert(eik->ftype != FTYPE_POINT_SOURCE);
  return eik->t_in[0];
}

dbl *eik3_get_t_out_ptr(eik3_s const *eik) {
  return eik->t_out[0];
}

void eik3_add_pt_src_BCs(eik3_s *eik, size_t l, jet3 jet) {
  assert(eik->ftype == FTYPE_POINT_SOURCE);
  eik3_add_trial(eik, l, jet);
}

void eik3_add_refl_BCs(eik3_s *eik, size_t const lf[3], jet3 const jet[3],
                       dbl33 const hess[3], dbl const t_in[3][3]) {
#if JMM_DEBUG
  assert(eik->ftype == FTYPE_REFLECTION);
  assert(mesh3_is_bdf(eik->mesh, lf));

  dbl nu[3];
  mesh3_get_face_normal(eik->mesh, lf, nu);
  for (size_t i = 0; i < 3; ++i)
    if (dbl3_isfinite(&jet[i].fx))
      assert(dbl3_dot(nu, &jet[i].fx) > 0);
#endif

  for (size_t i = 0; i < 3; ++i)
    if (!eik3_is_trial(eik, lf[i]))
      eik3_add_trial(eik, lf[i], jet[i]);

  /* Set the Hessian at each face. We assume the Hessian is already
   * transformed correctly when this function is called. */
  for (size_t i = 0; i < 3; ++i)
    memcpy(&eik->hess[lf[i]], hess[i], sizeof(dbl[3][3]));

  /* If any of the edges of `lf` are diffracting edges, we want to add
   * BCs for those diffracting edges now. */
  // TODO: refactor this and the same code in `eik3_add_valid_bde` out
  // into a separate function
  for (size_t i = 0, le[2]; i < 3; ++i) {
    size_t j = (i + 1) % 3;
    le[0] = lf[i];
    le[1] = lf[j];
    if (mesh3_is_diff_edge(eik->mesh, le)) {
      dbl f[2] = {jet[i].f, jet[j].f}, Df[2][3], x[2][3];
      dbl3_copy(&jet[i].fx, Df[0]);
      dbl3_copy(&jet[j].fx, Df[1]);
      mesh3_copy_vert(eik->mesh, le[0], x[0]);
      mesh3_copy_vert(eik->mesh, le[1], x[1]);

      bb31 bb;
      bb31_init_from_3d_data(&bb, f, Df, x);
      eik3_set_bde_bc(eik, le, &bb);
    }
  }

  for (size_t i = 0; i < 3; ++i) {
    if (dbl3_isfinite(eik->t_in[lf[i]]))
      assert(dbl3_dist(t_in[i], eik->t_in[lf[i]]) < 1e-14);
    dbl3_copy(t_in[i], eik->t_in[lf[i]]);
  }
}

void eik3_add_diff_edge_BCs(eik3_s *eik, size_t const le[2],
                            bb31 const *T, dbl const rho1[2],
                            dbl3 const t_in[2]) {
  assert(eik->ftype == FTYPE_EDGE_DIFFRACTION);

  assert(mesh3_is_diff_edge(eik->mesh, le));

  dbl f[2] = {T->c[0], T->c[3]};

  /* Add trial nodes at endpoints of `le` */
  for (size_t i = 0; i < 2; ++i)
    if (!eik3_is_trial(eik, le[i]))
      eik3_add_trial(eik, le[i], jet3_make_point_source(f[i]));

  /* Set the diff edge BCs using `T` */
  eik3_set_bde_bc(eik, le, T);

  /* Set `t_in` */
  for (size_t i = 0; i < 2; ++i)
    dbl3_copy(t_in[i], eik->t_in[le[i]]);

  /* Set `rho1` at the diffracting edge vertices */
  for (size_t i = 0; i < 2; ++i) {
    if (alist_contains(eik->rho1, &le[i])) {
      dbl rho1_;
      alist_get_by_key(eik->rho1, &le[i], &rho1_);
      assert(fabs(rho1[i] - rho1_) < 1e-14);
    } else {
      alist_append(eik->rho1, &le[i], &rho1[i]);
    }
  }
}

void eik3_set_bde_bc(eik3_s *eik, size_t const le[2], bb31 const *bb) {
  bde_bc_s bc = make_bde_bc(le, bb), *this_bc;

  /* Check if a boundary condition for the current edge exists
   * already, and if so, replace the boundary data.  */
  for (size_t i = 0; i < array_size(eik->bde_bc); ++i) {
    this_bc = array_get_ptr(eik->bde_bc, i);
    if (bc.le[0] == this_bc->le[0] && bc.le[1] == this_bc->le[1]) {
      this_bc->bb = *bb;
      return;
    }
  }

  /* Didn't find an existing BC for edge `le`, so add one. */
  array_append(eik->bde_bc, &bc);
}

bool eik3_get_bde_bc(eik3_s const *eik, size_t const le[2], bb31 *bb) {
  size_t le_sorted[2] = {le[0], le[1]};
  SORT2(le_sorted[0], le_sorted[1]);

  bool swapped = le_sorted[0] != le[0];

  /* Scan through `eik->bde_bc` for `le` and grab the corresponding
   * boundary, making sure to reverse `bb` in case the indices of `le`
   * weren't initially sorted.
   *
   * (This makes it so that from the perspective of the caller,
   * everything behaves correctly and transparently. Evaluating the
   * resulting `bb31` stored in `bb` at one of the edge endpoints just
   * does the right thing.) */
  bde_bc_s *bc;
  for (size_t i = 0; i < array_size(eik->bde_bc); ++i) {
    bc = array_get_ptr(eik->bde_bc, i);
    if (le_sorted[0] == bc->le[0] && le_sorted[1] == bc->le[1]) {
      *bb = bc->bb;
      if (swapped)
        bb31_reverse(bb);
      return true;
    }
  }

  /* We didn't find any boundary data for the edge indexed by `le`. */
  return false;

  // TODO: one little improvement we could make here is, if we *don't*
  // find any boundary data in `bde_bc`, but `le` *does* correspond to
  // a diffracting edge, we could just compute the boundary data and
  // return it now. This would make setting up edge diffraction
  // problems a little more streamlined. That is, after solving the
  // incident wave, grab use `eik3_get_bde_bc` to get the relevant
  // boundary data from an enveloped diffracting edge.
}

bool eik3_has_bde_bc(eik3_s const *eik, size_t const le[2]) {
  assert(le[0] != le[1]);
  size_t le_sorted[2] = {MIN(le[0], le[1]), MAX(le[0], le[1])};
  bde_bc_s *bc;
  for (size_t i = 0; i < array_size(eik->bde_bc); ++i) {
    bc = array_get_ptr(eik->bde_bc, i);
    if (bc->le[0] == le_sorted[0] && bc->le[1] == le_sorted[1]) {
      assert(mesh3_is_diff_edge(eik->mesh, bc->le));
      return true;
    }
  }
  return false;
}

ftype_e eik3_get_ftype(eik3_s const *eik) {
  return eik->ftype;
}

dbl eik3_get_slerp_tol(eik3_s const *eik) {
  return eik->slerp_tol;
}

bool eik3_has_BCs(eik3_s const *eik, size_t l) {
  return eik->num_BCs[l] > 0;
}

static void transport_dbl(eik3_s const *eik, size_t l0, dbl *values) {
  par3_s par = eik->par[l0];

  if (par3_is_empty(&par))
    return;

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  values[l0] = 0;
  for (size_t i = 0; i < num_active; ++i) {
    assert(isfinite(values[l[i]]));
    values[l0] += b[i]*values[l[i]];
  }
}

void eik3_transport_dbl(eik3_s const *eik, dbl *values, bool skip_filled) {
  assert(eik3_is_solved(eik));

  size_t nverts = mesh3_nverts(eik->mesh);
  for (size_t i = 0; i < nverts; ++i) {
    size_t l0 = eik->accepted[i];
    if (skip_filled && !isnan(values[l0]))
      continue;
    transport_dbl(eik, l0, values);
  }
}

static void transport_dblz(eik3_s const *eik, size_t l0, dblz *values) {
  par3_s par = eik->par[l0];

  if (par3_is_empty(&par))
    return;

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  values[l0] = 0;
  for (size_t i = 0; i < num_active; ++i) {
    assert(isfinite(creal(values[l[i]])));
    assert(isfinite(cimag(values[l[i]])));
    values[l0] += b[i]*values[l[i]];
  }
}

void eik3_transport_dblz(eik3_s const *eik, dblz *values, bool skip_filled) {
  assert(eik3_is_solved(eik));

  size_t nverts = mesh3_nverts(eik->mesh);
  for (size_t i = 0; i < nverts; ++i) {
    size_t l0 = eik->accepted[i];
    dblz z = values[l0];
    if (skip_filled && !isnan(creal(z)) && !isnan(cimag(z)))
      continue;
    transport_dblz(eik, l0, values);
  }
}

static void transport_curvature(eik3_s const *eik, size_t l0, dbl *kappa) {
  if (isfinite(kappa[l0]))
    return;

  par3_s par = eik->par[l0];

  if (par3_is_empty(&par))
    return;

  size_t num_active = par3_num_active(&par);
  size_t l[num_active];
  dbl b[num_active];
  par3_get_active(&par, l, b);

  bool finite[num_active];

  size_t num_finite = 0;
  for (size_t i = 0; i < num_active; ++i) {
    num_finite += finite[i] = isfinite(kappa[l[i]]);
    if (!finite[i])
      assert(isinf(kappa[l[i]]));
  }

  if (num_finite < num_active) {
    assert(num_finite == 1);
    assert(num_active == 2);

    dbl rho = 0;
    for (size_t i = 0; i < num_active; ++i)
      rho += b[i]/kappa[l[i]];
    kappa[l0] = 1/rho;

    return;
  }

  kappa[l0] = 0;
  for (size_t i = 0; i < num_active; ++i) {
    assert(isfinite(kappa[l[i]]));
    kappa[l0] += b[i]*kappa[l[i]];
  }
}

void eik3_transport_curvature(eik3_s const *eik, dbl *kappa, bool skip_filled) {
  assert(eik3_is_solved(eik));

  size_t nverts = mesh3_nverts(eik->mesh);
  for (size_t i = 0; i < nverts; ++i) {
    size_t l0 = eik->accepted[i];
    if (skip_filled && !isnan(kappa[l0]))
      continue;
    transport_curvature(eik, l0, kappa);
  }
}

dbl eik3_get_h(eik3_s const *eik) {
  return eik->h;
}

bool eik3_get_refl_bdf_inc_on_diff_edge(eik3_s const *eik, size_t const le[2],
                                        size_t lf[3]) {
  size_t nbdf = mesh3_get_num_bdf_inc_on_edge(eik->mesh, le);
  size_t (*bdf)[3] = malloc(nbdf*sizeof(size_t[3]));
  mesh3_get_bdf_inc_on_edge(eik->mesh, le, bdf);

  size_t i;

  for (i = 0; i < nbdf; ++i) {
    if (eik3_has_BCs(eik, bdf[i][0]) &&
        eik3_has_BCs(eik, bdf[i][1]) &&
        eik3_has_BCs(eik, bdf[i][2])) {
      memcpy(lf, bdf[i], sizeof(size_t[3]));
      break;
    }
  }

  bool found = i < nbdf;

#if JMM_DEBUG
  for (i = i + 1; i < nbdf; ++i) {
    if (eik3_has_BCs(eik, bdf[i][0]) &&
        eik3_has_BCs(eik, bdf[i][1]) &&
        eik3_has_BCs(eik, bdf[i][2])) {
      assert(false);
    }
  }
#endif

  free(bdf);

  return found;
}

size_t const *eik3_get_accepted_ptr(eik3_s const *eik) {
  return eik->accepted;
}
