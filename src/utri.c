#include <jmm/utri.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <jmm/array.h>
#include <jmm/bb.h>
#include <jmm/eik3.h>
#include <jmm/hybrid.h>
#include <jmm/mat.h>
#include <jmm/mesh3.h>
#include <jmm/slerp.h>

struct utri {
  eik3_s const *eik;

  dbl lam;
  dbl f;
  dbl Df;
  dbl x_minus_xb[3];
  dbl L;

  size_t l, l0, l1;
  dbl x[3];
  dbl x0[3];
  dbl x1[3];
  dbl x1_minus_x0[3];
  bb31 T;
};

void utri_alloc(utri_s **utri) {
  *utri = malloc(sizeof(utri_s));
}

void utri_dealloc(utri_s **utri) {
  free(*utri);
  *utri = NULL;
}

static void set_lambda_stype_constant(utri_s *utri, dbl lam) {
  utri->lam = lam;

  dbl xb[3];
  dbl3_saxpy(lam, utri->x1_minus_x0, utri->x0, xb);
  dbl3_sub(utri->x, xb, utri->x_minus_xb);
  utri->L = dbl3_norm(utri->x_minus_xb);

  dbl dL_dlam = -dbl3_dot(utri->x1_minus_x0, utri->x_minus_xb)/utri->L;

  dbl b[2] = {1 - lam, lam};
  dbl T = bb31_f(&utri->T, b);

  utri->f = T + utri->L;

  dbl a[2] = {-1, 1};
  dbl dT_dlam = bb31_df(&utri->T, b, a);

  utri->Df = dT_dlam + dL_dlam;
}

static void set_lambda(utri_s *utri, dbl lam) {
  switch(eik3_get_stype(utri->eik)) {
  case STYPE_CONSTANT:
    set_lambda_stype_constant(utri, lam);
    break;
  default:
    assert(false);
  }
}

static dbl hybrid_f(dbl lam, utri_s *utri) {
  set_lambda(utri, lam);
  return utri->Df;
}

static void
rotate_jet_for_diffraction(mesh3_s const *mesh,
                           size_t l_diff, size_t l_bdv,
                           jet31t *jet)
{
  assert(mesh3_vert_incident_on_diff_edge(mesh, l_diff));
  assert(!mesh3_vert_incident_on_diff_edge(mesh, l_bdv));
  assert(mesh3_bdv(mesh, l_bdv));

  /* Get the diffracting edges incident on `l_diff` */
  size_t num_inc_diff_edges = mesh3_get_num_inc_diff_edges(mesh, l_diff);
  assert(num_inc_diff_edges > 0);
  size_t (*le)[2] = malloc(num_inc_diff_edges*sizeof(size_t[2]));
  mesh3_get_inc_diff_edges(mesh, l_diff, le);

  /* Get an incident diffracting tangent edge */
  dbl3 te;
  mesh3_get_diff_edge_tangent(mesh, le[0], te);

  /* Get the diffraction vertex */
  dbl3 xe;
  mesh3_copy_vert(mesh, l_diff, xe);

  /* Get the boundary vertex */
  dbl3 xf;
  mesh3_copy_vert(mesh, l_bdv, xf);

  /* Project the boundary vertex onto the diffracting edge */
  dbl lam = (dbl3_dot(te, xf) - dbl3_dot(te, xe))/dbl3_dot(te, te);
  dbl3 xproj;
  dbl3_saxpy(lam, te, xe, xproj);

  /* Compute the face tangent vector */
  dbl3 tf;
  dbl3_sub(xf, xproj, tf);
  dbl3_normalize(tf);

  /* Compute the cosine of the angle between the jet and the
   * diffracting edge */
  dbl cos_beta = dbl3_dot(te, jet->Df);
  dbl sin_beta = sqrt(1 - cos_beta*cos_beta);

  /* Compute the diffracted tangent vector */
  for (size_t i = 0; i < 3; ++i)
    jet->Df[i] = cos_beta*te[i] + sin_beta*tf[i];

  /* Clean up */
  free(le);
}

void utri_init(utri_s *u, eik3_s const *eik, size_t lhat, size_t const l[2]) {
  u->eik = eik;

  mesh3_s const *mesh = eik3_get_mesh(eik);

  /* Initialize `u` */

  u->f = INFINITY;

  u->lam = u->f = u->Df = u->L = NAN;
  u->x_minus_xb[0] = u->x_minus_xb[1] = u->x_minus_xb[2] = NAN;

  u->l = lhat;
  u->l0 = l[0];
  u->l1 = l[1];

  mesh3_copy_vert(mesh, u->l, u->x);
  mesh3_copy_vert(mesh, u->l0, u->x0);
  mesh3_copy_vert(mesh, u->l1, u->x1);

  dbl3_sub(u->x1, u->x0, u->x1_minus_x0);

  jet31t jet[2];
  jet[0] = eik3_get_jet(eik, l[0]);
  jet[1] = eik3_get_jet(eik, l[1]);

  bool l0_on_diff_edge = mesh3_vert_incident_on_diff_edge(mesh, u->l0);
  bool l1_on_diff_edge = mesh3_vert_incident_on_diff_edge(mesh, u->l1);

  if (l0_on_diff_edge &&
      l1_on_diff_edge &&
      eik3_has_diff_bc(eik, l)) {
    /* If we're updating from a diffracting edge with BCs, then we
     * grab the cubic polynomial giving the BCs now */
    eik3_get_diff_bc(eik, l, &u->T);
  } else {
    /* If exactly one of x0 and x1 is incident on a diffracting edge,
     * this is a boundary triangle update, and we need to rotate the
     * tangent vector incident on the diffracting edge to account for
     * diffraction. */
    size_t which = (size_t)NO_INDEX;
    if (l0_on_diff_edge ^ l1_on_diff_edge)
      which = l0_on_diff_edge ? 0 : 1;
    if (which != (size_t)NO_INDEX) {
      size_t l_diff = which == 0 ? u->l0 : u->l1;
      size_t l_bdv = which == 0 ? u->l1 : u->l0;
      rotate_jet_for_diffraction(mesh, l_diff, l_bdv, &jet[which]);
    }

    dbl Xt[2][3];
    dbl3_copy(u->x0, Xt[0]);
    dbl3_copy(u->x1, Xt[1]);
    bb31_init_from_jets(&u->T, jet, Xt);
  }
}

void utri_solve(utri_s *utri) {
  dbl lam, f[2];

  if (hybrid((hybrid_cost_func_t)hybrid_f, 0, 1, utri, &lam))
    return;

  set_lambda(utri, 0);
  f[0] = utri->f;

  set_lambda(utri, 1);
  f[1] = utri->f;

  assert(f[0] != f[1]);

  if (f[0] < f[1])
    set_lambda(utri, 0);
  else
    set_lambda(utri, 1);
}

static void get_update_inds(utri_s const *utri, size_t l[2]) {
  l[0] = utri->l0;
  l[1] = utri->l1;
}

static void get_bary_coords(utri_s const *utri, dbl b[2]) {
  b[0] = 1 - utri->lam;
  b[1] = utri->lam;
}

par3_s utri_get_par(utri_s const *u) {
  par3_s par = {.l = {[2] = NO_PARENT}, .b = {[2] = NAN}};
  get_update_inds(u, par.l);
  get_bary_coords(u, par.b);
  return par;
}

dbl utri_get_value(utri_s const *u) {
  return u->f;
}

void utri_get_jet31t(utri_s const *utri, jet31t *jet) {
  jet->f = utri->f;
  dbl3_dbl_div(utri->x_minus_xb, utri->L, jet->Df);
}

static dbl get_lag_mult(utri_s const *utri) {
  dbl const atol = 1e-15;
  if (utri->lam < atol) {
    return utri->Df;
  } else if (utri->lam > 1 - atol) {
    return -utri->Df;
  } else {
    return 0;
  }
}

bool utri_ray_is_occluded(utri_s const *utri, eik3_s const *eik) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  par3_s par = utri_get_par(utri);
  return mesh3_local_ray_is_occluded(mesh, utri->l, &par);
}

bool utri_has_interior_point_solution(utri_s const *utri) {
  dbl const atol = 1e-14;
  return (atol < utri->lam && utri->lam < 1 - atol)
    || fabs(get_lag_mult(utri)) <= atol;
}

bool utri_active_vert_is_terminal_diff_vert(utri_s const *utri,
                                            eik3_s const *eik) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  size_t la = utri_get_active_ind(utri);
  return mesh3_vert_is_terminal_diff_edge_vert(mesh, la);
}

bool utri_is_backwards(utri_s const *utri, eik3_s const *eik) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  assert(mesh3_bdv(mesh, utri->l0) && mesh3_bdv(mesh, utri->l1));

  /* We need to allow diffracted triangle updates... */
  if (mesh3_is_diff_edge(mesh, (size_t[2]) {utri->l0, utri->l1}))
    return false;

  dbl3 xhat, x0, dx;
  mesh3_copy_vert(mesh, utri->l, xhat);
  mesh3_copy_vert(mesh, utri->l0, x0);
  dbl3_sub(mesh3_get_vert_ptr(mesh, utri->l1), x0, dx);

  dbl lam = (dbl3_dot(dx, xhat) - dbl3_dot(dx, x0))/dbl3_dot(dx, dx);
  dbl3 xproj;
  dbl3_saxpy(lam, dx, x0, xproj);

  dbl3 dxhat;
  dbl3_sub(xhat, xproj, dxhat);

  jet31t jet[2] = {
    eik3_get_jet(eik, utri->l0),
    eik3_get_jet(eik, utri->l1)
  };

  bool diff[2] = {
    mesh3_vert_incident_on_diff_edge(mesh, utri->l0),
    mesh3_vert_incident_on_diff_edge(mesh, utri->l1)
  };

  /* TODO: this can happen. Example: boundary triangle update on the
   * inside of a portal which is ~one triangle wedge. Both x0 and x1
   * could be incident on diffracting edges, but in this case they
   * would be different diffracting edges. Need to think about how
   * important it is to check this case carefully. It's doable, but
   * inefficient. */
  // assert(!(diff[0] && diff[1]));

  if (diff[0] ^ diff[1]) {
    size_t l_diff = diff[0] ? utri->l0 : utri->l1;
    size_t l_bdv = diff[0] ? utri->l1 : utri->l0;
    int which = diff[0] ? 0 : 1;
    rotate_jet_for_diffraction(mesh, l_diff, l_bdv, &jet[which]);
  }

  if (dbl3_dot(dxhat, jet[0].Df) <= 0)
    return true;

  if (dbl3_dot(dxhat, jet[1].Df) <= 0)
    return true;

  return false;
}

size_t utri_get_active_ind(utri_s const *utri) {
  dbl const atol = 1e-14;
  if (utri->lam < atol)
    return utri->l0;
  if (utri->lam > 1 - atol)
    return utri->l1;
  return (size_t)NO_INDEX;
}

size_t utri_get_inactive_ind(utri_s const *utri) {
  dbl const atol = 1e-15;
  if (utri->lam < atol)
    return utri->l1;
  if (utri->lam > 1 - atol)
    return utri->l0;
  return (size_t)NO_INDEX;
}

size_t utri_get_l(utri_s const *utri) {
  return utri->l;
}

/* Check whether `u`'s `l`, `l0`, and `l1` are collinear (i.e.,
 * whether the `utri` is "degenerate"). */
bool utri_is_degenerate(utri_s const *u) {
  dbl const atol = 1e-14;
  dbl dx0[3], dx1[3];
  dbl3_sub(u->x0, u->x, dx0);
  dbl3_sub(u->x1, u->x, dx1);
  dbl dot = dbl3_dot(dx0, dx1)/(dbl3_norm(dx0)*dbl3_norm(dx1));
  return fabs(1 - fabs(dot)) < atol;
}

bool utri_has_inds(utri_s const *u, size_t lhat, uint2 const l) {
  return u->l == lhat && u->l0 == l[0] && u->l1 == l[1];
}

/* Check whether two triangle updates have the same indices. The
 * target index shoudl be the same, but we allow the indices
 * parametrizing the base of the update to be out of order with
 * respect to one another. */
bool utris_have_same_inds(utri_s const *u1, utri_s const *u2) {
  return u1->l == u2-> l && (
    (u1->l0 == u2->l0 && u1->l1 == u2->l1) ||
    (u1->l0 == u2->l1 && u1->l1 == u2->l0));
}

static dbl get_lambda(utri_s const *u) {
  return u->lam;
}

bool utris_yield_same_update(utri_s const *utri1, utri_s const *utri2) {
  dbl const atol = 1e-14;

  if (utri1 == NULL || utri2 == NULL)
    return false;

  dbl lam1 = get_lambda(utri1);
  dbl lam2 = get_lambda(utri2);

  if (fabs(lam1 - lam2) > atol)
    return false;

  jet31t jet1, jet2;
  utri_get_jet31t(utri1, &jet1);
  utri_get_jet31t(utri2, &jet2);

  return jet31t_approx_eq(&jet1, &jet2, atol);
}

#if JMM_TEST
bool utri_is_causal(utri_s const *utri) {
  dbl dx0[3], dx1[3];
  dbl3_sub(utri->x0, utri->x, dx0);
  dbl3_sub(utri->x1, utri->x, dx1);

  dbl cos01 = dbl3_dot(dx0, dx1)/(dbl3_norm(dx0)*dbl3_norm(dx1));

  return cos01 >= 0;
}

dbl utri_get_lambda(utri_s const *utri) {
  return get_lambda(utri);
}
#endif
