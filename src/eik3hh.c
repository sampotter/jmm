#include "eik3hh.h"

#include <assert.h>
#include <stdio.h>

#include "bmesh.h"
#include "eik3.h"
#include "eik3_transport.h"
#include "mat.h"

struct eik3hh {
  dbl3 xsrc;
  dbl rfac;

  eik3_s *eik;

  /* Convenience variables: */
  mesh3_s const *mesh; /* == eik3_get_mesh(eik) */
  size_t nverts; /* == mesh3_nverts(mesh) */

  dbl33 *D2T;
  dbl *spread;
  dbl *origin;
};

void eik3hh_alloc(eik3hh_s **hh) {
  *hh = malloc(sizeof(eik3hh_s));
}

static void init(eik3hh_s *hh, mesh3_s const *mesh) {
  eik3_alloc(&hh->eik);
  eik3_init(hh->eik, mesh);

  /* set up convenience variables */
  hh->mesh = mesh;
  hh->nverts = mesh3_nverts(mesh);

  hh->D2T = malloc(sizeof(dbl33)*hh->nverts);
  for (size_t l = 0; l < hh->nverts; ++l)
    dbl33_nan(hh->D2T[l]);

  hh->spread = malloc(sizeof(dbl)*hh->nverts);
  for (size_t l = 0; l < hh->nverts; ++l)
    hh->spread[l] = NAN;

  hh->origin = malloc(sizeof(dbl)*hh->nverts);
  for (size_t l = 0; l < hh->nverts; ++l)
    hh->origin[l] = NAN;
}

static void add_pt_src_BCs(eik3hh_s *hh) {
  dbl3 x;

  /* iterate over all points inside the factoring radius and
   * initialize them with the exact jet */
  for (size_t l = 0; l < hh->nverts; ++l) {
    mesh3_copy_vert(hh->mesh, l, x);
    if (dbl3_dist(x, hh->xsrc) <= hh->rfac) {
      if (eik3_is_valid(hh->eik, l))
        continue;
      jet31t jet;
      jet.f = dbl3_dist(x, hh->xsrc);
      dbl3_sub(x, hh->xsrc, jet.Df);
      dbl3_dbl_div_inplace(jet.Df, jet.f);
      eik3_add_bc(hh->eik, l, jet);
    }
  }

  /* set all `FAR` neighbors of the now `VALID` BC nodes with exact
   * data and insert them as `TRIAL` nodes into the heap */
  for (size_t l = 0; l < hh->nverts; ++l) {
    if (!eik3_is_valid(hh->eik, l))
      continue;
    size_t nvv = mesh3_nvv(hh->mesh, l);
    size_t *vv = malloc(nvv*sizeof(size_t));
    mesh3_vv(hh->mesh, l, vv);
    for (size_t i = 0; i < nvv; ++i) {
      if (!eik3_is_far(hh->eik, vv[i]))
        continue;
      mesh3_copy_vert(hh->mesh, vv[i], x);
      jet31t jet;
      jet.f = dbl3_dist(x, hh->xsrc);
      dbl3_sub(x, hh->xsrc, jet.Df);
      dbl3_dbl_div_inplace(jet.Df, jet.f);
      eik3_add_trial(hh->eik, vv[i], jet);
    }
    free(vv);
  }
}

void eik3hh_init_pt_src(eik3hh_s *hh, mesh3_s const *mesh,
                        dbl3 const xsrc, dbl rfac) {
  init(hh, mesh);

  dbl3_copy(xsrc, hh->xsrc);
  hh->rfac = rfac;

  add_pt_src_BCs(hh);
}

void eik3hh_deinit(eik3hh_s *hh) {
  eik3_deinit(hh->eik);
  eik3_dealloc(&hh->eik);

  free(hh->D2T);
  hh->D2T = NULL;

  free(hh->spread);
  hh->spread = NULL;

  free(hh->origin);
  hh->origin = NULL;
}

void eik3hh_dealloc(eik3hh_s **hh) {
  free(*hh);
  *hh = NULL;
}

static void init_D2T_pt_src(eik3hh_s *hh) {
  dbl3 x;

  /* specially initialize D2T near the point source for the direct
   * eikonal */
  for (size_t l = 0; l < hh->nverts; ++l) {
    mesh3_copy_vert(hh->mesh, l, x);
    if (dbl3_dist(x, hh->xsrc) > hh->rfac)
      continue;
    jet31t jet = eik3_get_jet(hh->eik, l);
    dbl33 outer;
    dbl3_outer(jet.Df, jet.Df, outer);
    dbl33_eye(hh->D2T[l]);
    dbl33_sub_inplace(hh->D2T[l], outer);
    dbl33_dbl_div_inplace(hh->D2T[l], jet.f);
  }

  /* we also want to initialize D2T for points which are immediately
   * downwind of the diffracting edge */
  for (size_t l = 0; l < hh->nverts; ++l) {
    par3_s par = eik3_get_par(hh->eik, l);

    /* Get active update indices and skip if they don't correspond to
     * a diffracting edge */
    size_t la[3];
    size_t na = par3_get_active_inds(&par, la);
    if (na != 2 || !mesh3_is_diff_edge(hh->mesh, la))
      continue;

    /* Get the update point */
    dbl3 xhat;
    mesh3_copy_vert(hh->mesh, l, xhat);

    /* Unit tangent vector for diffracting edge */
    dbl3 te;
    mesh3_get_diff_edge_tangent(hh->mesh, la, te);

    /* Find the point of diffraction */
    dbl3 x0, xe;
    mesh3_copy_vert(hh->mesh, la[0], x0);
    dbl lam = (dbl3_dot(te, xhat) - dbl3_dot(te, x0))/dbl3_dot(te, te);
    dbl3_saxpy(lam, te, x0, xe);

    /* Compute unit vector pointing from `xe` to `xhat` */
    dbl3 tf;
    dbl3_sub(xhat, xe, tf);
    dbl rho = dbl3_normalize(tf); // cylindrical radius for `xhat`

    /* Get eikonal jet at `xhat` */
    jet31t J = eik3_get_jet(hh->eik, l);

    /* Get the ray direction */
    /* Compute unit vector orthogonal to `te` and `tf` (this vector
     * will be orthogonal to the ray direction) */
    dbl3 q1;
    dbl3_cross(te, tf, q1);
    assert(fabs(dbl3_dot(q1, J.Df)) < 1e-13);

    /* Compute unit vector orthogonal to `q1` and the ray
     * direction */
    dbl3 q2;
    dbl3_cross(J.Df, q1, q2);

    /* Compute first curvature outer product */
    dbl33 outer1;
    dbl3_outer(q1, q1, outer1);
    dbl33_dbl_div_inplace(outer1, rho);

    /* Compute second curvature outer product */
    dbl33 outer2;
    dbl3_outer(q2, q2, outer2);
    dbl33_dbl_div_inplace(outer2, J.f);

    /* Sum them up to get the Hessian */
    dbl33_add(outer1, outer2, hh->D2T[l]);
  }
}

static bool updated_from_diff_edge(eik3_s const *eik, size_t l) {
  mesh3_s const *mesh = eik3_get_mesh(eik);

  par3_s par = eik3_get_par(eik, l);

  size_t npar = 0;
  for (size_t i = 0; i < 3; ++i)
    npar += par.l[i] != NO_PARENT;

  if (npar == 0 || npar == 3)
    return false;
  else if (npar == 1)
    return mesh3_vert_incident_on_diff_edge(mesh, par.l[0]);
  else /* npar == 2 */
    return mesh3_is_diff_edge(mesh, par.l);
}

static bool
any_cell_vert_updated_from_diff_edge(eik3_s const *eik, size_t cv[4]) {
  for (size_t i = 0; i < 4; ++i)
    if (updated_from_diff_edge(eik, cv[i]))
      return true;
  return false;
}

static bool cell_incident_on_diff_edge(mesh3_s const *mesh, size_t cv[4]) {
  for (size_t i = 0; i < 4; ++i)
    if (mesh3_vert_incident_on_diff_edge(mesh, cv[i]))
      return true;
  return false;
}

/* Approximate the Hessian at each vertex, storing the result for
 * vertex `l` at `D2T[l]`. The user should have already allocated and
 * initialized `D2T`. Entries which are `NAN` which will be filled,
 * and those which are finite will be left alone and used to compute
 * other values. */
static void approx_D2T(eik3hh_s *hh) {
  eik3_s const *eik = hh->eik;
  mesh3_s const *mesh = eik3_get_mesh(eik);
  jet31t const *jet = eik3_get_jet_ptr(eik);

  dbl33 *D2T = hh->D2T;
  dbl33 *D2T_cell = malloc(4*mesh3_ncells(mesh)*sizeof(dbl33));

  bool *has_init = malloc(hh->nverts*sizeof(bool));
  for (size_t l = 0; l < hh->nverts; ++l)
    has_init[l] = dbl33_isfinite(hh->D2T[l]);

  /* first, compute the Hessian at each cell vertex */
  for (size_t lc = 0, lv[4]; lc < mesh3_ncells(mesh); ++lc) {
    mesh3_cv(mesh, lc, lv);

    /* copy in initial values of D2T */
    for (size_t i = 0; i < 4; ++i)
      if (has_init[lv[i]])
        dbl33_copy(D2T[lv[i]], D2T_cell[4*lc + i]);

    /* get T and DT */
    jet31t J[4];
    for (size_t i = 0; i < 4; ++i) {
      J[i] = jet[lv[i]];
    }

    /* set up A */
    dbl4 A[3];
    for (size_t i = 0; i < 3; ++i) {
      dbl4_zero(A[i]);
      A[i][i] = 1;
      A[i][3] = -1;
    }

    /* get cell verts */
    dbl43 X;
    for (size_t i = 0; i < 4; ++i)
      mesh3_copy_vert(mesh, lv[i], X[i]);

    /* set up dX */
    dbl33 dX;
    for (size_t i = 0; i < 3; ++i)
      dbl3_sub(X[i], X[3], dX[i]);

    dbl33 dXinv, dXinvT;
    dbl33_copy(dX, dXinv);
    dbl33_invert(dXinv);
    dbl33_transposed(dXinv, dXinvT);

    /* set up bb33 */
    bb33 bb;
    bb33_init_from_jets(&bb, J, X);

    /* compute Hessian at each vertex */
    for (size_t i = 0; i < 4; ++i) {
      if (has_init[lv[i]])
        continue;

      dbl4 b;
      dbl4_e(b, i);

      /* compute Hessian in affine coordinates */
      dbl33 D2T_affine;
      for (size_t p = 0; p < 3; ++p) {
        for (size_t q = 0; q < 3; ++q) {
          dbl4 a[2];
          dbl4_copy(A[p], a[0]); // blech
          dbl4_copy(A[q], a[1]); // blech
          D2T_affine[p][q] = bb33_d2f(&bb, b, a);
        }
      }

      /* transform back to Cartesian and store with cell vertex */
      dbl33 tmp;
      dbl33_mul(dXinv, D2T_affine, tmp);
      dbl33_mul(tmp, dXinvT, D2T_cell[4*lc + i]);
    }
  }

  /* zero out D2T */
  for (size_t l = 0; l < hh->nverts; ++l)
    if (!has_init[l])
      dbl33_zero(D2T[l]);

  /* number of terms in weighted average */
  size_t *N = calloc(hh->nverts, sizeof(size_t));

  /* accumulate each D2T_cell entry into D2T */
  for (size_t lc = 0, cv[4]; lc < mesh3_ncells(mesh); ++lc) {
    mesh3_cv(mesh, lc, cv);

    /* skip this cell if its data is invalid */
    if (dbl33_isnan(D2T_cell[4*lc]))
      continue;

    for (size_t i = 0; i < 4; ++i) {
      if (has_init[cv[i]])
        continue;

      /* If this vertex was updated from a diff edge, don't use data
       * from a cell which is incident on a diff edge... */
      if (updated_from_diff_edge(eik, cv[i]) &&
          cell_incident_on_diff_edge(mesh, cv))
        continue;

      /* If this vertex is incident on a diff edge, don't use data
       * from a cell which was updated from a diff edge */
      if (mesh3_vert_incident_on_diff_edge(mesh, cv[i]) &&
          any_cell_vert_updated_from_diff_edge(eik, cv))
        continue;

      dbl33_add_inplace(D2T[cv[i]], D2T_cell[4*lc + i]);
      ++N[cv[i]]; /* increment number of terms in average */
    }
  }

  /* normalize each entry by the number of incident cells */
  for (size_t lv = 0; lv < hh->nverts; ++lv) {
    if (has_init[lv])
      continue;
    size_t nvc = mesh3_nvc(mesh, lv);
    dbl33_dbl_div_inplace(D2T[lv], nvc);
  }

  free(N);
  free(has_init);
  free(D2T_cell);
}

static void init_spread_pt_src(eik3hh_s *hh) {
  for (size_t l = 0; l < hh->nverts; ++l) {
    par3_s par = eik3_get_par(hh->eik, l);
    if (!dbl3_all_nan(par.b)) continue;
    dbl3 x; mesh3_copy_vert(hh->mesh, l, x);
    hh->spread[l] = 1/dbl3_dist(x, hh->xsrc);
  }
}

static void prop_spread(eik3hh_s *hh) {
  eik3_s const *eik = hh->eik;
  mesh3_s const *mesh = eik3_get_mesh(eik);
  size_t const *accepted = eik3_get_accepted_ptr(eik);

  dbl33 const *D2T = hh->D2T;
  dbl *spread = hh->spread;

  for (size_t l = 0; l < hh->nverts; ++l) {
    size_t lhat = accepted[l];
    par3_s par = eik3_get_par(eik, lhat);

    if (par3_is_empty(&par))
      continue;

    dbl3 lam, abslam;
    size_t perm[3];
    dbl33_eigvals_sym(D2T[lhat], lam);
    dbl3_abs(lam, abslam);
    dbl3_argsort(abslam, perm);

    dbl kappa1 = lam[perm[2]], kappa2 = lam[perm[1]];

    dbl spread_b = 1;
    assert(isfinite(par.b[0]));
    for (size_t j = 0; j < 3; ++j) {
      if (isfinite(par.b[j])) {
        assert(isfinite(spread[par.l[j]]));
        spread_b *= pow(spread[par.l[j]], par.b[j]);
      }
    }

    dbl3 xlam = {0, 0, 0};
    for (size_t j = 0; j < 3; ++j) {
      if (isfinite(par.b[j])) {
        dbl3 x_;
        mesh3_copy_vert(mesh, par.l[j], x_);
        for (size_t k = 0; k < 3; ++k)
          xlam[k] += par.b[j]*x_[k];
      }
    }

    dbl3 xhat;
    mesh3_copy_vert(mesh, lhat, xhat);

    dbl L = dbl3_dist(xhat, xlam);

    spread[lhat] = spread_b*exp(-L*(kappa1 + kappa2)/2);
  }
}

static void init_origin_pt_src(eik3hh_s *hh) {
  /* Set all BC nodes' origin value to 1 initially */
  array_s const *bc_inds = eik3_get_bc_inds(hh->eik);
  for (size_t i = 0, l; i < array_size(bc_inds); ++i) {
    array_get(bc_inds, i, &l);
    hh->origin[l] = 1;
  }

  /* Set all initial trial nodes' origin value to 1 initially */
  array_s const *trial_inds = eik3_get_trial_inds(hh->eik);
  for (size_t i = 0, l; i < array_size(trial_inds); ++i) {
    array_get(trial_inds, i, &l);
    hh->origin[l] = 1;
  }

  /* Set all nodes without a value that are on a diffracting edge to 0
   * initially */
  for (size_t l = 0; l < hh->nverts; ++l)
    if (isnan(hh->origin[l]) && mesh3_vert_incident_on_diff_edge(hh->mesh, l))
      hh->origin[l] = 0;
}

static void prop_origin(eik3hh_s *hh) {
  /* Transport the origin value */
  eik3_transport_dbl(hh->eik, hh->origin, true);
}

void eik3hh_solve(eik3hh_s *hh) {
  eik3_solve(hh->eik);

  init_D2T_pt_src(hh);
  approx_D2T(hh);

  init_spread_pt_src(hh);
  prop_spread(hh);

  init_origin_pt_src(hh);
  prop_origin(hh);
}

size_t eik3hh_num_bc(eik3hh_s const *hh) {
  return eik3_num_bc(hh->eik);
}

static void dump_xy_T_slice(eik3hh_s const *hh,
                            grid2_to_mesh3_mapping_s const *mapping,
                            char const *path) {
  mesh3_s const *mesh = eik3_get_mesh(hh->eik);

  /* Set up the bmesh */
  bmesh33_s *bmesh;
  bmesh33_alloc(&bmesh);
  bmesh33_init_from_mesh3_and_jets(bmesh, mesh, eik3_get_jet_ptr(hh->eik));

  FILE *fp = fopen(path, "wb");

  for (size_t l = 0; l < grid2_nind(mapping->grid); ++l) {
    dbl value = isnan(mapping->b[l][0]) ?
      NAN :
      bb33_f(bmesh33_get_bb_ptr(bmesh, mapping->lc[l]), mapping->b[l]);
    fwrite(&value, sizeof(value), 1, fp);
  }

  fclose(fp);

  bmesh33_deinit(bmesh);
  bmesh33_dealloc(&bmesh);
}

static void
dump_xy_spreading_slice(eik3hh_s const *hh,
                        grid2_to_mesh3_mapping_s const *mapping,
                        char const *path) {
  FILE *fp = fopen(path, "wb");

  dbl spread_b;
  for (size_t l = 0; l < grid2_nind(mapping->grid); ++l) {
    if (isnan(mapping->b[l][0])) {
      spread_b = NAN;
    } else {
      spread_b = 1;
      for (size_t i = 0; i < 4; ++i)
        spread_b *= pow(hh->spread[mapping->cv[l][i]], mapping->b[l][i]);
    }
    fwrite(&spread_b, sizeof(spread_b), 1, fp);
  }

  fclose(fp);
}

static void
dump_xy_origin_slice(eik3hh_s const *hh,
                     grid2_to_mesh3_mapping_s const *mapping,
                     char const *path) {
  FILE *fp = fopen(path, "wb");

  dbl origin_b;
  for (size_t l = 0; l < grid2_nind(mapping->grid); ++l) {
    if (isnan(mapping->b[l][0])) {
      origin_b = NAN;
    } else {
      origin_b = 0;
      for (size_t i = 0; i < 4; ++i)
        origin_b += mapping->b[l][i]*hh->origin[mapping->cv[l][i]];
    }
    fwrite(&origin_b, sizeof(origin_b), 1, fp);
  }

  fclose(fp);
}

void eik3hh_dump_xy_slice(eik3hh_s const *hh,
                          grid2_to_mesh3_mapping_s const *mapping,
                          field_e field, char const *path) {
  switch (field) {
  case FIELD_T:
    dump_xy_T_slice(hh, mapping, path);
    break;
  case FIELD_SPREADING:
    dump_xy_spreading_slice(hh, mapping, path);
    break;
  case FIELD_ORIGIN:
    dump_xy_origin_slice(hh, mapping, path);
    break;
  default:
    assert(false);
  }
}

eik3_s const *eik3hh_get_eik_ptr(eik3hh_s const *hh) {
  return hh->eik;
}

jet31t const *eik3hh_get_jet_ptr(eik3hh_s const *hh) {
  return eik3_get_jet_ptr(hh->eik);
}
