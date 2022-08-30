#include <jmm/eik3hh_branch.h>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <jmm/array.h>
#include <jmm/bmesh.h>
#include <jmm/eik3.h>
#include <jmm/eik3hh.h>
#include <jmm/mat.h>
#include <jmm/mesh2.h>
#include <jmm/rtree.h>

#include <jmm/util.h>

struct eik3hh_branch {
  eik3hh_s const *hh;
  eik3_s *eik;
  eik3hh_branch_type_e type;
  size_t index;
  dbl33 *D2T;
  dbl *spread;
  dbl *origin;
  eik3hh_branch_s const *parent;
  array_s *children;
};

void eik3hh_branch_alloc(eik3hh_branch_s **branch) {
  *branch = malloc(sizeof(eik3hh_branch_s));
}

static void init(eik3hh_branch_s *branch, eik3hh_s const *hh,
                 eik3hh_branch_s const *parent, eik3hh_branch_type_e type) {
  mesh3_s const *mesh = eik3hh_get_mesh(hh);
  size_t nverts = mesh3_nverts(mesh);

  branch->hh = hh;

  eik3_alloc(&branch->eik);
  eik3_init(branch->eik, mesh);

  branch->type = type;

  branch->D2T = malloc(sizeof(dbl33)*nverts);
  for (size_t l = 0; l < nverts; ++l)
    dbl33_nan(branch->D2T[l]);

  branch->spread = malloc(sizeof(dbl)*nverts);
  for (size_t l = 0; l < nverts; ++l)
    branch->spread[l] = NAN;

  branch->origin = malloc(sizeof(dbl)*nverts);
  for (size_t l = 0; l < nverts; ++l)
    branch->origin[l] = NAN;

  /* The parent of a point source branch should always be NULL */
  assert(type != EIK3HH_BRANCH_TYPE_PT_SRC || parent == NULL);
  branch->parent = parent;

  array_alloc(&branch->children);
  array_init(branch->children, sizeof(eik3hh_branch_s *), ARRAY_DEFAULT_CAPACITY);
}

void eik3hh_branch_init_pt_src(eik3hh_branch_s *branch, eik3hh_s const *hh,
                               dbl3 const xsrc) {
  init(branch, hh, NULL, EIK3HH_BRANCH_TYPE_PT_SRC);

  /* Set the branch index to be the index of the point source */
  branch->index = mesh3_get_vert_index(eik3hh_get_mesh(hh), xsrc);
  assert(branch->index != (size_t)NO_INDEX);

  /* Set up the point source boundary conditions for eik */
  eik3_add_pt_src_bcs(branch->eik, xsrc, eik3hh_get_rfac(hh));
}

void eik3hh_branch_init_refl(eik3hh_branch_s *branch,
                             eik3hh_branch_s const *parent,
                             size_t refl_index) {
  eik3hh_s const *hh = parent->hh;

  init(branch, hh, parent, EIK3HH_BRANCH_TYPE_REFL);

  /* Set the branch index to be the index of the reflector */
  branch->index = refl_index;
  assert(branch->type != parent->type || branch->index != parent->index);

  /* Set up the reflection BCs */
  eik3_add_refl_bcs(branch->eik, parent->eik, refl_index);
}

void eik3hh_branch_deinit(eik3hh_branch_s *branch, bool free_children) {
  branch->hh = NULL;

  eik3_deinit(branch->eik);
  eik3_dealloc(&branch->eik);

  branch->type = EIK3HH_BRANCH_TYPE_UNINITIALIZED;
  branch->index = (size_t)NO_INDEX;

  free(branch->D2T);
  branch->D2T = NULL;

  free(branch->spread);
  branch->spread = NULL;

  free(branch->origin);
  branch->origin = NULL;

  /* Recursively free children */
  if (free_children) {
    for (size_t i = 0; i < array_size(branch->children); ++i) {
      eik3hh_branch_s *child = array_get_ptr(branch->children, i);
      eik3hh_branch_deinit(child, true);
      eik3hh_branch_dealloc(&child);
    }
  }
}

void eik3hh_branch_dealloc(eik3hh_branch_s **hh) {
  free(*hh);
  *hh = NULL;
}

static void init_D2T_downwind_from_diff_edges(eik3_s const *eik, dbl33 *D2T) {
  mesh3_s const *mesh = eik3_get_mesh(eik);
  size_t nverts = mesh3_nverts(mesh);

  for (size_t l = 0; l < nverts; ++l) {
    par3_s par = eik3_get_par(eik, l);

    /* Get active update indices and skip if they don't correspond to
     * a diffracting edge */
    size_t la[3];
    size_t na = par3_get_active_inds(&par, la);
    if (na != 2 || !mesh3_is_diff_edge(mesh, la))
      continue;

    /* Get the update point */
    dbl3 xhat;
    mesh3_copy_vert(mesh, l, xhat);

    /* Unit tangent vector for diffracting edge */
    dbl3 te;
    mesh3_get_diff_edge_tangent(mesh, la, te);

    /* Find the point of diffraction */
    dbl3 x0, xe;
    mesh3_copy_vert(mesh, la[0], x0);
    dbl lam = (dbl3_dot(te, xhat) - dbl3_dot(te, x0))/dbl3_dot(te, te);
    dbl3_saxpy(lam, te, x0, xe);

    /* Compute unit vector pointing from `xe` to `xhat` */
    dbl3 tf;
    dbl3_sub(xhat, xe, tf);
    dbl rho = dbl3_normalize(tf); // cylindrical radius for `xhat`

    /* Get eikonal jet at `xhat` */
    jet31t J = eik3_get_jet(eik, l);

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
    dbl33_add(outer1, outer2, D2T[l]);
  }
}

static void init_D2T_pt_src(eik3hh_branch_s *branch) {
  mesh3_s const *mesh = eik3hh_get_mesh(branch->hh);
  size_t nverts = mesh3_nverts(mesh);

  /* Initialize D2T analytically in the factoring ball around the
   * point source */
  size_t lsrc = branch->index;
  dbl3 xsrc; mesh3_copy_vert(mesh, lsrc, xsrc);
  dbl rfac = eik3hh_get_rfac(branch->hh);
  dbl3 x;
  for (size_t l = 0; l < nverts; ++l) {
    mesh3_copy_vert(mesh, l, x);
    if (dbl3_dist(x, xsrc) > rfac)
      continue;
    jet31t jet = eik3_get_jet(branch->eik, l);
    dbl33 outer;
    dbl3_outer(jet.Df, jet.Df, outer);
    dbl33_eye(branch->D2T[l]);
    dbl33_sub_inplace(branch->D2T[l], outer);
    dbl33_dbl_div_inplace(branch->D2T[l], jet.f);
  }

  init_D2T_downwind_from_diff_edges(branch->eik, branch->D2T);
}

static void init_D2T_refl(eik3hh_branch_s *branch, dbl33 const *D2T_in) {
  mesh3_s const *mesh = eik3hh_get_mesh(branch->hh);

  size_t refl_index = branch->index;
  size_t nf = mesh3_get_reflector_size(mesh, refl_index);
  uint3 *lf = malloc(nf*sizeof(uint3));
  mesh3_get_reflector(mesh, refl_index, lf);

  /* Reflect the incident Hessian across each face using the image
   * method (since the speed of sound is constant). */
  dbl33 *D2T = branch->D2T;
  for (size_t i = 0; i < nf; ++i) {
    dbl33 R; mesh3_get_R_for_face(mesh, lf[i], R);
    for (size_t j = 0; j < 3; ++j) {
      size_t l = lf[i][j];
      if (dbl33_isfinite(D2T[l]))
        continue;
      dbl33_conj(D2T_in[l], R, D2T[l]);
    }
  }

  free(lf);

  init_D2T_downwind_from_diff_edges(branch->eik, branch->D2T);
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

static bool cell_incident_on_diff_edge(mesh3_s const *mesh, size_t cv[4]) {
  for (size_t i = 0; i < 4; ++i)
    if (mesh3_vert_incident_on_diff_edge(mesh, cv[i]))
      return true;
  return false;
}

static bool
any_cell_vert_updated_from_diff_edge(eik3_s const *eik, size_t cv[4]) {
  for (size_t i = 0; i < 4; ++i)
    if (updated_from_diff_edge(eik, cv[i]))
      return true;
  return false;
}

/* Approximate the Hessian at each vertex, storing the result for
 * vertex `l` at `D2T[l]`. The user should have already allocated and
 * initialized `D2T`. Entries which are `NAN` which will be filled,
 * and those which are finite will be left alone and used to compute
 * other values. */
static void approx_D2T(eik3hh_branch_s *branch) {
  mesh3_s const *mesh = eik3hh_get_mesh(branch->hh);
  size_t nverts = mesh3_nverts(mesh);

  eik3_s const *eik = branch->eik;
  jet31t const *jet = eik3_get_jet_ptr(eik);

  dbl33 *D2T = branch->D2T;
  dbl33 *D2T_cell = malloc(4*mesh3_ncells(mesh)*sizeof(dbl33));

  bool *has_init = malloc(nverts*sizeof(bool));
  for (size_t l = 0; l < nverts; ++l)
    has_init[l] = dbl33_isfinite(branch->D2T[l]);

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
  for (size_t l = 0; l < nverts; ++l)
    if (!has_init[l])
      dbl33_zero(D2T[l]);

  /* number of terms in weighted average */
  size_t *N = calloc(nverts, sizeof(size_t));

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
  for (size_t lv = 0; lv < nverts; ++lv) {
    if (has_init[lv])
      continue;
    size_t nvc = mesh3_nvc(mesh, lv);
    dbl33_dbl_div_inplace(D2T[lv], nvc);
  }

  free(N);
  free(has_init);
  free(D2T_cell);
}

static void init_spread_pt_src(eik3hh_branch_s *branch) {
  mesh3_s const *mesh = eik3hh_get_mesh(branch->hh);
  size_t lsrc = branch->index;
  dbl3 xsrc;
  mesh3_copy_vert(mesh, lsrc, xsrc);
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l) {
    par3_s par = eik3_get_par(branch->eik, l);
    if (dbl3_all_nan(par.b) || par3_has_active_parent(&par, lsrc)) {
      dbl3 x;
      mesh3_copy_vert(mesh, l, x);
      branch->spread[l] = 1/dbl3_dist(x, xsrc);
    }
  }
}

static void init_spread_refl(eik3hh_branch_s *branch, dbl const *spread_in) {
  mesh3_s const *mesh = eik3hh_get_mesh(branch->hh);

  size_t refl_index = branch->index;
  size_t nf = mesh3_get_reflector_size(mesh, refl_index);
  uint3 *lf = malloc(nf*sizeof(uint3));
  mesh3_get_reflector(mesh, refl_index, lf);

  dbl *spread = branch->spread;
  for (size_t i = 0; i < nf; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      size_t l = lf[i][j];
      if (isnan(spread[l]))
        spread[l] = spread_in[l];
    }
  }
}

static void prop_spread(eik3hh_branch_s *branch) {
  mesh3_s const *mesh = eik3hh_get_mesh(branch->hh);
  size_t nverts = mesh3_nverts(mesh);

  eik3_s const *eik = branch->eik;
  size_t const *accepted = eik3_get_accepted_ptr(eik);
  dbl33 const *D2T = branch->D2T;
  dbl *spread = branch->spread;

  for (size_t l = 0; l < nverts; ++l) {
    size_t lhat = accepted[l];
    if (!isnan(spread[lhat]))
      continue;

    par3_s par = eik3_get_par(eik, lhat);
    if (par3_is_empty(&par))
      continue;

    dbl3 lam, abslam;
    size_t perm[3];
    dbl33_eigvals_sym(D2T[lhat], lam);
    assert(dbl3_isfinite(lam));
    dbl3_abs(lam, abslam);
    dbl3_argsort(abslam, perm);
    assert(abslam[perm[0]] <= abslam[perm[1]] && abslam[perm[1]] <= abslam[perm[2]]);

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

void eik3hh_branch_solve(eik3hh_branch_s *branch, bool verbose) {
  eik3_s *eik = branch->eik;
  mesh3_s const *mesh = eik3_get_mesh(eik);

  if (verbose)
    toc();

  if (verbose) {
    if (branch->type == EIK3HH_BRANCH_TYPE_PT_SRC)
      printf("Set up point source problem:\n");
    else if (branch->type == EIK3HH_BRANCH_TYPE_REFL)
      printf("Solving reflection problem (refl_index = %lu):\n", branch->index);
    printf("- number of points inside factoring radius: %lu\n",
           eik3_num_valid(eik));
    printf("- number of TRIAL points adjacent to factored nodes: %lu\n",
           eik3_num_trial(eik));
  }

  if (eik3_num_trial(eik) == 0) {
    printf("ERROR: didn't insert any TRIAL points!\n");
    exit(EXIT_FAILURE);
  }

  if (eik3_solve(eik) == JMM_ERROR_RUNTIME_ERROR) {
    size_t skipped = mesh3_nverts(mesh) - eik3_num_valid(eik);
    printf("- WARNING: didn't relax all points (skipped %lu)\n", skipped);
    if (!eik3_brute_force_remaining(eik))
      printf("- WARNING: failed to brute force skipped points\n");
  }

  if (branch->type == EIK3HH_BRANCH_TYPE_PT_SRC) {
    eik3_init_org_from_BCs(eik, branch->origin);
    init_D2T_pt_src(branch);
    init_spread_pt_src(branch);
  }

  else if (branch->type == EIK3HH_BRANCH_TYPE_REFL) {
    eik3hh_branch_s const *parent = branch->parent;
    assert(parent != NULL);
    eik3_init_org_for_refl(eik, branch->origin, branch->index, parent->origin);
    init_D2T_refl(branch, parent->D2T);
    init_spread_refl(branch, parent->spread);
  }

  eik3_prop_org(eik, branch->origin);
  approx_D2T(branch);
  prop_spread(branch);

  size_t num_viz_skipped = 0;
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    if (branch->origin[l] >= 0.5 && !eik3_is_valid(eik, l))
      ++num_viz_skipped;

  if (verbose) {
    if (num_viz_skipped > 0)
      printf("- WARNING: %lu visible nodes were skipped\n", num_viz_skipped);
    printf("- solved [%1.2gs]\n", toc());
  }
}

static void dump_xy_T_slice(eik3hh_branch_s const *branch,
                            grid2_to_mesh3_mapping_s const *mapping,
                            char const *path) {
  mesh3_s const *mesh = eik3hh_get_mesh(branch->hh);

  /* Set up the bmesh */
  bmesh33_s *bmesh;
  bmesh33_alloc(&bmesh);
  bmesh33_init_from_mesh3_and_jets(bmesh, mesh, eik3_get_jet_ptr(branch->eik));

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
dump_xy_spreading_slice(eik3hh_branch_s const *branch,
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
        spread_b *= pow(branch->spread[mapping->cv[l][i]], mapping->b[l][i]);
    }
    fwrite(&spread_b, sizeof(spread_b), 1, fp);
  }

  fclose(fp);
}

static void
dump_xy_origin_slice(eik3hh_branch_s const *branch,
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
        origin_b += mapping->b[l][i]*branch->origin[mapping->cv[l][i]];
    }
    fwrite(&origin_b, sizeof(origin_b), 1, fp);
  }

  fclose(fp);
}

eik3_s *eik3hh_branch_get_eik(eik3hh_branch_s *branch) {
  return branch->eik;
}

array_s *eik3hh_branch_get_children(eik3hh_branch_s *branch) {
  return branch->children;
}

dbl const *eik3hh_branch_get_spread(eik3hh_branch_s const *branch) {
  return branch->spread;
}

dbl const *eik3hh_branch_get_org(eik3hh_branch_s const *branch) {
  return branch->origin;
}

size_t eik3hh_branch_get_earliest_refl(eik3hh_branch_s const *branch) {
  mesh3_s const *mesh = eik3hh_get_mesh(branch->hh);
  size_t num_refl = mesh3_get_num_reflectors(mesh);
  size_t min_refl_index = (size_t)NO_INDEX;
  dbl min_T = INFINITY;
  for (size_t refl_index = 0; refl_index < num_refl; ++refl_index) {
    if (branch->type == EIK3HH_BRANCH_TYPE_REFL
        && branch->index == refl_index) {
      continue;
    }
    size_t nf = mesh3_get_reflector_size(mesh, refl_index);
    uint3 *lf = malloc(nf*sizeof(uint3));
    mesh3_get_reflector(mesh, refl_index, lf);
    for (size_t i = 0; i < nf; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        dbl T = eik3_get_T(branch->eik, lf[i][j]);
        if (T < min_T) {
          min_T = T;
          min_refl_index = refl_index;
        }
      }
    }
    free(lf);
  }
  return min_refl_index;
}

array_s *eik3hh_branch_get_visible_refls(eik3hh_branch_s const *branch) {
  mesh3_s const *mesh = eik3_get_mesh(branch->eik);

  dbl const *org_in = branch->origin;
#if JMM_DEBUG
  for (size_t l = 0; l < mesh3_nverts(mesh); ++l)
    assert(isfinite(org_in[l]));
#endif

  array_s *refl_inds;
  array_alloc(&refl_inds);
  array_init(refl_inds, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  size_t num_refl = mesh3_get_num_reflectors(mesh);
  for (size_t refl_ind = 0; refl_ind < num_refl; ++refl_ind) {
    /* Don't repeat reflections */
    if (branch->type == EIK3HH_BRANCH_TYPE_REFL && branch->index == refl_ind)
      continue;

    /* Get the faces of the current reflector */
    size_t nf = mesh3_get_reflector_size(mesh, refl_ind);
    uint3 *lf = malloc(nf*sizeof(uint3));
    mesh3_get_reflector(mesh, refl_ind, lf);

    /* This reflector is visible if any of the incident origins are
     * greater than 1/2. */
    bool visible = false;
    for (size_t i = 0; i < nf; ++i) {
      if (visible) break;
      for (size_t j = 0; j < 3; ++j) {
        if (org_in[lf[i][j]] > 0.5 && !array_contains(refl_inds, &refl_ind)) {
          array_append(refl_inds, &refl_ind);
          visible = true;
          break;
        }
      }
    }

    free(lf);
  }

  return refl_inds;
}

eik3hh_branch_s *
eik3hh_branch_add_refl(eik3hh_branch_s const *branch, size_t refl_index) {
  /* Make sure this isn't the same reflection */
  assert(branch->type != EIK3HH_BRANCH_TYPE_REFL
         || branch->index != refl_index);

  /* Make sure we haven't done this reflection already */
  for (size_t i = 0; i < array_size(branch->children); ++i) {
    eik3hh_branch_s const *child = array_get_ptr(branch->children, i);
    if (child->type == EIK3HH_BRANCH_TYPE_REFL
        && child->index == refl_index)
      assert(false);
  }

  eik3hh_branch_s *refl;
  eik3hh_branch_alloc(&refl);
  eik3hh_branch_init_refl(refl, branch, refl_index);

  array_append(branch->children, &refl);

  return refl;
}

void eik3hh_branch_dump_jet(eik3hh_branch_s const *branch, char const *path) {
  eik3_dump_jet(branch->eik, path);
}

void eik3hh_branch_dump_org(eik3hh_branch_s const *branch, char const *path) {
  size_t nverts = mesh3_nverts(eik3_get_mesh(branch->eik));

  FILE *fp = fopen(path, "wb");
  fwrite(branch->origin, sizeof(dbl), nverts, fp);
  fclose(fp);
}

void eik3hh_branch_dump_spread(eik3hh_branch_s const *branch, char const *path) {
  size_t nverts = mesh3_nverts(eik3_get_mesh(branch->eik));

  FILE *fp = fopen(path, "wb");
  fwrite(branch->spread, sizeof(dbl), nverts, fp);
  fclose(fp);
}

void eik3hh_branch_dump_xy_slice(eik3hh_branch_s const *branch,
                                 grid2_to_mesh3_mapping_s const *mapping,
                                 field_e field, char const *path) {
  switch (field) {
  case FIELD_T:
    dump_xy_T_slice(branch, mapping, path);
    break;
  case FIELD_SPREADING:
    dump_xy_spreading_slice(branch, mapping, path);
    break;
  case FIELD_ORIGIN:
    dump_xy_origin_slice(branch, mapping, path);
    break;
  default:
    assert(false);
  }
}

dbl squash(dbl x, size_t n) {
  dbl y = x;
  for (size_t i = 0; i < n; ++i)
    y = (3*y - 2*y*y)*y;
  return y;
}

void eik3hh_branch_render_frames(eik3hh_branch_s const *branch,
                                 camera_s const *camera,
                                 dbl T0, dbl T1, dbl frames_per_meter,
                                 bool verbose) {
  mesh3_s const *mesh = eik3_get_mesh(branch->eik);

  mesh2_s *surface_mesh = mesh3_get_surface_mesh(mesh);

  bmesh33_s *bmesh;
  bmesh33_alloc(&bmesh);
  bmesh33_init_from_mesh3_and_jets(bmesh, mesh, eik3_get_jet_ptr(branch->eik));

  size_t num_frames = floor(frames_per_meter*(T1 - T0));
  if (verbose)
    printf("rendering %lu frames\n", num_frames);

  dbl *T = malloc(num_frames*sizeof(dbl));
  for (size_t i = 0; i < num_frames; ++i)
    T[i] = T0 + i/frames_per_meter;

  for (size_t i = 0; i < num_frames; ++i) {
    if (verbose)
      printf("frame %lu/%lu (T = %g s)\n", i + 1, num_frames, T[i]);

    rtree_s *rtree;
    rtree_alloc(&rtree);
    rtree_init(rtree, 16, RTREE_SPLIT_STRATEGY_SURFACE_AREA);

    rtree_insert_mesh2(rtree, surface_mesh);

    bmesh33_s *level_bmesh = bmesh33_restrict_to_level(bmesh, T[i]);
    rtree_insert_bmesh33(rtree, level_bmesh);

    rtree_build(rtree);

    size_t npix = camera->dim[0]*camera->dim[1];

    dbl4 *img = malloc(npix*sizeof(dbl4));

    dbl3 surf_rgb = {0.54, 0.54, 0.54};
    dbl3 eik_rgb = {1.0, 1.0, 1.0};

    dbl surf_alpha = 0.5;
    dbl eik_alpha = 1;

    for (size_t i = 0, l = 0; i < camera->dim[0]; ++i) {
      for (size_t j = 0; j < camera->dim[1]; ++j, ++l) {
        ray3 ray = camera_get_ray_for_index(camera, i, j);

        isect isect;
        rtree_intersect(rtree, &ray, &isect, NULL);

        img[l][0] = 0;
        img[l][1] = 0;
        img[l][2] = 0;
        img[l][3] = isfinite(isect.t) ? 1 : 0;

        dbl transparency = 1;
        dbl const *rgb = NULL;
        dbl3 n;

        while (isfinite(isect.t)) {
          robj_type_e robj_type = robj_get_type(isect.obj);
          void const *robj_data = robj_get_data(isect.obj);

          dbl alpha = 1, scale = 1;

          /* Increment the distance along the ray */
          dbl3_saxpy_inplace(isect.t, ray.dir, ray.org);

          /* Update the current alpha and RGB value */
          switch (robj_type) {
          case ROBJ_MESH2_TRI:
            alpha *= surf_alpha;
            rgb = &surf_rgb[0];
            break;
          case ROBJ_BMESH33_CELL:
            alpha *= eik_alpha;
            rgb = &eik_rgb[0];
            break;
          default:
            assert(false);
          }

          if (robj_type == ROBJ_BMESH33_CELL) {
            bmesh33_cell_s const *bmesh33_cell = robj_data;
            dbl spread_interp, org_interp;
            if (bmesh33_cell->bmesh == level_bmesh) {
              spread_interp = mesh3_linterp(
                mesh, eik3hh_branch_get_spread(branch), ray.org);
              org_interp = mesh3_linterp(
                mesh, eik3hh_branch_get_org(branch), ray.org);
            } else {
              assert(false);
            }
            /* Convert the interpolated spreading factor to dB */
            dbl spread_dB = 20*log10(fmax(1e-16, spread_interp));
            /* Clamp and map the range [-60 dB, 0 dB] to [0, 1] for
             * use as a scaling factor */
            dbl spread_mapped = fmax(0, fmin(1, 1 - spread_dB/(-90)));
            alpha *= spread_mapped*squash(org_interp, 2);
          }

          /* Get the surface normal and dot it with the eye vector for
           * Lambertian shading */
          switch (robj_type) {
          case ROBJ_MESH2_TRI:
            mesh2_tri_s const *mesh2_tri = robj_data;
            mesh2_get_unit_surface_normal(surface_mesh, mesh2_tri->l, n);
            break;
          case ROBJ_BMESH33_CELL:
            bmesh33_cell_s const *bmesh33_cell = robj_data;
            bmesh33_cell_Df(bmesh33_cell, ray.org, n);
            dbl3_normalize(n);
            break;
          default:
            assert(false);
          }
          scale *= fabs(dbl3_dot(n, ray.dir));

          /* We're raymarching, so do backwards alpha blending */
          dbl3_saxpy_inplace(scale*alpha, rgb, img[l]);

          /* Update transparency for early stopping */
          transparency *= 1 - alpha;
          if (transparency < 1e-3)
            break;

          /* Advance the start of the ray and keep tracing.
           *
           * NOTE: if we have multiple overlapping intersections, we
           * might trip them repeatedly, so we need to skip any
           * intersections with a distance of zero here. */
          rtree_intersect(rtree, &ray, &isect, isect.obj);
          while (isect.t < EPS) {
            dbl3_saxpy_inplace(EPS, ray.dir, ray.org);
            rtree_intersect(rtree, &ray, &isect, isect.obj);
          }
        }
      }
    }

    char filename[128];
    snprintf(filename, 128, "image%04lu.bin", i);

    FILE *fp = fopen(filename, "wb");
    fwrite(img, sizeof(dbl4), npix, fp);
    fclose(fp);

    bmesh33_deinit(level_bmesh);
    bmesh33_dealloc(&level_bmesh);

    rtree_deinit(rtree);
    rtree_dealloc(&rtree);
  }

  mesh2_deinit(surface_mesh);
  mesh2_dealloc(&surface_mesh);

  bmesh33_deinit(bmesh);
  bmesh33_dealloc(&bmesh);
}
