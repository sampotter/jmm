#include "eik3hh.h"

#include <assert.h>
#include <stdio.h>

#include "bmesh.h"
#include "eik3.h"
#include "eik3_transport.h"
#include "mat.h"
#include "mesh2.h"
#include "rtree.h"

struct eik3hh {
  dbl c;

  size_t lsrc;
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

void eik3hh_init_pt_src(eik3hh_s *hh, mesh3_s const *mesh,
                        dbl3 const xsrc, dbl rfac, dbl c) {
  init(hh, mesh);

  hh->c = c;

  hh->lsrc = mesh3_get_vert_index(mesh, xsrc);
  dbl3_copy(xsrc, hh->xsrc);
  hh->rfac = rfac;

  eik3_add_pt_src_bcs(hh->eik, hh->xsrc, hh->rfac);
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
    if (dbl3_all_nan(par.b) || par3_has_active_parent(&par, hh->lsrc)) {
      dbl3 x;
      mesh3_copy_vert(hh->mesh, l, x);
      hh->spread[l] = 1/dbl3_dist(x, hh->xsrc);
    }
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
    if (!isnan(spread[lhat]))
      continue;

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

void eik3hh_solve(eik3hh_s *hh) {
  eik3_solve(hh->eik);

  init_D2T_pt_src(hh);
  approx_D2T(hh);

  init_spread_pt_src(hh);
  prop_spread(hh);

  eik3_get_org(hh->eik, hh->origin);
}

size_t eik3hh_num_trial(eik3hh_s const *hh) {
  return eik3_num_trial(hh->eik);
}

size_t eik3hh_num_valid(eik3hh_s const *hh) {
  return eik3_num_valid(hh->eik);
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

void eik3hh_render_frames(eik3hh_s const *hh, camera_s const *camera,
                          dbl t0, dbl t1, dbl frame_rate,
                          bool verbose) {
  bmesh33_s *bmesh;
  bmesh33_alloc(&bmesh);
  bmesh33_init_from_mesh3_and_jets(bmesh, hh->mesh, eik3hh_get_jet_ptr(hh));

  mesh2_s *surface_mesh = mesh3_get_surface_mesh(hh->mesh);

  size_t num_frames = floor(frame_rate*(t1 - t0));
  if (verbose)
    printf("rendering %lu frames\n", num_frames);

  dbl *t = malloc(num_frames*sizeof(dbl));
  for (size_t i = 0; i < num_frames; ++i)
    t[i] = t0 + i/frame_rate;

  for (size_t i = 0; i < num_frames; ++i) {
    if (verbose)
      printf("frame %lu/%lu (t = %g s)\n", i + 1, num_frames, t[i]);

    rtree_s *rtree;
    rtree_alloc(&rtree);
    rtree_init(rtree, 16, RTREE_SPLIT_STRATEGY_SURFACE_AREA);

    rtree_insert_mesh2(rtree, surface_mesh);

    bmesh33_s *level_bmesh = bmesh33_restrict_to_level(bmesh, hh->c*t[i]);
    rtree_insert_bmesh33(rtree, level_bmesh);

    rtree_build(rtree);

    size_t npix = camera->dim[0]*camera->dim[1];

    dbl4 *img = malloc(npix*sizeof(dbl4));

    dbl3 surf_rgb = {1.0, 1.0, 1.0};
    dbl3 eik_rgb = {1.0, 1.0, 1.0};

    dbl surf_alpha = 0.5;
    dbl eik_alpha = 0.95;

    for (size_t i = 0, l = 0; i < camera->dim[0]; ++i) {
      for (size_t j = 0; j < camera->dim[1]; ++j, ++l) {
        ray3 ray = camera_get_ray_for_index(camera, i, j);

        isect isect;
        rtree_intersect(rtree, &ray, &isect, NULL);

        img[l][0] = 0;
        img[l][1] = 0;
        img[l][2] = 0;
        img[l][3] = isfinite(isect.t) ? 1 : 0;

        dbl alpha = 1, scale;
        dbl const *rgb = NULL;
        dbl3 n;

        while (isfinite(isect.t)) {
          robj_type_e robj_type = robj_get_type(isect.obj);
          void const *robj_data = robj_get_data(isect.obj);

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
          scale = fabs(dbl3_dot(n, ray.dir));

          /* ... and accumulate */
          dbl3_saxpy_inplace(scale*alpha, rgb, img[l]);

          /* Advance the start of the ray and keep tracing */
          rtree_intersect(rtree, &ray, &isect, isect.obj);
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
