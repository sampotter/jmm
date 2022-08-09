#include "eik3hh.h"

#include <assert.h>
#include <stdio.h>

#include "bmesh.h"
#include "eik3.h"
#include "eik3hh_branch.h"
#include "mesh2.h"
#include "rtree.h"

struct eik3hh {
  mesh3_s const *mesh;
  dbl c;
  dbl rfac;
  eik3hh_branch_s *root;
};

void eik3hh_alloc(eik3hh_s **hh) {
  *hh = malloc(sizeof(eik3hh_s));
}

static void init(eik3hh_s *hh, mesh3_s const *mesh, dbl c, dbl rfac) {
  hh->mesh = mesh;
  hh->c = c;
  hh->rfac = rfac;
  hh->root = NULL;
}

void eik3hh_init_with_pt_src(eik3hh_s *hh, mesh3_s const *mesh, dbl c,
                             dbl rfac, dbl3 const xsrc) {
  init(hh, mesh, c, rfac);

  eik3hh_branch_alloc(&hh->root);
  eik3hh_branch_init_pt_src(hh->root, hh, xsrc);
}

void eik3hh_deinit(eik3hh_s *hh) {
  hh->mesh = NULL;
  hh->c = NAN;
  hh->rfac = NAN;

  if (hh->root != NULL) {
    eik3hh_branch_deinit(hh->root, /* free_children = */ true);
    eik3hh_branch_dealloc(&hh->root);
  }
}

void eik3hh_dealloc(eik3hh_s **hh) {
  free(*hh);
  *hh = NULL;
}

mesh3_s const *eik3hh_get_mesh(eik3hh_s const *hh) {
  return hh->mesh;
}

dbl eik3hh_get_rfac(eik3hh_s const *hh) {
  return hh->rfac;
}

eik3hh_branch_s *eik3hh_get_root_branch(eik3hh_s *hh) {
  return hh->root;
}

dbl xfer(dbl x, size_t n) {
  dbl y = x;
  for (size_t i = 0; i < n; ++i)
    y = (3*y - 2*y*y)*y;
  return y;
}

void eik3hh_render_frames(eik3hh_s const *hh, camera_s const *camera,
                          dbl t0, dbl t1, dbl frame_rate,
                          bool verbose) {
  mesh2_s *surface_mesh = mesh3_get_surface_mesh(hh->mesh);

  eik3hh_branch_s *branch_direct = hh->root;
  bmesh33_s *bmesh_direct;
  bmesh33_alloc(&bmesh_direct);
  bmesh33_init_from_mesh3_and_jets(
    bmesh_direct, hh->mesh, eik3_get_jet_ptr(eik3hh_branch_get_eik(branch_direct)));

  array_s *children = eik3hh_branch_get_children(branch_direct);
  eik3hh_branch_s *branch_refl;
  array_get(children, 0, &branch_refl);
  bmesh33_s *bmesh_refl;
  bmesh33_alloc(&bmesh_refl);
  bmesh33_init_from_mesh3_and_jets(
    bmesh_refl, hh->mesh, eik3_get_jet_ptr(eik3hh_branch_get_eik(branch_refl)));

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

    dbl T = hh->c*t[i];

    bmesh33_s *level_bmesh_direct = bmesh33_restrict_to_level(bmesh_direct, T);
    rtree_insert_bmesh33(rtree, level_bmesh_direct);

    bmesh33_s *level_bmesh_refl = bmesh33_restrict_to_level(bmesh_refl, T);
    rtree_insert_bmesh33(rtree, level_bmesh_refl);

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
            if (bmesh33_cell->bmesh == level_bmesh_direct) {
              spread_interp = mesh3_linterp(
                hh->mesh, eik3hh_branch_get_spread(branch_direct), ray.org);
              org_interp = mesh3_linterp(
                hh->mesh, eik3hh_branch_get_org(branch_direct), ray.org);
            } else if (bmesh33_cell->bmesh == level_bmesh_refl) {
              spread_interp = mesh3_linterp(
                hh->mesh, eik3hh_branch_get_spread(branch_refl), ray.org);
              org_interp = mesh3_linterp(
                hh->mesh, eik3hh_branch_get_org(branch_refl), ray.org);
            } else {
              assert(false);
            }
            /* Convert the interpolated spreading factor to dB */
            dbl spread_dB = 20*log10(fmax(1e-16, spread_interp));
            /* Clamp and map the range [-60 dB, 0 dB] to [0, 1] for
             * use as a scaling factor */
            dbl spread_mapped = fmax(0, fmin(1, 1 - spread_dB/(-90)));
            alpha *= spread_mapped*xfer(org_interp, 2);
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

    bmesh33_deinit(level_bmesh_direct);
    bmesh33_dealloc(&level_bmesh_direct);

    bmesh33_deinit(level_bmesh_refl);
    bmesh33_dealloc(&level_bmesh_refl);

    rtree_deinit(rtree);
    rtree_dealloc(&rtree);
  }

  mesh2_deinit(surface_mesh);
  mesh2_dealloc(&surface_mesh);

  bmesh33_deinit(bmesh_direct);
  bmesh33_dealloc(&bmesh_direct);

  bmesh33_deinit(bmesh_refl);
  bmesh33_dealloc(&bmesh_refl);
}
