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

void eik3hh_render_frames(eik3hh_s const *hh, camera_s const *camera,
                          dbl t0, dbl t1, dbl frame_rate,
                          bool verbose) {
  mesh2_s *surface_mesh = mesh3_get_surface_mesh(hh->mesh);

  eik3_s const *eik_direct = eik3hh_branch_get_eik(hh->root);
  jet31t const *jet_direct = eik3_get_jet_ptr(eik_direct);

  bmesh33_s *bmesh;
  bmesh33_alloc(&bmesh);
  bmesh33_init_from_mesh3_and_jets(bmesh, hh->mesh, jet_direct);

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
