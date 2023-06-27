#include <assert.h>
#include <argp.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <jmm/array.h>
#include <jmm/bmesh.h>
#include <jmm/camera.h>
#include <jmm/mesh2.h>
#include <jmm/mesh3.h>
#include <jmm/rtree.h>
#include <jmm/util.h>

const char *argp_program_version = "render 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "BRANCH_PATH ...";

enum long_opts {
  LONG_OPT_T0 = 777,
  LONG_OPT_T1,
  LONG_OPT_FRAMES_PER_METER,
  LONG_OPT_VERTS_PATH,
  LONG_OPT_CELLS_PATH,
};

static struct argp_option options[] = {
  {"verbose", 'v', 0, OPTION_ARG_OPTIONAL, "Verbose output", 0},
  {"speed", 'c', "SPEED", OPTION_ARG_OPTIONAL,
   "Speed of sound in m/s (default: 340.3 m/s)", 0},
  {"t0", LONG_OPT_T0, "T0", OPTION_ARG_OPTIONAL,
   "Initial time point for rendering", 0},
  {"t1", LONG_OPT_T1, "T1", OPTION_ARG_OPTIONAL,
   "Final time point for rendering", 0},
  {"frames_per_meter", LONG_OPT_FRAMES_PER_METER, "FRAMES_PER_METER",
   OPTION_ARG_OPTIONAL,
   "Rendering frame rate (default: 23.98 s^-1)", 0},
  {"verts_path", LONG_OPT_VERTS_PATH, "VERTS_PATH", OPTION_ARG_OPTIONAL,
   "Path to binary file with vertices (default: \"verts.bin\")", 0},
  {"cells_path", LONG_OPT_CELLS_PATH, "CELLS_PATH", OPTION_ARG_OPTIONAL,
   "Path to binary file with vertices (default: \"cells.bin\")", 0},
  {0}
};

typedef struct render_spec {
  bool verbose;
  dbl c;
  dbl t0;
  dbl t1;
  dbl frames_per_meter;
  char *verts_path;
  char *cells_path;
  array_s *branch_path_arr;
} render_spec_s;

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  render_spec_s *spec = state->input;
  size_t n;
  char *branch_path;
  switch (key) {
  case 'v':
    spec->verbose = true;
    break;
  case 'c':
    spec->c = atof(arg);
    break;
  case LONG_OPT_T0:
    spec->t0 = atof(arg);
    break;
  case LONG_OPT_T1:
    spec->t1 = atof(arg);
    break;
  case LONG_OPT_FRAMES_PER_METER:
    spec->frames_per_meter = atof(arg);
    break;
  case LONG_OPT_VERTS_PATH:
    n = strlen(arg);
    spec->verts_path = malloc(n + 1);
    strncpy(spec->verts_path, arg, n + 1);
    break;
  case LONG_OPT_CELLS_PATH:
    n = strlen(arg);
    spec->cells_path = malloc(n + 1);
    strncpy(spec->cells_path, arg, n + 1);
    break;
  case ARGP_KEY_ARG:
    n = strlen(arg);
    branch_path = malloc(n + 1);
    strncpy(branch_path, arg, n + 1);
    array_append(spec->branch_path_arr, &branch_path);
    break;
  case ARGP_KEY_END:
    if (array_is_empty(spec->branch_path_arr))
      argp_usage(state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return EXIT_SUCCESS;
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

int main(int argc, char *argv[]) {
  int code = EXIT_SUCCESS;

  render_spec_s spec = {
    .verbose = false,
    .c = 340.3, /* m/s */
    .t0 = NAN,
    .t1 = NAN,
    .frames_per_meter = 23.98,
    .verts_path = NULL,
    .cells_path = NULL
  };

  array_alloc(&spec.branch_path_arr);
  array_init(spec.branch_path_arr, sizeof(char *), ARRAY_DEFAULT_CAPACITY);

  argp_parse(&argp, argc, argv, 0, 0, &spec);

  if (spec.verbose) {
    printf("Render info:\n");
    printf("- t0 = %g\n", spec.t0);
    printf("- t1 = %g\n", spec.t1);
    printf("- frames_per_meter = %g\n", spec.frames_per_meter);
    printf("- verts_path = %s\n", spec.verts_path);
    printf("- cells_path = %s\n", spec.cells_path);
    printf("- branch_paths:\n");
    for (size_t i = 0; i < array_size(spec.branch_path_arr); ++i) {
      char *branch_path;
      array_get(spec.branch_path_arr, i, &branch_path);
      printf("  + %s\n", branch_path);
    }
  }

  mesh3_data_s data;
  mesh3_data_init_from_bin(
    &data,
    spec.verts_path ? spec.verts_path : "verts.bin",
    spec.cells_path ? spec.cells_path : "cells.bin");

  dbl eps = 1e-5;

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, &eps);

  array_s *bmesh_arr;
  array_alloc(&bmesh_arr);
  array_init(bmesh_arr, sizeof(bmesh33_s *), ARRAY_DEFAULT_CAPACITY);

  array_s *spread_arr;
  array_alloc(&spread_arr);
  array_init(spread_arr, sizeof(dbl *), ARRAY_DEFAULT_CAPACITY);

  array_s *org_arr;
  array_alloc(&org_arr);
  array_init(org_arr, sizeof(dbl *), ARRAY_DEFAULT_CAPACITY);

  for (size_t i = 0; i < array_size(spec.branch_path_arr); ++i) {
    char path[256];

    char *branch_path;
    array_get(spec.branch_path_arr, i, &branch_path);

    sprintf(path, "%s/jet.bin", branch_path);

    FILE *fp = fopen(path, "r");
    jet31t *jet = malloc(mesh3_nverts(mesh)*sizeof(jet31t));
    fread(jet, sizeof(jet31t), mesh3_nverts(mesh), fp);
    fclose(fp);

    bmesh33_s *bmesh;
    bmesh33_alloc(&bmesh);
    bmesh33_init_from_mesh3_and_jets(bmesh, mesh, jet);

    array_append(bmesh_arr, &bmesh);

    free(jet);

    sprintf(path, "%s/spread.bin", branch_path);
    fp = fopen(path, "r");
    dbl *spread = malloc(mesh3_nverts(mesh)*sizeof(dbl));
    fread(spread, sizeof(dbl), mesh3_nverts(mesh), fp);
    fclose(fp);
    array_append(spread_arr, &spread);

    sprintf(path, "%s/org.bin", branch_path);
    fp = fopen(path, "r");
    dbl *org = malloc(mesh3_nverts(mesh)*sizeof(dbl));
    fread(org, sizeof(dbl), mesh3_nverts(mesh), fp);
    fclose(fp);
    array_append(org_arr, &org);
  }

  mesh2_s *surface_mesh = mesh3_get_surface_mesh(mesh);

  /* TODO: read camera from file */
  camera_s camera = {
    // .type = CAMERA_TYPE_ORTHOGRAPHIC,
    .type = CAMERA_TYPE_PERSPECTIVE,
    .pos = {-6, 6, 15},
    .look = {0, 0, -1},
    .left = {-1, 0, 0},
    .up = {0, 1, 0},
    // .width = 22.0,
    // .height = 22.0/4,
    .fovy = 30,
    .aspect = 1,
    .dim = {512, 512}
  };

  dbl tau0 = spec.c*spec.t0;
  dbl tau1 = spec.c*spec.t1;

  size_t num_frames = floor(spec.frames_per_meter*(tau1 - tau0));
  if (spec.verbose)
    printf("Rendering %lu frames\n", num_frames);

  dbl *tau = malloc(num_frames*sizeof(dbl));
  for (size_t i = 0; i < num_frames; ++i)
    tau[i] = tau0 + i/spec.frames_per_meter;

  for (size_t i = 0; i < num_frames; ++i) {
    if (spec.verbose)
      printf("frame %lu/%lu (tau = %g m)\n", i + 1, num_frames, tau[i]);

    rtree_s *rtree;
    rtree_alloc(&rtree);
    rtree_init(rtree, 16, RTREE_SPLIT_STRATEGY_SURFACE_AREA);

    rtree_insert_mesh2(rtree, surface_mesh);

    array_s *level_bmesh_arr;
    array_alloc(&level_bmesh_arr);
    array_init(level_bmesh_arr, sizeof(bmesh33_s *), ARRAY_DEFAULT_CAPACITY);

    for (size_t j = 0; j < array_size(bmesh_arr); ++j) {
      bmesh33_s *bmesh;
      array_get(bmesh_arr, j, &bmesh);
      bmesh33_s *level_bmesh = bmesh33_restrict_to_level(bmesh, tau[i]);
      rtree_insert_bmesh33(rtree, level_bmesh);
      array_append(level_bmesh_arr, &level_bmesh);
    }

    rtree_build(rtree);

    size_t npix = camera.dim[0]*camera.dim[1];

    dbl4 *img = malloc(npix*sizeof(dbl4));

    dbl3 surf_rgb = {0.54, 0.54, 0.54};
    dbl3 eik_rgb = {1.0, 1.0, 1.0};

    dbl surf_alpha = 0.5;
    dbl eik_alpha = 1;

    for (size_t i = 0, l = 0; i < camera.dim[0]; ++i) {
      for (size_t j = 0; j < camera.dim[1]; ++j, ++l) {
        ray3 ray = camera_get_ray_for_index(&camera, i, j);

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

            /* Figure out which level mesh this cell belongs to :-( */
            size_t k;
            for (k = 0; k < array_size(level_bmesh_arr); ++k) {
              bmesh33_s *level_bmesh;
              array_get(level_bmesh_arr, k, &level_bmesh);
              if (bmesh33_cell->bmesh == level_bmesh)
                break;
            }
            assert(k < array_size(level_bmesh_arr));

            dbl const *spread;
            array_get(spread_arr, k, &spread);
            spread_interp = mesh3_linterp(mesh, spread, ray.org);

            dbl const *org;
            array_get(org_arr, k, &org);
            org_interp = mesh3_linterp(mesh, org, ray.org);

            // Convert the interpolated spreading factor to dB
            dbl spread_dB = 20*log10(fmax(1e-16, spread_interp));
            /* Clamp and map the range [-60 dB, 0 dB] to [0, 1] for
             * use as a scaling factor */
            dbl spread_mapped = fmax(0, fmin(1, 1 - spread_dB/(-90)));
            alpha *= spread_mapped*squash(org_interp, 2);
          }

          /* Get the surface normal and dot it with the eye vector for
           * Lambertian shading */
          if (robj_type == ROBJ_MESH2_TRI) {
            mesh2_tri_s const *mesh2_tri = robj_data;
            mesh2_get_unit_surface_normal(surface_mesh, mesh2_tri->l, n);
          } else if (robj_type == ROBJ_BMESH33_CELL) {
            bmesh33_cell_s const *bmesh33_cell = robj_data;
            bmesh33_cell_Df(bmesh33_cell, ray.org, n);
            dbl3_normalize(n);
          } else {
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

    for (size_t i = 0; i < array_size(level_bmesh_arr); ++i) {
      bmesh33_s *level_bmesh;
      array_get(level_bmesh_arr, i, &level_bmesh);
      bmesh33_deinit(level_bmesh);
      bmesh33_dealloc(&level_bmesh);
    }
    array_deinit(level_bmesh_arr);
    array_dealloc(&level_bmesh_arr);

    rtree_deinit(rtree);
    rtree_dealloc(&rtree);
  }

  mesh2_deinit(surface_mesh);
  mesh2_dealloc(&surface_mesh);

  for (size_t i = 0; i < array_size(bmesh_arr); ++i) {
    bmesh33_s *bmesh;
    array_get(bmesh_arr, i, &bmesh);
    bmesh33_deinit(bmesh);
    bmesh33_dealloc(&bmesh);
  }
  array_deinit(bmesh_arr);
  array_dealloc(&bmesh_arr);

  if (spec.verts_path != NULL) free(spec.verts_path);
  if (spec.cells_path != NULL) free(spec.cells_path);

  for (size_t i = 0; i < array_size(spec.branch_path_arr); ++i) {
    char *branch_path;
    array_get(spec.branch_path_arr, i, &branch_path);
    free(branch_path);
  }
  array_deinit(spec.branch_path_arr);
  array_dealloc(&spec.branch_path_arr);

  free(tau);

  return code;
}
