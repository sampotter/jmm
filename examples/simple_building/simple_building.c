#include <argp.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <bmesh.h>
#include <camera.h>
#include <eik3.h>
#include <eik3hh.h>
#include <mesh2.h>
#include <rtree.h>
#include <util.h>

const char *argp_program_version = "simple_building 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "EXAMPLE";

static struct argp_option options[] = {
  {"verbose", 'v', 0, OPTION_ARG_OPTIONAL, "Produce verbose output", 0},
  {"maxvol", 'a', "VOLUME", OPTION_ARG_OPTIONAL,
   "Maximum tetrahedron volume (default: 0.01)", 0},
  {"xsrc", 'x', "POINT", OPTION_ARG_OPTIONAL,
   "Sound-hard reflection coefficient (default: 1)", 0},
  {"rfac", 'r', "RADIUS", OPTION_ARG_OPTIONAL,
   "Factoring radius (default: 0.1)", 0},
  {"freq", 'f', "FREQUENCY", OPTION_ARG_OPTIONAL,
   "Angular frequency of wave (default: 1000)", 0},
  {"omega", 0, "FREQUENCY", OPTION_ALIAS, 0, 0},
  {0}
};

typedef struct problem_spec {
  bool verbose;
  dbl maxvol;
  dbl3 xsrc;
  dbl rfac;
  dbl omega;
  char *off_path;
} problem_spec_s;

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  problem_spec_s *spec = state->input;
  switch (key) {
  case 'v':
    spec->verbose = true;
    break;
  case 'a':
    spec->maxvol = atof(arg);
    break;
  case 'x':
    char *str = arg;
    spec->xsrc[0] = atof(strtok_r(str, ",", &str));
    spec->xsrc[1] = atof(strtok_r(str, ",", &str));
    spec->xsrc[2] = atof(strtok_r(str, ",", &str));
    break;
  case 'r':
    spec->rfac = atof(arg);
    break;
  case 'f':
    spec->omega = atof(arg);
    break;
  case ARGP_KEY_ARG:
    size_t n = strlen(arg);
    spec->off_path = malloc(n + 1);
    strncpy(spec->off_path, arg, n + 1);
    break;
  case ARGP_KEY_END:
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return EXIT_SUCCESS;
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

int main(int argc, char *argv[]) {
  int code = EXIT_SUCCESS;

  problem_spec_s spec = {
    .verbose = false,
    .maxvol = 0.01,
    .xsrc = {NAN, NAN, NAN},
    .rfac = 0.1,
    .omega = 100,
    .off_path = NULL
  };

  argp_parse(&argp, argc, argv, 0, 0, &spec);

  if (spec.verbose) {
    printf("Problem info:\n");
    printf("- maxvol: %g m^3\n", spec.maxvol);
    printf("- xsrc: (%g, %g, %g) m\n", spec.xsrc[0], spec.xsrc[1], spec.xsrc[2]);
    printf("- rfac: %g m\n", spec.rfac);
    printf("- omega: %g Hz\n", spec.omega);
    printf("- off_path: %s\n", spec.off_path);
  }

  toc();

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init_from_off_file(mesh, spec.off_path, spec.maxvol);

  rect3 bbox;
  mesh3_get_bbox(mesh, &bbox);

  printf("Built tetrahedron mesh [%1.2gs]:\n", toc());
  printf("- %lu vertices\n", mesh3_nverts(mesh));
  printf("- %lu cells.\n", mesh3_ncells(mesh));
  printf("- bounding box: [%g, %g] x [%g, %g] x [%g, %g]\n",
         bbox.min[0], bbox.max[0],
         bbox.min[1], bbox.max[1],
         bbox.min[2], bbox.max[2]);
  printf("- h: %g\n", mesh3_get_mean_edge_length(mesh));

  if (!mesh3_contains_ball(mesh, spec.xsrc, spec.rfac)) {
    fprintf(stderr, "ERROR: mesh doesn't fully contain factoring ball\n");
    code = EXIT_FAILURE;
    goto cleanup;
  }

  eik3hh_s *hh;
  eik3hh_alloc(&hh);
  eik3hh_init_pt_src(hh, mesh, spec.xsrc, spec.rfac);

  printf("Set up direct eikonal problem:\n");
  printf("- num. BC points: %lu\n", eik3hh_num_bc(hh));

  toc();
  eik3hh_solve(hh);
  printf("Computed direct eikonal [%1.2gs]\n", toc());

  /* Set up the image grid for making 2D plots */
  size_t ngrid = 1024 + 1;
  grid2_s img_grid = {
    .shape = {ngrid, ngrid},
    .xymin = {bbox.min[0], bbox.min[1]},
    .h = (bbox.max[0] - bbox.min[0])/(ngrid - 1),
    .order = ORDER_ROW_MAJOR
  };
  grid2_save(&img_grid, "img_grid.txt");

  /* Set up the mapping from the image grid to the tetrahedron mesh. */
  toc();
  grid2_to_mesh3_mapping_s mapping;
  grid2_to_mesh3_mapping_init_xy(&mapping, &img_grid, mesh, spec.xsrc[2]);
  printf("Set up grid-to-mesh mapping [%1.2gs]\n", toc());

  /* Save T slice to disk */
  toc();
  eik3hh_dump_xy_slice(hh, &mapping, FIELD_T, "T_slice.bin");
  printf("Saved xy-slice (z = %g) of the direct eikonal [%1.2gs]\n",
         spec.xsrc[2], toc());

  /* Save spreading factor slice to disk */
  toc();
  eik3hh_dump_xy_slice(hh, &mapping, FIELD_SPREADING, "spread_slice.bin");
  printf("Saved xy-slice (z = %g) of the direct spreading [%1.2gs]\n",
         spec.xsrc[2], toc());

  /* Save origin slice to disk */
  toc();
  eik3hh_dump_xy_slice(hh, &mapping, FIELD_ORIGIN, "origin_slice.bin");
  printf("Saved xy-slice (z = %g) of the direct origin [%1.2gs]\n",
         spec.xsrc[2], toc());

  /** <RAYTRACING> */

  camera_s camera = {
    // .type = CAMERA_TYPE_ORTHOGRAPHIC,
    .type = CAMERA_TYPE_PERSPECTIVE,
    .pos = {0, 0, 42},
    .look = {0, 0, -1},
    .left = {-1, 0, 0},
    .up = {0, 1, 0},
    // .width = 22.0,
    // .height = 22.0/4,
    .fovy = 30,
    .aspect = 1,
    .dim = {1024, 1024}
  };

  bmesh33_s *bmesh;
  bmesh33_alloc(&bmesh);
  bmesh33_init_from_mesh3_and_jets(bmesh, mesh, eik3hh_get_jet_ptr(hh));

  mesh2_s *surface_mesh = mesh3_get_surface_mesh(mesh);

  dbl T_max = eik3_get_max_T(eik3hh_get_eik_ptr(hh));

  dbl T0 = 0.2*T_max, T1 = 0.8*T_max;

  size_t num_T = 12;
  dbl *T = malloc(num_T*sizeof(dbl));
  for (size_t i = 0; i < num_T; ++i) {
    dbl t = i/(dbl)(num_T - 1);
    T[i] = (1 - t)*T0 + t*T1;
  }

  for (size_t i = 0; i < num_T; ++i) {
    printf("frame %lu/%lu (T = %g)\n", i, num_T - 1, T[i]);

    rtree_s *rtree;
    rtree_alloc(&rtree);
    rtree_init(rtree, 16, RTREE_SPLIT_STRATEGY_SURFACE_AREA);

    rtree_insert_mesh2(rtree, surface_mesh);

    bmesh33_s *level_bmesh = bmesh33_restrict_to_level(bmesh, T[i]);
    rtree_insert_bmesh33(rtree, level_bmesh);

    rtree_build(rtree);

    size_t npix = camera.dim[0]*camera.dim[1];

    dbl4 *img = malloc(npix*sizeof(dbl4));

    dbl3 surf_rgb = {1.0, 1.0, 1.0};
    dbl3 eik_rgb = {1.0, 1.0, 1.0};

    dbl surf_alpha = 0.5;
    dbl eik_alpha = 0.95;

    for (size_t i = 0, l = 0; i < camera.dim[0]; ++i) {
      for (size_t j = 0; j < camera.dim[1]; ++j, ++l) {
        ray3 ray = camera_get_ray_for_index(&camera, i, j);

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

  /** </RAYTRACING> */

cleanup:

  mesh2_deinit(surface_mesh);
  mesh2_dealloc(&surface_mesh);

  bmesh33_deinit(bmesh);
  bmesh33_dealloc(&bmesh);

  grid2_to_mesh3_mapping_deinit(&mapping);

  eik3hh_deinit(hh);
  eik3hh_dealloc(&hh);

  /* Free tetrahedron mesh */
  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  /* Free stuff related to problem specification */
  if (spec.off_path)
    free(spec.off_path);

  return code;
}
