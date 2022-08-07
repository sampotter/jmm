#include <argp.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <bmesh.h>
#include <camera.h>
#include <eik3.h>
#include <eik3hh.h>
#include <eik3hh_branch.h>
#include <util.h>

const char *argp_program_version = "simple_building 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "EXAMPLE";

enum long_opts {
  LONG_OPT_RENDER = 777,
  LONG_OPT_T0,
  LONG_OPT_T1,
  LONG_OPT_FRAMERATE
};

static struct argp_option options[] = {
  {"verbose", 'v', 0, OPTION_ARG_OPTIONAL, "Produce verbose output", 0},
  {"maxvol", 'a', "VOLUME", OPTION_ARG_OPTIONAL,
   "Maximum tetrahedron volume (default: 0.01)", 0},
  {"xsrc", 'x', "POINT", OPTION_ARG_OPTIONAL,
   "Point source location (format: \"x,y,z\")", 0},
  {"Rcoef", 'R', "RCOEF", OPTION_ARG_OPTIONAL,
   "Sound-hard reflection coefficient (default: 1)", 0},
  {"rfac", 'r', "RADIUS", OPTION_ARG_OPTIONAL,
   "Factoring radius (default: 0.1)", 0},
  {"speed", 'c', "SPEED", OPTION_ARG_OPTIONAL,
   "Speed of sound in m/s (default: 340.3 m/s)", 0},
  {"freq", 'f', "FREQUENCY", OPTION_ARG_OPTIONAL,
   "Angular frequency of wave (default: 1000)", 0},
  {"omega", 0, "FREQUENCY", OPTION_ALIAS, 0, 0},
  {"render", LONG_OPT_RENDER, 0, OPTION_ARG_OPTIONAL,
   "Render video of wavefront", 0},
  {"t0", LONG_OPT_T0, "T0", OPTION_ARG_OPTIONAL,
   "Initial time point for rendering", 0},
  {"t1", LONG_OPT_T1, "T1", OPTION_ARG_OPTIONAL,
   "Final time point for rendering", 0},
  {"framerate", LONG_OPT_FRAMERATE, "FRAMERATE", OPTION_ARG_OPTIONAL,
   "Rendering frame rate (default: 23.98 s^-1)", 0},
  {0}
};

typedef struct problem_spec {
  bool verbose;
  dbl maxvol;
  dbl3 xsrc;
  dbl rfac;
  dbl c;
  dbl omega;
  bool render;
  dbl t0;
  dbl t1;
  dbl frame_rate;
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
  case 'c':
    spec->c = atof(arg);
    break;
  case 'f':
    spec->omega = atof(arg);
    break;
  case LONG_OPT_RENDER:
    spec->render = true;
    break;
  case LONG_OPT_T0:
    spec->t0 = atof(arg);
    break;
  case LONG_OPT_T1:
    spec->t1 = atof(arg);
    break;
  case LONG_OPT_FRAMERATE:
    spec->frame_rate = atof(arg);
    break;
  case ARGP_KEY_ARG:
    size_t n = strlen(arg);
    spec->off_path = malloc(n + 1);
    strncpy(spec->off_path, arg, n + 1);
    break;
  case ARGP_KEY_END:
    if (spec->off_path == NULL)
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

  problem_spec_s spec = {
    .verbose = false,
    .maxvol = 0.01,
    .xsrc = {NAN, NAN, NAN},
    .rfac = 0.1,
    .c = 340.3,
    .omega = 100,
    .render = false,
    .t0 = NAN,
    .t1 = NAN,
    .frame_rate = 23.98,
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

  mesh3_data_s data;
  mesh3_data_from_off_file(&data, spec.off_path, spec.maxvol, spec.verbose);
  mesh3_data_insert_vert(&data, spec.xsrc, EPS);

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, NULL);

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
  eik3hh_init_with_pt_src(hh, mesh, spec.c, spec.rfac, spec.xsrc);

  eik3hh_branch_s *root_branch = eik3hh_get_root_branch(hh);
  eik3_s *eik_direct = eik3hh_branch_get_eik(root_branch);

  printf("Set up direct eikonal problem:\n");
  printf("- number of points inside factoring radius: %lu\n",
         eik3_num_valid(eik_direct));
  printf("- number of TRIAL points adjacent to factored ball: %lu\n",
         eik3_num_trial(eik_direct));
  if (eik3_num_trial(eik_direct) == 0) {
    printf("ERROR: didn't insert any TRIAL points!\n");
    exit(EXIT_FAILURE);
  }

  toc();
  eik3hh_branch_solve(root_branch);
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
  eik3hh_branch_dump_xy_slice(
    root_branch, &mapping, FIELD_T, "T_slice.bin");
  printf("Saved xy-slice (z = %g) of the direct eikonal [%1.2gs]\n",
         spec.xsrc[2], toc());

  /* Save spreading factor slice to disk */
  toc();
  eik3hh_branch_dump_xy_slice(
    root_branch, &mapping, FIELD_SPREADING, "spread_slice.bin");
  printf("Saved xy-slice (z = %g) of the direct spreading [%1.2gs]\n",
         spec.xsrc[2], toc());

  /* Save origin slice to disk */
  toc();
  eik3hh_branch_dump_xy_slice(
    root_branch, &mapping, FIELD_ORIGIN, "origin_slice.bin");
  printf("Saved xy-slice (z = %g) of the direct origin [%1.2gs]\n",
         spec.xsrc[2], toc());

  /** Render frames to disk */

  /* TODO: read camera from file */
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
  eik3hh_render_frames(hh, &camera, spec.t0, spec.t1, spec.frame_rate, true);

cleanup:

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
