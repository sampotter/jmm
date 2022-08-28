#include <argp.h>
#include <stdbool.h>
#include <stdlib.h>

#include "3d_wedge.h"

const char *argp_program_version = "3d_wedge 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "EXAMPLE";

static struct argp_option options[] = {
  {"verbose", 'v', 0, 0, "Produce verbose output", 0},
  {"visualize", 'z', 0, 0, "Make 3D visualization of solver running", 0},
  {"maxvol", 'a', "VOLUME", 0, "Maximum tetrahedron volume (default: 0.01)", 0},
  {0, 'n', "N", 0, "Wedge angle parameter (wedge angle is (2-n)pi, default: 1.75)", 0},
  {"width", 'w', "WIDTH", 0, "Domain bounding box width (default: 2)", 0},
  {"height", 'h', "HEIGHT", 0, "Domain bounding box height (default: 1)", 0},
  {"rho", 'R', "COEFFICIENT", 0, "Sound-hard reflection coefficient (default: 1)", 0},
  {"sp", 's', "DISTANCE", 0, "Distance of point source to edge (default: 1)", 0},
  {"phip", 'p', "ANGLE", 0, "Angle between source and o-face (default: pi/4)", 0},
  {"rfac", 'r', "DISTANCE", 0, "Radius of initialization ball around point source (default: 0)", 0},
  {"omega", 'o', "FREQUENCY", 0, "Angular frequency of wave (default: 1000)", 0},
  {0}
};

struct arguments {
  char *args[2];
  jmm_3d_wedge_spec_s *spec;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = state->input;

  jmm_3d_wedge_spec_s *spec = arguments->spec;

  switch (key) {
  case 'v':
    spec->verbose = true;
    break;
  case 'z':
    spec->visualize = true;
    break;
  case 'a':
    spec->maxvol = atof(arg); // TODO: switch to strtod
    break;
  case 'n':
    spec->n = atof(arg);
    break;
  case 'w':
    spec->w = atof(arg);
    break;
  case 'h':
    spec->h = atof(arg);
    break;
  case 'R':
    spec->R = atof(arg);
    break;
  case 'p':
    spec->phip = atof(arg);
    break;
  case 's':
    spec->sp = atof(arg);
    break;
  case 'r':
    spec->rfac = atof(arg);
    break;
  case 'o':
    spec->omega = atof(arg);
    break;
  case ARGP_KEY_ARG:
    arguments->args[state->arg_num] = arg;
    break;
  case ARGP_KEY_END:
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }

  return EXIT_SUCCESS;
}

static void print_params(jmm_3d_wedge_spec_s const *spec) {
  printf("wedge parameters:\n");
  printf("* maxvol: %g\n", spec->maxvol);
  printf("* sp: %g\n", spec->sp);
  printf("* phip: %g\n", spec->phip);
  printf("* rfac: %g\n", spec->rfac);
  printf("* omega: %g\n", spec->omega);
  printf("* n: %g\n", spec->n);
  printf("* w: %g\n", spec->w);
  printf("* h: %g\n", spec->h);
  printf("* R: %g\n", spec->R);
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

int main(int argc, char *argv[]) {
  jmm_3d_wedge_spec_s spec = {
    .verbose = false,
    .visualize = false,
    .maxvol = 0.01,
    .sp = 1,
    .phip = -JMM_PI/4,
    .rfac = 0,
    .omega = 1000,
    .n = 1.75,
    .w = 2,
    .h = 1,
    .R = 1
  };

  struct arguments arguments = {.spec = &spec};

  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  if (spec.verbose)
    print_params(&spec);

  jmm_3d_wedge_problem_s wedge;
  jmm_error_e error = jmm_3d_wedge_problem_init(&wedge, &spec);
  if (error != JMM_ERROR_NONE) {
    printf("ERROR: Failed to initialize wedge problem\n");
    goto cleanup;
  }

  error = jmm_3d_wedge_problem_solve(&wedge);
  if (error != JMM_ERROR_NONE) {
    printf("ERROR: Failed to solve wedge problem\n");
    goto cleanup;
  }

  jmm_3d_wedge_problem_dump(&wedge, ".", true, true, true);

  size_t ngrid = 256 + 1;
  grid2_s img_grid = {
    .shape = {ngrid, ngrid},
    .xymin = {-wedge.spec.w/2, -wedge.spec.w/2},
    .h = wedge.spec.w/(ngrid - 1),
    .order = ORDER_ROW_MAJOR
  };
  jmm_3d_wedge_problem_save_slice_plots(&wedge,".",true,true,true,&img_grid);

cleanup:
  exit(error == JMM_ERROR_NONE ? EXIT_SUCCESS : EXIT_FAILURE);
}
