#include <argp.h>
#include <stdbool.h>
#include <stdlib.h>

#include <pthread.h>

#include <GLFW/glfw3.h>

#include "3d_wedge.h"

const char *argp_program_version = "3d_wedge 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "EXAMPLE";

static struct argp_option options[] = {
  {"verbose", 'v', 0, 0, "Produce verbose output", 0},
  {"visualize", 'z', 0, 0, "Make 3D visualization of solver running", 0},
  {"maxvol", 'a', "VOLUME", 0, "Maximum tetrahedron volume (default: 0.01)", 0},
  {0, 'n', "N", 0, "Wedge angle parameter (wedge angle is (2-n)pi, default: 0.25)", 0},
  {"width", 'w', "WIDTH", 0, "Domain bounding box width (default: 2)", 0},
  {"height", 'h', "HEIGHT", 0, "Domain bounding box height (default: 1)", 0},
  {"rho", 'R', "COEFFICIENT", 0, "Sound-hard reflection coefficient (default: 1)", 0},
  {"sp", 's', "DISTANCE", 0, "Distance of point source to edge (default: 1)", 0},
  {"phip", 'p', "ANGLE", 0, "Angle between source and n-face (default: -pi/4)", 0},
  {"rfac", 'r', "DISTANCE", 0, "Radius of initialization ball around point source (default: 0)", 0},
  {"omega", 'o', "FREQUENCY", 0, "Angular frequency of wave (default: 1000)", 0},
  {0}
};

typedef struct {
  jmm_3d_wedge_problem_s wedge;

  // these correspond to the arguments of jmm_3d_wedge_problem_solve
  double phip;
  double sp;
  double rfac;
  double omega;

  jmm_error_e error;
  char *error_msg;
} wkspc_s;

struct arguments {
  char *args[2];

  wkspc_s wkspc;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = state->input;

  wkspc_s *wkspc = &arguments->wkspc;
  jmm_3d_wedge_spec_s *spec = &arguments->wkspc.wedge.spec;

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
    wkspc->phip = atof(arg);
    break;
  case 's':
    wkspc->sp = atof(arg);
    break;
  case 'r':
    wkspc->rfac = atof(arg);
    break;
  case 'o':
    wkspc->omega = atof(arg);
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

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

static void *solve_wedge_problem(void *vp) {
  wkspc_s *w = (wkspc_s *)vp;

  jmm_3d_wedge_problem_s *wedge = &w->wedge;
  jmm_3d_wedge_spec_s *spec = &wedge->spec;

  w->error = JMM_ERROR_NONE;

  w->error = jmm_3d_wedge_problem_init(wedge, spec);
  if (w->error != JMM_ERROR_NONE) {
    w->error_msg = "Failed to initialize wedge problem";
    pthread_exit(NULL);
  }

  w->error = jmm_3d_wedge_problem_solve(wedge, w->sp, w->phip, w->rfac, w->omega);
  if (w->error != JMM_ERROR_NONE) {
    w->error_msg = "Failed to solve wedge problem";
    pthread_exit(NULL);
  }

  pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
  struct arguments arguments = {
    .wkspc = {
      .wedge = {
        .spec = {
          .verbose = false,
          .visualize = false,
          .maxvol = 0.01,
          .n = 0.25,
          .w = 2,
          .h = 1,
          .R = 1
        }
      },
      .sp = 1,
      .phip = -JMM_PI/4,
      .rfac = 0,
      .omega = 1000
    }
  };

  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  wkspc_s *wkspc = &arguments.wkspc;
  jmm_3d_wedge_problem_s *wedge = &wkspc->wedge;
  jmm_3d_wedge_spec_s *spec = &wedge->spec;

  if (spec->verbose) {
    printf("wedge parameters:\n");
    printf("* maxvol: %g\n", spec->maxvol);
    printf("* n: %g\n", spec->n);
    printf("* w: %g\n", spec->w);
    printf("* h: %g\n", spec->h);
    printf("* R: %g\n", spec->R);
    printf("solve parameters:\n");
    printf("* sp: %g\n", wkspc->sp);
    printf("* phip: %g\n", wkspc->phip);
    printf("* rfac: %g\n", wkspc->rfac);
    printf("* omega: %g\n", wkspc->omega);
  }

  GLFWwindow *window = NULL;

  pthread_t thread;

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  pthread_create(&thread, &attr, solve_wedge_problem, wkspc);

  if (spec->visualize) {
    glfwInit();
    window = glfwCreateWindow(1280, 720, "Viz", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwDestroyWindow(window);
    glfwTerminate();
  }

  pthread_join(thread, NULL);

  if (wkspc->error != JMM_ERROR_NONE) {
    printf("%s: %s\n", jmm_error_to_string(wkspc->error), wkspc->error_msg);
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}
