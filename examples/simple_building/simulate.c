#include <argp.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <jmm/eik3.h>
#include <jmm/eik3hh.h>
#include <jmm/eik3hh_branch.h>
#include <jmm/util.h>

const char *argp_program_version = "simulate 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "EXAMPLE";

static struct argp_option options[] = {
  {"verbose", 'v', 0, OPTION_ARG_OPTIONAL, "Produce verbose output", 0},
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
  {0}
};

typedef struct problem_spec {
  bool verbose;
  dbl3 xsrc;
  dbl rfac;
  dbl c;
  dbl omega;
} problem_spec_s;

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  problem_spec_s *spec = state->input;
  switch (key) {
  case 'v':
    spec->verbose = true;
    break;
  case 'x': {
    char *str = arg;
    spec->xsrc[0] = atof(strtok_r(str, ",", &str));
    spec->xsrc[1] = atof(strtok_r(str, ",", &str));
    spec->xsrc[2] = atof(strtok_r(str, ",", &str));
    break;
  } case 'r':
    spec->rfac = atof(arg);
    break;
  case 'c':
    spec->c = atof(arg);
    break;
  case 'f':
    spec->omega = atof(arg);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return EXIT_SUCCESS;
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

void save_bin_files(problem_spec_s const *spec,
                    eik3hh_branch_s const *branch,
                    char const *name) {
  char path[256];

  struct stat st = {0};
  if (stat(name, &st) == -1)
    mkdir(name, 0700);

  sprintf(path, "%s/jet.bin", name);
  eik3hh_branch_dump_jet(branch, path);

  sprintf(path, "%s/org.bin", name);
  eik3hh_branch_dump_org(branch, path);

  sprintf(path, "%s/spread.bin", name);
  eik3hh_branch_dump_spread(branch, path);
}

int main(int argc, char *argv[]) {
  int code = EXIT_SUCCESS;

  problem_spec_s spec = {
    .verbose = false,
    .xsrc = {NAN, NAN, NAN},
    .rfac = 0.1,
    .c = 340.3,
    .omega = 100,
  };

  argp_parse(&argp, argc, argv, 0, 0, &spec);

  if (spec.verbose) {
    printf("Problem info:\n");
    printf("- xsrc: (%g, %g, %g) m\n",spec.xsrc[0],spec.xsrc[1],spec.xsrc[2]);
    printf("- rfac: %g m\n", spec.rfac);
    printf("- omega: %g Hz\n", spec.omega);
  }

  toc();

  dbl eps = 1e-5;

  mesh3_data_s data;
  mesh3_data_init_from_bin(&data, "verts.bin", "cells.bin");
  mesh3_data_insert_vert(&data, spec.xsrc, eps);

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, &eps);

  mesh3_dump_verts(mesh, "verts.bin");
  mesh3_dump_cells(mesh, "cells.bin");

  if (!mesh3_contains_ball(mesh, spec.xsrc, spec.rfac)) {
    fprintf(stderr, "ERROR: mesh doesn't fully contain factoring ball\n");
    code = EXIT_FAILURE;
    goto cleanup;
  }

  /* Set up the Helmholtz data structure to manage the different
   * branches as we solve eikonal problems. Initialize it with a point
   * source whose location is passed on the command line. */
  eik3hh_s *hh;
  eik3hh_alloc(&hh);
  eik3hh_init_with_pt_src(hh, mesh, spec.c, spec.rfac, spec.xsrc);

  /* Solve the "root" eikonal problem (the original point source
   * problem). */
  eik3hh_branch_s *root_branch = eik3hh_get_root_branch(hh);
  eik3hh_branch_solve(root_branch, spec.verbose);

  save_bin_files(&spec, root_branch, "direct");

  /* Find visible reflecting boundaries */
  array_s *refl_inds = eik3hh_branch_get_visible_refls(root_branch);
  printf("Found %lu visible reflections\n", array_size(refl_inds));

  /* Iterate over the reflecting boundaries and solve downwind eikonal
   * problems. Write some stuff to disk for plotting. */
  for (size_t i = 0; i < array_size(refl_inds); ++i) {
    size_t refl_ind;
    array_get(refl_inds, i, &refl_ind);

    eik3hh_branch_s *refl_branch = eik3hh_branch_add_refl(root_branch,refl_ind);
    eik3hh_branch_solve(refl_branch, spec.verbose);

    char name[128];
    sprintf(name, "refl%03d", (int)refl_ind);
    save_bin_files(&spec, refl_branch, name);
  }

  array_deinit(refl_inds);
  array_dealloc(&refl_inds);

cleanup:

  eik3hh_deinit(hh);
  eik3hh_dealloc(&hh);

  /* Free tetrahedron mesh */
  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  return code;
}
