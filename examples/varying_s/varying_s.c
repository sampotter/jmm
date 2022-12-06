#include <argp.h>
#include <stdio.h>
#include <string.h>

#include <jmm/eik3.h>
#include <jmm/mesh3.h>
#include <jmm/util.h>

const char *argp_program_version = "simulate 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "EXAMPLE";

static struct argp_option options[] = {
  {"verbose", 'v', 0, OPTION_ARG_OPTIONAL, "Produce verbose output", 0},
  {"xsrc", 'x', "POINT", OPTION_ARG_OPTIONAL,
   "Point source location (format: \"x,y,z\")", 0},
  {"rfac", 'r', "RADIUS", OPTION_ARG_OPTIONAL,
   "Factoring radius (default: 0.1)", 0},
  {0}
};

typedef struct problem_spec {
  bool verbose;
  dbl3 xsrc;
  dbl rfac;
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
    .xsrc = {NAN, NAN, NAN},
    .rfac = 0.1,
  };

  dbl eps = 1e-5;

  argp_parse(&argp, argc, argv, 0, 0, &spec);

  if (spec.verbose) {
    printf("Problem info:\n");
    printf("- xsrc: (%g, %g, %g) m\n",spec.xsrc[0],spec.xsrc[1],spec.xsrc[2]);
    printf("- rfac: %g m\n", spec.rfac);
  }

  toc();

  mesh3_data_s data;
  mesh3_data_init_from_bin(&data, "verts.bin", "cells.bin");
  mesh3_data_insert_vert(&data, spec.xsrc, eps);

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, &eps);

  if (!mesh3_contains_ball(mesh, spec.xsrc, spec.rfac)) {
    fprintf(stderr, "ERROR: mesh doesn't fully contain factoring ball\n");
    code = EXIT_FAILURE;
    goto cleanup;
  }

  /* Read the jet data for the slowness functions */
  jet31t *s_data = malloc(mesh3_nverts(mesh)*sizeof(jet31t));
  FILE *fp = fopen("s_data.bin", "r");
  fread(s_data, sizeof(jet31t), mesh3_nverts(mesh), fp);
  fclose(fp);

  /* Set up the slowness functions */
  sfunc_s s = {
    .stype = STYPE_JET31T,
    .data_jet31t = s_data
  };

  /* Solve the eikonal equation for the left ear */
  eik3_s *eik;
  eik3_alloc(&eik);
  eik3_init(eik, mesh, &s);
  eik3_add_pt_src_bcs(eik, spec.xsrc, spec.rfac);
  eik3_solve(eik);
  eik3_dump_jet(eik, "jet_T.bin");

cleanup:

  /* Free tetrahedron mesh */
  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  return code;
}
