#include <argp.h>
#include <stdio.h>
#include <string.h>

#include <jmm/bmesh.h>
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

// dbl const dtheta_dy = 1;

// dbl s(dbl3 x) {
//   return 1/(331.3 + 0.606*dtheta_dy*x[1]);
// }

// void Ds(dbl3 x, dbl3 Ds) {
//   // grad_c = 0.606 *grad_theta
//   // grad_theta = [0, dtheta_dy, 0]
//   // dtheta_dy = 1
//   // grad_s = -grad_c/c^2
//   dbl tmp = s(x);
//   Ds[0] = Ds[2] = 0;
//   Ds[1] = -0.606*dtheta_dy*tmp*tmp;
// }

// void D2s(dbl3 x, dbl33 D2s) {
//   dbl33_zero(D2s);
//   dbl tmp = -0.606*dtheta_dy*s(x);
//   D2s[1][1] = -2*tmp*tmp;
// }

static dbl const v0 = 1;
static dbl3 const v = {0.025, -0.025, 0.05};

static dbl const s0 = 1/v0;
static dbl const vnormsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

dbl c(dbl3 x) {
  return v0 + dbl3_dot(v, x);
}

dbl s(dbl3 x) {
  return 1/c(x);
}

void Ds(dbl3 x, dbl3 Ds) {
  dbl s_squared = s(x)*s(x);
  for (size_t i = 0; i < 3; ++i)
    Ds[i] = -s_squared*v[i];
}

void D2s(dbl3 x, dbl33 D2s) {
  dbl s_cubed = pow(s(x), 3);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      D2s[i][j] = -2*s_cubed*v[i]*v[j];
}

/* Analytic solution for a linear speed of sound */
dbl tau(dbl3 x) {
  return acosh(1 + s0*s(x)*vnormsq*dbl3_normsq(x)/2)/sqrt(vnormsq);
}

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

  printf("average edge length = %g\n", mesh3_get_mean_edge_length(mesh));

  if (!mesh3_contains_ball(mesh, spec.xsrc, spec.rfac)) {
    fprintf(stderr, "ERROR: mesh doesn't fully contain factoring ball\n");
    code = EXIT_FAILURE;
    goto cleanup;
  }

  /* Set up the slowness functions */
  sfunc_s sfunc = {
    .stype = STYPE_FUNC_PTR,
    .funcs = {.s = s, .Ds = Ds, .D2s = D2s}
  };

  /* Solve the eikonal equation */
  eik3_s *eik;
  eik3_alloc(&eik);
  eik3_init(eik, mesh, &sfunc);
  eik3_add_pt_src_bcs(eik, spec.xsrc, spec.rfac);
  eik3_solve(eik);
  eik3_dump_jet(eik, "jet_T.bin");

  printf("wrote jet_T.bin\n");

  /** Just dump the solution now for error checking... */

  jet31t const *jet = eik3_get_jet_ptr(eik);

  bmesh33_s *bmesh;
  bmesh33_alloc(&bmesh);
  bmesh33_init_from_mesh3_and_jets(bmesh, mesh, jet);

  size_t N = 256;

  dbl *value = malloc(N*N*sizeof(dbl));

  size_t k = 0;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      dbl3 x = {(2.0/N)*i - 1.0, (2.0/N)*j - 1.0, 0.5};
      value[k++] = bmesh33_f(bmesh, x);
    }
  }

  FILE *fp = fopen("T_grid.bin", "w");
  fwrite(value, sizeof(dbl), N*N, fp);
  fclose(fp);

  bmesh33_deinit(bmesh);
  bmesh33_dealloc(&bmesh);

  printf("wrote T_grid.bin\n");

  k = 0;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      dbl3 x = {(2.0/N)*i - 1.0, (2.0/N)*j - 1.0, 0.5};
      value[k++] = tau(x);
    }
  }

  fp = fopen("tau_grid.bin", "w");
  fwrite(value, sizeof(dbl), N*N, fp);
  fclose(fp);

  printf("wrote tau_grid.bin\n");

cleanup:

  /* Free tetrahedron mesh */
  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  return code;
}
