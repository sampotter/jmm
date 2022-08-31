#include <argp.h>
#include <string.h>

#include <jmm/mesh3.h>
#include <jmm/util.h>

const char *argp_program_version = "simulate 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] =
  "Read the boundary mesh contained in OFF_PATH and tetrahedralize "
  "the interior of the domain. The results are written to verts.bin and "
  "cells.bin in row-major order.\n"
  "\n"
  "(At the moment, this is just a thin CLI TetGen wrapper.)";

static char args_doc[] = "OFF_PATH";

static struct argp_option options[] = {
  {"verbose", 'v', 0, OPTION_ARG_OPTIONAL, "Produce verbose output", 0},
  {"maxvol", 'a', "VOLUME", OPTION_ARG_OPTIONAL,
   "Maximum tetrahedron volume (default: 0.01)", 0},
  {"eps", 'e', "EPS", OPTION_ARG_OPTIONAL,
   "Mesh epsilon (default: 1e-5)", 0},
  {0}
};

typedef struct mesh_spec {
  bool verbose;
  dbl maxvol;
  dbl eps;
  char *off_path;
} mesh_spec_s;

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  mesh_spec_s *spec = state->input;
  size_t n;
  switch (key) {
  case 'v':
    spec->verbose = true;
    break;
  case 'a':
    spec->maxvol = atof(arg);
    break;
  case 'e':
    spec->eps = atof(arg);
    break;
  case ARGP_KEY_ARG:
    n = strlen(arg);
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

  mesh_spec_s spec = {
    .verbose = false,
    .maxvol = 0.01,
    .eps = 1e-5,
    .off_path = NULL
  };

  argp_parse(&argp, argc, argv, 0, 0, &spec);

  if (spec.verbose) {
    printf("Mesh info:\n");
    printf("- maxvol: %g m^3\n", spec.maxvol);
    printf("- off_path: %s\n", spec.off_path);
  }

  toc();

  mesh3_data_s data;
  mesh3_data_init_from_off_file(&data, spec.off_path, spec.maxvol, spec.verbose);

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, &spec.eps);

  if (spec.verbose) {
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
    printf("- wrote mesh vertices to verts.bin\n");
    printf("- wrote mesh tetrahedra to cells.bins\n");
  }

  mesh3_dump_verts(mesh, "verts.bin");
  mesh3_dump_cells(mesh, "cells.bin");

  /* Free tetrahedron mesh */
  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  return code;
}
