#include <argp.h>
#include <string.h>

#include <jmm/mesh3.h>
#include <jmm/eik3hh_branch.h>
#include <jmm/util.h>

const char *argp_program_version = "slice 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = "";

static char args_doc[] = "BIN_PATH";

static struct argp_option options[] = {
  {"verbose", 'v', 0, OPTION_ARG_OPTIONAL, "Produce verbose output", 0},
  {NULL, 'z', 0, OPTION_ARG_OPTIONAL, "z-coord of slice (default value: 0)", 0},
  {"out", 'o', 0, OPTION_ARG_OPTIONAL,
   "Output directory (default value: current directory)", 0},
  {0}
};

typedef struct slice_spec {
  bool verbose;
  dbl z;
  char *out_path;
} slice_spec_s;

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  slice_spec_s *spec = state->input;
  size_t n;
  switch (key) {
  case 'v':
    spec->verbose = true;
    break;
  case 'z':
    spec->z = atof(arg);
    break;
  case 'o':
    n = strlen(arg);
    spec->out_path = malloc(n + 1);
    strncpy(spec->out_path, arg, n + 1);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return EXIT_SUCCESS;
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

int main(int argc, char *argv[]) {
  int code = EXIT_SUCCESS;

  slice_spec_s spec = {
    .verbose = false,
    .z = 0
  };

  argp_parse(&argp, argc, argv, 0, 0, &spec);

  if (spec.verbose) {
    printf("Problem info:\n");
    printf("- z: %g\n", spec.z);
  }

  mesh3_data_s data;
  mesh3_data_init_from_bin(&data, "verts.bin", "cells.bin");

  dbl eps = 1e-5;

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, &eps);

  rect3 bbox;
  mesh3_get_bbox(mesh, &bbox);

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
  grid2_to_mesh3_mapping_init_xy(&mapping, &img_grid, mesh, spec.z);
  if (spec.verbose) {
    printf("Set up grid-to-mesh mapping [%1.2gs]\n", toc());
  }

  // char path[256];

  printf("WARNING: not actually dumping a slice!\n");

  grid2_to_mesh3_mapping_deinit(&mapping);

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  return code;
}
