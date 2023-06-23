#include <argp.h>
#include <assert.h>

#include <sstream>

#define TETLIBRARY 1
#include <tetgen.h>

#include <jmm/eik3.h>
#include <jmm/mesh3.h>
#include <jmm/util.h>

const char *argp_program_version = "varying_s_convergence_test 0.0";
char const *argp_program_bug_address = "sfp@cims.nyu.edu";

static char doc[] = ""; // TODO: add some docs here when it makes sense

static char args_doc[] = "EXAMPLE";

static struct argp_option options[] = {
  {"verbose", 'v', 0, OPTION_ARG_OPTIONAL, "Produce verbose output", 0},
  {"rfac", 'r', "RADIUS", OPTION_ARG_OPTIONAL,
   "Factoring radius (default: 0.2)", 0},
  {"maxvol", 'a', "RADIUS", OPTION_ARG_OPTIONAL,
   "Factoring radius (default: 0.1)", 0},
  {0}
};

typedef struct problem_spec {
  bool verbose;
  dbl rfac;
  dbl maxvol;
} problem_spec_s;

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  problem_spec_s *spec = (problem_spec_s *)state->input;
  switch (key) {
  case 'v':
    spec->verbose = true;
    break;
  case 'r':
    spec->rfac = atof(arg);
    break;
  case 'a':
    spec->maxvol = atof(arg);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return EXIT_SUCCESS;
}

static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

#define X(i) in.pointlist[3*(i)]
#define Y(i) in.pointlist[3*(i) + 1]
#define Z(i) in.pointlist[3*(i) + 2]

static void set_tetgenio_for_cube(tetgenio &in) {
  in.numberofpoints = 8;
  in.pointlist = new REAL[3*in.numberofpoints];
  in.numberoffacets = 6;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

  X(0) = -1; Y(0) = -1; Z(0) = -1;
  X(1) = -1; Y(1) =  1; Z(1) = -1;
  X(2) =  1; Y(2) = -1; Z(2) = -1;
  X(3) =  1; Y(3) =  1; Z(3) = -1;
  X(4) = -1; Y(4) = -1; Z(4) =  1;
  X(5) = -1; Y(5) =  1; Z(5) =  1;
  X(6) =  1; Y(6) = -1; Z(6) =  1;
  X(7) =  1; Y(7) =  1; Z(7) =  1;

  auto set_facet = [&] (int i, int j0, int j1, int j2, int j3) {
    tetgenio::facet *f;
    tetgenio::polygon *p;

    f = &in.facetlist[i];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;

    p = &f->polygonlist[i];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = j0;
    p->vertexlist[1] = j1;
    p->vertexlist[2] = j2;
    p->vertexlist[3] = j3;
  };

  set_facet(0, 0, 1, 3, 2); // bottom face
  set_facet(1, 4, 6, 7, 5); // top face
  set_facet(2, 0, 2, 6, 4); // front face
  set_facet(3, 2, 3, 7, 6); // right face
  set_facet(4, 3, 1, 5, 7); // back face
  set_facet(5, 0, 4, 5, 1); // left face

  for (size_t i = 0; i < (size_t)in.numberoffacets; ++i)
    in.facetmarkerlist[i] = i;
}

mesh3_s *tetrahedralize_cube(dbl maxvol, dbl3 const xsrc) {
  tetgenio in, out, addin;

  /* set up flags to pass to tetgen */
  std::ostringstream oss;
  oss << "p" /* read a PLC */
      << "q1.414" /* quality bound == 1.414 */
      << "a" << maxvol /* max volume constraint */
      << "i";
  std::string switches = oss.str();

  set_tetgenio_for_cube(in);

  addin.numberofpoints = 1;
  addin.pointlist = new REAL[3];
  addin.pointlist[0] = addin.pointlist[1] = addin.pointlist[2] = 0;

  tetrahedralize((char *)switches.c_str(), &in, &out, &addin);

  mesh3_data_s data;
  data.nverts = out.numberofpoints;
  data.verts = (dbl3 *)malloc(data.nverts*sizeof(dbl3));
  data.ncells = out.numberoftetrahedra;
  data.cells = (uint4 *)malloc(data.ncells*sizeof(uint4));

  mesh3_data_insert_vert(&data, xsrc, 1e-10);

  memcpy(data.verts, out.pointlist, data.nverts*sizeof(dbl3));
  memcpy(data.cells, out.tetrahedronlist, data.ncells*sizeof(uint4));

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, &data, true, NULL);

  /* Make sure the point source is actually included in the mesh! */
  assert(mesh3_has_vertex(mesh, addin.pointlist));

  return mesh;
}

static dbl const v0 = 1;
static dbl3 const v = {0.025, -0.025, 0.05};

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

int main(int argc, char *argv[]) {
  int code = EXIT_SUCCESS;

  problem_spec_s spec = {
    .verbose = false,
    .rfac = 0.1,
    .maxvol = 0.2,
  };

  argp_parse(&argp, argc, argv, 0, 0, &spec);

  if (spec.verbose) {
    printf("Problem info:\n");
    printf("- rfac: %g m\n", spec.rfac);
    printf("- maxvol: %g m\n", spec.maxvol);
  }

  toc();

  dbl3 xsrc = {0, 0, 0};

  /* Set up the slowness functions */
  sfunc_s sfunc = {
    .stype = STYPE_FUNC_PTR,
    .funcs = {.s = s, .Ds = Ds, .D2s = D2s}
  };

  mesh3_s *mesh = tetrahedralize_cube(spec.maxvol, xsrc);

  if (!mesh3_contains_ball(mesh, xsrc, spec.rfac)) {
    fprintf(stderr, "ERROR: mesh doesn't fully contain factoring ball\n");
    code = EXIT_FAILURE;
    goto cleanup;
  }

  /* Solve the eikonal equation */
  eik3_s *eik;
  eik3_alloc(&eik);
  eik3_init(eik, mesh, &sfunc);
  eik3_add_pt_src_bcs(eik, xsrc, spec.rfac);
  eik3_solve(eik);
  eik3_dump_jet(eik, "jet_T.bin");

cleanup:

  /* Free tetrahedron mesh */
  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  return code;
}
