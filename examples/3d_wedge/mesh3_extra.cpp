#include "mesh3_extra.h"

#include "3d_wedge.h"

#include <sstream>

#define TETLIBRARY 1
#include "tetgen.h"

/* Functions for building the PLC describing the boundary of the
 * domain. There are five cases to handle, shown in the following
 * diagram:
 *
 *    n=3/4     n=1/4
 *      +---------+
 *      |\       /|
 *      | \  2  / |
 *      |  \   /  |
 *      |   \ / 1 |
 *      | 3  X----+ n=0
 *      |   / \ 5 |
 *      |  /   \  |
 *      | /  4  \ |
 *      |/       \|
 *      +---------+
 *    n=5/4     n=7/4
 *
 * The five cases correspond to the possible locations of the n-face
 * which result in different PLCs. The horizontal line between
 * regions 1 and 5 is the o-face. */

#define X(i) in.pointlist[3*(i)]
#define Y(i) in.pointlist[3*(i) + 1]
#define Z(i) in.pointlist[3*(i) + 2]

static void init_tetgenio_before(tetgenio & in,
                                 jmm_3d_wedge_spec_s const *spec,
                                 size_t num_poly_corners,
                                 dbl & xmin, dbl & ymin, dbl & zmin,
                                 dbl & xmax, dbl & ymax, dbl & zmax,
                                 dbl & theta)
{
  in.numberofpoints = 2*num_poly_corners;
  in.pointlist = new REAL[3*in.numberofpoints];

  xmin = -spec->w/2, xmax = spec->w/2;
  ymin = -spec->w/2, ymax = spec->w/2;
  zmin = -spec->h/2, zmax = spec->h/2;
  theta = JMM_PI*spec->n;

  for (size_t i = 0; i < num_poly_corners; ++i) {
    Z(i) = zmin;
    Z(i + num_poly_corners) = zmax;
  }
}

static void init_tetgenio_after(tetgenio & in,
                                size_t num_poly_corners,
                                dbl zmax) {
  /* The next num_poly_corners are the top polygon. Same as before but with
   * Z(i) = zmax instead of zmin. */

  for (size_t i = 0; i < num_poly_corners; ++i) {
    X(i + num_poly_corners) = X(i);
    Y(i + num_poly_corners) = Y(i);
    Z(i + num_poly_corners) = zmax;
  }

  /* The PLC has 7 (vertical) + 2 (horizontal) = 9 facets. Take care
   * of memory allocations for each facet now. */

  in.numberoffacets = num_poly_corners + 2;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

  tetgenio::facet *f;
  tetgenio::polygon *p;

  /* bottom facet */

  f = &in.facetlist[0];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;

  p = &f->polygonlist[0];
  p->numberofvertices = num_poly_corners;
  p->vertexlist = new int[p->numberofvertices];
  for (size_t i = 0; i < num_poly_corners; ++i)
    p->vertexlist[i] = i;

  /* top facet */

  f = &in.facetlist[1];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;

  p = &f->polygonlist[0];
  p->numberofvertices = num_poly_corners;
  p->vertexlist = new int[p->numberofvertices];
  for (size_t i = 0; i < num_poly_corners; ++i)
    p->vertexlist[i] = i + num_poly_corners;

  /* add the remaining facets. each facet is a vertical rectangle
   * connecting the ith and (i + 1)st corners on the bottom and top
   * facets */

  for (size_t i = 0; i < num_poly_corners; ++i) {
    f = &in.facetlist[i + 2];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;

    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = i;
    p->vertexlist[1] = (i + 1) % num_poly_corners;
    p->vertexlist[2] = p->vertexlist[1] + num_poly_corners;
    p->vertexlist[3] = p->vertexlist[0] + num_poly_corners;
  }

  /* set facet marker list... TODO: not 100% sure what this does */

  for (size_t i = 0; i < (size_t)in.numberoffacets; ++i)
    in.facetmarkerlist[i] = i;
}

static void init_tetgenio_in_case1(tetgenio & in, jmm_3d_wedge_spec_s const *spec) {
  size_t num_poly_corners = 7;
  dbl xmin, ymin, zmin, xmax, ymax, zmax, theta;

  init_tetgenio_before(
    in, spec, num_poly_corners, xmin, ymin, zmin, xmax, ymax, zmax, theta);

  X(0) = 0;     Y(0) = 0;
  X(1) = xmax;  Y(1) = ymax*sin(theta);
  X(2) = xmax;  Y(2) = ymax;
  X(3) = xmin;  Y(3) = ymax;
  X(4) = xmin;  Y(4) = ymin;
  X(5) = xmax;  Y(5) = ymin;
  X(6) = xmax;  Y(6) = 0;

  init_tetgenio_after(in, num_poly_corners, zmax);
}

static void init_tetgenio_in_case2(tetgenio & in, jmm_3d_wedge_spec_s const *spec) {
  size_t num_poly_corners = 6;
  dbl xmin, ymin, zmin, xmax, ymax, zmax, theta;

  init_tetgenio_before(
    in, spec, num_poly_corners, xmin, ymin, zmin, xmax, ymax, zmax, theta);

  X(0) = 0;                         Y(0) = 0;
  X(1) = xmax*sin(JMM_PI - theta);  Y(1) = ymax;
  X(2) = xmin;                      Y(2) = ymax;
  X(3) = xmin;                      Y(3) = ymin;
  X(4) = xmax;                      Y(4) = ymin;
  X(5) = xmax;                      Y(5) = 0;

  init_tetgenio_after(in, num_poly_corners, zmax);
}

static void init_tetgenio_in_case3(tetgenio & in, jmm_3d_wedge_spec_s const *spec) {
  size_t num_poly_corners = 5;
  dbl xmin, ymin, zmin, xmax, ymax, zmax, theta;

  init_tetgenio_before(
    in, spec, num_poly_corners, xmin, ymin, zmin, xmax, ymax, zmax, theta);

  X(0) = 0;     Y(0) = 0;
  X(1) = xmin;  Y(1) = ymax*sin(theta);
  X(2) = xmin;  Y(2) = ymin;
  X(3) = xmin;  Y(3) = ymax;
  X(4) = xmax;  Y(4) = 0;

  init_tetgenio_after(in, num_poly_corners, zmax);
}

static void init_tetgenio_in_case4(tetgenio & in, jmm_3d_wedge_spec_s const *spec) {
  size_t num_poly_corners = 4;
  dbl xmin, ymin, zmin, xmax, ymax, zmax, theta;

  init_tetgenio_before(
    in, spec, num_poly_corners, xmin, ymin, zmin, xmax, ymax, zmax, theta);

  X(0) = 0;                                 Y(0) = 0;
  X(1) = xmin*sin(theta - (5.0/4)*JMM_PI);  Y(1) = ymin;
  X(2) = xmax;                              Y(2) = ymin;
  X(3) = xmax;                              Y(3) = 0;

  init_tetgenio_after(in, num_poly_corners, zmax);
}

static void init_tetgenio_in_case5(tetgenio & in, jmm_3d_wedge_spec_s const *spec) {
  size_t num_poly_corners = 3;
  dbl xmin, ymin, zmin, xmax, ymax, zmax, theta;

  init_tetgenio_before(
    in, spec, num_poly_corners, xmin, ymin, zmin, xmax, ymax, zmax, theta);

  X(0) = 0;     Y(0) = 0;
  X(1) = xmax;  Y(1) = ymin*sin(-theta);
  X(2) = xmax;  Y(2) = 0;

  init_tetgenio_after(in, num_poly_corners, zmax);
}

#undef X
#undef Y
#undef Z

static std::string get_switches_str(jmm_3d_wedge_spec_s const *spec) {
  std::ostringstream oss;

  oss << "p" /* read a PLC */
      << "q1.414" /* "quality" mesh generation w/ quality bound = 1.414 */
      << "a" << spec->maxvol /* maximum volume constraint */
    ;

  return oss.str();
}

jmm_error_e mesh3_init_from_3d_wedge_spec(mesh3_s *mesh, jmm_3d_wedge_spec_s const *spec)
{
  tetgenio in, out;

  if (0 < spec->n && spec->n < 0.25) init_tetgenio_in_case1(in, spec);
  else if (0.25 <= spec->n && spec->n < 0.75) init_tetgenio_in_case2(in, spec);
  else if (0.75 <= spec->n && spec->n < 1.25) init_tetgenio_in_case3(in, spec);
  else if (1.25 <= spec->n && spec->n < 1.75) init_tetgenio_in_case4(in, spec);
  else if (1.75 <= spec->n && spec->n < 2) init_tetgenio_in_case5(in, spec);
  else return JMM_ERROR_BAD_ARGUMENTS;

  std::string switches = get_switches_str(spec);

  tetrahedralize((char *)switches.c_str(), &in, &out);

  size_t nverts = out.numberofpoints;
  double *verts = (double *)malloc(3*nverts*sizeof(double));
  for (size_t i = 0; i < 3*nverts; ++i)
    verts[i] = out.pointlist[i];

  size_t ncells = out.numberoftetrahedra;
  size_t *cells = (size_t *)malloc(4*ncells*sizeof(size_t));
  for (size_t i = 0; i < 4*ncells; ++i)
    cells[i] = out.tetrahedronlist[i];

  mesh3_init(mesh, verts, nverts, cells, ncells, true, NULL);

  rect3 bbox;
  mesh3_get_bbox(mesh, &bbox);

  printf("bounding box: [%g, %g] x [%g, %g] x [%g, %g]\n",
         bbox.min[0], bbox.max[0], bbox.min[1], bbox.max[1],
         bbox.min[2], bbox.max[2]);

  return JMM_ERROR_NONE;
}
