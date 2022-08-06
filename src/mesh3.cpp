#include "mesh3.h"

#include <assert.h>
#include <iostream>
#include <sstream>

#define TETLIBRARY 1
#include "../extra/tetgen1.6.0/tetgen.h"

void mesh3_data_from_off_file(mesh3_data_s *data, char const *path, dbl maxvol, bool verbose) {
  /* Set up string of command-line switches for TetGen */
  std::ostringstream oss;
  oss << "a" << maxvol
      << "p"
      << "q1.414"
    ;
  if (!verbose)
    oss << "Q";
  std::string switch_str = oss.str();

  /* Tetrahedralize the input OFF file */
  tetgenio in, out;
  in.load_plc((char *)path, (int)tetgenbehavior::OFF);
  tetrahedralize((char *)switch_str.c_str(), &in, &out);

  /* TODO: initializing a little inefficiently here by calling
   * mesh3_init because I copied this code from the 3d_wedge
   * example... should rewrite this to just set up the arrays directly
   * and *not* call mesh3_init */

  /* Copy over the vertices */
  size_t nverts = out.numberofpoints;
  double *verts = (double *)malloc(3*nverts*sizeof(double));
  for (size_t i = 0; i < 3*nverts; ++i)
    verts[i] = out.pointlist[i];

  /* Copy over the cells */
  size_t ncells = out.numberoftetrahedra;
  size_t *cells = (size_t *)malloc(4*ncells*sizeof(size_t));
  for (size_t i = 0; i < 4*ncells; ++i)
    cells[i] = out.tetrahedronlist[i];

  /* Initialize the mesh */
  mesh3_init(mesh, verts, nverts, cells, ncells, true, NULL);
}
