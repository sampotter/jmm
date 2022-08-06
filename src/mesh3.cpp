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
  data->nverts = out.numberofpoints;
  data->verts = (dbl3 *)malloc(data->nverts*sizeof(dbl3));
  memcpy(data->verts, out.pointlist, data->nverts*sizeof(dbl3));

  /* Copy over the cells */
  data->ncells = out.numberoftetrahedra;
  data->cells = (uint4 *)malloc(data->ncells*sizeof(uint4));
  for (size_t lc = 0; lc < data->ncells; ++lc)
    for (size_t i = 0; i < 4; ++i)
      data->cells[lc][i] = out.tetrahedronlist[4*lc + i];
}
