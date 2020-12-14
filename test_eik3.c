#include <cgreen/cgreen.h>

#include <math.h>

#include "eik3.h"
#include "mesh3.h"

Describe(costfunc);
BeforeEach(costfunc) {}
AfterEach(costfunc) {}

Ensure(costfunc, tetra) {
  dbl verts[12] = {
    1, 0, 0, // x0
    0, 1, 0, // x1
    0, 0, 1, // x2
    1, 1, 1  // xhat
  };

  size_t cells[4] = {0, 1, 2, 3};

  jet3 jet[3] = {
    {.f = 1, .fx = 1, .fy = 0, .fz = 0},
    {.f = 1, .fx = 0, .fy = 1, .fz = 0},
    {.f = 1, .fx = 0, .fy = 0, .fz = 1}
  };

  mesh3_s *mesh;
  mesh3_alloc(&mesh);
  mesh3_init(mesh, verts, 4, cells, 1);

  dbl lambda[2] = {1./3, 1./3};

  jet3 newjet;

  costfunc_s cf;
  costfunc_init(&cf, mesh, jet, 3, 0, 1, 2);
  costfunc_set_lambda(&cf, lambda);
  tetra(&cf, lambda, &newjet);

  assert_that(fabs(lambda[0] - 1./3) < 1e-15);
  assert_that(fabs(lambda[0] - 1./3) < 1e-15);
  assert_that(cf.niter, is_equal_to(0));

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);
}
