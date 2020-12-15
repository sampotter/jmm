#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <math.h>

#include "eik3.h"
#include "mesh3.h"

#define NUM_RANDOM_TRIALS 10

Describe(eik3);

BeforeEach(eik3) {
  significant_figures_for_assert_double_are(15);
}

AfterEach(eik3) {}

Ensure(eik3, tetra) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

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

  jet3 const jet_gt = {
    .f = 1.821367205045918,
    .fx = 0.57735026918962551,
    .fy = 0.57735026918962551,
    .fz = 0.57735026918962551
  };

  costfunc_s cf;
  costfunc_init(&cf, mesh, jet, 3, 0, 1, 2);

  costfunc_set_lambda(&cf, lambda);
  tetra(&cf, lambda, &newjet);
  assert_that_double(lambda[0], is_equal_to_double(1./3));
  assert_that_double(lambda[0], is_equal_to_double(1./3));
  assert_that_double(jet_gt.f, is_equal_to_double(newjet.f));
  assert_that_double(jet_gt.fx, is_equal_to_double(newjet.fx));
  assert_that_double(jet_gt.fy, is_equal_to_double(newjet.fy));
  assert_that_double(jet_gt.fz, is_equal_to_double(newjet.fz));
  assert_that(cf.niter, is_equal_to(0));

  lambda[0] = 0;
  lambda[1] = 0;
  costfunc_set_lambda(&cf, lambda);
  tetra(&cf, lambda, &newjet);
  assert_that_double(lambda[0], is_equal_to_double(1./3));
  assert_that_double(lambda[1], is_equal_to_double(1./3));
  assert_that_double(jet_gt.f, is_equal_to_double(newjet.f));
  assert_that_double(jet_gt.fx, is_equal_to_double(newjet.fx));
  assert_that_double(jet_gt.fy, is_equal_to_double(newjet.fy));
  assert_that_double(jet_gt.fz, is_equal_to_double(newjet.fz));
  assert_that(cf.niter, is_less_than(10));

  lambda[0] = 1;
  lambda[1] = 0;
  costfunc_set_lambda(&cf, lambda);
  tetra(&cf, lambda, &newjet);
  assert_that_double(lambda[0], is_equal_to_double(1./3));
  assert_that_double(lambda[1], is_equal_to_double(1./3));
  assert_that_double(jet_gt.f, is_equal_to_double(newjet.f));
  assert_that_double(jet_gt.fx, is_equal_to_double(newjet.fx));
  assert_that_double(jet_gt.fy, is_equal_to_double(newjet.fy));
  assert_that_double(jet_gt.fz, is_equal_to_double(newjet.fz));
  assert_that(cf.niter, is_less_than(10));

  lambda[0] = 0;
  lambda[1] = 1;
  costfunc_set_lambda(&cf, lambda);
  tetra(&cf, lambda, &newjet);
  assert_that_double(lambda[0], is_equal_to_double(1./3));
  assert_that_double(lambda[1], is_equal_to_double(1./3));
  assert_that_double(jet_gt.f, is_equal_to_double(newjet.f));
  assert_that_double(jet_gt.fx, is_equal_to_double(newjet.fx));
  assert_that_double(jet_gt.fy, is_equal_to_double(newjet.fy));
  assert_that_double(jet_gt.fz, is_equal_to_double(newjet.fz));
  assert_that(cf.niter, is_less_than(10));

  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    lambda[0] = gsl_ran_flat(rng, 0, 1);
    lambda[1] = gsl_ran_flat(rng, 0, 1);
    if (lambda[0] + lambda[1] > 1) {
      lambda[0] = 1 - lambda[0];
      lambda[1] = 1 - lambda[1];
    }

    costfunc_set_lambda(&cf, lambda);
    tetra(&cf, lambda, &newjet);

    assert_that_double(lambda[0], is_equal_to_double(1./3));
    assert_that_double(lambda[1], is_equal_to_double(1./3));

    assert_that(cf.niter, is_less_than(10));
  }

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  gsl_rng_free(rng);
}
