#include <cgreen/cgreen.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <math.h>
#include <string.h>

#include "eik3.h"
#include "mesh3.h"
#include "utetra.h"
#include "vec.h"

#define NUM_RANDOM_TRIALS 10

void get_gt_jet(dbl const xsrc[3], dbl const x[3], jet3 *jet) {
  dbl tmp[3];
  dbl3_sub(x, xsrc, tmp);
  dbl L = dbl3_norm(tmp);
  jet->f = L;
  dbl3_dbl_div(tmp, L, &jet->fx);
}

Describe(eik3);

BeforeEach(eik3) {
  double_absolute_tolerance_is(1e-15);
  double_relative_tolerance_is(1e-15);
}

AfterEach(eik3) {}

Ensure(eik3, tetra_works_for_olim18_122_update) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  dbl xsrc[3] = {0, 0, 0};

  dbl verts[12] = {
    1, 0, 0, // x0
    0, 1, 0, // x1
    1, 1, 0, // x2
    1, 1, 1  // xhat
  };

  jet3 newjet, jet[3];

  int perm[6][3] = {
    {0, 1, 2},
    {0, 2, 1},
    {1, 0, 2},
    {2, 0, 1},
    {1, 2, 0},
    {2, 1, 0}
  };

  dbl verts_perm[12];
  size_t cells_perm[4];

  // Last vertex stays fixed under permutation
  cells_perm[3] = 3;
  memcpy((void *)&verts_perm[9], (void *)&verts[9], sizeof(dbl)*3);

  dbl lam[2];

  dbl lam_gt[6][2] = {
    {0.5, 0.0},
    {0.0, 0.5},
    {0.5, 0.0},
    {0.5, 0.5},
    {0.0, 0.5},
    {0.5, 0.5}
  };

  dbl lam_verts[3][2] = {{0, 0}, {1, 0}, {0, 1}};

  mesh3_s *mesh;
  mesh3_alloc(&mesh);

  utetra_s *cf;
  utetra_alloc(&cf);

  utetra_spec_s spec;

  int p;
  for (int i = 0; i < 6; ++i) {
    /**
     * Get vertices for this permutation
     */
    for (int j = 0; j < 3; ++j) {
      cells_perm[j] = p = perm[i][j];
      memcpy((void *)&verts_perm[3*j], (void *)&verts[3*p], sizeof(dbl)*3);
    }

    /**
     * Create a mesh consisting of a single tetrahedron for this update
     */
    mesh3_init(mesh, verts_perm, 4, cells_perm, 1, true, NULL);

    /**
     * Get jets for vertex data
     */
    for (int j = 0; j < 3; ++j) {
      get_gt_jet(xsrc, &verts_perm[3*j],  &jet[j]);
    }

    spec = utetra_spec_from_ptrs(mesh, jet, 3, 0, 1, 2);
    spec.tol = 1e-15;
    utetra_init(cf, &spec);

    /**
     * Verify that cost function has correct nodal values
     */
    for (int j = 0; j < 3; ++j) {
      utetra_set_lambda(cf, lam_verts[perm[i][j]]);
      dbl L = dbl3_dist(&verts_perm[3*j], &verts[9]);
      dbl f = jet[j].f + L;
      assert_that_double(utetra_get_value(cf), is_nearly_double(f));
    }

    /**
     * Do tetrahedron updates starting from random initial iterates
     */
    for (int l = 0; l < NUM_RANDOM_TRIALS; ++l) {
      /**
       * Sample random a lambda
       */
      lam[0] = gsl_ran_flat(rng, 0, 1);
      lam[1] = gsl_ran_flat(rng, 0, 1);
      if (lam[0] + lam[1] > 1) {
        lam[0] = 1 - lam[0];
        lam[1] = 1 - lam[1];
      }

      utetra_solve(cf, lam);
      utetra_get_lambda(cf, lam);
      utetra_get_jet(cf, &newjet);

      if (fabs(lam[0]) < 1e-15) {
        assert_that(fabs(lam_gt[i][0]) < 1e-15);
      } else {
        assert_that_double(lam[0], is_nearly_double(lam_gt[i][0]));
      }

      if (fabs(lam[1]) < 1e-15) {
        assert_that(fabs(lam_gt[i][1]) < 1e-15);
      } else {
        assert_that_double(lam[1], is_nearly_double(lam_gt[i][1]));
      }
    }

    mesh3_deinit(mesh);
  }

  utetra_deinit(cf);
  utetra_dealloc(&cf);

  mesh3_dealloc(&mesh);

  gsl_rng_free(rng);
}

Ensure(eik3, tetra_works_for_olim18_222_update) {
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
  mesh3_init(mesh, verts, 4, cells, 1, true, NULL);

  dbl lambda[2] = {1./3, 1./3};

  jet3 newjet;

  jet3 const jet_gt = {
    .f = 1.821367205045918,
    .fx = 0.57735026918962551,
    .fy = 0.57735026918962551,
    .fz = 0.57735026918962551
  };

  utetra_s *cf;
  utetra_alloc(&cf);

  utetra_spec_s spec = utetra_spec_from_ptrs(mesh, jet, 3, 0, 1, 2);
  spec.tol = 1e-15;
  utetra_init(cf, &spec);

  utetra_solve(cf, lambda);
  utetra_get_jet(cf, &newjet);
  assert_that_double(lambda[0], is_nearly_double(1./3));
  assert_that_double(lambda[0], is_nearly_double(1./3));
  assert_that_double(jet_gt.f, is_nearly_double(newjet.f));
  assert_that_double(jet_gt.fx, is_nearly_double(newjet.fx));
  assert_that_double(jet_gt.fy, is_nearly_double(newjet.fy));
  assert_that_double(jet_gt.fz, is_nearly_double(newjet.fz));
  assert_that(utetra_get_num_iter(cf), is_equal_to(0));

  lambda[0] = 0;
  lambda[1] = 0;
  utetra_solve(cf, lambda);
  utetra_get_lambda(cf, lambda);
  utetra_get_jet(cf, &newjet);
  assert_that_double(lambda[0], is_nearly_double(1./3));
  assert_that_double(lambda[1], is_nearly_double(1./3));
  assert_that_double(jet_gt.f, is_nearly_double(newjet.f));
  assert_that_double(jet_gt.fx, is_nearly_double(newjet.fx));
  assert_that_double(jet_gt.fy, is_nearly_double(newjet.fy));
  assert_that_double(jet_gt.fz, is_nearly_double(newjet.fz));
  assert_that(utetra_get_num_iter(cf), is_less_than(10));

  lambda[0] = 1;
  lambda[1] = 0;
  utetra_solve(cf, lambda);
  utetra_get_lambda(cf, lambda);
  utetra_get_jet(cf, &newjet);
  assert_that_double(lambda[0], is_nearly_double(1./3));
  assert_that_double(lambda[1], is_nearly_double(1./3));
  assert_that_double(jet_gt.f, is_nearly_double(newjet.f));
  assert_that_double(jet_gt.fx, is_nearly_double(newjet.fx));
  assert_that_double(jet_gt.fy, is_nearly_double(newjet.fy));
  assert_that_double(jet_gt.fz, is_nearly_double(newjet.fz));
  assert_that(utetra_get_num_iter(cf), is_less_than(10));

  lambda[0] = 0;
  lambda[1] = 1;
  utetra_solve(cf, lambda);
  utetra_get_lambda(cf, lambda);
  utetra_get_jet(cf, &newjet);
  assert_that_double(lambda[0], is_nearly_double(1./3));
  assert_that_double(lambda[1], is_nearly_double(1./3));
  assert_that_double(jet_gt.f, is_nearly_double(newjet.f));
  assert_that_double(jet_gt.fx, is_nearly_double(newjet.fx));
  assert_that_double(jet_gt.fy, is_nearly_double(newjet.fy));
  assert_that_double(jet_gt.fz, is_nearly_double(newjet.fz));
  assert_that(utetra_get_num_iter(cf), is_less_than(10));

  for (int i = 0; i < NUM_RANDOM_TRIALS; ++i) {
    lambda[0] = gsl_ran_flat(rng, 0, 1);
    lambda[1] = gsl_ran_flat(rng, 0, 1);
    if (lambda[0] + lambda[1] > 1) {
      lambda[0] = 1 - lambda[0];
      lambda[1] = 1 - lambda[1];
    }

    utetra_solve(cf, lambda);
    utetra_get_lambda(cf, lambda);
    utetra_get_jet(cf, &newjet);

    assert_that_double(lambda[0], is_nearly_double(1./3));
    assert_that_double(lambda[1], is_nearly_double(1./3));

    assert_that(utetra_get_num_iter(cf), is_less_than(10));
  }

  utetra_deinit(cf);
  utetra_dealloc(&cf);

  mesh3_deinit(mesh);
  mesh3_dealloc(&mesh);

  gsl_rng_free(rng);
}

Ensure(eik3, olim18_222_is_symmetric) {
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
  mesh3_init(mesh, verts, 4, cells, 1, true, NULL);

  utetra_spec_s spec = utetra_spec_from_ptrs(mesh, jet, 3, 0, 1, 2);

  utetra_s *cf[2];
  for (int i = 0; i < 2; ++i) {
    utetra_alloc(&cf[i]);
    utetra_init(cf[i], &spec);
  }

  utetra_set_lambda(cf[0], (dbl[2]) {1, 0});
  utetra_set_lambda(cf[1], (dbl[2]) {0, 1});

  dbl lam[2][2];

  for (int k = 0; k < 10; ++k) {
    for (int i = 0; i < 2; ++i)
      utetra_get_lambda(cf[i], lam[i]);

    assert_that_double(lam[0][0], is_nearly_double(lam[1][1]));
    assert_that_double(lam[1][0], is_nearly_double(lam[0][1]));

    for (int i = 0; i < 2; ++i)
      utetra_step(cf[i]);
  }
}

Ensure(eik3, tetra_works_for_olim26_updates) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

  dbl xsrc[3] = {0, 0, 0};

  dbl verts[12] = {
    1, 1, 1, // x0
    2, 1, 1, // x1
    2, 2, 1, // x2
    2, 2, 2  // xhat
  };

  dbl3 x, x_gt = {1, 1, 1};

  jet3 newjet, jet[3];

  int perm[6][3] = {
    {0, 1, 2},
    {0, 2, 1},
    {1, 0, 2},
    {2, 0, 1},
    {1, 2, 0},
    {2, 1, 0}
  };

  dbl verts_perm[12];

  // Last vertex stays fixed under permutation
  memcpy((void *)&verts_perm[9], (void *)&verts[9], sizeof(dbl)*3);

  size_t cells[4] = {0, 1, 2, 3};

  dbl lam[2];

  dbl lam_gt[6][2] = {
    {0, 0},
    {0, 0},
    {1, 0},
    {1, 0},
    {0, 1},
    {0, 1}
  };

  mesh3_s *mesh;
  mesh3_alloc(&mesh);

  utetra_s *cf;
  utetra_alloc(&cf);

  for (int i = 0; i < 6; ++i) {
    /**
     * Get vertices for this permutation
     */
    for (int j = 0; j < 3; ++j)
      memcpy(&verts_perm[3*j], &verts[3*perm[i][j]], sizeof(dbl3));

    /**
     * Create a mesh consisting of a single tetrahedron for this update
     */
    mesh3_init(mesh, verts_perm, 4, cells, 1, true, NULL);

    /**
     * Get jets for vertex data
     */
    for (int j = 0; j < 3; ++j) {
      get_gt_jet(xsrc, &verts_perm[3*j],  &jet[j]);
    }

    utetra_spec_s spec = utetra_spec_from_ptrs(mesh, jet, 3, 0, 1, 2);
    spec.tol = 1e-15;
    utetra_init(cf, &spec);

    /**
     * Verify that cost function has correct nodal values
     */
    for (int j = 0; j < 3; ++j) {
      utetra_set_lambda(cf, lam_gt[2*j]);
      dbl L = dbl3_dist(&verts_perm[3*j], &verts[9]);
      dbl f = jet[j].f + L;
      assert_that_double(utetra_get_value(cf), is_nearly_double(f));
    }

    /**
     * Do tetrahedron updates starting from random initial iterates
     */
    for (int l = 0; l < NUM_RANDOM_TRIALS; ++l) {
      lam[0] = gsl_ran_flat(rng, 0, 1);
      lam[1] = gsl_ran_flat(rng, 0, 1);
      if (lam[0] + lam[1] > 1) {
        lam[0] = 1 - lam[0];
        lam[1] = 1 - lam[1];
      }

      utetra_solve(cf, lam);
      utetra_get_lambda(cf, lam);
      utetra_get_x(cf, x);
      utetra_get_jet(cf, &newjet);

      assert_that_double(x[0], is_nearly_double(x_gt[0]));
      assert_that_double(x[1], is_nearly_double(x_gt[1]));
      assert_that_double(x[2], is_nearly_double(x_gt[2]));

      if (fabs(lam[0]) < 1e-15) {
        assert_that_double(fabs(lam_gt[i][0]), is_less_than_double(1e-15));
      } else {
        assert_that_double(lam[0], is_nearly_double(lam_gt[i][0]));
      }

      if (fabs(lam[1]) < 1e-15) {
        assert_that_double(fabs(lam_gt[i][1]), is_less_than_double(1e-15));
      } else {
        assert_that_double(lam[1], is_nearly_double(lam_gt[i][1]));
      }
    }

    mesh3_deinit(mesh);
  }

  utetra_deinit(cf);
  utetra_dealloc(&cf);

  mesh3_dealloc(&mesh);

  gsl_rng_free(rng);
}
