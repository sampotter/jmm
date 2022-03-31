#include "3d_wedge.h"

#include <assert.h>
#include <stdio.h>

#include <error.h>
#include <mat.h>

#include "mesh3_extra.h"

void jmm_3d_wedge_problem_alloc(jmm_3d_wedge_problem_s **wedge) {
  *wedge = malloc(sizeof(jmm_3d_wedge_problem_s));
}

jmm_error_e jmm_3d_wedge_problem_init(jmm_3d_wedge_problem_s *wedge,
                                      jmm_3d_wedge_spec_s const *spec)
{
  wedge->spec = *spec; // just copy the spec over

  jmm_error_e error = JMM_ERROR_NONE;

  mesh3_alloc(&wedge->mesh);

  error = mesh3_init_from_3d_wedge_spec(wedge->mesh, spec);
  if (error != JMM_ERROR_NONE) {
    mesh3_dealloc(&wedge->mesh);
    return error;
  }

  // we allocate these here, but we'll initialize them later in solve
  eik3_alloc(&wedge->eik_direct);
  eik3_alloc(&wedge->eik_o_refl);
  eik3_alloc(&wedge->eik_n_refl);

  return error;
}

void jmm_3d_wedge_problem_deinit(jmm_3d_wedge_problem_s *wedge) {
  mesh3_deinit(wedge->mesh);
  mesh3_dealloc(&wedge->mesh);

  if (eik3_is_initialized(wedge->eik_direct))
    eik3_deinit(wedge->eik_direct);

  if (eik3_is_initialized(wedge->eik_o_refl))
    eik3_deinit(wedge->eik_o_refl);

  if (eik3_is_initialized(wedge->eik_n_refl))
    eik3_deinit(wedge->eik_n_refl);

  eik3_dealloc(&wedge->eik_direct);
  eik3_dealloc(&wedge->eik_o_refl);
  eik3_dealloc(&wedge->eik_n_refl);
}

void jmm_3d_wedge_problem_dealloc(jmm_3d_wedge_problem_s **wedge) {
  free(*wedge);
  *wedge = NULL;
}

jmm_error_e
jmm_3d_wedge_problem_solve(jmm_3d_wedge_problem_s *wedge, dbl sp, dbl phip,
                           dbl rfac, double omega) {
  dbl3 xsrc; // the point source
  dbl3 x; // a varying mesh vertex

  /* Reset the eikonal problems if we've called solve before */

  if (eik3_is_initialized(wedge->eik_direct))
    eik3_deinit(wedge->eik_direct);

  if (eik3_is_initialized(wedge->eik_n_refl))
    eik3_deinit(wedge->eik_n_refl);

  if (eik3_is_initialized(wedge->eik_o_refl))
    eik3_deinit(wedge->eik_o_refl);

  /* Compute the location of the point source */

  xsrc[0] = sp*cos(phip);
  xsrc[1] = -sp*sin(phip);
  xsrc[2] = 0;

  /* Set up and solve the direct eikonal problem */

  eik3_init(wedge->eik_direct, wedge->mesh, FTYPE_POINT_SOURCE);

  size_t num_initialized = 0;
  for (size_t i = 0; i < mesh3_ncells(wedge->mesh); ++i) {
    size_t cv[4];
    mesh3_cv(wedge->mesh, i, cv);

    /* first, check if any of the current cell's vertices intersection
     * the initialization ball... */
    bool found_intersection = false;
    for (size_t j = 0; j < 4; ++j) {
      mesh3_copy_vert(wedge->mesh, i, x);
      if (dbl3_dist(x, xsrc) <= rfac) {
        found_intersection = true;
        break;
      }
    }
    if (!found_intersection)
      continue;

    /* ... and if they do, initialize all of the vertices in the cell
     * with exact data */
    for (size_t j = 0; j < 4; ++j) {
      if (eik3_is_trial(wedge->eik_direct, cv[j]))
        continue;

      jet32t jet;

      mesh3_copy_vert(wedge->mesh, cv[j], x);
      jet.f = dbl3_dist(x, xsrc);

      dbl3_sub(x, xsrc, jet.Df);
      dbl3_dbl_div_inplace(jet.Df, jet.f);

      dbl33 eye, Df_otimes_Df;
      dbl33_eye(eye);
      dbl3_outer(jet.Df, jet.Df, Df_otimes_Df);
      dbl33_sub(eye, Df_otimes_Df, jet.D2f);
      dbl33_dbl_div_inplace(jet.D2f, jet.f);

      eik3_add_pt_src_BCs(wedge->eik_direct, cv[j], jet);
      ++num_initialized;
    }
  }

  if (num_initialized == 0) {
    return JMM_ERROR_BAD_ARGUMENTS;
  }

  if (wedge->spec.verbose) {
    printf("Number of vertices in the initialization ball: %lu\n",
           num_initialized);
  }

  eik3_solve(wedge->eik_direct);

  (void)omega;

  return JMM_ERROR_NONE;
}
