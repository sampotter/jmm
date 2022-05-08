#pragma once

#include <eik3.h>
#include <error.h>
#include <grid2.h>
#include <mesh3.h>

/* A structure containing parameters which specify a 3D wedge
 * scattering problem. */
typedef struct jmm_3d_wedge_spec {
  bool verbose;
  bool visualize;

  /* Maximum tetrahedron volume constraint for TetGen */
  double maxvol;

  /* UTD wedge angle parameter. Defined so that if theta is the
   * dihedral angle of the wedge, then n = theta/pi. We assume that 0
   * < n < 1, so that the wedge is acute. */
  double n;

  /* The width of the bounding box of the computational domain: the
   * diffracting edge is the central vertical axis of a box with
   * horizontal sides of length w. */
  double w;

  /* The height of the domain. */
  double h;

  /* The reflection coefficient, assuming sound-hard reflection. */
  double R;
} jmm_3d_wedge_spec_s;

typedef struct jmm_3d_wedge_problem {
  jmm_3d_wedge_spec_s spec;

  /* Tetrahedron mesh discretizing the domain */
  mesh3_s *mesh;

  /* The eikonal problem corresponding to the direct wavefront
   * originated by the point source */
  eik3_s *eik_direct;

  /* The eikonal problem for the "o-face" reflection. */
  eik3_s *eik_o_refl;

  /* The eikonal problem for the "n-face" reflection. */
  eik3_s *eik_n_refl;

  /* Cell-averaged Hessians: */
  dbl33 *D2T_direct;
  dbl33 *D2T_o_refl;
  dbl33 *D2T_n_refl;

  /* Approximate amplitudes: */
  dbl *A_direct;
  dbl *A_o_refl;
  dbl *A_n_refl;

  /* Groundtruth data: */
  jet32t *jet_direct_gt;
  jet32t *jet_o_refl_gt;
  jet32t *jet_n_refl_gt;

  /* Approximate origin: */
  dbl *origin_direct;
  dbl *origin_o_refl;
  dbl *origin_n_refl;

  /* t_in vectors: */
  dbl3 *t_in_direct;
  dbl3 *t_in_o_refl;
  dbl3 *t_in_n_refl;

  /* t_out vectors: */
  dbl3 *t_out_direct;
  dbl3 *t_out_o_refl;
  dbl3 *t_out_n_refl;
} jmm_3d_wedge_problem_s;

void jmm_3d_wedge_problem_alloc(jmm_3d_wedge_problem_s **wedge);
jmm_error_e jmm_3d_wedge_problem_init(jmm_3d_wedge_problem_s *wedge,
                                      jmm_3d_wedge_spec_s const *spec);
void jmm_3d_wedge_problem_deinit(jmm_3d_wedge_problem_s *wedge);
void jmm_3d_wedge_problem_dealloc(jmm_3d_wedge_problem_s **wedge);
jmm_error_e jmm_3d_wedge_problem_solve(jmm_3d_wedge_problem_s *wedge,
                                       dbl sp, dbl phip,
                                       dbl rfac, double omega);
void jmm_3d_wedge_problem_dump(jmm_3d_wedge_problem_s *wedge, char const *path,
                               bool dump_direct, bool dump_o_face,
                               bool dump_n_face);
void jmm_3d_wedge_problem_save_slice_plots(jmm_3d_wedge_problem_s const *wedge,
                                           char const *path,
                                           bool dump_direct,
                                           bool dump_o_face,
                                           bool dump_n_face,
                                           grid2_s const *grid);
