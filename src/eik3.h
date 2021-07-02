#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "bb.h"
#include "common.h"
#include "jet.h"
#include "par.h"

void eik3_alloc(eik3_s **eik);
void eik3_dealloc(eik3_s **eik);
void eik3_init(eik3_s *eik, mesh3_s *mesh, ftype_e ftype, eik3_s const *orig);
void eik3_deinit(eik3_s *eik);
size_t eik3_peek(eik3_s const *eik);
size_t eik3_step(eik3_s *eik);
void eik3_solve(eik3_s *eik);
bool eik3_is_solved(eik3_s const *eik);
void eik3_add_trial(eik3_s *eik, size_t ind, jet3 jet);
bool eik3_is_point_source(eik3_s const *eik, size_t ind);
bool eik3_is_far(eik3_s const *eik, size_t ind);
bool eik3_is_trial(eik3_s const *eik, size_t ind);
bool eik3_is_valid(eik3_s const *eik, size_t ind);
bool eik3_is_shadow(eik3_s const *eik, size_t l);
mesh3_s *eik3_get_mesh(eik3_s const *eik);
jet3 eik3_get_jet(eik3_s const *eik, size_t l);
void eik3_set_jet(eik3_s *eik, size_t l, jet3 jet);
jet3 *eik3_get_jet_ptr(eik3_s const *eik);
dbl33 *eik3_get_hess_ptr(eik3_s const *eik);
state_e *eik3_get_state_ptr(eik3_s const *eik);
par3_s eik3_get_par(eik3_s const *eik, size_t l);
void eik3_set_par(eik3_s *eik, size_t l, par3_s par);
bool eik3_has_par(eik3_s const *eik, size_t l);
void eik3_get_DT(eik3_s const *eik, size_t l, dbl DT[3]);
dbl const *eik3_get_DT_ptr(eik3_s const *eik, size_t l);
void eik3_get_D2T(eik3_s const *eik, size_t l, dbl D2T[3][3]);
dbl *eik3_get_t_in_ptr(eik3_s const *eik);
dbl *eik3_get_t_out_ptr(eik3_s const *eik);
void eik3_add_pt_src_BCs(eik3_s *eik, size_t l, jet3 jet);
void eik3_add_refl_BCs(eik3_s *eik, size_t const lf[3], jet3 const jet[3],
                       dbl33 const hess[3], dbl const t_in[3][3]);
void eik3_add_diff_edge_BCs(eik3_s *eik, size_t const le[2],
                            bb31 const *T, dbl const rho1[2],
                            dbl3 const t_in[2]);
void eik3_set_bde_bc(eik3_s *eik, size_t const le[2], bb31 const *bb);
bool eik3_get_bde_bc(eik3_s const *eik, size_t const le[2], bb31 *bb);
bool eik3_has_bde_bc(eik3_s const *eik, size_t const le[2]);
ftype_e eik3_get_ftype(eik3_s const *eik);
dbl eik3_get_slerp_tol(eik3_s const *eik);
bool eik3_has_BCs(eik3_s const *eik, size_t l);
void eik3_transport_dbl(eik3_s const *eik, dbl *values, bool skip_filled);
void eik3_transport_dblz(eik3_s const *eik, dblz *values, bool skip_filled);
void eik3_transport_curvature(eik3_s const *eik, dbl *kappa, bool skip_filled);
dbl eik3_get_h(eik3_s const *eik);
bool eik3_get_refl_bdf_inc_on_diff_edge(eik3_s const *eik, size_t const le[2],
                                        size_t lf[3]);
size_t const *eik3_get_accepted_ptr(eik3_s const *eik);

#ifdef __cplusplus
}
#endif
