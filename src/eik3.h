#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "bb.h"
#include "common.h"
#include "error.h"
#include "jet.h"
#include "par.h"

// TODO: put functions we want to go in the "public API" here

void eik3_alloc(eik3_s **eik);
void eik3_dealloc(eik3_s **eik);
void eik3_init(eik3_s *eik, mesh3_s *mesh, ftype_e ftype, eik3_s const *orig);
void eik3_deinit(eik3_s *eik);
bool eik3_is_initialized(eik3_s const *eik);

void eik3_add_pt_src_BCs(eik3_s *eik, size_t l, jet32t jet);
void eik3_add_refl_BCs(eik3_s *eik, size_t const lf[3], jet32t const jet[3]);

size_t eik3_peek(eik3_s const *eik);
jmm_error_e eik3_step(eik3_s *eik, size_t *l0);
jmm_error_e eik3_solve(eik3_s *eik);
bool eik3_is_solved(eik3_s const *eik);

bool eik3_is_far(eik3_s const *eik, size_t ind);
bool eik3_is_trial(eik3_s const *eik, size_t ind);
bool eik3_is_valid(eik3_s const *eik, size_t ind);

ftype_e eik3_get_ftype(eik3_s const *eik);
mesh3_s *eik3_get_mesh(eik3_s const *eik);
dbl eik3_get_slerp_tol(eik3_s const *eik);
dbl eik3_get_h(eik3_s const *eik);

jet32t eik3_get_jet(eik3_s const *eik, size_t l);
jet32t *eik3_get_jet_ptr(eik3_s const *eik);
state_e *eik3_get_state_ptr(eik3_s const *eik);
par3_s eik3_get_par(eik3_s const *eik, size_t l);
bool eik3_has_par(eik3_s const *eik, size_t l);
bool eik3_has_BCs(eik3_s const *eik, size_t l);
size_t const *eik3_get_accepted_ptr(eik3_s const *eik);

// TODO: "private API" functions go here...

void eik3_add_trial(eik3_s *eik, size_t l, jet32t jet);
bool eik3_is_point_source(eik3_s const *eik, size_t ind);
void eik3_set_jet(eik3_s *eik, size_t l, jet32t jet);
void eik3_set_par(eik3_s *eik, size_t l, par3_s par);
void eik3_add_diff_edge_BCs(eik3_s *eik, size_t const le[2],
                            bb31 const *T, dbl const rho1[2]);
void eik3_set_bde_bc(eik3_s *eik, size_t const le[2], bb31 const *bb);
bool eik3_get_bde_bc(eik3_s const *eik, size_t const le[2], bb31 *bb);
bool eik3_has_bde_bc(eik3_s const *eik, size_t const le[2]);
bool eik3_get_refl_bdf_inc_on_diff_edge(eik3_s const *eik, size_t const le[2],
                                        size_t lf[3]);

#ifdef __cplusplus
}
#endif
