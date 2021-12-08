#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "grid2.h"
#include "jet.h"
#include "par.h"

void eik2g1_alloc(eik2g1_s **eik);
void eik2g1_dealloc(eik2g1_s **eik);
void eik2g1_init(eik2g1_s *eik, grid2_s const *grid);
void eik2g1_deinit(eik2g1_s *eik);
size_t eik2g1_peek(eik2g1_s const *eik);
size_t eik2g1_step(eik2g1_s *eik);
void eik2g1_solve(eik2g1_s *eik);
bool eik2g1_is_solved(eik2g1_s const *eik);
void eik2g1_add_trial(eik2g1_s *eik, int2 const ind, jet22t jet);
void eik2g1_add_valid(eik2g1_s *eik, int2 const ind, jet22t jet);
bool eik2g1_is_valid(eik2g1_s const *eik, int2 const ind);
bool eik2g1_is_trial(eik2g1_s const *eik, int2 const ind);
bool eik2g1_is_far(eik2g1_s const *eik, int2 const ind);
state_e const *eik2g1_get_state_ptr(eik2g1_s const *eik);
jet22t const *eik2g1_get_jet_ptr(eik2g1_s const *eik);
bool eik2g1_has_par(eik2g1_s const *eik, int2 const ind);
par2_s eik2g1_get_par(eik2g1_s const *eik, int2 const ind);
par2_s const *eik2g1_get_par_ptr(eik2g1_s const *eik);

typedef struct eik2g1_sol_info {
  dbl lam_T, lam_T_check, lam_tau, lam_star;
  dbl FT_lamT, Ftau_lamT, Ftau_lamtau;
  dbl E0, That, E0_check0, E0_check1, E0_check2, E0_check3;
  dbl E0_term1, E0_term2, E0_term3;
  dbl E0_term1_v2;
  dbl T0_error, T1_error, DT0_error, DT1_error;
} eik2g1_sol_info_s;

eik2g1_sol_info_s eik2g1_get_sol_info(eik2g1_s const *eik, int2 const ind);
