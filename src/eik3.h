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
void eik3_init(eik3_s *eik, mesh3_s *mesh);
void eik3_deinit(eik3_s *eik);
size_t eik3_peek(eik3_s const *eik);
size_t eik3_step(eik3_s *eik);
void eik3_solve(eik3_s *eik);
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
state_e *eik3_get_state_ptr(eik3_s const *eik);
par3_s eik3_get_par(eik3_s const *eik, size_t l);
void eik3_set_par(eik3_s *eik, size_t l, par3_s par);
bool eik3_has_par(eik3_s const *eik, size_t l);
void eik3_get_DT(eik3_s const *eik, size_t l, dbl DT[3]);
bool eik3_is_refl_bdf(eik3_s const *eik, size_t const l[3]);
dbl *eik3_get_t0_ptr(eik3_s const *eik);
void eik3_add_valid_bde(eik3_s *eik, size_t const le[2], jet3 const jet[2]);
void eik3_set_bde_bc(eik3_s *eik, size_t const le[2], bb31 const *bb);
bool eik3_get_bde_bc(eik3_s const *eik, size_t const le[2], bb31 *bb);
bool eik3_has_bde_bc(eik3_s const *eik, size_t const le[2]);

#ifdef __cplusplus
}
#endif
