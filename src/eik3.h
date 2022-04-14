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
void eik3_init(eik3_s *eik, mesh3_s *mesh);
void eik3_deinit(eik3_s *eik);
bool eik3_is_initialized(eik3_s const *eik);

void eik3_dump_jet(eik3_s const *eik, char const *path);
void eik3_dump_state(eik3_s const *eik, char const *path);

size_t eik3_peek(eik3_s const *eik);
jmm_error_e eik3_step(eik3_s *eik, size_t *l0);
jmm_error_e eik3_solve(eik3_s *eik);
bool eik3_is_solved(eik3_s const *eik);

bool eik3_is_far(eik3_s const *eik, size_t ind);
bool eik3_is_trial(eik3_s const *eik, size_t ind);
bool eik3_is_valid(eik3_s const *eik, size_t ind);

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

void eik3_add_trial(eik3_s *eik, size_t l, jet32t jet);
void eik3_set_jet(eik3_s *eik, size_t l, jet32t jet);
void eik3_set_par(eik3_s *eik, size_t l, par3_s par);

#ifdef __cplusplus
}
#endif
