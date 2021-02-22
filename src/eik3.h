#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "jet.h"
#include "par.h"

typedef struct cutedge {
  /**
   * A cofficient in [0, 1] such that x(t) = (1 - t)*x0 + t*x1 is the
   * intersection between the shadow boundary and [x0, x1]. Here, x0
   * and x1 correspond to indices (l0, l1) where l0 < l1 (how edges
   * are stored internally).
   */
  dbl t;

  /**
   * The surface normal of the shadow boundary at x(t).
   */
  dbl n[3];
} cutedge_s;

void eik3_alloc(eik3_s **eik);
void eik3_dealloc(eik3_s **eik);
void eik3_init(eik3_s *eik, mesh3_s *mesh);
void eik3_deinit(eik3_s *eik);
size_t eik3_peek(eik3_s const *eik);
size_t eik3_step(eik3_s *eik);
void eik3_solve(eik3_s *eik);
void eik3_add_trial(eik3_s *eik, size_t ind, jet3 jet);
void eik3_add_valid(eik3_s *eik, size_t ind, jet3 jet);
bool eik3_is_point_source(eik3_s const *eik, size_t ind);
bool eik3_is_far(eik3_s const *eik, size_t ind);
bool eik3_is_trial(eik3_s const *eik, size_t ind);
bool eik3_is_valid(eik3_s const *eik, size_t ind);
bool eik3_is_shadow(eik3_s const *eik, size_t l);
mesh3_s *eik3_get_mesh(eik3_s const *eik);
jet3 eik3_get_jet(eik3_s const *eik, size_t l);
jet3 *eik3_get_jet_ptr(eik3_s const *eik);
state_e *eik3_get_state_ptr(eik3_s const *eik);
par3_s eik3_get_par(eik3_s const *eik, size_t l);
void eik3_get_DT(eik3_s const *eik, size_t l, dbl DT[3]);
edgemap_s const *eik3_get_cutset(eik3_s const *eik);

#ifdef __cplusplus
}
#endif
