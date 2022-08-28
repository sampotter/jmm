#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct dial3 dial3_s;

void dial3_alloc(dial3_s **dial);
error_e dial3_init(dial3_s *dial, stype_e stype, int const *shape, dbl h);
void dial3_deinit(dial3_s *dial);
void dial3_dealloc(dial3_s **dial);
void dial3_add_point_source(dial3_s *dial, int const *ind0, dbl T);
void dial3_add_boundary_points(dial3_s *dial, int const *inds, size_t n);
bool dial3_step(dial3_s *dial);
void dial3_solve(dial3_s *dial);
dbl dial3_get_T(dial3_s const *dial, int l);
void dial3_get_grad_T(dial3_s const *dial, int l, dbl *grad_T);
dbl *dial3_get_Toff_ptr(dial3_s const *dial);
dbl *dial3_get_xsrc_ptr(dial3_s const *dial);
state_e *dial3_get_state_ptr(dial3_s const *dial);

#ifdef __cplusplus
}
#endif
