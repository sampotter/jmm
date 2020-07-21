#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "vec.h"

typedef struct dial3 dial3_s;

void dial3_alloc(dial3_s **dial);
error_e dial3_init(dial3_s *dial, stype_e stype, ivec3 shape, dbl h);
void dial3_deinit(dial3_s *dial);
void dial3_dealloc(dial3_s **dial);
void dial3_add_trial(dial3_s *dial, ivec3 ind, dbl T, dvec3 grad_T);
void dial3_add_point_source_with_trial_nbs(dial3_s *dial, ivec3 ind0, dbl T0);
bool dial3_step(dial3_s *dial);
void dial3_solve(dial3_s *dial);
dbl dial3_get_T(dial3_s const *dial, int l);
dbl *dial3_get_T_ptr(dial3_s const *dial);
dvec3 dial3_get_grad_T(dial3_s const *dial, int l);

#ifdef __cplusplus
}
#endif
