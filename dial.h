#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"
#include "vec.h"

typedef struct dial3 dial3_s;

void dial3_add_trial(dial3_s *dial, ivec3 ind, dbl T, dvec3 grad_T);
bool dial3_step(dial3_s *dial);
void dial3_solve(dial3_s *dial);

#ifdef __cplusplus
}
#endif
