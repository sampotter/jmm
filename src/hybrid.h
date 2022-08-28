#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef dbl (*hybrid_cost_func_t)(dbl, void*);

bool hybrid(hybrid_cost_func_t f, dbl a, dbl b, void *context, dbl *t);

#ifdef __cplusplus
}
#endif
