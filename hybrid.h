#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef dbl (*hybrid_cost_func_t)(dbl, void*);

dbl hybrid(hybrid_cost_func_t f, dbl a, dbl b, void *context);

#ifdef __cplusplus
}
#endif
