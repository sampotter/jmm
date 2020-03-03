#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

dbl hybrid(dbl (*f)(dbl, void *), dbl a, dbl b, void *context);

#ifdef __cplusplus
}
#endif
