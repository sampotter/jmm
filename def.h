#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NDEBUG
#define SJS_EIK_DEBUG 1
#endif

#define EPS 1e-13
#define NO_INDEX -1
#define NO_PARENT -1
#define SQRT2 1.414213562373095
#define UNFACTORED -1

typedef enum state {FAR, TRIAL, VALID, BOUNDARY} state_e;

/**
 * Change this to float to use single precision instead.
 */
typedef double dbl;

#ifdef __cplusplus
}
#endif
