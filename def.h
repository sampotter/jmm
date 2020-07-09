#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NDEBUG
#define SJS_DEBUG 1
#endif

#define ROW_MAJOR_ORDERING 0
#define COLUMN_MAJOR_ORDERING 1
#define ORDERING ROW_MAJOR_ORDERING

#define EPS 1e-13
#define NO_INDEX -1
#define NO_PARENT -1
#define PI_OVER_FOUR 0.7853981633974483
#define SQRT2 1.414213562373095
#define UNFACTORED -1

typedef enum state {FAR, TRIAL, VALID, BOUNDARY, NEW_VALID} state_e;

/**
 * Change this to float to use single precision instead.
 */
typedef double dbl;

#ifdef __cplusplus
}
#endif
