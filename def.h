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
 * An enum encoding the "type" of slowness function to be
 * used. Specifically, this encodes the manner in which slowness data
 * is made available to a solver (e.g., whether the exact slowness
 * function is available, whether the gradient is available, whether
 * just the function values are provided on a grid, etc.).
 */
typedef enum stype {
  CONSTANT,
  NUM_STYPE
} stype_e;

/**
 * Change this to float to use single precision instead.
 */
typedef double dbl;

#ifdef __cplusplus
}
#endif
