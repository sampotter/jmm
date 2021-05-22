#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifndef NDEBUG
#define SJS_DEBUG 1
#endif

#define ROW_MAJOR_ORDERING 0
#define COLUMN_MAJOR_ORDERING 1
#define ORDERING ROW_MAJOR_ORDERING

#define EPS 1e-13
#define PI 3.141592653589793
#define PI_OVER_FOUR 0.7853981633974483
#define SQRT2 1.414213562373095
#define SQRT3 1.7320508075688772
#define SQRT5 2.236067977499789
#define SQRT10 3.1622776601683795
#define SQRT13 3.605551275463989
#define SQRT17 4.123105625617661

#define NO_INDEX -1
#define NO_LABEL SIZE_MAX
#define NO_PARENT SIZE_MAX

typedef enum state {
  FAR,
  TRIAL,
  VALID,
  BOUNDARY,
  ADJACENT_TO_BOUNDARY,
  NEW_VALID,
  UNKNOWN
} state_e;

typedef enum ftype {
  FTYPE_POINT_SOURCE,
  FTYPE_REFLECTION,
  FTYPE_EDGE_DIFFRACTION
} ftype_e;

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

typedef enum error {
  SUCCESS,
  BAD_ARGUMENT
} error_e;

typedef int (*compar_t)(void const *, void const *);

#ifdef __cplusplus
}
#endif
