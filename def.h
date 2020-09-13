#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stddef.h>

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
#define SQRT5 2.236067977499789
#define SQRT10 3.1622776601683795
#define SQRT13 3.605551275463989
#define SQRT17 4.123105625617661
#define UNFACTORED -1

// TODO: replace state_e with the commented out section of code
// below. Apart from saving a significant amount of memory, this will
// allow us to view the states as characters from Python by setting
// the memory view type appropriately (hopefully...)

typedef enum state {
  FAR,
  TRIAL,
  VALID,
  BOUNDARY,
  ADJACENT_TO_BOUNDARY,
  NEW_VALID,
} state_e;

// #define FAR 'f'
// #define TRIAL 't'
// #define VALID 'v'
// #define BOUNDARY 'b'
// #define NEW_VALID 'n'

// typedef char state;

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

#ifdef __cplusplus
}
#endif
