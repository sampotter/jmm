#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
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
#define JMM_PI 3.141592653589793
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

typedef enum policy {
  POLICY_INVALID,
  POLICY_COPY,
  POLICY_XFER,
  POLICY_VIEW
} policy_e;

typedef enum state {
  FAR,
  TRIAL,
  VALID,
  BOUNDARY,
  ADJACENT_TO_BOUNDARY,
  NEW_VALID,
  UNKNOWN
} state_e;

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

typedef enum field {
  FIELD_A,
  FIELD_T,
  FIELD_E_T,
  FIELD_ORIGIN,
  FIELD_SPREADING
} field_e;

typedef enum order {
  ORDER_ROW_MAJOR,
  ORDER_COLUMN_MAJOR,
  NUM_ORDERS
} order_e;

typedef double dbl;
typedef double _Complex dblz;

typedef dbl dbl2[2];
typedef dbl2 dbl22[2];

typedef dbl dbl3[3];
typedef dbl3 dbl43[4];
typedef dbl3 dbl33[3];

typedef dbl dbl4[4];
typedef dbl4 dbl34[3];
typedef dbl4 dbl44[4];

typedef int int2[2];

typedef int int3[3];

typedef size_t uint3[3];

typedef enum error {
  SUCCESS,
  BAD_ARGUMENT
} error_e;

typedef int (*compar_t)(void const *, void const *);

#ifdef __cplusplus
}
#endif
