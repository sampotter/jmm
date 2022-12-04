#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * An enum encoding the "type" of slowness function to be
 * used. Specifically, this encodes the manner in which slowness data
 * is made available to a solver (e.g., whether the exact slowness
 * function is available, whether the gradient is available, whether
 * just the function values are provided on a grid, etc.).
 */
typedef enum stype {
  CONSTANT,
  JET31T,
  NUM_STYPE
} stype_e;

typedef struct sfunc {
  stype_e stype;
  union {
    jet31t s_jet31t;
  };
} sfunc_s;

static sfunc_s const SFUNC_CONSTANT = {.stype = CONSTANT};

#ifdef __cplusplus
}
#endif
