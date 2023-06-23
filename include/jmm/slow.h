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
  STYPE_CONSTANT,
  STYPE_FUNC_PTR,
  STYPE_JET31T,
  STYPE_NUM_STYPE
} stype_e;

typedef struct sfunc {
  stype_e stype;
  struct {
    dbl (*s)(dbl3 x);
    void (*Ds)(dbl3 x, dbl3 Ds);
    void (*D2s)(dbl3 x, dbl33 D2s);
  } funcs;
  jet31t *data_jet31t;
} sfunc_s;

static sfunc_s const SFUNC_CONSTANT = {
  .stype = STYPE_CONSTANT,
  .funcs = {
    .s = NULL,
    .Ds = NULL,
    .D2s = NULL
  }
};

#ifdef __cplusplus
}
#endif
