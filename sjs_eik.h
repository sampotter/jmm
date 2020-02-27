#pragma once

#ifdef __cplusplus
extern "C" {
#endif

struct sjs;

#include "hermite.h"
#include "vec.h"

typedef dbl (*sfield)(dvec2);
typedef dvec2 (*vfield)(dvec2);

typedef struct sjs sjs_s;

/**
 * Allocate memory for `sjs`.
 */
void sjs_alloc(sjs_s **sjs);

/**
 * Deallocate memory used by `sjs`.
 *
 * TODO: remove after switching to GC
 */
void sjs_dealloc(sjs_s **sjs);

/**
 * Initialize the static jet scheme workspace `sjs`. Assumes that
 * `sjs` has already been allocated. The domain is a `shape.i`
 * $\times$ `shape.j` rectangular grid. The bottom left corner is
 * given by `xymin`, and the node spacing in each direction is
 * `h`. Additionally, the slowness and its gradient are provided by
 * `s` and `grad_s` respectively.
 */
void sjs_init(sjs_s *sjs, ivec2 shape, dvec2 xymin, dbl h, sfield s, vfield grad_s);

/**
 * Deinitialize `sjs`. This should be called before `sjs_dealloc` is
 * called.
 *
 * TODO: remove after switching to GC
 */
void sjs_deinit(sjs_s *sjs);

/**
 * Solve the problem specified by `sjs` with the given boundary
 * data. This assumes that `sjs_init` has been called, and boundary
 * data has been provided through at least one call to `sjs_add_bd`.
 */
void sjs_solve(sjs_s *sjs);

/**
 * Get $T(x, y)$ (the value of `T` at the possibly non-grid-aligned
 * point specified by `xy`). The value will be interpolated from
 * whatever cell `xy` happens to lie in.
 */
dbl sjs_T(sjs_s *sjs, dvec2 xy);

/**
 * Get $T_x(x, y)$ (see `sjs_T` for more details).
 */
dbl sjs_Tx(sjs_s *sjs, dvec2 xy);

/**
 * Get $T_y(x, y)$ (see `sjs_T` for more details).
 */
dbl sjs_Ty(sjs_s *sjs, dvec2 xy);

/**
 * Get $T_{xy}(x, y)$ (see `sjs_T` for more details).
 */
dbl sjs_Txy(sjs_s *sjs, dvec2 xy);

/**
 * Get the coefficients of the bicubic corresponding to the cell
 * indexed by `cind`.
 */
bicubic *sjs_bicubic(sjs_s *sjs, ivec2 cind);

#ifdef __cplusplus
}
#endif
