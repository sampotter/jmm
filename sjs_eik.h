#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "bicubic.h"
#include "jet.h"
#include "vec.h"

typedef dbl (*sfield_f)(void *, dvec2);
typedef dvec2 (*vfield_f)(void *, dvec2);

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
void sjs_init(sjs_s *sjs, ivec2 shape, dvec2 xymin, dbl h, sfield_f s,
              vfield_f grad_s, void *context);

/**
 * Deinitialize `sjs`. This should be called before `sjs_dealloc` is
 * called.
 *
 * TODO: remove after switching to GC
 */
void sjs_deinit(sjs_s *sjs);

/**
 * Run the solver for one step. See comments below about
 * `sjs_solve`. This function is just called repeatedly by `sjs_solve`
 * until there are no more `TRIAL` or `FAR` nodes.
 */
void sjs_step(sjs_s *sjs);

/**
 * Solve the problem specified by `sjs` with the given boundary
 * data. This assumes that `sjs_init` has been called, and boundary
 * data has been provided through at least one call to
 * `sjs_add_bd`. See also `sjs_step`.
 */
void sjs_solve(sjs_s *sjs);

/**
 * Set the jet at the node at `ind` to `jet`, set its state to
 * `TRIAL`, and insert it into the heap. This function should be
 * called to set the boundary data for the problem before solving it.
 */
void sjs_add_trial(sjs_s *sjs, ivec2 ind, jet_s jet);

/**
 * Set the state of the node at `ind` to `BOUNDARY`, effectively
 * removing it from the domain.
 */
void sjs_make_bd(sjs_s *sjs, ivec2 ind);

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
bicubic_s *sjs_bicubic(sjs_s *sjs, ivec2 cind);

#ifdef __cplusplus
}
#endif
