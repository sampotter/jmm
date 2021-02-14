#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/**
 * A simple linear memory pool, only supporting the acquisition of
 * variable-sized blocks of memory. The pattern of use is:
 * - allocate pool
 * - burn through many small blocks of memory
 * - release pool (& all blocks along with it)
 */
typedef struct pool pool_s;

void pool_alloc(pool_s **pool);
void pool_dealloc(pool_s **pool);
void pool_init(pool_s *pool, size_t capacity);
void pool_deinit(pool_s *pool);
void pool_grow(pool_s *pool, size_t capacity);
void *pool_get(pool_s *pool, size_t size);

#ifdef __cplusplus
}
#endif
