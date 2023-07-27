#pragma once

#include <stddef.h>

typedef struct pool pool_s;

void pool_alloc(pool_s **pool);
void pool_dealloc(pool_s **pool);
void pool_init(pool_s *pool, size_t initial_capacity);
void pool_deinit(pool_s *pool);
void *pool_get(pool_s *pool, size_t num_bytes);
