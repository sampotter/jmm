#include "pool.h"

#include <stdlib.h>

struct pool {
  void *mem;
  size_t size;
  size_t capacity;
};

void pool_alloc(pool_s **pool) {
  *pool = malloc(sizeof(pool_s));
}

void pool_dealloc(pool_s **pool) {
  free(*pool);
  *pool = NULL;
}

void pool_init(pool_s *pool, size_t capacity) {
  pool->mem = malloc(capacity);
  pool->size = 0;
  pool->capacity = capacity;
}

void pool_deinit(pool_s *pool) {
  free(pool->mem);
  pool->mem = NULL;
  pool->size = 0;
  pool->capacity = 0;
}

void pool_grow(pool_s *pool, size_t capacity) {
  if (pool->capacity >= capacity)
    return;
  pool->mem = realloc(pool->mem, capacity);
  pool->capacity = capacity;
}

void *pool_get(pool_s *pool, size_t size) {
  while (pool->size + size > pool->capacity)
    pool_grow(pool, 2*pool->capacity);
  void *ptr = pool->mem + pool->size;
  pool->size += size;
  return ptr;
}
