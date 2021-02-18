#include "pool.h"

#include <assert.h>
#include <stdlib.h>

#include "array.h"

typedef struct block {
  void *ptr;
  size_t size;
  size_t capacity;
} block_s;

void block_init(block_s *block, size_t capacity) {
  block->ptr = malloc(capacity);
  block->size = 0;
  block->capacity = capacity;
}

void block_deinit(block_s *block) {
  free(block->ptr);
  block->ptr = NULL;
}

size_t block_get_free_bytes(block_s const *block) {
  assert(block->size <= block->capacity);
  return block->capacity - block->size;
}

void *block_get(block_s *block, size_t num_bytes) {
  assert(num_bytes <= block_get_free_bytes(block));
  void *ptr = block->ptr + block->size;
  block->size += num_bytes;
  return ptr;
}

struct pool {
  array_s *blocks;
};

void pool_alloc(pool_s **pool) {
  *pool = malloc(sizeof(pool_s));
}

void pool_dealloc(pool_s **pool) {
  free(*pool);
  *pool = NULL;
}

static void pool_append_block(pool_s *pool, size_t capacity) {
  block_s block;
  block_init(&block, capacity);
  array_append(pool->blocks, &block);
}

void pool_init(pool_s *pool, size_t initial_capacity) {
  // Set up block list
  array_alloc(&pool->blocks);
  array_init(pool->blocks, sizeof(block_s), 1);

  // Set up first block
  pool_append_block(pool, initial_capacity);
}

void pool_deinit(pool_s *pool) {
  // Free blocks in block list
  for (size_t i = 0; i < array_size(pool->blocks); ++i)
    block_deinit(array_get_ptr(pool->blocks, i));

  // Free block list
  array_deinit(pool->blocks);
  array_dealloc(&pool->blocks);
}

void *pool_get(pool_s *pool, size_t num_bytes) {
  block_s *block;

  // First, traverse the block list and see if we can allocate from
  // any of the existing blocks
  for (size_t i = 0; i < array_size(pool->blocks); ++i) {
    block = array_get_ptr(pool->blocks, i);
    if (num_bytes <= block_get_free_bytes(block))
      return block_get(block, num_bytes);
  }

  // If we failed to do so, append new blocks to the end of free list,
  // doubling the capacity of each, until we're able to
  // allocate. Doubling the capacity each step allows us to keep
  // allocation time logarithmic, which should be fine for our
  // purposes.
  do {
    pool_append_block(pool, 2*block->capacity);
    block = array_get_ptr(pool->blocks, array_size(pool->blocks) - 1);
  } while (block->capacity < num_bytes);

  assert(block->capacity >= num_bytes);
  return block_get(block, num_bytes);
}
