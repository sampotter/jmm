#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "def.h"

typedef struct heap heap_s;

typedef dbl (*value_f)(void *, int);
typedef void (*setpos_f)(void *, int, int);

void heap_alloc(heap_s **heap);
void heap_dealloc(heap_s **heap);
void heap_init(heap_s *heap, int capacity, value_f value, setpos_f setpos,
               void *context);
void heap_deinit(heap_s *heap);
void heap_insert(heap_s *heap, int ind);
void heap_swim(heap_s *heap, int ind);
int heap_front(heap_s *heap);
void heap_pop(heap_s *heap);
int heap_size(heap_s *heap);

#ifdef __cplusplus
}
#endif
