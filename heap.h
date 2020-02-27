#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <Block.h>

#include "sjs_eik.h"

typedef struct heap heap_s;

typedef dbl (^value_b)(int);
typedef void (^setpos_b)(int, int);

void heap_alloc(heap_s **heap);
void heap_dealloc(heap_s **heap);
void heap_init(heap_s *heap, int capacity, value_b value, setpos_b setpos);
void heap_deinit(heap_s *heap);
void heap_insert(heap_s *heap, int ind);
void heap_swim(heap_s *heap, int ind);
int heap_front(heap_s *heap);
void heap_pop(heap_s *heap);
int heap_size(heap_s *heap);

#ifdef __cplusplus
}
#endif
