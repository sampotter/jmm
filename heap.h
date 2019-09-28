#pragma once

#include "sjs.h"

typedef struct heap heap_s;

void heap_init(heap_s *heap, int capacity);
void heap_insert(heap_s *heap, int ind);
void heap_swim(heap_s *heap, int ind);
int heap_front(heap_s *heap);
void heap_pop(heap_s *heap);
int heap_size(heap_s *heap);
