#include "heap.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

typedef struct heap {
  int capacity;
  int size;
  int* inds;
  dbl (*value)(int);
  void (*set_pos)(int, int);
} heap_s;

void heap_init(heap_s *heap, int capacity) {
  heap->capacity = capacity;
  heap->size = 0;
  heap->inds = malloc(heap->capacity*sizeof(int));
  assert(heap->inds != NULL);
#ifndef NDEBUG
  for (int i = 0; i < heap->capacity; ++i) {
    heap->inds[i] = NO_INDEX;
  }
#endif
}

void heap_grow(heap_s *heap) {
  heap->capacity *= 2;
  heap->inds = realloc(heap->inds, heap->capacity);
  assert(heap->inds != NULL);
#ifndef NDEBUG
  for (int i = heap->size; i < heap->capacity; ++i) {
    heap->inds[i] = NO_INDEX;
  }
#endif
}

int left(int pos) {
  return 2*pos + 1;
}

int right(int pos) {
  return 2*pos + 2;
}

int parent(int pos) {
  return (pos - 1)/2;
}

dbl value(heap_s *heap, int pos) {
  assert(pos >= 0);
  assert(pos < heap->size);

  int ind = heap->inds[pos];
  assert(ind != NO_INDEX);

  return heap->value(heap->inds[pos]);
  // return heap->sjs->jets[heap->inds[pos]].f;
}

void heap_set(heap_s *heap, int pos, int ind) {
  assert(pos >= 0);
  assert(pos < heap->size);

  heap->inds[pos] = ind;
  heap->set_pos(ind, pos);
}

void heap_swap(heap_s *heap, int pos1, int pos2) {
  assert(pos1 >= 0);
  assert(pos1 < heap->size);
  assert(pos2 >= 0);
  assert(pos2 < heap->size);

  int tmp = heap->inds[pos1];
  heap->inds[pos1] = heap->inds[pos2];
  heap->inds[pos2] = tmp;

  heap_set(heap, pos1, heap->inds[pos1]);
  heap_set(heap, pos2, heap->inds[pos2]);
}

// TODO: this calls `value` and `heap_set` about 2x as many times as
// necessary
void heap_swim(heap_s *heap, int pos) {
  assert(pos >= 0);
  assert(pos < heap->size);

  int par = parent(pos);
  while (pos > 0 && value(heap, par) > value(heap, pos)) {
    heap_swap(heap, par, pos);
    pos = par;
    par = parent(pos);
  }
}

void heap_insert(heap_s *heap, int ind) {
  if (heap->size == heap->capacity) {
    heap_grow(heap);
  }

  int pos = heap->size++;
  heap_set(heap, pos, ind);
  heap_swim(heap, pos);
}

int heap_front(heap_s *heap) {
#ifndef NDEBUG
  int ind = heap->inds[0];
  return ind;
#else
  return heap->inds[0];
#endif
}

void heap_sink(heap_s *heap, int pos) {
  assert(pos >= 0);
  assert(pos < heap->size);

  int ch = left(pos), next = ch + 1, n = heap->size;
  dbl cval, nval;
  while (ch < n) {
    cval = value(heap, ch);
    if (next < n) {
      nval = value(heap, next);
      if (cval > nval) {
        ch = next;
        cval = nval;
      }
    }
    if (value(heap, pos) > cval) {
      heap_swap(heap, pos, ch);
    }
    pos = ch;
    ch = left(pos);
    next = ch + 1;
  }
}

void heap_pop(heap_s *heap) {
#ifndef NDEBUG
  heap->set_pos(heap->inds[0], NO_INDEX);
  // heap->sjs->positions[heap->inds[0]] = NO_INDEX;
#endif
  heap_swap(heap, 0, heap->size - 1);
  if (--heap->size > 0) {
    heap_sink(heap, 0);
  }
}

int heap_size(heap_s *heap) {
  return heap->size;
}
