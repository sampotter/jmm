#include "heap.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

struct heap {
  int capacity;
  int size;
  int* inds;
  value_f value;
  setpos_f setpos;
  void *context;
};

void heap_alloc(heap_s **heap) {
  *heap = malloc(sizeof(heap_s));
  assert(*heap != NULL);
}

void heap_dealloc(heap_s **heap) {
  free(*heap);
  *heap = NULL;
}

void heap_init(heap_s *heap, int capacity, value_f value, setpos_f setpos,
               void *context) {
  heap->capacity = capacity;
  heap->size = 0;
  heap->inds = malloc(heap->capacity*sizeof(int));
  assert(heap->inds != NULL);
#if SJS_DEBUG
  for (int i = 0; i < heap->capacity; ++i) {
    heap->inds[i] = NO_INDEX;
  }
#endif
  heap->value = value;
  heap->setpos = setpos;
  heap->context = context;
}

void heap_deinit(heap_s *heap) {
  free(heap->inds);
  heap->inds = NULL;
}

void heap_grow(heap_s *heap) {
  heap->capacity *= 2;
  heap->inds = realloc(heap->inds, sizeof(int)*heap->capacity);
  assert(heap->inds != NULL);
#if SJS_DEBUG
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

#ifdef SJS_DEBUG
  int ind = heap->inds[pos];
  assert(ind != NO_INDEX);
#endif

  return heap->value(heap->context, heap->inds[pos]);
}

void heap_set(heap_s *heap, int pos, int ind) {
  assert(pos >= 0);
  assert(pos < heap->size);

  heap->inds[pos] = ind;
  heap->setpos(heap->context, ind, pos);
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
#if SJS_DEBUG
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
  size_t front_ind = heap->inds[0];
  heap_swap(heap, 0, heap->size - 1);
  if (--heap->size > 0) {
    heap_sink(heap, 0);
  }
  heap->setpos(heap->context, front_ind, NO_INDEX);
}

int heap_size(heap_s *heap) {
  return heap->size;
}
