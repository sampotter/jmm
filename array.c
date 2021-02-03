#include "array.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

struct array {
  char *data;
  size_t eltsize;
  size_t size;
  size_t capacity;
};

void array_alloc(array_s **arr) {
  *arr = malloc(sizeof(array_s));
}

void array_dealloc(array_s **arr) {
  assert(*arr != NULL);
  free(*arr);
  *arr = NULL;
}

void array_init(array_s *arr, size_t eltsize, size_t capacity) {
  arr->data = malloc(eltsize*capacity);
  arr->eltsize = eltsize;
  arr->size = 0;
  arr->capacity = capacity;
}

void array_deinit(array_s *arr) {
  assert(arr->data != NULL);
  free(arr->data);
  arr->data = NULL;
}

bool array_is_empty(array_s const *arr) {
  return arr->size == 0;
}

size_t array_size(array_s const *arr) {
  return arr->size;
}

size_t array_find(array_s const *arr, void const *elt) {
  char *ptr = arr->data;
  size_t i = 0;
  for (i = 0; i < arr->size; ++i) {
    if (!memcmp(ptr, elt, arr->eltsize)) {
      return i;
    }
    ptr += arr->eltsize;
  }
  assert(i == arr->size);
  return i;
}

bool array_contains(array_s const *arr, void const *elt) {
  return array_find(arr, elt) < arr->size;
}

static void grow_if_necessary(array_s *arr) {
  if (arr->size < arr->capacity) {
    return;
  }
  arr->capacity *= 2;
  arr->data = realloc(arr->data, arr->eltsize*arr->capacity);
}

void array_append(array_s *arr, void const *elt) {
  grow_if_necessary(arr);
  void *ptr = arr->data + arr->eltsize*arr->size;
  memcpy(ptr, elt, arr->eltsize);
  ++arr->size;
}

void array_get(array_s const *arr, size_t i, void *elt) {
  if (i >= arr->size) {
    return;
  }
  memcpy(elt, arr->data + arr->eltsize*i, arr->eltsize);
}

void *array_get_ptr(array_s const *arr, size_t i) {
  return i < arr->size ? arr->data + arr->eltsize*i : NULL;
}
