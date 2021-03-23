#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stddef.h>

#define ARRAY_DEFAULT_CAPACITY 32

typedef struct array array_s;

void array_alloc(array_s **arr);
void array_dealloc(array_s **arr);
void array_init(array_s *arr, size_t eltsize, size_t capacity);
void array_deinit(array_s *arr);
bool array_is_empty(array_s const *arr);
size_t array_size(array_s const *arr);
size_t array_find(array_s const *arr, void const *elt);
bool array_contains(array_s const *arr, void const *elt);
void array_append(array_s *arr, void const *elt);
void array_get(array_s const *arr, size_t i, void *elt);
void *array_get_ptr(array_s const *arr, size_t i);
void array_delete(array_s *arr, size_t i);
void array_pop_front(array_s *arr, void *elt);

#ifdef __cplusplus
}
#endif
