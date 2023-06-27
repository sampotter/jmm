#pragma once

#include <stdbool.h>
#include <stddef.h>

typedef struct alist alist_s;

void alist_alloc(alist_s **lst);
void alist_dealloc(alist_s **lst);
void alist_init(alist_s *lst, size_t keysize, size_t eltsize, size_t capacity);
void alist_deinit(alist_s *lst);
bool alist_is_empty(alist_s const *lst);
size_t alist_size(alist_s const *lst);
size_t alist_find(alist_s const *lst, void const *key);
bool alist_contains(alist_s const *lst, void const *key);
void alist_append(alist_s *lst, void const *key, void const *elt);
bool alist_get_by_index(alist_s const *lst, size_t i, void *elt);
bool alist_get_by_key(alist_s const *lst, void const *key, void *elt);
void alist_set_by_index(alist_s *lst, size_t i, void const *elt);
void alist_set_by_key(alist_s *lst, void const *key, void const *elt);
bool alist_remove_by_key(alist_s *lst, void const *key);
bool alist_get_key(alist_s const *lst, size_t i, void *key);
bool alist_get_pair(alist_s const *lst, size_t i, void *key, void *elt);
void alist_clear(alist_s *lst);
