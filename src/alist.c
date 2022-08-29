#include <jmm/alist.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

struct alist {
  char *data;
  size_t keysize;
  size_t eltsize;
  size_t nodesize;
  size_t size;
  size_t capacity;
};

void alist_alloc(alist_s **lst) {
  *lst = malloc(sizeof(alist_s));
}

void alist_dealloc(alist_s **lst) {
  assert(*lst != NULL);
  free(*lst);
  *lst = NULL;
}

void alist_init(alist_s *lst, size_t keysize, size_t eltsize, size_t capacity) {
  lst->nodesize = keysize + eltsize;
  lst->data = malloc(lst->nodesize*capacity);
  lst->keysize = keysize;
  lst->eltsize = eltsize;
  lst->size = 0;
  lst->capacity = capacity;
}

void alist_deinit(alist_s *lst) {
  assert(lst->data != NULL);
  free(lst->data);
  lst->data = NULL;
}

bool alist_is_empty(alist_s const *lst) {
  return lst->size == 0;
}

size_t alist_size(alist_s const *lst) {
  return lst->size;
}

size_t alist_find(alist_s const *lst, void const *key) {
  char *ptr = lst->data;
  size_t i = 0;
  for (i = 0; i < lst->size; ++i) {
    if (!memcmp(ptr, key, lst->keysize))
      return i;
    ptr += lst->nodesize;
  }
  assert(i == lst->size);
  return i;
}

bool alist_contains(alist_s const *lst, void const *key) {
  return alist_find(lst, key) < lst->size;
}

static void grow_if_necessary(alist_s *lst) {
  if (lst->size < lst->capacity)
    return;
  lst->capacity *= 2;
  lst->data = realloc(lst->data, lst->nodesize*lst->capacity);
}

void alist_append(alist_s *lst, void const *key, void const *elt) {
  grow_if_necessary(lst);
  char *ptr = lst->data + lst->nodesize*lst->size;
  memcpy(ptr, key, lst->keysize);
  memcpy(ptr + lst->keysize, elt, lst->eltsize);
  ++lst->size;
}

bool alist_get_by_index(alist_s const *lst, size_t i, void *elt) {
  if (i >= lst->size)
    return false;
  memcpy(elt, lst->data + lst->nodesize*i + lst->keysize, lst->eltsize);
  return true;
}

bool alist_get_by_key(alist_s const *lst, void const *key, void *elt) {
  size_t i = alist_find(lst, key);
  assert(i <= lst->size);
  if (i == lst->size)
    return false;
  memcpy(elt, lst->data + lst->nodesize*i + lst->keysize, lst->eltsize);
  return true;
}

void alist_set_by_index(alist_s *lst, size_t i, void const *elt) {
  assert(i < lst->size);
  memcpy(lst->data + lst->nodesize*i + lst->keysize, elt, lst->eltsize);
}

void alist_set_by_key(alist_s *lst, void const *key, void const *elt) {
  size_t i = alist_find(lst, key);
  assert(i <= lst->size);
  if (i == lst->size)
    alist_append(lst, key, elt);
  else
    memcpy(lst->data + lst->nodesize*i + lst->keysize, elt, lst->eltsize);
}

bool alist_remove_by_key(alist_s *lst, void const *key) {
  size_t i = alist_find(lst, key);
  assert(i <= lst->size);
  if (i == lst->size)
    return false;
  memmove(lst->data + lst->nodesize*i,
          lst->data + lst->nodesize*(i + 1),
          (lst->size - i - 1)*lst->nodesize);
  --lst->size;
  return true;
}

bool alist_get_key(alist_s const *lst, size_t i, void *key) {
  if (i >= lst->size)
    return false;
  memcpy(key, lst->data + i*lst->nodesize, lst->keysize);
  return true;
}

bool alist_get_pair(alist_s const *lst, size_t i, void *key, void *elt) {
  if (i >= lst->size)
    return false;
  if (key != NULL)
    memcpy(key, lst->data + i*lst->nodesize, lst->keysize);
  if (elt != NULL)
    memcpy(elt, lst->data + i*lst->nodesize + lst->keysize, lst->eltsize);
  return true;
}

void alist_clear(alist_s *lst) {
  lst->size = 0;
}
