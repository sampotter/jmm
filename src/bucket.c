#include <jmm/bucket.h>

#include <stdlib.h>

#define INIT_BUCKET_SIZE 16

/**
 * A bucket is a one-directional queue, implemented as a ring buffer,
 * storing (l)inear indices of nodes in a Dial-like solver. It's also
 * a node in a linked list of buckets.
 */
struct bucket {
  size_t size;
  size_t start;
  size_t stop;
  size_t capacity;
  int *l;
  bucket_s *next;
};

void bucket_alloc(bucket_s **bucket) {
  *bucket = malloc(sizeof(bucket_s));
}

void bucket_init(bucket_s *bucket) {
  bucket->size = 0;
  bucket->start = 0;
  bucket->stop = 0;
  bucket->capacity = INIT_BUCKET_SIZE;
  bucket->l = malloc(sizeof(int)*INIT_BUCKET_SIZE);
  bucket->next = NULL;
}

void bucket_deinit(bucket_s *bucket) {
  free(bucket->l);
}

void bucket_dealloc(bucket_s **bucket) {
  free(*bucket);
  *bucket = NULL;
}

void bucket_grow(bucket_s *bucket) {
  int *new_l = malloc(2*sizeof(int)*bucket->capacity);
  for (size_t i = 0, j = 0; i < bucket->size; ++i) {
    new_l[i] = bucket->l[j];
    j = (j + 1) % bucket->size;
  }
  free(bucket->l);
  bucket->l = new_l;

  // Update old parameters
  bucket->start = 0;
  bucket->stop = bucket->size;
  bucket->capacity *= 2;
}

void bucket_push(bucket_s *bucket, int l) {
  if (bucket->size == bucket->capacity) {
    bucket_grow(bucket);
  }
  bucket->l[bucket->stop] = l;
  bucket->stop = (bucket->stop + 1) % bucket->capacity;
  ++bucket->size;
}

int bucket_pop(bucket_s *bucket) {
  int l = bucket->l[bucket->start];
  bucket->start = (bucket->start + 1) % bucket->capacity;
  --bucket->size;
  return l;
}

bucket_s *bucket_get_next(bucket_s const *bucket) {
  return bucket->next;
}

void bucket_set_next(bucket_s *bucket, bucket_s *next) {
  bucket->next = next;
}

size_t bucket_get_size(bucket_s const *bucket) {
  return bucket->size;
}

bool bucket_is_empty(bucket_s const *bucket) {
  return bucket->size == 0;
}
