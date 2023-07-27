#pragma once

#include <stdbool.h>
#include <stddef.h>

typedef struct bucket bucket_s;

void bucket_alloc(bucket_s **bucket);
void bucket_init(bucket_s *bucket);
void bucket_deinit(bucket_s *bucket);
void bucket_dealloc(bucket_s **bucket);
void bucket_grow(bucket_s *bucket);
void bucket_push(bucket_s *bucket, int l);
int bucket_pop(bucket_s *bucket);
bucket_s *bucket_get_next(bucket_s const *bucket);
void bucket_set_next(bucket_s *bucket, bucket_s *next);
size_t bucket_get_size(bucket_s const *bucket);
bool bucket_is_empty(bucket_s const *bucket);
