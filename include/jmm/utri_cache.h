#pragma once

#include "utri.h"

typedef struct utri_cache utri_cache_s;

void utri_cache_alloc(utri_cache_s **cache);
void utri_cache_dealloc(utri_cache_s **cache);
void utri_cache_init(utri_cache_s *cache);
void utri_cache_deinit(utri_cache_s *cache);
bool utri_cache_contains_utri(utri_cache_s const *cache, utri_s const *utri);
bool utri_cache_contains_inds(utri_cache_s const *cache, size_t lhat, uint2 l);
utri_s *utri_cache_pop(utri_cache_s *cache, utri_s const *utri);
void utri_cache_purge(utri_cache_s *cache, size_t l);
bool utri_cache_try_add_unique(utri_cache_s *cache, utri_s *utri);
