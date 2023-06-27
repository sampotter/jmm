#pragma once

#include "array.h"
#include "utetra.h"

typedef struct utetra_cache utetra_cache_s;

void utetra_cache_alloc(utetra_cache_s **cache);
void utetra_cache_dealloc(utetra_cache_s **cache);
void utetra_cache_init(utetra_cache_s *cache);
void utetra_cache_deinit(utetra_cache_s *cache);
bool utetra_cache_contains_utetra(utetra_cache_s const *cache, utetra_s const *utetra);
bool utetra_cache_contains_inds(utetra_cache_s const *cache, size_t lhat, uint3 const l);
array_s *utetra_cache_pop_bracket(utetra_cache_s *cache, utetra_s const *utetra);
void utetra_cache_purge(utetra_cache_s *cache, size_t l);
bool utetra_cache_try_add_unique(utetra_cache_s *cache, utetra_s *utetra);
