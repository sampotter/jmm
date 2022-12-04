#include <jmm/utri_cache.h>

#include <jmm/array.h>

struct utri_cache {
  array_s *utri_arr;
};

void utri_cache_alloc(utri_cache_s **cache) {
  *cache = malloc(sizeof(utri_cache_s));
}

void utri_cache_dealloc(utri_cache_s **cache) {
  free(*cache);
  *cache = NULL;
}

void utri_cache_init(utri_cache_s *cache) {
  array_alloc(&cache->utri_arr);
  array_init(cache->utri_arr, sizeof(utri_s *), 16);
}

void utri_cache_deinit(utri_cache_s *cache) {
  utri_s *utri;
  for (size_t i = 0; i < array_size(cache->utri_arr); ++i) {
    array_get(cache->utri_arr, i, &utri);
    utri_dealloc(&utri);
  }

  array_deinit(cache->utri_arr);
  array_dealloc(&cache->utri_arr);
}

/* Check whether `utri` has been stored in the cache for
 * edge-diffracted updates already. */
bool utri_cache_contains_utri(utri_cache_s const *cache, utri_s const *utri) {
  utri_s const *utri_other = NULL;
  for (size_t i = 0; i < array_size(cache->utri_arr); ++i) {
    array_get(cache->utri_arr, i, &utri_other);
    if (utris_have_same_inds(utri, utri_other))
      return true;
  }
  return false;
}

bool utri_cache_contains_inds(utri_cache_s const *cache, size_t lhat, uint2 l) {
  for (size_t j = 0; j < array_size(cache->utri_arr); ++j) {
    utri_s *utri;
    array_get(cache->utri_arr, j, &utri);
    if (utri_has_inds(utri, lhat, l))
      return true;
  }
  return false;
}

/* Look through the cache of old edge-diffracted two-point updates for
 * an update which matches `utri` (is incident on it and which shares
 * the same active index). If we find one, delete it from the cache
 * and return it. */
utri_s *utri_cache_pop(utri_cache_s *cache, utri_s const *utri) {
  /* get indices of current `utri` */
  size_t l = utri_get_l(utri);
  size_t l_active = utri_get_active_ind(utri);
  size_t l_inactive = utri_get_inactive_ind(utri);

  utri_s *utri_other = NULL;

  /* iterate over the other `utri` in the cache... */
  for (size_t i = 0; i < array_size(cache->utri_arr); ++i) {
    array_get(cache->utri_arr, i, &utri_other);

    /* if this is a distinct `utri`, with the same target index and
     * the same active index (so, the inactive index must be
     * different!) ... */
    if (l == utri_get_l(utri_other) &&
        l_active == utri_get_active_ind(utri_other) &&
        l_inactive != utri_get_inactive_ind(utri_other)) {
      /* ... then delete it and break from the loop */
      array_delete(cache->utri_arr, i);
      break;
    }

    /* make sure to set `utri_other` back to `NULL` here---if we don't
     * end up breaking, on the last loop, we will return `NULL` to
     * signal that we didn't find a matching `utri` to delete */
    utri_other = NULL;
  }

  return utri_other;
}


/* Remove and free triangle updates targeting the node with index `l`
 * from `cache`. */
void utri_cache_purge(utri_cache_s *cache, size_t l) {
  utri_s *utri;
  for (size_t i = array_size(cache->utri_arr); i > 0; --i) {
    array_get(cache->utri_arr, i - 1, &utri);
    if (utri_get_l(utri) == l) {
      utri_dealloc(&utri);
      array_delete(cache->utri_arr, i - 1);
    }
  }
}

/* Try to add `utri` to `cache`. If `utri` is already contained, then
 * return `false` to signal failure. Otherwise, return `true`. */
bool utri_cache_try_add_unique(utri_cache_s *cache, utri_s *utri) {
  if (utri_cache_contains_utri(cache, utri))
    return false;
  array_append(cache->utri_arr, &utri);
  return true;
}
