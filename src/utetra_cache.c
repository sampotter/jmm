#include <jmm/utetra_cache.h>

struct utetra_cache {
  array_s *utetra_arr;
};

void utetra_cache_alloc(utetra_cache_s **cache) {
  *cache = malloc(sizeof(utetra_cache_s));
}

void utetra_cache_dealloc(utetra_cache_s **cache) {
  free(*cache);
  *cache = NULL;
}

void utetra_cache_init(utetra_cache_s *cache) {
  array_alloc(&cache->utetra_arr);
  array_init(cache->utetra_arr, sizeof(utetra_s *), 16);
}

void utetra_cache_deinit(utetra_cache_s *cache) {
  utetra_s *utetra;
  for (size_t i = 0; i < array_size(cache->utetra_arr); ++i) {
    array_get(cache->utetra_arr, i, &utetra);
    utetra_dealloc(&utetra);
  }

  array_deinit(cache->utetra_arr);
  array_dealloc(&cache->utetra_arr);
}

bool utetra_cache_contains_utetra(utetra_cache_s const *cache, utetra_s const *utetra) {
  utetra_s const *utetra_other = NULL;
  for (size_t i = 0; i < array_size(cache->utetra_arr); ++i) {
    array_get(cache->utetra_arr, i, &utetra_other);
    if (utetras_have_same_inds(utetra, utetra_other))
      return true;
  }
  return false;
}

bool utetra_cache_contains_inds(utetra_cache_s const *cache, size_t lhat, uint3 const l) {
  for (size_t j = 0; j < array_size(cache->utetra_arr); ++j) {
    utetra_s *utetra;
    array_get(cache->utetra_arr, j, &utetra);
    if (utetra_has_inds(utetra, lhat, l))
      return true;
  }
  return false;
}

array_s *utetra_cache_pop_bracket(utetra_cache_s *cache, utetra_s const *utetra) {
  size_t l = utetra_get_l(utetra);

  /* Array containing matched bracket utetra */
  array_s *utetras;
  array_alloc(&utetras);
  array_init(utetras, sizeof(utetra_s *), ARRAY_DEFAULT_CAPACITY);

  /* Array containing their indices */
  array_s *i_arr;
  array_alloc(&i_arr);
  array_init(i_arr, sizeof(size_t), ARRAY_DEFAULT_CAPACITY);

  /* First, find the indices of the cached utetra which share the same
   * target node and have the same active indices as `utetra`. */
  for (size_t i = 0; i < array_size(cache->utetra_arr); ++i) {
    utetra_s const *utetra_other;
    array_get(cache->utetra_arr, i, &utetra_other);
    if (l != utetra_get_l(utetra_other)
        || !utetras_have_same_minimizer(utetra, utetra_other))
      continue;
    array_append(i_arr, &i);
    array_append(utetras, &utetra_other);
  }

  /* If the utetras bracket the ray, we evict them from the cache
   * using the index array. If they don't, we clean up the utetra
   * array now. */
  if (utetra_is_bracketed_by_utetras(utetra, utetras)) {
    array_delete_all(cache->utetra_arr, i_arr);
  } else {
    array_deinit(utetras);
    array_dealloc(&utetras);
  }
  array_deinit(i_arr);
  array_dealloc(&i_arr);

  /* If we cleaned up utetras, it equals `NULL` now. So we return
   * `NULL` if there isn't a bracket, and the array of utetras
   * otherwise. */
  return utetras;
}

void utetra_cache_purge(utetra_cache_s *cache, size_t l) {
  utetra_s *utetra;
  for (size_t i = array_size(cache->utetra_arr); i > 0; --i) {
    array_get(cache->utetra_arr, i - 1, &utetra);
    if (utetra_get_l(utetra) == l) {
      utetra_dealloc(&utetra);
      array_delete(cache->utetra_arr, i - 1);
    }
  }
}

bool utetra_cache_try_add_unique(utetra_cache_s *cache, utetra_s *utetra) {
  if (utetra_cache_contains_utetra(cache, utetra))
    return false;
  array_append(cache->utetra_arr, &utetra);
  return true;
}
