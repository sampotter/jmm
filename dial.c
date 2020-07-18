#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dial.h"
#include "index.h"

#define INIT_BUCKET_SIZE 16

typedef struct bucket bucket_s;

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

void bucket_init(bucket_s *bucket) {
  bucket->size = 0;
  bucket->start = 0;
  bucket->stop = 0;
  bucket->capacity = INIT_BUCKET_SIZE;
  bucket->l = malloc(sizeof(int)*INIT_BUCKET_SIZE);
  bucket->next = NULL;
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

typedef dbl (*update_f)(dial3_s const *, int /* l0 */, int /* l */);

struct dial3 {
  stype_e stype;
  ivec3 shape;
  size_t size;
  dbl h;
  dbl gap;

  dbl *T;
  dvec3 *grad_T;

  /**
   * The state of each node.
   *
   * TODO: eventually, we should be able to just use a boolean for
   * this---or, better yet, delete this array of values entirely.
   */
  state_e *state;

  /**
   * The (l)inear (b)ucket index of each node. Either NO_INDEX, or a
   * bucket index, as computed by `dial3_bucket_T`.
   *
   * TODO: what invariants are enjoyed by this array?
   * - if node `l` is in bucket `lb`, then `nb[l] == lb`
   * - if it isn't in a bucket, we don't care?
   *
   * TODO: one issue here is that the semantics aren't 100%
   * well-defined, since a node can appear in multiple buckets. Is it
   * a problem that we only have one lb per node? Can we get away with
   * this, or will it be necessary to keep a list?
   *
   * TODO: a more provocative question than the previous TODO---is
   * this array of back-pointers necessary at all, in light of the
   * fact that we just re-add nodes to buckets?
   *
   * This is the Dial version of the "heap back-pointer" required by
   * Dijkstra-like algorithms.
   */
  int *lb;

  /**
   * The (l)inear (b)ucket index of the first bucket (pointed to by
   * `first`). The buckets form a linked list. This means that if a
   * node has index `lb` (where `lb >= lb0`), then we can retrieve the
   * corresponding bucket by traversing `lb - lb0` buckets starting
   * from `first`.
   *
   * After calling `dial3_init` and before adding any boundary data,
   * `lb0 == NO_INDEX` holds.
   */
  int lb0;

  /**
   * Linear offsets of six nearest neighbors in 3D.
   *
   * If `li` is `l`'s `i`th cardinal neighbor (with the ordering of
   * these neighbors left unspecified), then `l + nb_dl[i] = li`.
   */
  int nb_dl[6];

  /**
   * The first bucket to be processed. This is the bucket with
   * (l)inear (b)ucket index `lb0`.
   *
   * After calling `dial3_init` and before adding boundary data,
   * `first == NULL`.
   */
  bucket_s *first;

  update_f update;
};

dvec3 dial3_x(dial3_s const *dial, int l) {
  return ivec3_dbl_mul(l2ind3(dial->shape, l), dial->h);
}

/**
 * Compute a new value for the node at `l` from parent node `l0`.
 * This function modifies `dial->grad_T[l]`!
 *
 * Right now, the way this works is as follows:
 *
 * - a (possibly "virtual") point source, `xsrc` is located based on
 *   `dial->T[l0]`, `dial->grad_T[l0]`, and the location of node `l0`
 *
 * - the node `l0` is projected along the straight line connecting
 *   `l0` and `xsrc` onto the plane centered at `x0` and with normal
 *   vector `(x - x0)/h`, yielding `xs`
 *
 * - if the max norm distance between `x0` and `xs` is no more than `h`,
 *   new values of `dial->T[l]` and `dial->grad_T[l0]` are computed
 *
 * The idea here is to project the point along the characteristic onto
 * the "update box" surrounding `x`.
 *
 * Some caveats:
 *
 * 1. this does not take boundaries into consideration!
 * 2. we do *not* check if `xs` is incident on a set of `VALID` nodes!
 * 3. this only works for the "s = 1" case!
 *
 * for this experiment, we will definitely address point #1, and
 * *possibly* point #2 if proves to be necessary.
 *
 * We will not address #3 in this experiment, although this idea
 * should be able to be extended to more general slowness
 * functions. In that case, we would need to use e.g. RK integration
 * to track the characteristic back and find where it crosses the
 * plane, check what set of nodes it falls on, and compute a new value
 * of T from those nodes. There may also be some minimization
 * involved, since we no longer have a good initial guess for the
 * tangent vector... The goal would be to avoid this, too!
 *
 * TODO: looks like the l0 parameter could be replaced by x0 and xsrc
 */
dbl dial3_update_constant(dial3_s const *dial, int l0, int l) {
  // use gradient to compute spherical approximation to wavefront
  // passing through l0

  // TODO: these lines are being done for each neighborhod when they
  // should just be done once for all neighbors!
  dvec3 x0 = dial3_x(dial, l0);

// TODO: for s = 1, this will already be normalized! nice!
  dvec3 t0 = dvec3_normalized(dial->grad_T[l0]);

  dbl T0 = dial->T[l0];

  dvec3 xsrc = dvec3_saxpy(-T0, t0, x0);

  // do the projection
  dvec3 x = dial3_x(dial, l);
  dvec3 dx = dvec3_sub(x, x0);
  dvec3 t = dvec3_sub(x, xsrc);
  dbl s = dvec3_dot(dx, dvec3_sub(x0, xsrc))/dvec3_dot(dx, t);
  dvec3 xs = dvec3_saxpy(s, t, xsrc);

  // TODO: Right now we'll just try something really dumb---don't
  // check *where* xs lands, just check if it lands inside the box. If
  // it does, compute a new value and return it.

  dbl h = dial->h;
  if (dvec3_maxdist(xs, x0) > h + EPS) {
    return INFINITY;
  } else {
    dbl T = dvec3_norm(t);
    dial->grad_T[l] = dvec3_dbl_div(t, T);
    return T;
  }
}

update_f update_functions[NUM_STYPE] = {
  dial3_update_constant // CONSTANT
};

void dial3_init(dial3_s *dial, stype_e stype, ivec3 shape, dbl h) {
  dial->stype = stype;
  dial->shape = shape;
  dial->size = ivec3_prod(shape);
  dial->h = h;
  dial->gap = h*sqrt(3);

  dial->T = malloc(sizeof(dbl)*dial->size);
  for (size_t i = 0; i < dial->size; ++i) {
    dial->T[i] = INFINITY;
  }

  // grad_T is initialized with garbage
  dial->grad_T = malloc(sizeof(dvec3)*dial->size);

  dial->state = malloc(sizeof(state_e)*dial->size);
  for (size_t i = 0; i < dial->size; ++i) {
    dial->state[i] = FAR;
  }

  dial->lb = malloc(sizeof(int)*dial->size);
  for (size_t i = 0; i < dial->size; ++i) {
    dial->lb[i] = NO_INDEX;
  }

  dial->lb0 = NO_INDEX;

  // TODO: want to make sure `nb_dl` is in sorted order? (for cache
  // friendliness?)
#if ORDERING == ROW_MAJOR
  dial->nb_dl[0] = ind2l3(dial->shape, (ivec3) {.i = -1, .j =  0, .k =  0});
  dial->nb_dl[1] = ind2l3(dial->shape, (ivec3) {.i =  0, .j = -1, .k =  0});
  dial->nb_dl[2] = ind2l3(dial->shape, (ivec3) {.i =  0, .j =  0, .k = -1});
  dial->nb_dl[3] = ind2l3(dial->shape, (ivec3) {.i =  0, .j =  0, .k =  1});
  dial->nb_dl[4] = ind2l3(dial->shape, (ivec3) {.i =  0, .j =  1, .k =  0});
  dial->nb_dl[5] = ind2l3(dial->shape, (ivec3) {.i =  1, .j =  0, .k =  0});
#else
#  error not implemented yet
#endif

  dial->first = NULL;

  dial->update = update_functions[dial->stype];
}

int dial3_bucket_T(dial3_s const *dial, dbl T) {
  return T/dial->gap;
}

void dial3_prepend_buckets(dial3_s *dial, int lb) {
  while (lb < dial->lb0) {
    bucket_s *bucket = malloc(sizeof(bucket_s));
    bucket_init(bucket);
    bucket->next = dial->first;
    dial->first = bucket;
    --dial->lb0;
  }
}

bucket_s *dial3_find_bucket(dial3_s *dial, int lb) {
  if (lb < dial->lb0) {
    dial3_prepend_buckets(dial, lb);
  }
  bucket_s *bucket = dial->first;
  while (lb > dial->lb0) {
    if (bucket->next == NULL) {
      bucket->next = malloc(sizeof(bucket_s));
      bucket_init(bucket->next);
    }
    bucket = bucket->next;
    --lb;
  }
  assert(bucket != NULL);
  return bucket;
}

void dial3_insert(dial3_s *dial, int l, dbl T) {
  if (dial->first == NULL) {
    dial->lb0 = dial3_bucket_T(dial, T);
    dial->first = malloc(sizeof(bucket_s));
    bucket_init(dial->first);
    bucket_push(dial->first, l);
  } else {
    int lb = dial3_bucket_T(dial, T);
    bucket_s *bucket = dial3_find_bucket(dial, lb);
    bucket_push(bucket, l);
  }
}

void dial3_set_T(dial3_s *dial, int l, dbl T) {
  // TODO: for now, we just assume that T >= 0. In the future, we can
  // be more sophisticated about this, and introduce a `Tmin`
  // parameter. But then we'll need to adjust `dial3_bucket_T` to
  // ensure that the [Tmin, Tmin + gap) bucket corresponds corresponds
  // to `lb == 0` (so that `NO_INDEX == -1` still works).
  assert(T >= 0);
  dial->T[l] = T;
}

/**
 * TODO: this is very simple-minded for now... Later, we'll need to
 * account for different buckets, but for a point source this should
 * be fine
 */
void dial3_add_trial(dial3_s *dial, ivec3 ind, dbl T, dvec3 grad_T) {
  int l = ind2l3(dial->shape, ind);
  dial3_set_T(dial, l, T);
  dial->grad_T[l] = grad_T;
  dial3_insert(dial, l, T);
}

void dial3_update_nb(dial3_s *dial, int l0, int l) {
  dbl T = dial->update(dial, l0, l);
  // TODO: may want to add a little tolerance here to ensure we don't
  // mess with grad_T too much?
  if (T < dial->T[l]) {
    dial3_set_T(dial, l, T);
    int lb = dial3_bucket_T(dial, T);
    if (lb != dial->lb[l]) {
      assert(dial->lb[l] == NO_INDEX || (dial->lb0 <= lb && lb < dial->lb[l]));
      bucket_s *bucket = dial3_find_bucket(dial, lb);
      bucket_push(bucket, l);
      // TODO: should I actually do this?
      dial->lb[l] = lb;
    }
  }
}

void dial3_update_nbs(dial3_s *dial, int l0) {
  for (int b = 0; b < 6; ++b) {
    int l = l0 + dial->nb_dl[b];
    if (l < 0 || dial->size <= (size_t)l) continue;
    if (dial->state[l] == VALID) continue;
    dial3_update_nb(dial, l0, l);
  }
}

bool dial3_step(dial3_s *dial) {
  bucket_s *bucket = dial->first;
  if (bucket == NULL) {
    return false;
  }
  int l0;
  while (bucket->size > 0) {
    l0 = bucket_pop(bucket);
    // NOTE: a node can exist in multiple buckets
    if (dial->state[l0] == VALID) {
      continue;
    }
    dial->state[l0] = VALID;
    dial3_update_nbs(dial, l0);
  }
  return (dial->first = bucket->next) != NULL;
}

void dial3_solve(dial3_s *dial) {
  while (dial3_step(dial)) {}
}

void print_T(dial3_s const *dial) {
  dbl T;
  ivec3 ind, shape = dial->shape;
  for (ind.i = 0; ind.i < shape.i; ++ind.i) {
    for (ind.j = 0; ind.j < shape.j; ++ind.j) {
      for (ind.k = 0; ind.k < shape.k - 1; ++ind.k) {
        T = dial->T[ind2l3(shape, ind)];
        if (isinf(T)) {
          printf("   inf ");
        } else {
          printf("%0.4f ", T);
        }
      }
      T = dial->T[ind2l3(shape, ind)];
      if (isinf(T)) {
        printf("   inf\n");
      } else {
        printf("%0.4f\n", T);
      }
    }
    puts("");
  }
}

void print_bucket(bucket_s const *bucket) {
  printf("{");
  size_t i, j;
  for (i = 0, j = bucket->start; i < bucket->size - 1; ++i) {
    printf("%d, ", bucket->l[j]);
    j = (j + 1) % bucket->capacity;
  }
  printf("%d}\n", bucket->l[j]);
}

void print_buckets(dial3_s const *dial) {
  bucket_s *bucket = dial->first;
  while (bucket) {
    print_bucket(bucket);
    bucket = bucket->next;
  }
}

int main() {
  stype_e stype = CONSTANT;
  int nx = 5;
  int ny = 5;
  int nz = 5;
  dbl h = 2.0/(nx - 1);
  ivec3 shape = {.i = nx, .j = ny, .k = nz};

  dial3_s dial;
  dial3_init(&dial, stype, shape, h);

  ivec3 ind0 = ivec3_int_div(shape, 2);
  int i0 = ind0.i - 1, i1 = ind0.i + 1;
  int j0 = ind0.j - 1, j1 = ind0.j + 1;
  int k0 = ind0.k - 1, k1 = ind0.k + 1;
  ivec3 ind;
  dvec3 x0 = dial3_x(&dial, ind2l3(shape, ind0));
  for (ind.i = i0; ind.i <= i1; ++ind.i) {
    for (ind.j = j0; ind.j <= j1; ++ind.j) {
      for (ind.k = k0; ind.k <= k1; ++ind.k) {
        if (ivec3_equal(ind, ind0)) {
          continue;
        }
        dvec3 x = dial3_x(&dial, ind2l3(shape, ind));
        dvec3 grad_T = dvec3_sub(x, x0);
        dbl T = dvec3_norm(grad_T);
        grad_T = dvec3_dbl_div(grad_T, T);
        dial3_add_trial(&dial, ind, T, grad_T);
      }
    }
  }

  int l0 = ind2l3(shape, ind0);
  dial.T[l0] = 0;
  dial.grad_T[l0] = dvec3_nan();
  dial.state[l0] = VALID;

  print_T(&dial);

  dial3_solve(&dial);

  print_T(&dial);
}
