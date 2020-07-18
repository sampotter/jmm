#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "def.h"
#include "vec.h"

int ind2l3(ivec3 shape, ivec3 ind) {
#if ORDERING == ROW_MAJOR_ORDERING
  return ind.k + shape.k*(ind.j + shape.j*ind.i);
#else
#  error not implemented yet
#endif
}

ivec3 l2ind3(ivec3 shape, int l) {
#if ORDERING == ROW_MAJOR_ORDERING
  return (ivec3) {
    .i = l/(shape.k*shape.j),
    .j = l/shape.k % shape.j,
    .k = l % shape.k
  };
#else
#  error not implemented yet
#endif
}

typedef struct bucket bucket_s;

/**
 * A bucket is a one-directional queue, implemented as a ring buffer,
 * storing (l)inear indices of nodes in a Dial-like solver. It's also
 * a node in a linked list of buckets.
 */
struct bucket {
  int *l;
  size_t size;
  bucket_s *next;
};

typedef struct dial3 dial3_s;

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

dbl dial3_update_constant(dial3_s const *dial, int l0, int l) {
  // use gradient to compute spherical approximation to wavefront
  // passing through l0

  dvec3 x0 = dial3_x(dial, l0);
  dvec3 t0 = dvec3_normalized(dial->grad_T[l0]); // TODO: for s = 1, this will already be normalized! nice!
  dbl T0 = dial->T[l0];

  dvec3 xsrc = dvec3_saxpy(-T0, t0, x0);

  // project l onto the box of radius h surrounding x0

  // xsrc*dx/((xsrc - dx)*dx)

  dvec3 x = dial3_x(dial, l);
  dvec3 dx = dvec3_sub(x, x0);
  dvec3 t = dvec3_sub(x, xsrc);
  dbl s = dvec3_dot(dx, dvec3_sub(x0, xsrc))/dvec3_dot(dx, t);
  dvec3 xs = dvec3_saxpy(s, t, xsrc);

  // TODO: Right now we'll just try something really dumb---don't
  // check *where* xs lands, just check if it lands inside the box. If
  // it does, compute a new value and return it.

  dbl h = dial->h;
  if (dvec3_maxdist(xs, x0) > h) {
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
  dial->nb_dl[0] = ind2l3((ivec3) {.i = -1,  0,  0}, dial->shape);
  dial->nb_dl[1] = ind2l3((ivec3) {.i =  0, -1,  0}, dial->shape);
  dial->nb_dl[2] = ind2l3((ivec3) {.i =  0,  0, -1}, dial->shape);
  dial->nb_dl[3] = ind2l3((ivec3) {.i =  0,  0,  1}, dial->shape);
  dial->nb_dl[4] = ind2l3((ivec3) {.i =  0,  1,  0}, dial->shape);
  dial->nb_dl[5] = ind2l3((ivec3) {.i =  1,  0,  0}, dial->shape);
#else
#  error not implemented yet
#endif

  dial->first = NULL;

  dial->update = update_functions[dial->stype];
}

/**
 * TODO: this is very simple-minded for now... Later, we'll need to
 * account for different buckets, but for a point source this should
 * be fine
 */
void dial3_add_trial(dial3_s *dial, ivec3 ind, dbl T) {
  int l = ind2l3(dial->shape, ind);
  dial->T[l] = T;
  // TODO: add to bucket
}

int dial3_bucket_T(dial3_s const *dial, dbl T) {
  return T/dial->gap;
}

void dial3_update_nb(dial3_s *dial, int l0, int l) {
  dbl T = dial->update(dial, l0, l);
  if (T < dial->T[l]) {
    dial->T[l] = T;
    int lb = dial3_bucket_T(dial, T);
    if (lb != dial->lb[l]) {
      // TODO: update bucket index and move node...
    }
  }
}

void dial3_update_nbs(dial3_s *dial, int l0) {
  for (int b = 0; b < 6; ++b) {
    int l = l0 + dial->nb_dl[b];
    if (l < 0 || dial->size <= (size_t)l) continue;
    dial3_update_nb(dial, l0, l);
  }
}

void dial3_update_bucket_nodes(dial3_s *dial, int l0) {
  dial->state[l0] = VALID;
  // TODO: remove from bucket
  dial3_update_nbs(dial, l0);
}

// TODO: quick observations here:
// - after we set state(l0) to VALID, we want to evict it
// - when we update nodes, we may move them to a new bucket if their value changes enough
// - doing this in parallel will be tricky

bool dial3_step(dial3_s *dial) {
  bucket_s *bucket = dial->first;
  if (bucket == NULL) {
    return false;
  }
  for (size_t i = 0; i < bucket->size; ++i) {
    dial3_update_bucket_nodes(dial, bucket->l[i]);
  }
  // TODO: or "while dial->first is empty { ... }"
  if (bucket->size == 0) {
    dial->first = bucket->next;
  }
  return dial->first != NULL;
}

void dial3_solve(dial3_s *dial) {
  while (dial3_step(dial)) {}
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
  dial3_add_trial(&dial, ivec3_int_div(shape, 2), 0);
  dial3_solve(&dial);

  dbl T;
  ivec3 ind;
  for (ind.i = 0; ind.i < nx; ++ind.i) {
    for (ind.j = 0; ind.j < ny; ++ind.j) {
      for (ind.k = 0; ind.k < nz - 1; ++ind.k) {
        T = dial.T[ind2l3(shape, ind)];
        if (isinf(T)) {
          printf("   inf ");
        } else {
          printf("%0.4f ", T);
        }
      }
      T = dial.T[ind2l3(shape, ind)];
      if (isinf(T)) {
        printf("   inf\n");
      } else {
        printf("%0.4f\n", T);
      }
    }
    puts("");
  }
}
