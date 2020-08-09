#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include "bucket.h"
#include "dial.h"
#include "hybrid.h"
#include "index.h"

typedef enum update_status {
  CAUSAL,
  NONCAUSAL
} update_status_e;

typedef update_status_e (*update_f)(dial3_s const *, int /* l */, void * /* ptr */,
                                    dbl * /* T */, dvec3 * /* grad_T */);

struct dial3 {
  stype_e stype;
  ivec3 shape;
  size_t size;
  dbl h;
  dbl gap;

  dbl *T;
  dvec3 *grad_T;

  /**
   * The state of each node. For a Dial-like solver, `state` will
   * never be `TRIAL`. The only permissible states are `VALID`,
   * `FAR`, `BOUNDARY`, and `ADJACENT_TO_BOUNDARY`.
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

dvec3 get_x(ivec3 shape, dbl h, int l) {
#if ORDERING == ROW_MAJOR_ORDERING
  return (dvec3) {
    .data = {
      h*(l/(shape.data[2]*shape.data[1])),
      h*(l/shape.data[2] % shape.data[1]),
      h*(l % shape.data[2])
    }
  };
#else
#  error not implemented yet
#endif
}

typedef struct {
  dvec3 x1;
  dvec3 x;
  dvec3 xsrc;
  dvec3 d;
  dbl T;
  dvec3 grad_T;
} f_update_constant_context_s;

dbl f_update_constant(dbl q, void *context) {
  f_update_constant_context_s *ptr = (f_update_constant_context_s *)context;

  dvec3 xopt = dvec3_saxpy(q, ptr->d, ptr->x1);
  dvec3 D = dvec3_sub(ptr->x, xopt);
  dvec3 Dsrc = dvec3_sub(ptr->xsrc, xopt);

  dbl r = dvec3_norm(D);
  dbl rsrc = dvec3_norm(Dsrc);

  ptr->T = r + rsrc;
  ptr->grad_T = dvec3_dbl_div(D, r);

  return dvec3_dot(ptr->d, dvec3_add(dvec3_dbl_div(D, r), ptr->grad_T));
}

typedef struct update_constant_data {
  dvec3 x0;
  dvec3 xsrc;
  dvec3 x0_minus_xsrc;
} update_constant_data_s;

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
update_status_e
update_constant(dial3_s const *dial, int l, void *ptr, dbl *T, dvec3 *grad_T) {
  update_constant_data_s *data = (update_constant_data_s *)ptr;

  // do the projection
  dvec3 x = get_x(dial->shape, dial->h, l);
  dvec3 dx = dvec3_sub(x, data->x0);
  dvec3 t = dvec3_sub(x, data->xsrc);

  // Compute xs. While we're at it, checking xsrc is ahead of any
  // possible wavefront passing through x0. If it is, return INFINITY here.
  // This involves checking the angle between x - x0 and xsrc - x0.
  // Note that if we pass this check, then the denominator up ahead
  // involving the dot product between x - x0 and x - xsrc can't be zero.
  dbl s = dvec3_dot(dx, data->x0_minus_xsrc);
  if (s <= 0) {
    return NONCAUSAL;
  }
  s /= dvec3_dot(dx, t);
  dvec3 xs = dvec3_saxpy(s, t, data->xsrc);

  // TODO: Right now we'll just try something really dumb---don't
  // check *where* xs lands, just check if it lands inside the box. If
  // it does, compute a new value and return it.

  static dvec3 X1[3][4] = {
    {{.data = { 0,  1,  0}},
     {.data = { 0,  0,  1}},
     {.data = { 0, -1,  0}},
     {.data = { 0,  0, -1}}},
    {{.data = { 1,  0,  0}},
     {.data = { 0,  0,  1}},
     {.data = {-1,  0,  0}},
     {.data = { 0,  0, -1}}},
    {{.data = { 1,  0,  0}},
     {.data = { 0,  1,  0}},
     {.data = {-1,  0,  0}},
     {.data = { 0, -1,  0}}}
  };

  dbl h = dial->h;

  dvec3 xs_minus_x0 = dvec3_sub(xs, data->x0);
  dbl xs_minus_x0_maxnorm = dvec3_maxnorm(xs_minus_x0);
  if (xs_minus_x0_maxnorm > h + EPS) {
    if (dial->state[l] != ADJACENT_TO_BOUNDARY) {
      return NONCAUSAL;
    }

    // Find x1
    dvec3 x1;

    // Find the coordinate spanned by x - x0
    int i1;
    {
      dbl dx1 = x.data[1] - data->x0.data[1];
      dbl dx2 = x.data[2] - data->x0.data[2];
      // TODO: probably *do* need to use fabs for this one?
      i1 = (fabs(dx1) > EPS) + 2*(fabs(dx2) > EPS);
    }

    // Find x1 itself
    dbl dist_sq = INFINITY;
    for (int j = 0; j < 4; ++j) {
      dvec3 new_x1 = dvec3_add(data->x0, dvec3_dbl_mul(X1[i1][j], dial->h));
      dbl new_dist_sq = dvec3_dist_sq(new_x1, data->xsrc);
      if (new_dist_sq < dist_sq) {
        x1 = new_x1;
        dist_sq = new_dist_sq;
      }
    }

    int i2;
    { // TODO: should be able to just recover index when finding x1
      dbl dx1 = x1.data[1] - data->x0.data[1];
      dbl dx2 = x1.data[2] - data->x0.data[2];
      // TODO: probably don't need to use fabs for this one
      i2 = (fabs(dx1) > EPS) + 2*(fabs(dx2) > EPS);
    }

    // Find d
    dvec3 d = dvec3_one();
    d.data[i1] = 0;
    d.data[i2] = 0;

    f_update_constant_context_s context = {
      .x1 = x1,
      .x = x,
      .xsrc = data->xsrc,
      .d = d
    };
    hybrid(f_update_constant, -1, 1, (void *)&context);

    *T = context.T;
    *grad_T = context.grad_T;
  } else {
    *T = dvec3_norm(t);
    *grad_T = dvec3_dbl_div(t, *T);
  }

  return CAUSAL;
}

update_f update_functions[NUM_STYPE] = {
  update_constant // CONSTANT
};

void dial3_alloc(dial3_s **dial) {
  *dial = malloc(sizeof(dial3_s));
}

error_e dial3_init(dial3_s *dial, stype_e stype, int const *shape, dbl h) {
  if (stype != CONSTANT) {
    return BAD_ARGUMENT;
  }

  dial->stype = stype;
  dial->shape.data[0] = shape[0];
  dial->shape.data[1] = shape[1];
  dial->shape.data[2] = shape[2];
  dial->size = ivec3_prod(dial->shape);
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
  dial->nb_dl[0] = ind2l3(dial->shape, (ivec3) {.data = {-1,  0,  0}});
  dial->nb_dl[1] = ind2l3(dial->shape, (ivec3) {.data = { 0, -1,  0}});
  dial->nb_dl[2] = ind2l3(dial->shape, (ivec3) {.data = { 0,  0, -1}});
  dial->nb_dl[3] = ind2l3(dial->shape, (ivec3) {.data = { 0,  0,  1}});
  dial->nb_dl[4] = ind2l3(dial->shape, (ivec3) {.data = { 0,  1,  0}});
  dial->nb_dl[5] = ind2l3(dial->shape, (ivec3) {.data = { 1,  0,  0}});
#else
#  error not implemented yet
#endif

  dial->first = NULL;

  dial->update = update_functions[dial->stype];

  return SUCCESS;
}

void dial3_deinit(dial3_s *dial) {
  free(dial->T);
  free(dial->grad_T);
  free(dial->state);
  free(dial->lb);
}

void dial3_dealloc(dial3_s **dial) {
  free(*dial);
  *dial = NULL;
}

int get_lb(dial3_s const *dial, dbl T) {
  return T/dial->gap;
}

void prepend_buckets(dial3_s *dial, int lb) {
  while (lb < dial->lb0) {
    bucket_s *bucket;
    bucket_alloc(&bucket);
    bucket_init(bucket);
    bucket_set_next(bucket, dial->first);
    dial->first = bucket;
    --dial->lb0;
  }
}

bucket_s *find_bucket(dial3_s *dial, int lb) {
  if (lb < dial->lb0) {
    prepend_buckets(dial, lb);
  }
  bucket_s *bucket = dial->first;
  while (lb > dial->lb0) {
    bucket_s *next = bucket_get_next(bucket);
    if (next == NULL) {
      bucket_alloc(&next);
      bucket_init(next);
      bucket_set_next(bucket, next);
    }
    bucket = next;
    --lb;
  }
  assert(bucket != NULL);
  return bucket;
}

void dial3_insert(dial3_s *dial, int l, dbl T) {
  assert(dial->state[l] != BOUNDARY);
  if (dial->first == NULL) {
    dial->lb0 = get_lb(dial, T);
    bucket_alloc(&dial->first);
    bucket_init(dial->first);
    bucket_push(dial->first, l);
  } else {
    int lb = get_lb(dial, T);
    bucket_s *bucket = find_bucket(dial, lb);
    bucket_push(bucket, l);
  }
}

void dial3_set_T(dial3_s *dial, int l, dbl T) {
  // TODO: for now, we just assume that T >= 0. In the future, we can
  // be more sophisticated about this, and introduce a `Tmin`
  // parameter. But then we'll need to adjust `bucket_T` to
  // ensure that the [Tmin, Tmin + gap) bucket corresponds corresponds
  // to `lb == 0` (so that `NO_INDEX == -1` still works).
  assert(T >= 0);
  dial->T[l] = T;
}

void update_nb(dial3_s *dial, int l, void *ptr) {
  dbl T;
  dvec3 grad_T;
  update_status_e status = dial->update(dial, l, ptr, &T, &grad_T);

  // TODO: may want to add a little tolerance here to ensure we don't
  // mess with grad_T too much?
  if (status == CAUSAL && T < dial->T[l]) {
    dial3_set_T(dial, l, T);
    dial->grad_T[l] = grad_T;
    int lb = get_lb(dial, T);
    if (lb != dial->lb[l]) {
      assert(dial->lb[l] == NO_INDEX || (dial->lb0 <= lb && lb < dial->lb[l]));
      bucket_s *bucket = find_bucket(dial, lb);
      bucket_push(bucket, l);
      // TODO: should I actually do this?
      dial->lb[l] = lb;
    }
  }
}

bool dial3_can_update(dial3_s const *dial, int l) {
  return dial->state[l] == FAR || dial->state[l] == ADJACENT_TO_BOUNDARY;
}

void update_nbs(dial3_s *dial, int l0) {
  /**
   * Compute the local characteristic direction (tangent vector of the
   * ray) at the node l0. It is always the case that the gradient of T
   * points in the same direction. Additionally, for s = 1, the
   * gradient of T has unit magnitude and is already normalized
   * correctly.
   */
  dvec3 t0 = dial->grad_T[l0]; // Already normalized for s = 1!

  update_constant_data_s data;
  data.x0 = get_x(dial->shape, dial->h, l0);
  data.xsrc = dvec3_saxpy(-dial->T[l0], t0, data.x0);
  data.x0_minus_xsrc = dvec3_sub(data.x0, data.xsrc);

  for (int b = 0; b < 6; ++b) {
    int l = l0 + dial->nb_dl[b];
    if (l < 0 || dial->size <= (size_t)l) continue;
    if (!dial3_can_update(dial, l)) continue;
    update_nb(dial, l, (void *)&data);
  }
}

/**
 * TODO: check for nodes that are already VALID in the neighborhood of
 * the point source
 */
void dial3_add_point_source(dial3_s *dial, int const *ind0_data, dbl T0) {
  ivec3 shape = dial->shape;

  ivec3 ind0 = {.data = {ind0_data[0], ind0_data[1], ind0_data[2]}};
  int l0 = ind2l3(shape, ind0);
  if (dial->state[l0] != FAR && dial->state[l0] != ADJACENT_TO_BOUNDARY) {
    fprintf(stderr, "ERROR: tried to create a point source at a node "
            "that didn't have FAR or ADJACENT_TO_BOUNDARY state\n");
    exit(EXIT_FAILURE);
  }
  dvec3 x0 = ivec3_dbl_mul(ind0, dial->h);

  for (int i = 0; i < 26; ++i) {
    ivec3 ind = ivec3_add(ind0, NB_IND_OFFSET_26[i]);
    if (!inbounds(shape, ind)) {
      continue;
    }

    // TODO: want to think a little more carefully about what states
    // to skip on here...
    int l = ind2l3(dial->shape, ind);
    if (dial->state[l] == BOUNDARY) {
      continue;
    }

    dvec3 x = ivec3_dbl_mul(ind, dial->h);
    dvec3 grad_T = dvec3_sub(x, x0);
    dbl T = dvec3_norm(grad_T);
    grad_T = dvec3_dbl_div(grad_T, T);

    dial->T[l] = T;
    dial->grad_T[l] = grad_T;

    dial3_insert(dial, l, T);
  }

  dial->T[l0] = T0;
  dial->grad_T[l0] = dvec3_nan();
  dial->state[l0] = VALID;
}

void dial3_add_boundary_points(dial3_s *dial, int const *inds, size_t n) {
  ivec3 shape = dial->shape;

  int *l = malloc(sizeof(int)*n);
  for (size_t i = 0; i < n; ++i) {
    ivec3 ind = {.data = {inds[3*i], inds[3*i + 1], inds[3*i + 2]}};
    printf("%d %d %d\n", ind.data[0], ind.data[1], ind.data[2]);
    l[i] = ind2l3(shape, ind);
    for (size_t j = 0; j < 6; ++j) {
      ivec3 ind_ = ivec3_add(ind, NB_IND_OFFSET_6[j]);
      if (inbounds(shape, ind_)) {
        int l_ = ind2l3(shape, ind_);
        if (dial->state[l_] != BOUNDARY) {
          dial->state[l_] = ADJACENT_TO_BOUNDARY;
        }
      }
    }
  }
  for (size_t i = 0; i < n; ++i) {
    int l_ = l[i];
    dial->state[l_] = BOUNDARY;
    dial->T[l_] = NAN;
    dial->grad_T[l_] = dvec3_nan();
  }
  free(l);
}

bool dial3_step(dial3_s *dial) {
  bucket_s *bucket = dial->first;
  if (bucket == NULL) {
    return false;
  }
  int l0;
  while (bucket_get_size(bucket) > 0) {
    l0 = bucket_pop(bucket);
    assert(dial->state[l0] != BOUNDARY);
    // NOTE: a node can exist in multiple buckets
    if (dial->state[l0] == VALID) {
      continue;
    }
    dial->state[l0] = VALID;
    update_nbs(dial, l0);
  }
  assert(dial->first == bucket);
  do {
    dial->first = bucket_get_next(bucket);
    bucket_deinit(bucket);
    bucket_dealloc(&bucket);
    bucket = dial->first;
    ++dial->lb0;
  } while (bucket != NULL && bucket_is_empty(bucket));
  return bucket != NULL;
}

void dial3_solve(dial3_s *dial) {
  while (dial3_step(dial)) {}
}

dbl dial3_get_T(dial3_s const *dial, int l) {
  return dial->T[l];
}

dbl *dial3_get_T_ptr(dial3_s const *dial) {
  return dial->T;
}

void dial3_get_grad_T(dial3_s const *dial, int l, dbl *grad_T) {
  grad_T[0] = dial->grad_T[l].data[0];
  grad_T[1] = dial->grad_T[l].data[1];
  grad_T[2] = dial->grad_T[l].data[2];
}

dbl *dial3_get_grad_T_ptr(dial3_s const *dial) {
  return &dial->grad_T[0].data[0];
}

state_e *dial3_get_state_ptr(dial3_s const *dial) {
  return dial->state;
}
