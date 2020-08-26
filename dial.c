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

ivec3 NB_IND_OFFSET_6[6] = {
  {.data = {-1,  0,  0}},
  {.data = { 0, -1,  0}},
  {.data = { 0,  0, -1}},
  {.data = { 0,  0,  1}},
  {.data = { 0,  1,  0}},
  {.data = { 1,  0,  0}}
};

ivec3 NB_IND_OFFSET_26[26] = {
  {.data = {-1, -1, -1}},
  {.data = {-1, -1,  0}},
  {.data = {-1, -1,  1}},
  {.data = { 0, -1, -1}},
  {.data = { 0, -1,  0}},
  {.data = { 0, -1,  1}},
  {.data = { 1, -1, -1}},
  {.data = { 1, -1,  0}},
  {.data = { 1, -1,  1}},
  {.data = {-1,  0, -1}},
  {.data = {-1,  0,  0}},
  {.data = {-1,  0,  1}},
  {.data = { 0,  0, -1}},
  {.data = { 0,  0,  1}},
  {.data = { 1,  0, -1}},
  {.data = { 1,  0,  0}},
  {.data = { 1,  0,  1}},
  {.data = {-1,  1, -1}},
  {.data = {-1,  1,  0}},
  {.data = {-1,  1,  1}},
  {.data = { 0,  1, -1}},
  {.data = { 0,  1,  0}},
  {.data = { 0,  1,  1}},
  {.data = { 1,  1, -1}},
  {.data = { 1,  1,  0}},
  {.data = { 1,  1,  1}}
};

struct dial3 {
  stype_e stype;
  ivec3 shape;
  size_t size;
  dbl h;
  dbl gap;

  dbl *Toff;
  dvec3 *xsrc;

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
  dvec3 xsrc0;
  dvec3 d;
  dvec3 xsrc;
  dbl T0;
} f_update_constant_context_s;

dbl f_update_constant(dbl q, void *context) {
  f_update_constant_context_s *ptr = (f_update_constant_context_s *)context;
  ptr->xsrc = dvec3_saxpy(q, ptr->d, ptr->x1);
  dvec3 t = dvec3_sub(ptr->x, ptr->xsrc);
  ptr->T0 = dvec3_norm(t);
  t = dvec3_dbl_div(t, ptr->T0);
  dvec3 t0 = dvec3_normalized(dvec3_sub(ptr->xsrc, ptr->xsrc0));
  return dvec3_dot(ptr->d, dvec3_sub(t, t0));
}

//dbl f_update_constant(dbl q, void *context) {
//  f_update_constant_context_s *ptr = (f_update_constant_context_s *)context;
//
//  dvec3 xopt = dvec3_saxpy(q, ptr->d, ptr->x1);
//  dvec3 D = dvec3_sub(ptr->x, xopt);
//  dvec3 Dsrc = dvec3_sub(ptr->xsrc0, xopt);
//
//  dbl r = dvec3_norm(D);
//  dbl rsrc = dvec3_norm(Dsrc);
//
//  ptr->T = r + rsrc;
//  ptr->grad_T = dvec3_dbl_div(D, r);
//
//  return dvec3_dot(ptr->d, dvec3_add(dvec3_dbl_div(Dsrc, rsrc), ptr->grad_T));
//}

typedef struct update_constant_data {
  dvec3 x0;
  dvec3 xsrc0;
  dvec3 x0_minus_xsrc0;
  dbl Toff0;
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
update_constant(dial3_s const *dial, int l, void *ptr, dbl *Toff, dvec3 *xsrc) {
  update_constant_data_s *data = (update_constant_data_s *)ptr;

  // do the projection
  dvec3 x = get_x(dial->shape, dial->h, l);
  dvec3 dx = dvec3_sub(x, data->x0);
  dvec3 t = dvec3_sub(x, data->xsrc0);

  dbl h = dial->h;

  // Compute xs. While we're at it, checking xsrc is ahead of any
  // possible wavefront passing through x0. If it is, return INFINITY here.
  // This involves checking the angle between x - x0 and xsrc - x0.
  // Note that if we pass this check, then the denominator up ahead
  // involving the dot product between x - x0 and x - xsrc can't be zero.
  dvec3 xs;
  {
    dbl s = dvec3_dot(dx, data->x0_minus_xsrc0);
    if (s < -EPS) {
      return NONCAUSAL;
    } else if (fabs(s) <= EPS) {
      // If the dot product between dx and x0 - xsrc0 is zero,
      // then the projection is no longer well-defined.
      // Instead, we shrink x0 - xsrc onto the hit box.
      s = -h/dvec3_maxnorm(data->x0_minus_xsrc0);
      xs = dvec3_saxpy(s, data->x0_minus_xsrc0, data->x0);

      // We can compute Toff and xsrc immediately in this special case
      //
      // TODO: may not actually need to compute dist below... should be
      // able to just retrieve the distance from s when projecting onto
      // the hit box
      *Toff = data->Toff0 + dvec3_dist(xs, data->xsrc0);
      *xsrc = xs;
      return CAUSAL;
    } else {
      s /= dvec3_dot(dx, t);
      xs = dvec3_saxpy(s, t, data->xsrc0);
    }
  }

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

  dvec3 xs_minus_x0 = dvec3_sub(xs, data->x0);
  dbl xs_minus_x0_maxnorm = dvec3_maxnorm(xs_minus_x0);
  bool in_hit_box = xs_minus_x0_maxnorm < h + EPS;

  if (dial->state[l] != ADJACENT_TO_BOUNDARY) { // TODO: explain this shortcut!
    if (in_hit_box) {
      *Toff = data->Toff0;
      *xsrc = data->xsrc0;
      return CAUSAL;
    } else {
      return NONCAUSAL;
    }
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
    dbl new_dist_sq = dvec3_dist_sq(new_x1, data->xsrc0);
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
  //
  // TODO: explain what d is here
  dvec3 d = dvec3_one();
  d.data[i1] = 0;
  d.data[i2] = 0;

  // Get x1's state. If x1 is a BOUNDARY node, then we replace x1 with x0.
  // We do this *after* computing d, so that it is still the tangent vector of the original edge of the hit box.
  // TODO: document how (& why...?) this works.
  ivec3 ind1 = dvec3_to_ivec3(dvec3_dbl_div(x1, h));
  int l1 = ind2l3(dial->shape, ind1);
  if (dial->state[l1] == BOUNDARY) {
    x1 = data->x0;
  } else if (in_hit_box) {
    *Toff = data->Toff0;
    *xsrc = data->xsrc0;
    return CAUSAL;
  }

  f_update_constant_context_s context = {
    .x1 = x1,
    .x = x,
    .xsrc0 = data->xsrc0,
    .d = d
  };
  hybrid(f_update_constant, -1, 1, (void *)&context);

  *Toff = data->Toff0 + dvec3_dist(context.xsrc, data->xsrc0);
  *xsrc = context.xsrc;

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

  dial->Toff = malloc(sizeof(dbl)*dial->size);
  for (size_t i = 0; i < dial->size; ++i) {
    dial->Toff[i] = INFINITY;
  }

  // xsrc is initialized with garbage
  dial->xsrc = malloc(sizeof(dvec3)*dial->size);

  dial->state = malloc(sizeof(state_e)*dial->size);
  for (size_t i = 0; i < dial->size; ++i) {
    dial->state[i] = FAR;
  }

  dial->lb = malloc(sizeof(int)*dial->size);
  for (size_t i = 0; i < dial->size; ++i) {
    dial->lb[i] = NO_INDEX;
  }

  dial->lb0 = NO_INDEX;

  dial->first = NULL;

  dial->update = update_functions[dial->stype];

  return SUCCESS;
}

void dial3_deinit(dial3_s *dial) {
  free(dial->Toff);
  free(dial->xsrc);
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

void update_nb(dial3_s *dial, int l, void *ptr) {
  dbl Toff;
  dvec3 xsrc;
  update_status_e status = dial->update(dial, l, ptr, &Toff, &xsrc);

  dvec3 x = get_x(dial->shape, dial->h, l);
  dbl T = Toff + dvec3_dist(x, xsrc);

  // TODO: may want to add a little tolerance here to ensure we don't
  // mess with grad_T too much?
  if (status == CAUSAL && T < dial3_get_T(dial, l)) {
    dial->Toff[l] = Toff;
    dial->xsrc[l] = xsrc;
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

bool inbounds(ivec3 shape, ivec3 ind) {
  return 0 <= ind.data[0] && ind.data[0] < shape.data[0]
    && 0 <= ind.data[1] && ind.data[1] < shape.data[1]
    && 0 <= ind.data[2] && ind.data[2] < shape.data[2];
}

void update_nbs(dial3_s *dial, int l0) {
  ivec3 shape = dial->shape;

  // TODO: simplify get_x below now that we're computing this
  ivec3 ind0 = l2ind3(shape, l0);

  update_constant_data_s data;
  // data.x0 = get_x(shape, dial->h, l0);
  data.x0 = ivec3_dbl_mul(ind0, dial->h);
  data.xsrc0 = dial->xsrc[l0];
  data.x0_minus_xsrc0 = dvec3_sub(data.x0, data.xsrc0);
  data.Toff0 = dial->Toff[l0];

  for (int b = 0; b < 6; ++b) {
    ivec3 ind = ivec3_add(ind0, NB_IND_OFFSET_6[b]);
    int l = ind2l3(shape, ind);
    if (!inbounds(shape, ind) ||
        !dial3_can_update(dial, l)) {
      continue;
    }
    update_nb(dial, l, (void *)&data);
  }
}

/**
 * TODO: check for nodes that are already VALID in the neighborhood of
 * the point source
 */
void dial3_add_point_source(dial3_s *dial, int const *ind0_data, dbl Toff) {
  ivec3 shape = dial->shape;

  ivec3 ind0 = {.data = {ind0_data[0], ind0_data[1], ind0_data[2]}};
  int l0 = ind2l3(shape, ind0);
  if (dial->state[l0] != FAR && dial->state[l0] != ADJACENT_TO_BOUNDARY) {
    fprintf(stderr, "ERROR: tried to create a point source at a node "
            "that didn't have FAR or ADJACENT_TO_BOUNDARY state\n");
    exit(EXIT_FAILURE);
  }
  dvec3 xsrc = ivec3_dbl_mul(ind0, dial->h);

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

    dial->Toff[l] = Toff;
    dial->xsrc[l] = xsrc;

    // TODO: change dial3_insert interface so that it just takes dial and l?
    dvec3 x = ivec3_dbl_mul(ind, dial->h);
    dbl T = dvec3_dist(x, xsrc);
    dial3_insert(dial, l, T);
  }

  dial->Toff[l0] = Toff;
  dial->xsrc[l0] = xsrc;
  dial->state[l0] = VALID;
}

void dial3_add_boundary_points(dial3_s *dial, int const *inds, size_t n) {
  ivec3 shape = dial->shape;

  int *l = malloc(sizeof(int)*n);
  for (size_t i = 0; i < n; ++i) {
    ivec3 ind = {.data = {inds[3*i], inds[3*i + 1], inds[3*i + 2]}};
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
    dial->Toff[l_] = NAN;
    dial->xsrc[l_] = dvec3_nan();
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
  dvec3 x = get_x(dial->shape, dial->h, l);
  dbl dT = dvec3_dist(x, dial->xsrc[l]);
  return dial->Toff[l] + dT;
}

void dial3_get_grad_T(dial3_s const *dial, int l, dbl *grad_T) {
  dvec3 x = get_x(dial->shape, dial->h, l);
  dvec3 tmp = dvec3_sub(x, dial->xsrc[l]);
  tmp = dvec3_normalized(tmp);
  grad_T[0] = tmp.data[0];
  grad_T[1] = tmp.data[1];
  grad_T[2] = tmp.data[2];
}

dbl *dial3_get_Toff_ptr(dial3_s const *dial) {
  return dial->Toff;
}

dbl *dial3_get_xsrc_ptr(dial3_s const *dial) {
  return &dial->xsrc[0].data[0];
}

state_e *dial3_get_state_ptr(dial3_s const *dial) {
  return dial->state;
}
