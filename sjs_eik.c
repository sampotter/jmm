#include "sjs_eik.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "hybrid.h"
#include "index.h"
#include "jet.h"
#include "update.h"

#define NUM_CELL_VERTS 4
#define NUM_NB 8
#define NUM_NEARBY_CELLS 16

struct sjs {
  ivec2 shape;
  dvec2 xymin;
  dbl h;
  int nnodes, ncells;
  int nb_dl[NUM_NB + 1];
  int vert_dl[NUM_CELL_VERTS];
  int tri_dlc[NUM_NB];
  int nb_dlc[NUM_CELL_VERTS];
  int nearby_dlc[NUM_NEARBY_CELLS];
  bicubic_s *bicubics;
  jet_s *jets;
  state_e *states;
  int *positions;
  heap_s *heap;
};

static ivec2 offsets[NUM_NB + 1] = {
  {.i = -1, .j = -1},
  {.i = -1, .j =  0},
  {.i = -1, .j =  1},
  {.i =  0, .j =  1},
  {.i =  1, .j =  1},
  {.i =  1, .j =  0},
  {.i =  1, .j = -1},
  {.i =  0, .j = -1},
  {.i = -1, .j = -1}
};

static ivec2 cell_vert_offsets[NUM_CELL_VERTS] = {
  {.i = 0, .j = 0},
  {.i = 1, .j = 0},
  {.i = 0, .j = 1},
  {.i = 1, .j = 1}
};

/**
 * The array `tri_cell_offsets` contains the `ivec2` offsets that are
 * used to initialize the `tri_dlc` member `sjs_s`, which is in turn
 * used to grad cell-based interpolants used when doing triangle
 * updates. The idea here is that we only want to use cells that
 * belong to the valid part of the expanding solution. If the node `i`
 * below is being updated, then we should only look at cells that lie
 * on the other side of the ring of `i`'s eight nearest neighbors.
 *
 * In the diagram below, we label the cells with the offsets and order
 * them so that they match the ordering used in `tri_cell_offsets`.
 *
 * TODO: maybe try to come up with a slightly more compelling
 * explanation for why we want to do this.
 *
 *                0     1
 *	           +-----+-----+
 *			   |-2,-1|-2, 0|
 *	     +-----+-----+-----+-----+
 *	   7 |-1,-2|  X	 |	X  |-1, 1| 2
 *	     +-----+-----i-----+-----+
 *     6 | 0,-2|  X  |	X  | 0, 1| 3
 *	     +-----+-----+-----+-----+
 *			   | 1,-1| 1, 0|
 *	           +-----+-----+
 *                5     4
 */
static ivec2 tri_cell_offsets[NUM_NB] = {
  {.i = -2, .j = -1},
  {.i = -2, .j =  0},
  {.i = -1, .j =  1},
  {.i =  0, .j =  1},
  {.i =  1, .j =  0},
  {.i =  1, .j = -1},
  {.i =  0, .j = -2},
  {.i = -1, .j = -2}
};

static ivec2 nb_cell_offsets[NUM_CELL_VERTS] = {
  {.i = -1, .j = -1},
  {.i =  0, .j = -1},
  {.i = -1, .j =  0},
  {.i =  0, .j =  0}
};

static ivec2 nearby_cell_offsets[NUM_NEARBY_CELLS] = {
  {.i = -2, .j = -2}, {.i = -2, .j = -1}, {.i = -2, .j = 0}, {.i = -2, .j = 1},
  {.i = -1, .j = -2}, {.i = -1, .j = -1}, {.i = -1, .j = 0}, {.i = -1, .j = 1},
  {.i =  0, .j = -2}, {.i =  0, .j = -1}, {.i =  0, .j = 0}, {.i =  0, .j = 1},
  {.i =  1, .j = -2}, {.i =  1, .j = -1}, {.i =  1, .j = 0}, {.i =  1, .j = 1},
};

static bicubic_variable tri_bicubic_vars[NUM_NB] = {
  /* ic0 -> */ MU, MU, LAMBDA, LAMBDA, MU, MU, LAMBDA, LAMBDA
};

static int tri_edges[NUM_NB] = {
  /* ic0 -> */ 1, 1, 0, 0, 0, 0, 1, 1
};

static bool should_reverse_cubic[NUM_NB] = {
  /* ic0 -> */ true, false, true, false, false, true, false, true
};

#define _ NO_INDEX
static int dirty2updated[NUM_NEARBY_CELLS][NUM_CELL_VERTS - 1] = {
  {0, _, _}, {0, 1, _}, {1, 2, _}, {2, _, _},
  {0, 7, _}, {0, 1, 7}, {1, 2, 3}, {2, 3, _},
  {6, 7, _}, {5, 6, 7}, {3, 4, 5}, {3, 4, _},
  {6, _, _}, {5, 6, _}, {4, 5, _}, {4, _, _}
};
#undef _

static int updated2dirty[NUM_NB][NUM_CELL_VERTS] = {
  {0, 1, 4, 5},
  {1, 2, 5, 6},
  {2, 3, 6, 7},
  {6, 7, 10, 11},
  {10, 11, 14, 15},
  {9, 10, 13, 14},
  {8, 9, 12, 13},
  {4, 5, 8, 9}
};

static void set_nb_dl(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    sjs->nb_dl[i] = ind2l(sjs->shape, offsets[i]);
  }
}

static void set_vert_dl(sjs_s *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->vert_dl[i] = ind2l(sjs->shape, cell_vert_offsets[i]);
  }
}

static void set_tri_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB; ++i) {
    sjs->tri_dlc[i] = ind2lc(sjs->shape, tri_cell_offsets[i]);
  }
}

static void set_nb_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    sjs->nb_dlc[i] = ind2lc(sjs->shape, nb_cell_offsets[i]);
  }
}

static void set_nearby_dlc(sjs_s *sjs) {
  for (int i = 0; i < NUM_NEARBY_CELLS; ++i) {
    sjs->nearby_dlc[i] = ind2lc(sjs->shape, nearby_cell_offsets[i]);
  }
}

static dvec2 get_xy(sjs_s *sjs, int l) {
  ivec2 ind = l2ind(sjs->shape, l);
  dvec2 xy = {
    .x = sjs->h*ind.i + sjs->xymin.x,
    .y = sjs->h*ind.j + sjs->xymin.y
  };
  return xy;
}

static bool line(sjs_s *sjs, int l, int l0) {
  dvec2 xy = get_xy(sjs, l);
  dvec2 xy0 = get_xy(sjs, l0);
  dvec2 xy0_minus_xy = dvec2_sub(xy0, xy);
  dbl L = dvec2_norm(xy0_minus_xy);

  dbl T0 = sjs->jets[l0].f;
  dbl T = T0 + L;

  jet_s *J = &sjs->jets[l];
  bool updated = false;
  if (T < J->f) {
    updated = true;
    J->f = T;
    J->fx = -xy0_minus_xy.x/L;
    J->fy = -xy0_minus_xy.y/L;
  }
  return updated;
}

/**
 * In this function, `ic0` is used as an index to select a nearby
 * bicubic interpolant which will be used to approximate `T` locally.
 */
static bool tri(sjs_s *sjs, int l, int l0, int l1, int ic0) {
  assert(ic0 >= 0);
  assert(ic0 < NUM_NB);

  int lc = l2lc(sjs->shape, l) + sjs->tri_dlc[ic0];
  if (lc < 0 || sjs->ncells <= lc) {
    return false;
  }

  bicubic_s *bicubic = &sjs->bicubics[lc];
#if SJS_DEBUG
  assert(bicubic_valid(bicubic));
#endif

  /**
   * Get cubic along edge of interest.
   */
  bicubic_variable var = tri_bicubic_vars[ic0];
  int edge = tri_edges[ic0];
  cubic_s cubic = bicubic_restrict(bicubic, var, edge);
  if (should_reverse_cubic[ic0]) {
    cubic_reverse_on_unit_interval(&cubic);
  }

  update_data data = {
    .cubic = cubic,
    .xy = get_xy(sjs, l),
    .xy0 = get_xy(sjs, l0),
    .xy1 = get_xy(sjs, l1)
  };

  void *context = (void *)&data;
  dbl lam = hybrid(dF_dt, 0, 1, context);
  dbl T = F(lam, context);

  jet_s *J = &sjs->jets[l];

  bool updated = false;
  if (T < J->f) {
    updated = true;

    dvec2 xy = get_xy(sjs, l);
    dvec2 xylam = dvec2_ccomb(data.xy0, data.xy1, lam);
    dvec2 xylam_minus_xy = dvec2_sub(xylam, xy);
    dbl L = dvec2_norm(xylam_minus_xy);

    J->f = T;
    J->fx = -xylam_minus_xy.x/L;
    J->fy = -xylam_minus_xy.y/L;
  }
  return updated;
}

static bool can_build_cell(sjs_s const *sjs, int lc) {
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(sjs->shape, lc);
    state_e state = sjs->states[l + sjs->vert_dl[i]];
    if (state != TRIAL && state != VALID) {
      return false;
    }
  }
  return true;
}

static void est_Txy_values(sjs_s *sjs, int lc, dbl Txy[4]) {
  dbl fx[NUM_CELL_VERTS], fy[NUM_CELL_VERTS];

  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(sjs->shape, lc) + sjs->vert_dl[i];
    fx[i] = sjs->jets[l].fx;
    fy[i] = sjs->jets[l].fy;
  }

  dbl fxy[NUM_CELL_VERTS] = {
    (fy[1] - fy[0])/sjs->h, // left
    (fx[3] - fx[1])/sjs->h, // bottom
    (fx[2] - fx[0])/sjs->h, // top
    (fy[3] - fy[2])/sjs->h // right
  };

  static dbl lams[4] = {-0.5, 0.5, 0.5, 1.5};
  static dbl mus[4] = {0.5, -0.5, 1.5, 0.5};

  dbl lam, mu;

  for (int i = 0; i < 4; ++i) {
    lam = lams[i];
    mu = mus[i];
    Txy[i] = (1 - mu)*((1 - lam)*fxy[0] + lam*fxy[1]) +
      mu*((1 - lam)*fxy[2] + lam*fxy[3]);
  }
}

/**
 * Build the cell at index `lc`. This uses the jets at the cells
 * vertices to assemble the data matrix for the bicubic
 * interpolant. This function doesn't make any assumptions about the
 * state of the nodes, so it's assumed that the caller has already
 * made sure this is a reasonable thing to try to do.
 */
static void build_cell(sjs_s *sjs, int lc) {
  /* Get linear indices of cell vertices */
  int l[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    l[i] = lc2l(sjs->shape, lc) + sjs->vert_dl[i];
  }

  /* Get jet at each cell vertex */
  jet_s *J[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    J[i] = &sjs->jets[l[i]];
  }

  /* Precompute scaling factors for partial derivatives */
  dbl h = sjs->h, h_sq = h*h;

  /* Compute cell data from jets and scaling factors */
  dmat44 data;
  data.data[0][0] = J[0]->f;
  data.data[1][0] = J[1]->f;
  data.data[0][1] = J[2]->f;
  data.data[1][1] = J[3]->f;
  data.data[2][0] = h*J[0]->fx;
  data.data[3][0] = h*J[1]->fx;
  data.data[2][1] = h*J[2]->fx;
  data.data[3][1] = h*J[3]->fx;
  data.data[0][2] = h*J[0]->fy;
  data.data[1][2] = h*J[1]->fy;
  data.data[0][3] = h*J[2]->fy;
  data.data[1][3] = h*J[3]->fy;
  data.data[2][2] = h_sq*J[0]->fxy;
  data.data[3][2] = h_sq*J[1]->fxy;
  data.data[2][3] = h_sq*J[2]->fxy;
  data.data[3][3] = h_sq*J[3]->fxy;

  /* Set cell data */
  bicubic_set_data(&sjs->bicubics[lc], data);
}

static bool update(sjs_s *sjs, int l) {
  bool updated = false;

  for (int i0 = 1, l0, l1, ic0; i0 < 8; i0 += 2) {
    l0 = l + sjs->nb_dl[i0];
    if (sjs->states[l0] != VALID) {
      continue;
    }

    l1 = l + sjs->nb_dl[i0 - 1];
    if (sjs->states[l1] == VALID) {
      ic0 = i0 - 1;
      updated |= tri(sjs, l, l0, l1, ic0);
    }

    l1 = l + sjs->nb_dl[i0 + 1];
    if (sjs->states[l1] == VALID) {
      ic0 = i0;
      updated |= tri(sjs, l, l0, l1, ic0);
    }
  }

  for (int i0 = 0, l0; i0 < 8; ++i0) {
    l0 = l + sjs->nb_dl[i0];
    if (sjs->states[l0] == VALID) {
      updated |= line(sjs, l, l0);
    }
  }

  return updated;
}

static void adjust(sjs_s *sjs, int l0) {
  assert(sjs->states[l0] == TRIAL);
  assert(l0 >= 0);
  assert(l0 < sjs->nnodes);

  heap_swim(sjs->heap, sjs->positions[l0]);
}

/**
 * Checks if the linear index `l` is valid for `sjs`.
 */
static bool inbounds(sjs_s *sjs, int l) {
  return 0 <= l && l < sjs->nnodes;
}

static dbl value(void *vp, int l) {
  sjs_s *sjs = (sjs_s *)vp;
  assert(l >= 0);
  assert(l < sjs->nnodes);
  dbl T = sjs->jets[l].f;
  return T;
}

static void setpos(void *vp, int l, int pos) {
  sjs_s *sjs = (sjs_s *)vp;
  sjs->positions[l] = pos;
}

void sjs_alloc(sjs_s **sjs) {
  *sjs = malloc(sizeof(sjs_s));
  assert(*sjs != NULL);
}

void sjs_dealloc(sjs_s **sjs) {
  free(*sjs);
  *sjs = NULL;
}

// TODO: since the margins are BOUNDARY nodes, we actually don't need
// to allocate an extra margin of cells, since they will never be
// initialized (i.e., they will never have all of their vertex nodes
// become VALID because of the margin...)
void sjs_init(sjs_s *sjs, ivec2 shape, dvec2 xymin, dbl h) {
  sjs->shape = shape;
  sjs->ncells = (shape.i - 1)*(shape.j - 1);
  sjs->nnodes = shape.i*shape.j;
  sjs->xymin = xymin;
  sjs->h = h;
  sjs->bicubics = malloc(sjs->ncells*sizeof(bicubic_s));
  sjs->jets = malloc(sjs->nnodes*sizeof(jet_s));
  sjs->states = malloc(sjs->nnodes*sizeof(state_e));
  sjs->positions = malloc(sjs->nnodes*sizeof(int));

  assert(sjs->bicubics != NULL);
  assert(sjs->jets != NULL);
  assert(sjs->states != NULL);
  assert(sjs->positions != NULL);

#if SJS_DEBUG
  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->positions[l] = NO_INDEX;
  }
#endif

  heap_alloc(&sjs->heap);

  int capacity = (int) 3*sqrt(sjs->shape.i*sjs->shape.j);
  heap_init(sjs->heap, capacity, value, setpos, (void *)sjs);

  set_nb_dl(sjs);
  set_vert_dl(sjs);
  set_tri_dlc(sjs);
  set_nb_dlc(sjs);
  set_nearby_dlc(sjs);

#if SJS_DEBUG
  for (int lc = 0; lc < sjs->ncells; ++lc) {
    bicubic_invalidate(&sjs->bicubics[lc]);
  }
#endif

  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->jets[l].f = INFINITY;
    sjs->states[l] = FAR;
  }
}

void sjs_deinit(sjs_s *sjs) {
  free(sjs->bicubics);
  free(sjs->jets);
  free(sjs->states);
  free(sjs->positions);

  sjs->bicubics = NULL;
  sjs->jets = NULL;
  sjs->states = NULL;
  sjs->positions = NULL;

  heap_deinit(sjs->heap);
  heap_dealloc(&sjs->heap);
}

void sjs_step(sjs_s *sjs) {
  int l0 = heap_front(sjs->heap);
  assert(sjs->states[l0] == TRIAL);

  heap_pop(sjs->heap);

  // TODO: this is where we could try rebuilding the adjacent cells,
  // but I'm still not convinced this should be necessary if we do
  // things correctly below...

  // Set FAR nodes to TRIAL and insert them into the heap.
  sjs->states[l0] = VALID;
  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_dl[i];
    if (inbounds(sjs, l) && sjs->states[l] == FAR) {
      sjs->states[l] = TRIAL;
      heap_insert(sjs->heap, l);
    }
  }

  // Array keeping track of which nodes get updated.
  bool updated[NUM_NB];
  memset(updated, 0x0, NUM_NB*sizeof(bool));

  // Update neighboring nodes.
  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_dl[i];
    if (inbounds(sjs, l) && sjs->states[l] == TRIAL) {
      updated[i] = update(sjs, l);
      adjust(sjs, l);
    }
  }

  // Find adjacent "dirty" cells ("ready" cells which have at least
  // one updated vertex).
  bool dirty[NUM_NEARBY_CELLS];
  memset(dirty, 0x0, NUM_NEARBY_CELLS*sizeof(bool));
  for (int i = 0, lc; i < NUM_NEARBY_CELLS; ++i) {
    lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[i];
    if (!can_build_cell(sjs, lc)) {
      continue;
    }
    for (int j = 0, k; j < NUM_CELL_VERTS; ++j) {
      k = dirty2updated[i][j];
      dirty[i] |= k != NO_INDEX && updated[k];
    }
  }

  // Estimate Txy values at the vertices of each dirty cell.
  //
  // TODO: this is overkill---we won't actually use all these
  // values... but it may be more work to try and find the values we
  // don't need...
  dbl Txy[NUM_NEARBY_CELLS][4];
  for (int i = 0, lc; i < NUM_NEARBY_CELLS; ++i) {
    if (dirty[i]) {
      lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[i];
      est_Txy_values(sjs, lc, Txy[i]);
    }
  }

  // Set new Txy values at updated nodes to be the average of the Txy
  // values estimated at each dirty cell.
  dbl Txy_sum;
  int ndirty;
  for (int i = 0, l; i < NUM_NB; ++i) {
    Txy_sum = 0;
    ndirty = 0;
    for (int j = 0, k; j < NUM_CELL_VERTS; ++j) {
      k = updated2dirty[i][j];
      if (dirty[k]) {
        Txy_sum += Txy[k][j];
        ++ndirty;
      }
    }
    if (ndirty > 0) {
      l = l0 + sjs->nb_dl[i];
      sjs->jets[l].fxy = Txy_sum/ndirty;
    }
  }

  // Rebuild dirty cells.
  for (int i = 0, lc; i < NUM_NEARBY_CELLS; ++i) {
    if (dirty[i]) {
      lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[i];
      build_cell(sjs, lc);
    }
  }
}

void sjs_solve(sjs_s *sjs) {
  while (heap_size(sjs->heap) > 0) {
    sjs_step(sjs);
  }
}

void sjs_add_trial(sjs_s *sjs, ivec2 ind, jet_s jet) {
  int l = ind2l(sjs->shape, ind);
  sjs->jets[l] = jet;
  sjs->states[l] = TRIAL;
  heap_insert(sjs->heap, l);
}

void sjs_add_valid(sjs_s *sjs, ivec2 ind, jet_s jet) {
  int l = ind2l(sjs->shape, ind);
  sjs->jets[l] = jet;
  sjs->states[l] = VALID;
}

void sjs_make_bd(sjs_s *sjs, ivec2 ind) {
  int l = ind2l(sjs->shape, ind);
  sjs->states[l] = BOUNDARY;
}

ivec2 sjs_get_shape(sjs_s const *sjs) {
  return sjs->shape;
}

jet_s sjs_get_jet(sjs_s *sjs, ivec2 ind) {
  int l = ind2l(sjs->shape, ind);
  return sjs->jets[l];
}

state_e sjs_get_state(sjs_s const *sjs, ivec2 ind) {
  int l = ind2l(sjs->shape, ind);
  return sjs->states[l];
}

state_e *sjs_get_states_ptr(sjs_s const *sjs) {
  return sjs->states;
}

/**
 * The four functions below (`sjs_T`, `sjs_Tx`, `sjs_Ty`, and
 * `sjs_Txy`) are only intended to be used by people consuming this
 * API, not internally.
 *
 * TODO: at some point, we're going to want to add batch versions of
 * these functions to avoid calling `can_build_cell` over and over.
 */

dbl sjs_T(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = xy_to_lc_and_cc(sjs->shape, sjs->xymin, sjs->h, xy, &cc);
  if (!can_build_cell(sjs, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &sjs->bicubics[lc];
  return bicubic_f(bicubic, cc);
}

dbl sjs_Tx(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = xy_to_lc_and_cc(sjs->shape, sjs->xymin, sjs->h, xy, &cc);
  if (!can_build_cell(sjs, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &sjs->bicubics[lc];
  return bicubic_fx(bicubic, cc)/sjs->h;
}

dbl sjs_Ty(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = xy_to_lc_and_cc(sjs->shape, sjs->xymin, sjs->h, xy, &cc);
  if (!can_build_cell(sjs, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &sjs->bicubics[lc];
  return bicubic_fy(bicubic, cc)/sjs->h;
}

dbl sjs_Txy(sjs_s *sjs, dvec2 xy) {
  dvec2 cc;
  int lc = xy_to_lc_and_cc(sjs->shape, sjs->xymin, sjs->h, xy, &cc);
  if (!can_build_cell(sjs, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &sjs->bicubics[lc];
  return bicubic_fxy(bicubic, cc)/(sjs->h*sjs->h);
}

bool sjs_can_build_cell(sjs_s const *sjs, ivec2 indc) {
  int lc = indc2lc(sjs->shape, indc);
  return can_build_cell(sjs, lc);
}

void sjs_build_cells(sjs_s *sjs) {
  for (int lc = 0; lc < sjs->ncells; ++lc) {
    if (can_build_cell(sjs, lc)) {
      build_cell(sjs, lc);
    }
  }
}

bicubic_s sjs_get_bicubic(sjs_s const *sjs, ivec2 indc) {
  int lc = indc2lc(sjs->shape, indc);
  return sjs->bicubics[lc];
}

bicubic_s *sjs_get_bicubics_ptr(sjs_s const *sjs) {
  return sjs->bicubics;
}

heap_s *sjs_get_heap(sjs_s const *sjs) {
  return sjs->heap;
}
