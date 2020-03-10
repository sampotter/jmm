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
#define NUM_CELL_NB_VERTS 9
#define NUM_NB 8
#define NUM_NB_CELLS 4
#define NUM_NEARBY_CELLS 16

/**
 * TODO: add a few words about what `sjs` is and how it works
 *
 * TODO: add a few comments about how indexing works
 *
 * - row-major ordering is used
 * - l vs lc index spaces
 * - cell verts are in column major order
 */
struct sjs {
  ivec2 shape;
  dvec2 xymin;
  dbl h;
  int nnodes, ncells;
  int nb_dl[NUM_NB + 1];
  int cell_nb_verts_dl[NUM_CELL_NB_VERTS];
  int vert_dl[NUM_CELL_VERTS];
  int tri_dlc[NUM_NB];
  int nb_dlc[NUM_NB_CELLS];
  int nearby_dlc[NUM_NEARBY_CELLS];
  bicubic_s *bicubics;
  jet_s *jets;
  state_e *states;
  int *positions;
  heap_s *heap;
};

/**
 * This array maps a linear indexing of neighbors surrounding a grid
 * node to local (i, j) offsets. There are NUM_NB + 1 = 9 entries
 * instead of 8 so that it's not necessary to use modular arithmetic
 * to access the offsets at i +/- 1 in a for loop (this is nice for
 * triangle updates).
 */
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

/**
 * This array provides offsets that can be used to map an `lc`-index
 * to the `l`-indices of the four vertices of the corresponding
 * cell. In our indexing scheme, if we treat an `lc`-index as an `l`
 * index, we recover the upper-left corner of a cell. These offsets
 * can be used to map to any of the cell's vertices.
 *
 * TODO: this should be stored in row-major order instead of column
 * major order!
 *
 * TODO: add a picture
 */
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
 *             +-----+-----+
 *             |-2,-1|-2, 0|
 *       +-----+-----+-----+-----+
 *     7 |-1,-2|  X  |  X  |-1, 1| 2
 *       +-----+-----i-----+-----+
 *     6 | 0,-2|  X  |  X  | 0, 1| 3
 *       +-----+-----+-----+-----+
 *             | 1,-1| 1, 0|
 *             +-----+-----+
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

/**
 * TODO: document this
 *
 * - used to set nearby_dlc
 */
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
static int cell_verts_to_cell_nb_verts[NUM_NEARBY_CELLS][NUM_CELL_VERTS] = {
  {_, _, _, 0}, {_, 0, _, 1}, {_, 1, _, 2}, {_, 2, _, _},
  {_, _, 0, 3}, {0, 3, 1, 4}, {1, 4, 2, 5}, {2, 5, _, _},
  {_, _, 3, 6}, {3, 6, 4, 7}, {4, 7, 5, 8}, {5, 8, _, _},
  {_, _, 6, _}, {6, _, 7, _}, {7, _, 8, _}, {8, _, _, _}
};
#undef _

#define _ NO_INDEX
static int cell_nb_verts_to_nearby_cells[NUM_CELL_NB_VERTS][NUM_NB_CELLS] = {
  { 5,  1,  4,  0},
  { 6,  2,  5,  1},
  { 7,  3,  6,  2},
  { 9,  5,  8,  4},
  {10,  6,  9,  5},
  {11,  7, 10,  6},
  {13,  9, 12,  8},
  {14, 10, 13,  9},
  {15, 11, 14, 10}
};
#undef _

static bool nearby_cell_is_cell_nb[NUM_NEARBY_CELLS] = {
  false, false, false, false,
  false,  true,  true, false,
  false,  true,  true, false,
  false, false, false, false
};

static void set_nb_dl(sjs_s *sjs) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    sjs->nb_dl[i] = ind2l(sjs->shape, offsets[i]);
  }
}

static void set_cell_nb_verts_dl(sjs_s *sjs) {
  static ivec2 cell_nb_verts_offsets[NUM_CELL_NB_VERTS] = {
    {.i = -1, .j = -1},
    {.i = -1, .j =  0},
    {.i = -1, .j =  1},
    {.i =  0, .j = -1},
    {.i =  0, .j =  0},
    {.i =  0, .j =  1},
    {.i =  1, .j = -1},
    {.i =  1, .j =  0},
    {.i =  1, .j =  1}
  };
  for (int i = 0; i < NUM_CELL_NB_VERTS; ++i) {
    sjs->cell_nb_verts_dl[i] = ind2l(sjs->shape, cell_nb_verts_offsets[i]);
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
  /**
   * This array gives the offsets ind "indc" space from the linear
   * index of a grid node to its four neighboring cells.
   */
  static ivec2 nb_cell_offsets[NUM_NB_CELLS] = {
    {.i = -1, .j = -1},
    {.i =  0, .j = -1},
    {.i = -1, .j =  0},
    {.i =  0, .j =  0}
  };
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

static void line(sjs_s *sjs, int l, int l0) {
  dvec2 xy = get_xy(sjs, l);
  dvec2 xy0 = get_xy(sjs, l0);
  dvec2 xy0_minus_xy = dvec2_sub(xy0, xy);
  dbl L = dvec2_norm(xy0_minus_xy);

  dbl T0 = sjs->jets[l0].f;
  dbl T = T0 + L;

  // Check causality
  assert(T > sjs->jets[l0].f);

  jet_s *jet = &sjs->jets[l];
  if (T < jet->f) {
    jet->f = T;
    jet->fx = -xy0_minus_xy.x/L;
    jet->fy = -xy0_minus_xy.y/L;
  }
}

/**
 * In this function, `ic0` is used as an index to select a nearby
 * bicubic interpolant which will be used to approximate `T`
 * locally.
 *
 * If the cell being indexed by ic0 is out of bounds, or if the cell
 * is invalid, this function does nothing.
 */
static void tri(sjs_s *sjs, int l, int l0, int l1, int ic0) {
  assert(ic0 >= 0);
  assert(ic0 < NUM_NB);

  int lc = l2lc(sjs->shape, l) + sjs->tri_dlc[ic0];
  if (lc < 0 || sjs->ncells <= lc) {
    return;
  }

  bicubic_s *bicubic = &sjs->bicubics[lc];
  if (!bicubic_valid(bicubic)) {
    return;
  }

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

  // Check causality
  assert(T > sjs->jets[l0].f);
  assert(T > sjs->jets[l1].f);

  jet_s *jet = &sjs->jets[l];
  if (T < jet->f) {
    dvec2 xy = get_xy(sjs, l);
    dvec2 xylam = dvec2_ccomb(data.xy0, data.xy1, lam);
    dvec2 xylam_minus_xy = dvec2_sub(xylam, xy);
    dbl L = dvec2_norm(xylam_minus_xy);
    jet->f = T;
    jet->fx = -xylam_minus_xy.x/L;
    jet->fy = -xylam_minus_xy.y/L;
  }
}

/**
 * Returns whether the linear index `l` corresponds to a valid index
 * into `sjs`'s domain.
 */
static bool inbounds(sjs_s const *sjs, int l) {
  return 0 <= l && l < sjs->nnodes;
}

static bool can_build_cell(sjs_s const *sjs, int lc) {
  // TODO: do this using SIMD gathers
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(sjs->shape, lc) + sjs->vert_dl[i];
    if (!inbounds(sjs, l)) {
      return false;
    }
    state_e state = sjs->states[l];
    /**
     * TODO: we don't want to build cells that only have trial values,
     * I don't think...
     */
    if (state != VALID) {
      return false;
    }
  }
  return true;
}

/**
 * TODO: we should describe precisely what "valid" means for a
 * cell. Right now or definition is just "the cell is inbounds and all
 * of its incident vertices are VALID". If all of the incident
 * vertices are VALID, then T, Tx, and Ty should all be finite for
 * each incident vertex, which is enough to bilinearly interpolate Txy
 * values.
 */
static bool cell_is_valid(sjs_s const *sjs, int lc) {
  if (!(0 <= lc && lc <= sjs->ncells)) {
    return false;
  }
  int l = lc2l(sjs->shape, lc);
  for (int iv = 0, lv; iv < NUM_CELL_VERTS; ++iv) {
    lv = l + sjs->vert_dl[iv];
    if (sjs->states[lv] != VALID) {
      return false;
    }
  }
  return true;
}

static dvec4 interpolate_Txy_at_verts(sjs_s *sjs, int lc) {
  /**
   * TODO: this probably works, but we'll use the original thing just
   * to be safe... Once we're computing good `fxy` values we can turn
   * this back on
   */
  // dvec4 fx, fy;
  // int l = lc2l(sjs->shape, lc);
  // for (int iv = 0, lv; iv < NUM_CELL_VERTS; ++iv) {
  //   lv = l + sjs->vert_dl[iv];
  //   fx.data[iv] = sjs->jets[lv].fx;
  //   fy.data[iv] = sjs->jets[lv].fy;
  // }
  // return interpolate_fxy_at_verts(fx, fy, sjs->h);

  /**
   * ... we'll just use this for now since I have more confidence that
   * it *does* work
   */

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

  dvec4 Txy;

  for (int i = 0; i < 4; ++i) {
    lam = lams[i];
    mu = mus[i];
    Txy.data[i] = (1 - mu)*((1 - lam)*fxy[0] + lam*fxy[1]) +
      mu*((1 - lam)*fxy[2] + lam*fxy[3]);
  }

  return Txy;
}

static dvec4 get_cell_Txy_values(sjs_s const *sjs, int lc) {
  dvec4 Txy;
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(sjs->shape, lc) + sjs->vert_dl[i];
    Txy.data[i] = sjs->jets[l].fxy;
    assert(isfinite(Txy.data[i]));
  }
  return Txy;
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

static void update(sjs_s *sjs, int l) {
  for (int i0 = 1, l0, l1, ic0; i0 < 8; i0 += 2) {
    l0 = l + sjs->nb_dl[i0];
    if (!inbounds(sjs, l0) || sjs->states[l0] != VALID) {
      continue;
    }

    l1 = l + sjs->nb_dl[i0 - 1];
    if (inbounds(sjs, l1) && sjs->states[l1] == VALID) {
      ic0 = i0 - 1;
      tri(sjs, l, l0, l1, ic0);
    }

    l1 = l + sjs->nb_dl[i0 + 1];
    if (inbounds(sjs, l1) && sjs->states[l1] == VALID) {
      ic0 = i0;
      tri(sjs, l, l0, l1, ic0);
    }
  }

  for (int i0 = 0, l0; i0 < 8; ++i0) {
    l0 = l + sjs->nb_dl[i0];
    if (inbounds(sjs, l0) && sjs->states[l0] == VALID) {
      line(sjs, l, l0);
    }
  }
}

static void adjust(sjs_s *sjs, int l0) {
  assert(sjs->states[l0] == TRIAL);
  assert(l0 >= 0);
  assert(l0 < sjs->nnodes);

  heap_swim(sjs->heap, sjs->positions[l0]);
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
  set_cell_nb_verts_dl(sjs);
  set_vert_dl(sjs);
  set_tri_dlc(sjs);
  set_nb_dlc(sjs);
  set_nearby_dlc(sjs);

  for (int lc = 0; lc < sjs->ncells; ++lc) {
    bicubic_invalidate(&sjs->bicubics[lc]);
  }

  for (int l = 0; l < sjs->nnodes; ++l) {
    sjs->jets[l].f = INFINITY;
    sjs->jets[l].fx = NAN;
    sjs->jets[l].fy = NAN;
    sjs->jets[l].fxy = NAN;
  }

  for (int l = 0; l < sjs->nnodes; ++l) {
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

#if SJS_DEBUG
static void check_cell_consistency(sjs_s const *sjs, int l0) {
  dbl tol = 1e-10, h = sjs->h, h_sq = h*h, f, fx, fy, fxy;
  dvec2 cc[4] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};
  bicubic_s *bicubic;
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[ic];
    if (can_build_cell(sjs, lc)) {
      bicubic = &sjs->bicubics[lc];
      for (int jv = 0, l; jv < NUM_CELL_VERTS; ++jv) {
        l = lc2l(sjs->shape, lc) + sjs->vert_dl[jv];
        f = bicubic_f(bicubic, cc[jv]);
        fx = bicubic_fx(bicubic, cc[jv]);
        fy = bicubic_fy(bicubic, cc[jv]);
        fxy = bicubic_fxy(bicubic, cc[jv]);
        assert(fabs(f - sjs->jets[l].f) < tol);
        assert(fabs(fx - h*sjs->jets[l].fx) < tol);
        assert(fabs(fy - h*sjs->jets[l].fy) < tol);
        assert(fabs(fxy - h_sq*sjs->jets[l].fxy) < tol);
      }
    }
  }
}
#endif

void sjs_step(sjs_s *sjs) {
  int l0 = heap_front(sjs->heap);
  assert(sjs->states[l0] == TRIAL);
  heap_pop(sjs->heap);
  sjs->states[l0] = VALID;

  // Determine which of the cells surrounding l0 are now valid. It's
  // enough to check if any of the four nearest cells are valid: it's
  // impossible for them to have been valid (or built before), since
  // one of their vertices just became valid.
  bool valid_cell_nb[NUM_NB_CELLS];
  for (int ic = 0, lc; ic < NUM_NB_CELLS; ++ic) {
    lc = l2lc(sjs->shape, l0) + sjs->nb_dlc[ic];
    valid_cell_nb[ic] = cell_is_valid(sjs, lc);
  }

  // TODO: don't build this inline?
  bool nb_incident_on_valid_cell_nb[9] = {
    valid_cell_nb[0],
    valid_cell_nb[0] || valid_cell_nb[2],
    valid_cell_nb[2],
    valid_cell_nb[0] || valid_cell_nb[1],
    valid_cell_nb[0] || valid_cell_nb[1] || valid_cell_nb[2] || valid_cell_nb[3],
    valid_cell_nb[2] || valid_cell_nb[3],
    valid_cell_nb[1],
    valid_cell_nb[1] || valid_cell_nb[3],
    valid_cell_nb[3]
  };

  // Array used to indicated which nearby cells should be used to
  // compute average values for Txy at cell vertices.
  bool use_for_Txy_average[NUM_NEARBY_CELLS];
  memset(use_for_Txy_average, 0x0, sizeof(bool)*NUM_NEARBY_CELLS);

  // Traverse the nearby cells and figure out which should be used to
  // compute new Txy values.
  //
  // NOTE: `use_for_Txy_average` => the cell is inbounds
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    for (int j = 0, i; j < NUM_CELL_VERTS; ++j) {
      i = cell_verts_to_cell_nb_verts[ic][j];
      if (i == NO_INDEX) {
        continue;
      }
      use_for_Txy_average[ic] |= nb_incident_on_valid_cell_nb[i];
    }
    lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[ic];
    use_for_Txy_average[ic] &= cell_is_valid(sjs, lc);
  }

  // Compute new Txy values at the vertices of the cells that we
  // decided to use.
  dvec4 Txy[NUM_NEARBY_CELLS];
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    if (use_for_Txy_average[ic]) {
      lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[ic];
      // If the cell is one of `l0`'s neighbors, then we have to use
      // bilinear extrapolation to compute its Txy values. Otherwise,
      // we can just grab the cell's existing Txy values.
      Txy[ic] = nearby_cell_is_cell_nb[ic] ?
        interpolate_Txy_at_verts(sjs, lc) :
        get_cell_Txy_values(sjs, lc);
    }
  }

  // TODO: document this
  {
    dbl Txy_sum;
    int nterms;
    for (int i = 0, l; i < NUM_CELL_NB_VERTS; ++i) {
      Txy_sum = 0;
      nterms = 0;
      for (int j = 0, jc, lc; j < NUM_NB_CELLS; ++j) {
        jc = cell_nb_verts_to_nearby_cells[i][j];
        if (use_for_Txy_average[jc]) {
          lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[jc];
          // We effectively access Txy[jc] in reverse here to
          // efficiently get the Txy value for cell `jc` at the vertex
          // indexed by `i`. See cell_notes.pdf for a visual depiction
          // of why this works.
          Txy_sum += Txy[jc].data[j];
          ++nterms;
        }
      }
      if (nterms > 0) {
        l = l0 + sjs->cell_nb_verts_dl[i];
        sjs->jets[l].fxy = Txy_sum/nterms;
      }
    }
  }

  // Finally, rebuild the cells! Our criterion for whether we should
  // rebuild a cell is simply whether we recomputed one of its Txy
  // values, so we can just check `use_for_Txy_average` here.
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    if (use_for_Txy_average[ic]) {
      lc = l2lc(sjs->shape, l0) + sjs->nearby_dlc[ic];
      build_cell(sjs, lc);
    }
  }

#if SJS_DEBUG
  check_cell_consistency(sjs, l0);
#endif

  /**
   * The section below corresponds to what's done in a "normal"
   * Dijkstra-like algorithm for solving the eikonal equation.
   */

  // Set FAR nodes to TRIAL and insert them into the heap.
  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_dl[i];
    if (inbounds(sjs, l) && sjs->states[l] == FAR) {
      sjs->states[l] = TRIAL;
      heap_insert(sjs->heap, l);
    }
  }

  // Update neighboring nodes.
  for (int i = 0, l; i < NUM_NB; ++i) {
    l = l0 + sjs->nb_dl[i];
    if (inbounds(sjs, l) && sjs->states[l] == TRIAL) {
      update(sjs, l);
      adjust(sjs, l);
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
