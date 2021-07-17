#include "eik.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "eik_F3.h"
#include "eik_F4.h"
#include "eik_S4.h"
#include "hybrid.h"
#include "index.h"
#include "jet.h"

#define NUM_CELL_VERTS 4
#define NUM_CELL_NB_VERTS 9
#define NUM_NB 8
#define NUM_NB_CELLS 4
#define NUM_NEARBY_CELLS 16

/**
 * TODO: add a few words about what `eik` is and how it works
 *
 * TODO: add a few comments about how indexing works
 *
 * - row-major ordering is used
 * - l vs lc index spaces
 * - cell verts are in column major order
 */
struct eik {
  field2_s const *slow;
  int2 shape;
  dbl2 xymin;
  dbl h;
  int nnodes, ncells;
  int nb_dl[NUM_NB + 1];
  int cell_nb_verts_dl[NUM_CELL_NB_VERTS];
  int vert_dl[NUM_CELL_VERTS];
  int tri_dlc[NUM_NB];
  int nb_dlc[NUM_NB_CELLS];
  int nearby_dlc[NUM_NEARBY_CELLS];
  bicubic_s *bicubics;
  jet2 *jets;
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
static int2 offsets[NUM_NB + 1] = {
  {-1, -1},
  {-1,  0},
  {-1,  1},
  { 0,  1},
  { 1,  1},
  { 1,  0},
  { 1, -1},
  { 0, -1},
  {-1, -1}
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
static int2 cell_vert_offsets[NUM_CELL_VERTS] = {
  {0, 0},
  {1, 0},
  {0, 1},
  {1, 1}
};

/**
 * The array `tri_cell_offsets` contains the `ivec2` offsets that are
 * used to initialize the `tri_dlc` member `eik_s`, which is in turn
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
static int2 tri_cell_offsets[NUM_NB] = {
  {-2, -1},
  {-2,  0},
  {-1,  1},
  { 0,  1},
  { 1,  0},
  { 1, -1},
  { 0, -2},
  {-1, -2}
};

/**
 * TODO: document this
 *
 * - used to set nearby_dlc
 */
static int2 nearby_cell_offsets[NUM_NEARBY_CELLS] = {
  {-2, -2}, {-2, -1}, {-2, 0}, {-2, 1},
  {-1, -2}, {-1, -1}, {-1, 0}, {-1, 1},
  { 0, -2}, { 0, -1}, { 0, 0}, { 0, 1},
  { 1, -2}, { 1, -1}, { 1, 0}, { 1, 1},
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

/**
 * This array gives the offsets ind "indc" space from the linear
 * index of a grid node to its four neighboring cells.
 */
static int2 nb_cell_offsets[NUM_NB_CELLS] = {
    {-1, -1},
    { 0, -1},
    {-1,  0},
    { 0,  0}
};

static bool nearby_cell_is_cell_nb[NUM_NEARBY_CELLS] = {
  false, false, false, false,
  false,  true,  true, false,
  false,  true,  true, false,
  false, false, false, false
};

static void set_nb_dl(eik_s *eik) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    eik->nb_dl[i] = ind2l(eik->shape, offsets[i]);
  }
}

static void set_cell_nb_verts_dl(eik_s *eik) {
  static int2 cell_nb_verts_offsets[NUM_CELL_NB_VERTS] = {
    {-1, -1},
    {-1,  0},
    {-1,  1},
    { 0, -1},
    { 0,  0},
    { 0,  1},
    { 1, -1},
    { 1,  0},
    { 1,  1}
  };
  for (int i = 0; i < NUM_CELL_NB_VERTS; ++i) {
    eik->cell_nb_verts_dl[i] = ind2l(eik->shape, cell_nb_verts_offsets[i]);
  }
}

static void set_vert_dl(eik_s *eik) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    eik->vert_dl[i] = ind2l(eik->shape, cell_vert_offsets[i]);
  }
}

static void set_tri_dlc(eik_s *eik) {
  for (int i = 0; i < NUM_NB; ++i) {
    eik->tri_dlc[i] = ind2lc(eik->shape, tri_cell_offsets[i]);
  }
}

static void set_nb_dlc(eik_s *eik) {
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    eik->nb_dlc[i] = ind2lc(eik->shape, nb_cell_offsets[i]);
  }
}

static void set_nearby_dlc(eik_s *eik) {
  for (int i = 0; i < NUM_NEARBY_CELLS; ++i) {
    eik->nearby_dlc[i] = ind2lc(eik->shape, nearby_cell_offsets[i]);
  }
}

static void get_xy(eik_s *eik, int l, dbl2 xy) {
  int2 ind; l2ind(eik->shape, l, ind);
  xy[0] = eik->h*ind[0] + eik->xymin[0];
  xy[1] = eik->h*ind[1] + eik->xymin[1];
}

dbl S4_th(dbl th, void *data) {
  S4_context *context = (S4_context *)data;
  S4_compute(th, context);
  return context->S4_th;
}

static void line(eik_s *eik, int l, int l0) {
  jet2 const *J0 = &eik->jets[l0];
  dbl T0 = J0->f;
  dbl Tx0 = J0->fx;
  dbl Ty0 = J0->fy;

  dbl2 xy; get_xy(eik, l, xy);
  dbl2 xy0; get_xy(eik, l0, xy0);

  S4_context context;
  context.slow = eik->slow;
  context.s = field2_f(eik->slow, xy);
  context.s0 = field2_f(eik->slow, xy0);
  dbl2_sub(xy, xy0, context.lp);
  context.L = dbl2_norm(context.lp);
  dbl2_dbl_div_inplace(context.lp, context.L);
  dbl2_avg(xy, xy0, context.xy_xy0_avg);
  dbl2_normalized(DBL2(Tx0, Ty0), context.t0);

  dbl th = atan2(context.lp[1], context.lp[0]);
  {
    dbl th_min = th - PI_OVER_FOUR;
    dbl th_max = th + PI_OVER_FOUR;
    bool found = hybrid(S4_th, th_min, th_max, (void *)&context, &th);
#if JMM_DEBUG
    assert(found);
#else
    (void)found;
#endif
  }

  dbl T = T0 + context.L*context.S4;

  // Check causality
  assert(T > T0);

  jet2 *J = &eik->jets[l];
  if (T < J->f) {
    J->f = T;
    J->fx = context.s*cos(th);
    J->fy = context.s*sin(th);
  }
}

dbl F3_eta(dbl eta, void *data) {
  F3_context *context = (F3_context *)data;
  F3_compute(eta, context);
  return context->F3_eta;
}

/**
 * In this function, `ic0` is used as an index to select a nearby
 * bicubic interpolant which will be used to approximate `T`
 * locally.
 *
 * If the cell being indexed by ic0 is out of bounds, or if the cell
 * is invalid, this function does nothing.
 */
static void tri(eik_s *eik, int l, int l0, int l1, int ic0) {
  assert(ic0 >= 0);
  assert(ic0 < NUM_NB);

  int lc = l2lc(eik->shape, l) + eik->tri_dlc[ic0];
  if (lc < 0 || eik->ncells <= lc) {
    return;
  }

  bicubic_s *bicubic = &eik->bicubics[lc];
  if (!bicubic_valid(bicubic)) {
    return;
  }

  /**
   * Get cubic along edge of interest.
   */
  bicubic_variable var = tri_bicubic_vars[ic0];
  int edge = tri_edges[ic0];
  cubic_s T_cubic = bicubic_get_f_on_edge(bicubic, var, edge);
  cubic_s Tx_cubic = bicubic_get_fx_on_edge(bicubic, var, edge);
  cubic_s Ty_cubic = bicubic_get_fy_on_edge(bicubic, var, edge);
  if (should_reverse_cubic[ic0]) {
    cubic_reverse_on_unit_interval(&T_cubic);
    cubic_reverse_on_unit_interval(&Tx_cubic);
    cubic_reverse_on_unit_interval(&Ty_cubic);
  }

  dbl2 xy; get_xy(eik, l, xy);
  dbl2 xy0; get_xy(eik, l0, xy0);
  dbl2 xy1; get_xy(eik, l1, xy1);

  // TODO: try initializing from the mp0 minimizer since it's so cheap
  // to compute...

  /**
   * Compute initial guess for eta and theta by minimizing F3.
   */
  dbl eta, th;
  {
    F3_context context;
    context.T_cubic = T_cubic;
    dbl2_copy(xy, context.xy);
    dbl2_copy(xy0, context.xy0);
    dbl2_copy(xy1, context.xy1);
    context.slow = eik->slow;

    bool found = hybrid(F3_eta, 0, 1, (void *)&context, &eta);
#if JMM_DEBUG
    assert(found);
#else
    (void)found;
#endif
  }
  {
    dbl2 dxy; dbl2_sub(xy1, xy0, dxy);
    dbl2 dxy_times_eta; dbl2_dbl_mul(dxy, eta, dxy_times_eta);
    dbl2 xyeta; dbl2_add(xy0, dxy_times_eta, xyeta);
    dbl2 lp; dbl2_sub(xy, xyeta, lp);
    dbl2_normalize(lp);
    th = atan2(lp[1], lp[0]);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Minimize F4 in this section.
  //
  // TODO: this is still a little rough and unfinished (hence the
  // "construction signs" delimiting this section). We want to remove
  // the print statements below and decide on a robust error-handling
  // strategy so that `abort` is never called.

  dbl T = NAN;

  {
    F4_context context;
    context.T_cubic = T_cubic;
    context.Tx_cubic = Tx_cubic;
    context.Ty_cubic = Ty_cubic;
    dbl2_copy(xy, context.xy);
    dbl2_copy(xy0, context.xy0);
    dbl2_copy(xy1, context.xy1);
    context.slow = eik->slow;

    dbl2 xk, xk1, gk, gk1;
    dbl22 Hk, Hk1;
    F4_bfgs_init(eta, th, xk, gk, Hk, &context);
    dbl Tprev = context.F4;

    int iter = 0;
    while (F4_bfgs_step(xk, gk, Hk, xk1, gk1, Hk1, &context)) {
      if (xk1[0] < 0 || xk1[0] > 1) {
        printf("out of bounds: eta = %g\n", xk[0]);
        abort();
      }

      T = context.F4;

      ++iter;

      if (fabs(T - Tprev) <= EPS*fabs(fmax(T, Tprev)) + EPS) {
        break;
      }

      if (dbl2_maxnorm(gk1) <= EPS) {
        break;
      }

      if (iter >= 20) {
        printf("exceeded number of iterations\n");
        abort();
      }
      if (T > Tprev) {
        printf("no decrease in T: (%g > %g)\n", T, Tprev);
        abort();
      }

      Tprev = T;

      dbl2_copy(xk1, xk);
      dbl2_copy(gk1, gk);
      dbl22_copy(Hk1, Hk);
    }

    th = xk[1];
  }

  //////////////////////////////////////////////////////////////////////////////

  // Check causality
  assert(T > eik->jets[l0].f);
  assert(T > eik->jets[l1].f);

  /**
   * Commit new value if it's an improvement.
   */
  jet2 *jet = &eik->jets[l];
  if (T < jet->f) {
    jet->f = T;

    dbl s = field2_f(eik->slow, xy);
    jet->fx = s*cos(th);
    jet->fy = s*sin(th);
  }
}

/**
 * Returns whether the linear index `l` corresponds to a valid index
 * into `eik`'s domain.
 */
static bool inbounds(eik_s const *eik, int2 ind) {
  return 0 <= ind[0] && ind[0] < eik->shape[0] &&
         0 <= ind[1] && ind[1] < eik->shape[1];
}

static bool can_build_cell(eik_s const *eik, int lc) {
  // TODO: do this using SIMD gathers
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(eik->shape, lc) + eik->vert_dl[i];
    // TODO: this is probably wrong: need to double check what's in
    // vert_dl to be sure---it's probably wrapping around... although
    // this may not result in an actual bug?
    int2 ind; l2ind(eik->shape, l, ind);
    if (!inbounds(eik, ind)) {
      return false;
    }
    state_e state = eik->states[l];
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
static bool cell_is_valid(eik_s const *eik, int2 indc) {
  if (!(0 <= indc[0] && indc[0] < eik->shape[0] - 1 &&
		0 <= indc[1] && indc[1] < eik->shape[1] - 1)) {
    return false;
  }
  int l = indc2l(eik->shape, indc);
  for (int iv = 0, lv; iv < NUM_CELL_VERTS; ++iv) {
    int2 indv; int2_add(indc, cell_vert_offsets[iv], indv);
    if (!(0 <= indv[0] && indv[0] < eik->shape[0] &&
		  0 <= indv[1] && indv[1] < eik->shape[1])) {
      return false;
    }
    lv = l + eik->vert_dl[iv];
    if (eik->states[lv] != VALID) {
      return false;
    }
  }
  return true;
}

static void interpolate_Txy_at_verts(eik_s *eik, int lc, dbl4 Txy) {
  /**
   * TODO: this probably works, but we'll use the original thing just
   * to be safe... Once we're computing good `fxy` values we can turn
   * this back on
   */
  // dvec4 fx, fy;
  // int l = lc2l(eik->shape, lc);
  // for (int iv = 0, lv; iv < NUM_CELL_VERTS; ++iv) {
  //   lv = l + eik->vert_dl[iv];
  //   fx.data[iv] = eik->jets[lv].fx;
  //   fy.data[iv] = eik->jets[lv].fy;
  // }
  // return interpolate_fxy_at_verts(fx, fy, eik->h);

  /**
   * ... we'll just use this for now since I have more confidence that
   * it *does* work
   */

  dbl fx[NUM_CELL_VERTS], fy[NUM_CELL_VERTS];

  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(eik->shape, lc) + eik->vert_dl[i];
    fx[i] = eik->jets[l].fx;
    fy[i] = eik->jets[l].fy;
  }

  dbl fxy[NUM_CELL_VERTS] = {
    (fy[1] - fy[0])/eik->h, // left
    (fx[3] - fx[1])/eik->h, // bottom
    (fx[2] - fx[0])/eik->h, // top
    (fy[3] - fy[2])/eik->h // right
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

static void get_cell_Txy_values(eik_s const *eik, int lc, dbl4 Txy) {
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = lc2l(eik->shape, lc) + eik->vert_dl[i];
    Txy[i] = eik->jets[l].fxy;
    assert(isfinite(Txy[i]));
  }
}

/**
 * Build the cell at index `lc`. This uses the jets at the cells
 * vertices to assemble the data matrix for the bicubic
 * interpolant. This function doesn't make any assumptions about the
 * state of the nodes, so it's assumed that the caller has already
 * made sure this is a reasonable thing to try to do.
 */
static void build_cell(eik_s *eik, int lc) {
  /* Get linear indices of cell vertices */
  int l[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    l[i] = lc2l(eik->shape, lc) + eik->vert_dl[i];
  }

  /* Get jet at each cell vertex */
  jet2 *J[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    J[i] = &eik->jets[l[i]];
  }

  /* Precompute scaling factors for partial derivatives */
  dbl h = eik->h, h_sq = h*h;

  /* Compute cell data from jets and scaling factors */
  dbl44 data;
  data[0][0] = J[0]->f;
  data[1][0] = J[1]->f;
  data[0][1] = J[2]->f;
  data[1][1] = J[3]->f;
  data[2][0] = h*J[0]->fx;
  data[3][0] = h*J[1]->fx;
  data[2][1] = h*J[2]->fx;
  data[3][1] = h*J[3]->fx;
  data[0][2] = h*J[0]->fy;
  data[1][2] = h*J[1]->fy;
  data[0][3] = h*J[2]->fy;
  data[1][3] = h*J[3]->fy;
  data[2][2] = h_sq*J[0]->fxy;
  data[3][2] = h_sq*J[1]->fxy;
  data[2][3] = h_sq*J[2]->fxy;
  data[3][3] = h_sq*J[3]->fxy;

  /* Set cell data */
  bicubic_set_data(&eik->bicubics[lc], data);
}

static void update(eik_s *eik, int l) {
  /**
   * First, precompute whether each neighboring node is inbounds. We
   * need to do this in the (i, j) index space to avoid wrapping
   * errors.
   */
  bool inbounds_[9];
  {
    int2 ind; l2ind(eik->shape, l, ind);
    for (int i0 = 0; i0 < 9; ++i0) {
      int2 ind0; int2_add(ind, offsets[i0], ind0);
      inbounds_[i0] = inbounds(eik, ind0);
    }
  }

  for (int i0 = 1, l0, l1, ic0; i0 < 8; i0 += 2) {
    if (!inbounds_[i0]) {
      continue;
    }

    l0 = l + eik->nb_dl[i0];
    if (eik->states[l0] != VALID) {
      continue;
    }

    if (inbounds_[i0 - 1]) {
      l1 = l + eik->nb_dl[i0 - 1];
      if (eik->states[l1] == VALID) {
        ic0 = i0 - 1;
        tri(eik, l, l0, l1, ic0);
      }
    }

    if (inbounds_[i0 + 1]) {
      l1 = l + eik->nb_dl[i0 + 1];
      if (eik->states[l1] == VALID) {
        ic0 = i0;
        tri(eik, l, l0, l1, ic0);
      }
    }
  }

  for (int i0 = 0, l0; i0 < 8; ++i0) {
    if (inbounds_[i0]) {
      l0 = l + eik->nb_dl[i0];
      if (eik->states[l0] == VALID) {
        line(eik, l, l0);
      }
    }
  }
}

static void adjust(eik_s *eik, int l0) {
  assert(eik->states[l0] == TRIAL);
  assert(l0 >= 0);
  assert(l0 < eik->nnodes);

  heap_swim(eik->heap, eik->positions[l0]);
}

static dbl value(void *vp, int l) {
  eik_s *eik = (eik_s *)vp;
  assert(l >= 0);
  assert(l < eik->nnodes);
  dbl T = eik->jets[l].f;
  return T;
}

static void setpos(void *vp, int l, int pos) {
  eik_s *eik = (eik_s *)vp;
  eik->positions[l] = pos;
}

void eik_alloc(eik_s **eik) {
  *eik = malloc(sizeof(eik_s));
  assert(*eik != NULL);
}

void eik_dealloc(eik_s **eik) {
  free(*eik);
  *eik = NULL;
}

// TODO: since the margins are BOUNDARY nodes, we actually don't need
// to allocate an extra margin of cells, since they will never be
// initialized (i.e., they will never have all of their vertex nodes
// become VALID because of the margin...)
void eik_init(eik_s *eik, field2_s const *slow, int2 shape, dbl2 xymin, dbl h) {
  eik->slow = slow;
  eik->ncells = (shape[0] - 1)*(shape[1] - 1);
  eik->nnodes = shape[0]*shape[1];
  eik->h = h;
  eik->bicubics = malloc(eik->ncells*sizeof(bicubic_s));
  eik->jets = malloc(eik->nnodes*sizeof(jet2));
  eik->states = malloc(eik->nnodes*sizeof(state_e));
  eik->positions = malloc(eik->nnodes*sizeof(int));

  int2_copy(shape, eik->shape);
  dbl2_copy(xymin, eik->xymin);

  assert(eik->bicubics != NULL);
  assert(eik->jets != NULL);
  assert(eik->states != NULL);
  assert(eik->positions != NULL);

#if SJS_DEBUG
  for (int l = 0; l < eik->nnodes; ++l) {
    eik->positions[l] = NO_INDEX;
  }
#endif

  heap_alloc(&eik->heap);

  int capacity = (int) 3*sqrt(eik->shape[0]*eik->shape[1]);
  heap_init(eik->heap, capacity, value, setpos, (void *)eik);

  set_nb_dl(eik);
  set_cell_nb_verts_dl(eik);
  set_vert_dl(eik);
  set_tri_dlc(eik);
  set_nb_dlc(eik);
  set_nearby_dlc(eik);

  for (int lc = 0; lc < eik->ncells; ++lc) {
    bicubic_invalidate(&eik->bicubics[lc]);
  }

  for (int l = 0; l < eik->nnodes; ++l) {
    eik->jets[l].f = INFINITY;
    eik->jets[l].fx = NAN;
    eik->jets[l].fy = NAN;
    eik->jets[l].fxy = NAN;
  }

  for (int l = 0; l < eik->nnodes; ++l) {
    eik->states[l] = FAR;
  }
}

void eik_deinit(eik_s *eik) {
  eik->slow = NULL;

  free(eik->bicubics);
  free(eik->jets);
  free(eik->states);
  free(eik->positions);

  eik->bicubics = NULL;
  eik->jets = NULL;
  eik->states = NULL;
  eik->positions = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);
}

#if SJS_DEBUG
static void check_cell_consistency(eik_s const *eik, int l0) {
  dbl tol = 1e-10, h = eik->h, h_sq = h*h, f, fx, fy, fxy;
  dbl2 cc[4] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};
  bicubic_s *bicubic;
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    lc = l2lc(eik->shape, l0) + eik->nearby_dlc[ic];
    if (can_build_cell(eik, lc)) {
      bicubic = &eik->bicubics[lc];
      for (int jv = 0, l; jv < NUM_CELL_VERTS; ++jv) {
        l = lc2l(eik->shape, lc) + eik->vert_dl[jv];
        f = bicubic_f(bicubic, cc[jv]);
        fx = bicubic_fx(bicubic, cc[jv]);
        fy = bicubic_fy(bicubic, cc[jv]);
        fxy = bicubic_fxy(bicubic, cc[jv]);
        assert(fabs(f - eik->jets[l].f) < tol);
        assert(fabs(fx - h*eik->jets[l].fx) < tol);
        assert(fabs(fy - h*eik->jets[l].fy) < tol);
        assert(fabs(fxy - h_sq*eik->jets[l].fxy) < tol);
      }
    }
  }
}
#endif

void eik_step(eik_s *eik) {
  int l0 = heap_front(eik->heap);
  assert(eik->states[l0] == TRIAL);
  heap_pop(eik->heap);
  eik->states[l0] = VALID;

  int2 ind0; l2ind(eik->shape, l0, ind0);

  // Determine which of the cells surrounding l0 are now valid. It's
  // enough to check if any of the four nearest cells are valid: it's
  // impossible for them to have been valid (or built before), since
  // one of their vertices just became valid.
  bool valid_cell_nb[NUM_NB_CELLS];
  for (int ic = 0; ic < NUM_NB_CELLS; ++ic) {
    int2 indc; int2_add(ind0, nb_cell_offsets[ic], indc);
    valid_cell_nb[ic] = cell_is_valid(eik, indc);
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
  for (int ic = 0; ic < NUM_NEARBY_CELLS; ++ic) {
    for (int j = 0, i; j < NUM_CELL_VERTS; ++j) {
      i = cell_verts_to_cell_nb_verts[ic][j];
      if (i == NO_INDEX) {
        continue;
      }
      use_for_Txy_average[ic] |= nb_incident_on_valid_cell_nb[i];
    }
	{
      int2 indc; int2_add(ind0, nearby_cell_offsets[ic], indc);
      use_for_Txy_average[ic] &= cell_is_valid(eik, indc);
	}
  }

  // Compute new Txy values at the vertices of the cells that we
  // decided to use.
  dbl4 Txy[NUM_NEARBY_CELLS];
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    if (use_for_Txy_average[ic]) {
      lc = l2lc(eik->shape, l0) + eik->nearby_dlc[ic];
      // If the cell is one of `l0`'s neighbors, then we have to use
      // bilinear extrapolation to compute its Txy values. Otherwise,
      // we can just grab the cell's existing Txy values.
      if (nearby_cell_is_cell_nb[ic])
        interpolate_Txy_at_verts(eik, lc, Txy[ic]);
      else
        get_cell_Txy_values(eik, lc, Txy[ic]);
    }
  }

  // TODO: document this
  {
    dbl Txy_sum;
    int nterms;
    for (int i = 0, l; i < NUM_CELL_NB_VERTS; ++i) {
      Txy_sum = 0;
      nterms = 0;
      for (int j = 0, jc; j < NUM_NB_CELLS; ++j) {
        jc = cell_nb_verts_to_nearby_cells[i][j];
        if (use_for_Txy_average[jc]) {
          // We effectively access Txy[jc] in reverse here to
          // efficiently get the Txy value for cell `jc` at the vertex
          // indexed by `i`. See cell_notes.pdf for a visual depiction
          // of why this works.
          Txy_sum += Txy[jc][j];
          ++nterms;
        }
      }
      if (nterms > 0) {
        l = l0 + eik->cell_nb_verts_dl[i];
        eik->jets[l].fxy = Txy_sum/nterms;
      }
    }
  }

  // Finally, rebuild the cells! Our criterion for whether we should
  // rebuild a cell is simply whether we recomputed one of its Txy
  // values, so we can just check `use_for_Txy_average` here.
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    if (use_for_Txy_average[ic]) {
      lc = l2lc(eik->shape, l0) + eik->nearby_dlc[ic];
      build_cell(eik, lc);
    }
  }

#if SJS_DEBUG
  check_cell_consistency(eik, l0);
#endif

  /**
   * The section below corresponds to what's done in a "normal"
   * Dijkstra-like algorithm for solving the eikonal equation.
   */

  // Set FAR nodes to TRIAL and insert them into the heap.
  for (int i = 0, l; i < NUM_NB; ++i) {
    int2 ind; int2_add(ind0, offsets[i], ind);
    if (!inbounds(eik, ind)) {
      continue;
    }
    l = l0 + eik->nb_dl[i];
    if (eik->states[l] == FAR) {
      eik->states[l] = TRIAL;
      heap_insert(eik->heap, l);
    }
  }

  // Update neighboring nodes.
  for (int i = 0, l; i < NUM_NB; ++i) {
    int2 ind; int2_add(ind0, offsets[i], ind);
    if (!inbounds(eik, ind)) {
      continue;
    }
    l = l0 + eik->nb_dl[i];
    if (eik->states[l] == TRIAL) {
      update(eik, l);
      adjust(eik, l);
    }
  }
}

void eik_solve(eik_s *eik) {
  while (heap_size(eik->heap) > 0) {
    eik_step(eik);
  }
}

void eik_add_trial(eik_s *eik, int2 ind, jet2 jet) {
  int l = ind2l(eik->shape, ind);
  eik->jets[l] = jet;
  assert(eik->states[l] != TRIAL && eik->states[l] != VALID);
  eik->states[l] = TRIAL;
  heap_insert(eik->heap, l);
}

void eik_add_valid(eik_s *eik, int2 ind, jet2 jet) {
  int l = ind2l(eik->shape, ind);
  eik->jets[l] = jet;
  assert(eik->states[l] != TRIAL && eik->states[l] != VALID);
  eik->states[l] = VALID;
}

void eik_make_bd(eik_s *eik, int2 ind) {
  int l = ind2l(eik->shape, ind);
  eik->states[l] = BOUNDARY;
}

void eik_get_shape(eik_s const *eik, int2 shape) {
  int2_copy(eik->shape, shape);
}

jet2 eik_get_jet(eik_s *eik, int2 ind) {
  int l = ind2l(eik->shape, ind);
  return eik->jets[l];
}

jet2 *eik_get_jets_ptr(eik_s const *eik) {
  return eik->jets;
}

state_e eik_get_state(eik_s const *eik, int2 ind) {
  int l = ind2l(eik->shape, ind);
  return eik->states[l];
}

state_e *eik_get_states_ptr(eik_s const *eik) {
  return eik->states;
}

/**
 * The four functions below (`eik_T`, `eik_Tx`, `eik_Ty`, and
 * `eik_Txy`) are only intended to be used by people consuming this
 * API, not internally.
 *
 * TODO: at some point, we're going to want to add batch versions of
 * these functions to avoid calling `can_build_cell` over and over.
 */

dbl eik_T(eik_s *eik, dbl2 xy) {
  dbl2 cc;
  int lc = xy_to_lc_and_cc(eik->shape, eik->xymin, eik->h, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_f(bicubic, cc);
}

dbl eik_Tx(eik_s *eik, dbl2 xy) {
  dbl2 cc;
  int lc = xy_to_lc_and_cc(eik->shape, eik->xymin, eik->h, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_fx(bicubic, cc)/eik->h;
}

dbl eik_Ty(eik_s *eik, dbl2 xy) {
  dbl2 cc;
  int lc = xy_to_lc_and_cc(eik->shape, eik->xymin, eik->h, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_fy(bicubic, cc)/eik->h;
}

dbl eik_Txy(eik_s *eik, dbl2 xy) {
  dbl2 cc;
  int lc = xy_to_lc_and_cc(eik->shape, eik->xymin, eik->h, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_fxy(bicubic, cc)/(eik->h*eik->h);
}

bool eik_can_build_cell(eik_s const *eik, int2 indc) {
  int lc = indc2lc(eik->shape, indc);
  return can_build_cell(eik, lc);
}

void eik_build_cells(eik_s *eik) {
  for (int lc = 0; lc < eik->ncells; ++lc) {
    if (can_build_cell(eik, lc)) {
      build_cell(eik, lc);
    }
  }
}

bicubic_s eik_get_bicubic(eik_s const *eik, int2 indc) {
  int lc = indc2lc(eik->shape, indc);
  return eik->bicubics[lc];
}

bicubic_s *eik_get_bicubics_ptr(eik_s const *eik) {
  return eik->bicubics;
}

heap_s *eik_get_heap(eik_s const *eik) {
  return eik->heap;
}
