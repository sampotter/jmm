#include <jmm/eik.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <jmm/eik_F3.h>
#include <jmm/eik_F4.h>
#include <jmm/eik_S4.h>
#include <jmm/hybrid.h>

#define NUM_CELL_VERTS 4
#define NUM_CELL_NB_VERTS 9
#define NUM_NB 8
#define NUM_NB_CELLS 4
#define NUM_NEARBY_CELLS 16

/* See the documentation for `grid2_s` for more information about how
   indexing and the coordinate system works, depending on whether row
   or column major indexing is used. */
struct eik {
  field2_s const *slow;
  grid2_s const *grid;
  int nnodes, ncells;
  int nb_dl[NUM_NB + 1];
  int cell_nb_verts_dl[NUM_CELL_NB_VERTS];
  int vert_dl[NUM_CELL_VERTS];
  int tri_dlc[NUM_NB];
  int nb_dlc[NUM_NB_CELLS];
  int nearby_dlc[NUM_NEARBY_CELLS];
  bicubic_s *bicubics;
  jet21p *jets;
  state_e *states;
  int *positions;
  heap_s *heap;
  size_t num_accepted;
  size_t *accepted;
  par2_s *par;
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
static int2 cell_vert_offsets[NUM_ORDERS][NUM_CELL_VERTS] = {
  [ORDER_ROW_MAJOR] = {
    {0, 0},
    {0, 1},
    {1, 0},
    {1, 1}
  },
  [ORDER_COLUMN_MAJOR] = {
    {0, 0},
    {1, 0},
    {0, 1},
    {1, 1}
  }
};


static void get_cell_vert_offset(eik_s const *eik, size_t iv, int2 offset) {
  assert(iv <= NUM_CELL_VERTS);
  memcpy(offset, cell_vert_offsets[eik->grid->order][iv], sizeof(int2));
}

static dbl2 cell_vert_coefs[NUM_ORDERS][NUM_CELL_VERTS] = {
  [ORDER_ROW_MAJOR] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}},
  [ORDER_COLUMN_MAJOR] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}}
};

static void get_cell_vert_coefs(eik_s const *eik, size_t iv, dbl2 cc) {
  assert(iv <= NUM_CELL_VERTS);
  memcpy(cc, cell_vert_coefs[eik->grid->order][iv], sizeof(dbl2));
}

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
 * This looks like:
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
static int2 tri_cell_offsets[NUM_ORDERS][NUM_NB] = {
  [ORDER_ROW_MAJOR] = {
    {-2, -1},
    {-2,  0},
    {-1,  1},
    { 0,  1},
    { 1,  0},
    { 1, -1},
    { 0, -2},
    {-1, -2}
  }
};

static void get_tri_cell_offset(eik_s const *eik, size_t i, int2 ind) {
  assert(i < NUM_NB);
  if (eik->grid->order != ORDER_ROW_MAJOR)
    abort();
  memcpy(ind, tri_cell_offsets[eik->grid->order][i], sizeof(int2));
}

/**
 * TODO: document this
 *
 * - used to set nearby_dlc
 */
static int2 nearby_cell_offsets[NUM_ORDERS][NUM_NEARBY_CELLS] = {
  [ORDER_ROW_MAJOR] = {
    {-2, -2}, {-2, -1}, {-2, 0}, {-2, 1},
    {-1, -2}, {-1, -1}, {-1, 0}, {-1, 1},
    { 0, -2}, { 0, -1}, { 0, 0}, { 0, 1},
    { 1, -2}, { 1, -1}, { 1, 0}, { 1, 1},
  }
};

static void get_nearby_cell_offset(eik_s const *eik, size_t ic, int2 offset) {
  assert(ic < NUM_NEARBY_CELLS);
  if (eik->grid->order != ORDER_ROW_MAJOR)
    abort();
  memcpy(offset, nearby_cell_offsets[eik->grid->order][ic], sizeof(int2));
}

/* This maps from each cell neighboring the update zone to the
 * variable along its neighboring edge, where f = f(lambda, mu) is the
 * bicubic for that cell. The assumption is that the bicubics'
 * coordinate system matches the order of the indexing. */
static bicubic_variable tri_bicubic_vars[NUM_ORDERS][NUM_NB] = {
  [ORDER_ROW_MAJOR] = {
    /* ic -> */ MU, MU, LAMBDA, LAMBDA, MU, MU, LAMBDA, LAMBDA
  },
  [ORDER_COLUMN_MAJOR] = {
    /* ic -> */ LAMBDA, LAMBDA, MU, MU, LAMBDA, LAMBDA, MU, MU
  }
};

static bicubic_variable get_tri_bicubic_var(eik_s const *eik, size_t ic) {
  assert(ic <= NUM_NB);
  return tri_bicubic_vars[eik->grid->order][ic];
}

/* This gives the value of the coordinate not selected by
 * `tri_bicubic_vars`, restricting the bicubic to the edge incident on
 * the update zone. E.g., if `tri_bicubic_vars == LAMBDA`, and
 * `tri_edges == 1`, then this results in the cubic `f(lambda, 1)`. */
static int tri_edges[NUM_ORDERS][NUM_NB] = {
  [ORDER_ROW_MAJOR] = {
    /* ic -> */ 1, 1, 0, 0, 0, 0, 1, 1
  }
};

static int get_tri_edge(eik_s const *eik, size_t ic) {
  assert(ic <= NUM_NB);
  if (eik->grid->order != ORDER_ROW_MAJOR)
    abort();
  return tri_edges[eik->grid->order][ic];
}

/* This indicates whether the cubic selected by `tri_bicubic_vars` and
 * `tri_edges` above should be reversed so that the interval `[0, 1]`
 * spans the interval `[x[l0], x[l1]]`. */
static bool should_reverse_cubic[NUM_ORDERS][NUM_NB] = {
  [ORDER_ROW_MAJOR] = {
    /* ic -> */ true, false, true, false, false, true, false, true
  }
};

static bool get_should_reverse_cubic(eik_s const *eik, size_t ic) {
  assert(ic <= NUM_NB);
  if (eik->grid->order != ORDER_ROW_MAJOR)
    abort();
  return should_reverse_cubic[eik->grid->order][ic];
}

/* This lookup table provides a mapping from a nearby cell back to the
 * indices neighboring an update point. */
#define _ NO_INDEX
static int
cell_verts_to_cell_nb_verts[NUM_ORDERS][NUM_NEARBY_CELLS][NUM_CELL_VERTS] = {
  [ORDER_ROW_MAJOR] = {
    {_, _, _, 0}, {_, _, 0, 1}, {_, _, 1, 2}, {_, _, 2, _},
    {_, 0, _, 3}, {0, 1, 3, 4}, {1, 2, 4, 5}, {2, _, 5, _},
    {_, 3, _, 6}, {3, 4, 6, 7}, {4, 5, 7, 8}, {5, _, 8, _},
    {_, 6, _, _}, {6, 7, _, _}, {7, 8, _, _}, {8, _, _, _}
  },
  [ORDER_COLUMN_MAJOR] = {
    {_, _, _, 0}, {_, 0, _, 1}, {_, 1, _, 2}, {_, 2, _, _},
    {_, _, 0, 3}, {0, 3, 1, 4}, {1, 4, 2, 5}, {2, 5, _, _},
    {_, _, 3, 6}, {3, 6, 4, 7}, {4, 7, 5, 8}, {5, 8, _, _},
    {_, _, 6, _}, {6, _, 7, _}, {7, _, 8, _}, {8, _, _, _}
  }
};
#undef _

static int cell_vert_to_cell_nb_vert(eik_s const *eik, size_t ic, size_t jv) {
  return cell_verts_to_cell_nb_verts[eik->grid->order][ic][jv];
}

#define _ NO_INDEX
static int
cell_nb_verts_to_nearby_cells[NUM_ORDERS][NUM_CELL_NB_VERTS][NUM_NB_CELLS] = {
  [ORDER_ROW_MAJOR] = {
    { 5,  4,  1,  0},
    { 6,  5,  2,  1},
    { 7,  6,  3,  2},
    { 9,  8,  5,  4},
    {10,  9,  6,  5},
    {11, 10,  7,  6},
    {13, 12,  9,  8},
    {14, 13, 10,  9},
    {15, 14, 11, 10}
  }
};
#undef _

static int cell_nb_vert_to_nearby_cell(eik_s const *eik, size_t i, size_t j) {
  return cell_nb_verts_to_nearby_cells[eik->grid->order][i][j];
}

/**
 * This array gives the offsets in "indc" space from the linear
 * index of a grid node to its four neighboring cells.
 */
static int2 nb_cell_offsets[NUM_ORDERS][NUM_NB_CELLS] = {
  [ORDER_ROW_MAJOR] = {
    {-1, -1},
    {-1,  0},
    { 0, -1},
    { 0,  0}
  }
};

static void get_nb_cell_offset(eik_s const *eik, size_t ic, int2 indc) {
  assert(ic < NUM_NB_CELLS);
  memcpy(indc, nb_cell_offsets[eik->grid->order][ic], sizeof(int2));
}

static bool nearby_cell_is_cell_nb[NUM_NEARBY_CELLS] = {
  false, false, false, false,
  false,  true,  true, false,
  false,  true,  true, false,
  false, false, false, false
};

static void set_nb_dl(eik_s *eik) {
  for (int i = 0; i < NUM_NB + 1; ++i) {
    eik->nb_dl[i] = grid2_ind2l(eik->grid, offsets[i]);
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
    eik->cell_nb_verts_dl[i] = grid2_ind2l(eik->grid, cell_nb_verts_offsets[i]);
  }
}

static void set_vert_dl(eik_s *eik) {
  int2 offset;
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    get_cell_vert_offset(eik, i, offset);
    eik->vert_dl[i] = grid2_ind2l(eik->grid, offset);
  }
}

static void set_tri_dlc(eik_s *eik) {
  int2 offset;
  for (size_t i = 0; i < NUM_NB; ++i) {
    get_tri_cell_offset(eik, i, offset);
    eik->tri_dlc[i] = grid2_ind2lc(eik->grid, offset);
  }
}

static void set_nb_dlc(eik_s *eik) {
  int2 offset;
  for (size_t i = 0; i < NUM_CELL_VERTS; ++i) {
    get_nb_cell_offset(eik, i, offset);
    eik->nb_dlc[i] = grid2_ind2lc(eik->grid, offset);
  }
}

static void set_nearby_dlc(eik_s *eik) {
  int2 offset;
  for (size_t i = 0; i < NUM_NEARBY_CELLS; ++i) {
    get_nearby_cell_offset(eik, i, offset);
    eik->nearby_dlc[i] = grid2_ind2lc(eik->grid, offset);
  }
}

dbl S4_th(dbl th, void *data) {
  S4_context *context = (S4_context *)data;
  S4_compute(th, context);
  return context->S4_th;
}

dbl F3_eta(dbl eta, void *data) {
  F3_context *context = (F3_context *)data;
  F3_compute(eta, context);
  return context->F3_eta;
}

#if SJS_DEBUG
static void check_cell_consistency(eik_s const *eik, size_t lc) {
  dbl h = eik->grid->h, h_sq = h*h, f, fx, fy, fxy;

  dbl2 cc;
  bicubic_s const *bicubic = &eik->bicubics[lc];

  for (int iv = 0, l; iv < NUM_CELL_VERTS; ++iv) {
    get_cell_vert_coefs(eik, iv, cc);

    f = bicubic_f(bicubic, cc);
    fx = bicubic_fx(bicubic, cc);
    fy = bicubic_fy(bicubic, cc);
    fxy = bicubic_fxy(bicubic, cc);

    l = grid2_lc2l(eik->grid, lc) + eik->vert_dl[iv];

    assert(fabs(f - eik->jets[l].f) < EPS);
    assert(fabs(fx - h*eik->jets[l].Df[0]) < EPS);
    assert(fabs(fy - h*eik->jets[l].Df[1]) < EPS);
    assert(fabs(fxy - h_sq*eik->jets[l].fxy) < EPS);
  }
}
#endif

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

  int lc = grid2_l2lc(eik->grid, l) + eik->tri_dlc[ic0];
  if (lc < 0 || eik->ncells <= lc) {
    return;
  }

  bicubic_s *bicubic = &eik->bicubics[lc];
  if (!bicubic_valid(bicubic)) {
    return;
  }

#if SJS_DEBUG
  check_cell_consistency(eik, lc);
#endif

  /**
   * Get cubic along edge of interest.
   */
  bicubic_variable var = get_tri_bicubic_var(eik, ic0);
  int edge = get_tri_edge(eik, ic0);
  cubic_s T_cubic = bicubic_get_f_on_edge(bicubic, var, edge);
  cubic_s Tx_cubic = bicubic_get_fx_on_edge(bicubic, var, edge);
  cubic_s Ty_cubic = bicubic_get_fy_on_edge(bicubic, var, edge);
  if (get_should_reverse_cubic(eik, ic0)) {
    cubic_reverse_on_unit_interval(&T_cubic);
    cubic_reverse_on_unit_interval(&Tx_cubic);
    cubic_reverse_on_unit_interval(&Ty_cubic);
  }

  dbl h = eik->grid->h;

  assert(fabs(cubic_f(&T_cubic, 0) - eik->jets[l0].f) < EPS);
  assert(fabs(cubic_f(&T_cubic, 1) - eik->jets[l1].f) < EPS);
  assert(fabs(cubic_f(&Tx_cubic, 0)/h - eik->jets[l0].Df[0]) < EPS);
  assert(fabs(cubic_f(&Tx_cubic, 1)/h - eik->jets[l1].Df[0]) < EPS);
  assert(fabs(cubic_f(&Ty_cubic, 0)/h - eik->jets[l0].Df[1]) < EPS);
  assert(fabs(cubic_f(&Ty_cubic, 1)/h - eik->jets[l1].Df[1]) < EPS);

  dbl2 xy; grid2_l2xy(eik->grid, l, xy);
  dbl2 xy0; grid2_l2xy(eik->grid, l0, xy0);
  dbl2 xy1; grid2_l2xy(eik->grid, l1, xy1);

  if (var == LAMBDA) {
    assert(fabs(xy0[1] - xy1[1]) < EPS);
    assert(fabs(h - fabs(xy0[0] - xy1[0])) < EPS);
    assert(fabs(h - fabs(xy[1] - xy0[1])) < EPS);
    assert(fabs(h - fabs(xy[1] - xy0[1])) < EPS);
  } else if (var == MU) {
    assert(fabs(xy0[0] - xy1[0]) < EPS);
    assert(fabs(h - fabs(xy0[1] - xy1[1])) < EPS);
    assert(fabs(h - fabs(xy[0] - xy0[0])) < EPS);
    assert(fabs(h - fabs(xy[0] - xy1[0])) < EPS);
  }

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

    // If we failed to find an interior point minimizer, then find the
    // endpoint with the smaller value and use that as the minimizer.
    if (!found) {
      F3_compute(0, &context);
      dbl F3_eta0 = context.F3;
      F3_compute(1, &context);
      dbl F3_eta1 = context.F3;
      eta = F3_eta0 < F3_eta1 ? 0 : 1;
    }
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
    F4_context context = {
      .T_cubic = T_cubic,
      .Tx_cubic = Tx_cubic,
      .Ty_cubic = Ty_cubic,
      .slow = eik->slow
    };
    dbl2_copy(xy, context.xy);
    dbl2_copy(xy0, context.xy0);
    dbl2_copy(xy1, context.xy1);

    dbl2 xk, gk, xk1, gk1;
    dbl22 Hk, Hk1;
    F4_bfgs_init(eta, th, xk, gk, Hk, &context);

    dbl Tprev = context.F4;

    int iter = 0;
    while (dbl2_maxnorm(gk) > EPS &&
           F4_bfgs_step(xk, gk, Hk, xk1, gk1, Hk1, &context)) {
      if (xk1[0] < -EPS || xk1[0] > 1 + EPS) {
        printf("out of bounds: eta = %g\n", xk1[0]);
        abort();
      }

      xk1[0] = fmax(0, fmin(1, xk1[0]));

      T = context.F4;

      ++iter;

      if (fabs(T - Tprev) <= 1e1*EPS*(fabs(fmax(T, Tprev)) + 1)) {
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
        printf("no decrease in T: (%g > %g, num accepted: %lu)\n", T, Tprev,
               eik->num_accepted);
        abort();
      }

      Tprev = T;

      dbl2_copy(xk1, xk);
      dbl2_copy(gk1, gk);
      dbl22_copy(Hk1, Hk);
    }

    eta = xk[0];
    th = xk[1];
  }

  //////////////////////////////////////////////////////////////////////////////

  // Check causality
  assert(T > eik->jets[l0].f);
  assert(T > eik->jets[l1].f);

  /**
   * Commit new value if it's an improvement.
   */
  jet21p *jet = &eik->jets[l];
  if (T < jet->f) {
    jet->f = T;

    dbl s = field2_f(eik->slow, xy);
    jet->Df[0] = s*cos(th);
    jet->Df[1] = s*sin(th);

    // Update node's parent
    eik->par[l] = (par2_s) {.l = {l0, l1}, .b = {1 - eta, eta}};
  }
}

static bool can_build_cell(eik_s const *eik, int lc) {
  int2 ind;
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = grid2_lc2l(eik->grid, lc) + eik->vert_dl[i];
    grid2_l2ind(eik->grid, l, ind);
    if (!grid2_isind(eik->grid, ind) ||
        eik->states[l] != VALID ||
        jet21p_is_point_source(&eik->jets[l]))
      return false;
  }

  /* Once we know that we can build this cell, *then* we check that
   * all of its data is finite. If we try to check for finite jets in
   * the loop above, we'll run into trouble. A cell can have nodes
   * which are VALID and aren't point sources, but still haven't
   * received a value for Txy. All Txy values for a cell will be
   * finite *only if* all of the nodes are VALID. */
#if SJS_DEBUG
  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = grid2_lc2l(eik->grid, lc) + eik->vert_dl[i];
    assert(jet21p_is_finite(&eik->jets[l]));
  }
#endif

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
  if (!grid2_isindc(eik->grid, indc)) {
    return false;
  }
  int l = grid2_indc2l(eik->grid, indc);
  int2 offset, indv;
  for (int iv = 0, lv; iv < NUM_CELL_VERTS; ++iv) {
    get_cell_vert_offset(eik, iv, offset);
    int2_add(indc, offset, indv);
    if (!grid2_isind(eik->grid, indv)) {
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
  //   fx.data[iv] = eik->jets[lv].Df[0];
  //   fy.data[iv] = eik->jets[lv].Df[1];
  // }
  // return interpolate_fxy_at_verts(fx, fy, eik->grid->h);

  /**
   * ... we'll just use this for now since I have more confidence that
   * it *does* work
   */

  dbl fx[NUM_CELL_VERTS], fy[NUM_CELL_VERTS];

  for (int i = 0, l; i < NUM_CELL_VERTS; ++i) {
    l = grid2_lc2l(eik->grid, lc) + eik->vert_dl[i];
    fx[i] = eik->jets[l].Df[0];
    fy[i] = eik->jets[l].Df[1];
  }

  dbl fxy[NUM_CELL_VERTS] = {
    (fy[1] - fy[0])/eik->grid->h, // left
    (fx[3] - fx[1])/eik->grid->h, // bottom
    (fx[2] - fx[0])/eik->grid->h, // top
    (fy[3] - fy[2])/eik->grid->h // right
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
    l = grid2_lc2l(eik->grid, lc) + eik->vert_dl[i];
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
    l[i] = grid2_lc2l(eik->grid, lc) + eik->vert_dl[i];
  }

  /* Get jet at each cell vertex */
  jet21p *J[4];
  for (int i = 0; i < NUM_CELL_VERTS; ++i) {
    J[i] = &eik->jets[l[i]];
  }

  /* Precompute scaling factors for partial derivatives */
  dbl h = eik->grid->h, h_sq = h*h;

  /* Compute cell data from jets and scaling factors */
  dbl44 data;
  int2 offset;
  for (int i = 0, di, dj; i < NUM_CELL_VERTS; ++i) {
    get_cell_vert_offset(eik, i, offset);
    di = offset[0];
    dj = offset[1];
    data[di][dj] = J[i]->f;
    data[2 + di][dj] = h*J[i]->Df[0];
    data[di][2 + dj] = h*J[i]->Df[1];
    data[2 + di][2 + dj] = h_sq*J[i]->fxy;
  }

  /* Set cell data */
  bicubic_set_data(&eik->bicubics[lc], data);
}

/* Compute whether each neighboring node is inbounds. */
static void set_inbounds(eik_s const *eik, int l, bool inbounds[NUM_NB + 1]) {
  int2 ind, ind0;
  grid2_l2ind(eik->grid, l, ind);
  for (int i0 = 0; i0 < NUM_NB + 1; ++i0) {
    int2_add(ind, offsets[i0], ind0);
    inbounds[i0] = grid2_isind(eik->grid, ind0);
  }
}

static void update(eik_s *eik, int l) {
  bool inbounds[NUM_NB + 1];
  set_inbounds(eik, l, inbounds);

  for (int i0 = 1, l0, l1, ic0; i0 < 8; i0 += 2) {
    if (!inbounds[i0]) {
      continue;
    }

    l0 = l + eik->nb_dl[i0];
    if (eik->states[l0] != VALID) {
      continue;
    }

    if (inbounds[i0 - 1]) {
      l1 = l + eik->nb_dl[i0 - 1];
      if (eik->states[l1] == VALID) {
        ic0 = i0 - 1;
        tri(eik, l, l0, l1, ic0);
      }
    }

    if (inbounds[i0 + 1]) {
      l1 = l + eik->nb_dl[i0 + 1];
      if (eik->states[l1] == VALID) {
        ic0 = i0;
        tri(eik, l, l0, l1, ic0);
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
void eik_init(eik_s *eik, field2_s const *slow, grid2_s const *grid) {
  eik->slow = slow;
  eik->grid = grid;
  eik->ncells = grid2_nindc(grid);
  eik->nnodes = grid2_nind(grid);
  eik->bicubics = malloc(eik->ncells*sizeof(bicubic_s));
  eik->jets = malloc(eik->nnodes*sizeof(jet21p));
  eik->states = malloc(eik->nnodes*sizeof(state_e));
  eik->positions = malloc(eik->nnodes*sizeof(int));

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

  int capacity = (int) 3*sqrt(eik->nnodes);
  heap_init(eik->heap, capacity, value, setpos, (void *)eik);

  eik->num_accepted = 0;
  eik->accepted = malloc(eik->nnodes*sizeof(size_t));
  for (int l = 0; l < eik->nnodes; ++l)
    eik->accepted[l] = (size_t)NO_INDEX;

  eik->par = malloc(eik->nnodes*sizeof(par2_s));
  for (int l = 0; l < eik->nnodes; ++l)
    par2_init_empty(&eik->par[l]);

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
    eik->jets[l].Df[0] = NAN;
    eik->jets[l].Df[1] = NAN;
    eik->jets[l].fxy = NAN;
  }

  for (int l = 0; l < eik->nnodes; ++l) {
    eik->states[l] = FAR;
  }
}

void eik_deinit(eik_s *eik) {
  printf("eik_deinit()\n");

  eik->slow = NULL;

  free(eik->bicubics);
  free(eik->jets);
  free(eik->states);
  free(eik->positions);
  free(eik->accepted);
  free(eik->par);

  eik->bicubics = NULL;
  eik->jets = NULL;
  eik->states = NULL;
  eik->positions = NULL;
  eik->accepted = NULL;
  eik->par = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);
}

size_t eik_peek(eik_s const *eik) {
  return heap_front(eik->heap);
}

static bool recompute_nearby_cells(eik_s *eik, size_t l0) {
  int2 ind0;
  grid2_l2ind(eik->grid, l0, ind0);

  bool nb_inc_on_valid_cell[9];
  // Determine which of the cells surrounding l0 are now valid. It's
  // enough to check if any of the four nearest cells are valid: it's
  // impossible for them to have been valid (or built before), since
  // one of their vertices just became valid.
  bool valid_cell[NUM_NB_CELLS];
  for (int ic = 0; ic < NUM_NB_CELLS; ++ic) {
    int2 indc, offset;
    get_nb_cell_offset(eik, ic, offset);
    int2_add(ind0, offset, indc);
    valid_cell[ic] = cell_is_valid(eik, indc);
  }

  if (!valid_cell[0] && !valid_cell[1] && !valid_cell[2] && !valid_cell[3])
    return false;

  bool *inc = &nb_inc_on_valid_cell[0];
  if (eik->grid->order == ORDER_ROW_MAJOR) {
    inc[0] = valid_cell[0];
    inc[1] = valid_cell[0] || valid_cell[1];
    inc[2] =                  valid_cell[1];
    inc[3] = valid_cell[0]                  || valid_cell[2];
    inc[4] = valid_cell[0] || valid_cell[1] || valid_cell[2] || valid_cell[3];
    inc[5] =                  valid_cell[1]                  || valid_cell[3];
    inc[6] =                                   valid_cell[2];
    inc[7] =                                   valid_cell[2] || valid_cell[3];
    inc[8] =                                                    valid_cell[3];
  } else if (eik->grid->order == ORDER_COLUMN_MAJOR) {
    inc[0] = valid_cell[0];
    inc[1] = valid_cell[0]                  || valid_cell[2];
    inc[2] =                                   valid_cell[2];
    inc[3] = valid_cell[0] || valid_cell[1];
    inc[4] = valid_cell[0] || valid_cell[1] || valid_cell[2] || valid_cell[3];
    inc[5] =                                   valid_cell[2] || valid_cell[3];
    inc[6] =                  valid_cell[1];
    inc[7] =                  valid_cell[1]                  || valid_cell[3];
    inc[8] =                                                    valid_cell[3];
  } else {
    abort();
  }

  // Array used to indicated which nearby cells should be used to
  // compute average values for Txy at cell vertices.
  bool use_for_Txy_average[NUM_NEARBY_CELLS];
  memset(use_for_Txy_average, 0x0, sizeof(bool)*NUM_NEARBY_CELLS);

  // Traverse the nearby cells and figure out which should be used to
  // compute new Txy values.
  //
  // NOTE: `use_for_Txy_average` => the cell is inbounds
  for (int ic = 0; ic < NUM_NEARBY_CELLS; ++ic) {
    for (int jv = 0, i; jv < NUM_CELL_VERTS; ++jv) {
      i = cell_vert_to_cell_nb_vert(eik, ic, jv);
      if (i == NO_INDEX) {
        continue;
      }
      use_for_Txy_average[ic] |= nb_inc_on_valid_cell[i];
    }

    int2 offset, indc;
    get_nearby_cell_offset(eik, ic, offset);
    int2_add(ind0, offset, indc);
    use_for_Txy_average[ic] &= cell_is_valid(eik, indc);
  }

  // Compute new Txy values at the vertices of the cells that we
  // decided to use.
  dbl4 Txy[NUM_NEARBY_CELLS];
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    if (use_for_Txy_average[ic]) {
      lc = grid2_l2lc(eik->grid, l0) + eik->nearby_dlc[ic];
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
        jc = cell_nb_vert_to_nearby_cell(eik, i, j);
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
      lc = grid2_l2lc(eik->grid, l0) + eik->nearby_dlc[ic];
      build_cell(eik, lc);
    }
  }

  return true;
}

#if SJS_DEBUG
static void check_consistency_of_nearby_cells(eik_s const *eik, int l0) {
  for (int ic = 0, lc; ic < NUM_NEARBY_CELLS; ++ic) {
    lc = grid2_l2lc(eik->grid, l0) + eik->nearby_dlc[ic];
    if (can_build_cell(eik, lc))
      check_cell_consistency(eik, lc);
  }
}
#endif

void eik_step(eik_s *eik) {
  int l0 = heap_front(eik->heap);
  assert(eik->states[l0] == TRIAL);
  heap_pop(eik->heap);
  eik->states[l0] = VALID;

  eik->accepted[eik->num_accepted++] = l0;

  bool updated_cells = recompute_nearby_cells(eik, l0);
#if SJS_DEBUG
  if (updated_cells)
    check_consistency_of_nearby_cells(eik, l0);
#endif

  int2 ind0;
  grid2_l2ind(eik->grid, l0, ind0);

  // Set FAR nodes to TRIAL and insert them into the heap.
  for (int i = 0, l; i < NUM_NB; ++i) {
    int2 ind; int2_add(ind0, offsets[i], ind);
    if (!grid2_isind(eik->grid, ind)) {
      continue;
    }
    l = l0 + eik->nb_dl[i];
    if (eik->states[l] == FAR) {
      eik->states[l] = TRIAL;
      heap_insert(eik->heap, l);
    }
  }

  // If we didn't update any neighboring cells, we won't do any
  // updates, so return now.
  if (!updated_cells)
    return;

  // Update neighboring nodes.
  for (int i = 0, l; i < NUM_NB; ++i) {
    int2 ind; int2_add(ind0, offsets[i], ind);
    if (!grid2_isind(eik->grid, ind)) {
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

void eik_add_trial(eik_s *eik, int2 ind, jet21p jet) {
  int l = grid2_ind2l(eik->grid, ind);
  eik->jets[l] = jet;
  assert(eik->states[l] != TRIAL && eik->states[l] != VALID);
  eik->states[l] = TRIAL;
  heap_insert(eik->heap, l);
}

void eik_add_valid(eik_s *eik, int2 ind, jet21p jet) {
  int l = grid2_ind2l(eik->grid, ind);
  eik->jets[l] = jet;
  assert(eik->states[l] != TRIAL && eik->states[l] != VALID);
  eik->states[l] = VALID;
  eik->accepted[eik->num_accepted++] = l;
}

void eik_make_bd(eik_s *eik, int2 ind) {
  int l = grid2_ind2l(eik->grid, ind);
  eik->states[l] = BOUNDARY;
}

void eik_get_shape(eik_s const *eik, int2 shape) {
  int2_copy(eik->grid->shape, shape);
}

jet21p eik_get_jet(eik_s *eik, int2 ind) {
  int l = grid2_ind2l(eik->grid, ind);
  return eik->jets[l];
}

jet21p *eik_get_jets_ptr(eik_s const *eik) {
  return eik->jets;
}

state_e eik_get_state(eik_s const *eik, int2 ind) {
  int l = grid2_ind2l(eik->grid, ind);
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

dbl eik_T(eik_s const *eik, dbl2 xy) {
  dbl2 cc;
  int lc = grid2_xy2lc(eik->grid, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_f(bicubic, cc);
}

dbl eik_Tx(eik_s const *eik, dbl2 xy) {
  dbl2 cc;
  int lc = grid2_xy2lc(eik->grid, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_fx(bicubic, cc)/eik->grid->h;
}

dbl eik_Ty(eik_s const *eik, dbl2 xy) {
  dbl2 cc;
  int lc = grid2_xy2lc(eik->grid, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_fy(bicubic, cc)/eik->grid->h;
}

dbl eik_Txx(eik_s const *eik, dbl2 xy) {
  dbl2 cc;
  int lc = grid2_xy2lc(eik->grid, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_fxx(bicubic, cc)/(eik->grid->h*eik->grid->h);
}

dbl eik_Txy(eik_s const *eik, dbl2 xy) {
  dbl2 cc;
  int lc = grid2_xy2lc(eik->grid, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_fxy(bicubic, cc)/(eik->grid->h*eik->grid->h);
}

dbl eik_Tyy(eik_s const *eik, dbl2 xy) {
  dbl2 cc;
  int lc = grid2_xy2lc(eik->grid, xy, cc);
  if (!can_build_cell(eik, lc)) {
    return NAN;
  }
  bicubic_s *bicubic = &eik->bicubics[lc];
  return bicubic_fyy(bicubic, cc)/(eik->grid->h*eik->grid->h);
}

bool eik_can_build_cell(eik_s const *eik, int2 indc) {
  int lc = grid2_indc2lc(eik->grid, indc);
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
  int lc = grid2_indc2lc(eik->grid, indc);
  return eik->bicubics[lc];
}

bicubic_s *eik_get_bicubics_ptr(eik_s const *eik) {
  return eik->bicubics;
}

heap_s *eik_get_heap(eik_s const *eik) {
  return eik->heap;
}

par2_s eik_get_par(eik_s const *eik, int2 ind) {
  int l = grid2_ind2l(eik->grid, ind);
  return eik->par[l];
}

bool eik_has_par(eik_s const *eik, int2 ind) {
  int l = grid2_ind2l(eik->grid, ind);
  return !par2_is_empty(&eik->par[l]);
}

size_t const *eik_get_accepted_ptr(eik_s const *eik) {
  return eik->accepted;
}
