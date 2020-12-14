#include "eik3.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bb.h"
#include "heap.h"
#include "mat.h"
#include "mesh3.h"
#include "vec.h"

struct eik3 {
  mesh3_s const *mesh;
  jet3 *jet;
  state_e *state;
  int *pos;
  heap_s *heap;
};

void eik3_alloc(eik3_s **eik) {
  *eik = malloc(sizeof(eik3_s));
}

void eik3_dealloc(eik3_s **eik) {
  free(*eik);
  *eik = NULL;
}

static dbl value(void *ptr, int l) {
  eik3_s *eik = (eik3_s *)ptr;
  assert(l >= 0);
  assert(l < (int)mesh3_nverts(eik->mesh));
  dbl T = eik->jet[l].f;
  return T;
}

static void setpos(void *ptr, int l, int pos) {
  eik3_s *eik = (eik3_s *)ptr;
  eik->pos[l] = pos;
}

void eik3_init(eik3_s *eik, mesh3_s const *mesh) {
  eik->mesh = mesh;

  size_t nverts = mesh3_nverts(mesh);

  eik->jet = malloc(nverts*sizeof(jet3));
  for (size_t l = 0; l < nverts; ++l) {
    eik->jet[l] = (jet3) {.f = INFINITY, .fx = NAN, .fy = NAN, .fz = NAN};
  }

  eik->state = malloc(nverts*sizeof(state_e));
  for (size_t l = 0; l < nverts; ++l) {
    eik->state[l] = FAR;
  }

  eik->pos = malloc(nverts*sizeof(int));
  for (size_t l = 0; l < nverts; ++l) {
    eik->pos[l] = NO_INDEX;
  }

  /**
   * When we compute the initial heap capacity, we want to estimate
   * the number of nodes that could comprise the expanding numerical
   * front at any one time. We can't know this ahead of time, so we
   * set it to a constant multiple times (# nodes)^(1/d). In this
   * case, d=3. Even if this is an underestimate, well still reduce
   * the number of times the heap needs to be expanded at solve time.
   */
  int capacity = (int) 3*cbrt(nverts);

  heap_alloc(&eik->heap);
  heap_init(eik->heap, capacity, value, setpos, (void *)eik);
}

void eik3_deinit(eik3_s *eik) {
  free(eik->jet);
  eik->jet = NULL;

  free(eik->state);
  eik->state = NULL;

  free(eik->pos);
  eik->pos = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);
}

void costfunc_init(costfunc_s *cf, mesh3_s const *mesh, jet3 const *jet,
                   size_t l, size_t l0, size_t l1, size_t l2) {
  assert(jet3_is_finite(&jet[l0]));
  assert(jet3_is_finite(&jet[l1]));
  assert(jet3_is_finite(&jet[l2]));

  mesh3_get_vert(mesh, l, cf->x);

  mesh3_get_vert(mesh, l0, cf->Xt[0]);
  mesh3_get_vert(mesh, l1, cf->Xt[1]);
  mesh3_get_vert(mesh, l2, cf->Xt[2]);

  dbl33_transposed(cf->Xt, cf->X);
  dbl33_mul(cf->X, cf->Xt, cf->XXt);

  /**
   * Compute Bernstein-Bezier coefficients before transposing Xt and computing XXt
   */

  dbl T[3];
  dbl DT[3][3];

  T[0] = jet[l0].f;
  DT[0][0] = jet[l0].fx;
  DT[0][1] = jet[l0].fy;
  DT[0][2] = jet[l0].fz;

  T[1] = jet[l1].f;
  DT[1][0] = jet[l1].fx;
  DT[1][1] = jet[l1].fy;
  DT[1][2] = jet[l1].fz;

  T[2] = jet[l2].f;
  DT[2][0] = jet[l2].fx;
  DT[2][1] = jet[l2].fy;
  DT[2][2] = jet[l2].fz;

  bb3tri_interp3(T, &DT[0], cf->Xt, cf->Tc);
}

void costfunc_set_lambda(costfunc_s *cf, dbl const *lambda) {
  static dbl a1[3] = {-1, 1, 0};
  static dbl a2[3] = {-1, 0, 1};

  dbl b[3], xb[3], tmp1[3], tmp2[3][3], L, DL[2], D2L[2][2], DT[2], D2T[2][2];

  b[1] = lambda[0];
  b[2] = lambda[1];
  b[0] = 1 - b[1] - b[2];

  assert(b[0] >= 0);
  assert(b[1] >= 0);
  assert(b[2] >= 0);

  dbl33_dbl3_mul(cf->X, b, xb);
  dbl3_sub(cf->x, xb, cf->x_minus_xb);
  L = dbl3_norm(cf->x_minus_xb);

  dbl33_dbl3_mul(cf->Xt, cf->x_minus_xb, tmp1);
  dbl3_dbl_div(tmp1, L, tmp1);

  DL[0] = dbl3_dot(a1, tmp1);
  DL[1] = dbl3_dot(a2, tmp1);

  dbl3_outer(tmp1, tmp1, tmp2);
  dbl33_sub(cf->XXt, tmp2, tmp2);
  dbl33_dbl_div(tmp2, L, tmp2);

  dbl33_dbl3_mul(tmp2, a1, tmp1);
  D2L[0][0] = dbl3_dot(tmp1, a1);
  D2L[1][0] = D2L[0][1] = dbl3_dot(tmp1, a2);
  dbl33_dbl3_mul(tmp2, a2, tmp1);
  D2L[1][1] = dbl3_dot(tmp1, a2);

  DT[0] = dbb3tri(cf->Tc, b, a1);
  DT[1] = dbb3tri(cf->Tc, b, a2);

  D2T[0][0] = d2bb3tri(cf->Tc, b, a1, a1);
  D2T[1][0] = D2T[0][1] = d2bb3tri(cf->Tc, b, a1, a2);
  D2T[1][1] = d2bb3tri(cf->Tc, b, a2, a2);

  cf->f = L + bb3tri(cf->Tc, b);
  dbl2_add(DL, DT, cf->g);
  dbl22_add(D2L, D2T, cf->H);

  /**
   * Finally, compute Newton step, making sure to perturb the Hessian
   * if it's indefinite.
   */

  // Conditionally perturb the Hessian
  dbl tr = cf->H[0][0] + cf->H[1][1];
  dbl det = cf->H[0][0]*cf->H[1][1] - cf->H[0][1]*cf->H[1][0];
  dbl min_eig_doubled = tr - sqrt(tr*tr - 4*det);
  if (min_eig_doubled < 0) {
    cf->H[0][0] -= min_eig_doubled;
    cf->H[1][1] -= min_eig_doubled;
  }

  // Compute the Newton step
  dbl22_dbl2_solve(cf->H, cf->g, cf->p);
  dbl2_negate(cf->p);
}

void tetra(costfunc_s *cf, dbl const lam[2], jet3 *jet) {
  dbl const tscale = 0.5; // Step size scaling parameter
  dbl const c1 = 1e-4; // Constant for backtracking line search
  dbl const ftol = 1e-15, xtol = 1e-15;

  dbl lam1[2], dlam[2];
  dbl f = cf->f; // Current value of cost function
  dbl t = 1; // Initial step size
  dbl tc; // Breakpoint used to find Cauchy point
  dbl c1_times_g_dot_p;
  dbl Df; // Directional derivative used in Cauchy point calculation

  cf->niter = 0;

  /**
   * Newton iteration
   */
  while (true) {
    c1_times_g_dot_p = c1*dbl2_dot(cf->g, cf->p);
    assert(c1_times_g_dot_p < 0);

  cauchy:

    tc = (1 - lam[0] - lam[1])/dbl2_sum(cf->p);
    if (0 < tc && tc < 1) {
      dbl2_saxpy(tc, cf->p, lam, lam1);
      costfunc_set_lambda(cf, lam1);
      dbl2_sub(lam1, lam, dlam);
      Df = dbl2_dot(cf->g, dlam);
      if (Df >= 0 || cf->f >= f) {
        t = tc;
        goto backtrack;
      } else {
        cf->p[0] = (cf->p[0] - cf->p[1])/2;
        cf->p[1] = (cf->p[1] - cf->p[0])/2;
        lam = lam1;
        goto cauchy;
      }
    }

    tc = -lam[0]/cf->p[0];
    if (0 < tc && tc < 1) {
      dbl2_saxpy(tc, cf->p, lam, lam1);
      costfunc_set_lambda(cf, lam1);
      dbl2_sub(lam1, lam, dlam);
      Df = dbl2_dot(cf->g, dlam);
      if (Df >= 0 || cf->f >= f) {
        t = tc;
        goto backtrack;
      } else {
        cf->p[0] = 0;
        lam = lam1;
        goto cauchy;
      }
    }

    tc = -lam[1]/cf->p[1];
    if (0 < tc && tc < 1) {
      dbl2_saxpy(tc, cf->p, lam, lam1);
      costfunc_set_lambda(cf, lam1);
      dbl2_sub(lam1, lam, dlam);
      Df = dbl2_dot(cf->g, dlam);
      if (Df >= 0 || cf->f >= f) {
        t = tc;
        goto backtrack;
      } else {
        cf->p[1] = 0;
        lam = lam1;
        goto cauchy;
      }
    }

    // We didn't trip any breakpoints, compute lam1
    dbl2_saxpy(t, cf->p, lam, lam1);
    costfunc_set_lambda(cf, lam1);

  backtrack:
    while (cf->f > f + t*c1_times_g_dot_p) {
      t *= tscale;
      dbl2_saxpy(t, cf->p, lam, lam1);
      costfunc_set_lambda(cf, lam1);
    }

    // Check for convergence
    if (fabs(cf->f - f) <= ftol*fmax(cf->f, f) + ftol ||
        dbl2_maxdist(lam, lam1)
        <= xtol*fmax(dbl2_maxnorm(lam), dbl2_maxnorm(lam1)) + xtol) {
      break;
    }

    // Reset for next iteration
    t = 1;
    lam = lam1;
    f = cf->f;
    ++cf->niter;
  }

  dbl DT[3];
  memcpy((void *)DT, (void *)cf->x_minus_xb, sizeof(dbl)*3);
  dbl3_normalize(DT);

  jet->f = f;
  jet->fx = DT[0];
  jet->fy = DT[1];
  jet->fz = DT[2];
}

static void update(eik3_s *eik, size_t l, size_t l0) {
  /**
   * Update setup phase
   */

  // To find the updates incident on l and l0, first find the
  // tetrahedra that are incident on the edge (l0, l).
  int nec = mesh3_nec(eik->mesh, l0, l);
  size_t *ec = malloc(nec*sizeof(size_t));
  mesh3_ec(eik->mesh, l0, l, ec);

  // Allocate space for the l1 and l2 indices.
  size_t *l1 = malloc(nec*sizeof(size_t));
  size_t *l2 = malloc(nec*sizeof(size_t));

  int nup = 0;
  for (int i = 0; i < nec; ++i) {
    // Get the indices of cell ec[i]'s vertices
    size_t lv[4];
    mesh3_cv(eik->mesh, ec[i], lv);

    // Find the VALID vertices that *aren't* l0 or l
    int lnew[2];
    int k = 0;
    for (int j = 0; j < 4; ++j) {
      if (lv[j] != l0 && lv[j] != l && eik->state[lv[j]] == VALID) {
        lnew[k++] = lv[j];
      }
    }
    assert(k <= 2);

    // If we found a pair of VALID verts, add a new update, assigning
    // l1 and l2
    if (k == 2) {
      l1[nup] = lnew[0];
      l2[nup++] = lnew[1];
    }
  }

  /**
   * Update evaluation phase
   */

  costfunc_s cf;
  dbl lambda[2] = {0, 0};

  // Start by searching for an update tetrahedron that might have an
  // interior point solution
  for (int i = 0; i < nup; ++i) {
    costfunc_init(&cf, eik->mesh, eik->jet, l, l0, l1[i], l2[i]);
    costfunc_set_lambda(&cf, lambda);
    if (cf.g[0] > 0 || cf.g[0] > 0) {
      continue;
    } else if (dbl2_maxnorm(cf.g) <= EPS) {
      // Minimizing ray passes through x0
      eik->jet[l] = eik->jet[l0];
      eik->jet[l].f = cf.f;
      break;
    } else {
      tetra(&cf, lambda, &eik->jet[l]);
    }
  }

  /**
   * Cleanup
   */
  free(ec);
  free(l1);
  free(l2);
}

static void adjust(eik3_s *eik, size_t l) {
  assert(eik->state[l] == TRIAL);
  assert(l >= 0);
  assert(l < mesh3_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l]);
}

void eik3_step(eik3_s *eik) {
  size_t l, l0 = heap_front(eik->heap);
  assert(eik->state[l0] == TRIAL);
  heap_pop(eik->heap);
  eik->state[l0] = VALID;

  // Get i0's neighboring nodes.
  int nnb = mesh3_nvv(eik->mesh, l0);
  size_t *nb = malloc(nnb*sizeof(size_t));
  mesh3_vv(eik->mesh, l0, nb);

  // Set FAR nodes to TRIAL and insert them into the heap.
  for (int i = 0; i < nnb; ++i) {
    if (eik->state[l = nb[i]] == FAR) {
      eik->state[l] = TRIAL;
      heap_insert(eik->heap, l);
    }
  }

  // Update neighboring nodes.
  for (int i = 0; i < nnb; ++i) {
    if (eik->state[l = nb[i]] == TRIAL) {
      update(eik, l, l0);
      adjust(eik, l);
    }
  }

  free(nb);
}

void eik3_solve(eik3_s *eik) {
  while (heap_size(eik->heap) > 0) {
    eik3_step(eik);
  }
}

void eik3_add_trial(eik3_s *eik, size_t l, jet3 jet) {
  if (eik->state[l] == VALID) {
    return;
  } else if (eik->pos[l] == NO_INDEX) {
    assert(eik->state[l] == FAR);
    eik->jet[l] = jet;
    eik->state[l] = TRIAL;
    heap_insert(eik->heap, l);
  } else if (jet.f < eik->jet[l].f) {
    assert(eik->state[l] == TRIAL);
    eik->jet[l] = jet;
    adjust(eik, l);
  }
}

void eik3_add_valid(eik3_s *eik, size_t l, jet3 jet) {
  // TODO: need to decide what to do if the state is TRIAL
  if (eik->state[l] == TRIAL) {
    abort();
  }

  eik->jet[l] = jet;
  eik->state[l] = VALID;
}

bool eik3_is_valid(eik3_s const *eik, size_t l) {
  return eik->state[l] == VALID;
}

jet3 *eik3_get_jet_ptr(eik3_s const *eik) {
  return eik->jet;
}
