#include "eik3.h"

#include <assert.h>
#include <stdlib.h>

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

typedef struct costfunc {
  dbl f;
  dvec2 g;
  dmat22 H;

  dvec3 x; // x[l]
  dmat33 X; // X = [x[l0]'; x[l1]'; x[l2]']
  dmat33 Xt;

  // B-coefs for 9-point triangle interpolation T on base of update
  dbl Tc[10];

  dvec3 x_minus_xb;
} costfunc_s;

void costfunc_init(costfunc_s *cf, eik3_s *eik, size_t l,
                   size_t l0, size_t l1, size_t l2) {
  mesh3_get_vert(eik->mesh, l, cf->x.data);

  mesh3_get_vert(eik->mesh, l0, cf->Xt.rows[0].data);
  mesh3_get_vert(eik->mesh, l1, cf->Xt.rows[1].data);
  mesh3_get_vert(eik->mesh, l2, cf->Xt.rows[2].data);

  jet3 jet;
  dbl T[3];
  dvec3 DT[3];

  jet = eik->jet[l0];
  T[0] = jet.f;
  DT[0] = (dvec3) {.data = {jet.fx, jet.fy, jet.fz}};

  jet = eik->jet[l1];
  T[1] = jet.f;
  DT[1] = (dvec3) {.data = {jet.fx, jet.fy, jet.fz}};

  jet = eik->jet[l2];
  T[2] = jet.f;
  DT[2] = (dvec3) {.data = {jet.fx, jet.fy, jet.fz}};

  bb3tri_interp3(T, DT, cf->Xt.rows, cf->Tc);
}

void costfunc_set_lambda(costfunc_s *cf, dbl const *lambda) {
  static dvec3 const a1 = {.data = {-1, 1, 0}};
  static dvec3 const a2 = {.data = {-1, 0, 1}};

  dvec3 b;
  b.data[1] = lambda[0];
  b.data[2] = lambda[1];
  b.data[0] = 1 - b.data[1] - b.data[2];

  dvec3 xb = dmat33_dvec3_mul(cf->X, b);
  cf->x_minus_xb = dvec3_sub(cf->x, xb);
  dbl L = dvec3_norm(cf->x_minus_xb);

  dvec3 tmp1 = dmat33_dvec3_mul(cf->X, cf->x_minus_xb);
  tmp1 = dvec3_dbl_div(tmp1, L);

  dvec2 DL;
  DL.x = dvec3_dot(a1, tmp1);
  DL.y = dvec3_dot(a2, tmp1);

  dmat33 X = cf->Xt;
  dmat33_transpose(&X);

  dmat33 tmp2 = dmat33_mul(X, cf->Xt);
  tmp2 = dmat33_sub(tmp2, dvec3_outer(tmp1, tmp1));
  tmp2 = dmat33_dbl_div(tmp2, L);

  dmat22 D2L;
  D2L.data[0][0] = dvec3_dot(dmat33_dvec3_mul(tmp2, a1), a1);
  D2L.data[1][0] = D2L.data[0][1] = dvec3_dot(dmat33_dvec3_mul(tmp2, a1), a2);
  D2L.data[1][1] = dvec3_dot(dmat33_dvec3_mul(tmp2, a2), a2);

  dvec2 DT;
  DT.x = dbb3tri(cf->Tc, b.data, a1.data);
  DT.y = dbb3tri(cf->Tc, b.data, a2.data);

  dmat22 D2T;
  D2T.data[0][0] = d2bb3tri(cf->Tc, b.data, a1.data, a1.data);
  D2T.data[1][0] = D2T.data[0][1] = d2bb3tri(cf->Tc, b.data, a1.data, a2.data);
  D2T.data[1][1] = d2bb3tri(cf->Tc, b.data, a1.data, a2.data);

  cf->f = L + bb3tri(cf->Tc, b.data);
  cf->g = dvec2_add(DL, DT);
  cf->H = dmat22_add(D2L, D2T);
}

static dbl get_residual(costfunc_s const *cf, dbl lam0, dbl const *lam,
                        dbl const *alp, dbl t) {
  dvec2 tmp1 = cf->g;
  tmp1.x += alp[0] - alp[1];
  tmp1.y += alp[0] - alp[2];

  dbl tmp2 = 1/(t*alp[0]) - lam0;

  dvec2 tmp3;
  tmp3.x = 1/(t*alp[1]) - lam[0];
  tmp3.x = 1/(t*alp[2]) - lam[1];

  return sqrt(dvec2_norm_sq(tmp1) + tmp2*tmp2 + dvec2_norm_sq(tmp3));
}

static jet3 tetra(eik3_s *eik, size_t l, size_t l0, size_t l1, size_t l2) {
  static int const maxniter = 100;
  static dbl const eps = 1e-15;
  static dbl const beta = 0.5;
  static dbl const gamma = 1e-2;

  costfunc_s cf;
  costfunc_init(&cf, eik, l, l0, l1, l2);

  dvec2 lam = {.x = 1./3, .y = 1./3}, dlam;
  costfunc_set_lambda(&cf, &lam.x);

  dbl lam0 = 1./3;

  // Vector of Lagrange multipliers
  dvec3 alp = {.data = {3, 3, 3}}, dalp;

  static dbl const mu = 1e2;

  dbl t = mu, s, res, newres;
  int k = 1;
  dbl eta = 1;

  res = get_residual(&cf, lam0, &lam.x, alp.data, t);

  // primal-dual iteration
  while (res > eps || eta > 1 || k < maxniter) {
    /**
     * compute Newton step for modified KKT system
     */

    dlam.x = alp.data[1] - alp.data[0];
    dlam.y = alp.data[2] - alp.data[0];
    dlam = dvec2_sub(dlam, cf.g);

    dalp.data[0] = lam0 - 1/(t*alp.data[0]);
    dalp.data[1] = lam.x - 1/(t*alp.data[1]);
    dalp.data[2] = lam.y - 1/(t*alp.data[2]);

    dlam.x += lam0*dalp.data[0]/alp.data[0] - lam.x*dalp.data[1]/alp.data[1];
    dlam.y += lam0*dalp.data[0]/alp.data[0] - lam.y*dalp.data[2]/alp.data[2];

    dmat22 H = cf.H;
    dbl tmp = lam.x/alp.data[0];
    H.data[0][0] += tmp + lam.x/alp.data[1];
    H.data[1][0] = H.data[0][1] = tmp;
    H.data[1][1] += tmp + lam.y/alp.data[2];

    dlam = dmat22_dvec2_solve(H, dlam);
    dbl dlam0 = dlam.x + dlam.y;

    dalp.data[0] += lam0*dlam0/alp.data[0];
    dalp.data[1] -= lam.x*dlam.x/alp.data[1];
    dalp.data[2] -= lam.y*dlam.y/alp.data[2];

    /**
     * do backtracking line search
     */

    s = 1;
    s = fmin(s, dalp.data[0] < 0 ? -alp.data[0]/dalp.data[0] : 1);
    s = fmin(s, dalp.data[1] < 0 ? -alp.data[1]/dalp.data[1] : 1);
    s = fmin(s, dalp.data[2] < 0 ? -alp.data[2]/dalp.data[2] : 1);

    s = fmin(s, dlam0 < 0 ? lam0/dlam0 : 1);

    s = fmin(s, dlam.x < 0 ? -lam.x/dlam.x : 1);
    s = fmin(s, dlam.y < 0 ? -lam.y/dlam.y : 1);

    while (true) {
      lam = dvec2_saxpy(s, dlam, lam);
      alp = dvec3_saxpy(s, dalp, alp);
      costfunc_set_lambda(&cf, &lam.x);
      newres = get_residual(&cf, lam0, &lam.x, alp.data, t);
      if (res < (1 - gamma*s)*newres) {
        break;
      }
      s *= beta;
    }

    // prepare for next iteration
    res = newres;
    lam0 = 1 - lam.x - lam.y;
    eta = lam0*alp.data[0] + lam.x*alp.data[1] + lam.y*alp.data[2];
    t = 3*mu/eta;
    k += 1;
  }

  assert(k < maxniter);

  dvec3 DT = dvec3_normalized(cf.x_minus_xb);

  return (jet3) {
    .f = cf.f, .fx = DT.data[0], .fy = DT.data[1], .fz = DT.data[2]};
}

static void update(eik3_s *eik, size_t l, size_t l0) {
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

  // Don't need ec any more
  free(ec);

  // Do the updates
  jet3 jet;
  for (int i = 0; i < nup; ++i) {
    jet = tetra(eik, l, l0, l1[i], l2[i]);
    if (jet.f < eik->jet[l].f) {
      eik->jet[l] = jet;
    }
  }

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
