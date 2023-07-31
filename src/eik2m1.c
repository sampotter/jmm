#include <jmm/eik2m1.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <jmm/heap.h>
#include <jmm/log.h>
#include <jmm/mat.h>
#include <jmm/utri21.h>

struct eik2m1 {
  mesh22_s const *mesh;
  jet22t *jet;
  state_e *state;
  size_t *pos;
  heap_s *heap;
  par2_s *par;
  size_t nvalid;
};

static dbl value(eik2m1_s const *eik, int l) {
  assert(l >= 0);
  assert(l < (int)mesh22_nverts(eik->mesh));
  return eik->jet[l].f;
}

static void setpos(eik2m1_s const *eik, int l, int pos) {
  eik->pos[l] = pos;
}

void eik2m1_alloc(eik2m1_s **eik) {
  *eik = malloc(sizeof(eik2m1_s));
}

void eik2m1_dealloc(eik2m1_s **eik) {
  free(*eik);
}

void eik2m1_init(eik2m1_s *eik, mesh22_s const *mesh) {
  eik->mesh = mesh;

  size_t nverts = mesh22_nverts(eik->mesh);

  eik->jet = malloc(nverts*sizeof(jet22t));
  for (size_t l = 0; l < nverts; ++l)
    eik->jet[l] = jet22t_make_empty();

  eik->state = malloc(nverts*sizeof(jet22t));
  for (size_t l = 0; l < nverts; ++l)
    eik->state[l] = FAR;

  eik->pos = malloc(nverts*sizeof(jet22t));
  for (size_t l = 0; l < nverts; ++l)
    eik->pos[l] = NO_INDEX;

  heap_alloc(&eik->heap);
  heap_init(eik->heap, 3*sqrt(nverts), (value_f)value, (setpos_f)setpos, eik);

  eik->par = malloc(nverts*sizeof(par2_s));
  for (size_t l = 0; l < nverts; ++l)
    par2_init_empty(&eik->par[l]);

  eik->nvalid = 0;
}

void eik2m1_deinit(eik2m1_s *eik) {
  eik->mesh = NULL;

  free(eik->jet);
  eik->jet = NULL;

  free(eik->state);
  eik->state = NULL;

  free(eik->pos);
  eik->pos = NULL;

  heap_deinit(eik->heap);
  heap_dealloc(&eik->heap);

  free(eik->par);
  eik->par = NULL;
}

size_t eik2m1_peek(eik2m1_s const *eik) {
  return heap_front(eik->heap);
}

static void tri(eik2m1_s *eik, size_t l, size_t l0, size_t l1) {
  mesh22_s const *mesh = eik->mesh;

  dbl2 xhat, x[2];
  mesh22_get_vert(mesh, l, xhat);
  mesh22_get_vert(mesh, l0, x[0]);
  mesh22_get_vert(mesh, l1, x[1]);

  jet22t jet[2] = {eik->jet[l0], eik->jet[l1]};

  utri21_s utri;
  utri21_init(&utri, xhat, x, jet);

  dbl lam;
  if (!utri21_solve(&utri, &lam))
    return;

  // check if tangent vector of ray lies in cone spanned by DT0 and DT1

  dbl22 DT0_and_DT1;
  dbl2_copy(eik->jet[l0].Df, DT0_and_DT1[0]);
  dbl2_copy(eik->jet[l1].Df, DT0_and_DT1[1]);
  dbl22_transpose(DT0_and_DT1);

  dbl2 cone_coefs;
  dbl22_dbl2_solve(DT0_and_DT1, utri.jet.Df, cone_coefs);
  if (cone_coefs[0] < 0 || cone_coefs[1] < 0)
    return;

  // update jet
  eik->jet[l] = utri.jet;

  // set parent
  eik->par[l] = (par2_s) {.l = {l0, l1}, .b = {1 - lam, lam}};
}

static void get_op_edge(size_t const lf[3], size_t l, size_t le[2]) {
  size_t j = 0;
  for (size_t i = 0; i < 3; ++i)
    if (lf[i] != l)
      le[j++] = lf[i];
}

// static void update(eik2m1_s *eik, size_t l) {
//   log_debug("update(l = %lu)", l);

//   size_t nvf = mesh22_nvf(eik->mesh, l);
//   size_t *vf = malloc(nvf*sizeof(size_t));
//   mesh22_vf(eik->mesh, l, vf);

//   size_t lv[3], le[2];

//   for (size_t i = 0; i < nvf; ++i) {
//     mesh22_fv(eik->mesh, vf[i], lv);
//     get_op_edge(lv, l, le);
//     log_debug("l = %lu, le[0] = %lu, le[1] = %lu", l, le[0], le[1]);
//     if (eik->state[le[0]] == VALID && eik->state[le[1]] == VALID)
//       tri(eik, l, le[0], le[1]);
//   }

//   free(vf);
// }

static void adjust(eik2m1_s *eik, size_t l0) {
  assert(eik->state[l0] == TRIAL);
  assert(l0 < mesh22_nverts(eik->mesh));

  heap_swim(eik->heap, eik->pos[l0]);
}

size_t eik2m1_step(eik2m1_s *eik) {
  size_t l0 = heap_front(eik->heap);
  assert(eik->state[l0] == TRIAL);
  heap_pop(eik->heap);
  eik->state[l0] = VALID;

  if (!jet22t_is_finite(&eik->jet[l0])) {
    log_error("found non-finite jet!");

    dbl d_min = INFINITY;
    size_t d_argmin = (size_t)NO_INDEX;
    for (size_t l = 0; l < mesh22_nverts(eik->mesh); ++l) {
      dbl2 x;
      mesh22_get_vert(eik->mesh, l, x);
      if (eik->state[l] == TRIAL) {
        dbl d = dbl2_norm(x);
        if (d < d_min) {
          d_min = d;
          d_argmin = l;
        }
      }
    }
    printf("min d (= %g) obtained by l = %lu\n", d_min, d_argmin);

    // return NO_INDEX;
    abort();
  }

  size_t nvv = mesh22_nvv(eik->mesh, l0);
  size_t* vv = malloc(nvv*sizeof(size_t));
  mesh22_vv(eik->mesh, l0, vv);

  for (size_t i = 0; i < nvv; ++i) {
    size_t l = vv[i];
    if (eik->state[l] == FAR) {
      eik->state[l] = TRIAL;
      heap_insert(eik->heap, l);
    }
  }

  // Return early if all of l0's neighbors are VALID:

  size_t ntrial = 0;
  for (size_t i = 0; i < nvv; ++i)
    ntrial += eik->state[vv[i]] == TRIAL;
  if (ntrial == 0)
    goto cleanup;

  // First, find the edges on the VALID front that are incident on l0:

  size_t nvf = mesh22_nvf(eik->mesh, l0);
  size_t *vf = malloc(nvf*sizeof(size_t));
  mesh22_vf(eik->mesh, l0, vf);

  size_t lv[3], le[2], l1[2] = {NO_INDEX, NO_INDEX}, nvalid;
  for (size_t i = 0; i < nvf; ++i) {
    if (l1[0] != (size_t)NO_INDEX && l1[1] != (size_t)NO_INDEX)
      break;

    mesh22_fv(eik->mesh, vf[i], lv);
    get_op_edge(lv, l0, le);

    nvalid = (eik->state[le[0]] == VALID) + (eik->state[le[1]] == VALID);

    if (nvalid == 0 || nvalid == 2)
      continue;

    size_t l_valid = eik->state[le[0]] == VALID ? le[0] : le[1];

    if (l1[0] == (size_t)NO_INDEX)
      l1[0] = l_valid;
    if (l1[0] != l_valid && l1[1] == (size_t)NO_INDEX)
      l1[1] = l_valid;
  }

  free(vf);

  assert(l1[0] != (size_t)NO_INDEX || l1[1] != (size_t)NO_INDEX);

  // Update each TRIAL node l neighboring l0 from the edges on the
  // VALID front:

  for (size_t i = 0; i < nvv; ++i) {
    size_t l = vv[i];

    if (eik->state[l] != TRIAL)
      continue;

    if (l1[0] != (size_t)NO_INDEX)
      tri(eik, l, l0, l1[0]);

    if (l1[1] != (size_t)NO_INDEX)
      tri(eik, l, l0, l1[1]);

    adjust(eik, l);
  }

cleanup:
  free(vv);

  return l0;
}

void eik2m1_solve(eik2m1_s *eik) {
  while (heap_size(eik->heap) > 0) {
    // size_t l0 = eik2m1_step(eik);
    eik2m1_step(eik);
    ++eik->nvalid;
  }
}

void eik2m1_add_trial(eik2m1_s *eik, size_t l, jet22t jet) {
  eik->jet[l] = jet;
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = TRIAL;
  heap_insert(eik->heap, l);
}

void eik2m1_add_valid(eik2m1_s *eik, size_t l, jet22t jet) {
  eik->jet[l] = jet;
  assert(eik->state[l] != TRIAL && eik->state[l] != VALID);
  eik->state[l] = VALID;
  ++eik->nvalid;
}

bool eik2m1_is_valid(eik2m1_s const *eik, size_t l) {
  return eik->state[l] == VALID;
}

bool eik2m1_is_trial(eik2m1_s const *eik, size_t l) {
  return eik->state[l] == TRIAL;
}

bool eik2m1_is_far(eik2m1_s const *eik, size_t l) {
  return eik->state[l] == FAR;
}

state_e const *eik2m1_get_state_ptr(eik2m1_s const *eik) {
  return eik->state;
}

jet22t const *eik2m1_get_jet_ptr(eik2m1_s const *eik) {
  return eik->jet;
}

par2_s eik2m1_get_par(eik2m1_s const *eik, size_t l) {
  return eik->par[l];
}

par2_s const *eik2m1_get_par_ptr(eik2m1_s const *eik) {
  return eik->par;
}
