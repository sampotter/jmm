#include "eik.h"
#include "npy.h"

#include <stdlib.h>
#include <stdio.h>

#define MAX(x, y) x > y ? x : y
#define VX 0.133
#define VY -0.0933
#define VNORM sqrt(VX*VX + VY*VY);
#define VNORMSQ (VX*VX + VY*VY)

dbl vdot(dbl x, dbl y) {
  return VX*x + VY*y;
}

dbl normsq(dbl x, dbl y) {
  return x*x + y*y;
}

dbl s(dbl x, dbl y, void *context) {
  (void) context;
  return 1.0/(1.0 + vdot(x, y));
}

dbl sx(dbl x, dbl y) {
  return -VX*pow(s(x, y, NULL), 2.0);
}

dbl sy(dbl x, dbl y) {
  return -VY*pow(s(x, y, NULL), 2.0);
}

dbl sxy(dbl x, dbl y) {
  return 2*VX*VY*pow(s(x, y, NULL), 3.0);
}

dvec2 grad_s(dbl x, dbl y, void *context) {
  (void) context;
  return (dvec2) {.x = sx(x, y), .y = sy(x, y)};
}

// Note: below, `u` is the solution of |grad(tau)| = s, with s defined
// as above. We write `u` in terms of an auxiliary function we define
// below called `f`. This makes it simpler to write down its partial
// derivatives.

dbl f(dbl x, dbl y) {
  return 1 + s(x, y, NULL)*VNORMSQ*normsq(x, y)/2;
}

dbl fx(dbl x, dbl y) {
  return VNORMSQ*(sx(x, y)*normsq(x, y) + 2*s(x, y, NULL)*x)/2;
}

dbl fy(dbl x, dbl y) {
  return VNORMSQ*(sy(x, y)*normsq(x, y) + 2*s(x, y, NULL)*y)/2;
}

dbl fxy(dbl x, dbl y) {
  return VNORMSQ*(sxy(x, y)*normsq(x, y) + 2*(sx(x, y)*y + sy(x, y)*x))/2;
}

dbl dacosh(dbl z) {
  return pow(z - 1, -0.5)*pow(z + 1, -0.5);
}

dbl d2acosh(dbl z) {
  return -z*pow(z - 1, -1.5)*pow(z + 1, -1.5);
}

dbl u(dbl x, dbl y) {
  return acosh(f(x, y))/VNORM;
}

dbl ux(dbl x, dbl y) {
  return dacosh(f(x, y))*fx(x, y)/VNORM;
}

dbl uy(dbl x, dbl y) {
  return dacosh(f(x, y))*fy(x, y)/VNORM;
}

dbl uxy(dbl x, dbl y) {
  dbl tmp = f(x, y);
  return (d2acosh(tmp)*fx(x, y)*fy(x, y) + dacosh(tmp)*fxy(x, y))/VNORM;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("usage: %s <N>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  eik * scheme;
  eik_alloc(&scheme);

  field2_s slow = {
    .f = s,
    .grad_f = grad_s,
    .context = NULL
  };

  int N = atoi(argv[1]);
  int i0 = N/2;
  ivec2 shape = {N, N};
  dvec2 xymin = {-1, -1};
  dbl h = 2.0/(N-1);
  eik_init(scheme, &slow, shape, xymin, h);

  int R = N/20;
  if (R < 5) R = 5;

  /**
   * Initialize inside box.
   */
  // for (int i = 0; i < N; ++i) {
  //   int abs_di = abs(i - i0);
  //   dbl x = h*i + xymin.x;
  //   for (int j = 0; j < N; ++j) {
  //     int abs_dj = abs(j - i0);
  //     dbl y = h*j + xymin.y;
  //     int r = MAX(abs_di, abs_dj);
  //     if (r <= R) {
  //       ivec2 ind = {i, j};
  //       jet J = {u(x, y), ux(x, y), uy(x, y), uxy(x, y)};
  //       if (r < R) {
  //         eik_add_valid(scheme, ind, J);
  //       } else {
  //         eik_add_trial(scheme, ind, J);
  //       }
  //     }
  //   }
  // }

  /**
   * Initialize inside disk.
   */
  for (int i = 0; i < N; ++i) {
    int di = i - i0, di_sq = di*di;
    for (int j = 0; j < N; ++j) {
      int dj = j - i0, dj_sq = dj*dj;
      dbl r = sqrt(di_sq + dj_sq);
      if (r < R) {
        dbl x = h*i + xymin.x;
        dbl y = h*j + xymin.y;
        jet J = {u(x, y), ux(x, y), uy(x, y), uxy(x, y)};
        eik_add_valid(scheme, (ivec2) {i, j}, J);
      }
    }
  }
  {
    dbl di[4] = {1, 0, -1,  0};
    dbl dj[4] = {0, 1,  0, -1};
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        ivec2 ind = {i, j};
        if (eik_get_state(scheme, ind) != VALID) {
          continue;
        }
        for (int k = 0; k < 4; ++k) {
          int i_ = i + di[k];
          int j_ = j + dj[k];
          if (0 <= i_ && i_ < N && 0 <= j_ && j_ < N) {
            ivec2 ind_ = {i_, j_};
            state_e state = eik_get_state(scheme, ind_);
            if (state != VALID && state != TRIAL) {
              dbl x = h*i_ + xymin.x;
              dbl y = h*j_ + xymin.y;
              jet J = {u(x, y), ux(x, y), uy(x, y), uxy(x, y)};
              eik_add_trial(scheme, ind_, J);
            }
          }
        }
      }
    }
  }

  eik_build_cells(scheme);

  eik_solve(scheme);

  jet_s *jets = eik_get_jets_ptr(scheme);

  npy_write_2d_dbl_array("T.npy", &jets[0].f, N, N, sizeof(jet_s));
  npy_write_2d_dbl_array("Tx.npy", &jets[0].fx, N, N, sizeof(jet_s));
  npy_write_2d_dbl_array("Ty.npy", &jets[0].fy, N, N, sizeof(jet_s));
  npy_write_2d_dbl_array("Txy.npy", &jets[0].fxy, N, N, sizeof(jet_s));

  eik_deinit(scheme);
  eik_dealloc(&scheme);
}
