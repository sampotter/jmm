import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np

from jmm.grid import Grid2
from jmm.jet import Jet2

from util import *

if __name__ == '__main__':
    plt.ion()

    # Speed function
    x0, y0 = 0, 0
    c = np.array([2.0, 7.0, 5.0])

    # Discretization parameters
    h = 0.01
    xmin, ymin = 0, 0
    xmax, ymax = 1, 1
    M = int(np.round(xmax/h)) + 1
    N = int(np.round(ymax/h)) + 1
    shape = np.array([M, N], dtype=np.intc) # add an extra column
    xymin = np.array([xmin, ymin], dtype=np.float64)
    grid = Grid2(shape, xymin, h)
    rfac = 0.1

    x, y = np.meshgrid(np.linspace(xmin, xmax, M),
                       np.linspace(ymin, ymax, N),
                       indexing='ij')

    s, eik = get_eik_for_pt_src(grid, c, x0, y0, rfac)
    eik.solve()

    plt.figure()
    plt.contour(x, y, eik.T)
    plt.gca().set_aspect('equal')
    plt.show()

    # # FIXME: this doesn't work
    # eik_refl = Eik.from_s_and_grid(s, grid)
    # j = N - 1
    # for i in range(M):
    #     ind = np.array([i, j], dtype=np.intc)
    #     jet = Jet2(eik.T[i, j], eik.Tx[i, j], -eik.Ty[i, j], eik.Txy[i, j])
    #     eik_refl.add_trial(ind, jet)
    # eik_refl.build_cells()
    # eik_refl.solve()
