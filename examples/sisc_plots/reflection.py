import colorcet as cc
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

from jmm.eik import Eik
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
    shape = np.array([M, N], dtype=np.intc)
    xymin = np.array([xmin, ymin], dtype=np.float64)
    grid = Grid2(shape, xymin, h)
    rfac = 0.1

    x, y = np.meshgrid(np.linspace(xmin, xmax, M),
                       np.linspace(ymin, ymax, N),
                       indexing='ij')

    s, eik = get_eik_for_pt_src(grid, c, x0, y0, rfac)
    eik.solve()


    shape_refl = np.array([M, N + 2], dtype=np.intc) # add an extra column
    grid_refl = Grid2(shape_refl, xymin, h)

    eik_refl = Eik.from_s_and_grid(s, grid_refl)

    # add row of reflected valid jets in second ghost layer
    for i in range(M):
        ind = np.array([i, N+1], dtype=np.intc)
        jet = Jet2(eik.T[i, N-3], eik.Tx[i, N-3], -eik.Ty[i, N-3], eik.Txy[i, N-3])
        eik_refl.add_valid(ind, jet)

    # add row of reflected valid jets in first ghost layer
    for i in range(M):
        ind = np.array([i, N], dtype=np.intc)
        jet = Jet2(eik.T[i, N-2], eik.Tx[i, N-2], -eik.Ty[i, N-2], eik.Txy[i, N-2])
        eik_refl.add_valid(ind, jet)

    # add row of reflected trial jets on boundary
    for i in range(M):
        ind = np.array([i, N-1], dtype=np.intc)
        jet = Jet2(eik.T[i, N-1], eik.Tx[i, N-1], -eik.Ty[i, N-1], eik.Txy[i, N-1])
        eik_refl.add_trial(ind, jet)

    eik_refl.build_cells()

    eik_refl.solve()

    Tmin, Tmax = 0, eik_refl.T.max()
    levels = np.linspace(Tmin, Tmax, 19)

    plt.figure()
    cs = plt.contour(x, y, eik.T, levels, linewidths=1.5, linestyles='--',
                     zorder=1)
    plt.gca().clabel(cs, [cs.levels[7]], inline=True,
                     fmt={cs.levels[7]: r'$T_{in}$'}, fontsize=14, zorder=3)
    # plt.quiver(x[:, -1], y[:, -1], eik.Tx[:, -1], eik.Ty[:, -1])
    cs = plt.contour(x, y, eik_refl.T[:, :-2], levels, linewidths=1.5,
                     linestyles='-', zorder=2)
    plt.gca().clabel(cs, [cs.levels[-9]], inline=True,
                     fmt={cs.levels[-9]: r'$T_{refl}$'}, fontsize=14, zorder=4)
    # plt.quiver(x[:, -1], y[:, -1], eik_refl.Tx[:, -1], eik_refl.Ty[:, -1])
    norm = matplotlib.colors.Normalize(vmin=Tmin, vmax=Tmax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cs.cmap)
    sm.set_array([])
    plt.colorbar(sm, ticks=levels[::3])
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.show()
