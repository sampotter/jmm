import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np

from jmm.grid import Grid2

from util import *

if __name__ == '__main__':
    plt.ion()

    # Speed function
    x0, y0 = 0, 0
    x1, y1 = 0.8, 0
    c0 = np.array([2, 5, 20], dtype=np.float64)
    c1 = c0 + np.array([c0[1]*x1 + c0[2]*y1, 0, 0])

    # Discretization parameters
    h = 0.001
    xmin, ymin = 0, 0
    xmax, ymax = 1, 1
    M = int(np.round(xmax/h)) + 1
    N = int(np.round(ymax/h)) + 1
    shape = np.array([M, N], dtype=np.intc)
    xymin = np.array([xmin, ymin], dtype=np.float64)
    grid = Grid2(shape, xymin, h)
    rfac = 0.1

    s0, eik0 = get_eik_for_pt_src(grid, c0, x0, y0, rfac)
    eik0.solve()

    s1, eik1 = get_eik_for_pt_src(grid, c1, x1, y1, rfac)
    eik1.solve()

    T = np.minimum(eik0.T, eik1.T)

    x, y = np.meshgrid(np.linspace(xmin, xmax, M),
                       np.linspace(ymin, ymax, N),
                       indexing='ij')
    tau0 = np.array([s0.tau(*_) for _ in zip(x.ravel(), y.ravel())]).reshape(x.shape)
    tau1 = np.array([s1.tau(*_) for _ in zip(x.ravel(), y.ravel())]).reshape(x.shape)
    tau = np.minimum(tau0, tau1)

    extent = [xmin, xmax, ymin, ymax]

    E = abs(tau - np.minimum(eik0.T, eik1.T))
    plt.figure()
    plt.imshow(np.rot90(E), extent=extent, cmap=cc.cm.fire)
    plt.colorbar()
    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.imshow(np.rot90(T), extent=extent, cmap=cc.cm.rainbow)
    plt.colorbar()
    plt.contour(x, y, T, linewidths=1, linestyles='--', colors='k')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.tight_layout()
    plt.show()
