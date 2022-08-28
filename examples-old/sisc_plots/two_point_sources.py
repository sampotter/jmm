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
    h = 1/32
    xmin, ymin = 0, 0
    xmax, ymax = 1, 1
    M = int(np.round(xmax/h)) + 1
    N = int(np.round(ymax/h)) + 1
    shape = np.array([M, N], dtype=np.intc)
    xymin = np.array([xmin, ymin], dtype=np.float64)
    grid = Grid2(shape, xymin, h)
    rfac = 0.1

    s0 = LinearSpeedFunc2(c0, x0, y0)
    s1 = LinearSpeedFunc2(c1, x1, y1)

    eik = Eik.from_s_and_grid(s0, grid)
    add_exact_data_for_pt_src(eik, grid, s0, x0, y0, rfac)
    add_exact_data_for_pt_src(eik, grid, s1, x1, y1, rfac)
    eik.build_cells()
    eik.solve()

    X, Y = np.meshgrid(np.linspace(xmin, xmax, M),
                       np.linspace(ymin, ymax, N),
                       indexing='ij')
    tau0 = np.array([
        s0.tau(*_) for _ in zip(X.ravel(), Y.ravel())]).reshape(X.shape)
    tau1 = np.array([
        s1.tau(*_) for _ in zip(X.ravel(), Y.ravel())]).reshape(X.shape)

    argmin = np.argmin(np.array([tau0, tau1]), axis=0)

    tau = np.minimum(tau0, tau1)

    taux = np.array([
        [s0.taux(x, y), s1.taux(x, y)][i]
        for x, y, i in zip(X.ravel(), Y.ravel(), argmin.ravel())
    ]).reshape(X.shape)

    tauy = np.array([
        [s0.tauy(x, y), s1.tauy(x, y)][i]
        for x, y, i in zip(X.ravel(), Y.ravel(), argmin.ravel())
    ]).reshape(X.shape)

    tauxy = np.array([
        [s0.tauxy(x, y), s1.tauxy(x, y)][i]
        for x, y, i in zip(X.ravel(), Y.ravel(), argmin.ravel())
    ]).reshape(X.shape)

    # T = np.array([eik.get_T(np.array(_)) for _
    #               in zip(X.ravel(), Y.ravel())]).reshape(X.shape)

    # Tx = np.array([eik.get_Tx(np.array(_)) for _
    #               in zip(X.ravel(), Y.ravel())]).reshape(X.shape)

    # Ty = np.array([eik.get_Ty(np.array(_)) for _
    #               in zip(X.ravel(), Y.ravel())]).reshape(X.shape)

    T = eik.T
    Tx = eik.Tx
    Ty = eik.Ty

    extent = [xmin, xmax, ymin, ymax]

    plt.figure(figsize=(11, 6))

    # column #1: tau vs T

    plt.subplot(2, 3, 1)
    plt.imshow(np.rot90(T), extent=extent, cmap=cc.cm.rainbow)
    plt.colorbar()
    plt.contour(X, Y, T, linewidths=1, linestyles='--', colors='k')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title(r'$T$')

    E = tau - T
    vmin, vmax = -np.nanmax(abs(E)), np.nanmax(abs(E))
    plt.subplot(2, 3, 4)
    plt.imshow(np.rot90(E), extent=extent, vmin=vmin, vmax=vmax,
               cmap=cc.cm.coolwarm)
    cbar = plt.colorbar()
    cbar.formatter.set_powerlimits((0, 0))
    plt.contour(X, Y, T, linewidths=1, linestyles='--', colors='k')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title(r'$|\tau - T|$')

    # column #2: tau_x vs T_x

    plt.subplot(2, 3, 2)
    plt.imshow(np.rot90(Tx), extent=extent, cmap=cc.cm.rainbow)
    plt.colorbar()
    plt.contour(X, Y, T, linewidths=1, linestyles='--', colors='k')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title(r'$T_x$')

    E = taux - Tx
    vmin, vmax = -np.nanmax(abs(E)), np.nanmax(abs(E))
    plt.subplot(2, 3, 5)
    plt.imshow(np.rot90(E), extent=extent, vmin=vmin, vmax=vmax,
               cmap=cc.cm.coolwarm)
    cbar = plt.colorbar()
    cbar.formatter.set_powerlimits((0, 0))
    plt.contour(X, Y, T, linewidths=1, linestyles='--', colors='k')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title(r'$|\tau_x - T_x|$')

    # column #3: tau_y vs T_y

    plt.subplot(2, 3, 3)
    plt.imshow(np.rot90(Ty), extent=extent, cmap=cc.cm.rainbow)
    plt.colorbar()
    plt.contour(X, Y, T, linewidths=1, linestyles='--', colors='k')
    plt.contour(X, Y, T, linewidths=1, linestyles='--', colors='k')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title(r'$T_y$')

    E = tauy - Ty
    vmin, vmax = -np.nanmax(abs(E)), np.nanmax(abs(E))
    plt.subplot(2, 3, 6)
    plt.imshow(np.rot90(E), extent=extent, vmin=vmin, vmax=vmax,
               cmap=cc.cm.coolwarm)
    cbar = plt.colorbar()
    cbar.formatter.set_powerlimits((0, 0))
    plt.contour(X, Y, T, linewidths=1, linestyles='--', colors='k')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title(r'$|\tau_y - T_y|$')

    plt.tight_layout()

    plt.savefig('two-point-sources.pdf')

    plt.show()
