import colorcet as cc
import itertools as it
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

from jmm.eik import Eik
from jmm.grid import Grid2
from jmm.jet import Jet2

from util import *

def compute_J(eik, s, grid, X, Y, R=None, J=None):
    v0 = 1/s.s0
    v = s.v
    initial_J = J is not None
    if not initial_J:
        assert R is not None
        J = np.empty(eik.shape, dtype=np.float64)
        J[...] = np.nan
    for l in eik.accepted:
        ind = grid.l2ind(l)
        xhat = np.array([X[ind], Y[ind]])
        par = eik.get_par(np.array(ind, dtype=np.intc))
        if par is None:
            if not initial_J:
                J[ind] = np.linalg.norm(xhat)
            else:
                assert np.isfinite(J[ind])
        else:
            l0, l1 = par.l
            ind0, ind1 = grid.l2ind(l0), grid.l2ind(l1)
            if not np.isfinite(J[ind0]):
                import pdb; pdb.set_trace()
            if not np.isfinite(J[ind1]):
                import pdb; pdb.set_trace()
            x0, x1 = np.array([X[ind0], Y[ind0]]), np.array([X[ind1], Y[ind1]])
            b0, b1 = par.b
            xlam = b0*x0 + b1*x1
            L = np.linalg.norm(xhat - xlam)
            eps = L*(v0 + v@(xhat + xlam)/2)
            DeltaTlam = eik.get_Txx(xlam) + eik.get_Tyy(xlam)
            tlam = np.array([eik.get_Tx(xlam), eik.get_Ty(xlam)])
            tlam /= np.linalg.norm(tlam)
            gradslam = np.array([s.sx(*xlam), s.sy(*xlam)])
            Jlam = b0*J[ind0] + b1*J[ind1]
            J[ind] = abs(1 + eps*(DeltaTlam - tlam@gradslam))*Jlam
    return J

if __name__ == '__main__':
    plt.ion()

    omega = 1000

    # Speed function
    x0, y0 = 0, 0
    c = np.array([2, 5, 7], dtype=np.float64)

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

    C = np.array([1/s.s(*_) for _ in zip(x.ravel(), y.ravel())]).reshape(x.shape)

    J = compute_J(eik, s, grid, x, y, R=np.sqrt((x - x0)**2 + (y - y0)**2))
    A = np.exp(1j*np.pi/4)*np.sqrt(C/J)/(2*np.sqrt(2*np.pi*omega))

    shape_refl = np.array([M, N + 2], dtype=np.intc) # add an extra column
    grid_refl = Grid2(shape_refl, xymin, h)
    X_refl, Y_refl = np.meshgrid(np.linspace(xmin, xmax, M),
                                 np.linspace(ymin, ymax + 2*h, N + 2),
                                 indexing='ij')

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
    T_refl = eik_refl.T[:, :-2]

    J_refl = np.empty(shape_refl)
    J_refl[...] = np.nan
    J_refl[:, N+1] = J[:, N-3]
    J_refl[:, N] = J[:, N-2]
    J_refl[:, N-1] = J[:, N-1]
    J_refl = compute_J(eik_refl, s, grid_refl, X_refl, Y_refl, J=J_refl)
    A_refl = np.exp(1j*np.pi/4)*np.sqrt(C/J_refl[:, :-2])/(2*np.sqrt(2*np.pi*omega))

    W = A*np.exp(-1j*omega*eik.T)
    W_refl = A_refl*np.exp(-1j*omega*T_refl)
    U = W + W_refl

    Tmin, Tmax = 0, eik_refl.T.max()
    levels = np.linspace(Tmin, Tmax, 19)

    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    cs = plt.contour(x, y, eik.T, levels, linewidths=1.5, linestyles='--',
                     zorder=1, cmap=cc.cm.rainbow)
    plt.gca().clabel(cs, [cs.levels[7]], inline=True,
                     fmt={cs.levels[7]: r'$T_{in}$'}, fontsize=14, zorder=3)
    cs = plt.contour(x, y, T_refl, levels, linewidths=1.5,
                     linestyles='-', zorder=2, cmap=cc.cm.rainbow)
    plt.gca().clabel(cs, [cs.levels[-9]], inline=True,
                     fmt={cs.levels[-9]: r'$T_{refl}$'}, fontsize=14, zorder=4)
    norm = matplotlib.colors.Normalize(vmin=Tmin, vmax=Tmax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cs.cmap)
    sm.set_array([])
    plt.colorbar(sm, ticks=levels[::3])
    plt.gca().set_aspect('equal')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title('$T$ [s]')
    plt.subplot(1, 2, 2)
    U_real = np.real(U)
    vmax = abs(U_real).max()
    vmin = -vmax
    plt.imshow(np.rot90(U_real), extent=[xmin, xmax, ymin, ymax],
               vmin=vmin, vmax=vmax, cmap=cc.cm.rainbow)
    plt.colorbar()
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title(r'$\Re(U)$')
    plt.tight_layout()
    plt.savefig('reflection.pdf')
    plt.show()
