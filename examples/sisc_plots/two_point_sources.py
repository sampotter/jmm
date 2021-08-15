import itertools as it
import matplotlib.pyplot as plt
import numpy as np

from jmm.eik import Eik
from jmm.field import LinearSpeedFunc2
from jmm.grid import Grid2

def get_trial_mask(valid_mask):
    shape = valid_mask.shape
    trial_mask = np.zeros_like(valid_mask)
    for i, j in it.product(*(range(dim) for dim in shape)):
        if valid_mask[i, j]:
            continue
        for di, dj in it.product([-1, 0, 1], repeat=2):
            if 0 <= i + di < shape[0] and \
               0 <= j + dj < shape[1] and \
               valid_mask[i + di, j + dj]:
                trial_mask[i, j] = True
    return trial_mask

if __name__ == '__main__':
    # Speed function
    c = np.array([0.5, 5, 20])

    # Discretization parameters
    h = 0.01
    xmin, ymin = 0, 0
    xmax, ymax = 1, 0.5
    M = int(np.round(xmax/h)) + 1
    N = int(np.round(ymax/h)) + 1
    shape = np.array([M, N], dtype=np.intc)
    xymin = np.array([xmin, ymin], dtype=np.float64)
    grid = Grid2(shape, xymin, h)
    rfac = 0.1
    x, y = np.meshgrid(np.linspace(xmin, xmax, M),
                       np.linspace(ymin, ymax, N),
                       indexing='ij')

    def get_eik_for_pt_src(x0, y0):
        s = LinearSpeedFunc2(c, x0, y0)
        eik = Eik.from_s_and_grid(s, grid)
        valid_mask = (x0 - x)**2 + (y0 - y)**2 < rfac**2
        for i, j in zip(*np.where(valid_mask)):
            ind = np.array([i, j], dtype=np.intc)
            eik.add_valid(ind, s.get_jet(x[i, j], y[i, j]))
        for i, j in zip(*np.where(get_trial_mask(valid_mask))):
            ind = np.array([i, j], dtype=np.intc)
            eik.add_trial(ind, s.get_jet(x[i, j], y[i, j]))
        eik.build_cells()
        return eik

    print('eik0')
    eik0 = get_eik_for_pt_src(0, 0)
    eik0.solve()

    print('eik1')
    eik1 = get_eik_for_pt_src(0.8, 0)
    eik1.solve()

    ########################################################################
    # debugging

    s = LinearSpeedFunc2(c, 0, 0)

    T = np.empty_like(x)
    for i, j in it.product(range(65), repeat=2):
        T[i, j] = s.tau(x[i, j], y[i, j])

    Tx = np.empty_like(x)
    for i, j in it.product(range(65), repeat=2):
        Tx[i, j] = s.taux(x[i, j], y[i, j])

    Ty = np.empty_like(x)
    for i, j in it.product(range(65), repeat=2):
        Ty[i, j] = s.tauy(x[i, j], y[i, j])

    Txx = np.empty_like(x)
    for i, j in it.product(range(65), repeat=2):
        Txx[i, j] = s.tauxx(x[i, j], y[i, j])

    Txy = np.empty_like(x)
    for i, j in it.product(range(65), repeat=2):
        Txy[i, j] = s.tauxy(x[i, j], y[i, j])

    Tyy = np.empty_like(x)
    for i, j in it.product(range(65), repeat=2):
        Tyy[i, j] = s.tauyy(x[i, j], y[i, j])

    Tx_fd, Ty_fd = np.gradient(T, h, h, edge_order=2)
    Txx_fd, Txy_fd = np.gradient(Tx, h, h, edge_order=2)
    Tyx_fd, Tyy_fd = np.gradient(Ty, h, h, edge_order=2)

    plt.figure()
    plt.imshow(T, extent=[0, 1, 1, 0], vmin=(0, 1))
    # plt.imshow(np.log10(abs(eik0.T - T)), extent=[0, 1, 1, 0])
    plt.colorbar()
    plt.show()
