import itertools as it
import numpy as np

from jmm.eik import Eik2g1
from jmm.grid import Grid2
from jmm.jet import Jet22t

if __name__ == '__main__':
    n = 129
    h = 2/(n-1)
    shape = np.array([n, n], dtype=np.intc)
    xymin = np.array([-1, -1], dtype=np.float64)
    grid = Grid2(shape, xymin, h)
    eik = Eik2g1.from_grid(grid)

    rfac = 0.1
    X, Y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n), indexing='ij')

    #    r^2 = x^2 + y^2
    # => r rx = x
    # => rx = x/r
    #  & rx rx + r rxx = 1
    # => x^2/r^2 + r rxx = 1
    # => rxx = (1 - x^2/r^2)/r
    #  & r rxy + rx ry = 0
    # => rxy = -rx ry / r

    R = np.sqrt(X**2 + Y**2)
    Rx, Ry = X/R, Y/R
    Rxx = (1 - Rx**2)/R
    Ryx = -Ry*Rx/R
    Rxy = Ryx
    Ryy = (1 - Ry**2)/R

    valid_mask = R <= rfac

    trial_mask = np.zeros_like(valid_mask)
    shape = trial_mask.shape
    for i, j in it.product(*(range(dim) for dim in shape)):
        if valid_mask[i, j]:
            continue
        for di, dj in it.product([-1, 0, 1], repeat=2):
            if 0 <= i + di < shape[0] and \
               0 <= j + dj < shape[1] and \
               valid_mask[i + di, j + dj]:
                trial_mask[i, j] = True

    def get_jet_for_index(i, j):
        r = R[i, j]
        rx, ry = Rx[i, j], Ry[i, j]
        rxx, ryx, rxy, ryy = Rxx[i, j], Ryx[i, j], Rxy[i, j], Ryy[i, j]
        return Jet22t(r, rx, ry, rxx, ryx, rxy, ryy)

    for i, j in zip(*np.where(valid_mask)):
        ind = np.array([i, j], dtype=np.intc)
        jet = get_jet_for_index(i, j)
        eik.add_valid(ind, jet)

    for i, j in zip(*np.where(trial_mask)):
        ind = np.array([i, j], dtype=np.intc)
        jet = get_jet_for_index(i, j)
        eik.add_trial(ind, jet)

    eik.solve()
