import itertools as it
import numpy as np

from jmm.eik import Eik
from jmm.field import LinearSpeedFunc2

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

def get_eik_for_pt_src(grid, c, x0, y0, rfac):
    s = LinearSpeedFunc2(c, x0, y0)
    eik = Eik.from_s_and_grid(s, grid)
    M, N = grid.shape
    x, y = np.meshgrid(np.linspace(grid.xmin, grid.xmax, M),
                       np.linspace(grid.ymin, grid.ymax, N),
                       indexing='ij')
    valid_mask = (x - x0)**2 + (y - y0)**2 < rfac**2
    for i, j in zip(*np.where(valid_mask)):
        ind = np.array([i, j], dtype=np.intc)
        eik.add_valid(ind, s.get_jet(x[i, j], y[i, j]))
    for i, j in zip(*np.where(get_trial_mask(valid_mask))):
        ind = np.array([i, j], dtype=np.intc)
        eik.add_trial(ind, s.get_jet(x[i, j], y[i, j]))
    eik.build_cells()
    return s, eik
