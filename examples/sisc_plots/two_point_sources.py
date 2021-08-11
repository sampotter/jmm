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

if __name__ == '__main__':
    # Speed function
    c = np.array([0.5, 5, 20])

    # Discretization parameters
    N = 65
    shape = np.array([N, N], dtype=np.uintp)
    xymin = np.zeros(2)
    h = 1/(N - 1)
    rfac = 0.1
    x, y = np.meshgrid(np.linspace(0, 1, N), np.linspace(0, 1, N))

    def get_eik_for_pt_src(x0, y0):
        s = LinearSpeedFunc2(c, x0, y0)
        eik = Eik.from_s_and_grid(s, shape, xymin, h)
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
