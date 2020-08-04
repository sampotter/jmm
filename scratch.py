import itertools as it
import jmm
import numpy as np

n = 65
h = 2/(n-1)
shape = (n, n, n)
ind0 = (n//2, n//2, n//2)

dial = jmm.Dial(jmm.Stype.Constant, shape, h)
dial.add_point_source_with_trial_nbs(ind0, 0)

inds = np.array([_.flatten() for _ in np.meshgrid(range(48), range(48), range(48))]).T
dial.add_boundary_points(inds.astype(np.intc))
print(dial.state)

# dial.solve()

# L = np.linspace(-1, 1, n)
# X, Y, Z = np.meshgrid(L, L, L)
# tau = np.sqrt(X**2 + Y**2 + Z**2)

# print(abs(tau - dial.T).max())
