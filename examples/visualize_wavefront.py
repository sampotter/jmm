import jmm
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

from scipy.optimize import brentq

tetra_mesh_path = 'room/room.1.vtk'
tetra_mesh = pv.read(tetra_mesh_path)

surf_mesh_path = 'room.obj'
surf_mesh = pv.read(surf_mesh_path)

points = tetra_mesh.points.copy().astype(np.float64)
cells = tetra_mesh.cells.copy().astype(np.uint64)
cells = cells.reshape(cells.size//5, 5)[:, 1:]

num_points = points.shape[0]
num_cells = cells.shape[0]

print('%d points and %d cells' % (num_points, num_cells))

jets_path = 'room/jets.bin'
jets = np.fromfile(jets_path, np.float64, -1).reshape(num_points, 4)
T = jets[:, 0]
DT = jets[:, 1:]
del jets

Tslice = 10.0

# Find the tetrahedra that bracket the slice value

Tmin = T[cells].min(1)
Tmax = T[cells].max(1)
I = np.where((Tmin <= Tslice) & (Tslice <= Tmax))[0]

e = lambda n, i: np.eye(4)[i]

def get_tris_for_cell(Tslice, i):
    Tc = T[cells[i]]
    Tbb = jmm.Bb3Tet(Tc, DT[cells[i]], points[cells[i]])
    P = points[cells[i]]

    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    active = []
    Ps = []
    for i, j in edges:
        f = lambda s: Tbb.f((1 - s)*e(4, i) + s*e(4, j)) - Tslice
        active.append(np.sign(f(0)) != np.sign(f(1)))
        if (active[-1]):
            s = brentq(f, 0, 1)
            Ps.append((1 - s)*P[i] + s*P[j])
    active = np.array(active)
    Ps = np.array(Ps)

    if active.sum() == 4:
        F = np.array([
            [3, 0, 1, 2],
            [3, 1, 2, 3]
        ])
    else:
        assert False

    return pv.PolyData(Ps, F)

poly_data = get_tris_for_cell(Tslice, I[0])

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(poly_data, show_edges=True)
