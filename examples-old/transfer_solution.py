import colorcet as cc
import jmm
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

plt.ion()

tetra_mesh_path = 'room/room.1.vtk'
tetra_mesh = pv.read(tetra_mesh_path)

surf_mesh_path = 'room.obj'
surf_mesh = pv.read(surf_mesh_path)

points = tetra_mesh.points.copy().astype(np.float64)
cells = tetra_mesh.cells.astype(np.uint64)
cells = cells.reshape(cells.size//5, 5)[:, 1:].copy()

points = np.array(points)
cells = np.array(cells)

num_points = points.shape[0]
num_cells = cells.shape[0]

print('%d points and %d cells' % (num_points, num_cells))

mesh_eps = meshplex.MeshTetra(points, cells).edge_lengths.min()
mesh = jmm.Mesh3.from_verts_and_cells(points, cells, mesh_eps)
eik = jmm.Eik3(mesh)
eik.add_trial(3019, jmm.Jet3(0, np.nan, np.nan, np.nan))
eik.solve()

pmin = points.min(0) - 0.5
pmax = points.max(0) + 0.5
h = 0.15
dim = np.ceil((pmax - pmin)/h).astype(np.intc)
grid = jmm.Grid3(dim, pmin, h)

T = eik.transfer_solution_to_grid(grid)

plt.figure()
plt.imshow(T[:, 10, :], extent=[pmin[0], pmax[0], pmin[2], pmax[2]], cmap=cc.cm.rainbow)
plt.colorbar()
plt.gca().set_aspect('equal')
plt.show()
