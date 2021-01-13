import pyvista as pv
import pyvistaqt as pvqt

tet_mesh = 'room/room.1.vtk'
surf_mesh = pv.read('room.obj')

grid = pv.read(tet_mesh)
cells = grid.cells.reshape(-1, 5)[:, 1:]
points = grid.points.copy().astype(np.float64)
cells = grid.cells.reshape(-1, 5)[:, 1:].copy().astype(np.uint64)
centroids = grid.points[cells].mean(1)

mask = centroids[:, 1] < 2
inds = mask.nonzero()[0]
subgrid = grid.extract_cells(inds)

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(surf_mesh, 'r', 'wireframe')
plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
plotter.add_legend([[' Input Mesh ', 'r'],
                    [' Tesselated Mesh ', 'black']])
plotter.show()
