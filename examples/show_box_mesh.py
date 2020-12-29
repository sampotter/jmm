import pyvista as pv
import pyvistaqt as pvqt

root = 'box/0.001/box'

path = '%s.1.vtk' % root

surf_mesh = pv.read('%s.obj' % root)

grid = pv.read(path)

cells = grid.cells.reshape(-1, 5)[:, 1:]
centroids = grid.points[cells].mean(1)

mask = centroids[:, 0] < 0
inds = mask.nonzero()[0]
subgrid = grid.extract_cells(inds)

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(surf_mesh, 'r', 'wireframe')
plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
plotter.add_legend([[' Input Mesh ', 'r'],
                    [' Tesselated Mesh ', 'black']])
plotter.show()
