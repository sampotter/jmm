#!/usr/bin/env python

import jmm
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

tet_mesh = 'L.1.vtk'
surf_mesh = pv.read('L.obj')

grid = pv.read(tet_mesh)
points = grid.points.copy().astype(np.float64)
cells = grid.cells.reshape(-1, 5)[:, 1:].copy().astype(np.uint64)
centroids = grid.points[cells].mean(1)

mask = centroids[:, 2] < 0.5
inds = mask.nonzero()[0]
subgrid = grid.extract_cells(inds)

plotter = pvqt.BackgroundPlotter()
plotter.background_color = (1, 1, 1)
plotter.add_mesh(surf_mesh, 'r', 'wireframe')
plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
plotter.add_mesh(pv.Sphere(0.05, np.zeros(3)), color='blue')
plotter.add_legend([['Input mesh', 'r'],
                    ['Tetrahedron mesh', 'black'],
                    ['Point source', 'blue']])
