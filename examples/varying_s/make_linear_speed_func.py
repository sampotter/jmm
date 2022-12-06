#!/usr/bin/env python

import pyvista as pv
import pyvistaqt as pvqt
import numpy as np
import vtk

V = np.fromfile('verts.bin').reshape(-1, 3)
C = np.fromfile('cells.bin', dtype=np.uintp).reshape(-1, 4)

dtheta_dy = 1

theta = lambda x: dtheta_dy*x[:, 1]
grad_theta = np.array([0, dtheta_dy, 0])

T = lambda x: 273.15 + theta(x)
F = lambda x: 32 + 9*theta(x)/5

c = lambda x: 331.3 + 0.606*theta(x)
grad_c = 0.606*grad_theta

s = lambda x: 1/c(x)
grad_s = lambda x: -np.outer(1/c(x)**2, grad_c)

grad_s(V)

grid = pv.UnstructuredGrid({vtk.VTK_TETRA: C}, V)
grid['c'] = c(V)
grid['s'] = s(V)

# plotter = pvqt.BackgroundPlotter()
# plotter.add_mesh(grid, scalars='c')

s_data = np.concatenate([s(V).reshape(-1, 1), grad_s(V)], axis=1)
s_data.tofile('s_data.bin')
