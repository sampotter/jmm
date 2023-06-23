#!/usr/bin/env python

import pyvista as pv
import pyvistaqt as pvqt
import numpy as np
import vtk

V = np.fromfile('verts.bin').reshape(-1, 3)
C = np.fromfile('cells.bin', dtype=np.uintp).reshape(-1, 4)

v = np.array([0.025, -0.025, 0.05])

c = lambda x: 1 + v[0]*x[:, 0] + v[1]*x[:, 1] + v[2]*x[:, 2]
grad_c = v

s = lambda x: 1/c(x)
grad_s = lambda x: -np.outer(1/c(x)**2, grad_c)

grad_s(V)

grid = pv.UnstructuredGrid({vtk.VTK_TETRA: C}, V)
grid['c'] = c(V)
grid['s'] = s(V)

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(grid, scalars='c')

s_data = np.concatenate([s(V).reshape(-1, 1), grad_s(V)], axis=1)
s_data.tofile('s_data.bin')
