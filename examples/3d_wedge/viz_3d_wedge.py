import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

FAR, TRIAL, VALID = 0, 1, 2

verts = np.fromfile('verts.bin', dtype=np.float64).reshape(-1, 3)
cells = np.fromfile('cells.bin', dtype=np.uintp).reshape(-1, 4)

direct_jet = np.fromfile('direct_jet.bin', dtype=np.float64).reshape(-1, 13)
direct_T = direct_jet[:, 0]
direct_grad_T = direct_jet[:, 1:4]
direct_hess_T = direct_jet[:, 4:].reshape(-1, 3, 3)

direct_state = np.fromfile('direct_state.bin', dtype=np.intc)

mesh = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

verts_direct_far = pv.PolyData(verts[direct_state == FAR])
verts_direct_trial = pv.PolyData(verts[direct_state == TRIAL])
verts_direct_valid = pv.PolyData(verts[direct_state == VALID])

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(mesh, show_edges=True, opacity=0.25)
# plotter.add_mesh(verts_direct_far, color='red')
plotter.add_mesh(verts_direct_trial, color='yellow')
plotter.add_mesh(verts_direct_valid, color='green')
