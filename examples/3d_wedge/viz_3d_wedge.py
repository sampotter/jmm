import colorcet as cc
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

FAR, TRIAL, VALID = 0, 1, 2

verts = np.fromfile('verts.bin', dtype=np.float64).reshape(-1, 3)
cells = np.fromfile('cells.bin', dtype=np.uintp).reshape(-1, 4)
mesh = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

#
# load data for direct eikonal
#

direct_jet = np.fromfile('direct_jet.bin', dtype=np.float64).reshape(-1, 13)
direct_T = direct_jet[:, 0]
direct_grad_T = direct_jet[:, 1:4]
direct_hess_T = direct_jet[:, 4:].reshape(-1, 3, 3)
direct_lap_T = direct_hess_T[:, 0, 0] + direct_hess_T[:, 1, 1] \
    + direct_hess_T[:, 2, 2]

direct_state = np.fromfile('direct_state.bin', dtype=np.intc)
# verts_direct_far = pv.PolyData(verts[direct_state == FAR])
# verts_direct_trial = pv.PolyData(verts[direct_state == TRIAL])
# verts_direct_valid = pv.PolyData(verts[direct_state == VALID])

direct_jet_gt = np.fromfile('direct_jet_gt.bin', dtype=np.float64).reshape(-1, 13)
direct_T_gt = direct_jet_gt[:, 0]
direct_grad_T_gt = direct_jet_gt[:, 1:4]
direct_hess_T_gt = direct_jet_gt[:, 4:].reshape(-1, 3, 3)
direct_lap_T_gt = direct_hess_T_gt[:, 0, 0] + direct_hess_T_gt[:, 1, 1] \
    + direct_hess_T_gt[:, 2, 2]

direct_error_jet = direct_jet - direct_jet_gt
direct_error_T = direct_error_jet[:, 0]
direct_error_grad_T = direct_error_jet[:, 1:4]
direct_error_hess_T = direct_error_jet[:, 4:].reshape(-1, 3, 3)

#
# do plotting
#

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(mesh, show_edges=True, opacity=0.25)
# plotter.add_mesh(verts_direct_far, color='red')
# plotter.add_mesh(verts_direct_trial, color='yellow')
# plotter.add_mesh(verts_direct_valid, color='cyan')

# use this to plot errors
values = direct_error_hess_T[:, 2, 2]
vmax = abs(values).max()
vmin = -vmax
clim = (vmin, vmax)
# clim = (-1e-1, 1e-1)
cmap = cc.cm.coolwarm

# # use this to plot values
# values = direct_hess_T[:, 1, 1]
# clim = (values.min(), values.max())
# cmap = cc.cm.bmy

pts = pv.PolyData(verts)
pts['values'] = values
plotter.add_mesh(pts, scalars='values', clim=clim, cmap=cmap)
