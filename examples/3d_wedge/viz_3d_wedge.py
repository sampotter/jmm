import colorcet as cc
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

FAR, TRIAL, VALID = 0, 1, 2

verts = np.fromfile('verts.bin', dtype=np.float64).reshape(-1, 3)
cells = np.fromfile('cells.bin', dtype=np.uintp).reshape(-1, 4)
mesh = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

with open('spec.txt', 'r') as f:
    spec = {k: v.strip() for k, v in map(lambda s: s.split(':'), f)}

for k, Type in {
        'verbose': bool,
        'visualize': bool,
        'maxvol': float,
        'n': float,
        'w': float,
        'h': float,
        'R': float}.items():
    spec[k] = Type(spec[k])

print('problem specification (3d wedge):')
for k, v in spec.items():
    print(f'- {k}: {v}')

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
# load data for o-face reflection
#

o_refl_jet = np.fromfile('o_refl_jet.bin', dtype=np.float64).reshape(-1, 13)
o_refl_T = o_refl_jet[:, 0]
o_refl_grad_T = o_refl_jet[:, 1:4]
o_refl_hess_T = o_refl_jet[:, 4:].reshape(-1, 3, 3)
o_refl_lap_T = o_refl_hess_T[:, 0, 0] + o_refl_hess_T[:, 1, 1] \
    + o_refl_hess_T[:, 2, 2]

# o_refl_jet_gt = np.fromfile('o_refl_jet_gt.bin', dtype=np.float64).reshape(-1, 13)
# o_refl_T_gt = o_refl_jet_gt[:, 0]
# o_refl_grad_T_gt = o_refl_jet_gt[:, 1:4]
# o_refl_hess_T_gt = o_refl_jet_gt[:, 4:].reshape(-1, 3, 3)
# o_refl_lap_T_gt = o_refl_hess_T_gt[:, 0, 0] + o_refl_hess_T_gt[:, 1, 1] \
#     + o_refl_hess_T_gt[:, 2, 2]

# o_refl_error_jet = o_refl_jet - o_refl_jet_gt
# o_refl_error_T = o_refl_error_jet[:, 0]
# o_refl_error_grad_T = o_refl_error_jet[:, 1:4]
# o_refl_error_hess_T = o_refl_error_jet[:, 4:].reshape(-1, 3, 3)

#
# load data for n-face reflection
#

# n_refl_jet = np.fromfile('n_refl_jet.bin', dtype=np.float64).reshape(-1, 13)
# n_refl_T = n_refl_jet[:, 0]
# n_refl_grad_T = n_refl_jet[:, 1:4]
# n_refl_hess_T = n_refl_jet[:, 4:].reshape(-1, 3, 3)
# n_refl_lap_T = n_refl_hess_T[:, 0, 0] + n_refl_hess_T[:, 1, 1] \
#     + n_refl_hess_T[:, 2, 2]

# n_refl_jet_gt = np.fromfile('n_refl_jet_gt.bin', dtype=np.float64).reshape(-1, 13)
# n_refl_T_gt = n_refl_jet_gt[:, 0]
# n_refl_grad_T_gt = n_refl_jet_gt[:, 1:4]
# n_refl_hess_T_gt = n_refl_jet_gt[:, 4:].reshape(-1, 3, 3)
# n_refl_lap_T_gt = n_refl_hess_T_gt[:, 0, 0] + n_refl_hess_T_gt[:, 1, 1] \
#     + n_refl_hess_T_gt[:, 2, 2]

# n_refl_error_jet = n_refl_jet - n_refl_jet_gt
# n_refl_error_T = n_refl_error_jet[:, 0]
# n_refl_error_grad_T = n_refl_error_jet[:, 1:4]
# n_refl_error_hess_T = n_refl_error_jet[:, 4:].reshape(-1, 3, 3)

#
# do plotting
#

# # use this to plot errors
# values = direct_error_hess_T[:, 2, 2]
# vmax = abs(values).max()
# vmin = -vmax
# clim = (vmin, vmax)
# # clim = (-1e-1, 1e-1)
# cmap = cc.cm.coolwarm

# use this to plot values
values = o_refl_T.copy()
values[np.isinf(values)] = np.nan
clim = (np.nanmin(values), np.nanmax(values))
cmap = cc.cm.bmy

pts = pv.PolyData(verts)
pts['values'] = values

h = 0.025

l = 2313
l0 = 245
l1 = 25

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(mesh, show_edges=True, opacity=0.25)
plotter.add_mesh(pts, scalars='values', clim=clim, cmap=cmap)
plotter.add_mesh(pv.Sphere(h, verts[l]), 'chartreuse')
plotter.add_mesh(pv.Sphere(h, verts[l0]), 'blue')
plotter.add_mesh(pv.Sphere(h, verts[l1]), 'cyan')
