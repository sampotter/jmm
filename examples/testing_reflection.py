import colorcet as cc
import json
import numpy as np
import pyvista as pv
import time
import vtk

############################################################################

import jmm
import plotting

############################################################################

t0 = np.inf

def tic():
    global t0
    t0 = time.time()

def toc():
    global t0
    return time.time() - t0

############################################################################

lsrc = 0

should_plot_reflectors = False
should_plot_BCs = False
should_plot_refl_t_in = True
should_plot_refl_t_out = True
should_plot_refl_hists = True
should_plot_refl_angle_cond = True

############################################################################

tic()

surf_grid = pv.read(f'sethian_shadow/Building.obj')
grid = pv.read(f'sethian_shadow/Building.vtu')
with open('sethian_shadow/Building.json', 'r') as f:
    info = json.load(f)

print(f'- loaded scene [%1.2fs]' % toc())

############################################################################

cells = grid.cells.reshape(-1, 5)[:, 1:].astype(np.uintp)
points = grid.points.astype(np.float64)

num_dom_cells = info['num_dom_cells']
num_dom_points = info['num_dom_points']

dom_cells = cells[:num_dom_cells]
dom_points = points[:num_dom_points]

mesh_dom = jmm.Mesh3.from_verts_and_cells(dom_points, dom_cells)

mesh = jmm.Mesh3.from_verts_and_cells(points, cells)
for le in mesh_dom.get_diff_edges():
    mesh.set_boundary_edge(*le, True)

print(f'- set up meshes [%1.2fs]' % toc())

############################################################################

tic()

eik = jmm.Eik3(mesh)
eik.add_trial(lsrc, jmm.Jet3.make_point_source(0))
eik.solve()

print('- solved on extended domain [%1.2fs]' % toc())

tic()

eik_dom = jmm.Eik3(mesh_dom)
eik_dom.add_trial(lsrc, jmm.Jet3.make_point_source(0))
eik_dom.solve()

print('- solved on domain [%1.2fs]' % toc())

tic()

T = np.array([_[0] for _ in eik.jet[:num_dom_points]])
T_diff = np.array([_[0] for _ in eik_dom.jet])

T[~np.isfinite(T)] = -abs(T[np.isfinite(T)]).max()
T_diff[~np.isfinite(T_diff)] = -abs(T_diff[np.isfinite(T_diff)]).max()

D = T_diff - T
Z = D > 0.007

print('- found shadow [%1.2fs]' % toc())

############################################################################

if should_plot_reflectors:
    plotter = pv.Plotter()
    plotter.background_color = 'white'

    nrefl = mesh_dom.num_reflectors
    cmap = cc.cm.glasbey
    clim = (0, nrefl - 1)

    for i in range(nrefl):
        faces = mesh_dom.get_reflector(i)
        mask = (~Z[faces]).any(1)
        grid = pv.UnstructuredGrid({vtk.VTK_TRIANGLE: faces}, dom_points)
        grid.cell_arrays['label'] = i*np.ones(faces.shape[0])
        # grid.cell_arrays['label'][mask] *= 2
        # grid.cell_arrays['label'][mask] += 1
        plotter.add_mesh(grid, scalars='label', cmap=cmap, clim=clim,
                         show_scalar_bar=False)

    plotter.show()

    del plotter

############################################################################

tic()

i, faces = next(mesh_dom.get_active_reflectors(Z))
L = np.unique(faces)
L = L[~Z[L]]

x0, x1, x2 = points[faces[0]]
n = np.cross(x1 - x0, x2 - x0)
n /= np.linalg.norm(n)
R = np.eye(n.size) - 2*np.outer(n, n)

Tau = np.array([_[0] for _ in eik_dom.jet[L]])
DTau = np.array([(_[1], _[2], _[3]) for _ in eik_dom.jet[L]])
DTau = DTau@R

print('- computed BCs for reflection [%1.2fs]' % toc())

############################################################################

bcs_path = 'refl_bcs.txt'

with open(bcs_path, 'w') as f:
    for l, tau, dtau in zip(L, Tau, DTau):
        print('%lu %1.16g %1.16g %1.16g %1.16g' % (l, tau, *dtau), file=f)

print('- wrote BCs to %s' % bcs_path)

############################################################################

if should_plot_BCs:
    plotter = pv.BackgroundPlotter()
    plotter.background_color = 'white'

    plotter.add_mesh(
        pv.UnstructuredGrid({vtk.VTK_TRIANGLE: faces}, dom_points))
    for l_, dtau in zip(L, DTau):
        p = dom_points[l_]
        h = 0.1
        plotter.add_mesh(pv.Sphere(h, p), color='white')
        plotter.add_mesh(pv.Arrow(p, dtau), color='white')
    plotter.add_mesh(pv.Sphere(h, dom_points[lsrc]), color='red')

    plotter.show()

    del plotter

############################################################################

tic()
eik_refl, eik_dom_refl = jmm.Eik3(mesh), jmm.Eik3(mesh_dom)
for lf in faces:
    assert (lf < num_dom_points).all()
    if Z[lf].any():
        continue
    jets = [jmm.Jet3(*_) for _ in eik.jet[lf]]
    eik_refl.add_valid_bdf(*lf, *jets)
    eik_dom_refl.add_valid_bdf(*lf, *jets)
print('- set up BCs for reflection [%1.2fs]' % toc())

tic()
eik_refl.solve()
print('- computed reflection on extended domain [%1.2fs]' % toc())

plotter = pv.Plotter()
plotter.background_color = 'white'
plotting.plot_wavefront(plotter, eik_refl)
plotter.show()

tic()
eik_dom_refl.solve()
print('- computed reflection on original domain [%1.2fs]' % toc())

tic()

T = np.array([_[0] for _ in eik_refl.jet[:num_dom_points]])
T_diff = np.array([_[0] for _ in eik_dom_refl.jet])

T[~np.isfinite(T)] = -abs(T[np.isfinite(T)]).max()
T_diff[~np.isfinite(T_diff)] = -abs(T_diff[np.isfinite(T_diff)]).max()

D = T_diff - T
Z = D > 0.007

print('- found shadow for reflection [%1.2fs]' % toc())

dom_grid = pv.UnstructuredGrid({vtk.VTK_TETRA: dom_cells}, dom_points)
dom_grid.point_arrays['Z'] = Z

plotter = pv.Plotter()
plotter.background_color = 'white'
plotter.add_mesh(dom_grid, scalars='Z', cmap=cc.cm.gray_r,
                 show_scalar_bar=False)
plotter.show()

if should_plot_refl_t_in:
    plotter = pv.Plotter()
    plotter.background_color = 'white'
    dom_grid.vectors = (1 - Z).reshape(num_dom_points, 1)*eik_dom_refl.t_in
    dom_grid.vectors /= 2
    plotter.add_mesh(
        dom_grid.arrows,cmap=cc.cm.gray,lighting=False,show_scalar_bar=False)
    plotter.show()

if should_plot_refl_t_out:
    plotter = pv.Plotter()
    plotter.background_color = 'white'
    dom_grid.vectors = (1 - Z).reshape(num_dom_points, 1)*eik_dom_refl.t_out
    dom_grid.vectors /= 2
    plotter.add_mesh(
        dom_grid.arrows,cmap=cc.cm.gray,lighting=False,show_scalar_bar=False)
    plotter.show()
