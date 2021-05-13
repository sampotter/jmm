import colorcet as cc
import json
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import time
import vtk

plt.ion()

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

should_plot_diffractors = False
should_plot_wavefront = False
should_plot_BCs = False
should_plot_diff_shadow_mask = True
should_plot_diff_t_in = True
should_plot_diff_t_out = True
should_plot_diff_hists = True
should_plot_diff_angle_cond = True

############################################################################

h = mesh.min_edge_length
r = h/2

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

def plot_diffractor(plotter, i):
    edges = mesh_dom.get_diffractor(i)
    for l0, l1 in edges:
        x0, x1 = points[l0], points[l1]
        p, d = (x0 + x1)/2, x1 - x0
        H = np.linalg.norm(d)
        d /= H
        plotter.add_mesh(pv.Cylinder(p, d, r, H), color=colors[i])

if should_plot_diffractors:
    num_diff = mesh_dom.num_diffractors

    colors = np.random.choice(
        list(pv.colors.hexcolors.keys()), num_diff, replace=False)

    plotter = pv.Plotter()
    plotter.background_color = 'white'
    for i in range(num_diff):
        plot_diffractor(plotter, i)
    plotter.show()

    del plotter

############################################################################

tic()

i, edges = next(mesh_dom.get_active_diffractors(Z))
L = np.unique(edges)
L = L[~Z[L]]

Tau = np.array([_[0] for _ in eik.jet[L]])
DTau = np.array([(_[1], _[2], _[3]) for _ in eik.jet[L]])

print('- computed BCs for diffraction [%1.2fs]' % toc())

############################################################################

bcs_path = 'diff_bcs.txt'
with open(bcs_path, 'w') as f:
    for l, tau, t_in in zip(L, Tau, DTau):
        print('%lu %1.16g %1.16g %1.16g %1.16g' % (l, tau, *t_in), file=f)
print('- wrote BCs to %s' % bcs_path)

edges_path = 'diff_edges.txt'
with open(edges_path, 'w') as f:
    for l0, l1 in edges:
        print('%lu %lu' % (l0, l1), file=f)
print('- wrote diffracting edges to %s' % edges_path)

############################################################################

tic()
eik_diff, eik_dom_diff = jmm.Eik3(mesh), jmm.Eik3(mesh_dom)
for l0, l1 in edges:
    if Z[l0] or Z[l1]: continue
    jet0 = jmm.Jet3(*eik.jet[l0])
    jet1 = jmm.Jet3(*eik.jet[l1])
    eik_diff.add_valid_bde(l0, l1, jet0, jet1)
    eik_dom_diff.add_valid_bde(l0, l1, jet0, jet1)

print('- set up BCs for diffraction [%1.2fs]' % toc())

tic()
eik_diff.solve()
print('- computed diffraction on extended domain [%1.2fs]' % toc())

if should_plot_wavefront:
    plotter = pv.Plotter()
    plotter.background_color = 'white'
    plotting.plot_wavefront(plotter, eik_diff)
    plotter.show()

tic()
eik_dom_diff.solve()
print('- computed diffraction on original domain [%1.2fs]' % toc())

tic()

T = np.array([_[0] for _ in eik_diff.jet[:num_dom_points]])
T_diff = np.array([_[0] for _ in eik_dom_diff.jet])

T[~np.isfinite(T)] = -abs(T[np.isfinite(T)]).max()
T_diff[~np.isfinite(T_diff)] = -abs(T_diff[np.isfinite(T_diff)]).max()

D = T_diff - T
Z = D > 0.0033

print('- found shadow for diffraction [%1.2fs]' % toc())

dom_grid = pv.UnstructuredGrid({vtk.VTK_TETRA: dom_cells}, dom_points)
dom_grid.point_arrays['Z'] = Z

if should_plot_diff_shadow_mask:
    plotter = pv.Plotter()
    plotter.background_color = 'white'
    plotter.add_mesh(dom_grid, scalars='Z', cmap=cc.cm.gray_r,
                     show_scalar_bar=False)
    plot_diffractor(plotter, 0)
    plotter.show()

if should_plot_diff_t_in:
    plotter = pv.Plotter()
    plotter.background_color = 'white'
    dom_grid.vectors = (1 - Z).reshape(num_dom_points, 1)*eik_dom_diff.t_in
    dom_grid.vectors /= 2
    plotter.add_mesh(
        dom_grid.arrows,cmap=cc.cm.gray,lighting=False,show_scalar_bar=False)
    plot_diffractor(plotter, 0)
    plotter.show()

if should_plot_diff_t_out:
    plotter = pv.Plotter()
    plotter.background_color = 'white'
    dom_grid.vectors = (1 - Z).reshape(num_dom_points, 1)*eik_dom_diff.t_out
    dom_grid.vectors /= 2
    plotter.add_mesh(
        dom_grid.arrows,cmap=cc.cm.gray,lighting=False,show_scalar_bar=False)
    plot_diffractor(plotter, 0)
    plotter.show()

# Check whether diffraction angle condition holds
x0, x1 = points[edges[0]]
dx = x1 - x0
lhs = np.sqrt(1 - (eik_dom_diff.t_in@dx)**2)
rhs = np.sqrt(1 - (eik_dom_diff.t_out@dx)**2)

if should_plot_diff_hists:
    line_x = 2*h**3
    bins = 128 + 1
    plt.figure()
    plt.hist(lhs, bins=bins, label=r'$t_{in}^\top \delta x$', histtype='step')
    plt.hist(rhs, bins=bins, label=r'$t_{out}^\top \delta x$', histtype='step')
    plt.hist(lhs - rhs, bins=bins, label=r'$t_{in}^\top \delta x - t_{out}^\top \delta x$', histtype='step')
    plt.legend()
    plt.axvline(x=line_x, linestyle='--', linewidth=1, c='k')
    plt.axvline(x=-line_x, linestyle='--', linewidth=1, c='k')
    plt.show()

if should_plot_diff_angle_cond:
    dom_grid.point_arrays['diff_angle_cond'] = abs(lhs - rhs)
    dom_grid.point_arrays['diff_angle_cond'][L] = 0
    dom_grid.point_arrays['diff_angle_cond'][Z] = np.nan
    clipped = dom_grid.clip(origin=points[lsrc], normal=(0, 1, 0))
    plotter = pv.Plotter()
    plotter.add_mesh(pv.Sphere(2*r, points[lsrc]), color='cyan')
    plotter.add_mesh(clipped, scalars='diff_angle_cond', cmap=cc.cm.gray_r,
                     clim=(0, h**2), lighting=False)
    plot_diffractor(plotter, 0)
    plotter.show()
