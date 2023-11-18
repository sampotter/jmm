import sys; sys.path.insert(-1, '../../wrappers/python')

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

import scipy.interpolate
import scipy.spatial

import jmm

verbose = True

c = 343
# xsrc = np.array([-7.5, -7.5, 1.25]) # room.off
xsrc = np.array([3, 2, 2], dtype=np.float64) # room_small.off
# xtgt = np.array([3, 6, 2], dtype=np.float64) # room_small.off - direct
xtgt = np.array([0.6, 7, 2], dtype=np.float64) # room_small.off - indirect
rfac = 0.5
eps = 1e-5
maxvol = 1e-2
R = 0.9
off_path = '../../examples/data/off/room_small.off'

surf_mesh = pv.read(off_path)

# plotter = pvqt.BackgroundPlotter()
# plotter.background_color = 'white'
# plotter.add_mesh(surf_mesh, show_edges=False, opacity=0.15)
# plotter.add_mesh(pv.Sphere(0.5, xsrc), color='green')

mesh_data = jmm.Mesh3Data.from_off(off_path, maxvol, verbose)
mesh_data.insert_vert(xsrc, eps)

mesh = jmm.Mesh3(mesh_data, eps=eps)

hh = jmm.Eik3hh.new_with_pt_src(mesh, c, rfac, xsrc)

root_branch = hh.get_root_branch()
root_branch.solve(verbose)

print('hi')

branches = [(0, root_branch)]

for i in range(16):
    depth, parent_branch = branches[i]
    for refl_ind in parent_branch.get_visible_refls():
        refl_branch = parent_branch.add_refl(refl_ind)
        refl_branch.solve(verbose)
        branches.append((depth + 1, refl_branch))

taps = []

for d, branch in branches:
    eik = branch.get_eik()
    T = eik.build_T_bmesh()
    spread = branch.spread
    org = branch.org

    def get_bary_coords(x, y, z):
        lc = mesh.find_cell_containing_point((x, y, z))
        if lc >= mesh.ncells:
            return
        A = np.empty((4, 4))
        A[:3] = mesh.verts[mesh.cells[lc]].T
        A[3] = 1
        b = np.array([x, y, z, 1])
        theta = np.linalg.solve(A, b)
        return theta, lc

    def get_spread(x, y, z):
        _ = get_bary_coords(x, y, z)
        if _ is None:
            return np.nan
        theta, lc = _
        return theta@spread[mesh.cells[lc]]

    def get_org(x, y, z):
        _ = get_bary_coords(x, y, z)
        if _ is None:
            return np.nan
        theta, lc = _
        return theta@org[mesh.cells[lc]]

    rho = T(*xtgt)
    t = rho/c

    print('multiply by R^depth!')
    A = R**d*get_org(*xtgt)/rho

    taps.append((t, A))

taps = sorted(taps, key=lambda _: _[0])
taps = np.array(taps)

plt.figure()
plt.stem(*taps.T)
plt.xlim(0, 0.05)
plt.show()

assert False

############################################################################
# PYVISTA STUFF BELOW:

h = 1/8
xmin, ymin, zmin = mesh.verts.min(0)
xmax, ymax, zmax = mesh.verts.max(0)
dx, dy, dz = mesh.verts.ptp(0)
xlin = np.linspace(xmin, xmax, int(np.ceil(dx/h)))
ylin = np.linspace(ymin, ymax, int(np.ceil(dy/h)))
zlin = np.linspace(zmin, zmax, int(np.ceil(dz/h)))

# f = lambda x, y, z: 20*np.log10(get_spread(x, y, z))
f = lambda x, y, z: T(x, y, z)

X, Y, Z = np.meshgrid(xlin, ylin, [xsrc[2]], indexing='ij')
fxy = np.array([f(x, y, z) for x, y, z in zip(X.ravel(), Y.ravel(), Z.ravel())]).reshape(X.shape)
gridxy = pv.StructuredGrid(X, Y, Z)

X, Y, Z = np.meshgrid(xlin, [xsrc[1]], zlin, indexing='ij')
fxz = np.array([f(x, y, z) for x, y, z in zip(X.ravel(), Y.ravel(), Z.ravel())]).reshape(X.shape)
gridxz = pv.StructuredGrid(X, Y, Z)

X, Y, Z = np.meshgrid([xsrc[0]], ylin, zlin, indexing='ij')
fyz = np.array([f(x, y, z) for x, y, z in zip(X.ravel(), Y.ravel(), Z.ravel())]).reshape(X.shape)
gridyz = pv.StructuredGrid(X, Y, Z)

gridxy['f'] = fxy.T.ravel()
gridxz['f'] = fxz.T.ravel()
gridyz['f'] = fyz.T.ravel()

contours = np.linspace(0, max(np.nanmax(fxy), np.nanmax(fxz), np.nanmax(fyz)), 11)
contours += (contours[1] - contours[0])/2
contours = contours[:-1]

# grid = pv.UnstructuredGrid({vtk.VTK_TETRA: mesh.cells}, mesh.verts)
surf = pv.read(off_path)
surf = pv.PolyData(surf.points, surf.cells.reshape(-1, 4)[:, 1:])

plotter = pvqt.BackgroundPlotter()
# plotter.background_color = 'white'
# plotter.add_mesh(surf_mesh, show_edges=False, opacity=0.15)
plotter.add_mesh(gridxy, nan_opacity=0, opacity=1)
plotter.add_mesh(gridxz, nan_opacity=0, opacity=1)
plotter.add_mesh(gridyz, nan_opacity=0, opacity=1)
plotter.add_mesh(gridxy.contour(contours), color='white', line_width=2)
plotter.add_mesh(gridxz.contour(contours), color='white', line_width=2)
plotter.add_mesh(gridyz.contour(contours), color='white', line_width=2)
