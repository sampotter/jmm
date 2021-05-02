# Note: can sometimes be helpful to run this section first before
# running the rest of the script if running from python-mode in Emacs

import colorcet as cc
import jmm
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

import plotting

################################################################################
# parameters

vtu_path = None # 'room.vtu'

verts_path = 'Building_dom_verts.bin'
cells_path = 'Building_dom_cells.bin'

bc_path = 'refl_bcs.txt'

lsrc = 0 if bc_path is None else None
# l = [1657, 1384, 5658, 6228]
l = 1384
l0 = 2010
l1 = None
l2 = None
l3 = None
lbad = None

l_color = 'red'
l0_color = 'yellow'
l1_color = 'blue'
l2_color = 'green'
l3_color = 'black'
lbad_color = 'green'

plot_surf_tris = False
plot_wavefront = True
plot_ray_from_lsrc_to_l = False
plot_lsrc = False

################################################################################
# GEOMETRY SETUP

if vtu_path:
    grid = pv.read(vtu_path)

if vtu_path:
    verts = grid.points.astype(np.float64)
else:
    verts = np.fromfile(verts_path, dtype=np.float64)
    verts = verts.reshape(verts.size//3, 3)
num_verts = verts.shape[0]

if vtu_path:
    cells = grid.cells.reshape(-1, 5)[:, 1:].astype(np.uintp)
else:
    cells = np.fromfile(cells_path, dtype=np.uintp)
    cells = cells.reshape(cells.size//4, 4)

################################################################################
# COMPUTE FACES FOR SURFACE MESH

faces = set()
for C in cells:
    F = [tuple(sorted([C[0], C[1], C[2]])),
         tuple(sorted([C[0], C[1], C[3]])),
         tuple(sorted([C[0], C[2], C[3]])),
         tuple(sorted([C[1], C[2], C[3]]))]
    for f in F:
        if f in faces:
            faces.remove(f)
        else:
            faces.add(f)
faces = np.array(list(faces), dtype=np.uintp)

################################################################################
# SOLVE

mesh = jmm.Mesh3.from_verts_and_cells(verts, cells)

eik = jmm.Eik3(mesh)

# set up BCs
lsrcs = []
if bc_path is not None:
    with open(bc_path, 'r') as f:
        for line in f:
            strs = line.split()
            lsrc = int(strs[0])
            jet = jmm.Jet3(*(float(_) for _ in strs[1:]))
            eik.add_trial(lsrc, jet)
            lsrcs.append(lsrc)
        lsrc = None
else:
    eik.add_trial(lsrc, jmm.Jet3(0.0, np.nan, np.nan, np.nan))
    lsrcs.append(lsrc)

# if l0 is None:
#     eik.solve()
# else:
#     while eik.peek() != l0:
#         eik.step()

for _ in range(int(np.round(8*len(lsrcs)))):
    eik.step()

################################################################################
# HELPER FUNCTIONS FOR PLOTTING

plotter = pvqt.BackgroundPlotter()

def plot_tri(L):
    plotter.add_mesh(
        pv.make_tri_mesh(verts, np.array(L).reshape(1, 3)),
        opacity=0.95, color='white', show_edges=True)

def plot_cell(lc, **kwargs):
    plotter.add_mesh(
        pv.UnstructuredGrid(
            {vtk.VTK_TETRA: cells[lc].reshape(1, 4)},
            verts),
        **kwargs)

def plot_point(points, l, scale=1.25, color='white', opacity=1):
    plotter.add_mesh(
        pv.Sphere(scale*sphere_radius, points[l]),
        color=color, opacity=opacity)

def plot_jet(x, jet=None):
    if isinstance(jet, jmm.Jet3):
        d = np.array([jet.fx, jet.fy, jet.fz])
    elif len(jet) > 3:
        d = np.array([jet[1], jet[2], jet[3]])
    else:
        d = jet
    dtol = abs(1 - np.linalg.norm(d))
    if dtol >= 1e-13:
        print(dtol)
        assert False
    plotter.add_mesh(
        pv.Arrow(x, d, scale=0.1),
        color='white', opacity=1)

def plot_x(x, scale=1, **kwargs):
    plotter.add_mesh(
        pv.Sphere(scale*sphere_radius, x),
        **kwargs)

################################################################################
# MAKE PLOTS

surf_mesh = pv.PolyData(
    verts,
    np.concatenate([3*np.ones((faces.shape[0], 1), dtype=np.uintp), faces], axis=1)
)
plotter.add_mesh(surf_mesh, color='white', opacity=0.25, show_edges=plot_surf_tris)

highlight_inds = dict()
if lsrc is not None: highlight_inds[lsrc] = 'pink'
if l0 is not None:   highlight_inds[l0] = l0_color
if isinstance(l, list):
    for l_ in l:
        highlight_inds[l_] = l_color
else:
    highlight_inds[l] = l_color
if lbad is not None: highlight_inds[lbad] = lbad_color
if l1 is not None:   highlight_inds[l1] = l1_color
if l2 is not None:   highlight_inds[l2] = l2_color
if l3 is not None:   highlight_inds[l3] = l3_color

h = eik.mesh.min_tetra_alt/2
sphere_radius = h/6

if plot_wavefront:
    plotting.plot_wavefront(plotter, eik)

for l_, color in highlight_inds.items():
    plot_point(verts, l_, 1.5, color)

if plot_lsrc:
    for l_ in lsrcs:
        plot_point(verts, l_, 1, 'pink')

if plot_ray_from_lsrc_to_l:
    if l is None:
        raise Exception('l must be defined to plot ray connecting lsrc and l')
    x, xsrc = verts[l], verts[lsrc]
    xm, xd = (x + xsrc)/2, x - xsrc
    r, h = 0.5*sphere_radius, np.linalg.norm(xd)
    plotter.add_mesh(
        pv.Cylinder(xm, xd, r, h), color=highlight_inds[l], opacity=1)

############################################################################
