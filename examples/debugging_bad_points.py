# Note: can sometimes be helpful to run this section first before
# running the rest of the script if running from python-mode in Emacs

import colorcet as cc
import jmm
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

import plotting

############################################################################
# parameters

vtu_path = None # 'room.vtu'

verts_path = 'Building_verts.bin'
cells_path = 'Building_cells.bin'
bc_path = 'bcs.pickle'

dom_verts_path = 'Building_dom_verts.bin'
dom_cells_path = 'Building_dom_cells.bin'

lsrc = 0 if bc_path is None else None
l = None
l0 = 6121
l1 = None
l2 = None
l3 = None
lbad = None

l_color = 'red'
l0_color = 'purple'
l1_color = 'blue'
l2_color = 'green'
l3_color = 'black'
lbad_color = 'green'

plot_surf_tris = False
plot_wavefront = True
plot_ray_from_lsrc_to_l = False
plot_lsrc = False
plot_diffractors = True

############################################################################
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

############################################################################
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

############################################################################
# SOLVE

mesh = jmm.Mesh3.from_verts_and_cells(verts, cells)
for le in mesh.get_diff_edges():
    mesh.set_boundary_edge(*le, True)

dom_mesh = None
if dom_verts_path is not None:
    assert dom_cells_path is not None
    dom_verts = np.fromfile(dom_verts_path, dtype=np.float64)
    dom_verts = dom_verts.reshape(dom_verts.size//3, 3)
    dom_cells = np.fromfile(dom_cells_path, dtype=np.uintp)
    dom_cells = dom_cells.reshape(dom_cells.size//4, 4)
    dom_mesh = jmm.Mesh3.from_verts_and_cells(dom_verts, dom_cells)
    diff_edges = list(dom_mesh.get_diff_edges())
    assert len(diff_edges) > 0
    for le in diff_edges:
        mesh.set_boundary_edge(*le, True)

eik = jmm.Eik3(mesh)

# set up BCs
lsrcs = []
if bc_path is not None:
    ext = bc_path.split('.')[-1]
    if ext == 'txt':
        print('- bc file extension is "txt": reflection')
        with open(bc_path, 'r') as f:
            for line in f:
                strs = line.split()
                lsrc = int(strs[0])
                jet = jmm.Jet3(*(float(_) for _ in strs[1:]))
                eik.add_trial(lsrc, jet)
                lsrcs.append(lsrc)
            lsrc = None
    elif ext == 'pickle':
        print('- bc file extension is "pickle": diffraction')
        import pickle
        with open(bc_path, 'rb') as f:
            bcs = pickle.load(f)
        assert isinstance(bcs, dict)
        edges = list(bcs.keys())
        for le in edges:
            f, Df, x = bcs[le]
            jet0 = jmm.Jet3(f[0], *Df[0])
            jet1 = jmm.Jet3(f[1], *Df[1])
            eik.add_valid_bde(*le, jet0, jet1)
else:
    eik.add_trial(lsrc, jmm.Jet3(0.0, np.nan, np.nan, np.nan))
    lsrcs.append(lsrc)

# if l0 is None:
#     eik.solve()
# else:
#     while eik.peek() != l0:
#         eik.step()

for _ in range(1000):
    eik.step()

############################################################################
# HELPER FUNCTIONS FOR PLOTTING

h = eik.mesh.min_tetra_alt/2
sphere_radius = h/6

plotter = pvqt.BackgroundPlotter()
plotter.background_color = 'white'

def plot_tri(L, **kwargs):
    plotter.add_mesh(
        pv.make_tri_mesh(verts, np.array(L).reshape(1, 3)), **kwargs)

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

def plot_line(plotter, x, y, scale=1, color='white'):
    plot_x(x, scale=scale, color=color)
    plot_x(y, scale=scale, color=color)
    xm = (x + y)/2
    xd = x - y
    d = np.linalg.norm(xd)
    xd /= d
    r = scale*sphere_radius
    plotter.add_mesh(
        pv.Cylinder(xm, xd, scale*sphere_radius, d, capping=False),
        color=color)

############################################################################
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
elif l is not None:
    highlight_inds[l] = l_color
if lbad is not None: highlight_inds[lbad] = lbad_color
if l1 is not None:   highlight_inds[l1] = l1_color
if l2 is not None:   highlight_inds[l2] = l2_color
if l3 is not None:   highlight_inds[l3] = l3_color

if plot_wavefront:
    plotting.plot_wavefront(plotter, eik, opacity=1)

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

colors = np.random.choice(list(pv.colors.hexcolors.keys()), 28, replace=False)
if plot_diffractors:
    for i in range(mesh.num_diffractors):
        edges = mesh.get_diffractor(i)
        for le in edges:
            x0, x1 = verts[le]
            p, d = (x0 + x1)/2, x1 - x0
            r, h = sphere_radius, np.linalg.norm(d)
            d /= h
            plotter.add_mesh(pv.Cylinder(p, d, r, h), color=colors[i])

############################################################################

if dom_mesh is not None:
    for le in dom_mesh.get_diffractor(0):
        plot_line(plotter, *verts[le], color='black')

lvalid = []
for m in range(num_verts):
    if eik.is_valid(m):
        plot_point(verts, m, 0.75, 'green')
        plot_jet(verts[m], eik.jet[m])
        if m not in lvalid:
            lvalid.append(m)
lvalid = np.array(lvalid)

ltrial = []
for m0 in lvalid:
    for m1 in mesh.vv(m0):
        if eik.is_valid(m1):
            continue
        elif eik.is_trial(m1):
            plot_point(verts, m1, 0.75, 'yellow')
            if m1 not in ltrial:
                ltrial.append(m1)
        else:
            assert False
ltrial = np.array(ltrial)

D = []
for m in ltrial:
    dx = verts[np.unique(edges)] - verts[m]
    D.append(np.sqrt(np.sum(dx**2, axis=1)).min())

print(ltrial[np.argmin(D)])
