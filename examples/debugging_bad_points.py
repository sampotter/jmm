# Note: can sometimes be helpful to run this section first before
# running the rest of the script if running from python-mode in Emacs

import colorcet as cc
import jmm
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

plotter = pvqt.BackgroundPlotter()
# plotter = pv.Plotter()

################################################################################
# parameters

# verts_path = 'visualize_cutset/room_verts.bin'
# cells_path = 'visualize_cutset/room_cells.bin'
# surf_mesh_path = 'visualize_cutset/room.obj'

# lsrc = 0 # index of point source
# l = 3167
# l0 = 499
# l1 = 2138
# l2 = 4385
# l3 = None
# lbad = None

scene = 'L'

verts_path = 'visualize_cutset/%s_verts.bin' % scene
cells_path = 'visualize_cutset/%s_cells.bin' % scene

lsrc = 0 # index of point source
l = 93
l0 = 151
l1 = None
l2 = None
l3 = None
lbad = None

l_color = 'yellow'
l0_color = 'blue'
l1_color = 'black'
l2_color = 'blue'
l3_color = 'cyan'
lbad_color = 'green'

plot_surf_tris = False
plot_wavefront = False
plot_cutset = True
plot_ray_from_lsrc_to_l = False

################################################################################
# GEOMETRY SETUP

verts = np.fromfile(verts_path, dtype=np.float64)
verts = verts.reshape(verts.size//3, 3)

num_verts = verts.shape[0]

cells = np.fromfile(cells_path, dtype=np.uintp)
cells = cells.reshape(cells.size//4, 4)

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

mesh = jmm.Mesh3.from_verts_and_cells(verts, cells)

################################################################################
# SOLVE

eik = jmm.Eik3(mesh)
eik.add_trial(lsrc, jmm.Jet3(0.0, np.nan, np.nan, np.nan))

# l0 = eik.step()
# rounds = 3
# for _ in range(rounds):
#     while (eik.is_valid(l0) and verts[l0, 0] >= 1) or \
#           (eik.is_shadow(l0) and verts[l0, 0] < 1):
#         l0 = eik.step()
#     print(f'l0 = {l0}')
#     if _ < rounds - 1:
#         l0 = eik.step()

while eik.peek() != l0:
    eik.step()
# eik.step()

# eik.solve()

# if l is not None:
#     _ = verts[l] - verts[lsrc]
#     tau = np.linalg.norm(_)
#     Dtau = _/tau
#     del _
#     print(f'tau = {tau}')
#     print(f'T = {eik.jet[l][0]}')
#     print(f'T - tau = {tau - eik.jet[l][0]}')
#     print()
#     print(f'Dtau = {tuple(Dtau)}')
#     print(f'DT = ({eik.jet[l][1]}, {eik.jet[l][2]}, {eik.jet[l][3]})')
#     print(f'DT - Dtau = ({Dtau[0] - eik.jet[l][1]}, {Dtau[1] - eik.jet[l][2]}, {Dtau[2] - eik.jet[l][3]})')
#     print(f'angle(DT, Dtau) = {np.rad2deg(np.arccos(Dtau@[eik.jet[l][_] for _ in range(1, 4)]))} deg')

################################################################################
# HELPER FUNCTIONS FOR PLOTTING

def plot_tri(L):
    plotter.add_mesh(
        pv.make_tri_mesh(verts, np.array(L).reshape(1, 3)),
        opacity=0.95, color='white', show_edges=True)

def plot_point(points, l, scale=1.25, color='white', opacity=1):
    plotter.add_mesh(
        pv.Sphere(scale*sphere_radius, points[l]), color=color, opacity=opacity)

def plot_edge(l0, l1, color='white', opacity=1):
    x0, x1 = verts[l0], verts[l1]
    xm, xd = (x0 + x1)/2, x1 - x0
    r, h = 0.5*sphere_radius, np.linalg.norm(xd)
    plotter.add_mesh(
        pv.Cylinder(xm, xd, r, h), color=color, opacity=opacity)

def plot_cutset_edge(l0, l1):
    try:
        t = eik.shadow_cutset[min(l0, l1), max(l0, l1)].t
        plot_edge(l0, l1)
        xt = (1 - t)*verts[l0] + t*verts[l1]
        plotter.add_mesh(
            pv.Sphere(sphere_radius, xt), color='white', opacity=1)
    except:
        pass

def plot_jet(x, jet=None):
    if isinstance(jet, jmm.Jet3):
        d = np.array([jet.fx, jet.fy, jet.fz])
    else:
        d = np.array([jet[1], jet[2], jet[3]])
    assert(abs(1 - np.linalg.norm(d)) < 1e-13)
    plotter.add_mesh(
        pv.Arrow(x, d, scale=0.1),
        color='white', opacity=1)

################################################################################
# MAKE PLOTS

surf_mesh = pv.PolyData(
    verts,
    np.concatenate([3*np.ones((faces.shape[0], 1), dtype=np.uintp), faces], axis=1)
)
plotter.add_mesh(surf_mesh, color='white', opacity=0.25, show_edges=plot_surf_tris)

highlight_inds = {
    lsrc: 'pink',
    l0: l0_color,
}
if l is not None:
    highlight_inds[l] = l_color
if lbad is not None:
    highlight_inds[lbad] = lbad_color
if l1 is not None:
    highlight_inds[l1] = l1_color
if l2 is not None:
    highlight_inds[l2] = l2_color
if l3 is not None:
    highlight_inds[l3] = l3_color

h = eik.mesh.min_tetra_alt/2
sphere_radius = h/6

if plot_wavefront:

    # First, find the cells on the front---we initially take these to
    # be the cells which have exactly three VALID vertices.
    cells_on_front = \
        cells[(eik.state[cells] == jmm.State.Valid.value).sum(1) == 3]

    # Next, we want to filter out any cells on the front that contain
    # the same triangle. To do this, we count the corresponding
    # triangles on the front using a dictionary.
    tris_on_front = dict()
    for cell, state in zip(cells_on_front, eik.state[cells_on_front]):
        sorted_tri_inds = sorted(
            i for i, s in zip(cell, state) if s == jmm.State.Valid.value)
        key = tuple(sorted_tri_inds)
        if key not in tris_on_front:
            tris_on_front[key] = 1
        else:
            tris_on_front[key] += 1
    tris_on_front = np.array(
        [_ for _, count in tris_on_front.items() if count == 1],
        dtype=cells.dtype)

    # Now, find the unique indices so that we can look things up
    uniq_inds = np.unique(tris_on_front.ravel()).tolist()
    for i in range(tris_on_front.shape[0]):
        for j in range(3):
            tris_on_front[i, j] = uniq_inds.index(tris_on_front[i, j])

    # Pull out the jets and vertices corresponding to the unique
    # indices
    T = np.array([J[0] for J in eik.jet[uniq_inds]])
    DT = np.array([(J[1], J[2], J[3]) for J in eik.jet[uniq_inds]])

    # Put together the data for a PyVista PolyData instance containing
    # a triangle mesh representing the front
    if tris_on_front.size > 0:
        points = verts[uniq_inds]
        faces = np.concatenate([
            3*np.ones((tris_on_front.shape[0], 1), dtype=tris_on_front.dtype),
            tris_on_front
        ], axis=1)
        poly_data = pv.PolyData(points, faces)
        poly_data.point_arrays['T'] = T

        # Add it to the plotter and plot the values of the eikonal
        plotter.add_mesh(poly_data, scalars='T', cmap=cc.cm.rainbow,
                         opacity=1.0, show_edges=True)

    # Now, traverse each point and gradient, and add a colored arrow
    # for ind, p, d in zip(uniq_inds, points, DT):
    #     color = 'white'
    #     sphere = pv.Sphere(sphere_radius, p)
    #     plotter.add_mesh(sphere, color=color)
    #     arrow = pv.Arrow(p, d, tip_length=0.1, scale=h)
    #     plotter.add_mesh(arrow, color=color)

for ind, color in highlight_inds.items():
    p = verts[ind]
    sphere = pv.Sphere(1.5*sphere_radius, p)
    plotter.add_mesh(sphere, color=color)

if plot_ray_from_lsrc_to_l:
    if l is None:
        raise Exception('l must be defined to plot ray connecting lsrc and l')
    x, xsrc = verts[l], verts[lsrc]
    xm, xd = (x + xsrc)/2, x - xsrc
    r, h = 0.5*sphere_radius, np.linalg.norm(xd)
    plotter.add_mesh(
        pv.Cylinder(xm, xd, r, h), color=highlight_inds[l], opacity=1)

if plot_cutset:
    for (m0, m1), cutedge in eik.shadow_cutset.items():
        plot_point(verts, m0, color='purple' if eik.is_shadow(m0) else 'green')
        plot_point(verts, m1, color='purple' if eik.is_shadow(m1) else 'green')
        x0, x1 = verts[m0], verts[m1]
        xm, xd = (x0 + x1)/2, x1 - x0
        r, h = 0.5*sphere_radius, np.linalg.norm(xd)
        plotter.add_mesh(
            pv.Cylinder(xm, xd, r, h), color='white', opacity=1)
        t = cutedge.t
        xt = (1 - t)*verts[m0] + t*verts[m1]
        plotter.add_mesh(
            pv.Sphere(sphere_radius, xt), color='white', opacity=1)
        jet = cutedge.jet
        plot_jet(xt, jet)

################################################################################
# TMP

# xt = np.array([1, 1.6315787092401441, 0.5328554895202533])
# plotter.add_mesh(pv.Sphere(sphere_radius, xt), color='red')

# xb = eik.get_parent(l0).b@verts[eik.get_parent(l0).l]
# plotter.add_mesh(pv.Sphere(sphere_radius, xb), color='blue')

# x0, x1 = verts[135], verts[138]
# xm, xd = (x0 + x1)/2, x1 - x0
# r, h = 0.5*sphere_radius, np.linalg.norm(xd)
# plotter.add_mesh(
#     pv.Cylinder(xm, xd, r, h), color='cyan', opacity=1)

# for l_ in mesh.vv(72):
#     if not eik.is_valid(l_) and not eik.is_shadow(l_):
#         continue
#     c = 'cyan' if eik.is_valid(l_) else 'pink'
#     plot_point(verts, l_, color=c, scale=2)

# plot_edge(l0, l1)

# plot_point(verts, 124, color='red', scale=1.25)
# plot_point(verts, 11, color='magenta', scale=1.25)
# plot_point(verts, 98, color='orange', scale=1.25)

# for l_ in [48, 98, 130]:
#     c = 'yellow' if eik.is_valid(l_) else 'blue'
#     plot_point(verts, l_, color=c, scale=1.1)
#     plot_jet(verts[l_], eik.jet[l_])

for l_ in mesh.vv(l):
    state = jmm.State(eik.state[l_])
    print(l_, state)
    c = {
        jmm.State.Valid: 'cyan',
        jmm.State.Trial: 'yellow',
        jmm.State.Far: 'red',
        jmm.State.Shadow: 'brown'
    }[state]
    if l_ == l0:
        plot_point(verts, l_, color='white', scale=1.7)
    else:
        plot_point(verts, l_, color=c, scale=1.4)
