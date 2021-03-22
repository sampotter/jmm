# Note: can sometimes be helpful to run this section first before
# running the rest of the script if running from python-mode in Emacs

import colorcet as cc
import jmm
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

plotter = pvqt.BackgroundPlotter()

################################################################################
# parameters

verts_path = 'make_L_video/L_verts.bin'
cells_path = 'make_L_video/L_cells.bin'

lsrc = 0 # index of point source
l0 = 87 # index of first node to stop at


################################################################################

verts = np.fromfile(verts_path, dtype=np.float64)
verts = verts.reshape(verts.size//3, 3)

cells = np.fromfile(cells_path, dtype=np.uintp)
cells = cells.reshape(cells.size//4, 4)

mesh = jmm.Mesh3.from_verts_and_cells(verts, cells)

lsrc = 16

eik = jmm.Eik3(mesh)
eik.add_trial(lsrc, jmm.Jet3(0.0, np.nan, np.nan, np.nan))


l = 295

# l0 = 186
# l0 = 215
l0 = 160
# l0 = 319
while eik.peek() != l0:
    eik.step()


plotter = pvqt.BackgroundPlotter()

highlight_inds = [
    (lsrc, 'pink'),
    (l, 'magenta'),
    (l0, 'blue'),
]

h = eik.mesh.min_tetra_alt/2
sphere_radius = h/6

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
points = verts[uniq_inds]
faces = np.concatenate([
    3*np.ones((tris_on_front.shape[0], 1), dtype=tris_on_front.dtype),
    tris_on_front
], axis=1)
poly_data = pv.PolyData(points, faces)
poly_data.point_arrays['T'] = T

# Add it to the plotter and plot the values of the eikonal
plotter.add_mesh(poly_data, scalars='T', cmap=cc.cm.rainbow,
                 opacity=0.5, show_edges=True)

# Now, traverse each point and gradient, and add a colored arrow
# for ind, p, d in zip(uniq_inds, points, DT):
#     color = 'white'
#     sphere = pv.Sphere(sphere_radius, p)
#     plotter.add_mesh(sphere, color=color)
#     arrow = pv.Arrow(p, d, tip_length=0.1, scale=h)
#     plotter.add_mesh(arrow, color=color)

for ind, color in highlight_inds:
    p = verts[ind]
    sphere = pv.Sphere(1.25*sphere_radius, p)
    plotter.add_mesh(sphere, color=color)

state_to_color = {
    jmm.State.Valid.value: 'green',
    jmm.State.Trial.value: 'yellow',
    jmm.State.Far.value: 'white',
}

for l1 in mesh.vv(l):
    if eik.is_valid(l1):
        print('%d is VALID' % l1)
    plotter.add_mesh(pv.Sphere(1.1*sphere_radius, verts[l1]),
                     color=state_to_color[eik.state[l1]])

def plot_tri(L):
    plotter.add_mesh(
        pv.make_tri_mesh(verts, np.array(L).reshape(1, 3)),
        opacity=0.95, color='white', show_edges=True)
plot_tri([160, 186, 215])
