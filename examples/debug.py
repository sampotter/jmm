import jmm
import json
import numpy as np
import pickle
import pyvista as pv
import pyvistaqt as pvqt
import vtk

from plotting import plot_line, plot_wavefront

############################################################################
# parameters

should_plot = False

############################################################################
# load stuff from disk

grid = pv.read(f'sethian_shadow/Building.vtu')

with open('sethian_shadow/Building.json', 'r') as f:
    info = json.load(f)

cells = grid.cells.reshape(-1, 5)[:, 1:].astype(np.uintp)
points = grid.points.astype(np.float64)

num_dom_cells = info['num_dom_cells']
num_dom_points = info['num_dom_points']

dom_cells = cells[:num_dom_cells]
dom_points = points[:num_dom_points]

faces = set()
for C in dom_cells:
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
surf_mesh = pv.UnstructuredGrid({vtk.VTK_TRIANGLE: faces}, dom_points)

############################################################################
# set up problem to solve

mesh = jmm.Mesh3.from_verts_and_cells(dom_points, dom_cells)

with open('bcs.pickle', 'rb') as f:
    bcs = pickle.load(f)

eik = jmm.Eik3(mesh)

for le in bcs.keys():
    f, Df, x = bcs[le]
    bb = jmm.Bb31.from_3d_data(f, Df, x)
    eik.set_bde_bc(np.array(le), bb)

jets = dict()
for le in bcs.keys():
    f, Df, x = bcs[le]
    for i, l in enumerate(le):
        jets[l] = jmm.Jet3(f[i], *Df[i])

for l in jets.keys():
    tau = jets[l].f
    jet = jmm.Jet3.make_point_source(tau)
    eik.add_trial(l, jet)

############################################################################
# debugging

eik.solve()

if should_plot:
    p = pvqt.BackgroundPlotter()
    p.background_color = 'lightgray'
    p.add_mesh(surf_mesh, color='white', opacity=0.25, show_edges=False)
    for l0, l1 in bcs.keys():
        print(l0, l1)
        x0, x1 = points[l0], points[l1]
        plot_line(p, x0, x1, color='magenta')
    plot_wavefront(p, eik, opacity=0.5)
