import colorcet as cc
import logging
import pickle
import pyvista as pv
import pyvistaqt as pvqt
import time

from scipy.spatial import distance_matrix

import jmm.multiple_arrivals

from jmm.defs import State
from jmm.plot import *

PLOT_BAD_POINTS = True

logging.basicConfig(level=logging.INFO)

log = logging.getLogger('testing_field_pickling.py')

with open('field.pickle', 'rb') as f:
    field = pickle.load(f)

eik = field.extended_eik

# eik.solve()
# assert False

l0 = 2450
l = 9011
while eik.peek() != l0:
    eik.step()
eik.step()

mesh = field.domain.extended_mesh
verts = mesh.verts
cells = mesh.cells
bd_inds = np.unique(field.bd_inds)

r = mesh.min_edge_length

surf_mesh = mesh.get_surface_mesh()

plotter = pvqt.BackgroundPlotter()
plotter.background_color = 'white'

plot_mesh2(plotter, surf_mesh, color='white', opacity=0.25)
plot_eik3(plotter, eik, cmap=cc.cm.bmw, opacity=1)

points = pv.PolyData(verts[np.unique(field.bd_inds)])
plotter.add_mesh(points.glyph(geom=pv.Sphere(r/2)), color='grey',
                 show_scalar_bar=False)

plot_point(plotter, verts[l0], r, color='yellow')
plot_point(plotter, verts[l], r, color='blue')

assert False

valid_mask = eik.state == 2
valid_points = pv.PolyData(verts[valid_mask])
valid_points.vectors = eik.grad_T[valid_mask]
plotter.add_mesh(valid_points.glyph(orient=True, geom=pv.Arrow(scale=np.sqrt(r))),
                 color='grey', show_scalar_bar=False)

for m in mesh.vv(l):
    plot_point(plotter, verts[m], r/4, color='white')

if PLOT_BAD_POINTS:
    if isinstance(field, jmm.multiple_arrivals.PointSourceField):
        L = np.where(eik.state == State.Trial.value)[0]
        M = L[np.argsort(np.sqrt(np.sum((verts[1] - verts[L])**2, axis=1)))]
    elif isinstance(field, jmm.multiple_arrivals.DiffractedField):
        L = np.where(eik.state == State.Trial.value)[0]
        min_bd_dist = distance_matrix(verts[L], verts[bd_inds]).min(1)
        M = L[np.argsort(min_bd_dist)]
    else:
        assert False

    num_bad_points = 5
    M = M[:num_bad_points]

    for m in M:
        plot_point(plotter, verts[m], (2/3)*r, color='cyan')
    print(M)
