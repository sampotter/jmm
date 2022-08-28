import colorcet as cc
import logging
import numpy as np
import pickle
import pyvista as pv
import pyvistaqt as pvqt
import time

from scipy.spatial import distance_matrix

import jmm.multiple_arrivals

from jmm.defs import State
from jmm.plot import *

np.seterr(all='raise')

PLOT_BAD_POINTS = True

logging.basicConfig(level=logging.INFO)

log = logging.getLogger('testing_field_pickling.py')

with open('field.pickle', 'rb') as f:
    field = pickle.load(f)

lhat = 973
l0 = 723
# l1 = 3064
while field.eik.peek() != l0:
    field.eik.step()
# field.eik.step()

log.info('%s', field.parent_labels)

verts = field.domain.verts

mesh = field.domain.mesh
surf_mesh = mesh.get_surface_mesh()
r = mesh.min_edge_length

I = np.arange(verts.shape[0])

plotter = pvqt.BackgroundPlotter(shape=(2, 2))

plotter.subplot(0, 0)
plot_mesh2(plotter, surf_mesh, opacity=0.25, color='white')
F = field.parent.parent.parent; F.solve()
values = 20*np.log10(np.maximum(1e-3, abs(F.amplitude[I])))
points = pv.PolyData(verts[I])
points['values'] = values
plotter.add_mesh(points, cmap=cc.cm.fire, clim=(-60, 0))

plotter.subplot(1, 0)
plot_mesh2(plotter, surf_mesh, opacity=0.25, color='white')
F = field.parent.parent; F.solve()
values = 20*np.log10(np.maximum(1e-3, abs(F.amplitude[I])))
points = pv.PolyData(verts[I])
points['values'] = values
plotter.add_mesh(points, cmap=cc.cm.fire, clim=(-60, 0))

plotter.subplot(0, 1)
plot_mesh2(plotter, surf_mesh, opacity=0.25, color='white')
F = field.parent; F.solve()
values = 20*np.log10(np.maximum(1e-3, abs(F.amplitude[I])))
points = pv.PolyData(verts[I])
points['values'] = values
plotter.add_mesh(points, cmap=cc.cm.fire, clim=(-60, 0))

# plotter.subplot(1, 1)
# plot_mesh2(plotter, surf_mesh, opacity=0.25, color='white')
# F = field; F.solve()
