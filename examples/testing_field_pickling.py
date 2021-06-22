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

# field.solve()
l0 = 3064
while field.eik.peek() != l0:
    field.eik.step()

log.info('%s', field.parent_labels)

verts = field.domain.verts

mesh = field.domain.mesh
r = mesh.min_edge_length

plotter = pvqt.BackgroundPlotter()
# plotter.background_color = 'white'

plot_mesh2(plotter, field.domain.mesh.get_surface_mesh(), opacity=0.4,
           color='white', show_edges=False)

def plot_BCs(F, c, opacity=0.8):
    plot_field_BCs(plotter, F, color=c, r=0.75*r, opacity=opacity)
    if not isinstance(F, jmm.multiple_arrivals.PointSourceField):
        I = np.unique(F.bd_inds)
        plot_vector_field(
            plotter, verts[I], F.eik.t_in[I], 3*r, color=c, opacity=opacity)

plot_BCs(field, 'cyan')
plot_BCs(field.parent, 'orange')
plot_BCs(field.parent.parent, 'pink')
plot_BCs(field.parent.parent.parent, 'brown')
plot_BCs(field.parent.parent.parent.parent, 'purple')

# points = pv.PolyData(verts)
# points['_'] = 20*np.log10(np.maximum(1e-3, abs(field.parent.amplitude)))
# plotter.add_mesh(points, scalars='_', cmap=cc.cm.fire)
