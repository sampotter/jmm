import colorcet as cc
import logging
import json
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import sys

from jmm.mesh import Mesh3
from jmm.multiple_arrivals import Domain, Field, PointSourceField, MultipleArrivals
from jmm.plot import *

np.seterr(all='raise')

logging.basicConfig(level=logging.INFO)
log = logging.getLogger('testing_multiple_arrivals.py')

surf_grid = pv.read(f'sethian_shadow/Building.obj')
grid = pv.read(f'sethian_shadow/Building.vtu')
with open('sethian_shadow/Building.json', 'r') as f:
    info = json.load(f)

extended_cells = grid.cells.reshape(-1, 5)[:, 1:].astype(np.uintp)
extended_verts = grid.points.astype(np.float64)

num_dom_cells = info['num_dom_cells']
num_dom_verts = info['num_dom_points']

verts = extended_verts[:num_dom_verts]
cells = extended_cells[:num_dom_cells]

domain = Domain(verts, cells, refl_coef=0.7)

l_int = np.array([_ for _ in range(num_dom_verts) if not domain.mesh.bdv(_)])
src_index = l_int[int(sys.argv[1])]
num_arrivals = 5
omega = 3000 # Hz

log.info('computing %d arrivals starting from index %d (%dth interior vert)', num_arrivals, src_index, int(sys.argv[1]))

field = PointSourceField(domain, src_index, omega, r=1.5)

ma = MultipleArrivals(domain, field, num_arrivals)
ma.traverse(max_arrivals=30)

F = ma._fields[19]

mesh = F.domain.mesh
surf_mesh = mesh.get_surface_mesh()
verts = F.domain.mesh.verts
r = mesh.min_edge_length

P = pvqt.BackgroundPlotter()
plot_field_BCs(P, ma._fields[0], r=r, color='cyan')
plot_mesh2(P, surf_mesh, color='white', opacity=0.25)
plot_field_BCs(P, F, r=r/2, color='cyan')

values = 20*np.log10(np.maximum(Field.minimum_magnitude, abs(F.amplitude)))
plot_scalar_field(P, verts, values, clim=(-60, 0), cmap=cc.cm.fire)
