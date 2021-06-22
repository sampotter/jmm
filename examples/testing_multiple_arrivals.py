import colorcet as cc
import logging
import json
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

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
src_index = l_int[3]
num_arrivals = 5
omega = 1000 # Hz

log.info('computing %d arrivals starting from index %d', num_arrivals, src_index)

field = PointSourceField(domain, src_index, omega)

ma = MultipleArrivals(domain, field, num_arrivals)
ma.traverse()

############################################################################

f = ma._fields[2]

poly = pv.PolyData(domain.mesh.verts)
poly['time'] = f.time
poly['scale'] = f._scale
poly['amplitude'] = 20*np.log10(np.clip(f.amplitude, 1e-3, 1))

_ = pvqt.BackgroundPlotter()
plot_mesh2(_, domain.mesh.get_surface_mesh(), opacity=0.25, color='white')
plot_field_BCs(_, f, opacity=0.4, r=domain.h, show_edges=True, color='magenta')
# _.add_mesh(poly, scalars='time', cmap=cc.cm.rainbow)
# _.add_mesh(poly, scalars='amplitude', cmap=cc.cm.fire)
_.add_mesh(poly, scalars='scale', cmap=cc.cm.gray)
