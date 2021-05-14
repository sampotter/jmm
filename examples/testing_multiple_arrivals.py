import logging
import json
import numpy as np
import pyvista as pv

from jmm.mesh import Mesh3
from jmm.multiple_arrivals import PointSourceTree, MultipleArrivals

logging.basicConfig(level=logging.INFO)

log = logging.getLogger('testing_multiple_arrivals.py')

src_index = 0
num_arrivals = 5

surf_grid = pv.read(f'sethian_shadow/Building.obj')
grid = pv.read(f'sethian_shadow/Building.vtu')
with open('sethian_shadow/Building.json', 'r') as f:
    info = json.load(f)

extended_cells = grid.cells.reshape(-1, 5)[:, 1:].astype(np.uintp)
extended_points = grid.points.astype(np.float64)

num_cells = info['num_dom_cells']
cells = extended_cells[:num_cells]

num_points = info['num_dom_points']
points = extended_points[:num_points]

log.info('domain has %d points and %d cells', num_points, num_cells)
log.info('extended domain has %d points and %d cells', extended_points.shape[0],
         extended_cells.shape[0])

mesh = Mesh3.from_verts_and_cells(points, cells)
extended_mesh = Mesh3.from_verts_and_cells(extended_points, extended_cells)
for le in mesh.get_diff_edges():
    extended_mesh.set_boundary_edge(*le, True)

tree = PointSourceTree(mesh, extended_mesh, src_index)

ma = MultipleArrivals(tree, num_arrivals)
