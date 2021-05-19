import logging
import pickle
import pyvista as pv
import pyvistaqt as pvqt

from jmm.plot import *

logging.basicConfig(level=logging.INFO)

with open('field.pickle', 'rb') as f:
    field = pickle.load(f)
field.solve()

mesh = field.domain.mesh

surf_mesh = mesh.get_surface_mesh()

plotter = pvqt.BackgroundPlotter()
plotter.background_color = 'white'

plot_mesh2(plotter, surf_mesh, color='white', opacity=0.25)
plot_shadow(plotter, field, color='purple', opacity=0.8, show_edges=True)
plot_point(plotter, mesh.verts[field.bd_inds[0]], mesh.min_edge_length, color='red')
