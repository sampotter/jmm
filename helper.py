import jmm
import meshio
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

################################################################################
# Set up some test mesh stuff

# mesh_path = 'room.vtu'
mesh_path = 'sphere.vtu'

print(f'- loading {mesh_path}')
mesh = meshio.read(mesh_path)
verts = mesh.points.astype(np.float64).copy()
cells = mesh.cells_dict['tetra'].astype(np.uintp).copy()
del mesh

mesh = jmm.Mesh3(verts, cells)
