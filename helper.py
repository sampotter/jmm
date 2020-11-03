import colorcet as cc
import matplotlib.pyplot as plt
import meshio
import numpy as np
import jmm

plt.ion()

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

################################################################################
# Set up some test Bernstein polynomial stuff

r = lambda x: np.linalg.norm(x)
Dr = lambda x: x/r(x)

X = np.eye(3)

f = np.array([r(x) for x in X])
Df = np.array([Dr(x) for x in X])

bb = jmm.Bb3Tri(f, Df, X)