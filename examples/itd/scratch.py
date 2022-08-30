import pyvista as pv
import numpy as np

grid = pv.read('HUTUB_pp2_in_cube_50k.off')

V = grid.points
F = grid.cells.reshape(-1, 4)[:, 1:]

v0 = V[F][:, 0, :]
v1 = V[F][:, 1, :]
v2 = V[F][:, 2, :]

t1 = v1 - v0
t2 = v2 - v0

A = np.sqrt(np.sum(np.cross(t1, t2)**2, axis=1))/2
