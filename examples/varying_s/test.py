import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

norm = np.linalg.norm

V = np.fromfile('verts.bin', dtype=np.float64).reshape(-1, 3)
C = np.fromfile('cells.bin', dtype=np.uintp).reshape(-1, 4)

r = 0.1
scale = 0.5

xsrc = np.array([5, -5, 1])

lhat, l0, l1, l2 = 1847, 70, 1584, 0
x = V[lhat]
x0 = V[l0]
x1 = V[l1]
x2 = V[l2]

jet = np.fromfile('jet_T.bin', dtype=np.float64).reshape(-1, 4)

T = jet[:-1, 0]
DT = jet[:-1, 1:]

surf = pv.read('room.off')
grid = pv.UnstructuredGrid({vtk.VTK_TETRA: C}, V)

points = pv.PolyData(V[np.isfinite(T)])
points['T'] = T[np.isfinite(T)]
points['DT'] = DT[np.isfinite(T)]
points['scale'] = scale*np.ones(np.isfinite(T).sum())

tangents = points.glyph(scale='scale', orient='DT')

plotter = pvqt.BackgroundPlotter()
def plot_sphere(x, c):
    plotter.add_mesh(pv.Sphere(r, x), color=c)
def plot_arrow(x, d, c):
    plotter.add_mesh(pv.Arrow(x, d, scale=scale), color=c)
plotter.background_color = 'white'
plot_sphere(xsrc, 'red')
plotter.add_mesh(surf, opacity=0.1)
# plotter.add_mesh(points, scalars='T', cmap=cc.cm.gray_r)
plotter.add_mesh(tangents, scalars='T', cmap=cc.cm.colorwheel, show_scalar_bar=False)
#  # update:
# plot_sphere(x, 'yellow')
# plot_sphere(x0, 'cyan')
# plot_sphere(x1, 'blue')
# plot_sphere(x2, 'teal')
# plotter.add_mesh(pv.UnstructuredGrid({vtk.VTK_TRIANGLE: np.array([[0, 1, 2]], dtype=int)}, np.array([x0, x1, x2])), color='blue', opacity=0.5)
# # solution:
# plot_sphere([4.8418079764632918, -4.6681211690261026, 0.99000090880374736], 'green')
# plot_arrow(x, [-0.43162481465055197, 0.90106296411845488, 0.042255817018970993], 'yellow')
# # ...
# plotter.enable_point_picking()
plotter.show()
