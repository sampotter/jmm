import colorcet as cc
import jmm
import numpy as np
import os
import pyvista as pv
import pyvistaqt as pvqt
import skimage.morphology
import time

from PIL import Image

from skimage.measure import marching_cubes

# Load mesh, etc.

tetra_mesh_path = 'room/room.1.vtk'
tetra_mesh = pv.read(tetra_mesh_path)

surf_mesh_path = 'room.obj'
surf_mesh = pv.read(surf_mesh_path)

points = tetra_mesh.points.copy().astype(np.float64)
cells = tetra_mesh.cells.astype(np.uint64)
cells = cells.reshape(cells.size//5, 5)[:, 1:].copy()

points = np.array(points)
cells = np.array(cells)

num_points = points.shape[0]
num_cells = cells.shape[0]

print('%d points and %d cells' % (num_points, num_cells))

# Compute eikonal

mesh = jmm.Mesh3.from_verts_and_cells(points, cells)
eik = jmm.Eik3(mesh)
eik.add_trial(3019, jmm.Jet3(0, np.nan, np.nan, np.nan))
eik.solve()

# Transfer solution

pmin = points.min(0) - 0.5
pmax = points.max(0) + 0.5
h = 0.1
dim = np.ceil((pmax - pmin)/h).astype(np.intc)
grid = jmm.Grid3(dim, pmin, h)

T = eik.transfer_solution_to_grid(grid)

# Make movie

levels = np.linspace(0.75, 20, 30*24)[1:-1]

def get_level_set_poly_data(level):
    selem = np.ones((3, 3, 3), dtype=np.intc)
    mask = skimage.morphology.binary_erosion(np.isfinite(T), selem)
    V, F = marching_cubes(T, level=level, spacing=(h, h, h),
                          allow_degenerate=False, mask=mask)[:2]
    return pv.PolyData(
        V + pmin,
        np.concatenate([
            3*np.ones((F.shape[0], 1)), F], axis=1).astype(F.dtype))

plotter = pvqt.BackgroundPlotter()
plotter.set_position((-20, 35, -20))
plotter.set_viewup((0, 1, 0))
plotter.set_focus((0, 0, 0))


plotter.add_mesh(surf_mesh, color='white', opacity=0.3)
plotter.add_mesh(get_level_set_poly_data(levels[0]), color='white')

os.mkdir('frames')

i = 0

level = levels[i]
plotter.clear()
plotter.add_mesh(surf_mesh, color='white', opacity=0.3)
plotter.add_mesh(get_level_set_poly_data(level), color='white')
i += 1

Image.fromarray(plotter.image).save('frames/%03d.png' % i)




print('finished')
