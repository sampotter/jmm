import colorcet as cc
import json
import meshplex
import numpy as np
import os
import pyvista as pv
import pyvistaqt as pvqt
import skimage.morphology
import time
import vtk

from PIL import Image

from skimage.measure import marching_cubes

from jmm.grid import Grid3
from jmm.multiple_arrivals import Domain, PointSourceField, ReflectedField, \
    DiffractedField
from jmm.plot import *

class LevelSetPlotter:
    def __init__(self, field, pmin, pmax, grid_h):
        assert field.is_solved
        self.pmin = pmin
        self.pmax = max
        self.dim = np.ceil((pmax - pmin)/grid_h).astype(np.intc)
        self.grid_h = grid_h
        self.grid = Grid3(self.dim, self.pmin, self.grid_h)
        self.field = field
        self.T = self.field.eik.transfer_solution_to_grid(self.grid)

    def plot_level_set(self, plotter, level, **kwargs):
        mask = skimage.morphology.binary_erosion(
            np.isfinite(self.T),
            selem=np.ones((3, 3, 3), dtype=np.intc)
        )
        h = self.field.h
        verts, faces = marching_cubes(
            self.T,
            level=level,
            spacing=(h, h, h),
            allow_degenerate=False,
            mask=mask
        )[:2]
        grid = pv.UnstructuredGrid({vtk.VTK_TRIANGLE: faces}, verts + self.pmin)
        plotter.add_mesh(grid, **kwargs)

if __name__ == '__main__':

    # Load mesh, etc.

    surf_grid = pv.read(f'../sethian_shadow/Building.obj')

    grid = pv.read(f'../sethian_shadow/Building.vtu')
    with open('../sethian_shadow/Building.json', 'r') as f:
        info = json.load(f)

    num_points = info['num_dom_points']
    points = grid.points[:num_points].astype(np.float64)

    num_cells = info['num_dom_cells']
    cells = grid.cells.reshape(-1, 5)[:num_cells, 1:].astype(np.uintp)

    print('%d points and %d cells' % (num_points, num_cells))

    # Compute eikonal

    domain = Domain(points, cells, 1.0)

    omega = 3000

    l_int = np.array([_ for _ in range(num_points) if not domain.mesh.bdv(_)])
    src_index = l_int[0]

    root_field = PointSourceField(domain, src_index, omega, 1.5)
    root_field.solve()

    pmin = points.min(0) - 0.5
    pmax = points.max(0) + 0.5
    grid_h = 0.1
    lsp = LevelSetPlotter(root_field, pmin, pmax, grid_h)

    # Make movie

    levels = np.linspace(0.75, 20, 30*24)[1:-1]

    plotter = pvqt.BackgroundPlotter()
    plotter.set_position((-20, 35, -20))
    plotter.set_viewup((0, 1, 0))
    plotter.set_focus((0, 0, 0))

    surf_mesh = domain.mesh.get_surface_mesh()
    plot_mesh2(plotter, surf_mesh, color='white', opacity=0.5)

    lsp.plot_level_set(plotter, levels[100], color='white')

    # os.mkdir('frames')

    # i = 0

    # level = levels[i]
    # plotter.clear()
    # plotter.add_mesh(surf_mesh, color='white', opacity=0.3)
    # plotter.add_mesh(get_level_set_poly_data(level), color='white')
    # i += 1

    # Image.fromarray(plotter.image).save('frames/%03d.png' % i)




    # print('finished')
