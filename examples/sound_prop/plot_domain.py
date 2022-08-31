import numpy as np
import pyvista as pv

import vtk

from pathlib import Path

if __name__ == '__main__':
    path = Path('../../build/examples/simple_building')
    verts = np.fromfile(path/'verts.bin', dtype=np.float64).reshape(-1, 3)
    cells = np.fromfile(path/'cells.bin', dtype=np.uintp).reshape(-1, 4)
    grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

    centroids = grid.points[cells].mean(1)
    mask = np.where(centroids[:, 2] < 2)[0]
    subgrid = grid.extract_cells(mask)

    plotter = pv.Plotter(off_screen=True, polygon_smoothing=True)
    plotter.add_mesh(grid, color='white', show_edges=True)
    plotter.camera.focal_point = (0, 0, 0)
    plotter.camera.up = (0, 0, 1)
    plotter.camera.position = (-50, 25, 30)
    plotter.camera.zoom(2)
    plotter.window_size = (1280, 720)
    plotter.screenshot('simple-building.png', transparent_background=True)
    plotter.close()

    plotter = pv.Plotter(off_screen=True, polygon_smoothing=True)
    plotter.add_mesh(subgrid, color='white', show_edges=True)
    plotter.camera.focal_point = (0, 0, 0)
    plotter.camera.up = (0, 1, 0)
    plotter.camera.position = (0, 0, 40)
    plotter.camera.zoom(1)
    plotter.window_size = (1024, 1024)
    plotter.screenshot('simple-building-cutaway.png', transparent_background=True)
    plotter.close()
