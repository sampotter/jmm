import pyvista as pv
import pyvistaqt as pvqt
import numpy as np
import vtk

from pathlib import Path

itd_build_dir = Path('../../builddir/examples/itd')

grid = pv.read(itd_build_dir/'HUTUB_pp2_in_cube_50k.off')

V = grid.points
F = grid.cells.reshape(-1, 4)[:, 1:]

V0 = V[F][:, 0, :]
V1 = V[F][:, 1, :]
V2 = V[F][:, 2, :]

T1 = V1 - V0
T2 = V2 - V0

A = np.sqrt(np.sum(np.cross(T1, T2)**2, axis=1))/2

V_tet = np.fromfile(itd_build_dir/'verts.bin').reshape(-1, 3)
C_tet = np.fromfile(itd_build_dir/'cells.bin', dtype=np.uint64).reshape(-1, 4)

grid_tet = pv.UnstructuredGrid({vtk.VTK_TETRA: C_tet}, V_tet)

P_tet = grid_tet.cell_centers().points

cell_ind = (P_tet[:, 1] > 400).nonzero()[0]
subgrid = grid_tet.extract_cells(cell_ind)

bbox = pv.Box((-1e3, 1e3, -1e3, 1e3, -1e3, 1e3))

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
plotter.add_mesh(bbox, 'r', 'wireframe')
plotter.add_legend([[' Input Mesh ', 'r'],
                    [' Tessellated Mesh ', 'black']])
plotter.show()
