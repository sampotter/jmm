import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

import jmm

eps = 1e-5
xsrc = np.array([-7.5, -7.5, 1.25])
rfac = 0.5

mesh_data = jmm.Mesh3Data.from_off('../../examples/varying_s/room.off', 1e-1, False)
mesh_data.insert_vert(xsrc, eps)

mesh = jmm.Mesh3(mesh_data, eps=eps)

eik = jmm.Eik3(mesh, jmm.Sfunc.Constant)

grid = pv.UnstructuredGrid({vtk.VTK_TETRA: mesh.cells}, mesh.verts)

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(grid, show_edges=True)
