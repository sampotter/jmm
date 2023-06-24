import sys; sys.path.insert(-1, '../../wrappers/python')

import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

import jmm

verbose = True

c = 340.3
xsrc = np.array([-7.5, -7.5, 1.25])
rfac = 0.5

eps = 1e-5
off_path = '../../examples/varying_s/room.off'

mesh_data = jmm.Mesh3Data.from_off(off_path, 1e-1, verbose)
mesh_data.insert_vert(xsrc, eps)

mesh = jmm.Mesh3(mesh_data, eps=eps)

hh = jmm.Eik3hh.new_with_pt_src(mesh, c, rfac, xsrc)

root_branch = hh.get_root_branch()
root_branch.solve(verbose)

refl_inds = root_branch.get_visible_refls()

for refl_ind in refl_inds:
    refl_branch = root_branch.add_refl(refl_ind)
    refl_branch.solve(verbose)

# grid = pv.UnstructuredGrid({vtk.VTK_TETRA: mesh.cells}, mesh.verts)
grid = pv.read(off_path)

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(grid, show_edges=True)
