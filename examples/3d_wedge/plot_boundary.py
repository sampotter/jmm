import colorcet as cc
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

from pathlib import Path

path = Path('n0.25_a0.01_rfac0.2_phip0.7853981633974483_sp1.4142135623730951_w4_h2')

verts = np.fromfile(path/'verts.bin', dtype=np.float64).reshape(-1, 3)
cells = np.fromfile(path/'cells.bin', dtype=np.uintp).reshape(-1, 4)

origin = np.fromfile(path/'o_refl_origin.bin', dtype=np.float64)

edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

def get_level_set(verts, cells, values, level):
    cell_values = values[cells]
    mins = cell_values.min(1)
    maxs = cell_values.max(1)
    I = np.where((mins <= level) & (level <= maxs))[0]

    level_set_faces = []
    level_set_verts = []

    j0 = 0 # current vertex index offset

    for i in cells[I]:
        face_verts = []
        for e in edges:
            i0, i1 = i[e[0]], i[e[1]]
            v0, v1 = values[i0], values[i1]
            if min(v0, v1) <= level <= max(v0, v1):
                t = (level - v0)/(v1 - v0)
                x0, x1 = verts[i0], verts[i1]
                xt = x0 + t*(x1 - x0)
                face_verts.append(xt)

        num_face_verts = len(face_verts)

        if num_face_verts == 3:
            faces = [[j0 + j for j in range(num_face_verts)]]
        elif num_face_verts == 4:
            p0 = sum(face_verts)/num_face_verts
            dp = np.array([p - p0 for p in face_verts]).T
            t = np.linalg.svd(dp, full_matrices=False)[0][:, :-1]
            x, y = t.T@dp
            J = np.argsort(np.arctan2(y, x))
            face_verts = [face_verts[j] for j in J] + [p0]
            faces = [[j0,     j0 + 1, j0 + 4],
                    [j0 + 1, j0 + 2, j0 + 4],
                    [j0 + 2, j0 + 3, j0 + 4],
                    [j0 + 3, j0,     j0 + 4]]
        else:
            raise RuntimeError(f'???: num_face_verts == {num_face_verts}')

        level_set_faces.extend(faces)
        level_set_verts.extend(face_verts)

        j0 += len(face_verts)

    level_set_verts = np.array(level_set_verts, dtype=np.float64)
    level_set_faces = np.array(level_set_faces, dtype=np.uintp)

    return pv.UnstructuredGrid({vtk.VTK_TRIANGLE:level_set_faces},level_set_verts)

grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

points = pv.PolyData(verts)
points['origin'] = origin

shadow_boundary = get_level_set(verts, cells, origin, 0.5)

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(grid, show_edges=True, opacity=0.25)
# plotter.add_mesh(points, scalars='origin', cmap=cc.cm.fire)
plotter.add_mesh(shadow_boundary, color='purple', opacity=0.95)
