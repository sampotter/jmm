import colorcet as cc
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

from scipy.spatial.transform import Rotation

from pathlib import Path

path = Path('n0.25_a0.01_rfac0.2_phip0.7853981633974483_sp1.4142135623730951_w4_h2')

verts = np.fromfile(path/'verts.bin', dtype=np.float64).reshape(-1, 3)
cells = np.fromfile(path/'cells.bin', dtype=np.uintp).reshape(-1, 4)
num_verts, num_cells = verts.shape[0], cells.shape[0]

eik = 'direct'

jet = np.fromfile(path/f'{eik}_jet.bin', dtype=np.float64).reshape(-1, 4)

origin = np.fromfile(path/f'{eik}_origin.bin', dtype=np.float64)
par_l = np.fromfile(path/f'{eik}_par_l.bin', dtype=np.uintp).reshape(-1, 3)
par_b = np.fromfile(path/f'{eik}_par_b.bin', dtype=np.float64).reshape(-1, 3)
accepted = np.fromfile(path/f'{eik}_accepted.bin', dtype=np.uintp)
has_bc = np.fromfile(path/f'{eik}_has_bc.bin', dtype=np.bool8)

T = jet[:, 0]
grad_T = jet[:, 1:4]
# hess_T = jet[:, 4:].reshape(-1, 3, 3)

hess_T_bb_cell = np.fromfile(path/f'{eik}_hess_bb.bin', dtype=np.float64)
hess_T_bb_cell = hess_T_bb.reshape(-1, 4, 3, 3)

def hess_T_bb_for_vert(l):
    H = []
    for i, L in enumerate(cells):
        if l in L:
            j = np.where(L == l)[0][0]
            H.append(hess_T_bb_cell[i][j])
    return H

hess_T_bb = np.array([
    np.array(hess_T_bb_for_vert(l)).mean(0)
    for l in range(num_verts)
])

jet_gt = np.fromfile(path/f'{eik}_jet_gt.bin', dtype=np.float64).reshape(-1, 13)
T_gt = jet_gt[:, 0]
grad_T_gt = jet_gt[:, 1:4]
hess_T_gt = jet_gt[:, 4:].reshape(-1, 3, 3)

edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

def get_level_set(verts, cells, values, level, f=None):
    cell_values = values[cells]
    mins = cell_values.min(1)
    maxs = cell_values.max(1)
    I = np.where((mins <= level) & (level <= maxs))[0]

    level_set_faces = []
    level_set_verts = []
    if f is not None:
        f_values = []

    j0 = 0 # current vertex index offset

    for i in cells[I]:
        face_verts = []
        if f is not None:
            f_face = []
        for e in edges:
            i0, i1 = i[e[0]], i[e[1]]
            v0, v1 = values[i0], values[i1]
            if min(v0, v1) <= level <= max(v0, v1):
                t = (level - v0)/(v1 - v0)
                x0, x1 = verts[i0], verts[i1]
                xt = x0 + t*(x1 - x0)
                face_verts.append(xt)
                if f is not None:
                    # simple linear interpolation... blech
                    f0, f1 = f[i0], f[i1]
                    ft = f0 + t*(f1 - f0)
                    f_face.append(ft)

        num_face_verts = len(face_verts)

        if num_face_verts == 3:
            faces = [[j0 + j for j in range(num_face_verts)]]
        else:
            p0 = sum(face_verts)/num_face_verts
            dp = np.array([p - p0 for p in face_verts]).T
            t = np.linalg.svd(dp, full_matrices=False)[0][:, :-1]
            x, y = t.T@dp
            J = np.argsort(np.arctan2(y, x))
            face_verts = [face_verts[j] for j in J] + [p0]
            if f is not None:
                f_face = [f_face[j] for j in J] + [sum(f_face)/num_face_verts]
            faces = []
            for j in range(num_face_verts):
                faces.append([
                    j0 + j,
                    j0 + ((j + 1) % num_face_verts),
                    j0 + num_face_verts
                ])

        level_set_faces.extend(faces)
        level_set_verts.extend(face_verts)
        if f is not None:
            f_values.extend(f_face)

        j0 += len(face_verts)

    level_set_verts = np.array(level_set_verts, dtype=np.float64)
    level_set_faces = np.array(level_set_faces, dtype=np.uintp)
    if f is not None:
        f_values = np.array(f_values)

    level_set_grid = pv.UnstructuredGrid(
        {vtk.VTK_TRIANGLE:level_set_faces},
        level_set_verts
    )

    if f is not None:
        return level_set_grid, f_values
    else:
        return level_set_grid

############################################################################
# PROPAGATING D2T

xsrc = np.array([1, -1, 0])

hess_T = np.empty((num_verts, 3, 3), dtype=np.float64)
hess_T[...] = np.nan

for l in accepted:
    if has_bc[l]:
        x = verts[l]
        r = x - xsrc
        tau = np.linalg.norm(r)
        hess_T[l] = (np.eye(r.size) - np.outer(r, r)/np.dot(r, r))/tau
    else:
        npar = np.isfinite(par_b[l]).sum()
        Ds, Qs, qs, ws = [], [], [], []
        for i in range(npar):
            Q, D = np.linalg.svd(hess_T[par_l[l][i]])[:2]
            if Q[:, -1]@grad_T[par_l[l][i]] < 0:
                Q[:, -1] *= -1
            Qs.append(Q)
            Ds.append(D)
            q = Rotation.from_matrix(Q).as_quat()
            qs.append(q)
            ws.append(par_b[l][i])
        qlam = sum(w*q for w, q in zip(ws, qs))
        qlam /= np.linalg.norm(qlam)
        # Qlam = Rotation.from_quat(qlam).as_matrix()
        Qlam = sum(w*Q for w, Q in zip(ws, Qs))
        Qlam = np.linalg.qr(Qlam)[0]
        Dlam = sum(w*D for w, D in zip(ws, Ds))
        xlam = sum(w*x for w, x in zip(ws, verts[par_l[l]]))
        Llam = np.linalg.norm(verts[l] - xlam)
        hess_T[l] = Qlam@np.diag(Dlam/(1 + Llam*Dlam))@Qlam.T

############################################################################
# PLOTTING

grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

def xfer(x):
    return x**2*(3 - 2*x)

h = 0.025

points = pv.PolyData(verts)
points['origin'] = origin # xfer(xfer(xfer(origin)))

shadow_boundary = get_level_set(verts, cells, origin, 0.5)

# Fhat = T
# F = T_gt
# f = abs(F - Fhat) # /abs(F)
f = hess_T_bb[:, 2, 2]
# f = hess_T_gt[:, 2, 2]
hplane, f_interp = get_level_set(verts, cells, verts[:, 2], 0, f)
hplane['f'] = f_interp

plotter = pvqt.BackgroundPlotter()
# plotter.background_color = 'white'
plotter.add_mesh(grid, show_edges=False, opacity=0.25)
# plotter.add_mesh(points, scalars='origin', cmap=cc.cm.fire)
# plotter.add_mesh(shadow_boundary, color='purple', opacity=0.95)
plotter.add_mesh(hplane, scalars='f', cmap=cc.cm.bmy,
                 show_edges=False, interpolate_before_map=True, clim=(0, 4.5))
