import colorcet as cc
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

from scipy.spatial.transform import Rotation

from pathlib import Path

path = Path('n0.25_a0.001_rfac0.2_phip0.7853981633974483_sp1.4142135623730951_w4_h2')

verts = np.fromfile(path/'verts.bin', dtype=np.float64).reshape(-1, 3)
cells = np.fromfile(path/'cells.bin', dtype=np.uintp).reshape(-1, 4)
num_verts, num_cells = verts.shape[0], cells.shape[0]

eik = 'direct'

jet = np.fromfile(path/f'{eik}_jet.bin', dtype=np.float64).reshape(-1, 4)
T = jet[:, 0]
grad_T = jet[:, 1:4]
hess_T = np.fromfile(path/f'{eik}_hess.bin', dtype=np.float64).reshape(-1, 3, 3)

jet_gt = np.fromfile(path/f'{eik}_jet_gt.bin', dtype=np.float64).reshape(-1, 13)
T_gt = jet_gt[:, 0]
grad_T_gt = jet_gt[:, 1:4]
hess_T_gt = jet_gt[:, 4:].reshape(-1, 3, 3)

origin = np.fromfile(path/f'{eik}_origin.bin', dtype=np.float64)
par_l = np.fromfile(path/f'{eik}_par_l.bin', dtype=np.uintp).reshape(-1, 3)
par_b = np.fromfile(path/f'{eik}_par_b.bin', dtype=np.float64).reshape(-1, 3)
accepted = np.fromfile(path/f'{eik}_accepted.bin', dtype=np.uintp)
has_bc = np.fromfile(path/f'{eik}_has_bc.bin', dtype=np.bool8)

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
            if v0 != v1 and min(v0, v1) <= level <= max(v0, v1):
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
# estimate h

dX = verts[cells][:, 1:, :] - verts[cells][:, 0, :].reshape(-1, 1, 3)
H = np.sqrt(np.sum(dX**2, axis=1))
h = H.mean()

############################################################################
# compute Hessian determinants

hess_det = np.array([np.linalg.det(_) for _ in hess_T])
hess_det_gt = np.array([np.linalg.det(_) for _ in hess_T_gt])

############################################################################
# computing the geometric part of the amplitude

om = 5

xsrc = np.array([-1, 1, 0])

# initialize the amplitude

def get_amp(g, H):
    amp = np.empty(num_verts, np.complex128)
    amp[...] = np.nan
    for l in range(num_verts):
        if np.isnan(par_b[l]).all():
            xhat = verts[l]
            r = np.linalg.norm(xhat - xsrc)
            amp[l] = 1
    for lhat in accepted:
        if np.isnan(amp[lhat]):
            b, l = par_b[lhat], par_l[lhat]
            npar = np.isfinite(b).sum()
            b, l = b[:npar], l[:npar]
            assert np.isfinite(amp[l]).all()
            amplam = np.product(amp[l]**b)
            xlam, xhat = b@verts[l], verts[lhat]
            L = np.linalg.norm(xhat - xlam)
            DT, D2T = g[lhat], H[lhat]
            q1, q2 = np.linalg.svd(np.eye(3) - np.outer(DT, DT))[0][:, :2].T
            k1, k2 = -q1@D2T@q1, -q2@D2T@q2
            amp[lhat] = amplam*np.exp(L*(k1 + k2)/2)
    return amp

amp = get_amp(grad_T, hess_T)
amp_gt = get_amp(grad_T_gt, hess_T_gt) # not actually the "true" amplitude!

############################################################################
# PLOTTING

grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

def xfer(x):
    return x**2*(3 - 2*x)

points = pv.PolyData(verts)
points['origin'] = origin # xfer(xfer(xfer(origin)))

shadow_boundary = get_level_set(verts, cells, origin, 0.5)

i, j = 1, 0
# Fhat, F, emax = T, T_gt, h**2
Fhat, F, emax = grad_T[:, i], grad_T_gt[:, i], h**2
# Fhat, F, emax = hess_T[:, i, j], hess_T_gt[:, i, j], h

f, clim, cmap = origin, (0, 1), cc.cm.CET_D1A
# f, clim, cmap = abs(F - Fhat)/np.maximum(1, abs(F)), (0, emax), cc.cm.gouldian
# f, clim, cmap = F, (-abs(F).max(), abs(F).max()), cc.cm.CET_D13
# f, clim, cmap = Fhat, (-abs(Fhat).max(), abs(Fhat).max()), cc.cm.CET_D13
# f, clim, cmap = 20*np.log10(np.real(amp)), (-60, 0), cc.cm.gouldian
# f, clim, cmap = 20*np.log10(np.real(amp_gt)), (-60, 0), cc.cm.gouldian
# f, clim, cmap = abs(20*np.log10(np.real(amp_gt)) - 20*np.log10(np.real(amp))), None, cc.cm.gouldian
# f, clim, cmap = np.real(amp*np.exp(1j*om*T)), (-1, 1), cc.cm.CET_D13
# f, clim, cmap = np.real(amp_gt*np.exp(1j*om*T)), (-1, 1), cc.cm.CET_D13
# f, clim, cmap = np.log(np.maximum(1e-16, abs(hess_det_gt))), None, cc.cm.gouldian

hplane, f_interp = get_level_set(verts, cells, verts[:, 2], 0, f)
hplane['f'] = f_interp

plotter = pvqt.BackgroundPlotter()
plotter.background_color = 'white'
plotter.add_mesh(grid, show_edges=False, opacity=0.25)
# plotter.add_mesh(points, scalars='origin', cmap=cc.cm.fire)
# plotter.add_mesh(shadow_boundary, color='purple', opacity=0.95)
plotter.add_mesh(hplane, scalars='f', cmap=cmap, clim=clim,
                 show_edges=False, interpolate_before_map=True)
#                 clim=(0, emax))

# plt.figure()
# plt.hist(f[f < emax], bins=129)
# plt.show()
