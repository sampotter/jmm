import colorcet as cc
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

from scipy.spatial.transform import Rotation

from pathlib import Path

path = Path('n1.75_a0.01_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2')
# path = Path('../../build/examples/3d_wedge')

verts = np.fromfile(path/'verts.bin', dtype=np.float64).reshape(-1, 3)
cells = np.fromfile(path/'cells.bin', dtype=np.uintp).reshape(-1, 4)
num_verts, num_cells = verts.shape[0], cells.shape[0]

# estimate h
dX = verts[cells][:, 1:, :] - verts[cells][:, 0, :].reshape(-1, 1, 3)
H = np.sqrt(np.sum(dX**2, axis=1))
h = H.mean()

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
t_in = np.fromfile(path/f'{eik}_t_in.bin', dtype=np.float64).reshape(-1, 3)
t_out = np.fromfile(path/f'{eik}_t_out.bin', dtype=np.float64).reshape(-1, 3)
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
# compute Hessian determinants

hess_det = np.array([np.linalg.det(_) for _ in hess_T])
hess_det_gt = np.array([np.linalg.det(_) for _ in hess_T_gt])

############################################################################
# checking logarithmic vs linear...

X, Y = np.meshgrid(np.linspace(0, 1), np.linspace(0, 1), indexing='xy')
x0, y0 = 1, 1
mask = X + Y > 1

Rinv = 1/np.sqrt((X + x0)**2 + (Y + y0)**2)
Rinv[mask] = np.nan

r0, r1, r2 = np.sqrt(x0**2 + y0**2), np.sqrt((x0 + 1)**2 + y0**2), np.sqrt(x0**2 + (y0 + 1)**2)
Rinv_lin = (1 - X - Y)*(1/r0) + X*(1/r1) + Y*(1/r2)
Rinv_lin[mask] = np.nan

# Rinv_linlog = np.exp((1 - X - Y)*np.log(1/r0) + X*np.log(1/r1) + Y*np.log(1/r2))
Rinv_linlog = ((1/r0)**(1 - X - Y))*((1/r1)**X)*((1/r2)**Y)
Rinv_linlog[mask] = np.nan

plt.figure(figsize=(12, 8))
plt.subplot(2, 3, 1)
plt.contourf(X, Y, Rinv)
plt.colorbar()
plt.gca().set_aspect('equal')
plt.subplot(2, 3, 2)
plt.contourf(X, Y, Rinv_lin)
plt.colorbar()
plt.gca().set_aspect('equal')
plt.subplot(2, 3, 3)
plt.contourf(X, Y, Rinv_linlog)
plt.colorbar()
plt.gca().set_aspect('equal')
plt.subplot(2, 3, 5)
plt.contourf(X, Y, abs(Rinv_lin - Rinv)/abs(Rinv))
plt.colorbar()
plt.gca().set_aspect('equal')
plt.subplot(2, 3, 6)
plt.contourf(X, Y, abs(Rinv_linlog - Rinv)/abs(Rinv))
plt.colorbar()
plt.gca().set_aspect('equal')
plt.show()

############################################################################
# computing the geometric part of the amplitude (spreading factors)

xsrc = np.array([-1, 1, 0])

# initialize the amplitude

def get_amp(g, H):
    amp = np.empty(num_verts, np.complex128)
    amp[...] = np.nan
    for l in range(num_verts):
        if np.isnan(par_b[l]).all():
            xhat = verts[l]
            r = np.linalg.norm(xhat - xsrc)
            amp[l] = 1/r
    for lhat in accepted:
        if np.isnan(amp[lhat]):
            b, l = par_b[lhat], par_l[lhat]
            npar = np.isfinite(b).sum()
            b, l = b[:npar], l[:npar]
            assert np.isfinite(amp[l]).all()
            amplam = np.product(amp[l]**b)
            xlam, xhat = b@verts[l], verts[lhat]
            L = np.linalg.norm(xhat - xlam)
            D2T = H[lhat]
            k1, k2 = -np.linalg.svd(D2T, compute_uv=False)[:2]
            amp[lhat] = amplam*np.exp(L*(k1 + k2)/2)
    return amp

amp = get_amp(grad_T, hess_T)
amp_gt = get_amp(grad_T_gt, hess_T_gt) # not actually the "true" amplitude!

############################################################################
# UTD coefs

# physical & geometric parameters

om = 10
c = 1
k = om/c
n = 7/8

n_o = np.array([0, -1, 0])
t_o = np.array([-1, 0, 0])
t_e = np.array([0, 0, 1])

############################################################################
# generic propagate function

def prop(BCs, init):
    value = np.empty(num_verts)
    value[...] = np.nan

    for l in BCs:
        init(value, l)

    queue = set()
    for l in np.where(np.isfinite(value))[0]:
        queue.add(l)
    while queue:
        l = queue.pop()
        init(value, l)
        b, L = par_b[l], par_l[l]
        npar = np.isfinite(b).sum()
        for i in range(npar):
            if np.isnan(value[L[i]]):
                queue.add(L[i])

    for l in accepted:
        if np.isfinite(value[l]):
            continue
        b, L = par_b[l], par_l[l]
        npar = np.isfinite(b).sum()
        if npar == 0:
            continue
        if np.isfinite(value[L[:npar]]).all():
            value[l] = b[:npar]@value[L[:npar]]

    return value

L_diff_edge = np.where(np.sqrt(np.sum(verts[:, :2]**2, axis=1)) < 1e-10)[0]

# compute s (distance along ray from diffracting edge)

def init_s(s, l):
    s[l] = T[l]

s = prop(L_diff_edge, init_s)
s = T - s

# compute rho_e

def init_rho_e(rho_e, l):
    t_in = grad_T[l]
    q_e = t_e - (t_in@t_e)*t_in
    q_e /= np.linalg.norm(q_e)
    rho_e[l] = 1/(q_e@hess_T[l]@q_e)

rho_e = prop(L_diff_edge, init_rho_e)

# compute rho_1

def init_rho_1(rho_1, l):
    t_in = grad_T[l]
    q_1 = n_o - (t_in@n_o)*t_in
    q_1 /= np.linalg.norm(q_1)
    rho_1[l] = 1/(q_1@hess_T[l]@q_1)

rho_1 = prop(L_diff_edge, init_rho_1)

# compute rho_2

def init_rho_2(rho_2, l):
    t_in = grad_T[l]
    q_1 = n_o - (t_in@n_o)*t_in
    q_1 /= np.linalg.norm(q_1)
    q_2 = np.cross(t_in, q_1)
    rho_2[l] = 1/(q_2@hess_T[l]@q_2)

rho_2 = prop(L_diff_edge, init_rho_2)

############################################################################
# compute groundtruth amplitude parameters



############################################################################
# compute D and D_gt for sound-hard boundary

# D = utd.D(1, n_o, t_o, t_e, k, n, t_in, t_out, hess_T, s)
# D_gt = utd.D(1, n_o, t_o, t_e, k, n, t_in_gt, t_out_gt, hess_T_gt, s_gt)

############################################################################
# PLOTTING

grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

def xfer(x):
    return x**2*(3 - 2*x)

points = pv.PolyData(verts)
points['origin'] = origin # xfer(xfer(xfer(origin)))

shadow_boundary = get_level_set(verts, cells, origin, 0.5)

i, j = 0, 0
# Fhat, F, emax = T, T_gt, h**2
# Fhat, F, emax = grad_T[:, i], grad_T_gt[:, i], h**2
Fhat, F, emax = hess_T[:, i, j], hess_T_gt[:, i, j], h

# f, clim, cmap = origin, (0, 1), cc.cm.CET_D1A
# f, clim, cmap = abs(F - Fhat)/np.maximum(1, abs(F)), (0, emax), cc.cm.gouldian
# f, clim, cmap = F, (-abs(F).max(), abs(F).max()), cc.cm.CET_D13
# f, clim, cmap = F, (-abs(F).max(), abs(F).max()), cc.cm.CET_D13
# f, clim, cmap = Fhat, (-abs(Fhat).max(), abs(Fhat).max()), cc.cm.CET_D13
# f, clim, cmap = 20*np.log10(np.real(amp)), (-60, 0), cc.cm.gouldian
# f, clim, cmap = 20*np.log10(np.real(amp_gt)), (-60, 0), cc.cm.gouldian
f, clim, cmap = abs(20*np.log10(np.real(amp_gt)) - 20*np.log10(np.real(amp))), None, cc.cm.gouldian
# f, clim, cmap = np.real(amp*np.exp(1j*om*T)), (-1, 1), cc.cm.CET_D13
# f, clim, cmap = np.real(amp_gt*np.exp(1j*om*T)), (-1, 1), cc.cm.CET_D13
# f, clim, cmap = np.log10(np.maximum(1e-16, abs(amp*np.exp(1j*om*T) - amp_gt*np.exp(1j*om*T_gt)))), (-4, 0), cc.cm.gouldian
# f, clim, cmap = np.log(np.maximum(1e-16, abs(hess_det_gt))), None, cc.cm.gouldian
# f, clim, cmap = hess_det_gt, None, cc.cm.gouldian

Z = [0]
hplanes = []
for z in Z:
    hplane, f_interp = get_level_set(verts, cells, verts[:, 2], z, f)
    hplane['f'] = f_interp
    hplanes.append(hplane)

points['T'] = T

points['t_in'] = t_in
glyph_t_in = points.glyph(orient='t_in', factor=0.1, geom=pv.Arrow())

points['t_out'] = t_out
glyph_t_out = points.glyph(orient='t_out', factor=0.1, geom=pv.Arrow())

points['s'] = rho_2
plotter = pvqt.BackgroundPlotter()
# plotter.background_color = 'white'
plotter.add_mesh(grid, show_edges=False, opacity=0.25)
# plotter.add_mesh(glyph_t_in, color='blue')
# plotter.add_mesh(glyph_t_out, color='red')
# plotter.add_mesh(points, scalars='s', clim=(0, 3), cmap=cc.cm.fire, nan_opacity=0)
# plotter.add_mesh(shadow_boundary, color='white', opacity=0.95)

for hplane in hplanes:
    plotter.add_mesh(hplane, scalars='f', cmap=cmap, clim=clim,
                     show_edges=False, interpolate_before_map=True)
    #                 clim=(0, emax))

# plt.figure()
# plt.hist(f[f < emax], bins=129)
# plt.show()
