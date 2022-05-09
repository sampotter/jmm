import colorcet as cc
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import pyvistaqt as pvqt
import vtk

plt.ion()

from pathlib import Path

xsrc = [1, 1]

paths = [
    Path('n1.75_a0.01_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2'),
    Path('n1.75_a0.0031622776601683794_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2'),
    Path('n1.75_a0.001_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2'),
]

def read_img_grid(path):
    with open(path, 'r') as f:
        # read in the shape
        s = f.readline()
        shape = tuple(int(_) for _ in s.split(':')[1].strip().split(','))

        # read in the bottom left corner of the image
        s = f.readline()
        xymin = tuple(float(_) for _ in s.split(':')[1].strip().split(','))

        # read in the value of h
        s = f.readline()
        h = float(s.split(':')[1].strip())

        # read in the order
        s = f.readline()
        order = s.split(':')[1].strip()

    return shape, xymin, h, order

def get_grid(shape, xymin, h, order):
    x0, y0 = xymin
    x1, y1 = x0 + shape[0]*h, y0 + shape[1]*h
    xgrid = np.linspace(x0, x1, shape[0])
    ygrid = np.linspace(y0, y1, shape[1])
    x, y = np.meshgrid(xgrid, ygrid, indexing='ij')
    return x, y

def plot_slice(shape, xymin, h, order, path, levels=21, cmap=None,
               vmin=None, vmax=None, diverging=True, **kwargs):
    assert order == 'row major'

    x, y = get_grid(shape, xymin, h, order)
    f = np.fromfile(path).reshape(shape)

    if diverging:
        if cmap is None:
            cmap = cc.cm.CET_D1A
        if vmax is None:
            vmax = np.nanmax(abs(f))
        if vmin is None:
            vmin = -vmax
    else:
        if cmap is None:
            cmap = cc.cm.CET_L1
        if vmax is None:
            vmax = np.nanmax(f)
        if vmin is None:
            vmin = np.nanmin(f)

    if isinstance(levels, int):
        levels = np.linspace(vmin, vmax, levels)

    plt.contourf(x, y, f, levels=levels, cmap=cmap, **kwargs)

    return x, y, f

def get_h(path):
    V = np.fromfile(path/'verts.bin', np.float64).reshape(-1, 3)
    C = np.fromfile(path/'cells.bin', np.uintp).reshape(-1, 4)
    dV = V[C[:, 1:]] - V[C[:, 0]].reshape(-1, 1, 3)
    return np.mean(np.sqrt(np.sum(dV**2, axis=1)))

def add_plot(path, field_path):
    h = get_h(path)
    img_grid = read_img_grid(path/'img_grid.txt')
    x, y, f = plot_slice(*img_grid, path/field_path, zorder=1)
    plt.colorbar()
    # plt.contour(x, y, f, [h**3/10, h**2/10], colors=['black'], linestyles=['--'])
    plt.plot([10, 0, 10], [0, 0, -10], linewidth=3, c='k', zorder=3)
    plt.plot([1, -10], [1, -10], linewidth=2, c='k', linestyle='--', zorder=3)
    plt.scatter(*xsrc, s=25, c='k', marker='o', linewidth=2, zorder=3)
    plt.text(*xsrc, r'$x_{\mathrm{src}}$')
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.gca().set_aspect('equal')

plt.figure(figsize=(12, 9))
for i, path in enumerate(paths):
    plt.subplot(len(paths), 3, 3*i + 1)
    add_plot(path, 'slice_direct_E_T.bin')
    plt.subplot(len(paths), 3, 3*i + 2)
    add_plot(path, 'slice_o_refl_E_T.bin')
    plt.subplot(len(paths), 3, 3*i + 3)
    add_plot(path, 'slice_n_refl_E_T.bin')
plt.tight_layout()
plt.savefig('error-field-plots.png')
plt.show()

############################################################################
# plotting total field...

path = paths[-1]

img_grid = read_img_grid(path/'img_grid.txt')
x, y = get_grid(*img_grid)
shape = img_grid[0]

k = 25

Td = np.fromfile(path/'slice_direct_T.bin').reshape(shape)
Ad = np.fromfile(path/'slice_direct_A.bin').reshape(shape)
Od = np.fromfile(path/'slice_direct_origin.bin').reshape(shape)
ud = Ad*np.exp(-1j*k*Td)

To = np.fromfile(path/'slice_o_refl_T.bin').reshape(shape)
Ao = np.fromfile(path/'slice_o_refl_A.bin').reshape(shape)
Oo = np.fromfile(path/'slice_o_refl_origin.bin').reshape(shape)
uo = Ao*np.exp(-1j*k*To)

Tn = np.fromfile(path/'slice_n_refl_T.bin').reshape(shape)
An = np.fromfile(path/'slice_n_refl_A.bin').reshape(shape)
On = np.fromfile(path/'slice_n_refl_origin.bin').reshape(shape)
un = An*np.exp(-1j*k*Tn)

def make_plot(field, clim, cmap):
    vmin, vmax = clim
    plt.imshow(np.flipud(field.T), vmin=vmin, vmax=vmax, extent=[-2, 2, -2, 2],
               cmap=cmap, zorder=1)
    plt.gca().set_aspect('equal')
    plt.plot([4, 0, 4], [0, 0, -4], linewidth=1, c='k', zorder=2)
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)

Tmax = max(np.nanmax(Td), np.nanmax(To), np.nanmax(Tn))
Amax = max(np.nanmax(Ad), np.nanmax(Ao), np.nanmax(An))
umax = Amax

plt.figure(figsize=(12, 10))
plt.subplot(3, 3, 1)
make_plot(Td, (0, Tmax), cc.cm.gouldian)
plt.subplot(3, 3, 2)
make_plot(To, (0, Tmax), cc.cm.gouldian)
plt.subplot(3, 3, 3)
make_plot(Tn, (0, Tmax), cc.cm.gouldian)
plt.colorbar()
plt.subplot(3, 3, 4)
make_plot(Ad, (0, Amax), cc.cm.gouldian)
plt.subplot(3, 3, 5)
make_plot(Ao, (0, Amax), cc.cm.gouldian)
plt.subplot(3, 3, 6)
make_plot(An, (0, Amax), cc.cm.gouldian)
plt.colorbar()
plt.subplot(3, 3, 7)
make_plot(np.real(ud), (-umax, umax), cc.cm.CET_D13)
plt.subplot(3, 3, 8)
make_plot(np.real(uo), (-umax, umax), cc.cm.CET_D13)
plt.subplot(3, 3, 9)
make_plot(np.real(un), (-umax, umax), cc.cm.CET_D13)
plt.colorbar()
plt.tight_layout()
plt.savefig('individual-fields.png')
plt.show()

u = ud + uo + un

plt.figure(figsize=(8, 8))
make_plot(np.real(u), (-umax, umax), cc.cm.CET_D13)
plt.tight_layout()
plt.savefig('total-field.png')
plt.show()

plt.figure(figsize=(14, 4))
for i, A in enumerate([Ad, Ao, An]):
    plt.subplot(1, 3, i + 1)
    plt.contourf(x, y, 20*np.log10(A/Amax), 11, cmap=cc.cm.gouldian)
    plt.gca().set_aspect('equal')
    plt.plot([4, 0, 4], [0, 0, -4], linewidth=1, c='k', zorder=2)
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.colorbar()
plt.tight_layout()
plt.savefig('amp-contours.png')
plt.show()
