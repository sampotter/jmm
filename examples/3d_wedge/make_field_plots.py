import colorcet as cc
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import pyvistaqt as pvqt
import vtk

from matplotlib.ticker import FuncFormatter

plt.ion()

from pathlib import Path

xsrc = [1, 1]

paths = [
    Path('n1.75_a0.01_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2'),
    Path('n1.75_a0.0031622776601683794_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2'),
    Path('n1.75_a0.001_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2'),
    Path('n1.75_a0.00031622776601683794_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2'),
    Path('n1.75_a0.0001_rfac0.25_phip0.7853981633974483_sp1.4142135623730951_w4_h2'),
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

def plot_slice(shape, xymin, h, order, path, origin_path, levels=21, cmap=None,
               vmin=None, vmax=None, diverging=True, **kwargs):
    assert order == 'row major'

    x, y = get_grid(shape, xymin, h, order)
    f = np.fromfile(path).reshape(shape)
    origin = np.fromfile(origin_path).reshape(shape)

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
        num_levels = levels
        levels = np.linspace(vmin, vmax, num_levels)
        levels[num_levels//2] = 0

    plt.contourf(x, y, f, levels=levels, cmap=cmap, **kwargs)
    fmt = lambda value, pos: f'{value:1.2e}'
    cbar = plt.colorbar(format=FuncFormatter(fmt))

    if (np.nanmin(origin) < 0.5 and 0.5 < np.nanmax(origin)):
        plt.contour(x, y, origin, levels=[0.5], colors='k', linestyles='--')

    return x, y, f, cbar

def get_h(path):
    V = np.fromfile(path/'verts.bin', np.float64).reshape(-1, 3)
    C = np.fromfile(path/'cells.bin', np.uintp).reshape(-1, 4)
    dV = V[C[:, 1:]] - V[C[:, 0]].reshape(-1, 1, 3)
    return np.mean(np.sqrt(np.sum(dV**2, axis=1)))

def get_N(path):
    return np.fromfile(path/'verts.bin', np.float64).reshape(-1, 3).shape[0]

def add_plot(path, field_path, origin_path):
    h = get_h(path)
    img_grid = read_img_grid(path/'img_grid.txt')
    x, y, f, cbar = plot_slice(*img_grid, path/field_path, path/origin_path, zorder=1)
    plt.plot([10, 0, 10], [0, 0, -10], linewidth=2, c='k', zorder=3)
    plt.scatter(*xsrc, s=3, c='k', marker='o', linewidth=2, zorder=3)
    plt.text(*xsrc, r'$x_{\mathrm{src}}$')
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.gca().set_aspect('equal')
    return cbar

plt.figure(figsize=(9.5, 12))
for i, path in enumerate(paths):
    plt.subplot(len(paths), 3, 3*i + 1)
    add_plot(path, 'slice_direct_E_T.bin', 'slice_direct_origin.bin')
    if i == 0: plt.title('Direct')
    if i == len(paths) - 1: plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.subplot(len(paths), 3, 3*i + 2)
    add_plot(path, 'slice_o_refl_E_T.bin', 'slice_o_refl_origin.bin')
    if i == 0: plt.title('Reflection ($o$-face)')
    if i == len(paths) - 1: plt.xlabel('$x$')
    plt.subplot(len(paths), 3, 3*i + 3)
    cbar = add_plot(path, 'slice_n_refl_E_T.bin', 'slice_n_refl_origin.bin')
    if i == 0: plt.title('Reflection ($n$-face)')
    if i == len(paths) - 1: plt.xlabel('$x$')
    h, N = get_h(path), get_N(path)
    cbar.set_label(f'$h = {h:.03g}$ ($N = {N:,}$)')
plt.tight_layout()
plt.savefig('error-field-plots.pdf')
plt.show()

############################################################################
# plotting total field...

path = paths[-1]

img_grid = read_img_grid(path/'img_grid.txt')
x, y = get_grid(*img_grid)
shape = img_grid[0]

k = 50

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
    plt.plot([4, 0, 4], [0, 0, -4], linewidth=2, c='k', zorder=2)
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)

Tmax = max(np.nanmax(Td), np.nanmax(To), np.nanmax(Tn))
Amax = max(np.nanmax(Ad), np.nanmax(Ao), np.nanmax(An))
umax = 3

m, n = 4, 3
plt.figure(figsize=(9, 10))
plt.subplot(m, n, 1)
make_plot(Td, (0, Tmax), cc.cm.gouldian)
plt.ylabel('$y$')
plt.title('$T$ (direct)')
plt.subplot(m, n, 2)
make_plot(To, (0, Tmax), cc.cm.gouldian)
plt.title('$T$ ($o$-face)')
plt.subplot(m, n, 3)
make_plot(Tn, (0, Tmax), cc.cm.gouldian)
plt.title('$T$ ($n$-face)')
plt.colorbar()
plt.subplot(m, n, 4)
make_plot(Ad, (0, 10), cc.cm.gouldian)
plt.title('Spreading factor (direct)')
plt.ylabel('$y$')
plt.subplot(m, n, 5)
make_plot(Ao, (0, Amax), cc.cm.gouldian)
plt.title('Spreading factor ($o$-face)')
plt.subplot(m, n, 6)
make_plot(An, (0, Amax), cc.cm.gouldian)
plt.title('Spreading factor ($n$-face)')
plt.colorbar()
plt.subplot(m, n, 7)
make_plot(np.real(ud), (-umax, umax), cc.cm.CET_D13)
plt.title('$u$ (direct)')
plt.ylabel('$y$')
plt.subplot(m, n, 8)
make_plot(np.real(uo), (-umax, umax), cc.cm.CET_D13)
plt.title('$u$ ($o$-face)')
plt.subplot(m, n, 9)
make_plot(np.real(un), (-umax, umax), cc.cm.CET_D13)
plt.title('$u$ ($n$-face)')
plt.colorbar()
plt.subplot(m, n, 10)
make_plot((Od >= 0.5).astype(float)*np.real(ud), (-umax, umax), cc.cm.CET_D13)
plt.title('$[\mathtt{origin} \geq 1/2] \cdot u$ (direct)')
plt.ylabel('$y$')
plt.subplot(m, n, 11)
make_plot((Oo >= 0.5).astype(float)*np.real(uo), (-umax, umax), cc.cm.CET_D13)
plt.title('$[\mathtt{origin} \geq 1/2] \cdot u$ ($o$-face)')
plt.subplot(m, n, 12)
make_plot(np.real(un), (-umax, umax), cc.cm.CET_D13)
plt.title('$D \cdot u$ ($n$-face)')
plt.colorbar()
plt.tight_layout()
plt.savefig('individual-fields.pdf')
plt.show()

plt.figure()
plt.imshow(Oo)
plt.show()

# u = ud + uo + un
u = ud*Od + uo*Oo + un
# u = ud*(Od >= 0.5).astype(float) + uo*(Oo >= 0.5).astype(float)
plt.figure(figsize=(9, 7))
make_plot(np.real(u), (-umax, umax), cc.cm.CET_D13)
plt.colorbar()
plt.tight_layout()
plt.ylabel('$y$')
plt.xlabel('$x$')
plt.savefig('total-field.pdf')
plt.tight_layout()
plt.show()

xfer = lambda x: (1 - x)*x**2 + x*(1 - (1 - x)**2)

plt.figure(figsize=(10, 6))
levels = np.linspace(-60, 0, 2*12 + 1)
for i, (A, O) in enumerate(zip([Ad, Ao, An], [Od, Oo, On])):
    plt.subplot(2, 3, i + 1)
    A_dB = 20*np.log10(np.maximum(1e-16, A/Amax))
    plt.contourf(x, y, A_dB, levels, cmap=cc.cm.gouldian)
    plt.gca().set_aspect('equal')
    plt.title(['Direct', r'Reflection ($o$-face)', r'Reflection ($n$-face)'][i])
    plt.plot([4, 0, 4], [0, 0, -4], linewidth=1, c='k', zorder=2)
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    if i == 0: plt.ylabel('$y$')
    O = xfer(xfer(xfer(O)))
    # O = (O >= 0.5).astype(np.float64)
    plt.subplot(2, 3, 3 + i + 1)
    AO_dB = np.maximum(-60, 20*np.log10(np.maximum(1e-16, O*A/Amax)))
    im = plt.contourf(x, y, AO_dB, levels, cmap=cc.cm.gouldian)
    plt.gca().set_aspect('equal')
    plt.plot([4, 0, 4], [0, 0, -4], linewidth=1, c='k', zorder=2)
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.xlabel('$x$')
    if i == 0: plt.ylabel('$y$')
plt.tight_layout()
plt.subplots_adjust(right=0.9)
cax = plt.gcf().add_axes([0.9, 0.09, 0.025, 0.85])
plt.colorbar(im, cax=cax)
plt.savefig('amp-contours.pdf')
plt.show()

plt.figure(figsize=(14, 4))
plt.tight_layout()
plt.savefig('amp-contours-weighted.png')
plt.show()
