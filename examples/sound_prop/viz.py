import colorcet as cc
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import vtk

plt.ion()

from scipy.spatial.transform import Rotation

from pathlib import Path

def get_level_set(verts, cells, values, level, f=None):
    cell_values = values[cells]
    mins = cell_values.min(1)
    maxs = cell_values.max(1)
    I = np.where((mins <= level) & (level <= maxs))[0]

    level_set_faces = []
    level_set_verts = []
    if f is not None:
        f_values = []

    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

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
        if num_face_verts == 0:
            continue

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

    if level_set_verts.size == 0:
        return

    level_set_grid = pv.UnstructuredGrid(
        {vtk.VTK_TRIANGLE:level_set_faces},
        level_set_verts
    )

    if f is not None:
        return level_set_grid, f_values
    else:
        return level_set_grid


def draw_boundary(linewidth=1.75):
    xbd = [-6.5, -10,   -10, -2,  -2,    -5, # lower left room
           -5, # portal
           -2, -2, -1.75, -1.75, 10, 10,  # lower right room
           2.5, # portal
           2.5, 10, 10, -6.75, -6.75, 1, # upper right room
           1, # portal
           -1.75, -1.75, # lower right room
           -2, # portal
           -2, -7, -7, -10, -10, -6.5, # left room
           -6.5, # portal
           ]
    ybd = [-4.25, -4.25, -10, -10, -4.25, -4.25, # lower left room
           -4, # portal
           -4, 0.25, 0.25, -10, -10, 5, # lower right room
           5, # portal
           5.25, 5.25, 10, 10, 5.25, 5.25, # upper right room
           5, # portal
           5, 2, # lower right room
           2, # portal
           5, 5, 10, 10, -4, -4, # left room
           -4.25, # portal
           ]
    plt.plot(xbd, ybd, linewidth=linewidth, c='k')
    xcol = np.array([0, 0.5, 0.5, 0, 0])
    ycol = np.array([0, 0, 0.5, 0.5, 0])
    for xoff, yoff in [(-0.5, 3), (-0.5, -1), (-0.5, -5), (-0.5, -8.45),
                       (3.25, -8.45), (7, -8.45),
                       (3.25, 3), (7, 3)]:
        plt.plot(xcol + xoff, ycol + yoff, linewidth=linewidth, c='k')

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
    xgrid = np.linspace(x0, x1, shape[0], endpoint=False)
    ygrid = np.linspace(y0, y1, shape[1], endpoint=False)
    x, y = np.meshgrid(xgrid, ygrid, indexing='ij')
    return x, y

def xfer(x):
    return 3*x**2 - 2*x**3

if __name__ == '__main__':
    path = Path('../../build/examples/simple_building')
    # name = 'direct'
    name = 'refl081'

    if name == 'direct':
        plotname = name.capitalize()
    elif 'refl' in name:
        plotname = f'Refl. #{int(name[4:])}'

    verts = np.fromfile(path/'verts.bin', dtype=np.float64).reshape(-1, 3)
    cells = np.fromfile(path/'cells.bin', dtype=np.uintp).reshape(-1, 4)
    num_verts, num_cells = verts.shape[0], cells.shape[0]

    jet = np.fromfile(path/f'{name}_jet.bin', dtype=np.float64).reshape(-1, 4)
    jet[~np.isfinite(jet)] = np.nan
    T = jet[:, 0]
    grad_T = jet[:, 1:4]

    img_grid = read_img_grid(path/'img_grid.txt')
    x, y = get_grid(*img_grid)
    shape = img_grid[0]
    h = img_grid[2]

    freq = 1000
    c = 340

    T = np.fromfile(path/f'{name}_T_slice.bin').reshape(shape)
    spread = np.fromfile(path/f'{name}_spread_slice.bin').reshape(shape)
    spread_dB = 20*np.log10(np.maximum(1e-16, spread))
    origin = np.fromfile(path/f'{name}_origin_slice.bin').reshape(shape)
    u = xfer(xfer(origin))*spread*np.exp(2*np.pi*1j*freq*T/c)

    margin = 0.5

    figsize = (5, 4)

    plt.figure(figsize=figsize)
    plt.contourf(x, y, T/c, levels=21, cmap=cc.cm.gouldian)
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.xlim(x.min() - margin, x.max() + margin)
    plt.ylim(y.min() - margin, y.max() + margin)
    plt.title(f'$T$ [s] ({plotname})')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    draw_boundary()
    plt.show()
    plt.tight_layout()
    plt.savefig(f'{name}_T.pdf')

    plt.figure(figsize=figsize)
    plt.imshow(np.flipud((spread_dB).T),
               extent=[x.min() + h/2, x.max() + h/2,
                       y.min(), y.max()],
               clim=(-60, 0), cmap=cc.cm.gouldian)
    plt.colorbar()
    plt.xlim(x.min() - margin, x.max() + margin)
    plt.ylim(y.min() - margin, y.max() + margin)
    plt.title(f'Spreading factor [dB] ({plotname})')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    draw_boundary()
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{name}_spread_dB.pdf')

    vmax = np.nanmax(abs(u[np.isfinite(u)]))
    vmin = -vmax
    plt.figure(figsize=figsize)
    plt.imshow(np.flipud(np.real(u).T),
               extent=[x.min() + h/2, x.max() + h/2, y.min(), y.max()],
               clim=(vmin, vmax),
               cmap=cc.cm.gouldian)
    plt.colorbar()
    plt.xlim(x.min() - margin, x.max() + margin)
    plt.ylim(y.min() - margin, y.max() + margin)
    plt.title(f'$u$ ({plotname})')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    draw_boundary()
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{name}_u.pdf')

    plt.figure(figsize=figsize)
    plt.contour(x, y, origin, [0.5], colors=['k'], linewidths=[1.5],
                linestyles=['--'])
    plt.imshow(np.flipud(origin.T), extent=[x.min(), x.max(), y.min(), y.max()],
               cmap=cc.cm.gouldian)
    plt.colorbar()
    plt.xlim(x.min() - margin, x.max() + margin)
    plt.ylim(y.min() - margin, y.max() + margin)
    plt.title(r'$\mathtt{origin}$' + f'({plotname})')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    draw_boundary()
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{name}_origin.pdf')
