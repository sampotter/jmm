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

if __name__ == '__main__':
    path = Path('../../build/examples/simple_building')

    verts = np.fromfile(path/'verts.bin', dtype=np.float64).reshape(-1, 3)
    cells = np.fromfile(path/'cells.bin', dtype=np.uintp).reshape(-1, 4)
    num_verts, num_cells = verts.shape[0], cells.shape[0]

    jet = np.fromfile(path/f'jet.bin', dtype=np.float64).reshape(-1, 4)
    jet[~np.isfinite(jet)] = np.nan
    T = jet[:, 0]
    grad_T = jet[:, 1:4]

    # make 3D plot

    grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)
    points = pv.PolyData(verts)
    points['T'] = T
    plotter = pvqt.BackgroundPlotter()
    plotter.background_color = 'white'
    plotter.add_mesh(grid, color='white', show_edges=False, opacity=0.1)
    plotter.add_mesh(points, scalars='T', cmap=cc.cm.gouldian, nan_opacity=0)

    # plot 2D slice

    img_grid = read_img_grid(path/'img_grid.txt')
    x, y = get_grid(*img_grid)
    shape = img_grid[0]
    h = img_grid[2]

    freq = 3000
    c = 340

    T = np.fromfile(path/'T_slice.bin').reshape(shape)
    spread = np.fromfile(path/'spread_slice.bin').reshape(shape)
    spread_dB = 20*np.log10(np.maximum(1e-16, spread))
    origin = np.fromfile(path/'origin_slice.bin').reshape(shape)
    u = origin*spread*np.exp(1j*freq*T/c)

    margin = 0.5

    figsize = (5, 4)

    plt.figure(figsize=figsize)
    plt.contourf(x, y, T/c, levels=21, cmap=cc.cm.gouldian)
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.xlim(x.min() - margin, x.max() + margin)
    plt.ylim(y.min() - margin, y.max() + margin)
    plt.title('$T$ [s] (Direct)')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    draw_boundary()
    plt.show()
    plt.tight_layout()
    plt.savefig('T.pdf')

    plt.figure(figsize=figsize)
    plt.imshow(np.flipud((spread_dB).T),
               extent=[x.min() + h/2, x.max() + h/2,
                       y.min(), y.max()],
               clim=(-60, 0), cmap=cc.cm.gouldian)
    plt.colorbar()
    plt.xlim(x.min() - margin, x.max() + margin)
    plt.ylim(y.min() - margin, y.max() + margin)
    plt.title('Spreading factor [dB] (Direct)')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    draw_boundary()
    plt.tight_layout()
    plt.show()
    plt.savefig('spread_dB.pdf')

    plt.figure(figsize=figsize)
    plt.imshow(np.flipud(np.real(u).T),
               extent=[x.min() + h/2, x.max() + h/2, y.min(), y.max()],
               cmap=cc.cm.gouldian)
    plt.colorbar()
    plt.xlim(x.min() - margin, x.max() + margin)
    plt.ylim(y.min() - margin, y.max() + margin)
    plt.title('$u$ (Direct)')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    draw_boundary()
    plt.tight_layout()
    plt.show()
    plt.savefig('u.pdf')

    plt.figure(figsize=figsize)
    plt.contour(x, y, origin, [0.5], colors=['k'], linewidths=[1.5],
                linestyles=['--'])
    plt.imshow(np.flipud(origin.T), extent=[x.min(), x.max(), y.min(), y.max()],
               cmap=cc.cm.gouldian)
    plt.colorbar()
    plt.xlim(x.min() - margin, x.max() + margin)
    plt.ylim(y.min() - margin, y.max() + margin)
    plt.title(r'$\mathtt{origin}$ (Direct)')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    draw_boundary()
    plt.tight_layout()
    plt.show()
    plt.savefig('origin.pdf')
