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

def plot_sphere(p, x, r, c):
    p.add_mesh(pv.Sphere(r, x), color=c)

if __name__ == '__main__':
    path = Path('../../build/examples/simple_building')

    verts = np.fromfile(path/'verts.bin', dtype=np.float64).reshape(-1, 3)
    cells = np.fromfile(path/'cells.bin', dtype=np.uintp).reshape(-1, 4)
    num_verts, num_cells = verts.shape[0], cells.shape[0]


    org = np.fromfile(path/'org.bin')
    grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)
    points = pv.PolyData(verts)
    points['org'] = org
    # points['T'] = T
    p = pvqt.BackgroundPlotter()
    p.add_mesh(grid, color='white', show_edges=True, opacity=0.25)
    p.add_mesh(points, scalars='org', cmap=cc.cm.fire, nan_opacity=0, clim=(0, 1))
    # p.add_mesh(points, scalars='T', cmap=cc.cm.rainbow)
    p.add_mesh(pv.Sphere(0.1, (5, -5, 1)), 'red')
    p.enable_point_picking()
    plot_sphere(p, verts[1416], 0.1, 'cyan')
    plot_sphere(p, verts[1260], 0.1, 'orange')

    jet = np.fromfile(path/f'jet.bin', dtype=np.float64).reshape(-1, 4)
    jet[~np.isfinite(jet)] = np.nan

    T = jet[:, 0]
    grad_T = jet[:, 1:4]

    # make 3D plot

    grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, verts)

    points = pv.PolyData(verts)
    points['T'] = T
    points['grad_T'] = grad_T

    arrows = points.glyph(factor=0.25, orient='grad_T')

    T_cmap = cc.cm.rainbow

    p = pvqt.BackgroundPlotter()
    p.add_mesh(grid, color='white', show_edges=False, opacity=0.1)
    p.add_mesh(points, scalars='T', cmap=T_cmap, nan_opacity=0)
    p.add_mesh(arrows, scalars='T', cmap=T_cmap, nan_opacity=0, opacity=1)

    lhat = 3150
    # l0 = 5732
    # l1 = 181
    # l2 = 7080
    # # l3 = 638

    plot_sphere(p, verts[lhat], 0.1, 'red')
    # plot_sphere(p, verts[l0], 0.1, 'blue')
    # plot_sphere(p, verts[l1], 0.1, 'cyan')
    # plot_sphere(p, verts[l2], 0.1, 'teal')
    # # plot_sphere(p, verts[l3], 0.1, 'pink')
    # for l in np.unique(cells[(cells == lhat).any(1)]):
    #     plot_sphere(p, verts[l], 0.025, 'black')

    # lf = np.unique([91, 128, 1703, 91, 1702, 1703, 126, 130, 14320, 126, 181, 771, 126, 181, 14320, 126, 771, 772, 130, 182, 1908, 130, 182, 14320, 130, 673, 1908, 131, 183, 11104, 131, 11104, 14232, 225, 227, 7505, 225, 227, 14813, 226, 228, 1096, 226, 228, 11097, 226, 771, 11097, 227, 774, 7505, 227, 1495, 14813, 228, 774, 1096, 239, 240, 2697, 239, 240, 14825, 239, 1526, 14825, 240, 673, 2697, 402, 403, 1494, 402, 403, 1702, 402, 1494, 1495, 402, 1702, 1703, 669, 670, 1526, 669, 670, 11104, 669, 1526, 14825, 669, 11104, 14232, 672, 673, 1908, 672, 673, 2697, 771, 772, 11097, 774, 775, 1096, 774, 775, 7505, 1494, 1495, 14813])
    # for l in lf:
    #     plot_sphere(p, verts[l], 0.035, 'orange')

    #     # 4653, 181, 7080

    # plot_sphere(p, verts[5732], 0.11, 'chartreuse')

    plot_sphere(p, verts[448], 0.1, 'yellow')
    plot_sphere(p, verts[125], 0.1, 'yellow')

    p.add_mesh(pv.Arrow(verts[132], [-0.16025809732166316, -0.0058130399910618297, 0.98705802808593912]), color='magenta')

#     plot_sphere(p, verts[3010], 0.1, 'orange')
    # p.add_mesh(pv.Arrow(verts[132], [-0.16025809732166316, -0.0058130399910618297, 0.98705802808593912]), color='magenta')
