#!/usr/bin/env python

import argparse
import colorcet as cc
import jmm
import json
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
import sys
import time
import vtk

t0 = np.inf

def tic():
    global t0
    t0 = time.time()

def toc():
    global t0
    return time.time() - t0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('scene', type=str)
    parser.add_argument('lsrc', type=int)
    args = parser.parse_args()

    tic()

    try:
        surf_grid = pv.read(f'{args.scene}.obj')
    except:
        print(f'! failed to read surface mesh from {args.scene}.obj')
        sys.exit(1)

    try:
        grid = pv.read(f'{args.scene}.vtu')
    except:
        print(f'! failed to read tetrahedron mesh from {args.scene}.vtu')
        sys.exit(1)

    try:
        with open(f'{args.scene}.json', 'r') as f:
            info = json.load(f)
    except:
        print(f'! failed to read scene info from {args.scene}.json')
        sys.exit(1)

    cells = grid.cells.reshape(-1, 5)[:, 1:].astype(np.uintp)
    points = grid.points.astype(np.float64)

    num_dom_cells = info['num_dom_cells']
    num_dom_points = info['num_dom_points']

    dom_cells = cells[:num_dom_cells]
    dom_points = points[:num_dom_points]

    mesh_dom = jmm.Mesh3.from_verts_and_cells(dom_points, dom_cells)
    mesh = jmm.Mesh3.from_verts_and_cells(points, cells)

    print(f'- set up meshes [%1.2fs]' % toc())

    tic()

    eik = jmm.Eik3(mesh)
    eik.add_trial(args.lsrc, jmm.Jet3(0, np.nan, np.nan, np.nan))
    eik.solve()

    print('- solved on extended domain [%1.2fs]' % toc())

    tic()

    eik_dom = jmm.Eik3(mesh_dom)
    eik_dom.add_trial(args.lsrc, jmm.Jet3(0, np.nan, np.nan, np.nan))
    eik_dom.solve()

    print('- solved on domain [%1.2fs]' % toc())

    T = np.array([_[0] for _ in eik.jet[:num_dom_points]])
    T_diff = np.array([_[0] for _ in eik_dom.jet])

    T[~np.isfinite(T)] = -abs(T[np.isfinite(T)]).max()
    T_diff[~np.isfinite(T_diff)] = -abs(T_diff[np.isfinite(T_diff)]).max()

    D = T_diff - T
    Z = D > 0.007

    dom_grid = pv.UnstructuredGrid({vtk.VTK_TETRA: dom_cells}, dom_points)
    dom_grid.point_arrays['Z'] = Z
    dom_grid.point_arrays['T'] = T
    dom_grid.point_arrays['T_diff'] = T_diff

    if points[args.lsrc, 1] in {0, 4.5}:
        plotter = pvqt.BackgroundPlotter()
        plotter.add_mesh(dom_grid, scalars='Z', cmap=cc.cm.gray_r, show_scalar_bar=False)
        # plotter.add_mesh(dom_grid, scalars='T', cmap=cc.cm.coolwarm)
        # plotter.add_mesh(dom_grid, scalars='T_diff', cmap=cc.cm.coolwarm)
        plotter.add_mesh(pv.Sphere(0.1, points[args.lsrc]), color='pink')
    else:
        clipped_dom_grid = dom_grid.clip('y', origin=(0, points[args.lsrc, 1], 0))
        plotter = pvqt.BackgroundPlotter()
        plotter.add_mesh(clipped_dom_grid, scalars='Z', cmap=cc.cm.gray_r, show_scalar_bar=False)
        # plotter.add_mesh(clipped_dom_grid, scalars='T', cmap=cc.cm.coolwarm)
        # plotter.add_mesh(clipped_dom_grid, scalars='T_diff', cmap=cc.cm.coolwarm)
        plotter.add_mesh(pv.Sphere(0.1, points[args.lsrc]), color='pink')
