#!/usr/bin/env python

import argparse
import json
import numpy as np
import pyvista as pv
import trimesh as tm
import vtk

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bd_path', type=str)
    parser.add_argument('surf_path', type=str)
    parser.add_argument('out_path', type=str)
    parser.add_argument('json_path', type=str)
    parser.add_argument('--plot', default=False, action='store_true')
    args = parser.parse_args()

    try:
        bd_grid = pv.read(args.bd_path)
    except:
        print(f'! failed to read tetrahedron mesh from {args.bd_path}')
        exit(1)

    try:
        surf_grid = pv.read(args.surf_path)
    except:
        print(f'! failed to read surfmesh from {args.surf_path}')
        exit(1)

    surf_mesh = tm.Trimesh(surf_grid.points, surf_grid.faces.reshape(-1, 4)[:, 1:])

    points = bd_grid.points.astype(np.float64).copy()
    num_points = points.shape[0]

    cells = bd_grid.cells.reshape(-1, 5)[:, 1:].astype(np.uintp).copy()
    cell_means = points[cells].mean(1)
    num_cells = cells.shape[0]

    dom_cell_mask = surf_mesh.contains(cell_means)

    num_dom_cells = dom_cell_mask.sum()
    print(f'- found {num_dom_cells} domain cells (out of {num_cells})')

    dom_point_mask = np.zeros(num_points, dtype=bool)
    dom_point_mask[np.unique(cells[dom_cell_mask].ravel())] = True
    num_dom_points = dom_point_mask.sum()
    print(f'- found {num_dom_points} domain points (out of {num_points})')

    dom_cell_perm = np.argsort(dom_cell_mask)[::-1]

    dom_point_perm = np.argsort(dom_point_mask)[::-1]
    dom_point_inv_perm = np.argsort(dom_point_perm)

    cells = cells[dom_cell_perm]
    points = points[dom_point_perm]

    # Adjust cell indices after permuting points
    cells = dom_point_inv_perm[cells]

    # Quick sanity check
    assert((cells >= num_dom_points).any(1).nonzero()[0][0] == num_dom_cells)

    out_grid = pv.UnstructuredGrid({vtk.VTK_TETRA: cells}, points)
    out_grid.save(args.out_path)
    print(f'- wrote sorted mesh to {args.out_path}')

    stats = {
        'num_dom_cells': int(num_dom_cells),
        'num_dom_points': int(num_dom_points)
    }

    with open(args.json_path, 'w') as f:
        json.dump(stats, f)
    print(f'- wrote stats to {args.json_path}')

    if args.plot:
        dom_points = points[:num_dom_points]
        dom_cells = cells[:num_dom_cells]
        dom_grid = pv.UnstructuredGrid({vtk.VTK_TETRA: dom_cells}, dom_points)
        dom_grid.plot(show_edges=True)
