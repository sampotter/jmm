#!/usr/bin/env python

import argparse
import itertools as it
import meshio
import numpy as np
import os
import pyvista as pv

from make_box_mesh import faces as box_faces

MARGIN_PERCENT = 0.1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str)
    args = parser.parse_args()

    try:
        grid = pv.read(args.path)
    except:
        print(f'! failed to read {args.path} (unknown file format?)')
        exit(1)
    print(f'- loaded mesh from {args.path}')
    print(f'- bounds: x = (%g, %g), y = (%g, %g), z = (%g, %g)' %
          tuple(grid.bounds))

    xm, ym, zm = (grid.points.min(0) + grid.points.max(0))/2
    dx, dy, dz = grid.points.ptp(0)/2
    p = MARGIN_PERCENT
    bounding_box = np.array([
        (
            (1 + p)*dx*(-1)**i + xm,
            (1 + p)*dy*(-1)**j + ym,
            (1 + p)*dz*(-1)**k + zm
        )
        for i, j, k in it.product(range(1, 3), repeat=3)
    ], dtype=np.float64)

    bounds = np.empty(6)
    bounds[::2] = bounding_box.min(0)
    bounds[1::2] = bounding_box.max(0)
    print(f'- new bounds: x = (%g, %g), y = (%g, %g), z = (%g, %g)' %
          tuple(bounds))

    # NOTE: this is assumes that the verts in `box_verts` are ordered
    # in the same way as `bounding_box`. This is true right now since
    # `box_verts` was generated using `it.product(range(2),
    # repeat=3)`. But this is a dependency that could break...

    verts = grid.points.astype(np.float64)

    faces = grid.faces.reshape(grid.faces.size//4, 4)
    faces = faces[:, 1:].astype(np.uintp)

    faces = np.concatenate([faces, box_faces.astype(np.uintp) + verts.shape[0]],
                           dtype=faces.dtype, axis=0)
    verts = np.concatenate([verts, bounding_box],
                           dtype=verts.dtype, axis=0)

    filename = os.path.splitext(args.path)[0]
    outpath = os.path.join(filename + '_bd.off')

    cells = [('triangle', faces)]
    mesh = meshio.Mesh(verts, cells)
    mesh.write(outpath)
    print(f'- wrote mesh with bounding box to {outpath}')
