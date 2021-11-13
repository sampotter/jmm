#!/usr/bin/env python

import argparse
import meshio
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str)

    args = parser.parse_args()

    mesh = meshio.read(args.path)
    print('- read %s' % args.path)

    scene_name = args.path.split('.')[0]
    print('- using scene name "%s"' % scene_name)

    verts_bin_path = '%s_verts.bin' % scene_name
    verts = mesh.points.astype(np.float64)
    with open(verts_bin_path, 'wb') as f:
        verts.tofile(f)
    print('- wrote %s' % verts_bin_path)

    cells_bin_path = '%s_cells.bin' % scene_name
    cells = mesh.cells_dict['tetra'].astype(np.uintp)
    with open(cells_bin_path, 'wb') as f:
        cells.tofile(f)
    print('- wrote %s' % cells_bin_path)
