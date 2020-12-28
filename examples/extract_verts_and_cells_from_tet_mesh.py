#!/usr/bin/env python

import argparse
import meshio
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=str)
    args = parser.parse_args()

    root = args.root

    mesh_path = '%s.1.vtk' % root
    mesh = meshio.read(mesh_path)
    print('- read %s' % mesh_path)

    verts_bin_path = '%s_verts.bin' % root
    verts = mesh.points.astype(np.float64)
    with open(verts_bin_path, 'wb') as f:
        verts.tofile(f)
    print('- wrote %s' % verts_bin_path)

    cells_bin_path = '%s_cells.bin' % root
    cells = mesh.cells_dict['tetra'].astype(np.uintp)
    with open(cells_bin_path, 'wb') as f:
        cells.tofile(f)
    print('- wrote %s' % cells_bin_path)
