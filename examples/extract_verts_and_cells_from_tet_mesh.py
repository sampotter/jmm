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
    indsrc = np.argmin(np.sqrt(np.sum(verts**2, axis=1)))
    if np.linalg.norm(verts[indsrc]) > np.finfo(np.float64).resolution:
        print('- rounding point closest to origin to (0, 0, 0)')
        verts[indsrc, :] = 0
        mesh.points[indsrc, :] = 0
        mesh.write(mesh_path)
        print('- wrote corrected %s' % mesh_path)
    print('- point source index: %d' % indsrc)
    with open(verts_bin_path, 'wb') as f:
        verts.tofile(f)
    print('- wrote %s' % verts_bin_path)

    cells_bin_path = '%s_cells.bin' % root
    cells = mesh.cells_dict['tetra'].astype(np.uintp)
    with open(cells_bin_path, 'wb') as f:
        cells.tofile(f)
    print('- wrote %s' % cells_bin_path)
