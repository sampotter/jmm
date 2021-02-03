#!/usr/bin/env python

import argparse
import meshio
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=str)
    parser.add_argument('--xsrc', type=str, default="(0, 0, 0)")
    args = parser.parse_args()

    root = args.root

    mesh_path = '%s.1.vtk' % root
    mesh = meshio.read(mesh_path)
    print('- read %s' % mesh_path)

    xsrc = eval(args.xsrc)
    print('- using %s for point source location' % (xsrc,))
    xsrc = np.array(xsrc, dtype=np.float64)

    verts_bin_path = '%s_verts.bin' % root
    verts = mesh.points.astype(np.float64)
    indsrc = np.argmin(np.sqrt(np.sum((xsrc - verts)**2, axis=1)))
    if np.linalg.norm(xsrc - verts[indsrc]) > np.finfo(np.float64).resolution:
        print('- rounding point closest to origin to point source')
        verts[indsrc, :] = xsrc
        mesh.points[indsrc, :] = xsrc
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
