#!/usr/bin/env python

import colorcet as cc
import itertools as it
import jmm
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.optimize

from common import get_camera_basis, get_view_direction

LEVEL = 0.5

if __name__ == '__main__':
    path = '.'
    scene = 'box'
    indsrc = 16

    verts_bin_path = os.path.join(path, scene + '_verts.bin')
    verts = np.fromfile(verts_bin_path, np.float64)
    verts = verts.reshape(verts.size//3, 3)
    print('- read verts from %s' % verts_bin_path)

    cells_bin_path = os.path.join(path, scene + '_cells.bin')
    cells = np.fromfile(cells_bin_path, np.uintp)
    cells = cells.reshape(cells.size//4, 4)
    print('- reading cells from %s' % cells_bin_path)

    mesh = jmm.Mesh3.from_verts_and_cells(verts, cells)
    eik = jmm.Eik3(mesh)
    eik.add_trial(indsrc, jmm.Jet3(0, np.nan, np.nan, np.nan))
    eik.solve()

    bmesh = jmm.Bmesh33.from_eik3(eik)

    # Initialize an R-tree (don't build the BVH yet!) with the
    # triangles in the surface mesh of our domain. For each level in
    # our video, we'll copy this R-tree, insert the tetrahedra for the
    # current level, and build the R-tree before raytracing and
    # visualizing.
    surface_mesh = mesh.get_surface_mesh()
    base_rtree = jmm.Rtree()
    base_rtree.insert(surface_mesh)

    level_bmesh = bmesh.get_level_bmesh(LEVEL)

    level_rtree = base_rtree.copy()
    level_rtree.insert(level_bmesh.mesh)
    level_rtree.build()
