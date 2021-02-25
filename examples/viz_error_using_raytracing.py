#!/usr/bin/env python

import colorcet as cc
import itertools as it
import jmm
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.optimize

from common import get_camera_basis, get_view_direction

plt.ion()

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

    # Compute camera frame (front x left x up)
    origin = np.array([-4, 0, 0], dtype=np.float64)
    target = np.array([0, 0, 0], dtype=np.float64)
    up = np.array([0, 0, 1])
    left, front, up = get_camera_basis(origin, target, up)

    # Image size
    w, h = 512, 512
    num_rays = w*h

    # Camera parameters
    s = w/h # aspect ratio
    fov_y = 45 # vertical field of view

    Theta = np.deg2rad(np.linspace(-s*fov_y/2, s*fov_y/2, w))
    Phi = np.deg2rad(np.linspace(-fov_y/2, fov_y/2, h))

    orgs = np.outer(np.ones(num_rays), origin)
    dirs = np.array([get_view_direction(left, front, up, phi, theta)
                     for phi, theta in it.product(Phi, Theta)])

    T = np.empty((num_rays,), dtype=np.float64)
    I = np.empty((num_rays,), dtype=np.float64)
    for l, (org, dir_) in enumerate(zip(orgs, dirs)):
        isect = level_rtree.intersect(org, dir_)
        T[l] = isect.t if isect.hit else np.nan
        I[l] = isect.obj.astype(jmm.Mesh2Tri).index if isect.hit else np.nan
    T = T.reshape(h, w)
    I = I.reshape(h, w)

    plt.figure(figsize=(9, 4))

    plt.subplot(1, 2, 1)
    plt.imshow(I, cmap=cc.cm.rainbow, interpolation='none')
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title("index")

    plt.subplot(1, 2, 2)
    plt.imshow(T, cmap=cc.cm.bmw)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title("$t$")

    plt.tight_layout()

    plt.show()

    f = np.array([0, 1, 1, 1], dtype=np.float64)
    Df = np.zeros((4, 3), dtype=np.float64); Df[1, 0] = 1; Df[2, 1] = 1; Df[3, 2] = 1;
    X = Df.copy()
    bb = jmm.Bb33.from_3d_data(f, Df, X)
    print('- set up test Bezier tetrahedron')
