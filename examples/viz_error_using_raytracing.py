#!/usr/bin/env python

import colorcet as cc
import itertools as it
import jmm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.optimize

from common import get_camera_basis, get_view_direction
from PIL import Image

matplotlib.use('agg')
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
    print('- solved point source problem')

    bmesh = jmm.Bmesh33.from_eik3(eik)
    print('- assembled Bernstein-Bezier mesh for solution')

    # Initialize an R-tree (don't build the BVH yet!) with the
    # triangles in the surface mesh of our domain. For each level in
    # our video, we'll copy this R-tree, insert the tetrahedra for the
    # current level, and build the R-tree before raytracing and
    # visualizing.
    surface_mesh = mesh.get_surface_mesh()
    base_rtree = jmm.Rtree()
    base_rtree.insert(surface_mesh)
    print('- got surface mesh and inserted into BVH')

    level_bmesh = bmesh.get_level_bmesh(LEVEL)
    print('- got Bernstein-Bezier mesh for level (T = %f)' % LEVEL)

    level_rtree = base_rtree.copy()
    level_rtree.insert(level_bmesh.mesh)
    level_rtree.build()
    print('- built BVH for level set')

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
    print('- set up camera rays')

    T = np.empty((num_rays,), dtype=np.float64)
    I = np.empty((num_rays,), dtype=np.float64)
    Type = np.empty((num_rays,), dtype=np.float64)
    for l, (org, dir_) in enumerate(zip(orgs, dirs)):
        assert abs(1 - np.linalg.norm(dir_)) < 1e-15
        isect = level_rtree.intersect(org, dir_)
        T[l] = isect.t if isect.hit else np.nan
        I[l] = isect.obj.astype(jmm.Mesh2Tri).index if isect.hit else np.nan
        Type[l] = isect.obj.type.value if isect.hit else np.nan
    T = T.reshape(h, w)
    I = I.reshape(h, w)
    Type = Type.reshape(h, w)
    print('- finished raytracing')

    plt.figure(figsize=(12, 4))

    plt.subplot(1, 3, 1)
    plt.imshow(I, cmap=cc.cm.rainbow, interpolation='none')
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title("index")

    plt.subplot(1, 3, 2)
    plt.imshow(T, cmap=cc.cm.bmw)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title("$t$")

    plt.subplot(1, 3, 3)
    plt.imshow(Type, cmap=cc.cm.rainbow)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title('Type')

    plt.tight_layout()

    plt.show()

    print('- finished making plot')

    tet_inds = set()
    TriColor = np.array([0.5, 0.5, 1.0])
    TriAlpha = 0.1
    Img = np.empty((h, w, 3), dtype=np.float64)
    Img[:] = 1
    for (i, j), org, dir_ in zip(it.product(range(h), range(w)), orgs.copy(), dirs.copy()):
        isect = level_rtree.intersect(org, dir_)
        while isect.hit:
            if isect.obj.type.value == 0: # Mesh2Tri
                Img[i, j] *= 1 - TriAlpha
                Img[i, j] += TriAlpha*TriColor
            elif isect.obj.type.value == 1:
                tetra = isect.obj.astype(jmm.Mesh3Tetra)
                tet_inds.add(tetra.index) # For debugging...
                cell = level_bmesh.get_cell(tetra.index)
                ray = jmm.Ray3(org, dir_)
                print(f'cell.ray_intersects_level(l = {tetra.index})')
                b = cell.ray_intersects_level(ray, LEVEL)
                if b is not None:
                    print(b)
            org += (isect.t + 1e-10)*dir_
            isect = level_rtree.intersect(org, dir_)
    print('- finished raytracing for image')
    plt.figure()
    plt.imshow(Img, cmap=cc.cm.gray, vmin=0, vmax=1, interpolation='none')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()
    print('- plotted image')

    f = np.array([0, 1, 1, 1], dtype=np.float64)
    Df = np.zeros((4, 3), dtype=np.float64); Df[1, 0] = 1; Df[2, 1] = 1; Df[3, 2] = 1;
    X = Df.copy()
    bb = jmm.Bb33.from_3d_data(f, Df, X)
    print('- set up test Bezier tetrahedron')
