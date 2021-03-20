#!/usr/bin/env python

import colorcet as cc
import itertools as it
import jmm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import scipy.optimize

from common import get_camera_basis, get_view_direction
from PIL import Image

plt.ion()

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

    errors = []
    for vert, jet in zip(verts, eik.jet):
        tau = np.linalg.norm(vert)
        errors.append(jet[0] - tau)
    errors = np.array(errors)

    # plt.figure()
    # plt.scatter(np.sqrt(np.sum(verts**2, axis=1)), errors, s=3, c='k')
    # plt.show()

    e_max = abs(errors).max()

    print('- computed errors (max abs error: %g)' % e_max)

    bmesh = jmm.Bmesh33.from_eik3(eik)

    print('- assembled Bernstein-Bezier mesh for solution')

    # Set up camera
    pos = (-2, 0, 0)
    look = (1, 0, 0)
    up = (0, 0, 1)
    left = np.cross(up, look)
    shape = (512, 512)
    width = 2.02
    height = 2.02
    camera = jmm.Camera.make_orthographic(pos, look, left, up, shape, width, height)
    print('- set up camera')

    # Define extent for imshow
    extent = [pos[1] - width/2, pos[1] + width/2,
              pos[2] - height/2, pos[2] + height/2]

    # Set the clim based on the two-sigma value of the absolute errors
    vmax = 2*abs(errors).std()

    levels = np.linspace(0.0, 1.8, 18*30 + 1)

    for i, LEVEL in enumerate(levels):
        print('level = %g (%d/%d)' % (LEVEL, i + 1, levels.size))

        level_bmesh = bmesh.restrict_to_level(LEVEL)
        print('- got Bernstein-Bezier mesh for level (T = %f)' % LEVEL)

        level_rtree = jmm.Rtree()
        level_rtree.insert(level_bmesh)
        level_rtree.build()
        print('- built BVH for level set')

        T = np.empty((camera.num_rays,), dtype=np.float64)
        I = np.empty((camera.num_rays,), dtype=np.float64)
        # Type = np.empty((camera.num_rays,), dtype=np.float64)
        Error = np.empty((camera.num_rays,), dtype=np.float64)

        # Set default values
        T[...] = np.nan
        I[...] = np.nan
        # Type[...] = np.nan
        Error[...] = np.nan

        for l, ray in enumerate(camera.get_image_rays()):
            isect = level_rtree.intersect(ray)

            if not isect.hit:
                continue

            T[l] = isect.t
            I[l] = isect.obj.astype(jmm.Bmesh33Cell).index
            # Type[l] = isect.obj.type.value
            Error[l] = LEVEL - np.linalg.norm(ray.get_point(T[l]))


        T = T.reshape(shape)
        I = I.reshape(shape)
        # Type = Type.reshape(shape)
        Error = Error.reshape(shape)

        print('- finished raytracing')

        fig = plt.figure(figsize=(11, 6))

        gs = gridspec.GridSpec(2, 3, figure=fig)

        ax = fig.add_subplot(gs[0, 0])
        im = ax.imshow(I, extent=extent, cmap=cc.cm.rainbow,
                       interpolation='none')
        # fig.colorbar(im, ax=ax)
        ax.set_aspect('equal')
        ax.set_title("cell index")

        ax = fig.add_subplot(gs[1, 0])
        im = ax.imshow(T, vmin=1, vmax=3, extent=extent,
                       cmap=cc.cm.rainbow)
        fig.colorbar(im, ax=ax)
        ax.set_aspect('equal')
        ax.set_title("$t$")

        ax = fig.add_subplot(gs[:, 1:])
        im = ax.imshow(Error, vmin=-vmax, vmax=vmax, extent=extent,
                       cmap=cc.cm.coolwarm)
        fig.colorbar(im, ax=ax)
        ax.set_aspect('equal')
        ax.set_title(r'$e(x) = T(x) - \tau(x)$')

        fig.suptitle('level = %1.4f' % LEVEL)

        fig.tight_layout()

        fig.savefig('frames/frame%03d.png' % i)

        plt.close(fig)

        print('- finished making "info" plot')
