#!/usr/bin/env python

import argparse
import colorcet as cc
import itertools as it
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import scipy.optimize

import jmm

from PIL import Image
from pathlib import Path

norm = np.linalg.norm

matplotlib.use('Agg')

def get_tau_diff(x, xsrc):
    x0 = np.array([1, 1, 0])
    dx = np.array([0, 0, 1])

    def f(t):
        xt = x0 + t*dx
        return norm(x - xt) + norm(xt - xsrc)

    def ft(t):
        xt = x0 + t*dx
        y, ysrc = xt - x, xt - xsrc
        return np.dot(dx, y/norm(y) + ysrc/norm(ysrc))

    t = scipy.optimize.brentq(ft, 0, 1)

    return f(t)

def get_tau(x, xsrc):
    assert xsrc[0] >= 1 or xsrc[1] >= 1
    if x[0] < 1 and x[1] < 1:
        return np.nan
    elif xsrc[0] < 1:
        if x[0] < 1:
            return norm(x - xsrc)
        m = (xsrc[1] - 1)/(xsrc[0] - 1)
        y0 = m*(x[0] - 1) + 1
        if x[1] > y0:
            return norm(x - xsrc)
        return get_tau_diff(x, xsrc)
    elif xsrc[1] < 1:
        if x[1] < 1:
            return norm(x - xsrc)
        m = (xsrc[0] - 1)/(xsrc[1] - 1)
        b = 1 - m
        x0 = (1 - b)*x[1] + b
        if x0 < x[0]:
            return norm(x - xsrc)
        return get_tau_diff(x, xsrc)
    else:
        return norm(x - xsrc)

def get_default_camera_array_dict(camera):
    arrs = dict()
    for key in camera.keys():
        arrs[key] = np.empty((camera[key].num_rays,), dtype=np.float64)
        arrs[key].fill(np.nan)
    return arrs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
Solve a point source in an L-shaped room for a given location extract the
numerical level sets and make plots from several different views showing
the raytracing information and error rendered on the level sets.
''')
    parser.add_argument('-l', '--indsrc', type=int, default=36,
                        help='index of vertex to use for source location',)
    parser.add_argument('-m', '--minval', type=float, default=0.0,
                        help='minimum eikonal value')
    parser.add_argument('-M', '--maxval', type=float, default=1.0,
                        help='minimum eikonal value ')
    parser.add_argument('-n', '--levels', type=int, default=11,
                        help='number of level sets to plot')
    parser.add_argument('-p', '--outpath', type=str, default='frames',
                        help='path to which to write output frames')

    args = parser.parse_args()

    path = '.'
    scene = 'L'
    indsrc = args.indsrc
    levels = np.linspace(args.minval, args.maxval, args.levels)

    verts_bin_path = os.path.join(path, scene + '_verts.bin')
    verts = np.fromfile(verts_bin_path, np.float64)
    verts = verts.reshape(verts.size//3, 3)

    xsrc = verts[indsrc]

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

    errors = np.empty(verts.shape[0])
    for l, (vert, jet) in enumerate(zip(verts, eik.jet)):
        if l == indsrc:
            errors[l] = 0
        else:
            tau = get_tau(vert, xsrc)
            errors[l] = np.nan if np.isinf(jet[0]) else jet[0] - tau
    e_max = np.nanmax(abs(errors))

    print('- computed errors (max abs error: %g)' % e_max)

    bmesh = jmm.Bmesh33.from_eik3(eik)

    print('- assembled Bernstein-Bezier mesh for solution')

    camera = {
        'front': jmm.Camera.make_orthographic(
            pos=(-1, 1, 0.5),
            look=(1, 0, 0),
            left=(0, 1, 0),
            up=(0, 0, 1),
            shape=(512, 256),
            width=2.02,
            height=1.02
        ),
        'side': jmm.Camera.make_orthographic(
            pos=(1, -1, 0.5),
            look=(0, 1, 0),
            left=(-1, 0, 0),
            up=(0, 0, 1),
            shape=(512, 256),
            width=2.02,
            height=1.02
        ),
        'bottom': jmm.Camera.make_orthographic(
            pos=(1, 1, -1),
            look=(0, 0, 1),
            left=(-1, 0, 0),
            up=(0, 1, 0),
            shape=(512, 512),
            width=2.02,
            height=2.02
        ),
        'top': jmm.Camera.make_orthographic(
            pos=(1, 1, 2),
            look=(0, 0, -1),
            left=(1, 0, 0),
            up=(0, -1, 0),
            shape=(512, 512),
            width=2.02,
            height=2.02
        )
    }

    print('- set up cameras')

    # Set the clim based on the two-sigma value of the absolute errors
    # vmax = 2*np.nanstd(abs(errors))
    vmax = 0.01

    print(f'- using error clim: ({-vmax}, {vmax})')

    outpath = Path(args.outpath)
    outpath.mkdir(parents=True, exist_ok=True)

    print(f'- writing frames to {outpath}')

    for i, LEVEL in enumerate(levels):
        print('level = %g (%d/%d)' % (LEVEL, i + 1, len(levels)))

        level_bmesh = bmesh.restrict_to_level(LEVEL)
        print('- got Bernstein-Bezier mesh for level (T = %f)' % LEVEL)

        level_rtree = jmm.Rtree()
        level_rtree.insert(level_bmesh)
        level_rtree.build()
        print('- built BVH for level set')

        T = get_default_camera_array_dict(camera)
        I = get_default_camera_array_dict(camera)
        Error = get_default_camera_array_dict(camera)
        for key in camera.keys():
            for l, ray in enumerate(camera[key].get_image_rays()):
                isect = level_rtree.intersect(ray)
                if not isect.hit:
                    continue
                T[key][l] = isect.t
                I[key][l] = isect.obj.astype(jmm.Bmesh33Cell).index
                Error[key][l] = LEVEL - get_tau(ray.get_point(T[key][l]), xsrc)
            shape = camera[key].shape
            T[key] = T[key].reshape(shape)
            I[key] = I[key].reshape(shape)
            Error[key] = Error[key].reshape(shape)
        print('- finished raytracing')

        fig = plt.figure(figsize=(16, 12))

        gs = gridspec.GridSpec(4, 4, figure=fig)

        # Plot of hit index and hit distance for "front" camera

        ax = fig.add_subplot(gs[0, 0])
        im = ax.imshow(I['front'], extent=camera['front'].extent,
                       cmap=cc.cm.rainbow, interpolation='none')
        ax.set_aspect('equal')
        ax.set_title("cell index (front)")

        ax = fig.add_subplot(gs[0, 1])
        im = ax.imshow(T['front'], extent=camera['front'].extent,
                       cmap=cc.cm.rainbow, interpolation='none')
        ax.set_aspect('equal')
        ax.set_title("hit distance (front)")

        # Plot of hit index and hit distance for "side" camera

        ax = fig.add_subplot(gs[1, 0])
        im = ax.imshow(I['side'], extent=camera['side'].extent,
                       cmap=cc.cm.rainbow, interpolation='none')
        ax.set_aspect('equal')
        ax.set_title("cell index (side)")

        ax = fig.add_subplot(gs[1, 1])
        im = ax.imshow(T['side'], extent=camera['side'].extent,
                       cmap=cc.cm.rainbow, interpolation='none')
        ax.set_aspect('equal')
        ax.set_title("hit distance (side)")

        # Views of error

        ax = fig.add_subplot(gs[0, 2:])
        im = ax.imshow(Error['front'], extent=camera['front'].extent,
                       vmin=-vmax, vmax=vmax, cmap=cc.cm.coolwarm)
        fig.colorbar(im, ax=ax)
        ax.set_aspect('equal')
        ax.set_title(r'$e(x) = T(x) - \tau(x)$ (front)')

        ax = fig.add_subplot(gs[1, 2:])
        im = ax.imshow(Error['side'], extent=camera['side'].extent,
                       vmin=-vmax, vmax=vmax, cmap=cc.cm.coolwarm)
        fig.colorbar(im, ax=ax)
        ax.set_aspect('equal')
        ax.set_title(r'$e(x) = T(x) - \tau(x)$ (side)')

        ax = fig.add_subplot(gs[2:, :2])
        im = ax.imshow(Error['bottom'], extent=camera['bottom'].extent,
                       vmin=-vmax, vmax=vmax, cmap=cc.cm.coolwarm)
        fig.colorbar(im, ax=ax)
        ax.set_aspect('equal')
        ax.set_title(r'$e(x) = T(x) - \tau(x)$ (bottom)')

        ax = fig.add_subplot(gs[2:, 2:])
        im = ax.imshow(Error['top'], extent=camera['top'].extent,
                       vmin=-vmax, vmax=vmax, cmap=cc.cm.coolwarm)
        fig.colorbar(im, ax=ax)
        ax.set_aspect('equal')
        ax.set_title(r'$e(x) = T(x) - \tau(x)$ (top)')

        fig.suptitle('level = %1.4f' % LEVEL)

        fig.tight_layout()

        fig.savefig(outpath/('frame%03d.png' % i))

        plt.close(fig)

        print('- finished making "info" plot')
