#!/usr/bin/env python

import colorcet as cc
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import sys

from glob import glob

def max_area_from_path(path):
    return float('.'.join(path.split('.')[:-1]).split('_')[-1])

if __name__ == '__main__':
    verts_paths = glob('out/verts_*.bin')
    verts_paths = sorted(verts_paths, key=max_area_from_path)

    faces_paths = glob('out/faces_*.bin')
    faces_paths = sorted(faces_paths, key=max_area_from_path)

    jet_gt_paths = glob('out/eik2m1_jets_gt_*.bin')
    jet_gt_paths = sorted(jet_gt_paths, key=max_area_from_path)

    jet_paths = glob('out/eik2m1_jets_*.bin')
    jet_paths = [_ for _ in jet_paths if 'gt' not in _]
    jet_paths = sorted(jet_paths, key=max_area_from_path)

    max_areas_gt = sorted(max_area_from_path(path) for path in jet_gt_paths)
    max_areas = sorted(max_area_from_path(path) for path in jet_paths)

    if max_areas != max_areas_gt:
        raise RuntimeError('found mismatched number of jet and gt jet files')

    verts = {
        max_area: np.fromfile(path, dtype=np.float64).reshape(-1, 2)
        for max_area, path in zip(max_areas, verts_paths)
    }

    faces = {
        max_area: np.fromfile(path, dtype=np.uintp).reshape(-1, 3)
        for max_area, path in zip(max_areas, faces_paths)
    }

    jet_gt = {
        max_area: np.fromfile(path, dtype=np.float64).reshape(-1, 7)
        for max_area, path in zip(max_areas, jet_gt_paths)
    }

    jet = {
        max_area: np.fromfile(path, dtype=np.float64).reshape(-1, 7)
        for max_area, path in zip(max_areas, jet_paths)
    }

    ########################################################################
    # make pointwise error plots

    print('making pointwise error plots')

    for max_area in max_areas:
        print(f'- max area: {max_area}')

        V, F = verts[max_area], faces[max_area]

        E = jet[max_area] - jet_gt[max_area]

        vmax = np.nanmax(abs(E), axis=0)
        vmin = -vmax

        def get_kwargs(i):
            return {
                'vmin': vmin[i],
                'vmax': vmax[i],
                'cmap': cc.cm.coolwarm
            }

        def make_subplot(loc, V, F, X, i, title):
            tri = mtri.Triangulation(*V.T, F)

            Y = X.copy()
            Y[~np.isfinite(Y)] = 0

            plt.subplot(*loc)
            plt.tricontourf(tri, Y, **get_kwargs(i))
            plt.colorbar()
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.gca().set_aspect('equal')
            plt.title(title)

        plt.figure(figsize=(14, 4))
        make_subplot((2, 4, 1), V, F, E[:, 0], 0, r'$T - \tau$')
        make_subplot((2, 4, 2), V, F, E[:, 1], 1, r'$T_x - \tau_x$')
        make_subplot((2, 4, 6), V, F, E[:, 2], 2, r'$T_y - \tau_y$')
        make_subplot((2, 4, 3), V, F, E[:, 3], 3, r'$T_{xx} - \tau_{xx}$')
        make_subplot((2, 4, 4), V, F, E[:, 4], 4, r'$T_{yx} - \tau_{yx}$')
        make_subplot((2, 4, 7), V, F, E[:, 5], 5, r'$T_{xy} - \tau_{xy}$')
        make_subplot((2, 4, 8), V, F, E[:, 6], 6, r'$T_{yy} - \tau_{yy}$')
        # plt.tight_layout()
        plt.savefig('out/eik2m1_E_%s.png' % max_area)
        plt.close()

    ########################################################################
    # make pointwise jet gt plots

    print('making pointwise groundtruth jet plots')

    for max_area in max_areas:
        print(f'- max_area: {max_area}')

        V, F = verts[max_area], faces[max_area]

        J = jet_gt[max_area]

        vmax = np.nanmax(abs(J), axis=0)
        vmin = -vmax
        vmin[0] = 0

        def get_kwargs(i):
            return {
                'vmin': vmin[i],
                'vmax': vmax[i],
                'cmap': cc.cm.blues if i == 0 else cc.cm.coolwarm
            }

        def make_subplot(loc, V, F, X, i, title):
            tri = mtri.Triangulation(*V.T, F)

            Y = X.copy()
            Y[~np.isfinite(Y)] = 0

            plt.subplot(*loc)
            plt.tricontourf(tri, Y, **get_kwargs(i))
            plt.colorbar()
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.gca().set_aspect('equal')
            plt.title(title)

        plt.figure(figsize=(14, 4))
        make_subplot((2, 4, 1), V, F, J[:, 0], 0, r'$\tau$')
        make_subplot((2, 4, 2), V, F, J[:, 1], 1, r'$\tau_x$')
        make_subplot((2, 4, 6), V, F, J[:, 2], 2, r'$\tau_y$')
        make_subplot((2, 4, 3), V, F, J[:, 3], 3, r'$\tau_{xx}$')
        make_subplot((2, 4, 4), V, F, J[:, 4], 4, r'$\tau_{yx}$')
        make_subplot((2, 4, 7), V, F, J[:, 5], 5, r'$\tau_{xy}$')
        make_subplot((2, 4, 8), V, F, J[:, 6], 6, r'$\tau_{yy}$')
        # plt.tight_layout()
        plt.savefig('out/eik2g1_jet_gt_%s.png' % max_area)
        plt.close()

    ########################################################################
    # make pointwise jet plots

    print('making pointwise jet plots')

    for max_area in max_areas:
        print(f'- max area: {max_area}')

        V, F = verts[max_area], faces[max_area]

        J = jet[max_area]

        vmax = np.nanmax(abs(J), axis=0)
        vmin = -vmax
        vmin[0] = 0

        def get_kwargs(i):
            return {
                'extent': [-1, 1, -1, 1],
                'vmin': vmin[i],
                'vmax': vmax[i],
                'cmap': cc.cm.blues if i == 0 else cc.cm.coolwarm
            }

        def make_subplot(loc, V, F, X, i, title):
            tri = mtri.Triangulation(*V.T, F)

            Y = X.copy()
            Y[~np.isfinite(Y)] = 0

            plt.subplot(*loc)
            plt.tricontourf(tri, Y, **get_kwargs(i))
            plt.colorbar()
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.gca().set_aspect('equal')
            plt.title(title)

        plt.figure(figsize=(14, 4))
        make_subplot((2, 4, 1), V, F, J[:, 0], 0, r'$T$')
        make_subplot((2, 4, 2), V, F, J[:, 1], 1, r'$T_x$')
        make_subplot((2, 4, 6), V, F, J[:, 2], 2, r'$T_y$')
        make_subplot((2, 4, 3), V, F, J[:, 3], 3, r'$T_{xx}$')
        make_subplot((2, 4, 4), V, F, J[:, 4], 4, r'$T_{yx}$')
        make_subplot((2, 4, 7), V, F, J[:, 5], 5, r'$T_{xy}$')
        make_subplot((2, 4, 8), V, F, J[:, 6], 6, r'$T_{yy}$')
        plt.tight_layout()
        plt.savefig('out/eik2g1_jet_%s.png' % max_area)
        plt.close()

    ########################################################################
    # make convergence plots

    print('making convergence plots')

    error_rel2_J = []
    for max_area in max_areas:
        I = np.arange(1)
        e = (jet[max_area][:, I] - jet_gt[max_area][:, I]).ravel()
        mask = np.isfinite(e)
        error_rel2_J.append(np.linalg.norm(e[mask])/np.linalg.norm(jet[max_area][:, I].ravel()[mask]))

    error_rel2_DJ = []
    for max_area in max_areas:
        I = np.arange(1, 3)
        e = (jet[max_area][:, I] - jet_gt[max_area][:, I]).ravel()
        mask = np.isfinite(e)
        error_rel2_DJ.append(np.linalg.norm(e[mask])/np.linalg.norm(jet[max_area][:, I].ravel()[mask]))

    error_rel2_D2J = []
    for max_area in max_areas:
        I = np.arange(3, 7)
        e = (jet[max_area][:, I] - jet_gt[max_area][:, I]).ravel()
        mask = np.isfinite(e)
        error_rel2_D2J.append(np.linalg.norm(e[mask])/np.linalg.norm(jet[max_area][:, I].ravel()[mask]))

    A = np.array([float(_) for _ in max_areas])
    H = np.sqrt(np.array(A))

    plt.figure(figsize=(4, 4))
    plt.loglog(H, error_rel2_J, label=r'$\tau - T$', marker='*')
    plt.loglog(H, error_rel2_DJ, label=r'$\nabla\tau - \nabla T$', marker='*')
    plt.loglog(H, error_rel2_D2J, label=r'$\nabla^2\tau - \nabla^2 T$', marker='*')
    plt.legend()
    plt.xlabel(r'$h = \sqrt{A_{\max}}$')
    plt.ylabel(r'Relative $\ell_2$ error')
    plt.tight_layout()
    plt.savefig('out/eik2m1_error_rel2.png')
    plt.close()

    ########################################################################
    # fit Ch^p to each set of errors

    print('writing convergence rates to out/eik2m1_convergence_rates.txt')

    p_J, C_J = np.polyfit(np.log(H), np.log(error_rel2_J), 1)
    p_DJ, C_DJ = np.polyfit(np.log(H), np.log(error_rel2_DJ), 1)
    p_D2J, C_D2J = np.polyfit(np.log(H), np.log(error_rel2_D2J), 1)

    with open('out/eik2m1_convergence_rates.txt', 'w') as f:
        print('J: p = %1.3f (C = %1.3f)' % (p_J, C_J), file=f)
        print('DJ: p = %1.3f (C = %1.3f)' % (p_DJ, C_DJ), file=f)
        print('D2J: p = %1.3f (C = %1.3f)' % (p_D2J, C_D2J), file=f)
