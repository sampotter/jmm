#!/usr/bin/env python

import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import sys

from glob import glob

def n_from_path(path):
    return int(path.split('.')[0].split('_')[-1])

if __name__ == '__main__':
    jet_gt_paths = glob('out/jet_gt_*.bin')
    jet_gt_paths = sorted(jet_gt_paths, key=n_from_path)

    jet_paths = glob('out/jet_*.bin')
    jet_paths = [_ for _ in jet_paths if 'gt' not in _]
    jet_paths = sorted(jet_paths, key=n_from_path)

    N_gt = sorted(n_from_path(path) for path in jet_gt_paths)
    N = sorted(n_from_path(path) for path in jet_paths)

    if N != N_gt:
        raise RuntimeError('found mismatched number of jet and gt jet files')

    jet_gt = {
        n: np.fromfile(path, dtype=np.float64).reshape(-1, 7)
        for n, path in zip(N, jet_gt_paths)
    }
    for n, J in jet_gt.items():
        if J.shape[0] != n**2:
            raise RuntimeError('found jets with weird shape')

    jet = {
        n: np.fromfile(path, dtype=np.float64).reshape(-1, 7)
        for n, path in zip(N, jet_paths)
    }
    for n, J in jet.items():
        if J.shape[0] != n**2:
            raise RuntimeError('found jets with weird shape')

    ########################################################################
    # make pointwise error plots

    print('making pointwise error plots')

    for n in N:
        print(f'- n = {n}')

        E = jet[n] - jet_gt[n]

        vmax = np.nanmax(abs(E), axis=0)
        vmin = -vmax

        def get_kwargs(i):
            return {
                'extent': [-1, 1, -1, 1],
                'vmin': vmin[i],
                'vmax': vmax[i],
                'cmap': cc.cm.coolwarm
            }

        def make_subplot(loc, X, i, title):
            plt.subplot(*loc)
            plt.imshow(np.rot90(X), **get_kwargs(i))
            plt.colorbar()
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.title(title)

        plt.figure(figsize=(14, 4))
        make_subplot((2, 4, 1), E[:, 0].reshape(n, n), 0, r'$T - \tau$')
        make_subplot((2, 4, 2), E[:, 1].reshape(n, n), 1, r'$T_x - \tau_x$')
        make_subplot((2, 4, 6), E[:, 2].reshape(n, n), 2, r'$T_y - \tau_y$')
        make_subplot((2, 4, 3), E[:, 3].reshape(n, n), 3, r'$T_{xx} - \tau_{xx}$')
        make_subplot((2, 4, 4), E[:, 4].reshape(n, n), 4, r'$T_{yx} - \tau_{yx}$')
        make_subplot((2, 4, 7), E[:, 5].reshape(n, n), 5, r'$T_{xy} - \tau_{xy}$')
        make_subplot((2, 4, 8), E[:, 6].reshape(n, n), 6, r'$T_{yy} - \tau_{yy}$')
        plt.tight_layout()
        # plt.show()
        plt.savefig('out/eik2g1_E_%d.png' % n)
        plt.close()

    ########################################################################
    # make pointwise jet gt plots

    print('making pointwise groundtruth jet plots')

    for n in N:
        print(f'- n = {n}')

        J = jet_gt[n]

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

        def make_subplot(loc, X, i, title):
            plt.subplot(*loc)
            plt.imshow(np.rot90(X), **get_kwargs(i))
            plt.colorbar()
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.title(title)

        plt.figure(figsize=(14, 4))
        make_subplot((2, 4, 1), J[:, 0].reshape(n, n), 0, r'$\tau$')
        make_subplot((2, 4, 2), J[:, 1].reshape(n, n), 1, r'$\tau_x$')
        make_subplot((2, 4, 6), J[:, 2].reshape(n, n), 2, r'$\tau_y$')
        make_subplot((2, 4, 3), J[:, 3].reshape(n, n), 3, r'$\tau_{xx}$')
        make_subplot((2, 4, 4), J[:, 4].reshape(n, n), 4, r'$\tau_{yx}$')
        make_subplot((2, 4, 7), J[:, 5].reshape(n, n), 5, r'$\tau_{xy}$')
        make_subplot((2, 4, 8), J[:, 6].reshape(n, n), 6, r'$\tau_{yy}$')
        plt.tight_layout()
        # plt.show()
        plt.savefig('out/eik2g1_jet_gt_%d.png' % n)
        plt.close()

    print('making pointwise jet plots')

    for n in N:
        print(f'- n = {n}')

        J = jet[n]

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

        def make_subplot(loc, X, i, title):
            plt.subplot(*loc)
            plt.imshow(np.rot90(X), **get_kwargs(i))
            plt.colorbar()
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.title(title)

        plt.figure(figsize=(14, 4))
        make_subplot((2, 4, 1), J[:, 0].reshape(n, n), 0, r'$T$')
        make_subplot((2, 4, 2), J[:, 1].reshape(n, n), 1, r'$T_x$')
        make_subplot((2, 4, 6), J[:, 2].reshape(n, n), 2, r'$T_y$')
        make_subplot((2, 4, 3), J[:, 3].reshape(n, n), 3, r'$T_{xx}$')
        make_subplot((2, 4, 4), J[:, 4].reshape(n, n), 4, r'$T_{yx}$')
        make_subplot((2, 4, 7), J[:, 5].reshape(n, n), 5, r'$T_{xy}$')
        make_subplot((2, 4, 8), J[:, 6].reshape(n, n), 6, r'$T_{yy}$')
        plt.tight_layout()
        # plt.show()
        plt.savefig('out/eik2g1_jet_%d.png' % n)
        plt.close()

    ########################################################################
    # make convergence plots

    print('making convergence plots')

    error_rel2_J = []
    for n in N:
        I = np.arange(1)
        e = (jet[n][:, I] - jet_gt[n][:, I]).ravel()
        mask = np.isfinite(e)
        error_rel2_J.append(np.linalg.norm(e[mask])/np.linalg.norm(jet[n][:, I].ravel()[mask]))

    error_rel2_DJ = []
    for n in N:
        I = np.arange(1, 3)
        e = (jet[n][:, I] - jet_gt[n][:, I]).ravel()
        mask = np.isfinite(e)
        error_rel2_DJ.append(np.linalg.norm(e[mask])/np.linalg.norm(jet[n][:, I].ravel()[mask]))

    error_rel2_D2J = []
    for n in N:
        I = np.arange(3, 7)
        e = (jet[n][:, I] - jet_gt[n][:, I]).ravel()
        mask = np.isfinite(e)
        error_rel2_D2J.append(np.linalg.norm(e[mask])/np.linalg.norm(jet[n][:, I].ravel()[mask]))

    plt.figure(figsize=(8, 8))
    plt.loglog(N, error_rel2_J, label=r'$\tau - T$', marker='*')
    plt.loglog(N, error_rel2_DJ, label=r'$\nabla\tau - \nabla T$', marker='*')
    plt.loglog(N, error_rel2_D2J, label=r'$\nabla^2\tau - \nabla^2 T$', marker='*')
    plt.legend()
    plt.xlabel('$n$')
    plt.ylabel(r'Relative $\ell_2$ error')
    plt.tight_layout()
    plt.savefig('out/error_rel2.png')
    plt.close()

    ########################################################################
    # fit Ch^p to each set of errors

    print('writing convergence rates to out/convergence_rates.txt')

    H = 1/np.array(N)

    p_J, C_J = np.polyfit(np.log(H), np.log(error_rel2_J), 1)
    p_DJ, C_DJ = np.polyfit(np.log(H), np.log(error_rel2_DJ), 1)
    p_D2J, C_D2J = np.polyfit(np.log(H), np.log(error_rel2_D2J), 1)

    with open('out/convergence_rates.txt', 'w') as f:
        print('J: p = %1.3f (C = %1.3f)' % (p_J, C_J), file=f)
        print('DJ: p = %1.3f (C = %1.3f)' % (p_DJ, C_DJ), file=f)
        print('D2J: p = %1.3f (C = %1.3f)' % (p_D2J, C_D2J), file=f)
