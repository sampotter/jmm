#!/usr/bin/env python


'''This is a version of basic_solve.c written in Python, intended to
demonstrate the Python API's use.'''


import argparse
import colorcet as cc
import jmm
import matplotlib.pyplot as plt
import numpy as np
import os
import pyvista as pv
import pyvistaqt as pvqt
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, help='''
path to directory from which tetrahedron mesh will be read and to
which output will be written''')
    parser.add_argument('--scene', type=str)
    parser.add_argument('--mask', type=str, help='''
mask to use when coloring scatter plots (only "boundary" now)''')
    parser.add_argument('--indsrc', type=int, help='''
index to vertex with location of point source''')
    parser.add_argument('--r0', type=float, help='''
factoring radius (TODO: not currently being used)''', default=None)
    args = parser.parse_args()

    path = args.path
    scene = args.scene
    indsrc = args.indsrc
    r0 = args.r0

    if args.mask is not None:
        if args.mask in {'boundary'}:
            mask = args.mask
        else:
            raise Exception(
                'mask should be "boundary" (got "%s")' % args.mask)
    else:
        mask = args.mask

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

    if indsrc is None:
        print('- placing point source at vertex closest to (0, 0, 0)')
        xsrc = np.zeros(3)
        R = verts - xsrc
        tau = np.sqrt(np.sum(R**2, axis=1))
        indsrc = np.argmin(tau)
    else:
        xsrc = verts[indsrc]
        R = verts - xsrc
        tau = np.sqrt(np.sum(R**2, axis=1))
    print('- point source is at index %d: %s' % (indsrc, tuple(xsrc)))

    Dtau = np.empty((tau.size, 3), dtype=tau.dtype)
    nz = tau != 0
    Dtau[nz] = R[nz]/tau[nz].reshape(tau[nz].shape[0], 1)
    Dtau[~nz] = np.nan
    del nz

    if r0 is None:
        eik.add_trial(indsrc, jmm.Jet3(0, np.nan, np.nan, np.nan))
    else:
        num_fac = 0
        for l in np.where(tau < r0)[0]:
            eik.add_valid(l, jmm.Jet3(tau[l], *Dtau[l]))
            num_fac += 1
        for l in np.where(eik.state == jmm.State.Valid.value)[0]:
            for m in mesh.vv(l):
                if eik.state[m] == jmm.State.Valid.value:
                    continue
                eik.add_trial(m, jmm.Jet3(tau[m], *Dtau[m]))
                num_fac += 1
        print('- number of nodes inside factoring radius: %d' % num_fac)

    raise Exception()

    eik.solve()

    T = np.array([jet[0] for jet in eik.jet])
    DT = np.array([(jet[1], jet[2], jet[3]) for jet in eik.jet])

    E = T - tau
    nn = ~np.isnan(DT).any(1)
    ED = np.empty_like(E)
    ED[nn] = np.sqrt(np.sum((DT[nn] - Dtau[nn])**2, axis=1))
    ED[~nn] = np.nan
    angle = np.empty_like(E)
    dot = (DT[nn]*Dtau[nn]).sum(1)
    if (abs(dot).max() > 1 + 1e-15).any():
        print('* WARNING: found weird dot product: %g' % abs(dot).max())
    dot = np.maximum(-1, np.minimum(1, dot))
    angle[nn] = np.rad2deg(np.arccos(dot))
    angle[~nn] = np.nan
    eT = np.linalg.norm(E)/np.linalg.norm(T)
    print('- l2 error (p = 0): %g' % eT)

    if mask == 'boundary':
        import meshplex
        mesh_tetra = meshplex.MeshTetra(verts, cells)
        mesh_tetra.mark_boundary()
        mask = mesh_tetra.is_boundary_point
    else:
        mask = None

    plt.figure(figsize=(12, 4))

    plt.subplot(1, 3, 1)
    if mask is None:
        plt.scatter(T, E, s=1, c='k')
    else:
        plt.scatter(T[~mask], E[~mask], s=1, c='k', zorder=1)
        plt.scatter(T[mask], E[mask], s=1, c='r', zorder=2)
    plt.xlabel(r'$\tau(x) = \|x\|$')
    plt.ylabel(r'$T(x) - \tau(x)$')

    plt.subplot(1, 3, 2)
    if mask is None:
        plt.scatter(T, ED, s=1, c='k')
    else:
        plt.scatter(T[~mask], ED[~mask], s=1, c='k', zorder=1)
        plt.scatter(T[mask], ED[mask], s=1, c='r', zorder=2)
    plt.xlabel(r'$\tau(x) = \|x\|$')
    plt.ylabel(r'$\|\nabla T(x) - \nabla \tau(x)\|$')

    plt.subplot(1, 3, 3)
    if mask is None:
        plt.scatter(T, angle, s=1, c='k')
    else:
        plt.scatter(T[~mask], angle[~mask], s=1, c='k', zorder=1)
        plt.scatter(T[mask], angle[mask], s=1, c='r', zorder=2)
    plt.xlabel(r'$\tau(x) = \|x\|$')
    plt.ylabel(r'$\angle (\nabla T, \nabla \tau)$ [Deg.]')
    plt.tight_layout()
    plt.savefig('%s/hists.pdf' % path)

    plt.figure(figsize=(10, 8))
    if mask is None:
        plt.scatter(T, np.maximum(2.2e-16, abs(E)), s=0.25, c='k')
    else:
        plt.scatter(T[~mask], np.maximum(2.2e-16, abs(E[~mask])), s=0.25, c='k', zorder=1)
        plt.scatter(T[mask], np.maximum(2.2e-16, abs(E[mask])), s=0.25, c='r', zorder=2)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig('%s/log_abs_error_hist.pdf' % path)
