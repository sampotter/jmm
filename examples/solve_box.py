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


r0 = 0.1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=str)
    parser.add_argument('--mask', type=str, help='''
mask to use when coloring scatter plots (one of: "full", "boundary")''')
    args = parser.parse_args()

    root = args.root

    if args.mask is not None:
        if args.mask in {'full', 'boundary'}:
            mask = args.mask
        else:
            raise Exception(
                'mask should be "full" or "boundary" (got "%s")' % args.mask)

    verts_bin_path = root + '_verts.bin'
    cells_bin_path = root + '_cells.bin'

    print('- reading verts from %s' % verts_bin_path)
    print('- reading cells from %s' % cells_bin_path)

    verts = np.fromfile(verts_bin_path, np.float64)
    verts = verts.reshape(verts.size//3, 3)

    cells = np.fromfile(cells_bin_path, np.uintp)
    cells = cells.reshape(cells.size//4, 4)

    mesh = jmm.Mesh3(verts, cells)

    eik = jmm.Eik3(mesh)

    xsrc = np.zeros(3)
    R = verts - xsrc
    tau = np.sqrt(np.sum(R**2, axis=1))

    Dtau = np.empty((tau.size, 3), dtype=tau.dtype)
    nz = tau != 0
    Dtau[nz] = R[nz]/tau[nz].reshape(tau[nz].shape[0], 1)
    Dtau[~nz] = np.nan
    del nz

    indsrc = np.argmin(tau)
    print('- point source is at index %d: %s' % (indsrc, tuple(xsrc)))

    eik.add_trial(indsrc, jmm.Jet3(0, np.nan, np.nan, np.nan))

    eik.solve()

    print('- number of full updates: %d' % eik.num_full_updates)

    T = np.array([jet[0] for jet in eik.jet])
    DT = np.array([(jet[1], jet[2], jet[3]) for jet in eik.jet])

    E = T - tau
    nn = ~np.isnan(DT).any(1)
    ED = np.empty_like(E)
    ED[nn] = np.sqrt(np.sum((DT[nn] - Dtau[nn])**2, axis=1))
    ED[~nn] = np.nan
    angle = np.empty_like(E)
    angle[nn] = np.rad2deg(np.arccos((DT[nn]*Dtau[nn]).sum(1)))
    angle[~nn] = np.nan

    eT = np.linalg.norm(E)/np.linalg.norm(T)
    print('- l2 error (p = 0): %g' % eT)

    if mask == 'full':
        mask = eik.full_update
    elif mask == 'boundary':
        import meshplex
        mesh_tetra = meshplex.MeshTetra(verts, cells)
        mesh_tetra.mark_boundary()
        mask = mesh_tetra.is_boundary_point
    else:
        mask = None

    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    if mask is None:
        plt.scatter(T, angle, s=1, c='k')
    else:
        plt.scatter(T[~mask], angle[~mask], s=1, c='k', zorder=1)
        plt.scatter(T[mask], angle[mask], s=1, c='r', zorder=2)
    plt.xlabel(r'$\tau(x) = \|x\|$')
    plt.ylabel(r'$\angle (\nabla T, \nabla \tau)$ [Deg.]')

    plt.subplot(1, 2, 2)
    if mask is None:
        plt.scatter(T, angle, s=1, c='k')
    else:
        plt.scatter(T[~mask], E[~mask], s=1, c='k', zorder=1)
        plt.scatter(T[mask], E[mask], s=1, c='r', zorder=2)
    plt.xlabel(r'$\tau(x) = \|x\|$')
    plt.ylabel(r'$T(x) - \tau(x)$')

    plt.tight_layout()
    plt.savefig('%s_hists.pdf' % root)
