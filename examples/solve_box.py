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
    args = parser.parse_args()

    root = args.root

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

    valid_inds = [indsrc]
    eik.add_valid(indsrc, jmm.Jet3(tau[indsrc], *Dtau[indsrc]))

    # Start by making all neighbors of the point source VALID.
    for i in mesh.vv(indsrc):
        eik.add_valid(i, jmm.Jet3(tau[i], *Dtau[i]))
        if i not in valid_inds:
            valid_inds.append(i)

    # Make all remaining nodes in the factored ball VALID.
    for i in np.where(np.sqrt(np.sum(verts**2, axis=1)) <= r0)[0]:
        if eik.is_valid(i):
            continue
        eik.add_valid(i, jmm.Jet3(Dtau[i], *Dtau[i]))
        if i not in valid_inds:
            valid_inds.append(i)

    # Make all FAR neighbors of VALID nodes TRIAL.
    for i in valid_inds:
        for j in mesh.vv(i):
            if eik.is_far(j):
                eik.add_trial(j, jmm.Jet3(tau[i], *Dtau[i]))

    eik.solve()

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

    plt.figure(figsize=(9, 5))
    plt.subplot(1, 2, 1)
    plt.scatter(T, angle, s=2, c='k')
    plt.xlabel(r'$\tau(x) = \|x\|$')
    plt.ylabel(r'$\angle (\nabla T, \nabla \tau)$ [Deg.]')
    plt.subplot(1, 2, 2)
    plt.scatter(T, E, s=2, c='k')
    plt.xlabel(r'$\tau(x) = \|x\|$')
    plt.ylabel(r'$|T(x) - \tau(x)|$')
    plt.tight_layout()
    plt.savefig('%s_hists.pdf' % root)
