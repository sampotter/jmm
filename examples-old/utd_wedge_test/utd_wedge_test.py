import colorcet as cc
import logging
import matplotlib.pyplot as plt
import numpy as np
import pyvistaqt as pvqt

from meshplex import MeshTetra
from pathlib import Path

import jmm.bmesh

from jmm.plot import *

from utd_wedge_problem import UtdWedgeProblem

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    pv.set_plot_theme('dark')

    plt.ion()

    maxvol = 0.1
    maxvol *= 1.0
    n = 2/3
    sp = 5/np.sqrt(2)
    phip = 2*np.pi/3
    w = 10
    h = 2
    R = 1
    r = 1.5
    omega = 10000
    prob = UtdWedgeProblem(maxvol, n, sp, phip, w, h, R, r, omega)

    edge_lengths = MeshTetra(prob.verts, prob.cells).edge_lengths

    print('max vol = %g' % maxvol)

    print('avg. edge length = %g' % edge_lengths.mean())

    print('max edge length = %g' % edge_lengths.max())

    print('RMS error tau = %g (pt src)' %
          np.sqrt((np.sum(prob.error_pt_src_tau**2)/prob.domain.verts.size)))

    print('RMS error tau = %g (near refl)' %
          np.sqrt((np.sum(prob.error_near_refl_tau**2)/prob.domain.verts.size)))

    print('RMS error tau = %g (far refl)' %
          np.sqrt((np.sum(prob.error_far_refl_tau**2)/prob.domain.verts.size)))

    print('RMS error tau = %g (diff)' %
          np.sqrt((np.sum(prob.error_diff_tau**2)/prob.domain.verts.size)))

    ########################################################################
    # plot settings

    savefig = True
    outpath = Path(f'./wedge_plots/maxvol_{maxvol}')
    outpath.mkdir(parents=True, exist_ok=True)

    ########################################################################
    # debugging

    style = 'classic'
    figsize = (10, 4)

    # point source T and SPL
    with plt.style.context(style):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 2, 1)
        prob.plot_pt_src_T(fig, ax)
        ax = fig.add_subplot(1, 2, 2)
        prob.plot_pt_src_SPL(fig, ax)
        fig.tight_layout()
        fig.savefig(outpath/Path('pt_src_T_and_SPL.png'))

    # near refl T and SPL
    with plt.style.context(style):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 2, 1)
        prob.plot_near_refl_T(fig, ax)
        ax = fig.add_subplot(1, 2, 2)
        prob.plot_near_refl_SPL(fig, ax)
        fig.tight_layout()
        fig.savefig(outpath/Path('near_refl_T_and_SPL.png'))

    # far refl T and SPL
    with plt.style.context(style):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 2, 1)
        prob.plot_far_refl_T(fig, ax)
        ax = fig.add_subplot(1, 2, 2)
        prob.plot_far_refl_SPL(fig, ax)
        fig.tight_layout()
        fig.savefig(outpath/Path('far_refl_T_and_SPL.png'))

    # diff T and SPL
    with plt.style.context(style):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 2, 1)
        prob.plot_diff_T(fig, ax)
        ax = fig.add_subplot(1, 2, 2)
        prob.plot_diff_SPL(fig, ax)
        fig.tight_layout()
        fig.savefig(outpath/Path('diff_T_and_SPL.png'))

    # pt src T and gradT error
    with plt.style.context(style):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 2, 1)
        prob.plot_error_pt_src_tau(fig, ax)
        ax = fig.add_subplot(1, 2, 2)
        prob.plot_error_pt_src_grad_tau(fig, ax)
        fig.tight_layout()
        fig.savefig(outpath/Path('pt_src_T_and_grad_T_error.png'))

    # near refl T and gradT error
    with plt.style.context(style):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 2, 1)
        prob.plot_error_near_refl_tau(fig, ax)
        ax = fig.add_subplot(1, 2, 2)
        prob.plot_error_near_refl_grad_tau(fig, ax)
        fig.tight_layout()
        fig.savefig(outpath/Path('near_refl_T_and_grad_T_error.png'))

    # far refl T and gradT error
    with plt.style.context(style):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 2, 1)
        prob.plot_error_far_refl_tau(fig, ax)
        ax = fig.add_subplot(1, 2, 2)
        prob.plot_error_far_refl_grad_tau(fig, ax)
        fig.tight_layout()
        fig.savefig(outpath/Path('far_refl_T_and_grad_T_error.png'))

    # diff refl T and gradT error
    with plt.style.context(style):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 2, 1)
        prob.plot_error_diff_tau(fig, ax)
        ax = fig.add_subplot(1, 2, 2)
        prob.plot_error_diff_grad_tau(fig, ax)
        fig.tight_layout()
        fig.savefig(outpath/Path('diff_T_and_grad_T_error.png'))

    ########################################################################
    # NOTE: this isn't ready for prime-time... some issues:
    #
    # * It's slow. We need to spatially hash the meshgrid points and
    #   batch together the evaluations for each cell. If we wanted to
    #   speed it up even more past that point, we would need to write
    #   vectorized evaluation functions for the Bezier tetras.
    #
    # * It's wrong. Near a point source, this will end up computing
    #   the wrong values. We need to incorporate our local model of
    #   the point source into this since that's the function we're
    #   actually computing.
    #
    # * It's a little inflexible. We'll need to design a different
    #   version of this for each of the different quantities we're
    #   interested in investigating.
    #
    # * We're not getting the most out of it that we can. Once we
    #   define a continuous, finite element version of each of these
    #   functions, we should really be evaluating integrals to compute
    #   errors instead of just sums. And then we should really be
    #   using the error bounds described in M.J. Lai's book to see how
    #   well we're doing.

    # bmesh = jmm.bmesh.Bmesh33.from_eik3(prob.pt_src_field.eik)

    # padding = 1.1
    # ngrid = 101
    # z = 0
    # X, Y = np.meshgrid(
    #     np.linspace(-padding*w/2, padding*w/2, ngrid),
    #     np.linspace(-padding*w/2, padding*w/2, ngrid)
    # )

    # F = np.array([
    #     bmesh(np.array([x, y, z])) for x, y in zip(X.ravel(), Y.ravel())
    # ]).reshape(X.shape)

    # plt.figure()
    # plt.imshow(F)
    # plt.show()
