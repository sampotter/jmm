import colorcet as cc
import logging
import itertools as it
import json
import meshplex
import numpy as np
import os
import pyvista as pv
import pyvistaqt as pvqt
import skimage.morphology
import time

from PIL import Image

from skimage.measure import marching_cubes

from jmm.grid import Grid3
from utd_wedge_problem import UtdWedgeProblem

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

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

    # Transfer solution

    padding_scale = 1.1
    nw = 51
    x, y, z = np.meshgrid(
        np.linspace(-padding_scale*w/2, padding_scale*w/2, nw),
        np.linspace(-padding_scale*w/2, padding_scale*w/2, nw),
        np.linspace(-padding_scale*h/2, padding_scale*h/2, int(np.ceil((h/w)*nw)))
    )

    def grid_eik3(eik):
        bmesh = jmm.bmesh.Bmesh33.from_eik3(eik)
        return np.array([
            bmesh(np.array([x[i, j, k], y[i, j, k], z[i, j, k]]))
            for i, j, k in it.product(*(range(_) for _ in x.shape))
        ]).reshape(x.shape)

    # Make movie

    levels = np.linspace(0.75, 20, 30*24)[1:-1]

    def get_level_set_poly_data(level):
        selem = np.ones((3, 3, 3), dtype=np.intc)
        mask = skimage.morphology.binary_erosion(np.isfinite(T), selem)
        V, F = marching_cubes(T, level=level, spacing=(h, h, h),
                              allow_degenerate=False, mask=mask)[:2]
        return pv.PolyData(
            V + pmin,
            np.concatenate([
                3*np.ones((F.shape[0], 1)), F], axis=1).astype(F.dtype))

    plotter = pvqt.BackgroundPlotter()
    plotter.set_position((-20, 35, -20))
    plotter.set_viewup((0, 1, 0))
    plotter.set_focus((0, 0, 0))


    plotter.add_mesh(surf_mesh, color='white', opacity=0.3)
    plotter.add_mesh(get_level_set_poly_data(levels[0]), color='white')

    os.mkdir('frames')

    i = 0

    level = levels[i]
    plotter.clear()
    plotter.add_mesh(surf_mesh, color='white', opacity=0.3)
    plotter.add_mesh(get_level_set_poly_data(level), color='white')
    i += 1

    Image.fromarray(plotter.image).save('frames/%03d.png' % i)




    print('finished')
