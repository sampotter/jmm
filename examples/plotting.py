import colorcet as cc
import numpy as np
import pyvista as pv

import jmm.defs

def plot_x(plotter, x, scale=1, sphere_radius=0.05, **kwargs):
    plotter.add_mesh(
        pv.Sphere(scale*sphere_radius, x),
        **kwargs)

def plot_line(plotter, x, y, scale=1, color='white', sphere_radius=0.05):
    plot_x(plotter, x, scale=scale, color=color)
    plot_x(plotter, y, scale=scale, color=color)
    xm = (x + y)/2
    xd = x - y
    d = np.linalg.norm(xd)
    xd /= d
    r = scale*sphere_radius
    plotter.add_mesh(
        pv.Cylinder(xm, xd, scale*sphere_radius, d, capping=False),
        color=color)
