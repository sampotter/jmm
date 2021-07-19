import numpy as np
import pyvistaqt as pvqt

from utd_wedge_problem import UtdWedgeProblem

if __name__ == '__main__':
    maxvol = 0.1
    n = 0.1
    sp = 5
    phip = np.pi/4
    w = 10
    h = 2
    prob = UtdWedgeProblem(maxvol, n, sp, phip, w, h)

    plotter = pvqt.BackgroundPlotter()
    prob.plot_mesh(plotter, show_edges=True)
