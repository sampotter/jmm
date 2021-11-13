import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import vtk

import jmm.utd

def get_D(omega, nplot=1000):
    t0 = 10000
    xsrc = t0*np.array([1, -1, 0])*np.sqrt(2)

    c = 340.3

    k = omega/c

    alpha = np.pi/2
    no = np.array([0, -1, 0], dtype=np.float64)
    e = np.array([0, 0, -1], dtype=np.float64)

    theta = np.linspace(0, 3*np.pi/2, nplot)
    x = np.cos(-theta)
    y = np.sin(-theta)
    z = np.zeros(nplot)

    s = np.array([x, y, z]).T
    thetap = np.pi/4
    sp = np.outer(np.ones(nplot), [np.cos(thetap), np.sin(thetap), 0])

    t = 1*np.ones(nplot)

    r = np.linalg.norm(xsrc)
    grad = -xsrc/r
    hess = (np.eye(3) - np.outer(grad, grad))/r
    hess = np.outer(np.ones(nplot), hess).reshape(-1, 3, 3)

    D = jmm.utd.D_from_geometry(k, alpha, no, e, s, sp, t, hess)

    return theta, D

if __name__ == '__main__':
    plt.ion()

    pv.set_plot_theme('document')

    theta, D_1_000 = get_D(1000, nplot=5000)
    theta, D_3_000 = get_D(3000, nplot=5000)
    theta, D_10_000 = get_D(10000, nplot=5000)

    plt.figure()
    plt.plot(theta, 20*np.log10(abs(D_1_000)), label=r'$\omega = 1\mathrm{kHz}$')
    plt.plot(theta, 20*np.log10(abs(D_3_000)), label=r'$\omega = 3\mathrm{kHz}$')
    plt.plot(theta, 20*np.log10(abs(D_10_000)), label=r'$\omega = 10\mathrm{kHz}$')
    plt.xlabel(r'$\theta_{\mathrm{out}}$')
    plt.ylabel(r'$|D| \mathrm{dB}$')
    plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2])
    plt.gca().set_xticklabels(['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$'])
    plt.legend()
    plt.plot()

    verts = [[0, 10, -5],
             [0, 10, 5],
             [0, 0, -5],
             [0, 0, 5],
             [10, 0, -5],
             [10, 0, 5]]
    quads = [[0, 1, 3, 2],
             [2, 3, 5, 4]]
