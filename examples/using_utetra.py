import jmm
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt


def F(x, Xt, T, DT):
    '''Evaluate the solution mapping for a tetrahedron update assuming a
constant (c = 1 = s) speed of sound.

    - x:    A length 3 array indicating the location of the evaluation point.

    - Xt:   A 3x3 array where each row (i.e., Xt[i] for i = 0, 1, 2) is
            one of the vertices of the base of the update tetrahedron.

    - T:    A length 3 array where T[i] is the eikonal at Xt[i].

    - DT:   A 3x3 array where row DT[i] is the gradient of the eikonal
            at point Xt[i].

    Returns a 3-tuple where:

    - The first entry is a length 2 array containing the minimizing
      argument (lam[0], lam[1]) such that Xt[0] + lam@(Xt[1:] - Xt[0])
      is the starting point of the optimal ray leading into x.

    - The second entry is the optimal eikonal value.

    - The third entry is a length 3 array with the optimal eikonal
      gradient.

    '''
    lam0 = np.array([1, 1])/3
    utetra = jmm.UpdateTetra(x, Xt, T, DT)
    utetra.set_lambda(lam0)
    utetra.solve()
    lam = utetra.get_lambda()
    jet = utetra.get_jet()
    return lam, jet.f, np.array([jet.fx, jet.fy, jet.fz])


if __name__ == '__main__':
    x = np.ones(3)
    Xt = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]], dtype=np.float64)
    T = np.sqrt(2)*np.ones(3)
    DT = Xt/np.sqrt(np.sum(Xt**2, axis=1)).reshape(3, 1)

    # For this test problem, should have lam == (1/3, 1/3), Tx ~=
    # sqrt(3), and DTx == (1, 1, 1)/sqrt(3). Because of approximation,
    # Tx will be a little bigger than sqrt(3).
    lam, Tx, DTx = F(x, Xt, T, DT)

    xlam = Xt[0] + lam@(Xt[1:] - Xt[0])

    # Visualize the result using PyVista...
    r = 0.04
    plotter = pvqt.BackgroundPlotter()
    for i in range(3):
        plotter.add_mesh(pv.Sphere(r, Xt[i]), color='white')
        dx = Xt[(i+1)%3] - Xt[i]
        xm = (Xt[(i+1)%3] + Xt[i])/2
        h = np.linalg.norm(dx)
        dx /= h
        plotter.add_mesh(pv.Cylinder(xm, dx, r/3, h), color='white')
    update_base = pv.PolyData(Xt, np.array([3, 0, 1, 2]))
    plotter.add_mesh(pv.Sphere(r, x), color='red')
    plotter.add_mesh(pv.Sphere(r, xlam), color='red')
    plotter.add_mesh(update_base, color='white', opacity=0.5)
    hlam = np.linalg.norm(x - xlam)
    plotter.add_mesh(pv.Cylinder((x + xlam)/2, (x - xlam)/hlam, r/3, hlam),
                     color='red')
    plotter.show()
