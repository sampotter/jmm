import jmm
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt


def F(x, Xt, T, DT):
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

    # Visualize the result...
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
