import itertools as it
import numpy as np
import jmm
import unittest

class TestBicubic(unittest.TestCase):

    def test_ctor(self):
        for _ in range(10):
            data = np.random.randn(4, 4)
            bicubic = jmm.Bicubic(data)
            for (i, lam), (j, mu) in it.product(enumerate([0.0,1.0]), repeat=2):
                self.assertAlmostEqual(bicubic.f(lam, mu), data[i, j])
                self.assertAlmostEqual(bicubic.fx(lam, mu), data[2 + i, j])
                self.assertAlmostEqual(bicubic.fy(lam, mu), data[i, 2 + j])
                self.assertAlmostEqual(bicubic.fxy(lam, mu), data[2 + i, 2 + j])

    def test_get_f_on_edge(self):
        for _ in range(10):
            data = np.random.randn(4, 4)
            bicubic = jmm.Bicubic(data)

            cubic = bicubic.get_f_on_edge(jmm.BicubicVariable.Lambda, 0)
            self.assertAlmostEqual(cubic.f(0), data[0, 0])
            self.assertAlmostEqual(cubic.f(1), data[1, 0])
            self.assertAlmostEqual(cubic.df(0), data[2, 0])
            self.assertAlmostEqual(cubic.df(1), data[3, 0])
            for _ in range(10):
                lam = np.random.rand()
                self.assertAlmostEqual(cubic.f(lam), bicubic.f(lam, 0))

            cubic = bicubic.get_f_on_edge(jmm.BicubicVariable.Lambda, 1)
            self.assertAlmostEqual(cubic.f(0), data[0, 1])
            self.assertAlmostEqual(cubic.f(1), data[1, 1])
            self.assertAlmostEqual(cubic.df(0), data[2, 1])
            self.assertAlmostEqual(cubic.df(1), data[3, 1])
            for _ in range(10):
                lam = np.random.rand()
                self.assertAlmostEqual(cubic.f(lam), bicubic.f(lam, 1))

            cubic = bicubic.get_f_on_edge(jmm.BicubicVariable.Mu, 0)
            self.assertAlmostEqual(cubic.f(0), data[0, 0])
            self.assertAlmostEqual(cubic.f(1), data[0, 1])
            self.assertAlmostEqual(cubic.df(0), data[0, 2])
            self.assertAlmostEqual(cubic.df(1), data[0, 3])
            for _ in range(10):
                mu = np.random.rand()
                self.assertAlmostEqual(cubic.f(mu), bicubic.f(0, mu))

            cubic = bicubic.get_f_on_edge(jmm.BicubicVariable.Mu, 1)
            self.assertAlmostEqual(cubic.f(0), data[1, 0])
            self.assertAlmostEqual(cubic.f(1), data[1, 1])
            self.assertAlmostEqual(cubic.df(0), data[1, 2])
            self.assertAlmostEqual(cubic.df(1), data[1, 3])
            for _ in range(10):
                mu = np.random.rand()
                self.assertAlmostEqual(cubic.f(mu), bicubic.f(1, mu))

    def test_get_fx_on_edge(self):
        for _ in range(10):
            data = np.random.randn(4, 4)
            bicubic = jmm.Bicubic(data)

            cubic = bicubic.get_fx_on_edge(jmm.BicubicVariable.Lambda, 0)
            self.assertAlmostEqual(cubic.f(0), data[2, 0])
            self.assertAlmostEqual(cubic.f(1), data[3, 0])
            for _ in range(10):
                lam = np.random.rand()
                self.assertAlmostEqual(cubic.f(lam), bicubic.fx(lam, 0))

            cubic = bicubic.get_fx_on_edge(jmm.BicubicVariable.Lambda, 1)
            self.assertAlmostEqual(cubic.f(0), data[2, 1])
            self.assertAlmostEqual(cubic.f(1), data[3, 1])
            for _ in range(10):
                lam = np.random.rand()
                self.assertAlmostEqual(cubic.f(lam), bicubic.fx(lam, 1))

            cubic = bicubic.get_fx_on_edge(jmm.BicubicVariable.Mu, 0)
            self.assertAlmostEqual(cubic.f(0), data[2, 0])
            self.assertAlmostEqual(cubic.f(1), data[2, 1])
            for _ in range(10):
                mu = np.random.rand()
                self.assertAlmostEqual(cubic.f(mu), bicubic.fx(0, mu))

            cubic = bicubic.get_fx_on_edge(jmm.BicubicVariable.Mu, 1)
            self.assertAlmostEqual(cubic.f(0), data[3, 0])
            self.assertAlmostEqual(cubic.f(1), data[3, 1])
            for _ in range(10):
                mu = np.random.rand()
                self.assertAlmostEqual(cubic.f(mu), bicubic.fx(1, mu))

    def test_get_fy_on_edge(self):
        for _ in range(10):
            data = np.random.randn(4, 4)
            bicubic = jmm.Bicubic(data)

            cubic = bicubic.get_fy_on_edge(jmm.BicubicVariable.Lambda, 0)
            self.assertAlmostEqual(cubic.f(0), data[0, 2])
            self.assertAlmostEqual(cubic.f(1), data[1, 2])
            for _ in range(10):
                lam = np.random.rand()
                self.assertAlmostEqual(cubic.f(lam), bicubic.fy(lam, 0))

            cubic = bicubic.get_fy_on_edge(jmm.BicubicVariable.Lambda, 1)
            self.assertAlmostEqual(cubic.f(0), data[0, 3])
            self.assertAlmostEqual(cubic.f(1), data[1, 3])
            for _ in range(10):
                lam = np.random.rand()
                self.assertAlmostEqual(cubic.f(lam), bicubic.fy(lam, 1))

            cubic = bicubic.get_fy_on_edge(jmm.BicubicVariable.Mu, 0)
            self.assertAlmostEqual(cubic.f(0), data[0, 2])
            self.assertAlmostEqual(cubic.f(1), data[0, 3])
            for _ in range(10):
                mu = np.random.rand()
                self.assertAlmostEqual(cubic.f(mu), bicubic.fy(0, mu))

            cubic = bicubic.get_fy_on_edge(jmm.BicubicVariable.Mu, 1)
            self.assertAlmostEqual(cubic.f(0), data[1, 2])
            self.assertAlmostEqual(cubic.f(1), data[1, 3])
            for _ in range(10):
                mu = np.random.rand()
                self.assertAlmostEqual(cubic.f(mu), bicubic.fy(1, mu))

    def test_interpolate_fxy_at_verts(self):
        fx = np.zeros(4)
        fy = np.zeros(4)
        h = 1
        fxy = jmm.interpolate_fxy_at_verts(fx, fy, h)
        self.assertAlmostEqual(fxy[0], 0)
        self.assertAlmostEqual(fxy[1], 0)
        self.assertAlmostEqual(fxy[2], 0)
        self.assertAlmostEqual(fxy[3], 0)

        A = np.array([
            [ 3, -1,  3, -1],
            [ 3,  3, -1, -1],
            [-1, -1,  3,  3],
            [-1,  3, -1,  3]
        ])/4
        Bx = np.array([
            [ 0, 0,  0, 0],
            [ 0, 0, -1, 1],
            [-1, 1,  0, 0],
            [ 0, 0,  0, 0]
        ])
        By = np.array([
            [-1, 1,  0, 0],
            [ 0, 0,  0, 0],
            [ 0, 0,  0, 0],
            [ 0, 0, -1, 1]
        ])

        for _ in range(10):
            fx = np.random.randn(4)
            fy = np.random.randn(4)
            h = np.random.rand()
            fxy_gt = (A@Bx@fx + A@By@fy)/h
            fxy = jmm.interpolate_fxy_at_verts(fx, fy, h)
            self.assertAlmostEqual(fxy[0], fxy_gt[0])
            self.assertAlmostEqual(fxy[1], fxy_gt[1])
            self.assertAlmostEqual(fxy[2], fxy_gt[2])
            self.assertAlmostEqual(fxy[3], fxy_gt[3])

        for _ in range(10):
            a = np.random.randn(2)
            h = np.random.rand()
            f = np.array([a[0], a[0], a[1], a[1]])
            fx = (a[1] - a[0])*np.ones(4)
            fy = np.zeros(4)
            fxy = jmm.interpolate_fxy_at_verts(fx, fy, h)
            self.assertAlmostEqual(fxy[0], 0.0)
            self.assertAlmostEqual(fxy[1], 0.0)
            self.assertAlmostEqual(fxy[2], 0.0)
            self.assertAlmostEqual(fxy[3], 0.0)

        for _ in range(10):
            b = np.random.randn(2)
            h = np.random.rand()
            f = np.array([b[0], b[1], b[0], b[1]])
            fx = np.zeros(4)
            fy = (b[1] - b[0])*np.ones(4)
            fxy = jmm.interpolate_fxy_at_verts(fx, fy, h)
            self.assertAlmostEqual(fxy[0], 0.0)
            self.assertAlmostEqual(fxy[1], 0.0)
            self.assertAlmostEqual(fxy[2], 0.0)
            self.assertAlmostEqual(fxy[3], 0.0)

        for _ in range(10):
            a = np.random.randn(3)
            h = np.random.rand()
            f = np.array([
                a[2],
                a[0] + a[2],
                a[1] + a[2],
                a[0] + a[1] + a[2],
            ])
            fx = a[0]*np.ones(4)
            fy = a[1]*np.ones(4)
            fxy = jmm.interpolate_fxy_at_verts(fx, fy, h)
            self.assertAlmostEqual(fxy[0], 0.0)
            self.assertAlmostEqual(fxy[1], 0.0)
            self.assertAlmostEqual(fxy[2], 0.0)
            self.assertAlmostEqual(fxy[3], 0.0)

        # for _ in range(10):
        #     f = np.random.randn(4)
        #     h = np.random.rand()
        #     # fx = np.array([
        #     #     f[1] - f[0],
        #     #     f[1] - f[0],
        #     #     f[3] - f[2],
        #     #     f[3] - f[2]
        #     # ])/h
        #     # fy = np.array([
        #     #     f[2] - f[0],
        #     #     f[3] - f[1],
        #     #     f[2] - f[0],
        #     #     f[3] - f[1]
        #     # ])/h
        #     fx = h*np.array([
        #         f[1] - f[0],
        #         f[1] - f[0],
        #         (1 - h)*(f[1] - f[0]) - h*(f[3] - f[2]),
        #         (1 - h)*(f[1] - f[0]) - h*(f[3] - f[2])
        #     ])
        #     fy = h*np.array([
        #         f[2] - f[0],
        #         (1 - h)*(f[2] - f[0]) - h*(f[3] - f[1]),
        #         f[2] - f[0],
        #         (1 - h)*(f[2] - f[0]) - h*(f[3] - f[1])
        #     ])
        #     fxy = jmm.interpolate_fxy_at_verts(fx, fy, h)
        #     fxy_gt = h**2*(f[0] - f[1] - f[2] + f[3])
        #     print(fxy, h)
        #     self.assertAlmostEqual(fxy[0], fxy_gt)
        #     self.assertAlmostEqual(fxy[1], fxy_gt)
        #     self.assertAlmostEqual(fxy[2], fxy_gt)
        #     self.assertAlmostEqual(fxy[3], fxy_gt)

if __name__ == '__main__':
    unittest.main()
