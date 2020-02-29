import itertools as it
import numpy as np
import sjs_eik as sjs
import unittest

class TestBicubic(unittest.TestCase):

    def test_ctor(self):
        for _ in range(10):
            data = np.random.randn(4, 4)
            bicubic = sjs.Bicubic(data)
            for (i, lam), (j, mu) in it.product(enumerate([0.0,1.0]), repeat=2):
                self.assertAlmostEqual(bicubic.f(lam, mu), data[i, j])
                self.assertAlmostEqual(bicubic.fx(lam, mu), data[2 + i, j])
                self.assertAlmostEqual(bicubic.fy(lam, mu), data[i, 2 + j])
                self.assertAlmostEqual(bicubic.fxy(lam, mu), data[2 + i, 2 + j])

    def test_restrict(self):
        for _ in range(10):
            data = np.random.randn(4, 4)
            bicubic = sjs.Bicubic(data)

            cubic = bicubic.restrict(sjs.BicubicVariable.Lambda, 0)
            self.assertAlmostEqual(cubic.f(0), data[0, 0])
            self.assertAlmostEqual(cubic.f(1), data[1, 0])
            self.assertAlmostEqual(cubic.df(0), data[2, 0])
            self.assertAlmostEqual(cubic.df(1), data[3, 0])

            cubic = bicubic.restrict(sjs.BicubicVariable.Lambda, 1)
            self.assertAlmostEqual(cubic.f(0), data[0, 1])
            self.assertAlmostEqual(cubic.f(1), data[1, 1])
            self.assertAlmostEqual(cubic.df(0), data[2, 1])
            self.assertAlmostEqual(cubic.df(1), data[3, 1])

            cubic = bicubic.restrict(sjs.BicubicVariable.Mu, 0)
            self.assertAlmostEqual(cubic.f(0), data[0, 0])
            self.assertAlmostEqual(cubic.f(1), data[0, 1])
            self.assertAlmostEqual(cubic.df(0), data[0, 2])
            self.assertAlmostEqual(cubic.df(1), data[0, 3])

            cubic = bicubic.restrict(sjs.BicubicVariable.Mu, 1)
            self.assertAlmostEqual(cubic.f(0), data[1, 0])
            self.assertAlmostEqual(cubic.f(1), data[1, 1])
            self.assertAlmostEqual(cubic.df(0), data[1, 2])
            self.assertAlmostEqual(cubic.df(1), data[1, 3])

if __name__ == '__main__':
    unittest.main()
