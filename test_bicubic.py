import itertools as it
import numpy as np
import sjs_eik
import unittest

class TestBicubic(unittest.TestCase):

    def test_ctor(self):
        for _ in range(10):
            data = np.random.randn(4, 4)
            bicubic = sjs_eik.Bicubic(data)
            for (i, lam), (j, mu) in it.product(enumerate([0.0,1.0]), repeat=2):
                self.assertTrue(np.isclose(bicubic.f(lam, mu), data[i, j]))
                self.assertTrue(np.isclose(bicubic.fx(lam, mu), data[2 + i, j]))
                self.assertTrue(np.isclose(bicubic.fy(lam, mu), data[i, 2 + j]))
                self.assertTrue(np.isclose(bicubic.fxy(lam, mu), data[2 + i, 2 + j]))

if __name__ == '__main__':
    unittest.main()
