import numpy as np
import jmm
import unittest

class TestCubic(unittest.TestCase):

    def test_ctor(self):
        for _ in range(10):
            data = np.random.randn(4)
            cubic = jmm.Cubic(data)
            self.assertAlmostEqual(cubic.f(0), data[0])
            self.assertAlmostEqual(cubic.f(1), data[1])
            self.assertAlmostEqual(cubic.df(0), data[2])
            self.assertAlmostEqual(cubic.df(1), data[3])

    def test_reverse_on_unit_interval(self):
        for _ in range(10):
            data = np.random.randn(4)
            cubic = jmm.Cubic(data)
            reversed_cubic = jmm.Cubic(data)
            reversed_cubic.reverse_on_unit_interval()
            for _ in range(10):
                lam = np.random.rand()
                self.assertAlmostEqual(cubic.f(1 - lam), reversed_cubic.f(lam))
                self.assertAlmostEqual(cubic.df(1 - lam), -reversed_cubic.df(lam))

if __name__ == '__main__':
    unittest.main()
