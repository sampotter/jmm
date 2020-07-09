import numpy as np
import jmm
import unittest

class TestHybrid(unittest.TestCase):

    def test_line(self):
        for _ in range(10):
            t0 = np.random.rand()
            f = lambda t: t - t0
            t = jmm.hybrid(f, 0, 1)
            self.assertAlmostEqual(t, t0)

    def test_cubic(self):
        for _ in range(10):
            t0 = np.random.rand()
            f = lambda t: (t - t0)**3
            t = jmm.hybrid(f, 0, 1)
            self.assertAlmostEqual(t, t0)

    def test_quadratic(self):
        for _ in range(10):
            t0 = np.random.rand()
            f = lambda t: (t/t0)**2 - 1
            t = jmm.hybrid(f, 0, 1)
            self.assertAlmostEqual(t, t0)
            t = jmm.hybrid(f, -1, 0)
            self.assertAlmostEqual(t, -t0)

if __name__ == '__main__':
    main.unittest()
