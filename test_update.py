import numpy as np
import sjs
import unittest

class TestF3(unittest.TestCase):

    def test_basic(self):
        cubic = sjs.Cubic([0, 0, 0, 0])
        xy0 = sjs.Dvec2(1, 0)
        xy1 = sjs.Dvec2(0, 1)
        xy = sjs.Dvec2(0, 0)
        h = 1
        data = sjs.F3Context(cubic, xy, xy0, xy1)
        self.assertAlmostEqual(data.F3(0), 1.0)
        self.assertAlmostEqual(data.F3(1), 1.0)
        self.assertAlmostEqual(data.F3(1/2), 1.0/np.sqrt(2))
        self.assertAlmostEqual(data.dF3_deta(0), -1.0)
        self.assertAlmostEqual(data.dF3_deta(1), 1.0)

    def test_hybrid(self):
        cubic = sjs.Cubic([0, 0, 0, 0])

        xy0 = sjs.Dvec2(1, 0)
        xy1 = sjs.Dvec2(0, 1)
        xy = sjs.Dvec2(0, 0)
        h = 1
        data = sjs.F3Context(cubic, xy, xy0, xy1)

        lam = sjs.hybrid(lambda lam: data.dF3_deta(lam), 0, 1)
        self.assertAlmostEqual(lam, 0.5)
        self.assertAlmostEqual(data.F3(lam), 1.0/np.sqrt(2))
        self.assertAlmostEqual(data.dF3_deta(lam), 0.0)

        data.cubic = sjs.Cubic([0, 1, 1, 1])
        lam = sjs.hybrid(lambda lam: data.dF3_deta(lam), 0, 1)
        self.assertAlmostEqual(lam, 0.0)
        self.assertAlmostEqual(data.F3(lam), 1.0)
        self.assertAlmostEqual(data.dF3_deta(lam), 0.0)

        data.cubic = sjs.Cubic([1, 0, -1, -1])
        lam = sjs.hybrid(lambda lam: data.dF3_deta(lam), 0, 1)
        self.assertAlmostEqual(lam, 1.0)
        self.assertAlmostEqual(data.F3(lam), 1.0)
        self.assertAlmostEqual(data.dF3_deta(lam), 0.0)

        data.cubic = sjs.Cubic([1, 1, -1, 1])
        lam = sjs.hybrid(lambda lam: data.dF3_deta(lam), 0, 1)
        self.assertAlmostEqual(lam, 0.5)
        self.assertAlmostEqual(data.dF3_deta(lam), 0.0)

if __name__ == '__main__':
    unittest.main()
