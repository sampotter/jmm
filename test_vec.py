import numpy as np
import sjs_eik as sjs
import unittest

class TestDvec2(unittest.TestCase):

    def test_ctor(self):
        v = sjs.Dvec2(1, 2)
        self.assertEqual(v.x, 1)
        self.assertEqual(v.y, 2)

    def test_ccomb(self):
        u, v = sjs.Dvec2(2, 0), sjs.Dvec2(1, 3)
        w = sjs.ccomb(u, v, 0.1)
        self.assertAlmostEqual(w.x, 1.9)
        self.assertAlmostEqual(w.y, 0.3)

    def test_dist(self):
        for _ in range(10):
            u1, u2, v1, v2 = np.random.randn(4)
            u, v = sjs.Dvec2(u1, u2), sjs.Dvec2(v1, v2)
            d = np.sqrt((v1 - u1)**2 + (v2 - u2)**2)
            self.assertAlmostEqual(sjs.dist(u, v), d)

if __name__ == '__main__':
    unittest.main()
