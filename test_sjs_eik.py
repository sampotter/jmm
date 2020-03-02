import numpy as np
import sjs_eik as sjs
import unittest

class TestStaticJetScheme(unittest.TestCase):

    def test_set_jet(self):
        shape = (2, 2)
        xymin = (0, 0)
        h = 1
        s = lambda x, y: 1.0 + 0.0*x + 0.0*y
        grad_s = lambda x, y: 0.0 + 0.0*x + 0.0*y
        scheme = sjs.StaticJetScheme(shape, xymin, h, s, grad_s)
        scheme.add_trial(0, 0, sjs.Jet(1, 0, 0, 0))
        scheme.add_trial(1, 0, sjs.Jet(1, 0, 0, 0))
        scheme.add_trial(0, 1, sjs.Jet(1, 0, 0, 0))
        scheme.add_trial(1, 1, sjs.Jet(1, 0, 0, 0))
        self.assertTrue(scheme.can_build_cell(0, 0))
        scheme.build_cells()
        for _ in range(100):
            x, y = np.random.rand(2)
            self.assertAlmostEqual(scheme.T(x, y), 1.0)
            self.assertAlmostEqual(scheme.Tx(x, y), 0.0)
            self.assertAlmostEqual(scheme.Ty(x, y), 0.0)
            self.assertAlmostEqual(scheme.Txy(x, y), 0.0)

if __name__ == '__main__':
    unittest.main()
