import numpy as np
import sjs
import unittest

class TestF3(unittest.TestCase):

    def test_evaluate(self):
        cubic = sjs.Cubic([0, 0, 0, 0])
        xy0 = sjs.Dvec2(1, 0)
        xy1 = sjs.Dvec2(0, 1)
        xy = sjs.Dvec2(0, 0)
        h = 1
        context = sjs.F3Context(cubic, xy, xy0, xy1, sjs.slow1)

        context.compute(0)
        self.assertAlmostEqual(context.F3, 1.0)
        self.assertAlmostEqual(context.F3_eta, -1.0)

        context.compute(1)
        self.assertAlmostEqual(context.F3, 1.0)
        self.assertAlmostEqual(context.F3_eta, 1.0)

        context.compute(1/2)
        self.assertAlmostEqual(context.F3, 1.0/np.sqrt(2))

    def test_hybrid(self):
        cubic = sjs.Cubic([0, 0, 0, 0])

        xy0 = sjs.Dvec2(1, 0)
        xy1 = sjs.Dvec2(0, 1)
        xy = sjs.Dvec2(0, 0)
        h = 1
        context = sjs.F3Context(cubic, xy, xy0, xy1, sjs.slow1)

        def func(lam):
            context.compute(lam)
            return context.F3_eta

        lam = sjs.hybrid(func, 0, 1)
        self.assertAlmostEqual(lam, 0.5)
        self.assertAlmostEqual(context.F3, 1.0/np.sqrt(2))
        self.assertAlmostEqual(context.F3_eta, 0.0)

        context.T_cubic = sjs.Cubic([0, 1, 1, 1])
        lam = sjs.hybrid(func, 0, 1)
        self.assertAlmostEqual(lam, 0.0)
        self.assertAlmostEqual(context.F3, 1.0)
        self.assertAlmostEqual(context.F3_eta, 0.0)

        context.T_cubic = sjs.Cubic([1, 0, -1, -1])
        lam = sjs.hybrid(func, 0, 1)
        self.assertAlmostEqual(lam, 1.0)
        self.assertAlmostEqual(context.F3, 1.0)
        self.assertAlmostEqual(context.F3_eta, 0.0)

        context.T_cubic = sjs.Cubic([1, 1, -1, 1])
        lam = sjs.hybrid(func, 0, 1)
        self.assertAlmostEqual(lam, 0.5)
        self.assertAlmostEqual(context.F3_eta, 0.0)
