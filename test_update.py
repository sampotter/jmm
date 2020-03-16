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
        context = sjs.F3Context(cubic, xy, xy0, xy1)
        self.assertAlmostEqual(context.F3(0), 1.0)
        self.assertAlmostEqual(context.F3(1), 1.0)
        self.assertAlmostEqual(context.F3(1/2), 1.0/np.sqrt(2))
        self.assertAlmostEqual(context.dF3_deta(0), -1.0)
        self.assertAlmostEqual(context.dF3_deta(1), 1.0)

    def test_hybrid(self):
        cubic = sjs.Cubic([0, 0, 0, 0])

        xy0 = sjs.Dvec2(1, 0)
        xy1 = sjs.Dvec2(0, 1)
        xy = sjs.Dvec2(0, 0)
        h = 1
        context = sjs.F3Context(cubic, xy, xy0, xy1)

        lam = sjs.hybrid(lambda lam: context.dF3_deta(lam), 0, 1)
        self.assertAlmostEqual(lam, 0.5)
        self.assertAlmostEqual(context.F3(lam), 1.0/np.sqrt(2))
        self.assertAlmostEqual(context.dF3_deta(lam), 0.0)

        context.cubic = sjs.Cubic([0, 1, 1, 1])
        lam = sjs.hybrid(lambda lam: context.dF3_deta(lam), 0, 1)
        self.assertAlmostEqual(lam, 0.0)
        self.assertAlmostEqual(context.F3(lam), 1.0)
        self.assertAlmostEqual(context.dF3_deta(lam), 0.0)

        context.cubic = sjs.Cubic([1, 0, -1, -1])
        lam = sjs.hybrid(lambda lam: context.dF3_deta(lam), 0, 1)
        self.assertAlmostEqual(lam, 1.0)
        self.assertAlmostEqual(context.F3(lam), 1.0)
        self.assertAlmostEqual(context.dF3_deta(lam), 0.0)

        context.cubic = sjs.Cubic([1, 1, -1, 1])
        lam = sjs.hybrid(lambda lam: context.dF3_deta(lam), 0, 1)
        self.assertAlmostEqual(lam, 0.5)
        self.assertAlmostEqual(context.dF3_deta(lam), 0.0)

class TestF4(unittest.TestCase):

    def test_hybrid(self):
        T = sjs.Cubic([
            0.31622776601683794,
            -0.031622776601683805,
            0.014562255152853733,
            0.0008327554319921132
        ])
        Tx = sjs.Cubic([
            0.031622776601683805,
            -0.029124510305707466,
            -0.0024982662959763396,
            0.0
        ])
        Ty = sjs.Cubic([
            0.09486832980505126,
            0.009486832980505164,
            -0.003578655376164433,
            -0.0007765074093920993
        ])
        xy0 = sjs.Dvec2(0.1, 0.3)
        xy1 = sjs.Dvec2(0.0, 0.3)
        xy = sjs.Dvec2(0.1, 0.4)
        context = sjs.F4Context(T, Tx, Ty, xy, xy0, xy1)

        argeta3 = 0.24944078587475857
        argth3 = 1.3263440520581946

        self.assertAlmostEqual(context.F4(argeta3, argth3), 0.4123228462738079)
        self.assertAlmostEqual(context.dF4_deta(argeta3, argth3), 1.145686251245516e-06)
        self.assertAlmostEqual(context.dF4_dth(argeta3, argth3), -2.8057630904607244e-07)

if __name__ == '__main__':
    unittest.main()
