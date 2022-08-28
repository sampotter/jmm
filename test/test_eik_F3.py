import autograd
import autograd.numpy as np
import jmm
import unittest

from scipy.optimize import brentq
from test_util import get_linear_speed_s, get_linear_speed_tau

class F3:
    def __init__(self, s, a_T, p, p0, p1):
        self._s = s
        self.s1 = s(*p)
        self.a_T = a_T
        self.p = p
        self.p0 = p0
        self.dp = p1 - p0
        self._T_eta = autograd.grad(lambda eta: self.T(eta))
        self._L_eta = autograd.grad(lambda eta: self.L(eta))
        self._grad_F3 = autograd.grad(lambda eta: self.F3(eta))

    def T(self, eta):
        return self.a_T@[eta**p for p in range(4)]

    def T_eta(self, eta):
        return self._T_eta(eta)

    def L(self, eta):
        return np.linalg.norm(self.p - self.p0 - eta*self.dp)

    def L_eta(self, eta):
        return self._L_eta(eta)

    def s0(self, eta):
        peta = self.p0 + eta*self.dp
        return self._s(*peta)

    def S1(self, eta):
        return (self.s1 + self.s0(eta))/2

    def F3(self, eta):
        return self.T(eta) + self.L(eta)*self.S1(eta)

    def grad_F3(self, eta):
        return self._grad_F3(eta)

class TestF3(unittest.TestCase):

    def test_evaluate_constant_slowness(self):
        cubic = jmm.Cubic([0, 0, 0, 0])
        xy0 = jmm.Dvec2(1, 0)
        xy1 = jmm.Dvec2(0, 1)
        xy = jmm.Dvec2(0, 0)
        h = 1
        slow = jmm.get_constant_slowness_field2()
        context = jmm.F3Context(cubic, xy, xy0, xy1, slow)

        context.compute(0)
        self.assertAlmostEqual(context.F3, 1.0)
        self.assertAlmostEqual(context.F3_eta, -1.0)

        context.compute(1)
        self.assertAlmostEqual(context.F3, 1.0)
        self.assertAlmostEqual(context.F3_eta, 1.0)

        context.compute(1/2)
        self.assertAlmostEqual(context.F3, 1.0/np.sqrt(2))

    def test_evaluate_linear_speed(self):
        for _ in range(10):
            h = 0.1

            x0, y0 = np.random.uniform(-1, 1 - h, (2,))
            x1, y1 = x0 + h, y0 + h
            x, y = x0 - h, y0

            xy = jmm.Dvec2(x, y)
            xy0 = jmm.Dvec2(x0, y0)
            xy1 = jmm.Dvec2(x0, y1)

            vx, vy = np.random.uniform(-0.05, 0.05, (2,))

            tau_xy = get_linear_speed_tau(vx, vy)
            tau = lambda eta: tau_xy(x0, y0 + eta*(y1 - y0))
            tau_eta = autograd.grad(tau)

            s_gt = get_linear_speed_s(vx, vy)
            slow = jmm.get_linear_speed_field2(vx, vy)

            cubic = jmm.Cubic([tau(0.0), tau(1.0), tau_eta(0.0), tau_eta(1.0)])

            context = jmm.F3Context(cubic, xy, xy0, xy1, slow)

            a_T = np.array([cubic.a[0], cubic.a[1], cubic.a[2], cubic.a[3]])
            p, p0, p1 = np.array([x, y]), np.array([x0, y0]), np.array([x0, y1])
            context_gt = F3(s_gt, a_T, p, p0, p1)

            for _ in range(10):
                eta = np.random.random()

                T, T_gt = cubic.f(eta), context_gt.T(eta)
                T_eta, T_eta_gt = cubic.df(eta), context_gt.T_eta(eta)

                self.assertAlmostEqual(T, T_gt)
                self.assertAlmostEqual(T_eta, T_eta_gt)

                context.compute(eta)

                f3, f3_eta = context.F3, context.F3_eta
                f3_gt, f3_eta_gt = context_gt.F3(eta), context_gt.grad_F3(eta)

                self.assertAlmostEqual(f3, f3_gt)
                self.assertAlmostEqual(f3_eta, f3_eta_gt)

    def test_hybrid_constant_slowness(self):
        cubic = jmm.Cubic([0, 0, 0, 0])

        xy0 = jmm.Dvec2(1, 0)
        xy1 = jmm.Dvec2(0, 1)
        xy = jmm.Dvec2(0, 0)
        h = 1
        slow = jmm.get_constant_slowness_field2()
        context = jmm.F3Context(cubic, xy, xy0, xy1, slow)

        def func(lam):
            context.compute(lam)
            return context.F3_eta

        lam = jmm.hybrid(func, 0, 1)
        self.assertAlmostEqual(lam, 0.5)
        self.assertAlmostEqual(context.F3, 1.0/np.sqrt(2))
        self.assertAlmostEqual(context.F3_eta, 0.0)

        context.T_cubic = jmm.Cubic([0, 1, 1, 1])
        lam = jmm.hybrid(func, 0, 1)
        self.assertAlmostEqual(lam, 0.0)
        self.assertAlmostEqual(context.F3, 1.0)
        self.assertAlmostEqual(context.F3_eta, 0.0)

        context.T_cubic = jmm.Cubic([1, 0, -1, -1])
        lam = jmm.hybrid(func, 0, 1)
        self.assertAlmostEqual(lam, 1.0)
        self.assertAlmostEqual(context.F3, 1.0)
        self.assertAlmostEqual(context.F3_eta, 0.0)

        context.T_cubic = jmm.Cubic([1, 1, -1, 1])
        lam = jmm.hybrid(func, 0, 1)
        self.assertAlmostEqual(lam, 0.5)
        self.assertAlmostEqual(context.F3_eta, 0.0)

if __name__ == '__main__':
    unittest.main()
