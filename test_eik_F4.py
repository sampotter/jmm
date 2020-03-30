import autograd
import autograd.numpy as np
import sjs
import unittest

from scipy.optimize import brentq
from test_util import get_linear_speed_s, get_linear_speed_tau

class F4:
    def __init__(self, s, a_T, a_Tx, a_Ty, p, p0, p1):
        self._s = s
        self._grad_s = autograd.grad(lambda args: self._s(args))
        self.a_T = a_T
        self.a_Tx = a_Tx
        self.a_Ty = a_Ty
        self.p = p
        self.p0 = p0
        self.p1 = p1
        self._grad_F4 = autograd.grad(lambda args: self.F4(args))
        self._hess_F4 = autograd.hessian(lambda args: self.F4(args))

    def s(self, x):
        return self._s(*x)

    def grad_s(self, x):
        return self._grad_s(x)

    def T(self, lam):
        return self.a_T@[lam**p for p in range(4)]

    def Tx(self, lam):
        return self.a_Tx@[lam**p for p in range(4)]

    def Ty(self, lam):
        return self.a_Ty@[lam**p for p in range(4)]

    def plam(self, lam):
        return (1 - lam)*self.p0 + lam*self.p1

    def L(self, lam):
        return np.linalg.norm(self.plam(lam) - self.p)

    def lprime(self, lam):
        return (self.p - self.plam(lam))/self.L(lam)

    def gradTlam(self, lam):
        return np.array([self.Tx(lam), self.Ty(lam)])

    def t0(self, lam):
        return self.gradTlam(lam)/np.linalg.norm(self.gradTlam(lam))

    def t1(self, th):
        return np.array([np.cos(th), np.sin(th)])

    def s0(self, lam):
        return self.s(self.plam(lam))

    def pm(self, args):
        lam, th = args
        return \
            (self.p + self.plam(lam))/2 \
            - self.L(lam)*(self.t1(th) - self.t0(lam))/8

    def tm(self, args):
        lp = self.lprime(args[0])
        t0 = self.t0(args[0])
        t1 = self.t1(args[1])
        return 3*lp/2 - (t0 + t1)/4

    def sm(self, args):
        return self.s(self.pm(args))

    def S4(self, args):
        tm = self.tm(args)
        tmnorm = np.linalg.norm(tm)
        s0 = self.s0(args[0])
        s1 = self.s(self.p)
        sm = self.sm(args)
        return (s0 + s1 + 4*sm*tmnorm)/6

    def F4(self, args):
        T = self.T(args[0])
        L = self.L(args[0])
        S = self.S4(args)
        return T + L*S

    def grad_F4(self, args):
        return self._grad_F4(args)

    def hess_F4(self, args):
        return self._hess_F4(args)

class TestF4(unittest.TestCase):

    def test_evaluate(self):
        eps = 1e-7

        for _ in range(10):
            vx, vy = np.random.uniform(-0.05, 0.05, (2,))
            s_gt = get_linear_speed_s(vx, vy)
            slow = sjs.get_linear_speed_field2(vx, vy)

            data = np.random.randn(4, 4)
            h = np.random.random()
            H = np.diag([1, 1, h, h])
            data = H@data@H

            p = np.random.randn(2)
            p0 = p + h*np.random.randn(2)
            p1 = p + h*np.random.randn(2)

            bicubic = sjs.Bicubic(data)
            T = bicubic.get_f_on_edge(sjs.BicubicVariable.Lambda, 0)
            Tx = bicubic.get_fx_on_edge(sjs.BicubicVariable.Lambda, 0)
            Ty = bicubic.get_fy_on_edge(sjs.BicubicVariable.Lambda, 0)

            a_T = np.array([T.a[i] for i in range(4)])
            a_Tx = np.array([Tx.a[i] for i in range(4)])
            a_Ty = np.array([Ty.a[i] for i in range(4)])

            context_gt = F4(s_gt, a_T, a_Tx, a_Ty, p, p0, p1)

            xy = sjs.Dvec2(*p)
            xy0 = sjs.Dvec2(*p0)
            xy1 = sjs.Dvec2(*p1)

            context = sjs.F4Context(T, Tx, Ty, xy, xy0, xy1, slow)

            for _ in range(10):
                args = np.random.random(2)
                args[1] *= 2*np.pi

                context.compute(*args)

                f4_gt = context_gt.F4(args)
                f4_eta_gt, f4_th_gt = context_gt.grad_F4(args)

                self.assertAlmostEqual(f4_gt, context.F4)
                self.assertAlmostEqual(f4_eta_gt, context.F4_eta)
                self.assertAlmostEqual(f4_th_gt, context.F4_th)

    def test_bfgs_linear_speed(self):
        for _ in range(10):
            h = 0.1

            x0, y0 = np.random.uniform(-1, 1 - h, (2,))
            x1, y1 = x0 + h, y0 + h
            x, y = x0 - h, y0

            vx, vy = np.random.uniform(-0.05, 0.05, (2,))
            s_gt = get_linear_speed_s(vx, vy)
            slow = sjs.get_linear_speed_field2(vx, vy)

            tau_Omega = get_linear_speed_tau(vx, vy)
            tau = lambda lam, mu: tau_Omega(
                (1 - lam)*x0 + lam*x1,
                (1 -  mu)*y0 +  mu*y1
            )
            grad_tau = autograd.grad(lambda args: tau(args[0], args[1]))
            hess_tau = autograd.hessian(lambda args: tau(args[0], args[1]))
            tau_x = lambda lam, mu: grad_tau(np.array([lam, mu]))[0]
            tau_y = lambda lam, mu: grad_tau(np.array([lam, mu]))[1]
            tau_xy = lambda lam, mu: hess_tau(np.array([lam, mu]))[1][0]

            data = np.array([
                [  tau(0., 0.),    tau(0., 1.),  tau_y(0., 0.),  tau_y(0., 1.)],
                [  tau(1., 0.),    tau(1., 1.),  tau_y(1., 0.),  tau_y(1., 1.)],
                [tau_x(0., 0.),  tau_x(0., 1.), tau_xy(0., 0.), tau_xy(0., 1.)],
                [tau_x(1., 0.),  tau_x(1., 1.), tau_xy(1., 0.), tau_xy(1., 1.)],
            ])
            bicubic = sjs.Bicubic(data)
            T = bicubic.get_f_on_edge(sjs.BicubicVariable.Lambda, 0)
            Tx = bicubic.get_fx_on_edge(sjs.BicubicVariable.Lambda, 0)
            Ty = bicubic.get_fy_on_edge(sjs.BicubicVariable.Lambda, 0)

            a_T = np.array([T.a[i] for i in range(4)])
            a_Tx = np.array([Tx.a[i] for i in range(4)])
            a_Ty = np.array([Ty.a[i] for i in range(4)])

            p, p0, p1 = np.array([x, y]), np.array([x0, y0]), np.array([x0, y1])
            context_gt = F4(s_gt, a_T, a_Tx, a_Ty, p, p0, p1)

            xy = sjs.Dvec2(*p)
            xy0 = sjs.Dvec2(*p0)
            xy1 = sjs.Dvec2(*p1)

            context = sjs.F4Context(T, Tx, Ty, xy, xy0, xy1, slow)

            context3 = sjs.F3Context(T, xy, xy0, xy1, slow)
            def F3(eta):
                context3.compute(eta)
                return context3.F3
            def F3_eta(eta):
                context3.compute(eta)
                return context3.F3_eta
            if np.sign(F3_eta(0.)) == np.sign(F3_eta(1.)):
                argeta3 = 0. if F3(0.) < F3(1.) else 1.
            else:
                argeta3 = brentq(F3_eta, 0, 1)
            lp = (p - p0 - argeta3*(p1 - p0))
            lp /= np.linalg.norm(lp)
            argth3 = np.arctan2(*reversed(lp))

            xk_gt = np.array([argeta3, argth3])

            eps = 1e-7

            hess_gt = context_gt.hess_F4(xk_gt)
            hess_fd = context.hess_fd(argeta3, argth3, eps)

            self.assertTrue(abs(hess_gt[0, 0] - hess_fd[0, 0]) < eps)
            self.assertTrue(abs(hess_gt[1, 0] - hess_fd[1, 0]) < eps)
            self.assertTrue(abs(hess_gt[0, 1] - hess_fd[0, 1]) < eps)
            self.assertTrue(abs(hess_gt[1, 1] - hess_fd[1, 1]) < eps)

            gk_gt = context_gt.grad_F4(xk_gt)
            Hk_gt = np.linalg.inv(hess_gt)

            xk, gk, Hk = context.bfgs_init(*xk_gt)

            eps = 1e-6
            self.assertTrue(abs(Hk_gt[0, 0] - Hk[0, 0]) < eps)
            self.assertTrue(abs(Hk_gt[1, 0] - Hk[1, 0]) < eps)
            self.assertTrue(abs(Hk_gt[0, 1] - Hk[0, 1]) < eps)
            self.assertTrue(abs(Hk_gt[1, 1] - Hk[1, 1]) < eps)

if __name__ == '__main__':
    unittest.main()
