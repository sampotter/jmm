import autograd
import autograd.numpy as np
import sjs
import unittest

from scipy.optimize import brentq

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
        return self._s(x)

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
            v0 = (1 + np.random.random())/2
            vx, vy = np.random.randn(2)/4
            while max(abs(vx), abs(vy)) > 1/3:
                vx, vy = np.random.randn(2)/4
            s = lambda x, y: 0.399 + vx*x + vy*y
            slow = sjs.Field2(s, lambda x, y: (vx, vy))

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

            context_gt = F4(lambda args: s(*args), a_T, a_Tx, a_Ty, p, p0, p1)

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

                hess_gt = context_gt.hess_F4(args)
                hess_fd = context.hess_fd(*args, eps)

    def test_bfgs(self):
        eps = 1e-7

        for _ in range(10):
            v0 = (1 + np.random.random())/2
            vx, vy = np.random.randn(2)/4
            while max(abs(vx), abs(vy)) > 1/3:
                vx, vy = np.random.randn(2)/4
            s = lambda x, y: 0.399 + vx*x + vy*y
            slow = sjs.Field2(s, lambda x, y: (vx, vy))

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

            context_gt = F4(lambda args: s(*args), a_T, a_Tx, a_Ty, p, p0, p1)

            xy = sjs.Dvec2(*p)
            xy0 = sjs.Dvec2(*p0)
            xy1 = sjs.Dvec2(*p1)

            context = sjs.F4Context(T, Tx, Ty, xy, xy0, xy1, slow)



            self.assertTrue(abs(hess_gt[0, 0] - hess_fd[0, 0]) < eps)
            self.assertTrue(abs(hess_gt[1, 0] - hess_fd[1, 0]) < eps)
            self.assertTrue(abs(hess_gt[0, 1] - hess_fd[0, 1]) < eps)
            self.assertTrue(abs(hess_gt[1, 1] - hess_fd[1, 1]) < eps)

            xk_gt = args
            gk_gt = context_gt.grad_F4(xk_gt)
            Hk_gt = np.linalg.inv(context_gt.hess_F4(xk_gt))

            xk, gk, Hk = context.bfgs_init(*args)

            eps = 1e-5

            print(Hk_gt[1, 1] - Hk[1, 1])

            self.assertTrue(abs(Hk_gt[0, 0] - Hk[0, 0]) < eps)
            self.assertTrue(abs(Hk_gt[1, 0] - Hk[1, 0]) < eps)
            self.assertTrue(abs(Hk_gt[0, 1] - Hk[0, 1]) < eps)
            self.assertTrue(abs(Hk_gt[1, 1] - Hk[1, 1]) < eps)

if __name__ == '__main__':
    unittest.main()
