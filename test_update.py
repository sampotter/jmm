import autograd
import autograd.numpy as np
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

class F4:

    def __init__(self, a_T, a_Tx, a_Ty, p, p0, p1):
        self.a_T = a_T
        self.a_Tx = a_Tx
        self.a_Ty = a_Ty
        self.p = p
        self.p0 = p0
        self.p1 = p1
        self._grad_F4 = autograd.grad(lambda args: self.F4(args))
        self._hess_F4 = autograd.hessian(lambda args: self.F4(args))

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

    def S4(self, args):
        lp = self.lprime(args[0])
        t0 = self.t0(args[0])
        t1 = self.t1(args[1])
        t = 3*lp/2 - (t0 + t1)/4
        return (1 + 2*np.linalg.norm(t))/3

    def F4(self, args):
        T = self.T(args[0])
        L = self.L(args[0])
        S = self.S4(args)
        return T + L*S

    def grad_F4(self, args):
        return self._grad_F4(args)

class TestF4(unittest.TestCase):

    def test_evaluate(self):
        for _ in range(10):
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

            context_gt = F4(a_T, a_Tx, a_Ty, p, p0, p1)

            xy = sjs.Dvec2(*p)
            xy0 = sjs.Dvec2(*p0)
            xy1 = sjs.Dvec2(*p1)

            context = sjs.F4Context(T, Tx, Ty, xy, xy0, xy1)

            for _ in range(10):
                args = np.random.random(2)
                args[1] *= 2*np.pi

                f4_gt = context_gt.F4(args)
                f4 = context.F4(*args)
                self.assertAlmostEqual(f4_gt, f4)

                grad_gt = context_gt.grad_F4(args)
                grad = np.array([context.dF4_deta(*args), context.dF4_dth(*args)])
                self.assertAlmostEqual(grad_gt[0], grad[0])
                self.assertAlmostEqual(grad_gt[1], grad[1])

if __name__ == '__main__':
    unittest.main()
