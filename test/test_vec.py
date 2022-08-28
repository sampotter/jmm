import numpy as np
import jmm
import unittest

class TestDvec2(unittest.TestCase):

    def test_ctor(self):
        v = jmm.Dvec2(1, 2)
        self.assertEqual(v.x, 1)
        self.assertEqual(v.y, 2)

    def test_ccomb(self):
        u, v = jmm.Dvec2(2, 0), jmm.Dvec2(1, 3)
        w = jmm.ccomb(u, v, 0.1)
        self.assertAlmostEqual(w.x, 1.9)
        self.assertAlmostEqual(w.y, 0.3)

    def test_dist(self):
        for _ in range(10):
            u1, u2, v1, v2 = np.random.randn(4)
            u, v = jmm.Dvec2(u1, u2), jmm.Dvec2(v1, v2)
            d = np.sqrt((v1 - u1)**2 + (v2 - u2)**2)
            self.assertAlmostEqual(jmm.dist(u, v), d)

    def test_norm(self):
        for _ in range(10):
            u1, u2 = np.random.randn(2)
            u = jmm.Dvec2(u1, u2)
            d = np.sqrt(u1**2 + u2**2)
            self.assertAlmostEqual(u.norm(), d)

    def test_sub(self):
        for _ in range(10):
            u1, u2, v1, v2 = np.random.randn(4)
            u, v = jmm.Dvec2(u1, u2), jmm.Dvec2(v1, v2)
            self.assertAlmostEqual(u*v, u1*v1 + u2*v2)

    def test_dot(self):
        for _ in range(10):
            u1, u2, v1, v2 = np.random.randn(4)
            u, v = jmm.Dvec2(u1, u2), jmm.Dvec2(v1, v2)
            self.assertAlmostEqual(u*v, u1*v1 + u2*v2)

    def test_add(self):
        for _ in range(10):
            u1, u2, v1, v2 = np.random.randn(4)
            u, v = jmm.Dvec2(u1, u2), jmm.Dvec2(v1, v2)
            w = u + v
            self.assertAlmostEqual(w.x, u1 + v1)
            self.assertAlmostEqual(w.y, u2 + v2)

    def test_sub(self):
        for _ in range(10):
            u1, u2, v1, v2 = np.random.randn(4)
            u, v = jmm.Dvec2(u1, u2), jmm.Dvec2(v1, v2)
            w = u - v
            self.assertAlmostEqual(w.x, u1 - v1)
            self.assertAlmostEqual(w.y, u2 - v2)

    def test_dbl_div(self):
        for _ in range(10):
            u1, u2, a = np.random.randn(3)
            u = jmm.Dvec2(u1, u2)
            v = u/a
            self.assertAlmostEqual(v.x, u1/a)
            self.assertAlmostEqual(v.y, u2/a)

    def test_dbl_mul(self):
        for _ in range(10):
            u1, u2, a = np.random.randn(3)
            u = jmm.Dvec2(u1, u2)
            v = u*a
            self.assertAlmostEqual(v.x, a*u1)
            self.assertAlmostEqual(v.y, a*u2)
            u *= a
            self.assertAlmostEqual(v.x, u.x)
            self.assertAlmostEqual(v.y, u.y)

    def test_floor(self):
        for _ in range(10):
            u1, u2 = np.random.randn(2)
            u = jmm.Dvec2(u1, u2)
            v = u.floor()
            self.assertAlmostEqual(v.x, np.floor(u1))
            self.assertAlmostEqual(v.y, np.floor(u2))

    def test_normalize(self):
        for _ in range(10):
            u = jmm.Dvec2(*np.random.randn(2))
            u.normalize()
            self.assertAlmostEqual(u.norm(), 1.0)

    def test_cproj(self):
        for _ in range(10):
            u, v = np.random.randn(2, 2)
            w = jmm.cproj(jmm.Dvec2(*u), jmm.Dvec2(*v))
            w_gt = (np.eye(2) - np.outer(u, u))@v
            self.assertAlmostEqual(w[0], w_gt[0])
            self.assertAlmostEqual(w[1], w_gt[1])

class TestDvec4(unittest.TestCase):

    def test_basic(self):
        u = jmm.Dvec4(list(range(4)))
        for i in range(4):
            self.assertEqual(u[i], i)

    def test_dot(self):
        for _ in range(10):
            u_data, v_data = np.random.randn(2, 4)
            u, v = jmm.Dvec4(u_data), jmm.Dvec4(v_data)
            self.assertAlmostEqual(jmm.dot(u, v), u_data@v_data)
            self.assertAlmostEqual(u.dot(v), u_data@v_data)

    def test_sum(self):
        for _ in range(10):
            u_data = np.random.randn(4)
            u = jmm.Dvec4(u_data)
            self.assertAlmostEqual(u.sum(), u_data.sum())
            self.assertAlmostEqual(jmm.sum(u), u_data.sum())

    def test_add(self):
        for _ in range(10):
            u_data, v_data = np.random.randn(2, 4)
            u, v = jmm.Dvec4(u_data), jmm.Dvec4(v_data)
            w = u + v
            w_gt = u_data + v_data
            self.assertAlmostEqual(w[0], w_gt[0])
            self.assertAlmostEqual(w[1], w_gt[1])
            self.assertAlmostEqual(w[2], w_gt[2])
            self.assertAlmostEqual(w[3], w_gt[3])

    def test_div(self):
        for _ in range(10):
            u_data = np.random.randn(4)
            u = jmm.Dvec4(u_data)
            a = np.random.randn()
            v = u/a
            v_gt = u_data/a
            self.assertAlmostEqual(v[0], v_gt[0])
            self.assertAlmostEqual(v[1], v_gt[1])
            self.assertAlmostEqual(v[2], v_gt[2])
            self.assertAlmostEqual(v[3], v_gt[3])

    def test_m(self):
        for _ in range(10):
            x = np.random.randn()
            m = jmm.Dvec4.m(x)
            for i in range(4):
                self.assertAlmostEqual(m[i], x**i);

    def test_dm(self):
        for _ in range(10):
            x = np.random.randn()
            dm = jmm.Dvec4.dm(x)
            for i in range(4):
                self.assertAlmostEqual(dm[i], i*x**(i - 1))

    def test_e1(self):
        e1 = jmm.Dvec4.e1()
        self.assertEqual(e1[0], 1)
        for i in range(1, 4):
            self.assertEqual(e1[i], 0)

    def test_one(self):
        one = jmm.Dvec4.one()
        for i in range(4):
            self.assertEqual(one[i], 1)

    def test_iota(self):
        iota = jmm.Dvec4.iota()
        for i in range(4):
            self.assertEqual(iota[i], i)

class TestIvec2(unittest.TestCase):

    def test_dvec2_to_ivec2(self):
        for _ in range(10):
            u1, u2 = np.random.randn(2)
            u = jmm.Dvec2(u1, u2)
            v = jmm.Ivec2(u)
            self.assertAlmostEqual(v.i, int(u1))
            self.assertAlmostEqual(v.j, int(u2))

if __name__ == '__main__':
    unittest.main()
