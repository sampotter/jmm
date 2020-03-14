import numpy as np
import sjs
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

    def test_norm(self):
        for _ in range(10):
            u1, u2 = np.random.randn(2)
            u = sjs.Dvec2(u1, u2)
            d = np.sqrt(u1**2 + u2**2)
            self.assertAlmostEqual(u.norm(), d)

    def test_sub(self):
        for _ in range(10):
            u1, u2, v1, v2 = np.random.randn(4)
            u, v = sjs.Dvec2(u1, u2), sjs.Dvec2(v1, v2)
            self.assertAlmostEqual(u*v, u1*v1 + u2*v2)

    def test_dot(self):
        for _ in range(10):
            u1, u2, v1, v2 = np.random.randn(4)
            u, v = sjs.Dvec2(u1, u2), sjs.Dvec2(v1, v2)
            w = u - v
            self.assertAlmostEqual(w.x, u1 - v1)
            self.assertAlmostEqual(w.y, u2 - v2)

    def test_dbl_div(self):
        for _ in range(10):
            u1, u2, a = np.random.randn(3)
            u = sjs.Dvec2(u1, u2)
            v = u/a
            self.assertAlmostEqual(v.x, u1/a)
            self.assertAlmostEqual(v.y, u2/a)

    def test_floor(self):
        for _ in range(10):
            u1, u2 = np.random.randn(2)
            u = sjs.Dvec2(u1, u2)
            v = u.floor()
            self.assertAlmostEqual(v.x, np.floor(u1))
            self.assertAlmostEqual(v.y, np.floor(u2))

class TestDvec4(unittest.TestCase):

    def test_basic(self):
        u = sjs.Dvec4(list(range(4)))
        for i in range(4):
            self.assertEqual(u[i], i)

    def test_dot(self):
        for _ in range(10):
            u_data, v_data = np.random.randn(2, 4)
            u, v = sjs.Dvec4(u_data), sjs.Dvec4(v_data)
            self.assertAlmostEqual(sjs.dot(u, v), u_data@v_data)
            self.assertAlmostEqual(u.dot(v), u_data@v_data)

    def test_sum(self):
        for _ in range(10):
            u_data = np.random.randn(4)
            u = sjs.Dvec4(u_data)
            self.assertAlmostEqual(u.sum(), u_data.sum())
            self.assertAlmostEqual(sjs.sum(u), u_data.sum())

    def test_add(self):
        for _ in range(10):
            u_data, v_data = np.random.randn(2, 4)
            u, v = sjs.Dvec4(u_data), sjs.Dvec4(v_data)
            w = u + v
            w_gt = u_data + v_data
            self.assertAlmostEqual(w[0], w_gt[0])
            self.assertAlmostEqual(w[1], w_gt[1])
            self.assertAlmostEqual(w[2], w_gt[2])
            self.assertAlmostEqual(w[3], w_gt[3])

    def test_div(self):
        for _ in range(10):
            u_data = np.random.randn(4)
            u = sjs.Dvec4(u_data)
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
            m = sjs.Dvec4.m(x)
            for i in range(4):
                self.assertAlmostEqual(m[i], x**i);

    def test_dm(self):
        for _ in range(10):
            x = np.random.randn()
            dm = sjs.Dvec4.dm(x)
            for i in range(4):
                self.assertAlmostEqual(dm[i], i*x**(i - 1))

    def test_e1(self):
        e1 = sjs.Dvec4.e1()
        self.assertEqual(e1[0], 1)
        for i in range(1, 4):
            self.assertEqual(e1[i], 0)

    def test_one(self):
        one = sjs.Dvec4.one()
        for i in range(4):
            self.assertEqual(one[i], 1)

    def test_iota(self):
        iota = sjs.Dvec4.iota()
        for i in range(4):
            self.assertEqual(iota[i], i)

class TestIvec2(unittest.TestCase):

    def test_dvec2_to_ivec2(self):
        for _ in range(10):
            u1, u2 = np.random.randn(2)
            u = sjs.Dvec2(u1, u2)
            v = sjs.Ivec2(u)
            self.assertAlmostEqual(v.i, int(u1))
            self.assertAlmostEqual(v.j, int(u2))

if __name__ == '__main__':
    unittest.main()
