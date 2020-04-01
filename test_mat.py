import itertools as it
import numpy as np
import sjs
import unittest

class TestDmat22(unittest.TestCase):

    def test_add(self):
        for _ in range(10):
            A_data, B_data = np.random.randn(2, 2, 2)
            C_gt = A_data + B_data
            A = sjs.Dmat22(A_data)
            B = sjs.Dmat22(B_data)
            C = A + B
            for i, j in it.product(range(2), repeat=2):
                self.assertAlmostEqual(C[i, j], C_gt[i, j])

    def test_sub(self):
        for _ in range(10):
            A_data, B_data = np.random.randn(2, 2, 2)
            C_gt = A_data - B_data
            A = sjs.Dmat22(A_data)
            B = sjs.Dmat22(B_data)
            C = A - B
            for i, j in it.product(range(2), repeat=2):
                self.assertAlmostEqual(C[i, j], C_gt[i, j])

    def test_mul(self):
        for _ in range(10):
            A_data, B_data = np.random.randn(2, 2, 2)
            C_gt = A_data@B_data
            A = sjs.Dmat22(A_data)
            B = sjs.Dmat22(B_data)
            C = A@B
            for i, j in it.product(range(2), repeat=2):
                self.assertAlmostEqual(C[i, j], C_gt[i, j])
            x_data = np.random.randn(2)
            x = sjs.Dvec2(*x_data)
            y, y_gt = A@x, A_data@x_data
            self.assertAlmostEqual(y[0], y_gt[0])
            self.assertAlmostEqual(y[1], y_gt[1])

    def test_dbl_mul(self):
        for _ in range(10):
            A_data = np.random.randn(2, 2)
            a = np.random.randn()
            C_gt = a*A_data
            A = sjs.Dmat22(A_data)
            C = A*a
            for i, j in it.product(range(2), repeat=2):
                self.assertAlmostEqual(C[i, j], C_gt[i, j])

    def test_dbl_mul(self):
        for _ in range(10):
            A_data = np.random.randn(2, 2)
            a = np.random.randn()
            C_gt = A_data/a
            A = sjs.Dmat22(A_data)
            C = A/a
            for i, j in it.product(range(2), repeat=2):
                self.assertAlmostEqual(C[i, j], C_gt[i, j])

    def test_solve(self):
        for _ in range(10):
            A_gt = np.random.randn(2, 2)
            b_gt = np.random.randn(2)
            A, b = sjs.Dmat22(A_gt), sjs.Dvec2(*b_gt)
            x, x_gt = A.solve(b), np.linalg.solve(A_gt, b_gt)
            self.assertAlmostEqual(x[0], x_gt[0])
            self.assertAlmostEqual(x[1], x_gt[1])

    def test_outer(self):
        for _ in range(10):
            u_gt, v_gt = np.random.randn(2, 2)
            u, v = sjs.Dvec2(*u_gt), sjs.Dvec2(*v_gt)
            uv, uv_gt = sjs.outer(u, v), np.outer(u_gt, v_gt)
            for i, j in it.product(range(2), repeat=2):
                self.assertAlmostEqual(uv[i, j], uv_gt[i, j])

    def test_invert(self):
        for _ in range(10):
            A_gt = np.random.randn(2, 2)
            A = sjs.Dmat22(A_gt)
            A.invert()
            A_gt = np.linalg.inv(A_gt)
            for i, j in it.product(range(2), repeat=2):
                self.assertAlmostEqual(A[i, j], A_gt[i, j])

    def test_trace(self):
        for _ in range(10):
            A_gt = np.random.randn(2, 2)
            A = sjs.Dmat22(A_gt)
            trace, trace_gt = A.trace(), np.trace(A_gt)
            self.assertAlmostEqual(trace, trace_gt)

    def test_det(self):
        for _ in range(10):
            A_gt = np.random.randn(2, 2)
            A = sjs.Dmat22(A_gt)
            det, det_gt = A.det(), np.linalg.det(A_gt)
            self.assertAlmostEqual(det, det_gt)

    def test_eigvals(self):
        for _ in range(10):
            A_gt = np.random.randn(2, 2)
            A_gt = (A_gt + A_gt.T)/2
            A = sjs.Dmat22(A_gt)
            lam1, lam2 = A.eigvals()
            lam1_gt, lam2_gt = np.linalg.eig(A_gt)[0]
            lam1_gt, lam2_gt = max(lam1_gt, lam2_gt), min(lam1_gt, lam2_gt)
            self.assertAlmostEqual(lam1, lam1_gt)
            self.assertAlmostEqual(lam2, lam2_gt)

    def test_transpose(self):
        for _ in range(10):
            A_gt = np.random.randn(2, 2)
            A = sjs.Dmat22(A_gt)
            A.transpose()
            A_gt = A_gt.T
            for i, j in it.product(range(2), repeat=2):
                self.assertAlmostEqual(A[i, j], A_gt[i, j])

class TestDmat44(unittest.TestCase):

    def test_ctor(self):
        for _ in range(10):
            A_data = np.random.randn(4, 4)
            A = sjs.Dmat44(A_data)
            for i, j in it.product(range(4), repeat=2):
                self.assertEqual(A[i, j], A_data[i, j])

    def test_vec_mat_mul(self):
        for _ in range(10):
            x_data = np.random.randn(4)
            A_data = np.random.randn(4, 4)
            x, A = sjs.Dvec4(x_data), sjs.Dmat44(A_data)
            y = x*A
            y_gt = x_data@A_data
            for i in range(4):
                self.assertAlmostEqual(y[i], y_gt[i])

    def test_mat_vec_mul(self):
        for _ in range(10):
            A_data = np.random.randn(4, 4)
            x_data = np.random.randn(4)
            A, x = sjs.Dmat44(A_data), sjs.Dvec4(x_data)
            y = A*x
            y_gt = A_data@x_data
            for i in range(4):
                self.assertAlmostEqual(y[i], y_gt[i])

    def test_mat_mat_mul(self):
        for _ in range(10):
            A_data = np.random.randn(4, 4)
            B_data = np.random.randn(4, 4)
            A, B = sjs.Dmat44(A_data), sjs.Dmat44(B_data)
            C = A*B
            C_gt = A_data@B_data
            for i, j in it.product(range(4), repeat=2):
                self.assertAlmostEqual(C[i, j], C_gt[i, j])

    def test_col(self):
        for _ in range(10):
            A_data = np.random.randn(4, 4)
            A = sjs.Dmat44(A_data)
            for j in range(4):
                col = sjs.col(A, j)
                for i in range(4):
                    self.assertEqual(col[i], A_data[i, j])

if __name__ == '__main__':
    unittest.main()
