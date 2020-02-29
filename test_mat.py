import itertools as it
import numpy as np
import sjs_eik as sjs
import unittest

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
