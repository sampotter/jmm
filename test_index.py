import itertools as it
import numpy as np
import jmm
import unittest

class TestIndex(unittest.TestCase):

    def test_ind2l(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            ind = (np.random.randint(shape[0]), np.random.randint(shape[1]))
            l = jmm._ind2l(shape, ind)
            l_gt = np.ravel_multi_index(ind, shape, order='C')
            self.assertEqual(l, l_gt)

    def test_ind2lc(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            ind = (
                np.random.randint(shape[0] - 1),
                np.random.randint(shape[1] - 1)
            )
            l = jmm._ind2lc(shape, ind)
            l_gt = np.ravel_multi_index(ind, shape - 1, order='C')
            self.assertEqual(l, l_gt)

    def test_indc2l(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            indc = (
                np.random.randint(shape[0] - 1),
                np.random.randint(shape[1] - 1)
            )
            l = jmm._ind2l(shape, indc)
            l_gt = np.ravel_multi_index(indc, shape, order='C')
            self.assertEqual(l, l_gt)

    def test_indc2lc(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            indc = (
                np.random.randint(shape[0] - 1),
                np.random.randint(shape[1] - 1)
            )
            lc = jmm._indc2lc(shape, indc)
            lc_gt = np.ravel_multi_index(indc, shape - 1, order='C')
            self.assertEqual(lc, lc_gt)

    def test_l2ind(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            l = np.random.randint(np.product(shape))
            ind = jmm._l2ind(shape, l)
            ind_gt = np.unravel_index(l, shape)
            self.assertEqual(ind.i, ind_gt[0])
            self.assertEqual(ind.j, ind_gt[1])

    def test_l2indc(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            l = np.random.randint(np.product(shape))
            indc = jmm._l2indc(shape, l)
            indc_gt = np.unravel_index(l, shape)
            self.assertEqual(indc.i, indc_gt[0])
            self.assertEqual(indc.j, indc_gt[1])

    def test_lc2ind(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            lc = np.random.randint((shape[0] - 1)*(shape[1] - 1))
            ind = jmm._lc2ind(shape, lc)
            ind_gt = np.unravel_index(lc, (shape[0] - 1, shape[1] - 1))
            self.assertEqual(ind.i, ind_gt[0])
            self.assertEqual(ind.j, ind_gt[1])

    def test_lc2indc(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            lc = np.random.randint((shape[0] - 1)*(shape[1] - 1))
            indc = jmm._lc2indc(shape, lc)
            indc_gt = np.unravel_index(lc, (shape[0] - 1, shape[1] - 1))
            self.assertEqual(indc.i, indc_gt[0])
            self.assertEqual(indc.j, indc_gt[1])

    def test_l2lc(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            l = np.random.randint(np.product(shape))
            lc = jmm._l2lc(shape, l)
            ind_gt = jmm._l2ind(shape, l)
            lc_gt = jmm._ind2lc(shape, (ind_gt.i, ind_gt.j))
            self.assertEqual(lc, lc_gt)

    def test_lc2l(self):
        for _ in range(10):
            shape = np.random.randint((2, 1), 1000)
            # TODO: this fails sometimes! I wonder if it has to do
            # with the bleeding at the edge of the domain?
            lc = np.random.randint((shape[0] - 1)*(shape[1] - 1))
            l = jmm._lc2l(shape, lc)
            indc_gt = jmm._lc2indc(shape, lc)
            l_gt = jmm._indc2l(shape, (indc_gt.i, indc_gt.j))
            self.assertEqual(l, l_gt)

    def test_xy_to_lc_and_cc(self):
        shape = (2, 2)
        xymin = (0, 0)
        h = 1

        for x, y in it.product(range(2), repeat=2):
            lc, cc = jmm._xy_to_lc_and_cc(shape, xymin, h, (x, y))
            self.assertEqual(lc, 0)
            self.assertEqual(x, cc.x)
            self.assertEqual(y, cc.y)

        for _ in range(10):
            x, y = np.random.rand(2)
            lc, cc = jmm._xy_to_lc_and_cc(shape, xymin, h, (x, y))
            self.assertEqual(lc, 0)
            self.assertEqual(x, cc.x)
            self.assertEqual(y, cc.y)

        for _ in range(10):
            shape = np.random.randint((2, 1), 99) + 1
            xymin = np.random.randn(2)
            h = np.random.rand()

            xy = h*np.random.randint(shape) + xymin
            assert(all(xy[i] >= xymin[i] for i in range(2)))
            assert(all(xy[i] <= xymin[i] + h*shape[i] for i in range(2)))

            cc_gt = (xy - xymin)/h
            ind_gt = np.floor(cc_gt).astype(int)
            cc_gt -= ind_gt
            for i in range(2):
                if ind_gt[i] < 0:
                    ind_gt[i] = 0
                    cc_gt[i] = 0.0
                if ind_gt[i] >= shape[i] - 1:
                    ind_gt[i] -= 1
                    cc_gt[i] = 1.0
            lc_gt = np.ravel_multi_index(ind_gt, shape - 1, order='C')

            lc, cc = jmm._xy_to_lc_and_cc(shape, xymin, h, xy)

            self.assertEqual(lc, lc_gt)
            self.assertEqual(cc.x, cc_gt[0])
            self.assertEqual(cc.y, cc_gt[1])

if __name__ == '__main__':
    unittest.main()
