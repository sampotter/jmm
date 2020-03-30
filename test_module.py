import autograd
import autograd.numpy as np
import sjs
import unittest

from test_util import get_linear_speed_s

class TestLinearSpeed(unittest.TestCase):
    def test(self):
        vx, vy = np.random.uniform(-0.05, 0.05, (2,))

        s_gt = get_linear_speed_s(vx, vy)
        grad_s_gt = autograd.grad(lambda args: s_gt(args[0], args[1]))

        slow = sjs.get_linear_speed_field2(vx, vy)

        for _ in range(10):
            x, y = np.random.uniform(-1, 1, (2,))

        sx_gt, sy_gt = grad_s_gt(np.array([x, y]))
        tmp = slow.grad_s(x, y)
        sx, sy = tmp.x, tmp.y

        self.assertAlmostEqual(s_gt(x, y), slow.s(x, y))
        self.assertAlmostEqual(sx_gt, sx)
        self.assertAlmostEqual(sy_gt, sy)

if __name__ == '__main__':
    unittest.main()
