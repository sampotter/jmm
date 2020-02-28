import numpy as np
import sjs_eik
import unittest

class TestCubic(unittest.TestCase):

    def test_basic(self):
        for _ in range(10):
            data = np.random.randn(4)
            cubic = sjs_eik.Cubic(data)
            self.assertTrue(np.isclose(cubic.f(0), data[0]))
            self.assertTrue(np.isclose(cubic.f(1), data[1]))
            self.assertTrue(np.isclose(cubic.df(0), data[2]))
            self.assertTrue(np.isclose(cubic.df(1), data[3]))

class TestBicubic(unittest.TestCase):

    def test_basic(self):
        pass

if __name__ == '__main__':
    unittest.main()
