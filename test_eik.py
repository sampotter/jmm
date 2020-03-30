import numpy as np
import sjs
import unittest

# TODO: definitely need to add some more tests here!

class TestEik(unittest.TestCase):
    def test_set_jet(self):
        shape = (2, 2)
        xymin = (0, 0)
        h = 1
        slow = sjs.get_constant_slowness_field2()
        eik = sjs.Eik(slow, shape, xymin, h)
        eik.add_trial(0, 0, sjs.Jet(1, 0, 0, 0))
        eik.add_trial(1, 0, sjs.Jet(1, 0, 0, 0))
        eik.add_trial(0, 1, sjs.Jet(1, 0, 0, 0))
        eik.add_trial(1, 1, sjs.Jet(1, 0, 0, 0))
        self.assertFalse(eik.can_build_cell(0, 0))

if __name__ == '__main__':
    unittest.main()
