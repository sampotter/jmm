import sjs_eik
import unittest

class TestCubic(unittest.TestCase):

    def test_basic(self):
        values = [1, 2, 3, 4]
        heap_pos = 4 * [-1]

        value = lambda i: values[i]

        def setpos(i, pos):
            heap_pos[i] = pos

        heap = sjs_eik.Heap(4, value, setpos)

        self.assertTrue(heap.size == 0)
        self.assertIsNone(heap.front)

if __name__ == '__main__':
    unittest.main()
