import sjs_eik

values = [1, 2, 3, 4]
heap_pos = 4 * [-1]

def value(i):
    return values[i]

def setpos(i, pos):
    heap_pos[i] = pos

heap = sjs_eik.Heap(4, value, setpos)

heap.insert(0)
