#!/usr/bin/env python

import meshio
import numpy as np

verts = np.array([
    [-1, -1, -1],
    [-1, -1,  1],
    [-1,  1, -1],
    [-1,  1,  1],
    [ 1, -1, -1],
    [ 1, -1,  1],
    [ 1,  1, -1],
    [ 1,  1,  1]
])

faces = np.array([
    [0, 1, 2],
    [2, 1, 3],
    [0, 2, 4],
    [4, 2, 6],
    [0, 1, 4],
    [4, 1, 5],
    [7, 6, 5],
    [5, 6, 4],
    [7, 5, 3],
    [3, 5, 1],
    [7, 6, 3],
    [3, 6, 2]
])

if __name__ == '__main__':
    meshio.write_points_cells('box.obj', verts, [(('triangle'), faces)])
    meshio.write_points_cells('box.off', verts, [(('triangle'), faces)])
