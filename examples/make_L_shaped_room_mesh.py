#!/usr/bin/env python

import meshio
import numpy as np

verts = np.array([
    [1, 0, 0],
    [2, 0, 0],
    [2, 2, 0],
    [0, 2, 0],
    [0, 1, 0],
    [1, 1, 0],
    [1, 0, 1],
    [2, 0, 1],
    [2, 2, 1],
    [0, 2, 1],
    [0, 1, 1],
    [1, 1, 1],
])

faces = np.array([
    # Walls
    [0, 1, 7],
    [0, 7, 6],
    [1, 2, 8],
    [1, 8, 7],
    [2, 3, 9],
    [2, 9, 8],
    [3, 4, 10],
    [3, 10, 9],
    [4, 5, 11],
    [4, 11, 10],
    [5, 0, 6],
    [5, 6, 11],
    # Floor
    [0, 1, 5],
    [1, 2, 5],
    [4, 5, 2],
    [4, 2, 3],
    [6, 7, 11],
    [7, 8, 11],
    [10, 11, 8],
    [10, 8, 9]
])

meshio.write_points_cells('L.off', verts, [(('triangle'), faces)])
meshio.write_points_cells('L.obj', verts, [(('triangle'), faces)])
