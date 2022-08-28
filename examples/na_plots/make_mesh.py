#!/usr/bin/env python

import matplotlib.pyplot as plt
import meshpy.triangle as triangle
import numpy as np
import sys

if __name__ == '__main__':
    points = [(1, 1), (-1, 1), (-1, -1), (1, -1)]
    facets = [(0, 1), (1, 2), (2, 3), (3, 0)]

    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets)

    max_area = float(sys.argv[1])

    mesh = triangle.build(info, max_volume=max_area)

    V = np.array(mesh.points).astype(np.float64)
    V.tofile(sys.argv[2])

    F = np.array(mesh.elements).astype(np.uintp)
    F.tofile(sys.argv[3])

    plt.figure(figsize=(8, 8))
    plt.triplot(*V.T, triangles=F, linewidth=1, c='k')
    plt.title('Max area: %s' % sys.argv[1])
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(sys.argv[4])
