import colorcet as cc
import itertools as it
import jmm
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

n = 33
h = 2/(n-1)
shape = (1, n, n)
ind0 = (0, n - 1, 0)

dial = jmm.Dial(jmm.Stype.Constant, shape, h)

I, J, K = (
    _.flatten() for _ in np.meshgrid(
        np.arange(1, dtype=np.intc),
        np.arange(3*(n//4) + 1, dtype=np.intc),
        np.arange(3*(n//4) + 1, dtype=np.intc)
    )
)
inds = np.array([I, J, K], order='F')

dial.add_boundary_points(inds)
print(dial.state)
dial.add_point_source(ind0, 0)
dial.solve()
print(dial.state)
print(dial.T)
# dial.solve()
# print(dial.state)
# print(dial.T)
# print(dial.grad_T)

Y, Z = np.meshgrid(
    np.linspace(0, 2, n),
    np.linspace(0, 2, n),
    indexing='ij'
)

plt.figure()
plt.imshow(dial.T[0, :, :] - np.sqrt((Y - 2)**2 + Z**2))
plt.colorbar()
plt.show()

# OLD PYVISTA STUFF

# X, Y, Z = np.meshgrid(L, L, L)
# tau = np.sqrt(X**2 + Y**2 + Z**2)
# print(abs(tau - dial.T).max())
# def next():
#     dial.step()
#     I, J, K = np.where(dial.state == jmm.State.Valid.value)
#     grid = pv.StructuredGrid(I, J, K)
#     cell = h*np.array([
#         [0, 0, 0],
#         [1, 0, 0],
#         [1, 1, 0],
#         [0, 1, 0],
#         [0, 0, 1],
#         [1, 0, 1],
#         [1, 1, 1],
#         [0, 1, 1],
#     ])
#     cells = np.concatenate([
#         np.concatenate([[8], 8*i + np.arange(8)])
#         for i in range(len(I))
#     ])
#     cell_type = pv.vtk.VTK_HEXAHEDRON*np.ones(len(I), dtype=np.intc)
#     points = np.vstack(
#         tuple(
#             cell + h*np.array(ind)
#             for ind in zip(I, J, K)
#         )
#     )
#     grid = pv.UnstructuredGrid(cells, cell_type, points)
#     grid.cell_arrays['T'] = dial.T[dial.state == jmm.State.Valid.value]
#     grid.plot(interactive=True, show_edges=True)
