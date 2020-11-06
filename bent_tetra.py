import jmm
import matplotlib.pyplot as plt
import meshio
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

################################################################################
# Set up some test mesh stuff

# mesh_path = 'room.vtu'
mesh_path = 'sphere.vtu'

print(f'- loading {mesh_path}')
mesh = meshio.read(mesh_path)
verts = mesh.points.astype(np.float64).copy()
cells = mesh.cells_dict['tetra'].astype(np.uintp).copy()
del mesh

mesh = jmm.Mesh3(verts, cells)

def is_causal(i, f):
    X = verts[f] - verts[i]
    return (X@X.T >= 0).all()

def check(i):
    F = mesh.vf(i)
    for f in F:
        print(is_causal(i, f))

################################################################################
# Testing... do a quick plot

F = np.array([3, 0, 1, 2,
              3, 1, 2, 3,
              3, 2, 3, 0,
              3, 3, 0, 1])

i = 0
check(i)
h = np.sqrt(np.sum((verts[np.unique(mesh.vf(i)).ravel()] - verts[i])**2, axis=1)).mean()
plotter = pvqt.BackgroundPlotter()
for f in mesh.vf(i):
    V = verts[[i] + list(f.astype(int))]
    poly_data = pv.PolyData(V, F)
    color = 'white' if is_causal(i, f) else 'red'
    plotter.add_mesh(poly_data, opacity=0.1, color=color)
plotter.add_mesh(pv.Sphere(radius=h/10, center=verts[i]), color='red')

################################################################################

def count_noncausal_updates(i):
    num_noncausal = 0
    for f in mesh.vf(i):
        if not is_causal(i, f):
            num_noncausal += 1
    return num_noncausal

R = np.sqrt(np.sum(verts**2, axis=1))
noncausal = np.array(
    [count_noncausal_updates(i) for i in range(verts.shape[0])])

def find_opposite_cell(i, f):
    for j in f:
        for k in mesh.vc(j):
            C = cells[k]
            if i not in C and set(f).issubset(C):
                return k

def find_opposite_vert(i, f):
    for j in f:
        for k in mesh.vc(j):
            C = cells[k]
            if i not in C and set(f).issubset(C):
                return np.setdiff1d(C, f)[0]


i = 0
for f in mesh.vf(i):
    if not is_causal(i, f):
        j = find_opposite_vert(i, f)
        print(i, f, j)
