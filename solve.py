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
# States

VALID = 0
TRIAL = 1
FAR = 2

################################################################################

N = verts.shape[0]

T = np.empty((N,), dtype=np.float64)
T[...] = np.inf

gradT = np.empty((N, 3), dtype=np.float64)
gradT[...] = np.nan

R = np.sqrt(np.sum(verts**2, axis=1))
I = np.where(R < 0.3)[0]
T[I] = R[I]
gradT[I] = verts[I]/R[I].reshape(I.size, 1)

State = np.empty((N,), dtype=np.uint8)
State[...] = FAR
State[I] = VALID

I_trial = []
for i in I:
    nb = mesh.vv(i)
    far_nb = any(State[j] == FAR for j in nb)
    valid_nb = any(State[j] == VALID for j in nb)
    if far_nb and valid_nb:
        I_trial.append(i)
for i in I_trial:
    State[i] = TRIAL

################################################################################

A = np.array([
    [ 1, -1,  0],
    [ 0,  1, -1],
    [-1,  0,  1]
], dtype=np.float64)

def update(i, f):
    x = verts[i]
    X = verts[f]

    r = lambda y: np.linalg.norm(y - x)
    Dr = lambda y: (y - x)/r(y)
    D2r = lambda y: (np.eye(x.size) - np.outer(Dr(y), Dr(y)))/r(y)

    T = jmm.Bb3Tri(T[f], gradT[f], X)

    def F(th):
        y = th@X
        return bb.f(y) + r(y)

    def DF(th):
        DT = A@np.array([bb.DF(th, a) for a in A], dtype=np.float64)
        return DT + Dr(th@X)

    def D2F(th):
        D2T = np.array([[bb.D2f(x, a1, a2) for a1 in A] for a2 in A])
        D2T = A@D2T@A.T
        return D2T + A@D2r(th@X

################################################################################

I = np.where(State == TRIAL)[0]
i1 = I[np.argmin(T[I])]
State[i1] = VALID
for i in mesh.vv(i1):
    if State[i] == FAR:
        State[i] = TRIAL
for i in mesh.vv(i1):
    if State[i] == VALID:
        continue
    for f in mesh.vf(i):
        if i1 in f and (State[f] == VALID).all():
            update(i, f)
