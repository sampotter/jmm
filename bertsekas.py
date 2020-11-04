import colorcet as cc
import matplotlib.pyplot as plt
import meshio
import numpy as np
import jmm

plt.ion()

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

################################################################################
# Set up some test Bernstein polynomial stuff

r = lambda x: np.linalg.norm(x)
Dr = lambda x: x/r(x)

X = np.eye(3)

f = np.array([r(x) for x in X])
Df = np.array([Dr(x) for x in X])

bb = jmm.Bb3Tri(f, Df, X)

################################################################################
# Set up cost functions and derivatives

xhat = np.array([2, 1, 1.5], dtype=np.float64)

def F(x):
    return bb.f(x) + np.linalg.norm(xhat - x)

A = np.array([
    [ 1, -1,  0],
    [ 0,  1, -1],
    [-1,  0,  1]
], dtype=np.float64)

def DF(x):
    g = A@np.array([bb.Df(x, a) for a in A], dtype=np.float64)
    return g + (x - xhat)/np.linalg.norm(x - xhat)

def D2F(x):
    H1 = np.array([
        [bb.D2f(x, a1, a2) for a1 in A]
        for a2 in A
    ], dtype=np.float64)
    H1 = A@H1@A.T
    r = np.linalg.norm(x - xhat)
    g = (x - xhat)/r
    H2 = (np.eye(x.size) - np.outer(g, g))/r
    return H1 + H2

################################################################################
# Test out Bertsekas's optimization algorithm

# Optimization algorithm parameters
beta = 0.25
sigma = 0.1
eps = 0.1
ftol = 5*np.finfo(np.float64).eps
xtol = 5*np.finfo(np.float64).eps

# Initial iterate
x0 = np.random.random((2,))
x0 = np.concatenate([x0, [1 - x0.sum()]])
while any(_ < 0 for _ in x0):
    x0 = np.random.random((2,))
    x0 = np.concatenate([x0, [1 - x0.sum()]])

# Do one step of Bertsekas's iteration
def step(xk):
    ik = np.argmax(xk)

    Tk = np.eye(xk.size)
    Tk[ik, :] = 1
    Tk = np.linalg.inv(Tk)

    yk = np.linalg.solve(Tk, xk)

    fk = F(xk)
    gk = Tk.T@DF(xk)
    Hk = Tk.T@D2F(xk)@Tk

    pk = np.maximum(0, yk - gk/np.diag(Hk))

    I = np.setdiff1d(np.arange(xk.size), [ik])
    wk = np.linalg.norm((yk - pk)[I])

    epsk = min(eps, wk)

    Ikp = np.where((0 <= yk) & (yk <= epsk) & (gk > 0))[0]

    for i in Ikp:
        for j in Ikp:
            if i == j:
                continue
            Hk[i, j] = 0

    pk = np.linalg.solve(Hk, gk)

    def get_yk1(alpha):
        yk1 = np.maximum(0, yk - alpha*pk)
        yk1[ik] = 1
        return yk1

    Ikpc = np.setdiff1d(np.arange(xk.size), Ikp)

    m = 0
    yk1 = get_yk1(beta**m)
    fk1 = F(Tk@yk1)
    lhs = fk - fk1
    rhs = sigma*(gk[Ikpc]@pk[Ikpc]*beta**m + gk[Ikp]@(yk - yk1)[Ikp])
    while abs(fk - fk1)/abs(fk) > np.finfo(np.dtype(fk)).eps and \
          ((Tk@yk1 < 0).any() or lhs < rhs):
        m += 1
        yk1 = get_yk1(beta**m)
        fk1 = F(Tk@yk1)
        lhs = fk - fk1
        rhs = sigma*(gk[Ikpc]@pk[Ikpc]*beta**m + gk[Ikp]@(yk - yk1)[Ikp])

    xk1 = Tk@yk1

    return xk1, m

# Define a function that carries out one step of the main loop of the
# optimization algorithm
xs, fs, ms = [x0], [F(x0)], []
def iterate():
    x, m = step(xs[-1])
    xs.append(x)
    ms.append(m)
    fs.append(F(x))
    df, dx = fs[-2] - fs[-1], xs[-2] - xs[-1]
    print('|f| = %g, |dx| = %g, m = %d' % (df, np.linalg.norm(dx), m))
    return df, dx

# Solve the optimization problem to the prescribed tolerance
df, dx = iterate()
while abs(df) > ftol and abs(dx).max() > xtol:
    df, dx = iterate()
