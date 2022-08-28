import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

plt.ion()

def random_unit_vector():
    p = np.random.randn(3)
    return p/np.linalg.norm(p)

p1, p2, p3 = (random_unit_vector() for _ in range(3))
while p1@p2 <= 0 or p1@p3 <= 0 or p2@p3 <= 0:
    p1, p2, p3 = (random_unit_vector() for _ in range(3))
P = np.array([p1, p2, p3])

b = np.random.random(3)
b /= np.sum(b)
print(b)

q0 = b@P
q0 /= np.linalg.norm(q0)

def f(q):
    return np.dot(b, np.arccos(P@q)**2)/2

def get_next_q(q):
    dots = P@q
    angles = np.arccos(dots)
    dists = np.sqrt(1 - dots**2)
    q1 = q + (angles*b/dists)@P
    q1 /= np.linalg.norm(q1)
    return q1

c1, c2, c3 = 'blue', 'purple', 'magenta'

Q = [q0]
while True:
    q = Q[-1]
    q1 = get_next_q(q)
    Q.append(q1)
    if np.arccos(q@q1) < 1e-13:
        break

# plotter = pvqt.BackgroundPlotter()
# plotter.add_mesh(pv.Sphere(1, (0, 0, 0)), color='white', opacity=0.25)
# for p, c in zip([p1, p2, p3], [c1, c2, c3]):
#     plotter.add_mesh(pv.Sphere(0.025, p), color=c)
# for i, q in enumerate(Q):
#     c = cc.cm.gouldian(i/(len(Q) - 1))
#     plotter.add_mesh(pv.Sphere(0.025, q0), color=c)

arclength = [np.arccos(q@q1) for q, q1 in zip(Q[:-1], Q[1:])]
plt.figure()
plt.semilogy(np.arange(len(arclength)) + 1, arclength, marker='*')
plt.xlim(0, len(arclength))
plt.ylim(1e-16, 1e-1)
plt.show()

plt.figure()
plt.plot([f(q) for q in Q])
plt.show()
