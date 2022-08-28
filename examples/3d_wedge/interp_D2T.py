import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

from scipy.spatial.transform import Rotation

lam = np.array([0.075234048705463547, 0.17988900839685265])

H0 = np.array([[5.382782547084771, -1.6356222306139783, -2.3545392512490397],
               [-1.6356222306139783, 5.1576542146682858, -2.5221488227448852],
               [-2.3545392512490397, -2.5221488227448852, 3.2789820205623581]])
H1 = np.array([[7.8424779580662776, -1.442848627275954, 2.9280069470290777],
               [-1.442848627275954, 7.6688231558145672, 3.1095047487649903],
               [2.9280069470290777, 3.1095047487649903, 2.8909176112911799]])
H2 = np.array([[13.669125819660826, 4.5416000576547129, -2.0866093875754546],
               [4.5416000576547129, 4.2099619402316222, 5.1855775896009915],
               [-2.0866093875754546, 5.1855775896009915, 13.114126186393639]])

xsrc = np.array([1, -1, 0])

Xt = np.array([[1.0680329717582246, -0.92712405217419214, 0.10490764702032128],
               [0.95823711084613838, -1.0443516373066821, 0.090003829709602898],
               [0.97783990970831791, -0.94492842393774246, -0.025302286890365996]])

dXt = Xt[1:] - Xt[0]

R = Xt - xsrc
Tau = np.sqrt(np.sum(R**2, axis=1))
tau0, tau1, tau2 = Tau
T = np.diag(1/Tau)@R
t0, t1, t2 = T

x = np.array([1.1357843269562427, -0.80921907234928969, 0.27403999490762704])
r = x - xsrc
tau = np.linalg.norm(r)
t = r/tau

def get_H(y):
    r = y - xsrc
    return (np.eye(r.size) - np.outer(r, r)/np.dot(r, r))/np.linalg.norm(r)

xlam = Xt[0] + lam@dXt
Llam = np.linalg.norm(x - xlam)

rlam = xlam - xsrc
tlam = rlam/np.linalg.norm(rlam)

Hlam_gt = get_H(xlam)
Hlam = np.array([[7.0584574659747963, -0.50990469986772569, -1.9089142757827773],
                 [-0.50990469986772569, 5.1761001975803627, -0.7119214623424831],
                 [-1.9089142757827773, -0.7119214623424831, 5.0190206953258727]])

H = np.array([[2.4864662292277204, -0.10889688102695839, -0.32671008840866256],
              [-0.10889688102695834, 2.2336641427251931, -0.16169913509496675],
              [-0.32671008840866256, -0.16169913509496675, 2.1399973207655236]])
H_gt = get_H(x)

Q0, D0 = np.linalg.svd(H0)[:2]
Q1, D1 = np.linalg.svd(H1)[:2]
Q2, D2 = np.linalg.svd(H2)[:2]

if Q0[:, 2]@t0 < 0: Q0[:, 2] *= -1
if Q1[:, 2]@t1 < 0: Q1[:, 2] *= -1
if Q2[:, 2]@t2 < 0: Q2[:, 2] *= -1

Qlam_gt, Dlam_gt = np.linalg.svd(Hlam_gt)[:2]
Qlam, Dlam = np.linalg.svd(Hlam)[:2]
Q, D = np.linalg.svd(H)[:2]
Q_gt, D_gt = np.linalg.svd(H_gt)[:2]

# try new method

Dlam_new = D0 + lam[0]*(D1 - D0) + lam[1]*(D2 - D0)

q0 = Rotation.from_matrix(Q0).as_quat()
q1 = Rotation.from_matrix(Q1).as_quat()
q2 = Rotation.from_matrix(Q2).as_quat()

def qinterp(qs, ws):
    q0 = sum(w*q for w, q in zip(qs, ws))
    q0 /= np.linalg.norm(q0)
    return q0

qs = [q0, q1, q2]
ws = [1 - lam[0] - lam[1], lam[0], lam[1]]

qlam = qinterp(qs, ws)

Qlam_new = Rotation.from_quat(qlam).as_matrix()

Hlam_new = Qlam_new@np.diag(Dlam_new)@Qlam_new.T

D_new = Dlam_new/(1 + Llam*Dlam_new)
H_new = Qlam_new@np.diag(D_new)@Qlam_new.T

h = 0.005
plt = pvqt.BackgroundPlotter()
plt.add_mesh(pv.Sphere(h, xsrc), color='white')
plt.add_mesh(pv.Sphere(h, xlam), color='green')
plt.add_mesh(pv.Arrow(xlam, tlam, scale=10*h), color='cyan')
plt.add_mesh(pv.Sphere(h, x), color='pink')
for x_, Q_, in [(x0, Q0), (x1, Q1), (x2, Q2)]:
    for c, q_ in zip(['red', 'orange', 'yellow'], Q_.T):
        plt.add_mesh(pv.Arrow(x_, q_, scale=10*h), color=c)
# Qlam
for c, q_ in zip(['blue', 'purple', 'cyan'], Qlam_new.T):
    plt.add_mesh(pv.Arrow(xlam, q_, scale=7*h), color=c)
for c, q_ in zip(['black', 'gray', 'white'], Qlam.T):
    plt.add_mesh(pv.Arrow(xlam, q_, scale=7*h), color=c)
for c, q_ in zip(['green', 'chartreuse', 'yellow'], Qlam_gt.T):
    plt.add_mesh(pv.Arrow(xlam, q_, scale=7*h), color=c)
# Q
for c, q_ in zip(['blue', 'purple', 'cyan'], Qlam_new.T):
    plt.add_mesh(pv.Arrow(x, q_, scale=7*h), color=c)
for c, q_ in zip(['black', 'gray', 'white'], Q.T):
    plt.add_mesh(pv.Arrow(x, q_, scale=7*h), color=c)
for c, q_ in zip(['green', 'chartreuse', 'yellow'], Q_gt.T):
    plt.add_mesh(pv.Arrow(x, q_, scale=7*h), color=c)
