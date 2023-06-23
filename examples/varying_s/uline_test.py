import colorcet as cc
import itertools as it
import matplotlib.pyplot as plt
import numpy as np

norm = np.linalg.norm

v = np.array([0.1, -0.1])
s = lambda x: 1/(1 + v@x)
Ds = lambda x: -v/((1 + v@x)**2)

############################################################################
# WORKING LINE UPDATE FOLLOWS!
#

T0 = 0

x0 = np.array([1, 1])
xhat = np.array([0, 0])

s0 = s(x0)
shat = s(xhat)

L0 = lambda sigma: 1 - 3*sigma/L + 2*sigma**2/L**2
L1 = lambda sigma: 4*sigma/L - 4*sigma**2/L**2
L2 = lambda sigma: -sigma/L + 2*sigma**2/L**2

L0p = lambda sigma: -3/L + 4*sigma/L**2
L1p = lambda sigma: 4/L - 8*sigma/L**2
L2p = lambda sigma: -1/L + 4*sigma/L**2

phi = lambda sigma, xm: x0*L0(sigma) + xm*L1(sigma) + xhat*L2(sigma)
phip = lambda sigma, xm: x0*L0p(sigma) + xm*L1p(sigma) + xhat*L2p(sigma)
t = lambda sigma, xm: phip(sigma, xm)/norm(phip(sigma, xm))

L = np.linalg.norm(xhat - x0)

f = lambda xm: T0 + L*(s0*norm(phip(0, xm)) + 4*s(phi(L/2, xm))*norm(phip(L/2, xm)) + shat*norm(phip(L, xm)))/6
Df = lambda xm: (2/3)*(s0*t(0, xm) + L*Ds(phi(L/2, xm)) - shat*t(L, xm))

xm = (x0 + xhat)/2
niter = 0
while True:
    Dxm = Df(xm)
    alpha = 1
    niter_ = 0
    while alpha > 1e-4 and f(xm - alpha*Dxm) >= f(xm):
        alpha /= 2
        niter_ += 1
    xm -= alpha*Dxm
    if norm(Dxm) < 1e-12:
        break
    niter += 1
print(niter, norm(Dxm))

ngrid = 101

Sigma = np.linspace(0, L, ngrid)
Phi = np.array([phi(sigma, xm) for sigma in Sigma])

X, Y = np.meshgrid(np.linspace(0, 1, ngrid), np.linspace(0, 1, ngrid), indexing='xy')
F = np.empty_like(X)

for i, j in it.product(*(range(_) for _ in X.shape)):
    F[i, j] = f(np.array([X[i, j], Y[i, j]]))

plt.figure()
plt.contour(X, Y, F, levels=201, cmap=cc.cm.colorwheel, zorder=1, linewidths=1)
plt.scatter(*x0, s=50, c='k', zorder=2)
plt.scatter(*xhat, s=50, c='k', zorder=2)
plt.scatter(*xm, s=50, facecolor='white', edgecolor='k', zorder=2)
plt.plot(*Phi.T, zorder=2, linewidth=2, c='k')
plt.gca().set_aspect('equal')
plt.show()
