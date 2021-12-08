import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import brentq

plt.ion()

T = np.linspace(0, 1, 101)

H0 = lambda t: 1 - 3*t**2 + 2*t**3
H1 = lambda t: 3*t**2 - 2*t**3
K0 = lambda t: t - 2*t**2 + t**3
K1 = lambda t: -t**2 + t**3

H0p = lambda t: -6*t + 6*t**2
H1p = lambda t: 6*t - 6*t**2
K0p = lambda t: 1 - 4*t + 3*t**2
K1p = lambda t: -2*t + 3*t**2

def get_exact_deviation(xhat, yhat, x0, y0, x1, y1):
    dx, dy = x1 - x0, y1 - y0
    tau0 = np.sqrt(x0**2 + y0**2)
    tau1 = np.sqrt(x1**2 + y1**2)
    gradtau0 = np.array([x0, y0])/tau0
    gradtau1 = np.array([x1, y1])/tau1
    Dtau0 = np.array([dx, dy])@gradtau0
    Dtau1 = np.array([dx, dy])@gradtau1
    Htau = lambda t: tau0*H0(t) + tau1*H1(t) + Dtau0*K0(t) + Dtau1*K1(t)
    Htaup = lambda t: tau0*H0p(t) + tau1*H1p(t) + Dtau0*K0p(t) + Dtau1*K1p(t)
    L = lambda t: np.sqrt((xhat - x0 - t*dx)**2 + (yhat - y0 - t*dy)**2)
    Lp = lambda t: -np.array([dx, dy])@np.array([xhat - x0 - t*dx, yhat - y0 - t*dy])/L(t)
    f = lambda t: Htau(t) + L(t)
    fp = lambda t: Htaup(t) + Lp(t)
    if np.sign(fp(0)) == np.sign(fp(1)):
        t = 0 if f(0) < f(1) else 1
    else:
        t = brentq(fp, 0, 1, xtol=1e-15, rtol=1e-15)
    return Htau(t) - np.sqrt((x0 + t*dx)**2 + (y0 + t*dy)**2), t

def get_max_f4(x, y, h):
    f4 = -(3*h**4)/(y**2+(x+T*h)**2)**(3/2) \
        +(18*h**4*(x+T*h)**2)/(y**2+(x+T*h)**2)**(5/2) - \
        (15*h**4*(x+T*h)**4)/(y**2+(x+T*h)**2)**(7/2)
    return abs(f4).max()

X = np.linspace(-0.1, 0.1, 2001)
y = -0.1
h = 0.0005

exact_dev = np.array([
    get_exact_deviation(x + h, y + h, x, y, x + h, y)
    if x >= 0 else
    get_exact_deviation(x - h, y + h, x, y, x - h, y)
    for x in X
])

max_f4 = np.array([get_max_f4(x, y, h if x >= 0 else -h) for x in X])

plt.figure()
plt.axhline(y=h**3, linestyle='--', linewidth=1, c='r', zorder=1)
plt.axhline(y=h**4, linestyle='--', linewidth=1, c='b', zorder=1)
plt.semilogy(X, max_f4, c='k', zorder=2)
plt.semilogy(X, abs(exact_dev[:, 0]), c='gray', zorder=2)
plt.ylabel(r'$\max_{0 \leq \lambda \leq 1} |f^{(4)}(\lambda)|$')
plt.xlabel(r'$x$')
plt.xlim(X.min(), X.max())
plt.show()
