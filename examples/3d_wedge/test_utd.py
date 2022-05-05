import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import brentq
from scipy.special import modfresnelm

plt.ion()

ngrid = 201

n = 1.9
xsrc = np.array([0.5, 1, 0])
n_o = np.array([0, 1, 0])
t_o = np.array([1, 0, 0])
t_e = np.array([0, 0, 1])
k = 10

# image source for reflection
ximg = xsrc.copy()
ximg[1] *= -1 # reflect over o-face

alpha = (2 - n)*np.pi

def get_phi(x):
    phi = np.arctan2(x[1], x[0])
    phi = np.mod(phi, 2*np.pi)
    if phi < n*np.pi:
        return phi
    else:
        return np.nan

sp = np.linalg.norm(xsrc)
phip = get_phi(xsrc)

def tau(x):
    phi = get_phi(x)
    if phi < np.pi + phip:
        return np.linalg.norm(x - xsrc)
    elif phi < n*np.pi:
        return np.linalg.norm(x) + np.linalg.norm(xsrc)
    else:
        return np.nan

def tau_r(x):
    phi = get_phi(x)
    if phi < np.pi - phip:
        return np.linalg.norm(x - ximg)
    elif phi < n*np.pi:
        return np.linalg.norm(x) + np.linalg.norm(ximg)
    else:
        return np.nan

def s(x):
    if get_phi(x) < n*np.pi:
        return np.linalg.norm(x)
    else:
        return np.nan

def get_beta(sign, x):
    return get_phi(x) + sign*phip

def N(sign1, sign2, x):
    phi = get_phi(x)
    if phi < n*np.pi:
        beta = get_beta(sign2, x)
        return np.round((beta + sign1*np.pi)/(2*np.pi*n))
    else:
        return np.nan

def a(sign1, sign2, x):
    phi = get_phi(x)
    beta = get_beta(sign2, x)
    if phi < n*np.pi:
        return 2*np.cos((2*np.pi*n*N(sign1, sign2, x) - beta)/2)**2
    else:
        return np.nan

def get_L(x):
    phi = get_phi(x)
    if phi < n*np.pi:
        return sp*s(x)/(sp + s(x))
    else:
        return np.nan

def get_F(x):
    sqrtx = np.sqrt(x)
    return 2j*sqrtx*np.exp(1j*x)*modfresnelm(sqrtx)[0]

def F(sign1, sign2, x):
    phi = get_phi(x)
    if phi < n*np.pi:
        return get_F(k*get_L(x)*a(sign1, sign2, x))
    else:
        return np.nan

def D_old(sign1, sign2, x):
    phi = get_phi(x)
    beta = get_beta(sign2, x)
    if phi < n*np.pi:
        _  = -np.exp(-1j*np.pi/4)/(2*n*np.sqrt(2*np.pi*k))
        _ /=  np.tan((np.pi + sign1*beta)/(2*n))
        _ *=  F(sign1, sign2, x)
        return _
    else:
        return np.nan

def D(sign1, sign2, x):
    phi = get_phi(x)
    if phi < n*np.pi:
        beta = phi + sign2*phip
        s = np.linalg.norm(x)
        L = s*sp/(s + sp)
        N = np.round((beta + sign1*np.pi)/(2*np.pi*n))
        a = 2*np.cos((2*np.pi*n*N - beta)/2)**2
        _  = -np.exp(-1j*np.pi/4)/(2*n*np.sqrt(2*np.pi*k))
        _ /=  np.tan((np.pi + sign1*beta)/(2*n))
        _ *=  get_F(k*L*a)
        return _
    else:
        return np.nan

X, Y = np.meshgrid(np.linspace(-2, 2, ngrid), np.linspace(-2, 2, ngrid), indexing='xy')
x = np.array([(x, y, 0) for x, y in zip(X.ravel(), Y.ravel())])

Tau = np.array([tau(_) for _ in x]).reshape(X.shape)
Tau_r = np.array([tau_r(_) for _ in x]).reshape(X.shape)
S = np.array([s(_) for _ in x]).reshape(X.shape)

L = np.array([get_L(_) for _ in x]).reshape(X.shape)

Phi = np.array([get_phi(_) for _ in x]).reshape(X.shape)
Beta_p = np.array([get_beta( 1, _) for _ in x]).reshape(X.shape)
Beta_m = np.array([get_beta(-1, _) for _ in x]).reshape(X.shape)

Valid = (k*L > 1).astype(np.float64)
Valid[np.isnan(L)] = np.nan

N1 = np.array([N( 1, -1, _) for _ in x]).reshape(X.shape)
N2 = np.array([N(-1, -1, _) for _ in x]).reshape(X.shape)
N3 = np.array([N( 1,  1, _) for _ in x]).reshape(X.shape)
N4 = np.array([N(-1,  1, _) for _ in x]).reshape(X.shape)

a1 = np.array([a( 1, -1, _) for _ in x]).reshape(X.shape)
a2 = np.array([a(-1, -1, _) for _ in x]).reshape(X.shape)
a3 = np.array([a( 1,  1, _) for _ in x]).reshape(X.shape)
a4 = np.array([a(-1,  1, _) for _ in x]).reshape(X.shape)

F1 = np.array([F( 1, -1, _) for _ in x]).reshape(X.shape)
F2 = np.array([F(-1, -1, _) for _ in x]).reshape(X.shape)
F3 = np.array([F( 1,  1, _) for _ in x]).reshape(X.shape)
F4 = np.array([F(-1,  1, _) for _ in x]).reshape(X.shape)

D1 = np.array([D( 1, -1, _) for _ in x]).reshape(X.shape)
D2 = np.array([D(-1, -1, _) for _ in x]).reshape(X.shape)
D3 = np.array([D( 1,  1, _) for _ in x]).reshape(X.shape)
D4 = np.array([D(-1,  1, _) for _ in x]).reshape(X.shape)

Dh = D1 + D2 + D3 + D4

ui = (np.exp(-1j*k*Tau)/Tau)*(Phi < np.pi + phip).astype(np.float)
ur = (np.exp(-1j*k*Tau_r)/Tau_r)*(Phi < np.pi - phip).astype(np.float)
ud = Dh*np.exp(-1j*k*(sp + S))/(np.sqrt(S)*sp)
u = ui + ur + ud

inches_per_subplot = 1.75
subplot_n_rows = 7
subplot_n_cols = 4

def add_subplot(i, j, field, title, cmap=cc.cm.gouldian, vmin=None, vmax=None):
    plt.subplot(subplot_n_rows, subplot_n_cols, subplot_n_cols*i + j + 1)
    plt.imshow(field, extent=[-2, 2, -2, 2],
               interpolation='none', zorder=1, cmap=cmap,
               vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.scatter(xsrc[0], -xsrc[1], s=50, edgecolor='k',
                facecolor='red', zorder=3)
    plt.plot([0, 10*np.cos(phip + np.pi)], [0, -10*np.sin(np.pi - phip)],
             c='red', linestyle='--', linewidth=1, zorder=2)
    plt.plot([0, 10*np.cos(phip + np.pi)], [0, -10*np.sin(np.pi + phip)],
             c='red', linestyle='--', linewidth=1, zorder=2)
    plt.plot([xsrc[0], 0], [-xsrc[1], 0], linewidth=1, linestyle='--',
             c='r', zorder=2)
    plt.plot([20, 0, 20*np.cos(alpha)], [0, 0, 20*np.sin(alpha)],
             c='black', linestyle='-', linewidth=1, zorder=3)
    plt.scatter([0], [0], s=50, c='k', zorder=3)
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect('equal')

plt.figure(figsize=(4*inches_per_subplot + 2, 7*inches_per_subplot))
add_subplot(0, 0, Tau, r'$\tau_{\mathrm{i}}$')
add_subplot(0, 1, Tau_r, r'$\tau_{\mathrm{r}}$')
add_subplot(0, 2, S, r'$s$')
add_subplot(0, 3, L, r'$L$')
add_subplot(1, 0, Phi, r'$\phi$')
add_subplot(1, 1, Beta_p, r'$\beta^+$')
add_subplot(1, 2, Beta_m, r'$\beta^-$')
add_subplot(1, 3, Valid, r'$kL > 1$')
add_subplot(2, 0, N1, r'$N^+(\beta^-)$', vmin=-1, vmax=1)
add_subplot(2, 1, N2, r'$N^-(\beta^-)$', vmin=-1, vmax=1)
add_subplot(2, 2, N1, r'$N^+(\beta^+)$', vmin=-1, vmax=1)
add_subplot(2, 3, N2, r'$N^-(\beta^+)$', vmin=-1, vmax=1)
add_subplot(3, 0, a1, r'$a^+(\beta^-)$', vmin=0, vmax=2)
add_subplot(3, 1, a2, r'$a^-(\beta^-)$', vmin=0, vmax=2)
add_subplot(3, 2, a3, r'$a^+(\beta^+)$', vmin=0, vmax=2)
add_subplot(3, 3, a4, r'$a^-(\beta^+)$', vmin=0, vmax=2)
add_subplot(4, 0, np.real(F1), r'$\Re F^+(\beta^-)$')
add_subplot(4, 1, np.real(F2), r'$\Re F^-(\beta^-)$')
add_subplot(4, 2, np.real(F3), r'$\Re F^+(\beta^+)$')
add_subplot(4, 3, np.real(F4), r'$\Re F^-(\beta^+)$')
add_subplot(5, 0, np.real(D1), r'$\Re D^+(\beta^-)$')
add_subplot(5, 1, np.real(D2), r'$\Re D^-(\beta^-)$')
add_subplot(5, 2, np.real(D3), r'$\Re D^+(\beta^+)$')
add_subplot(5, 3, np.real(D4), r'$\Re D^-(\beta^+)$')
add_subplot(6, 0, np.real(ui), r'$\Re u_{\mathrm{i}}$', vmin=-1, vmax=1)
add_subplot(6, 1, np.real(ur), r'$\Re u_{\mathrm{r}}$', vmin=-1, vmax=1)
add_subplot(6, 2, np.real(ud), r'$\Re u_{\mathrm{d}}$', vmin=-1, vmax=1)
add_subplot(6, 3, np.real(u), r'$\Re u$', vmin=-1, vmax=1)
plt.tight_layout()
plt.savefig('utd.png', dpi=300)
plt.show()
