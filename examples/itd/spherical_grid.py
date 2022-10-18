import colorcet as cc
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import scipy.special

m = 20
p = 10
n = m*p

phi_lin = np.linspace(0, 2*np.pi, 2*n + 1)
# theta_lin = scipy.special.roots_legendre(n - 1)[0]
# theta_lin = np.concatenate([[-1], theta_lin, [1]])
# theta_lin = np.arcsin(theta_lin)
theta_lin = np.arcsin(np.linspace(-1, 1, 101))

phi, theta = np.meshgrid(phi_lin, theta_lin, indexing='ij')

x = np.cos(phi)*np.sin(theta)
y = np.sin(phi)*np.sin(theta)
z = np.cos(theta)

f = x**2 + 2*y**2 - 0.5*z**2 + 1.75*x*y*z - 4*y*z + x

vmax = abs(f).max()
vmin = -vmax
levels = np.linspace(vmin, vmax, 11)

plt.figure(figsize=(12, 3.5))
plt.contourf(phi, theta, f, levels, cmap=cc.cm.coolwarm)
plt.colorbar()
plt.xticks(phi_lin[::m])
plt.yticks(theta_lin[::m])
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\theta$')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.show()
