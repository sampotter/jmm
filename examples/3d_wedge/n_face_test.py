import matplotlib.pyplot as plt
import numpy as np

n = 0.25
phip = 4/5*np.pi
dphip = (n - 1)*np.pi + phip
phii = phip - 2*dphip
phivalid = n*np.pi + dphip

xp = np.cos(-phip)
yp = np.sin(-phip)

xi = np.cos(-phii)
yi = np.sin(-phii)

plt.figure()
plt.plot([0, xp], [0, yp], marker='.', linewidth=1, c='k', linestyle='--')
plt.text(xp, yp, r'$x_p$')
plt.plot([4*np.cos(n*np.pi), 0, 4], [4*np.sin(n*np.pi), 0, 0], linewidth=3, c='k', linestyle='-')
plt.plot([0, 4*np.cos((n - 1)*np.pi)], [0, 4*np.sin((n - 1)*np.pi)], linewidth=1, c='k', linestyle=':')
plt.plot([0, xi], [0, yi], marker='.', linewidth=1, c='k', linestyle='--')
plt.text(xi, yi, r'$x_i$')
plt.plot([0, 4*np.cos(phivalid)], [0, 4*np.sin(phivalid)], linewidth=1, c='blue', linestyle=':')
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.gca().set_aspect('equal')
plt.show()
