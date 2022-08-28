import matplotlib.pyplot as plt
import numpy as np

plt.ion()

T = np.linspace(0, 1)

H0 = 1 - 3*T**2 + 2*T**3
H1 = 3*T**2 - 2*T**3
K0 = T - 2*T**2 + T**3
K1 = -T**2 + T**3

plt.figure(figsize=(6.5, 3.5))
plt.plot(T, H0, label=r'$H_0(\lambda)$')
plt.plot(T, H1, label=r'$H_1(\lambda)$')
plt.plot(T, K0, label=r'$K_0(\lambda)$')
plt.plot(T, K1, label=r'$K_1(\lambda)$')
plt.xlabel(r'$\lambda$')
plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig('hermite_polys.pdf')
