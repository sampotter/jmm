import matplotlib.pyplot as plt; plt.ion()
import numpy as np

norm = np.linalg.norm

H0 = lambda sigma, L: 1 - 3*(sigma/L)**2 + 2*(sigma/L)**3
H1 = lambda sigma, L: 3*(sigma/L)**2 - 2*(sigma/L)**3
K0 = lambda sigma, L: sigma - 2*sigma**2/L + sigma**3/L**2
K1 = lambda sigma, L: -sigma**2/L + sigma**3/L**2

H0p = lambda sigma, L: -6*sigma/L**2 + 6*sigma**2/L**3
H1p = lambda sigma, L: 6*sigma/L**2 - 6*sigma**2/L**3
K0p = lambda sigma, L: 1 - 4*sigma/L + 3*sigma**2/L**2
K1p = lambda sigma, L: -2*sigma/L + 3*sigma**2/L**2

L = 0.1

Sigma = np.linspace(0, L)

plt.figure()
plt.plot(Sigma, H0(Sigma, L), label='H_0')
plt.plot(Sigma, H1(Sigma, L), label='H_1')
plt.plot(Sigma, K0(Sigma, L), label='K_0')
plt.plot(Sigma, K1(Sigma, L), label='K_1')
plt.xlim(0, L)
plt.legend()
plt.show()

plt.figure()
plt.plot(Sigma, H0p(Sigma, L), label="H_0'")
plt.plot(Sigma, H1p(Sigma, L), label="H_1'")
plt.plot(Sigma, K0p(Sigma, L), label="K_0'")
plt.plot(Sigma, K1p(Sigma, L), label="K_1'")
plt.xlim(0, L)
plt.legend()
plt.show()

phi = lambda sigma, xlam, tlam, that: xlam*H0(sigma) + xhat*H1(sigma) + tlam*K0(sigma) + that*K1(sigma)
phip = lambda sigma, xlam, tlam, that: xlam*H0p(sigma) + xhat*H1p(sigma) + tlam*K0p(sigma) + that*K1p(sigma)
