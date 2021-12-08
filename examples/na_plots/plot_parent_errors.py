import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np

plt.ion()

N = 257
h = 2/(N-1)
i0, j0 = N//2, N//2

rfac = 0.3
dirpath = f'out_rfac{rfac}'
path = f'{dirpath}/sol_infos_for_par_{N}.txt'

with open(path, 'r') as f:
    lines = f.readlines()

data = np.array([[float(_) for _ in line.split()] for line in lines])
# data = data[np.isfinite(data).all(1)]
D, I, J, LamT, Lamtau, Lamstar, E0, E0_term1, E0_term2, E0_term3, has_parent = data.T
D, I, J = D.astype(int), I.astype(int), J.astype(int)
has_parent = has_parent.astype(bool)

############################################################################
# make parent error plots

def make_plot(E_values, title):
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')

    for d in range(D.max()):
        mask = D == d
        Y = abs(E_values[mask])
        Y = np.maximum(np.finfo(np.float64).eps, Y)
        X = abs(h*np.sqrt((I[mask] - i0)**2 + (J[mask] - j0)**2))
        c = cc.cm.colorwheel(float(d/D.max()))
        plt.scatter(X, Y, s=5, color=c, zorder=d+2)
        plt.scatter(X[~has_parent[mask]], Y[~has_parent[mask]], s=15,
                    edgecolors='black', facecolors='none', zorder=d+3)
        if d == 0:
            plt.scatter(X, Y, marker='x', color='k', zorder=D.max() + 10)
    plt.title(r'%s' % (title))
    plt.xlabel(r'$\sqrt{x^2 + y^2}$')

    ax1 = plt.gca()

    plt.twinx()
    plt.gca().set_yscale('log')
    plt.axhline(y=h, linestyle='--', c='k', linewidth=1, zorder=1)
    plt.axhline(y=h**2, linestyle='--', c='k', linewidth=1, zorder=1)
    plt.axhline(y=h**3, linestyle='--', c='k', linewidth=1, zorder=1)
    plt.axhline(y=h**4, linestyle='--', c='k', linewidth=1, zorder=1)
    plt.yticks([h, h**2, h**3, h**4])
    plt.gca().set_yticklabels([r'$h$', r'$h^2$', r'$h^3$', r'$h^4$'])
    plt.xticks(np.linspace(0.1, 1.01, 10))
    plt.axvline(x=rfac, linestyle='-', c='gray', linewidth=1, zorder=1)

    ax2 = plt.gca()

    ymin1, ymax1 = ax1.get_ylim()
    ymin2, ymax2 = ax2.get_ylim()

    ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)

    ax1.set_ylim(1e-17, 1e2)
    ax2.set_ylim(1e-17, 1e2)

    ax1.set_xlim(0.9e-1, 1)
    ax2.set_xlim(0.9e-1, 1)

plt.figure(figsize=(12, 8))
plt.subplot(2, 2, 1); make_plot(E0, r'$\epsilon^0$')
plt.subplot(2, 2, 2); make_plot(E0_term1, r'$|F_T(\lambda_T) - F_\tau(\lambda_\tau)|$')
plt.subplot(2, 2, 3); make_plot(E0_term2, r'$|F_\tau(\lambda_T) - F_\tau(\lambda_\tau)|$')
plt.subplot(2, 2, 4); make_plot(E0_term3, r'$|F_\tau(\lambda_\tau) - f(\lambda_{\mathsf{opt}})|$')
plt.tight_layout()
# plt.show()
plt.savefig(f'{dirpath}/parent_plot_{N}.pdf')

# make domain of dependence plot

x, y = np.meshgrid(np.linspace(-1, 1, 201), np.linspace(-1, 1, 201), indexing='ij')
z = np.sqrt(x**2 + y**2)

plt.figure(figsize=(6.5, 6.5))
plt.axhline(y=0, linestyle='-', linewidth=1, c='gray', zorder=1)
plt.axvline(x=0, linestyle='-', linewidth=1, c='gray', zorder=1)
plt.plot([-1, 1], [1, -1], linestyle='--', linewidth=1, c='gray', zorder=1)
plt.plot([1, -1], [1, -1], linestyle='--', linewidth=1, c='gray', zorder=1)
plt.contour(x, y, z, zorder=2, linewidths=1)
plt.scatter(h*(I - i0), h*(J - j0), s=5, c=D, cmap=cc.cm.colorwheel, zorder=3)
plt.scatter(h*(I[~has_parent] - i0), h*(J[~has_parent] - j0), s=15,
            edgecolors='k', facecolors='none', zorder=4)
plt.scatter(h*(I[D == 0] - i0), h*(J[D == 0] - j0), marker='x', c='k', zorder=4)
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.xlim(-1, 0.1)
plt.ylim(-0.1, 1)
plt.gca().set_aspect('equal')
plt.tight_layout()
# plt.show()
plt.savefig(f'{dirpath}/domain_of_dependence_{N}.pdf')
