import itertools as it
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

paths = list(Path('.').glob('n*_a*_rfac*_phip*_sp*_w*_h*'))

def a_from_path(path):
    a_str = str(path).split('_')[1]
    return float(a_str[1:])

A = [a_from_path(path) for path in paths]

def load_data(path, field, selected_origin=None):
    # load numerical data
    jet = np.fromfile(path/f'{field}_jet.bin').reshape(-1, 4)
    T = jet[:, 0]
    DT = jet[:, 1:]
    D2T = np.fromfile(path/f'{field}_hess.bin').reshape(-1, 3, 3)

    # load groundtruth data
    jet_gt = np.fromfile(path/f'{field}_jet_gt.bin').reshape(-1, 13)
    T_gt = jet_gt[:, 0]
    DT_gt = jet_gt[:, 1:4]
    D2T_gt = jet_gt[:, 4:].reshape(-1, 3, 3)

    if selected_origin is not None:
        origin = np.fromfile(path/f'{field}_origin.bin')
        return T, T_gt, DT, DT_gt, D2T, D2T_gt, origin
    else:
        return T, T_gt, DT, DT_gt, D2T, D2T_gt

def get_h_and_N_for_each_mesh():
    h = []
    N = []
    for path in paths:
        V = np.fromfile(path/'verts.bin', np.float64).reshape(-1, 3)
        C = np.fromfile(path/'cells.bin', np.uintp).reshape(-1, 4)
        dV = V[C[:, 1:]] - V[C[:, 0]].reshape(-1, 1, 3)
        h.append(np.mean(np.sqrt(np.sum(dV**2, axis=1))))
        N.append(V.shape[0])
    return np.array(h), np.array(N)

def compute_rel_lp_errors_wrt_a(field, selected_origin=None, p=1):
    norm = lambda x: np.linalg.norm(x, ord=p)
    E_T, E_DT, E_D2T = [], [], []
    for path in paths:
        if selected_origin is not None:
            if selected_origin not in {0, 1}:
                raise RuntimeError('selected_origin must equal 0 or 1 if passed')
            T, T_gt, DT, DT_gt, D2T, D2T_gt, origin = \
                load_data(path, field, selected_origin)
            mask = origin < 0.5 if selected_origin == 0 else origin >= 0.5
        else:
            T, T_gt, DT, DT_gt, D2T, D2T_gt = load_data(path, field)
            mask = np.ones(T.shape, dtype=np.bool_)
        if mask.sum() == 0:
            break
        E_T.append(np.nanmean(abs(T[mask] - T_gt[mask]))/np.nanmean(abs(T_gt[mask])))
        E_DT.append(np.nanmean(np.apply_along_axis(norm, 1, DT[mask] - DT_gt[mask]))/
                    np.nanmean(np.apply_along_axis(norm, 1, DT_gt[mask])))
        E_D2T.append(
            np.nanmean([norm(H - H_gt) for H, H_gt in zip(D2T[mask], D2T_gt[mask])])/
            np.nanmean([norm(H_gt) for H_gt in D2T_gt[mask]]))
    return np.array(E_T), np.array(E_DT), np.array(E_D2T)

plt.figure(figsize=(10, 10))

p = dict()

field_tex_str = {
    'direct': 'Direct',
    'o_refl': r'$o$-refl.',
    'n_refl': r'$n$-refl.'
}

origin_tex_str = {
    None: '',
    0: ' ($\mathtt{origin} < 1/2$)',
    1: ' ($\mathtt{origin} \geq 1/2$)',
}

for i, (selected_origin, field) in enumerate(
        it.product([None, 0, 1], ['direct', 'o_refl', 'n_refl'])):
    print(f'origin: {selected_origin}, field: {field}')

    h, N = get_h_and_N_for_each_mesh()
    E_T, E_DT, E_D2T = compute_rel_lp_errors_wrt_a(field, selected_origin)
    if E_T.size == 0 or E_DT.size == 0 or E_D2T.size == 0:
        continue

    I = np.argsort(h)
    h = h[I]
    N = N[I]
    E_T = E_T[I]
    E_DT = E_DT[I]
    E_D2T = E_D2T[I]

    if np.isnan(E_T).any() or np.isnan(E_DT).any() or np.isnan(E_D2T).any():
        continue

    fit_T = np.poly1d(np.polyfit(np.log(h), np.log(E_T), 1))
    fit_DT = np.poly1d(np.polyfit(np.log(h), np.log(E_DT), 1))
    fit_D2T = np.poly1d(np.polyfit(np.log(h), np.log(E_D2T), 1))

    p[selected_origin, field, 'T'] = fit_T[1]
    p[selected_origin, field, 'DT'] = fit_DT[1]
    p[selected_origin, field, 'D2T'] = fit_D2T[1]

    print(f'- T: {fit_T[1]:0.2f}, DT: {fit_DT[1]:0.2f}, D2T: {fit_D2T[1]:0.2f}')

    plt.subplot(3, 3, i + 1)

    plt.loglog(h, E_T, c='k', marker='.', label=r'$\|\tau - T\|_1/\|\tau\|_1$')
    plt.loglog(h, np.exp(fit_T(np.log(h))), c='k', marker='.', linestyle='--')
    plt.loglog(h, E_DT, c='b', marker='.', label=r'$\|\nabla\tau - T\|_1/\|\nabla\tau\|_1$')
    plt.loglog(h, np.exp(fit_DT(np.log(h))), c='b', marker='.', linestyle='--')
    plt.loglog(h, E_D2T, c='g', marker='.', label=r'$\|\nabla^2\tau - T\|_1/\|\nabla^2\tau\|_1$')
    plt.loglog(h, np.exp(fit_D2T(np.log(h))), c='g', marker='.', linestyle='--')
    plt.title(f'{field_tex_str[field]}{origin_tex_str[selected_origin]}')
    if i == 0:
        plt.legend(fontsize=8, loc='upper left')
    if i == 6:
        plt.xlabel(r'$h$')

plt.tight_layout()
plt.savefig('error-plot.pdf')
plt.show()

# print table of least squares fits coefficients
f = open('error-table.tex', 'w')
print(r'\centering', file=f)
print(r'''\begin{tabular}{c|ccc|ccc|ccc}
  & \multicolumn{3}{c|}{All} & \multicolumn{3}{c|}{Direct} & \multicolumn{3}{c}{Diffracted} \\
  & $T$ & $\nabla T$ & $\nabla^2 T$ & $T$ & $\nabla T$ & $\nabla^2 T$ & $T$ & $\nabla T$ & $\nabla^2 T$ \\
  \midrule''', file=f)
for field in ['direct', 'o_refl', 'n_refl']:
    field_str = {
        'direct': 'Direct',
        'o_refl': 'Reflection ($o$-face)',
        'n_refl': 'Reflection ($n$-face)'
    }[field]
    print(rf'  {field_str}', file=f, end='')
    for selected_origin, deriv in it.product([None, 1, 0], ['T', 'DT', 'D2T']):
        if (selected_origin, field, deriv) in p:
            print(rf' & {p[selected_origin, field, deriv]:0.2f}', file=f, end='')
        else:
            print(rf' &', file=f, end='')
    print(r' \\', file=f)
print(r'\end{tabular}', file=f)
h_tex_str = ', '.join([f'{_:0.03g}' for _ in h])
N_tex_str = ', '.join([str(_) for _ in N])
print(r"\caption{Least squares fits for the relative $\ell_1$ norm error for $\Eik$, $\grad\Eik$, and $\hess\Eik$, measured for $h = {" + h_tex_str + r"}$ ($N = {" + N_tex_str + r"}$). First column (``All''): the error measured in all of $\domain_h$. Second column (``Direct''): for $\mx \in \calV_h$ with $\mathtt{origin}(\mx) > 0.5$. Third column (``Diffracted''): for $\mx \in \calV_h$ with $\mathtt{origin}(\mx) \leq 0.5$.}\label{table:errors}", file=f)
f.close()
