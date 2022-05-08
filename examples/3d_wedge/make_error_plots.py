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

def get_h_for_each_mesh():
    h = []
    for path in paths:
        V = np.fromfile(path/'verts.bin', np.float64).reshape(-1, 3)
        C = np.fromfile(path/'cells.bin', np.uintp).reshape(-1, 4)
        dV = V[C[:, 1:]] - V[C[:, 0]].reshape(-1, 1, 3)
        h.append(np.mean(np.sqrt(np.sum(dV**2, axis=1))))
    return np.array(h)

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
        E_T.append(np.mean(abs(T[mask] - T_gt[mask]))/np.mean(abs(T_gt[mask])))
        E_DT.append(np.mean(np.apply_along_axis(norm, 1, DT[mask] - DT_gt[mask]))/
                    np.mean(np.apply_along_axis(norm, 1, DT_gt[mask])))
        E_D2T.append(
            np.mean([norm(H - H_gt) for H, H_gt in zip(D2T[mask], D2T_gt[mask])])/
            np.mean([norm(H_gt) for H_gt in D2T_gt[mask]]))
    return np.array(E_T), np.array(E_DT), np.array(E_D2T)

plt.figure(figsize=(9, 9))

p = dict()

for i, (selected_origin, field) in enumerate(
        it.product([None, 0, 1], ['direct', 'o_refl', 'n_refl'])):
    print(f'origin: {selected_origin}, field: {field}')

    h = get_h_for_each_mesh()
    E_T, E_DT, E_D2T = compute_rel_lp_errors_wrt_a(field, selected_origin)
    if E_T.size == 0 or E_DT.size == 0 or E_D2T.size == 0:
        continue

    I = np.argsort(h)
    h = h[I]
    E_T = E_T[I]
    E_DT = E_DT[I]
    E_D2T = E_D2T[I]

    fit_T = np.poly1d(np.polyfit(np.log(h), np.log(E_T), 1))
    fit_DT = np.poly1d(np.polyfit(np.log(h), np.log(E_DT), 1))
    fit_D2T = np.poly1d(np.polyfit(np.log(h), np.log(E_D2T), 1))

    p[selected_origin, field, 'T'] = fit_T[1]
    p[selected_origin, field, 'DT'] = fit_DT[1]
    p[selected_origin, field, 'D2T'] = fit_D2T[1]

    print(f'- T: {fit_T[1]:0.2f}, DT: {fit_DT[1]:0.2f}, D2T: {fit_D2T[1]:0.2f}')

    plt.subplot(3, 3, i + 1)

    plt.loglog(h, E_T, c='k')
    plt.loglog(h, np.exp(fit_T(np.log(h))), c='k', linestyle='--')
    plt.loglog(h, E_DT, c='b')
    plt.loglog(h, np.exp(fit_DT(np.log(h))), c='b', linestyle='--')
    plt.loglog(h, E_D2T, c='g')
    plt.loglog(h, np.exp(fit_D2T(np.log(h))), c='g', linestyle='--')
    plt.title(f'{field} (origin = {selected_origin})')

plt.tight_layout()
plt.savefig('error_plot.pdf')
plt.show()

# print table of least squares fits coefficients
f = open('error_table.tex', 'w')
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
f.close()
