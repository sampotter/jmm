import matplotlib.pyplot as plt

eik.solve()
R = verts - verts[indsrc]
tau = np.sqrt(np.sum(R**2, axis=1))
Dtau = R/tau.reshape(tau.shape[0], 1)
T = np.array([jet[0] for jet in eik.jet])
DT = np.array([(jet[1], jet[2], jet[3]) for jet in eik.jet])
E = T - tau
ED = np.sqrt(np.sum((DT - Dtau)**2, axis=1))
angle = np.rad2deg(np.arccos((DT*Dtau).sum(1)))

plt.figure(figsize=(9, 5))
plt.subplot(1, 2, 1)
plt.scatter(T, angle, s=2, c='k')
plt.xlabel(r'$\tau(x) = \|x\|$')
plt.ylabel(r'$\angle (\nabla T, \nabla \tau)$ [Deg.]')
plt.subplot(1, 2, 2)
plt.scatter(T, E, s=2, c='k')
plt.xlabel(r'$\tau(x) = \|x\|$')
plt.ylabel(r'$|T(x) - \tau(x)|$')
plt.tight_layout()
plt.show()

# Find the first VALID node with greater than 1 degree error in its
# gradient
while True:
    ind = eik.step()
    if angle[ind] > 1:
        print(ind)
        break

# Find noncausal vertices among bad_inds
noncausal = np.empty(bad_inds.shape, dtype=bool)
for i, l in enumerate(bad_inds):
    v = verts[l]
    for l0, l1 in mesh.ve(l):
        v0, v1 = verts[l0], verts[l1]
        if (v0 - v)@(v1 - v) < 0:
            noncausal[i] = True
            break

# Find all noncausal vertices
update_angles_by_vertex = dict()
noncausal = np.empty(nverts, dtype=bool)
for l, v in enumerate(verts):
    update_angles_by_vertex[l] = list()
    for l0, l1 in mesh.ve(l):
        v0, v1 = verts[l0], verts[l1]
        d0 = v0 - v
        d0 /= np.linalg.norm(d0)
        d1 = v1 - v
        d1 /= np.linalg.norm(d1)
        update_angle = np.rad2deg(np.arccos(d0@d1))
        update_angles_by_vertex[l].append(update_angle)
        if update_angle > 90:
            noncausal[l] = True
            break

def get_worst_vertex(update_angles_by_vertex):
    worst_update_angle = 0
    worst_l = None
    for l in range(nverts):
        max_update_angle = max(update_angles_by_vertex[l])
        if max_update_angle > worst_update_angle:
            worst_update_angle = max_update_angle
            worst_l = l
    return worst_l

update_angles = []
for l in range(nverts):
    update_angles.extend(update_angles_by_vertex[l])
update_angles = np.array(update_angles)

plt.figure()
plt.hist(update_angles, bins=129)
plt.xlim(0, 180)
plt.show()


def plot_hists(E, angle):
    plt.figure(figsize=(8, 5))
    plt.subplot(1, 2, 1)
    plt.hist(E[np.isfinite(E)], bins=257)
    xmax = abs(E[np.isfinite(E)]).max()
    plt.xlim(-xmax, xmax)
    plt.yscale('log')
    plt.xlabel(r'$T - \tau$')
    plt.subplot(1, 2, 2)
    plt.hist(angle[np.isfinite(angle)], bins=257)
    plt.xlim(0, 90)
    plt.yscale('log')
    plt.xlabel(r'$\angle (\nabla T, \nabla \tau)$ [Deg.]')
    plt.title('Error Histograms')
    plt.show()

# This algorithm pulls out two arrays, L1 and L2, containing pairs of
# indices for tetrahedron updates running around the edge [l, l0]
l = 1000
l0 = mesh.vv(l)[0]
ec = mesh.ec(l, l0)
L1 = []
L2 = []
for cell in cells[ec]:
    for j, l1 in enumerate(cell):
        if l1 == l or l1 == l0:
            continue
        L1.append(l1)
        break
    for l2 in cell[j+1:]:
        if l2 == l or l2 == l0:
            continue
        L2.append(l2)
        break

# This algorithm sorts L1 and L2 so that the tetrahedron updates are
# sequential, and so that they all oriented the same way; i.e.,
# l2[(i+1)%nec] == l1[i%nec]. This will be important for weeding out
# spurious minimizers, since I'm pretty sure we an have a pair of
# adjacent updates that each have minimizers on the boundary with a
# nonzero Lagrange multiplier. We need to be able to look at
# neighboring minimizers and check if they're the same point---and
# check if the computed jet is the same! If it's not, we might want to
# average... This may improve accuracy!
for i in range(ec.size):
    l1 = L1[i]
    for j in range(i + 1, ec.size):
        if L1[j] == l1:
            L1[j], L2[j] = L2[j], L1[j]
        if L2[j] == l1:
            L1[i + 1], L1[j] = L1[j], L1[i + 1]
            L2[i + 1], L2[j] = L2[j], L2[i + 1]

l = 2812
vc = set()
for l0 in mesh.vv(l):
    for lc in mesh.vc(l0):
        vc.add(lc)
vc = np.array(list(vc))
C = cells[vc]
plotter = pvqt.BackgroundPlotter()
cells_ = np.concatenate([4*np.ones((vc.size, 1)), C], axis=1).astype(np.uint64)
poly_data = pv.PolyData(verts, cells_)
plotter.add_mesh(poly_data, color='white', opacity=0.3, show_edges=True)
plotter.add_mesh(pv.Sphere(0.008, verts[l]), color='blue')

def get_vc(l):
    return np.where((cells == l).any(1))[0]

def get_vv(l):
    return np.setdiff1d(np.unique(cells[get_vc(l)]), l)

################################################################################
# Testing out 4-pt updates

from scipy.optimize import LinearConstraint, minimize

l0 = eik.peek()
l = next(l for l in mesh.vv(l0) if eik.state[l] == jmm.State.Trial.value)
l1, l2 = [
    _ for _ in mesh.vv(l)
    if _ != l0 and eik.state[_] == jmm.State.Valid.value
][:2]

Ls = []
for C in cells[mesh.vc(l)]:
    L = np.array([_ for _ in C if _ != l])
    if (eik.state[L] == jmm.State.Valid.value).all():
        Ls.append(L)
if not Ls:
    print("didn't find any valid update faces!")
Ls = np.array(Ls)

def get_opposite_vertex(l, L):
    return int(np.setdiff1d(cells[
        (cells != l).all(1) &
        (cells == L[0]).any(1) &
        (cells == L[1]).any(1) &
        (cells == L[2]).any(1)
    ], L))

xhat = verts[l]

def do_4point_update(L):
    x = verts[L]
    dx = x[1:] - x[0]

    f = np.array([eik.jet[_][0] for _ in L])
    Df = np.array([(eik.jet[_][1], eik.jet[_][2], eik.jet[_][3]) for _ in L])

    bb3tet = jmm.Bb3Tet(f, Df, x)

    def a(i, j):
        assert i != j
        tmp = np.zeros((4,), dtype=np.float64)
        tmp[i] = 1
        tmp[j] = -1
        return tmp

    def e(i):
        return np.eye(4)[i]

    df_dlam = np.array([bb3tet.Df(e(0), a(i, 0)) for i in [1, 2, 3]])

    def xlam(lam):
        return x[0] + lam@dx

    def L(lam):
        return np.linalg.norm(xlam(lam) - xhat)

    def dL_dlam(lam):
        return dx@xlam(lam)/L(lam)

    def F(lam):
        assert(lam.size == 3)
        assert(all(lam >= 0))
        assert(lam.sum() <= 1)
        b = np.empty((4,))
        b[0] = 1 - lam.sum()
        b[1:] = lam
        return bb3tet.f(b) + L(lam)

    def dF_dlam(lam):
        assert(lam.size == 3)
        assert(all(lam >= 0))
        assert(lam.sum() <= 1)
        b = np.empty((4,))
        b[0] = 1 - lam.sum()
        b[1:] = lam
        return np.array([bb3tet.Df(e(0), a(i, 0)) for i in [1, 2, 3]]) \
            + dL_dlam(lam)

    lb = np.array([0, 0, 0, -np.inf])
    ub = np.array([np.inf, np.inf, np.inf, 1])
    A = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 1]
    ])
    simplex_constraint = LinearConstraint(A, lb, ub, keep_feasible=True)

    lam0 = np.array([1, 1, 1])/4

    res = minimize(F, lam0, method='trust-constr', jac=dF_dlam,
                   constraints=[simplex_constraint], tol=1e-12)

    bopt = np.empty((4,))
    bopt[1:] = res.x
    bopt[0] = 1 - res.x.sum()
    print('optimal b: (%1.4f, %1.4f, %1.4f, %1.4f)' % tuple(bopt))

    Topt = res.fun
    print('optimal T: %1.4f' % Topt)

    DTopt = xhat - bopt@x
    DTopt /= np.linalg.norm(DTopt)
    print('optimal DT: (%1.4f, %1.4f, %1.4f)' % tuple(DTopt))

    Terror = Topt - np.linalg.norm(xhat)
    print('error in T: %1.4f' % Terror)

    DTerror = DTopt - xhat/np.linalg.norm(xhat)
    print('error in DT: (%1.4f, %1.4f, %1.4f)' % tuple(DTerror))

print('correct T: %1.4f' % np.linalg.norm(xhat))
print('correct DT: (%1.4f, %1.4f, %1.4f)' % tuple(xhat/np.linalg.norm(xhat)))
for L in Ls:
    l3 = get_opposite_vertex(l, L)
    L = np.concatenate([L, [l3]]).astype(int)
    print('L = [%d, %d, %d, %d]' % tuple(L))
    do_4point_update(L)


mesh = pv.read('boxfac2/0.00333/box.1.vtk')
V = mesh.points
C = mesh.cells.reshape(mesh.cells.size//5, 5)[:, 1:]
dV = np.concatenate([
    V[C][:, 1, :] - V[C][:, 0, :],
    V[C][:, 2, :] - V[C][:, 0, :],
    V[C][:, 3, :] - V[C][:, 0, :]
], axis=0)
H = np.sqrt(np.sum(dV**2, axis=1))
h = H.mean()
print('avg. edge length: h = %f' % h)

hs = np.array([0.270614, 0.179521, 0.124047, 0.082954, 0.057469, 0.038464, 0.026660])
l2 = np.array([0.000773176, 0.000433334, 0.000339295, 0.000223934, 0.000129865, 5.09184e-05, 2.28107e-05])

dlogh = np.log10(hs[1:]) - np.log10(hs[:-1])
dloge = np.log10(l2[1:]) - np.log10(l2[:-1])

print(dloge/dlogh)

plt.figure(figsize=(3, 3))
plt.loglog(hs, l2, marker='*')
plt.ylabel(r'Error (Relative $\ell_2$)')
plt.xlabel(r'$h$ (Average edge length)')
plt.tight_layout()
plt.show()



# l, l0, l1, l2 = 33, 30, 11, 77
l, l0, l1, l2 = 79, 40, 34, 82
r = 0.05
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(surf_mesh, 'r', 'wireframe')
plotter.add_mesh(pv.Sphere(0.75*r, points[indsrc]), color='white')
if l: plotter.add_mesh(pv.Sphere(r, points[l]), color='red')
if l0: plotter.add_mesh(pv.Sphere(r, points[l0]), color='blue')
if l1: plotter.add_mesh(pv.Sphere(r, points[l1]), color='green')
if l2: plotter.add_mesh(pv.Sphere(r, points[l2]), color='yellow')

l = 81
r = 0.05
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(surf_mesh, 'r', 'wireframe')
plotter.add_mesh(pv.Sphere(r, points[indsrc]), color='white')
plotter.add_mesh(pv.Sphere(r, points[l]), color='blue')
centers, labels = [points[l]], [str(l)]
for l0 in mesh.vv(l):
    plotter.add_mesh(pv.Sphere(0.75*r, points[l0]), color='cyan')
    centers.append(points[l0])
    labels.append(str(l0))
poly_data = pv.PolyData(centers + r*np.ones_like(centers))
poly_data['labels'] = np.array(labels, dtype=str)
plotter.add_point_labels(poly_data, 'labels', point_size=0, font_size=12)

r = 0.05
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(surf_mesh, 'r', 'wireframe')
plotter.add_mesh(pv.Sphere(r, points[indsrc]), color='red', opacity=0.5)
for l in range(points.shape[0]):
    color = 'blue' if not eik.is_diff(l) else 'yellow'
    plotter.add_mesh(pv.Sphere(r/2, points[l]), color=color)
plotter.add_mesh(pv.Sphere(r, points[81]), color='white', opacity=0.35)
plotter.add_mesh(pv.Sphere(r, points[33]), color='pink', opacity=0.35)

plotter = pvqt.BackgroundPlotter()

l0 = 85
while eik.peek() != l0:
    eik.step()

r = 0.05
plotter.clear()
plotter.background_color = (0.3, 0.3, 0.3)
plotter.add_mesh(surf_mesh, 'r', 'wireframe')
plotter.add_mesh(pv.Sphere(r, points[indsrc]), color='red', opacity=0.5)
state_to_color = {
    0: 'red', 1: 'yellow', 2: 'green', 3: None, 4: None, 5: None, 6: 'purple'}
for l in range(points.shape[0]):
    if eik.state[l] == 0:
        continue
    color = state_to_color[eik.state[l]]
    plotter.add_mesh(pv.Sphere(r/2, points[l]), color=color)
plotter.add_mesh(pv.Sphere(r, points[l0]), color='white', opacity=0.5)


L = np.array([l for l, (_, y, _) in enumerate(points)
              if y < 1 and not eik.is_shadow(l)])
print(L[np.argsort([_[0] for _ in eik.jet[L]])])

while True:
    l = eik.step()
    if points[l, 1] < 1 and eik.state[l] == jmm.State.Valid.value:
        print('bad valid point')
        break
    if points[l, 1] >= 1 and eik.state[l] == jmm.State.Shadow.value:
        print('bad shadow point')
        break
print(l)

diff_edges = []
for l0 in range(100):
    for l1 in range(l0):
        if mesh.is_diff_edge(l0, l1):
            diff_edges.append((l0, l1))
r = 0.05
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(surf_mesh, 'r', 'wireframe')
plotter.add_mesh(pv.Sphere(r, points[indsrc]), color='red', opacity=0.75)
for l0, l1 in diff_edges:
    plotter.add_mesh(pv.Sphere(r/2, points[l0]), color='white', opacity=0.5)
    plotter.add_mesh(pv.Sphere(r/2, points[l1]), color='white', opacity=0.5)
    plotter.add_mesh(
        pv.Cylinder((points[l0] + points[l1])/2, points[l1] - points[l0], r/4,
                    np.linalg.norm(points[l1] - points[l0])),
        color='white', opacity=0.5)


def min_tet_alt(i):
    P = points[cells[i]]
    J = np.arange(4)
    h = np.inf
    for j in J:
        I = np.setdiff1d(J, j)
        pj, PI = P[j], P[I]
        n = np.cross(PI[1] - PI[0], PI[2] - PI[0])
        n /= np.linalg.norm(n)
        hj = abs(n@(pj - PI[0]))
        h = min(h, hj)
    return h

min_tet_alts = np.array([min_tet_alt(i) for i in range(cells.shape[0])])
