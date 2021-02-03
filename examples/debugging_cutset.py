tet_mesh = 'L/L.1.vtk'
surf_mesh = pv.read('L.obj')

grid = pv.read(tet_mesh)
points = grid.points.copy().astype(np.float64)
cells = grid.cells.reshape(-1, 5)[:, 1:].copy().astype(np.uint64)

plotter = pvqt.BackgroundPlotter()

# l0 = 71 # 14
# while eik.peek() != l0:
#     eik.step()
eik.solve()

def plot_point(l, scale=1, color='white', opacity=1):
    plotter.add_mesh(pv.Sphere(scale*r, points[l]),
                     color=color, opacity=opacity)

r = 0.05
plotter.clear()
plotter.background_color = (0.3, 0.3, 0.3)
plotter.add_mesh(surf_mesh, 'r', 'wireframe')

plot_point(indsrc, color='red', opacity=0.5)

state_to_color = {
    0: 'red', 1: 'yellow', 2: 'green', 3: None, 4: None, 5: None, 6: 'purple'}
for l in range(points.shape[0]):
    if eik.state[l] == 0:
        continue
    color = state_to_color[eik.state[l]]
    plot_point(l, scale=1/2, color=color)

h = 0.1

def plot_jet(l):
    x = points[l]
    jet = eik.jet[l]
    DT = np.array([jet[1], jet[2], jet[3]])
    plotter.add_mesh(pv.Arrow(x, DT, scale=h), color='white', opacity=0.6)

# plot_point(l0, opacity=0.5, scale=1.1)

for l, state in enumerate(eik.state):
    if state == jmm.State.Valid.value:
        plot_jet(l)

def plot_cutset(plot_grad=True):
    for (m0, m1), cutedge in eik.shadow_cutset.items():
        t = cutedge.t
        p0, p1 = points[m0], points[m1]
        plotter.add_mesh(
            pv.Cylinder((p0 + p1)/2, p1 - p0, 0.125*r, np.linalg.norm(p1 - p0)),
            opacity=0.5, color='white')
        pt = (1 - t)*p0 + t*p1
        plotter.add_mesh(pv.Sphere(r/3, pt), color='white', opacity=0.6)
        if plot_grad:
            plotter.add_mesh(pv.Arrow(pt, cutedge.n, scale=h),
                             color='white', opacity=0.6)
plot_cutset(False)

def plot_update(l0):
    par = eik.get_parent(l0)
    if par.size == 1:
        raise Exception('blah')
    elif par.size == 2:
        p0, p1 = points[par.l[0]], points[par.l[1]]
        plotter.add_mesh(
            pv.Cylinder((p0 + p1)/2, p1 - p0, 0.1*r, np.linalg.norm(p1 - p0)),
            opacity=0.5, color='white')
    else:
        plotter.add_mesh(
            pv.make_tri_mesh(points, np.array(L).reshape(1, 3)),
            opacity=0.5, color='white')
    plotter.add_mesh(
        pv.Sphere(r/3, eik.get_parent(l0).b@points[eik.get_parent(l0).l]),
        color='red', opacity=0.6)
# plot_update(l0)
