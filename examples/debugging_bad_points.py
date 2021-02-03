h = 5*0.00666


def plot_solution(plotter, verts, cells, eik, highlight_ind=None, opacity=1):
    # First, find the cells on the front---we initially take these to
    # be the cells which have exactly three VALID vertices.
    cells_on_front = \
        cells[(eik.state[cells] == jmm.State.Valid.value).sum(1) == 3]
    # Next, we want to filter out any cells on the front that contain
    # the same triangle. To do this, we count the corresponding
    # triangles on the front using a dictionary.
    tris_on_front = dict()
    for cell, state in zip(cells_on_front, eik.state[cells_on_front]):
        sorted_tri_inds = sorted(
            i for i, s in zip(cell, state) if s == jmm.State.Valid.value)
        key = tuple(sorted_tri_inds)
        if key not in tris_on_front:
            tris_on_front[key] = 1
        else:
            tris_on_front[key] += 1
    tris_on_front = np.array(
        [_ for _, count in tris_on_front.items() if count == 1],
        dtype=cells.dtype)
    # Now, find the unique indices so that we can look things up
    uniq_inds = np.unique(tris_on_front.ravel()).tolist()
    for i in range(tris_on_front.shape[0]):
        for j in range(3):
            tris_on_front[i, j] = uniq_inds.index(tris_on_front[i, j])
    # Pull out the jets and vertices corresponding to the unique
    # indices
    T = np.array([J[0] for J in eik.jet[uniq_inds]])
    DT = np.array([(J[1], J[2], J[3]) for J in eik.jet[uniq_inds]])
    # Put together the data for a PyVista PolyData instance containing
    # a triangle mesh representing the front
    points = verts[uniq_inds]
    faces = np.concatenate([
        3*np.ones((tris_on_front.shape[0], 1), dtype=tris_on_front.dtype),
        tris_on_front
    ], axis=1)
    poly_data = pv.PolyData(points, faces)
    poly_data.point_arrays['T'] = T
    # Add it to the plotter and plot the values of the eikonal
    plotter.add_mesh(poly_data, scalars='T', cmap=cc.cm.rainbow,
                     opacity=opacity, show_edges=True)
    # Now, traverse each point and gradient, and add a colored arrow
    if highlight_ind is None:
        for ind, p, d in zip(uniq_inds, points, DT):
            Dtau = p/np.linalg.norm(p)
            alpha = 256*(np.dot(d, Dtau) + 1)/2
            color = cc.cm.coolwarm_r(alpha)[:3]
            sphere = pv.Sphere(h/12, p)
            plotter.add_mesh(sphere, color=color)
            arrow = pv.Arrow(p, d, tip_length=0.1, scale=h)
            plotter.add_mesh(arrow, color=color)
    else:
        for ind, p, d in zip(uniq_inds, points, DT):
            if ind == highlight_ind:
                color = 'blue'
            else:
                color = 'white'
            sphere = pv.Sphere(h/12, p)
            plotter.add_mesh(sphere, color=color)
            arrow = pv.Arrow(p, d, tip_length=0.1, scale=h)
            plotter.add_mesh(arrow, color=color)


def plot_edge(plotter, l0, l1, opacity=0.6):
    plotter.add_mesh(pv.Sphere(0.0125, verts[l0]), color='white', opacity=opacity)
    plotter.add_mesh(pv.Sphere(0.0125, verts[l1]), color='white', opacity=opacity)
    plotter.add_mesh(pv.Cylinder((verts[l0] + verts[l1])/2, verts[l1] - verts[l0],
                                 0.0025, np.linalg.norm(verts[l1] - verts[l0])),
                     color='white', opacity=opacity)


while eik.peek() != 2198:
    eik.step()
eik.step()
l = 2627
l0 = 2198


while eik.peek() != 769:
    eik.step()
eik.step()
l = 2627
l0 = 769


while eik.peek() != 1896:
    eik.step()
eik.step()
l = 2627
l0 = 1896
l1 = 769

l2 = 1927

alpha1 = np.array([0, 0, -0.00023929887561292829])
lam1 = np.array([0.50476302928905725, 0])
T1 = 1.3891790361956593
DT1 = np.array([0.32450279153558675, -0.6138109399381555, -0.71967636358147102])

l2 = 2198
alpha2 = np.array([0, 0, -0.000034241291624229064])
lam2 = np.array([0.50476302928905725, 0])
T2 = 1.3891790361956593
DT2 = np.array([0.32450279153558675, -0.6138109399381555, -0.71967636358147102])


plotter = pvqt.BackgroundPlotter()
plot_solution(plotter, verts, cells, eik, highlight_ind=l, opacity=1)
plotter.add_mesh(pv.Sphere(0.01, np.zeros(3)), color='red')
plotter.add_mesh(pv.Sphere(0.01, verts[l]), color='red')
plotter.add_mesh(pv.Cylinder(verts[l]/2, verts[l], 0.0025, np.linalg.norm(verts[l])), color='red')

lvalid = [_ for _ in mesh.vv(l) if eik.state[_] == jmm.State.Valid.value]
for _ in lvalid:
    plotter.add_mesh(pv.Sphere(0.01, verts[_]), color='blue')
ltrial = np.setdiff1d(mesh.vv(l), lvalid)
for _ in ltrial:
    plotter.add_mesh(pv.Sphere(0.01, verts[_]), color='orange')

plot_edge(plotter, l0, l1)
plot_edge(plotter, l1, l2)
plot_edge(plotter, l2, l0)
