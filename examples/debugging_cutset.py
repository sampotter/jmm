tet_mesh = 'L/L.1.vtk'
surf_mesh = pv.read('L.obj')

grid = pv.read(tet_mesh)
points = grid.points.copy().astype(np.float64)
cells = grid.cells.reshape(-1, 5)[:, 1:].copy().astype(np.uint64)

plotter = pvqt.BackgroundPlotter()

l0 = 29
while eik.peek() != l0:
    eik.step()
l1 = 60

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

plotter.add_mesh(pv.Sphere(r, points[l0]), color='blue', opacity=0.5)
plotter.add_mesh(pv.Sphere(r, points[l1]), color='white', opacity=0.5)

# for (m0, m1), t in eik.shadow_cutset.items():
#     p0, p1 = points[m0], points[m1]
#     plotter.add_mesh(
#         pv.Cylinder((p0 + p1)/2, p1 - p0, 0.125*r, np.linalg.norm(p1 - p0)),
#         opacity=0.5, color='white')

# def plot_tri(L):
#     plotter.add_mesh(
#         pv.make_tri_mesh(points, np.array(L).reshape(1, 3)),
#         opacity=0.5, color='white')
# plot_tri([20, 5, 85])
