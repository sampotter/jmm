import jmm
import meshio
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt

################################################################################
# Set up some test mesh stuff

# mesh_path = 'room.vtu'
mesh_path = 'sphere.vtu'

print(f'- loading {mesh_path}')
mesh = meshio.read(mesh_path)
verts = mesh.points.astype(np.float64).copy()
cells = mesh.cells_dict['tetra'].astype(np.uintp).copy()
del mesh

mesh = jmm.Mesh3(verts, cells)

def baryc(i, p):
    def V(x, y, z, w):
        return (x - w)@np.cross(y - w, z - w)/6
    a, b, c, d = verts[cells[i]]
    v = V(a, b, c, d)
    return np.array([
        V(p, b, c, d),
        V(a, p, c, d),
        V(a, b, p, d),
        V(a, b, c, p)
    ])/v

i = 0
j = np.random.choice(mesh.cc(i))

Pi = verts[cells[i]]
Pj = verts[cells[j]]

F = np.array([3, 0, 1, 2,
              3, 1, 2, 3,
              3, 2, 3, 0,
              3, 3, 0, 1])

Si = pv.PolyData(Pi, F)
Sj = pv.PolyData(Pj, F)

h = np.sum((Pi[1:] - Pi[0])**2, axis=1).mean()

p = verts[np.setxor1d(cells[i], cells[j])].mean(0)
S = pv.Sphere(radius=h/10, center=p)

plt = pvqt.BackgroundPlotter()
plt.add_mesh(Si, style='wireframe')
plt.add_mesh(Sj, style='wireframe')
plt.add_mesh(S, color='red')

print(baryc(i, p))
print(baryc(j, p))
