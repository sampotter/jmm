import colorcet as cc
import pyvista as pv
import pyvistaqt as pvqt
import matplotlib.pyplot as plt
import numpy as np
import vtk

from pathlib import Path

plt.ion()

itd_build_dir = Path('../../builddir/examples/itd')

grid = pv.read(itd_build_dir/'HUTUB_pp2_in_cube_50k.off')

V = grid.points
F = grid.cells.reshape(-1, 4)[:, 1:]

V0 = V[F][:, 0, :]
V1 = V[F][:, 1, :]
V2 = V[F][:, 2, :]

T1 = V1 - V0
T2 = V2 - V0

A = np.sqrt(np.sum(np.cross(T1, T2)**2, axis=1))/2

V_tet = np.fromfile(itd_build_dir/'verts.bin').reshape(-1, 3)
C_tet = np.fromfile(itd_build_dir/'cells.bin', dtype=np.uint64).reshape(-1, 4)

grid_tet = pv.UnstructuredGrid({vtk.VTK_TETRA: C_tet}, V_tet)

P_tet = grid_tet.cell_centers().points

cell_ind = (P_tet[:, 1] < 0).nonzero()[0]
# subgrid = grid_tet.extract_cells(cell_ind)

bbox = pv.Box(grid_tet.bounds)

# plotter = pvqt.BackgroundPlotter()
# # plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
# plotter.add_mesh(grid_tet, style='wireframe')
# plotter.add_mesh(pv.Sphere(4, (V_tet[-1, 0], V_tet[-1, 1], V_tet[-1, 2])), color='red', opacity=0.5)
# # plotter.add_mesh(bbox, 'r', 'wireframe')
# # plotter.add_legend([[' Input Mesh ', 'r'],
# #                     [' Tessellated Mesh ', 'black']])
# plotter.show()

jet_L = np.fromfile(itd_build_dir/'jet_L.bin').reshape(-1, 4)
Tau_L = jet_L[:, 0]
GradTau_L = jet_L[:, 1:]

jet_R = np.fromfile(itd_build_dir/'jet_R.bin').reshape(-1, 4)
Tau_R = jet_R[:, 0]
GradTau_R = jet_R[:, 1:]

c = 340.3e3

T_L = Tau_L/c
T_R = Tau_R/c
itd = T_R - T_L
itd_abs_max = abs(itd).max()
itd_contours = np.linspace(-0.5e-3, 0.5e-3, 11)

grid_tet['T_L'] = T_L
grid_tet['T_R'] = T_R
grid_tet['itd'] = itd

bust = pv.read('HUTUB_pp2_in_cube_50k_origin_fixed_bust_only_BROKEN.off')

clipped = grid_tet.clip('z', origin=(0, 0, 0))
contours = clipped.contour(scalars='T_L')
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(bbox, 'r', 'wireframe')
plotter.add_mesh(contours)
plotter.add_mesh(bust)

clipped = grid_tet.clip('z', origin=(0, 0, 0))
contours = clipped.contour(scalars='T_R')
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(bbox, 'r', 'wireframe')
plotter.add_mesh(contours)
plotter.add_mesh(bust)

clipped = grid_tet.clip('z', origin=(0, 0, 0))
contours = clipped.contour(itd_contours, scalars='itd')
plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(bbox, 'r', 'wireframe')
plotter.add_mesh(contours)
plotter.add_mesh(bust)

El, Az = [], []
with open('fliege64.txt', 'r') as f:
    for line in f:
        el, az = line.strip().split()
        El.append(float(el))
        Az.append(float(az))
El, Az = np.array(El), np.array(Az)

r = 500
x = r*np.cos(Az)*np.sin(El)
y = r*np.sin(Az)*np.sin(El)
z = r*np.cos(El)

X_fliege = np.fromfile(itd_build_dir/'fliege64.bin').reshape(-1, 3)
grid_fliege = pv.PolyData(X_fliege)

itd_fliege = np.fromfile(itd_build_dir/'itd.bin')
grid_fliege['itd'] = itd_fliege

plotter = pvqt.BackgroundPlotter()
plotter.add_mesh(bust)
plotter.add_mesh(grid_fliege)

plt.figure(figsize=(8, 4))
plt.scatter(Az, El, s=15, c=itd_fliege, cmap=cc.cm.gouldian)
plt.colorbar()
plt.xlim(-np.pi, np.pi)
plt.ylim(0, np.pi)
plt.title('ITD [s]')
plt.xlabel('Azimuth')
plt.ylabel('Elevation')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig('itd.png')

m_uniform = 256
n_uniform = 512

itd_uniform = np.fromfile(itd_build_dir/'itd_uniform.bin').reshape(
    m_uniform, n_uniform)

itd_uniform_abs_max = np.nanmax(abs(itd_uniform))
itd_uniform_levels = np.linspace(-itd_uniform_abs_max, itd_uniform_abs_max, 21)

Az_uniform, El_uniform = np.meshgrid(
    np.linspace(0, 2*np.pi, n_uniform),
    np.linspace(0, np.pi, m_uniform),
)

plt.figure()
# plt.imshow(itd_uniform, extent=[0, 2*np.pi, 0, np.pi], cmap=cc.cm.blues)
plt.contourf(Az_uniform, El_uniform, itd_uniform, itd_uniform_levels,
             cmap=cc.cm.gouldian)
plt.xticks(np.linspace(0, 2*np.pi, 9, endpoint=True))
plt.yticks(np.linspace(0, np.pi, 5, endpoint=True))
plt.gca().invert_yaxis()
plt.gca().set_yticklabels(['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
plt.gca().set_xticklabels(['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$',
                           r'$\pi$', r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$',
                           r'$2\pi$'])
plt.colorbar()
plt.gca().set_aspect('equal')
plt.show()
