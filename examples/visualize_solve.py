'''This is a version of basic_solve.c written in Python, intended to
demonstrate the Python API's use.'''


import colorcet as cc
import jmm
import numpy as np
import os
import pyvista as pv
import pyvistaqt as pvqt
import sys


# Factoring radius
R0 = 0.1


h = 0.025


def plot_solution(plotter, verts, cells, eik, highlight_ind=None):
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
    plotter.add_mesh(poly_data, scalars='T', cmap=cc.cm.rainbow)
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


if __name__ == '__main__':
    print(sys.argv)

    if len(sys.argv) == 1:
        print('- No command-line arguments: using defaults')
        verts_bin_path = 'sphere_verts.bin'
        cells_bin_path = 'sphere_cells.bin'
        indsrc = 830
        jets_bin_path = 'sphere_jets.bin'
    elif len(sys.argv) != 5:
        print('usage: %s <verts.bin> <cells.bin> <indsrc> <jets.bin>')
        exit(1)
    else:
        verts_bin_path = sys.argv[1]
        cells_bin_path = sys.argv[2]
        indsrc = int(sys.argv[3])
        jets_bin_path = sys.argv[4]

    print(f'- Reading vertices from {verts_bin_path}')
    print(f'- Reading cells from {cells_bin_path}')
    print(f'- Making vertex {indsrc} a point source')
    print(f'- Writing jets to {jets_bin_path}')

    verts = np.fromfile(verts_bin_path, np.float64)
    verts = verts.reshape(verts.size//3, 3)

    cells = np.fromfile(cells_bin_path, np.uintp)
    cells = cells.reshape(cells.size//4, 4)

    mesh = jmm.Mesh3(verts, cells)

    eik = jmm.Eik3(mesh)

    nverts = verts.shape[0]
    xsrc = verts[indsrc]

    valid_inds = [indsrc]

    eik.add_valid(indsrc, jmm.Jet3(0, np.nan, np.nan, np.nan))

    # Start by making all neighbors of the point source VALID.
    for i in mesh.vv(indsrc):
        dx = verts[i] - xsrc
        T = np.linalg.norm(dx)
        DT = dx/T
        eik.add_valid(i, jmm.Jet3(T, *DT))
        if i not in valid_inds:
            valid_inds.append(i)

    # Make all FAR neighbors of VALID nodes TRIAL.
    for i in valid_inds:
        for j in mesh.vv(i):
            if eik.is_far(j):
                dx = verts[j] - xsrc
                T = np.linalg.norm(dx)
                DT = dx/T
                eik.add_trial(j, jmm.Jet3(T, *DT))

    # Call plot_next to step once and plot the VALID front, highlight
    # the new VALID node.
    def plot_next():
        old_valid = np.where(eik.state == 2)[0]
        eik.step()
        new_valid_ind = np.setdiff1d(
            np.where(eik.state == 2)[0], old_valid)[0]
        plotter = pvqt.BackgroundPlotter()
        plot_solution(
            plotter, verts, cells, eik, highlight_ind=new_valid_ind)

    # plot_next()

    # L = []
    # l0 = eik.peek()
    # while l0 != 987:
    #     print(l0)
    #     L.append(l0)
    #     assert(eik.step() == l0)
    #     l0 = eik.peek()
