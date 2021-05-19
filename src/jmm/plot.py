import numpy as np
import pyvista as pv
import vtk

from scipy.optimize import brentq

import jmm.defs

def plot_point(plotter, point, radius, **kwargs):
    '''Uses `plotter` to plot `mesh`.'''
    plotter.add_mesh(pv.Sphere(radius, point), **kwargs)

def plot_arrow(plotter, origin, direction, length, **kwargs):
    plotter.add_mesh(pv.Arrow(origin, direction, scale=length), **kwargs)

def plot_mesh2(plotter, mesh, **kwargs):
    '''Uses `plotter` to plot `mesh`. The keyword arguments stored in
    `kwargs` are forwarded to the call `plotter.add_mesh`.

    '''
    grid = pv.UnstructuredGrid({vtk.VTK_TRIANGLE: mesh.faces}, mesh.verts)
    plotter.add_mesh(grid, **kwargs)

def plot_eik3(plotter, eik, **kwargs):

    '''Plot the wavefront of an instance of `jmm.eik.Eik3`. The keyword
    arguments passed in `kwargs` will be forwarded to the call to
    `pyvista.Plotter.add_mesh` which plots the triangulated wavefront.

    '''

    verts = eik.mesh.verts
    cells = eik.mesh.cells

    # First, find the cells on the front---we initially take these to
    # be the cells which have exactly three VALID vertices.
    cells_on_front = \
        cells[(eik.state[cells] == jmm.defs.State.Valid.value).sum(1) == 3]

    # Next, we want to filter out any cells on the front that contain
    # the same triangle. To do this, we count the corresponding
    # triangles on the front using a dictionary.
    tris_on_front = dict()
    for cell, state in zip(cells_on_front, eik.state[cells_on_front]):
        sorted_tri_inds = sorted(
            i for i, s in zip(cell, state) if s == jmm.defs.State.Valid.value)
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
    if tris_on_front.size > 0:
        points = verts[uniq_inds]
        faces = np.concatenate([
            3*np.ones((tris_on_front.shape[0], 1), dtype=tris_on_front.dtype),
            tris_on_front
        ], axis=1)
        poly_data = pv.PolyData(points, faces)
        poly_data.point_arrays['T'] = T

        # Add it to the plotter and plot the values of the eikonal
        plotter.add_mesh(poly_data, scalars='T', **kwargs)

def plot_shadow(plotter, field, **kwargs):
    '''Contour `field`'s shadow boundary and plot it using `plotter`,
    forwarding `**kwargs` to the call to `plotter.add_mesh`.

    '''

    mesh = field.domain.mesh

    cell_masks = field.shadow_mask[mesh.cells]
    num_shadow = cell_masks.sum(1)
    bracket_mask = (0 < num_shadow) & (num_shadow < 4)
    bracket_inds = np.where(bracket_mask)[0]
    bracket = mesh.cells[bracket_mask]

    zverts, zcells = [], []

    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    e = np.eye(4)

    for lc, cell, mask in zip(bracket_inds, bracket, cell_masks[bracket_inds]):
        X = mesh.verts[cell]

        bb = field.eik.get_bezier_tetra(lc)
        extended_bb = field.extended_eik.get_bezier_tetra(lc)
        z = lambda b: bb.f(b) - extended_bb.f(b) - field._shadow_mask_threshold

        cut_edges = [(i, j) for i, j in edges if mask[i] != mask[j]]
        assert len(cut_edges) in {3, 4}

        if len(cut_edges) == 3:
            zcell = []
            for i, j in cut_edges:
                b = lambda t: (1 - t)*e[i] + t*e[j]
                t = brentq(lambda t: z(b(t)), 0, 1)
                xt = b(t)@X
                zcell.append(len(zverts))
                zverts.append(xt)
            zcells.append(zcell)

        if len(cut_edges) == 4:
            offset = len(zverts)

            ################################################################
            # First, find all the vertices

            for i, j in cut_edges:
                b = lambda t: (1 - t)*e[i] + t*e[j]
                t = brentq(lambda t: z(b(t)), 0, 1)
                xt = b(t)@X
                zverts.append(xt)

            (i0, i1), (j0, j1) = set(edges) - set(cut_edges)

            ei = (e[i0] + e[i1])/2
            ej = (e[j0] + e[j1])/2
            b = lambda t: (1 - t)*ei + t*ej

            # try to find the vertex in the middle---this may fail
            try:
                t = brentq(lambda t: z(b(t)), 0, 1)
            except:
                print('found weird cell')
                continue
            xt = b(t)@X
            zverts.append(xt)

            ################################################################
            # Triangulate the five vertices

            i = 0 # pick an arbitrary starting index
            I = {i} # visited vertices

            # This next loop is truly horrible! Can definitely clean this up
            for j in range(4):

                # find the indices of the cutedges containing i, and use them
                # to index the edge of a triangle which will be on the
                # boundary of the tetrahedron
                zcell = [_ for _, e in enumerate(cut_edges) if i in e] + [4]
                zcells.append(offset + np.array(zcell))
                if j == 3:
                    break

                # find the next unvisited vertex that has the opposite state of i
                I_next = [_ for _ in cut_edges[zcell[0]] if _ not in I]
                if not I_next:
                    I_next = [_ for _ in cut_edges[zcell[1]] if _ not in I]
                assert mask[i] != mask[I_next[0]]
                i = I_next[0]
                I.add(i)

    zverts, zcells = np.array(zverts), np.array(zcells)

    grid = pv.UnstructuredGrid({vtk.VTK_TRIANGLE: zcells}, zverts)

    plotter.add_mesh(grid, **kwargs)
