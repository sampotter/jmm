import colorcet as cc
import jmm
import numpy as np
import pyvista as pv

def plot_wavefront(plotter, eik):
    verts = eik.mesh.verts
    cells = eik.mesh.cells

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
    if tris_on_front.size > 0:
        points = verts[uniq_inds]
        faces = np.concatenate([
            3*np.ones((tris_on_front.shape[0], 1), dtype=tris_on_front.dtype),
            tris_on_front
        ], axis=1)
        poly_data = pv.PolyData(points, faces)
        poly_data.point_arrays['T'] = T

        # Add it to the plotter and plot the values of the eikonal
        plotter.add_mesh(poly_data, scalars='T', cmap=cc.cm.rainbow,
                         opacity=1.0, show_edges=True)
