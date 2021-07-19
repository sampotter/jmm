import meshpy.tet
import numpy as np
import pyvista as pv
import vtk

class UtdWedgeProblem:
    def _get_mesh_plc(self):
        rw, rh = self.w/2, self.h/2

        points = [
            (0, 0, rh), (0, 0, -rh), # (0, 1): center
            (rw, 0, rh), (rw, 0, -rh), # (2, 3): right midpoint
            (rw, -rw, rh), (rw, -rw, -rh), # (4, 5): bottom right corner
            (-rw, -rw, rh), (-rw, -rw, -rh), # (6, 7): bottom left corner
        ]

        facets = [
            (0, 1, 3, 2), # fixed wedge face
            (2, 3, 5, 4), # right face below wedge
            (4, 5, 7, 6), # bottom face
            (0, 2, 4, 6), # ceiling bottom
            (1, 3, 5, 7), # floor bottom
        ]

        if 0 < self.n < 0.25:
            y = rw*np.tan(np.pi*self.n)
            points.extend([
                (-rw, rw, rh), (-rw, rw, -rh), # (8, 9): top left corner
                (rw, rw, rh), (rw, rw, -rh), # (10, 11): top right corner
                (rw, y, rh), (rw, y, -rh), # (12, 13): wedge/boundary intersection
            ])
            facets.extend([
                (6, 7, 9, 8), # left face
                (8, 9, 11, 10), # top face
                (10, 11, 13, 12), # right face above wedge
                (12, 13, 1, 0), # far wedge face
                (6, 8, 0), # ceiling left
                (7, 9, 1), # floor left
                (0, 8, 10, 12), # ceiling top right
                (1, 9, 11, 13), # floor top right
            ])
        elif 0.25 <= self.n < 0.75:
            x = rw/np.tan(np.pi*self.n)
            points.extend([
                (-rw, rw, rh), (-rw, rw, -rh), # (8, 9): top left corner
                (x, rw, rh), (x, rw, -rh), # (10, 11): wedge/boundary intersection
            ])
            facets.extend([
                (6, 7, 9, 8), # left face
                (8, 9, 11, 10), # top face
                (10, 11, 1, 0), # far wedge face
                (6, 8, 0), # ceiling left
                (7, 9, 1), # floor left
                (0, 8, 10), # ceiling top
                (1, 9, 11), # floor top
            ])
        elif 0.75 <= self.n < 1:
            y = -rw*np.tan(np.pi*self.n)
            points.extend([
                (-rw, y, rh), (-rw, y, -rh), # (8, 9): wedge/boundary intersection
            ])
            facets.extend([
                (6, 7, 9, 8), # left face
                (8, 9, 1, 0), # far wedge face
                (6, 8, 0), # ceiling left
                (7, 9, 1), # floor left
            ])
        else:
            raise RuntimeError('got n = %g (should have 0 < n < 1)' % self.n)

        points.append(self.xsrc)

        return points, facets

    def _set_up_mesh_info(self):
        '''Create the instance of `meshpy.MeshInfo`. This involves putting
        together a polygonal mesh which describes the boundary of the
        computational domain, which is a little tricky because of the
        different possible values of `n`.

        '''
        mesh_info = meshpy.tet.MeshInfo()

        points, facets = self._get_mesh_plc()

        mesh_info.set_points(points)
        mesh_info.set_facets(facets)

        return mesh_info

    def __init__(self, maxvol, n, sp, phip, w, h):
        '''Set up a wedge problem for testing the accuracy of `jmm.Eik3` and
        the contents of `jmm.multiple_arrivals` and `jmm.utd`.

        This will set up a wedge with a particular angle, and place a
        point source in its vicinity.

        Since this is an exterior problem, we just mesh a rectangular
        prism containing part of the wedge. The  has dimensions
        `[-w/2, w/2] x [-w/2, w/2] x [-h/2, h/2]`.

        Parameters
        ----------
        maxvol : float
            The maximum volume of a tetrahedron in the mesh.
        n : float
            The UTD wedge angle parameter. Defined so that if theta is
            the dihedral angle of the wedge, then n = theta/pi. For this
            problem, 0 < n < 1 is required.
        sp : float
            The distance from the point source to the diffracting edge.
        phip : float
            The angle between the wedge and the vector pointing to the
            point source from the closest point on the diffracting
            edge.
        w : float
            The width the computational domain bounding box.
        h : float
            The height of the computational domain bounding box.

        '''
        if not 0 < n < 1:
            raise ValueError('n should satisfy 0 < n < 1')

        self.maxvol = maxvol
        self.n = n
        self.sp = sp
        self.phip = phip
        self.w = w
        self.h = h

        self.xsrc = (sp*np.cos(-phip), sp*np.sin(-phip), 0)

        self._mesh_info = self._set_up_mesh_info()
        self._mesh = meshpy.tet.build(self._mesh_info,
                                      max_volume=self.maxvol)

        if self.verts.size == 0:
            raise RuntimeError('something weird happened---mesh is empty!')

        dist_to_xsrc = np.sqrt(np.sum((self.xsrc - self.verts)**2, axis=1))
        if min(dist_to_xsrc) > 0:
            raise RuntimeError("point source wasn't added to mesh")

        self.xsrc_index = np.argmin(dist_to_xsrc)

    @property
    def verts(self):
        return np.array(self._mesh.points, dtype=np.float64)

    @property
    def cells(self):
        return np.array(self._mesh.elements, dtype=np.uintp)

    def plot_mesh(self, plotter, **kwargs):
        '''Plot the tetrahedron mesh discretizing the computational
        domain. Keyword arguments will be forwarded like
        `plotter.add_mesh(..., **kwargs)` so that plotting can be
        controlled.

        '''
        grid = pv.UnstructuredGrid({vtk.VTK_TETRA: self.cells}, self.verts)
        plotter.add_mesh(grid, **kwargs)

if __name__ == '__main__':
    pass
