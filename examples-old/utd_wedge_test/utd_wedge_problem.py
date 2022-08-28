import colorcet as cc
import logging
import meshpy.tet
import numpy as np
import pyvista as pv
import scipy.optimize
import vtk

from cached_property import cached_property

from meshplex import MeshTetra

from jmm.multiple_arrivals import DiffractedField, Domain, PointSourceField, \
    ReflectedField

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

    def _get_refl_indices(self):
        near_index, far_index = None, None
        for i, faces in enumerate(self.domain.mesh.reflectors):
            x, y, z = self.verts[np.unique(faces)].T
            if (y == 0).all():
                near_index = i
            elif ((x == 0) & (y == 0) & (-self.h/2 < z) & (z < self.h/2)).any():
                far_index = i
        if near_index is None or far_index is None:
            raise RuntimeError("couldn't figure out wedge's reflector indices")
        return near_index, far_index

    def __init__(self, maxvol, n, sp, phip, w, h, R, r, omega,
                 min_dB=-60):
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
        R : float
            The reflection coefficient.
        r : float
            The factoring radius for the point source.
        omega : float
            The frequency of the wave.
        min_dB : float
            The minimum decibel value, below which clamping occurs (used for
            plots).

        '''
        self.log = logging.getLogger(type(self).__name__)

        if not 0 < n < 1:
            raise ValueError('n should satisfy 0 < n < 1')

        self.maxvol = maxvol
        self.n = n
        self.sp = sp
        self.phip = phip
        self.w = w
        self.h = h
        self.R = R
        self.r = r
        self.omega = omega
        self.min_dB = min_dB

        self.xsrc = (sp*np.cos(-phip), sp*np.sin(-phip), 0)

        self._mesh_info = self._set_up_mesh_info()
        self._mesh = meshpy.tet.build(self._mesh_info,
                                      max_volume=self.maxvol)

        if self.verts.size == 0:
            raise RuntimeError('something weird happened---mesh is empty!')

        if min(self.dist2xsrc) > 0:
            raise RuntimeError("point source wasn't added to mesh")

        self.xsrc_index = np.argmin(self.dist2xsrc)

        self.domain = Domain(self.verts, self.cells, self.R)

        self.edge_lengths = MeshTetra(self.verts, self.cells).edge_lengths
        self.max_edge_length = self.edge_lengths.max()
        self.avg_edge_length = np.mean(self.edge_lengths)

        self.pt_src_field = PointSourceField(
            self.domain, self.xsrc_index, self.omega, self.r)
        self.pt_src_field.solve()
        self.log.info('computed field for point source')

        near_refl_index, far_refl_index = self._get_refl_indices()
        self.near_refl_index = near_refl_index
        self.far_refl_index = far_refl_index

        self.near_refl_field = None
        if self.should_do_near_reflection:
            self.near_refl_faces = \
                self.domain.mesh.get_reflector(self.near_refl_index)
            self.near_refl_BCs = \
                self.pt_src_field._get_reflection_BCs(self.near_refl_faces)
            if not self.near_refl_BCs:
                raise RuntimeError("couldn't generate BCs for near reflection")
            self.near_refl_field = ReflectedField(
                self.domain, self.near_refl_index, self.omega,
                *self.near_refl_BCs, parent=self.pt_src_field)
            self.near_refl_field.solve()
            self.log.info('computed near reflection')

        self.far_refl_field = None
        if self.should_do_far_reflection:
            self.far_refl_faces = \
                self.domain.mesh.get_reflector(self.far_refl_index)
            self.far_refl_BCs = \
                self.pt_src_field._get_reflection_BCs(self.far_refl_faces)
            if not self.far_refl_BCs:
                raise RuntimeError("couldn't generate BCs for far reflection")
            self.far_refl_field = ReflectedField(
                self.domain, self.far_refl_index, self.omega,
                *self.far_refl_BCs, parent=self.pt_src_field)
            self.far_refl_field.solve()
            self.log.info('computed far reflection')

        diffractors = list(self.domain.mesh.diffractors)
        if len(diffractors) != 1:
            raise RuntimeError('domain contains more than one diffractor')

        self.diff_edges = diffractors[0]
        self.diff_BCs = self.pt_src_field._get_diffraction_BCs(self.diff_edges)
        if not self.diff_BCs:
            raise RuntimeError("couldn't generate BCs for diffraction")
        self.diff_field = DiffractedField(
            self.domain, 0, self.omega, *self.diff_BCs,
            parent=self.pt_src_field)
        self.diff_field.solve()
        self.log.info('computed diffracted field')

    @property
    def should_do_near_reflection(self):
        return self.phip < np.pi

    @property
    def should_do_far_reflection(self):
        return (2 - self.n)*np.pi - self.phip < np.pi

    @cached_property
    def verts(self):
        return np.array(self._mesh.points, dtype=np.float64)

    @cached_property
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

    def compute_RMSE(self, errors):
        return np.sqrt(np.sum(errors**2)/self.verts.shape[0])

    def compute_MAE(self, errors):
        return abs(errors).sum()/self.verts.shape[0]

    ########################################################################
    # Groundtruth geometric things

    @cached_property
    def dist2xsrc(self):
        return np.sqrt(np.sum((self.xsrc - self.verts)**2, axis=1))

    @cached_property
    def dist2xsrc_plus_dist2edge(self):
        e = np.array([0, 0, 1])
        def get_dist(x):
            if x[0] == 0 and x[1] == 0:
                return np.linalg.norm(x - self.xsrc)
            def d(z):
                return np.add(
                    np.dot(e, z*e - self.xsrc)/np.linalg.norm(z*e - self.xsrc),
                    np.dot(e, z*e - x)/np.linalg.norm(z*e - x)
                )
            z = scipy.optimize.brentq(d, -self.h/2, self.h/2)
            return np.linalg.norm(z*e - self.xsrc) + np.linalg.norm(z*e - x)
        return np.array([get_dist(x) for x in self.verts])

    @cached_property
    def grad_of_dist2xsrc_plus_dist2edge(self):
        e = np.array([0, 0, 1])
        def get_grad(x):
            def d(z):
                return np.add(
                    np.dot(e, z*e - self.xsrc)/np.linalg.norm(z*e - self.xsrc),
                    np.dot(e, z*e - x)/np.linalg.norm(z*e - x)
                )
            z = scipy.optimize.brentq(d, -self.h/2, self.h/2)
            return (x - z*e)/np.linalg.norm(x - z*e)
        return np.array([get_grad(x) for x in self.verts])

    @cached_property
    def angle(self):
        arctan2 = np.arctan2(self.verts[:, 1], self.verts[:, 0])
        return np.mod(2*np.pi - arctan2, 2*np.pi)

    def reflect_over_near_reflector(self, x):
        if not self.should_do_near_reflection:
            raise RuntimeError('not doing near reflection')
        x0 = self.verts[self.near_refl_faces[0, 0]]
        nu = self.domain.mesh.get_face_normal(*self.near_refl_faces[0])
        refl = np.eye(3) - 2*np.outer(nu, nu)
        return refl@(x - x0) + x0

    @cached_property
    def near_refl_ximg(self):
        return self.reflect_over_near_reflector(self.xsrc)

    @cached_property
    def near_refl_dist2img(self):
        if not self.should_do_near_reflection:
            raise RuntimeError('not doing near reflection')
        return np.sqrt(np.sum((self.near_refl_ximg - self.verts)**2, axis=1))

    @cached_property
    def near_refl_grad_of_dist2img(self):
        if not self.should_do_near_reflection:
            raise RuntimeError('not doing near reflection')
        xsrc_refl = self.reflect_over_near_reflector(self.xsrc)
        return (self.verts - xsrc_refl)/self.near_refl_dist2img.reshape(-1, 1)

    def reflect_over_far_reflector(self, x):
        if not self.should_do_far_reflection:
            raise RuntimeError('not doing far reflection')
        x0 = self.verts[self.far_refl_faces[0, 0]]
        nu = self.domain.mesh.get_face_normal(*self.far_refl_faces[0])
        refl = np.eye(3) - 2*np.outer(nu, nu)
        return refl@(x - x0) + x0

    @cached_property
    def far_refl_ximg(self):
        return self.reflect_over_far_reflector(self.xsrc)

    @cached_property
    def far_refl_dist2img(self):
        if not self.should_do_far_reflection:
            raise RuntimeError('not doing far reflection')
        xsrc_refl = self.reflect_over_far_reflector(self.xsrc)
        return np.sqrt(np.sum((self.far_refl_ximg - self.verts)**2, axis=1))

    @cached_property
    def far_refl_grad_of_dist2img(self):
        if not self.should_do_far_reflection:
            raise RuntimeError('not doing far reflection')
        xsrc_refl = self.reflect_over_far_reflector(self.xsrc)
        return (self.verts - xsrc_refl)/self.far_refl_dist2img.reshape(-1, 1)

    ########################################################################
    # Masks for each field

    @cached_property
    def pt_src_mask(self):
        return self.angle < self.phip + np.pi

    @cached_property
    def near_refl_mask(self):
        if not self.should_do_near_reflection:
            raise RuntimeError('not doing near reflection')
        return self.angle < np.pi - self.phip

    @cached_property
    def far_refl_mask(self):
        if not self.should_do_far_reflection:
            raise RuntimeError('not doing far reflection')
        # TODO: can probably simplify the next line, but the meaning of it is:
        # - with the angle increasing from 0 as we rotate from the
        #   near face to the far face...
        # - ... find the positive angle from the far reflector extended to other
        #   half of the domain to point source location... (this is the
        #   "phip - ((2 - n)*np.pi - np.pi)" term)
        # - ... and then subtract that angle from the angle of the far reflector
        #   to get the angle of the reflection zone boundary
        far_refl_angle = (2 - self.n)*np.pi
        return np.logical_or(
            (self.verts[:, 0] == 0) & (self.verts[:, 1] == 0),
            self.angle > far_refl_angle - (self.phip - (far_refl_angle - np.pi))
        )

    ########################################################################
    # Plotting utilities

    @property
    def padding_scale(self):
        return 1.1

    @cached_property
    def xlim(self):
        return -self.padding_scale*self.w/2, self.padding_scale*self.w/2

    @cached_property
    def ylim(self):
        return -self.padding_scale*self.w/2, self.padding_scale*self.w/2

    @cached_property
    def dB_cmap(self):
        return cc.cm.fire

    def make_scatter_plot(self, fig, ax, values, clim, title, cmap, unit_str,
                          plot_pt_src=False,
                          plot_near_img=False,
                          plot_far_img=False,
                          plot_diff_src=False,
                          plot_near_refl_bd=False,
                          plot_far_refl_bd=False):
        vmin, vmax = clim
        ax.set_title(title)
        zorder = 1
        if plot_near_refl_bd:
            zorder += 1
            r = self.w*np.sqrt(2)
            th = np.pi - self.phip
            x, y = r*np.cos(th), -r*np.sin(th)
            ax.plot([0, x], [0, y], linewidth=2, linestyle='--', c='white')
        if plot_far_refl_bd:
            zorder += 1
            r = self.w*np.sqrt(2)
            dth = (2*np.pi - self.phip) - (np.pi*self.n + np.pi/2)
            th = self.phip + 2*dth
            x, y = r*np.cos(th), r*np.sin(th)
            ax.plot([0, x], [0, y], linewidth=2, linestyle='--', c='white')
        zorder += 1
        sc = ax.scatter(*self.verts[:, :2].T, s=5, c=values,
                        vmin=vmin, vmax=vmax, cmap=cmap, zorder=zorder)
        if plot_pt_src:
            zorder += 1
            ax.scatter(*self.xsrc[:2], s=15, facecolor='red', edgecolor='pink',
                       zorder=zorder)
        if plot_near_img:
            zorder += 1
            ax.scatter(*self.near_refl_ximg[:2], s=15, facecolor='orange',
                       edgecolor='yellow', zorder=zorder)
        if plot_far_img:
            zorder += 1
            ax.scatter(*self.far_refl_ximg[:2], s=15,
                       edgecolor='green', facecolor='cyan',
                       zorder=zorder)
        if plot_diff_src:
            zorder += 1
            ax.scatter([0], [0], s=15, edgecolor='purple', facecolor='magenta',
                       zorder=zorder)
        ax.set_aspect('equal')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_xlim(*self.xlim)
        ax.set_ylim(*self.ylim)
        fig.colorbar(sc, ax=ax)

    ########################################################################
    # Point source fields

    # eikonal

    @cached_property
    def pt_src_T(self):
        return self.pt_src_field.eik.T

    def plot_pt_src_T(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.pt_src_T,
            clim=(np.nanmin(self.pt_src_T), np.nanmax(self.pt_src_T)),
            title=r'$T$ [m] (point source)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_pt_src=True,
        )

    @cached_property
    def pt_src_tau(self):
        mask = self.pt_src_mask
        tau = np.empty(self.verts.shape[0], dtype=self.verts.dtype)
        tau[mask] = self.dist2xsrc[mask]
        tau[~mask] = self.dist2xsrc_plus_dist2edge[~mask]
        return tau

    @cached_property
    def error_pt_src_tau(self):
        return abs(self.pt_src_T - self.pt_src_tau)

    @cached_property
    def clim_pt_src_tau(self):
        cmax = abs(self.error_pt_src_tau).max()
        return 0, cmax

    def plot_error_pt_src_tau(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.error_pt_src_tau,
            clim=self.clim_pt_src_tau,
            title=r'$T - \tau$ [m] (point source)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_pt_src=True,
        )

    # grad(eikonal)

    @cached_property
    def pt_src_grad_T(self):
        return self.pt_src_field.eik.grad_T

    @cached_property
    def pt_src_grad_tau(self):
        mask = self.pt_src_mask
        grad_tau = np.empty_like(self.verts)
        grad_tau[mask] = self.verts[mask] - self.xsrc
        grad_tau[mask] /= self.dist2xsrc[mask].reshape(-1, 1)
        grad_tau[~mask] = self.grad_of_dist2xsrc_plus_dist2edge[~mask]
        return grad_tau

    @cached_property
    def error_pt_src_grad_tau(self):
        prod = self.pt_src_grad_T*self.pt_src_grad_tau
        return np.arccos(np.clip(np.sum(prod, axis=1), -1, 1))

    @cached_property
    def clim_pt_src_grad_tau(self):
        return (0, self.error_pt_src_grad_tau.max())

    def plot_error_pt_src_grad_tau(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.error_pt_src_grad_tau,
            clim=self.clim_pt_src_grad_tau,
            title=r'$\cos^{-1}(\nabla T \cdot \nabla \tau)$ [m] (point source)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_pt_src=True
        )

    # amplitude

    @cached_property
    def pt_src_amp(self):
        return abs(self.pt_src_field.amplitude)

    # TODO: plot_pt_src_amp

    @cached_property
    def pt_src_amp_gt(self):
        amp = 1/self.pt_src_tau
        amp[self.xsrc_index] = np.nan
        return amp

    @cached_property
    def error_pt_src_amp(self):
        return abs(self.pt_src_amp - self.pt_src_amp_gt)

    # TODO: plot_error_pt_src_amp

    # sound pressure level

    @cached_property
    def pt_src_SPL(self):
        return 20*np.log10(abs(self.pt_src_amp))

    def plot_pt_src_SPL(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.pt_src_SPL,
            clim=(self.min_dB, np.nanmax(self.pt_src_SPL)),
            title=r'SPL [dB] (point source)',
            cmap=self.dB_cmap,
            unit_str='dB',
            plot_pt_src=True
        )

    @cached_property
    def pt_src_SPL_gt(self):
        return 20*np.log10(abs(self.pt_src_amp_gt))

    @cached_property
    def error_pt_src_SPL(self):
        return abs(self.pt_src_SPL - self.pt_src_SPL_gt)

    # TODO: plot_error_pt_src_SPL

    ########################################################################
    # Near reflection fields

    # eikonal

    @cached_property
    def near_refl_T(self):
        return self.near_refl_field.eik.T

    def plot_near_refl_T(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.near_refl_T,
            clim=(np.nanmin(self.near_refl_T), np.nanmax(self.near_refl_T)),
            title=r'$T$ [m] (near reflection)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_pt_src=True,
            plot_near_img=True,
            plot_near_refl_bd=True
        )

    @cached_property
    def near_refl_tau(self):
        if not self.should_do_near_reflection:
            raise RuntimeError('not doing near reflection')
        mask = self.near_refl_mask
        tau = np.empty(self.verts.shape[0], dtype=self.verts.dtype)
        tau[mask] = self.near_refl_dist2img[mask]
        tau[~mask] = self.dist2xsrc_plus_dist2edge[~mask]
        return tau

    @cached_property
    def error_near_refl_tau(self):
        if not self.should_do_near_reflection:
            raise RuntimeError('not doing near reflection')
        return abs(self.near_refl_T - self.near_refl_tau)

    @cached_property
    def clim_near_refl_tau(self):
        cmax = abs(self.error_near_refl_tau).max()
        return 0, cmax

    def plot_error_near_refl_tau(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.error_near_refl_tau,
            clim=self.clim_near_refl_tau,
            title=r'$T - \tau$ [m] (near reflection)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_pt_src=True,
            plot_near_img=True,
            plot_near_refl_bd=True
        )

    # grad(eikonal)

    @cached_property
    def near_refl_grad_tau(self):
        if not self.should_do_near_reflection:
            raise RuntimeError('not doing near reflection')
        mask = self.near_refl_mask
        grad_tau = np.empty_like(self.verts)
        xsrc_refl = self.reflect_over_near_reflector(self.xsrc)
        grad_tau[mask] = self.near_refl_grad_of_dist2img[mask]
        grad_tau[~mask] = self.grad_of_dist2xsrc_plus_dist2edge[~mask]
        return grad_tau

    @cached_property
    def error_near_refl_grad_tau(self):
        if not self.should_do_near_reflection:
            raise RuntimeError('not doing near reflection')
        prod = self.near_refl_field.eik.grad_T*self.near_refl_grad_tau
        return np.arccos(np.clip(np.sum(prod, axis=1), -1, 1))

    @cached_property
    def clim_near_refl_grad_tau(self):
        return (0, self.error_near_refl_grad_tau.max())

    def plot_error_near_refl_grad_tau(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.error_near_refl_grad_tau,
            clim=self.clim_near_refl_grad_tau,
            title=r'$\cos^{-1}(\nabla T\cdot\nabla\tau)$ [m] (near reflection)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_pt_src=True,
            plot_near_img=True,
            plot_near_refl_bd=True
        )

    # amplitude

    @cached_property
    def near_refl_amp(self):
        return abs(self.near_refl_field.amplitude)

    # TODO: plot_near_refl_amp

    # TODO: near_refl_amp_gt

    # TODO: error_near_refl_amp

    # TODO: plot_error_near_refl_amp

    # sound pressure level

    @cached_property
    def near_refl_SPL(self):
        return 20*np.log10(abs(self.near_refl_amp))

    def plot_near_refl_SPL(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.near_refl_SPL,
            clim=(self.min_dB, np.nanmax(self.near_refl_SPL)),
            title=r'SPL (near reflection)',
            cmap=self.dB_cmap,
            unit_str='dB',
            plot_pt_src=True,
            plot_near_img=True,
            plot_near_refl_bd=True
        )

    # TODO: near_refl_SPL_gt

    # TODO: error_near_refl_SPL

    # TODO: plot_error_near_refl_SPL

    ########################################################################
    # Far reflection fields

    # eikonal

    @cached_property
    def far_refl_T(self):
        return self.far_refl_field.eik.T

    def plot_far_refl_T(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.far_refl_T,
            clim=(np.nanmin(self.far_refl_T), np.nanmax(self.far_refl_T)),
            title=r'$T$ [m] (far reflection)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters'
        )

    @cached_property
    def far_refl_tau(self):
        if not self.should_do_far_reflection:
            raise RuntimeError('not doing far reflection')
        mask = self.far_refl_mask
        tau = np.empty(self.verts.shape[0], dtype=self.verts.dtype)
        tau[mask] = self.far_refl_dist2img[mask]
        tau[~mask] = self.dist2xsrc_plus_dist2edge[~mask]
        return tau

    @cached_property
    def error_far_refl_tau(self):
        if not self.should_do_far_reflection:
            raise RuntimeError('not doing far reflection')
        return abs(self.far_refl_tau - self.far_refl_field.eik.T)

    @cached_property
    def clim_far_refl_tau(self):
        cmax = abs(self.error_far_refl_tau).max()
        return 0, cmax

    def plot_error_far_refl_tau(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.error_far_refl_tau,
            clim=self.clim_far_refl_tau,
            title=r'$T - \tau$ (far reflection)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_pt_src=True,
            plot_far_img=True,
            plot_far_refl_bd=True
        )

    # grad(eikonal)

    @cached_property
    def far_refl_grad_tau(self):
        if not self.should_do_far_reflection:
            raise RuntimeError('not doing far reflection')
        mask = self.far_refl_mask
        grad_tau = np.empty_like(self.verts)
        xsrc_refl = self.reflect_over_far_reflector(self.xsrc)
        grad_tau[mask] = self.far_refl_grad_of_dist2img[mask]
        grad_tau[~mask] = self.grad_of_dist2xsrc_plus_dist2edge[~mask]
        return grad_tau

    @cached_property
    def error_far_refl_grad_tau(self):
        if not self.should_do_far_reflection:
            raise RuntimeError('not doing far reflection')
        prod = self.far_refl_field.eik.grad_T*self.far_refl_grad_tau
        return np.arccos(np.clip(np.sum(prod, axis=1), -1, 1))

    @cached_property
    def clim_far_refl_grad_tau(self):
        return (0, self.error_far_refl_grad_tau.max())

    def plot_error_far_refl_grad_tau(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.error_far_refl_grad_tau,
            clim=self.clim_far_refl_grad_tau,
            title=r'$\cos^{-1}(\nabla T \cdot \nabla \tau)$ (far reflection)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_pt_src=True,
            plot_far_img=True,
            plot_far_refl_bd=True
        )

    # amplitude

    @cached_property
    def far_refl_amp(self):
        return abs(self.far_refl_field.amplitude)

    # TODO: plot_far_refl_amp

    # TODO: far_refl_amp_gt

    # TODO: error_far_refl_amp

    # TODO: plot_error_far_refl_amp

    # sound pressure level

    @cached_property
    def far_refl_SPL(self):
        return 20*np.log10(abs(self.far_refl_amp))

    def plot_far_refl_SPL(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.far_refl_SPL,
            clim=(self.min_dB, np.nanmax(self.far_refl_SPL)),
            title=r'SPL (far reflection)',
            cmap=self.dB_cmap,
            unit_str='dB',
            plot_pt_src=True,
            plot_far_img=True,
            plot_far_refl_bd=True
        )

    # TODO: far_refl_SPL_gt

    # TODO: error_far_refl_SPL

    # TODO: plot_error_far_refl_SPL

    ########################################################################
    # Diffracted fields

    # eikonal

    @cached_property
    def diff_T(self):
        return self.diff_field.eik.T

    def plot_diff_T(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.diff_T,
            clim=(np.nanmin(self.diff_T), np.nanmax(self.diff_T)),
            title=r'$T$ [m] (diffracted field)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_diff_src=True,
            plot_near_refl_bd=True,
            plot_far_refl_bd=True,
        )

    @cached_property
    def diff_tau(self):
        return self.dist2xsrc_plus_dist2edge

    @cached_property
    def error_diff_tau(self):
        return abs(self.diff_tau - self.diff_field.eik.T)

    @cached_property
    def clim_diff_tau(self):
        cmax = abs(self.error_diff_tau).max()
        return 0, cmax

    def plot_error_diff_tau(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.error_diff_tau,
            clim=self.clim_diff_tau,
            title=r'$T - \tau$ (diffracted field)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_diff_src=True,
            plot_near_refl_bd=True,
            plot_far_refl_bd=True,
        )

    # grad(eikonal)

    @cached_property
    def diff_grad_tau(self):
        return self.grad_of_dist2xsrc_plus_dist2edge

    @cached_property
    def error_diff_grad_tau(self):
        prod = self.diff_field.eik.grad_T*self.diff_grad_tau
        return np.arccos(np.clip(np.sum(prod, axis=1), -1, 1))

    @cached_property
    def clim_diff_grad_tau(self):
        return (0, self.error_diff_grad_tau.max())

    def plot_error_diff_grad_tau(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.error_diff_grad_tau,
            clim=self.clim_diff_grad_tau,
            title=r'$\cos^{-1}(\nabla T\cdot\nabla\tau)$ (diffracted field)',
            cmap=cc.cm.CET_L18,
            unit_str='Meters',
            plot_diff_src=True,
            plot_near_refl_bd=True,
            plot_far_refl_bd=True,
        )

    # amplitude

    @cached_property
    def diff_amp(self):
        return abs(self.diff_field.amplitude)

    # TODO: plot_diff_amp

    # TODO: diff_amp_gt

    # TODO: error_diff_amp

    # TODO: plot_error_diff_amp

    # sound pressure level

    @cached_property
    def diff_SPL(self):
        return 20*np.log10(abs(self.diff_amp))

    def plot_diff_SPL(self, fig, ax):
        self.make_scatter_plot(
            fig, ax,
            values=self.diff_SPL,
            clim=(self.min_dB, np.nanmax(self.diff_SPL)),
            title=r'SPL (diffracted field)',
            cmap=self.dB_cmap,
            unit_str='dB',
            plot_diff_src=True,
            plot_near_refl_bd=True,
            plot_far_refl_bd=True,
        )

    # TODO: diff_SPL_gt

    # TODO: error_diff_SPL

    # TODO: plot_error_diff_SPL
