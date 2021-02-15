#!/usr/bin/env python

import colorcet as cc
import embree
import itertools as it
import jmm
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.optimize

from PIL import Image

EPS = 1e-5

plt.ion()

def get_bary_coords(P, x):
    aug = np.empty((4, 4))
    aug[0] = 1
    aug[1:] = P.T
    vol = np.linalg.det(aug)/6
    coords = np.empty(4)
    for i in range(4):
        aug[1:, i] = x
        coords[i] = np.linalg.det(aug)/(6*vol)
        aug[1:, i] = P[i]
    return coords

class BezierTetraMesh:
    def __init__(self, eik, value=None):
        self._eik = eik
        self._device = embree.Device()
        self._value = None
        self.value = value

    @property
    def mesh(self):
        return self._eik.mesh

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        if self._value == value:
            return
        self._value = value
        self._build_bvh()

    def _build_bvh(self):
        cells = self.mesh.cells

        # Mask the tetrahedra which bracket the slice value
        #
        # TODO: actually, we need to do this in terms of the convex
        # hull of the control points of the Bezier tetrahedra... but
        # let's do it this simplified way now until we finish making a
        # video
        T = np.array([_[0] for _ in self._eik.jet])
        mask = (T[cells] <= self.value).any(1) \
             & (self.value <= T[cells]).any(1)

        masked_cells = cells[mask]
        num_masked_cells = masked_cells.shape[0]

        num_faces = 4*num_masked_cells
        faces = np.empty((num_faces, 3), dtype=cells.dtype)
        faces[0::4] = masked_cells[:, [0, 1, 2]]
        faces[1::4] = masked_cells[:, [0, 1, 3]]
        faces[2::4] = masked_cells[:, [0, 2, 3]]
        faces[3::4] = masked_cells[:, [1, 2, 3]]

        self._scene = self._device.make_scene()
        geometry = self._device.make_geometry(embree.GeometryType.Triangle)
        vertex_buffer = geometry.set_new_buffer(
            embree.BufferType.Vertex, 0, embree.Format.Float3,
            3*np.dtype('float32').itemsize, self.mesh.num_verts)
        vertex_buffer[... ] = self.mesh.verts[...]
        index_buffer = geometry.set_new_buffer(
            embree.BufferType.Index, 0, embree.Format.Uint3,
            3*np.dtype('uint32').itemsize, num_faces)
        index_buffer[...] = faces[...]
        geometry.commit()
        self._scene.attach_geometry(geometry)
        geometry.release()
        self._scene.commit()

        self._bezier_tetra = np.array([
            self._eik.get_bezier_tetra(l) for l in np.where(mask)[0]])

    def get_intersecting_cells(self, orgs, dirs):
        if not self._value:
            raise Exception('level set (`value`) is unset')

        num_rays = orgs.shape[0]

        isect_lists = np.array([None for _ in range(num_rays)], dtype=object)

        rayhit = embree.RayHit1M(num_rays)
        context = embree.IntersectContext()
        rayhit.org[...] = orgs
        rayhit.dir[...] = dirs
        rayhit.tnear[...] = 0
        rayhit.tfar[...] = np.inf
        rayhit.flags[...] = 0
        rayhit.geom_id[...] = embree.INVALID_GEOMETRY_ID
        self._scene.intersect1M(context, rayhit)

        for i in range(num_rays):
            if ~np.isfinite(rayhit.tfar[i]):
                continue
            prim_id = rayhit.prim_id[i]//4
            tfar = rayhit.tfar[i]
            isect = (prim_id, tfar)
            if not isect_lists[i]:
                isect_lists[i] = [isect]
            else:
                isect_lists[i].append(isect)

        return isect_lists

    def get_values(self, orgs, dirs):
        if not self._value:
            raise Exception('level set (`value`) is unset')

        num_rays = orgs.shape[0]

        isect_lists = np.array([None for _ in range(num_rays)], dtype=object)

        rayhit = embree.RayHit1M(num_rays)
        context = embree.IntersectContext()
        rayhit.org[...] = orgs
        rayhit.dir[...] = dirs
        rayhit.tnear[...] = 0
        rayhit.tfar[...] = np.inf
        rayhit.flags[...] = 0
        rayhit.geom_id[...] = embree.INVALID_GEOMETRY_ID
        self._scene.intersect1M(context, rayhit)

        values = np.empty(num_rays, dtype=np.float64)
        values[...] = np.nan

        X = np.empty((num_rays, 3), dtype=np.float64)
        X[... ] = np.nan

        for i in range(num_rays):
            if ~np.isfinite(rayhit.tfar[i]):
                continue
            prim_id = rayhit.prim_id[i]//4
            t0 = rayhit.tfar[i]
            x0 = orgs[i] + (t0 + EPS)*dirs[i]

            rayhit_ = embree.RayHit()
            rayhit_.org = x0
            rayhit_.dir = dirs[i]
            rayhit_.tnear = 0
            rayhit_.tfar = np.inf
            rayhit_.geom_id = embree.INVALID_GEOMETRY_ID
            context_ = embree.IntersectContext()
            self._scene.intersect1(context_, rayhit_)
            t1 = rayhit_.tfar/(1 - EPS) # correct for offset
            if not np.isfinite(t1):
                raise Exception(f'bad value: t1 = {t1}')

            P = self.mesh.verts[self.mesh.cells[prim_id]]
            x1 = orgs[i] + t1*dirs[i]
            bb = self._bezier_tetra[prim_id]
            value = self.value
            f = lambda t: bb.f(get_bary_coords(P, (1 - t)*x0 + t*x1)) - value
            try:
                topt = scipy.optimize.brentq(f, 0, 1)
            except:
                continue
            values[i] = f(topt) + value
            X[i] = (1 - topt)*x0 + topt*x1

        return values, X

if __name__ == '__main__':
    path = '.'
    scene = 'box'
    indsrc = 16
    Tslice = 0.5

    verts_bin_path = os.path.join(path, scene + '_verts.bin')
    verts = np.fromfile(verts_bin_path, np.float64)
    verts = verts.reshape(verts.size//3, 3)
    print('- read verts from %s' % verts_bin_path)

    cells_bin_path = os.path.join(path, scene + '_cells.bin')
    cells = np.fromfile(cells_bin_path, np.uintp)
    cells = cells.reshape(cells.size//4, 4)
    print('- reading cells from %s' % cells_bin_path)

    mesh = jmm.Mesh3.from_verts_and_cells(verts, cells)
    eik = jmm.Eik3(mesh)
    eik.add_trial(indsrc, jmm.Jet3(0, np.nan, np.nan, np.nan))
    eik.solve()

    btm = BezierTetraMesh(eik, value=Tslice)

    width = 512
    height = 512
    num_rays = width*height
    x, z = np.meshgrid(np.linspace(-1.1, 1.1, width),
                       np.linspace(-1.1, 1.1, height))
    orgs = np.array([x.ravel(), -1.1*np.ones(num_rays), z.ravel()]).T
    dirs = np.outer(np.ones(num_rays), np.array([0, 1, 0]))

    # isect_lists = btm.get_intersecting_cells(orgs, dirs).reshape(width, height)
    # img = Image.new('RGB', (width, height), 'white')
    # for i, j in it.product(range(width), range(height)):
    #     lst = isect_lists[i, j]
    #     if not lst:
    #         continue
    #     ind, t = lst[0]
    #     rgba = tuple(
    #         int(np.round(255*c))
    #         for c in cc.cm.rainbow(ind/mesh.num_cells))
    #     img.putpixel((i, j), rgba)
    # img.show()

    values, X = btm.get_values(orgs, dirs)
    values = values.reshape(width, height)
    values_gt = np.sqrt(np.sum(X**2, axis=1)).reshape(width, height)

    plt.figure()
    plt.imshow(values - values_gt, cmap=cc.cm.fire_r)
    plt.gca().set_aspect('equal')
    plt.show()
