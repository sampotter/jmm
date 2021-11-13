import colorcet as cc
import itertools as it
import jmm
import matplotlib.pyplot as plt
import meshio
import meshzoo
import numpy as np
import embree
import time

from PIL import Image

from common import get_camera_basis, get_view_direction

if __name__ == '__main__':
    plt.ion()

    verts, faces = meshzoo.icosa_sphere(20)
    faces = faces.astype(np.uintp)
    num_faces = faces.shape[0]

    mesh = jmm.Mesh2.from_verts_and_faces(verts, faces)
    rtree = jmm.Rtree.from_mesh2(mesh)

    print('R-tree bounding box: %s' % rtree.bounding_box)

    # Compute camera frame (front x left x up)
    origin = np.array([-4, 0, 0], dtype=np.float64)
    target = np.array([0, 0, 0], dtype=np.float64)
    up = np.array([0, 0, 1])
    left, front, up = get_camera_basis(origin, target, up)

    # Image size
    w, h = 640, 480
    num_rays = w*h

    # Camera parameters
    s = w/h # aspect ratio
    fov_y = 30 # vertical field of view

    Theta = np.deg2rad(np.linspace(-s*fov_y/2, s*fov_y/2, w))
    Phi = np.deg2rad(np.linspace(-fov_y/2, fov_y/2, h))

    orgs = np.outer(np.ones(num_rays), origin)
    dirs = np.array([get_view_direction(left, front, up, phi, theta)
                     for phi, theta in it.product(Phi, Theta)])

    t0 = time.time()

    T = np.empty((num_rays,), dtype=np.float64)
    I = np.empty((num_rays,), dtype=np.float64)
    for l, (org, dir_) in enumerate(zip(orgs, dirs)):
        isect = rtree.intersect(org, dir_)
        T[l] = isect.t if isect.hit else np.nan
        if isect.hit:
            I[l] = isect.obj.astype(jmm.Mesh2Tri).index
        else:
            I[l] = np.nan
    T = T.reshape(h, w)
    I = I.reshape(h, w)

    print('Our BVH elapsed time: %1.2f' % (time.time() - t0,))

    # Compare against python-embree

    t0 = time.time()

    device = embree.Device()
    geometry = device.make_geometry(embree.GeometryType.Triangle)
    scene = device.make_scene()
    vertex_buffer = geometry.set_new_buffer(
        embree.BufferType.Vertex, 0, embree.Format.Float3,
        3*np.dtype('float32').itemsize, verts.shape[0])
    vertex_buffer[...] = verts[...]
    index_buffer = geometry.set_new_buffer(
        embree.BufferType.Index, 0, embree.Format.Uint3,
        3*np.dtype('uint32').itemsize, faces.shape[0])
    index_buffer[...] = faces[...]
    geometry.commit()
    scene.attach_geometry(geometry)
    geometry.release()
    scene.commit()

    rayhit = embree.RayHit1M(num_rays)
    context = embree.IntersectContext()
    rayhit.org[...] = orgs
    rayhit.dir[...] = dirs
    rayhit.tnear[...] = 0
    rayhit.tfar[...] = np.inf
    rayhit.flags[...] = 0
    rayhit.geom_id[...] = embree.INVALID_GEOMETRY_ID
    scene.intersect1M(context, rayhit)

    I_gt = rayhit.prim_id.reshape(h, w).astype(np.float64)
    T_gt = rayhit.tfar.reshape(h, w)
    T_gt[~np.isfinite(T_gt)] = np.nan
    I_gt[~np.isfinite(T_gt)] = np.nan

    print('Embree elapsed time: %1.2f' % (time.time() - t0,))

    T_min = min(np.nanmin(T_gt), np.nanmin(T))
    T_max = max(np.nanmax(T_gt), np.nanmax(T))

    plt.figure(figsize=(13, 6))

    plt.subplot(2, 3, 1)
    plt.imshow(I_gt, cmap=cc.cm.glasbey, interpolation='none')
    plt.gca().set_aspect('equal')
    plt.title('Hit index (Embree)')

    plt.subplot(2, 3, 4)
    plt.imshow(T_gt, cmap=cc.cm.bmw, vmin=T_min, vmax=T_max)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title('Hit distance (Embree)')

    plt.subplot(2, 3, 2)
    plt.imshow(I, cmap=cc.cm.glasbey, interpolation='none')
    plt.gca().set_aspect('equal')
    plt.title("Hit index (jmm's R-tree)")

    plt.subplot(2, 3, 5)
    plt.imshow(T, cmap=cc.cm.bmw, vmin=T_min, vmax=T_max)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title("Hit distance (jmm's R-tree)")

    plt.subplot(2, 3, 3)
    E_I = I - I_gt
    plt.imshow((E_I != 0).astype(np.float64), cmap=cc.cm.fire_r,
               interpolation='none')
    plt.gca().set_aspect('equal')
    plt.title('Hit indices equal? (mask)')

    plt.subplot(2, 3, 6)
    E_T = T - T_gt
    T_vmax = np.nanmax(abs(E_T))
    plt.imshow(E_T, cmap=cc.cm.coolwarm, vmin=-T_vmax, vmax=T_vmax)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title('Hit distance error (jmm - Embree)')

    plt.tight_layout()

    plt.show()
