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

def get_view_direction(left, front, up, phi, theta):
    ray = np.array([0, 1, 0])
    # Rotate by phi around z/up axis
    ray = np.array([
        [ np.cos(phi), np.sin(phi), 0],
        [-np.sin(phi), np.cos(phi), 0],
        [          0,            0, 1]
    ])@ray
    ray = np.array([
        [1,              0,             0],
        [0,  np.cos(theta), np.sin(theta)],
        [0, -np.sin(theta), np.cos(theta)],
    ])@ray
    return np.array([left, front, up]).T@ray

if __name__ == '__main__':
    plt.ion()

    verts, faces = meshzoo.icosa_sphere(8)
    faces = faces.astype(np.uintp)

    mesh = jmm.Mesh2(verts, faces)
    rtree = jmm.Rtree.from_mesh2(mesh)

    print('R-tree bounding box: %s' % rtree.bounding_box)

    # Compute camera frame (front x left x up)
    origin = np.array([-4, 0, 0], dtype=np.float64)
    target = np.array([0, 0, 0], dtype=np.float64)
    front = target - origin
    front /= np.linalg.norm(front)
    up = np.array([0, 0, 1])
    left = np.cross(up, front)
    left /= np.linalg.norm(left)
    up = np.cross(front, left)
    up /= np.linalg.norm(up)

    isect = rtree.intersect(origin, front)
    print(isect)

    print('left', left)
    print('front', front)
    print('up', up)

    width = 256
    height = 256
    num_rays = width*height

    Phi = np.deg2rad(np.linspace(-15, 15, width))
    Theta = np.deg2rad(np.linspace(-15, 15, height))

    orgs = np.outer(np.ones(num_rays), origin)
    dirs = np.array([get_view_direction(left, front, up, phi, theta)
                     for phi, theta in it.product(Phi, Theta)])

    t0 = time.time()

    T, I = rtree.intersectN(orgs, dirs)
    T = T.reshape(width, height)
    I = I.reshape(width, height).astype(np.float64)
    I[np.isnan(T)] = np.nan

    print('- our BVH elapsed time: %1.2f' % (time.time() - t0,))

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

    I_gt = rayhit.prim_id.reshape(width, height).astype(np.float64)
    T_gt = rayhit.tfar.reshape(width, height)
    I_gt[~np.isfinite(T_gt)] = np.nan

    print('- Embree elapsed time: %1.2f' % (time.time() - t0,))

    plt.figure(figsize=(10, 7))
    plt.subplot(2, 3, 1)
    plt.imshow(I_gt, cmap=cc.cm.rainbow)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.subplot(2, 3, 4)
    plt.imshow(T_gt, cmap=cc.cm.fire)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.subplot(2, 3, 2)
    plt.imshow(I, cmap=cc.cm.rainbow)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.subplot(2, 3, 5)
    plt.imshow(T, cmap=cc.cm.fire)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.subplot(2, 3, 3)
    E_I = I - I_gt
    I_vmax = np.nanmax(abs(E_I))
    plt.imshow(E_I, cmap=cc.cm.coolwarm, vmin=-I_vmax, vmax=I_vmax)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.subplot(2, 3, 6)
    E_T = T - T_gt
    T_vmax = np.nanmax(abs(E_T))
    plt.imshow(E_T, cmap=cc.cm.coolwarm, vmin=-T_vmax, vmax=T_vmax)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.tight_layout()
    plt.show()
