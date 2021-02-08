import colorcet as cc
import itertools as it
import jmm
import meshio
import meshzoo
import numpy as np

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
    # trimesh = meshio.read('L.obj')
    # verts = trimesh.points.astype(np.float64).copy()
    # faces = trimesh.cells_dict['triangle'].astype(np.uintp).copy()

    verts, faces = meshzoo.icosa_sphere(2)
    faces = faces.astype(np.uintp)

    mesh = jmm.Mesh2(verts, faces)
    rtree = jmm.Rtree(mesh)

    print('R-tree bounding box: %s' % rtree.bounding_box)

    # Compute camera frame (front x left x up)
    origin = np.array([-10, 0, 0], dtype=np.float64)
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

    T, I = rtree.intersectN(orgs, dirs)
    T = T.reshape(width, height)
    I = I.reshape(width, height)

    img = Image.new('RGB', (width, height), 'white')

    # Plot intersection distance
    tmax = np.nanmax(T)
    for i, j in it.product(range(width), range(height)):
        t = T[i, j]
        if np.isnan(t):
            continue
        rgba = tuple(int(np.round(255*_)) for _ in cc.cm.rainbow(t/tmax))
        img.putpixel((i, j), rgba)

    # # Plot index of intersected triangle
    # imax = I[~np.isnan(T)].max()
    # for i, j in it.product(range(width), range(height)):
    #     i_ = I[i, j]
    #     if np.isnan(T[i, j]):
    #         continue
    #     rgba = tuple(int(np.round(255*_)) for _ in cc.cm.rainbow(i_/imax))
    #     img.putpixel((i, j), rgba)

    img.show()
