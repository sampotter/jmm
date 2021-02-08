import itertools as it
import jmm
import meshio
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
    trimesh = meshio.read('L.obj')
    verts = trimesh.points.astype(np.float64).copy()
    faces = trimesh.cells_dict['triangle'].astype(np.uintp).copy()
    mesh = jmm.Mesh2(verts, faces)
    rtree = jmm.Rtree(mesh)

    print('R-tree bounding box: %s' % rtree.bounding_box)

    # Compute camera frame (front x left x up)
    origin = np.array([1, 1, 5], dtype=np.float64)
    target = np.array([0, 0, 0], dtype=np.float64)
    front = target - origin
    front /= np.linalg.norm(front)
    up = np.array([0, 0, 1])
    left = np.cross(up, front)
    left /= np.linalg.norm(left)
    up = np.cross(front, left)
    up /= np.linalg.norm(up)

    ray = jmm.Ray3(origin, front)
    isect = rtree.intersect(ray)
    print(isect)

    print('left', left)
    print('front', front)
    print('up', up)

    width = 640
    height = 480

    Phi_deg = np.linspace(-15, 15, width)
    Theta_deg = np.linspace(-7.5, 7.5, height)

    img = Image.new('RGB', (width, height), 'white')

    for (i, phi_deg), (j, theta_deg) in it.product(enumerate(Phi_deg),
                                                   enumerate(Theta_deg)):
        phi, theta = np.deg2rad(phi_deg), np.deg2rad(theta_deg)
        direction = get_view_direction(left, front, up, phi, theta)
        ray = jmm.Ray3(origin, direction)
        isect = rtree.intersect(ray)
        if isect is not None:
            img.putpixel((i, j), (0, 0, 0, 0))

    img.show()
