import numpy as np

# TODO: I don't think the conventions I've chosen for the camera basis
# are actually very conventional...

def get_camera_basis(origin, target, up):
    '''Get a triple of a vectors (left, front, up), indicating an
    orthonormal basis for the camera.

    '''
    front = target - origin
    front /= np.linalg.norm(front)
    left = np.cross(up, front)
    left /= np.linalg.norm(left)
    up = np.cross(front, left)
    up /= np.linalg.norm(up)
    return left, front, up

def get_view_direction(left, front, up, phi, theta):
    '''Get the view direction in the camera basis expressed by (left,
    front, up) corresponding to the spherical angles phi and
    theta. This is written so that:

    `get_view_direction(*get_camera_basis(origin, target, up), phi, theta)`

    works as expected.

    '''

    ray = np.array([0, 1, 0])
    # Rotate by phi around z/up axis
    ray = np.array([
        [ np.cos(phi), np.sin(phi), 0],
        [-np.sin(phi), np.cos(phi), 0],
        [          0,            0, 1]
    ])@ray
    # Rotate around the x/left axis.
    ray = np.array([
        [1,              0,             0],
        [0,  np.cos(theta), np.sin(theta)],
        [0, -np.sin(theta), np.cos(theta)],
    ])@ray
    return np.array([left, front, up]).T@ray
