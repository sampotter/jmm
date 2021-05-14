import numpy as np

from enum import Enum

from jmm.geom cimport ray3, Ray3

class CameraType(Enum):
    Uninitialized = 0
    Orthographic = 1
    Perspective = 2

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

cdef class Camera:
    @staticmethod
    def make_orthographic(pos, look, left, up, shape, width, height):
        camera = Camera()
        camera._camera.type = CAMERA_TYPE_ORTHOGRAPHIC
        cdef int i
        for i in range(3):
            camera._camera.pos[i] = pos[i]
            camera._camera.look[i] = look[i]
            camera._camera.left[i] = left[i]
            camera._camera.up[i] = up[i]
        camera._camera.dim[0] = shape[0]
        camera._camera.dim[1] = shape[1]
        camera._camera.width = width
        camera._camera.height = height
        return camera

    @property
    def shape(self):
        return (self._camera.dim[0], self._camera.dim[1])

    @property
    def camera_type(self):
        return CameraType(self._camera.type)

    @property
    def pos(self):
        return (self._camera.pos[0], self._camera.pos[1], self._camera.pos[2])

    @property
    def look(self):
        return np.asarray(self._camera.look)

    @property
    def left(self):
        return np.asarray(self._camera.left)

    @property
    def up(self):
        return np.asarray(self._camera.up)

    @property
    def width(self):
        if self.camera_type != CameraType.Orthographic:
            raise Exception('tried to access "width" of %s camera' % (
                self._get_camera_type_name(),))
        return self._camera.width

    @property
    def height(self):
        if self.camera_type != CameraType.Orthographic:
            raise Exception('tried to access "height" of %s camera' % (
                self._get_camera_type_name(),))
        return self._camera.height

    @property
    def extent(self):
        '''A guess for the Matplotlib `imshow` plotting function's `extent`
        keyword argument. The extent is computed as the axes of the
        image plane.

        '''
        if self.camera_type == CameraType.Uninitialized:
            raise Exception("can't compute extent of uninitialized camera")
        elif self.camera_type == CameraType.Orthographic:
            w, h = self.width, self.height
            return (-w/2, w/2, -h/2, h/2)
        elif self.camera_type == CameraType.Perspective:
            raise Exception("extent not implemented for perspective camera yet")
        else:
            assert False

    @property
    def num_rays(self):
        return np.product(self._camera.dim)

    def get_image_rays(self):
        cdef ray3 ray
        cdef size_t i, j
        for i in range(self._camera.dim[0]):
            for j in range(self._camera.dim[1]):
                ray = camera_get_ray_for_index(&self._camera, i, j)
                yield Ray3.from_ptr(&ray)

    def _get_camera_type_name(self):
        return str(self.camera_type).split('.')[1].lower()
