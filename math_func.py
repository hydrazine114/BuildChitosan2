from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import cosine
import numpy as np


vec_norm = lambda x, y, z: (x ** 2 + y ** 2 + z ** 2) ** 0.5


def vec_ang(vec1, vec2, degree=False):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = np.dot(vec1, vec2)
    sinang = np.linalg.norm(np.cross(vec1, vec2))
    if degree:
        return np.arctan2(sinang, cosang) * 180 / np.pi
    return np.arctan2(sinang, cosang)


def create_rotvec(vec1, vec2):
    """ Create rotation vector to rotate vector 1 to vector 2"""
    angle = vec_ang(vec1, vec2)
    vec3 = np.cross(vec1, vec2)
    vec3 *= angle / vec_norm(*vec3)
    return vec3


def rotate_structure(structure, str_vector, oriented_vector):
    rot_vec = create_rotvec(np.array(str_vector, dtype='float64'), np.array(oriented_vector, dtype='float64'))
    rotation = R.from_rotvec(rot_vec)
    for i, vector in enumerate(structure):
        structure[i] = rotation.apply(vector)
    return structure


def rotate_struc_euler(structure, anglex=0, anglez=0, angley=0):
    rotation = R.from_euler('x', anglex, degrees=True)
    for i, vector in enumerate(structure):
        structure[i] = rotation.apply(vector)

    rotation = R.from_euler('y', angley, degrees=True)
    for i, vector in enumerate(structure):
        structure[i] = rotation.apply(vector)

    rotation = R.from_euler('z', anglez, degrees=True)
    for i, vector in enumerate(structure):
        structure[i] = rotation.apply(vector)

    return structure


if __name__ == '__main__':
    structure = [[1, 1, 0],
                 [1, -1, 0],
                 [-1, -1, 0],
                 [-1, 1, 0]]
    rotate_structure(structure, structure[0], [1, 0, 0])
    for i in structure:
        print(i)
