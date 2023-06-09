#!/usr/bin/env python3
# coding=utf-8
"""
Math functions. Mostly for vector operations.
@Author: Mads Jeppesen
@Date: 9/21/22
"""
import math
import numpy as np
import warnings

def rotation_matrix_from_vector_to_vector(v1, v2):
    """Generates a rotation matrix that rotates v1 to v2

    :param ndarray v1: vector to rotate from. Shape(3,).
    :param ndarray v2: vector to rotate to. Shape(3,).
    """
    # normalize:
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    # if v1 and v2 are identical, return an identity matrix
    if np.isclose(v1, v2).all():
        warnings.warn("v1 and v2 are similar and an identity rotation matrix is returned!")
        return np.identity(3)
    # actual calculation
    axis = np.cross(v1,v2)
    angle = vector_angle(v1, v2)
    return rotation_matrix(axis, angle)

def vector_angle(vec1, vec2, degrees=True):
    """Calculates the angle between two vectors.

    :param ndarray vec1: The first vector.
    :param ndarray vec2: The second vector.
    :param degrees: Return the results in degrees. If false, returns in radians
    :return float: The angle in degrees.
    """
    denom = (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    # if the denominator is zero this will return we should not perform the division and we return 0.0.
    # This will happen without this if statement but it will throw a RuntimeWarning everytime it happens
    if math.isclose(denom, 0.0):
        return 0.0
    else:
        cos_angle = np.dot(vec1, vec2) / denom
        val = np.arccos(np.clip(cos_angle, -1, 1))
    if degrees:
        return math.degrees(val)
    else: # radians
        return val

def rotation_matrix(axis, angle):
    """Generates a rotation matrix for a rotation about an arbitrary axis.

    :param ndarray axis: Axis to rotate about. Shape(3,).
    :param float angle: The angle of the rotation in degrees.
    :return ndarray: The rotation matrix. Shape(3,3).

    """
    unit = axis / np.linalg.norm(axis)
    cos_theta = math.cos(math.radians(angle))
    one_minus_cos_theta = 1.0 - cos_theta
    sin_theta = math.sin(math.radians(angle))
    xx = cos_theta + unit[0] * unit[0] * one_minus_cos_theta
    xy = unit[0] * unit[1] * one_minus_cos_theta - unit[2] * sin_theta
    xz = unit[0] * unit[2] * one_minus_cos_theta + unit[1] * sin_theta
    yx = unit[1] * unit[0] * one_minus_cos_theta + unit[2] * sin_theta
    yy = cos_theta + unit[1] * unit[1] * one_minus_cos_theta
    yz = unit[1] * unit[2] * one_minus_cos_theta - unit[0] * sin_theta
    zx = unit[2] * unit[0] * one_minus_cos_theta - unit[1] * sin_theta
    zy = unit[2] * unit[1] * one_minus_cos_theta + unit[0] * sin_theta
    zz = cos_theta + unit[2] * unit[2] * one_minus_cos_theta
    rot = np.array([[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]])
    return rot

def rotate(vec, R):
    """Rotates a 1D vector trough a left multiplication.

    :param ndarray vec: Vector to rotate.
    :param ndarray R: Rotation matrix.
    :return: ndarray: Rotated vector.

    """
    return np.dot(R, vec)

def rotate_right_multiply(vec, R):
    """Rotates a 1D vector trough a right multiplication.

    :param ndarray vec: Vector to rotate.
    :param ndarray R: Rotation matrix.
    :return: ndarray: Rotated vector.

    """
    return np.dot(vec, R)

def scalar_projection(vec1, vec2):
    """Scalar projection of a vector (vec1) onto another vector (vec2).

    Ref: Linear Algebra with Applications, international edition, 2014 - p 248.

    :param ndarray vec1: Vector to project.
    :param ndarray vec2: Vector to project onto.
    :return float: The scalar projection.

    """
    inner_product = np.dot(vec1, vec2)
    norm = np.linalg.norm(vec2)
    return inner_product / norm


def vector_projection(vec1, vec2):
    """Vector projection of a vector (vec1) onto another vector (vec2).

    Ref: Linear Algebra with Applications, international edition, 2014 - p 248.

    :param ndarray vec1: Vector to project.
    :param ndarray vec2: Vector to project onto.
    :return ndarray: Projected vector.

    """
    return scalar_projection(vec1, vec2) * normalize(vec2)

def vector_projection_on_subspace(vec1, *vectors, atol=1e-6):
    """Vector projection of a vector (vec1) onto another subspace spanned by *vectors.

    Ref: My imagination.

    The projection is a sum of projections onto all subspace vectors

    :param np.ndarray vec1: Vector to project.
    :param vectors: Vector(s) that spans the subspace.
    :return np.ndarray: Projected vector.

    """
    assert all(np.isclose(np.dot(a, b), 0, atol=atol) for na, a in enumerate(vectors) for nb, b in enumerate(vectors) if na != nb), "vectors must be orthogonal"
    projection = np.zeros((len(vec1)))
    for vector in vectors:
        projection += scalar_projection(vec1, vector) * normalize(vector)
    return projection

def normalize(vec):
    """Normalize a vector.

    :param ndarray vec: Vector to normalize.
    :return ndarray: The normalized vector.

    """
    norm = np.linalg.norm(vec)
    return vec / norm

def create_orthonormal_plane(vec1, vec2):
    """Create a set of 2 orthonormal vectors spanning the plane defined by vec1 and vec2."""
    raise NotImplementedError
    # make them ortogonal
    cp = np.vector(vec1, vec2)

def rotate_z(angle):
    theta = np.radians(90)
    cos, sin = np.cos(math.radians(angle)), np.sin(math.radians(angle))
    return np.matrix([[cos, -sin, 0],
                      [sin, cos, 0],
                      [0, 0, 1]])
