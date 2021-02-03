#!/usr/bin/env python3
# coding=utf-8

import math
import numpy as np

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
    rot = [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]]
    return rot


def rotate(vec, R):
    """Rotates a 1D vector.

    :param ndarray vec: Vector to rotate.
    :param ndarray R: Rotation matrix.
    :return: ndarray: Rotated vector.

    """
    return np.dot(vec, R)
