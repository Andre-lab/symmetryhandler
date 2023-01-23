#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CoordinateFrame class
@Author: Mads Jeppesen
@Date: 9/21/22
"""
import math
import numpy as np
from symmetryhandler.mathfunctions import rotate_right_multiply, rotate
from symmetryhandler.exceptions import UndefinedParameters

class CoordinateFrame:
    """A coordinate system containing 4 virtual residues. It is used when modelling symmetry in Rosetta."""

    def __init__(self, name, x=None, y=None, z=None, orig=None, generate_z_from_xy=False):
        self.name = name

        if type(x) is list:
            self._vrt_x = np.array(x)
        else:
            self._vrt_x = x

        if type(y) is list:
            self._vrt_y = np.array(y)
        else:
            self._vrt_y = y

        if generate_z_from_xy:
            self.generate_z()
        else:
            if type(z) is list:
                self._vrt_z = np.array(z)
            else:
                self._vrt_z = z

        if type(orig) is list:
            self._vrt_orig = np.array(orig)
        else:
            self._vrt_orig = orig

    def generate_z(self):
        self.vrt_z = np.cross(self.vrt_x, self.vrt_y)

    def rotate_right_multiply(self, R, return_self=False):
        """Rotates the coordinate frame with the rotation matrix R"""
        self._vrt_orig = rotate_right_multiply(self._vrt_orig, R)
        self._vrt_x = rotate_right_multiply(self._vrt_x, R)
        self._vrt_y = rotate_right_multiply(self._vrt_y, R)
        self._vrt_z = rotate_right_multiply(self._vrt_z, R)
        if return_self:
            return self

    def rotate(self, R, return_self=False):
        """Rotates the coordinate frame with the rotation matrix R"""
        self._vrt_orig = rotate(self._vrt_orig, R)
        self._vrt_x = rotate(self._vrt_x, R)
        self._vrt_y = rotate(self._vrt_y, R)
        self._vrt_z = rotate(self._vrt_z, R)
        if return_self:
            return self

    @property
    def vrt_x(self):
        return self._vrt_x

    @vrt_x.setter
    def vrt_x(self, vrt_x):
        self._vrt_x = vrt_x

    @property
    def vrt_y(self):
        return self._vrt_y

    @vrt_y.setter
    def vrt_y(self, vrt_y):
        self._vrt_y = vrt_y

    @property
    def vrt_x(self):
        return self._vrt_x

    @vrt_x.setter
    def vrt_x(self, vrt_x):
        self._vrt_x = vrt_x

    @property
    def vrt_z(self):
        return self._vrt_z

    @vrt_z.setter
    def vrt_z(self, vrt_z):
        self._vrt_z = vrt_z

    @property
    def vrt_orig(self):
        return self._vrt_orig

    @vrt_orig.setter
    def vrt_orig(self, vrt_orig):
        self._vrt_orig = vrt_orig

    def __str__(self):
        str_rosetta_style = "xyz " + self.name
        for vrt in (self._vrt_x, self._vrt_y, self._vrt_orig):
            if vrt is None:  # do not print anything if these are not defined
                raise UndefinedParameters("parameters of the CoodinateFrame have not been defined")
            str_rosetta_style += " "
            for i, index in enumerate(vrt, 1):
                str_rosetta_style += "{:.6f}".format(index)
                if i != 3:
                    str_rosetta_style += ","
        return str_rosetta_style

    def identical_x(self, coordinateframe: 'CoordinateFrame'):
        """Returns True if the x axis of the coordinateframe is the same between this and the parsed CoordinateFrame."""
        return np.all(np.isclose(coordinateframe._vrt_x, self._vrt_x))

    def identical_y(self, coordinateframe: 'CoordinateFrame'):
        """Returns True if the y axis of the coordinateframe is the same between this and the parsed CoordinateFrame."""
        return np.all(np.isclose(coordinateframe._vrt_y, self._vrt_y))

    def identical_z(self, coordinateframe: 'CoordinateFrame'):
        """Returns True if the z axis of the coordinateframe is the same between this and the parsed CoordinateFrame."""
        return np.all(np.isclose(coordinateframe._vrt_z, self._vrt_z))

    def identical_coordinate_axes(self, coordinateframe: 'CoordinateFrame'):
        """Returns True if the coordinate axes (x, y, z) are the same between this and the parsed CoordinateFrame."""
        same = True
        same *= np.all(np.isclose(coordinateframe._vrt_x, self._vrt_x))
        same *= np.all(np.isclose(coordinateframe._vrt_y, self._vrt_y))
        same *= np.all(np.isclose(coordinateframe._vrt_z, self._vrt_z))
        return same


