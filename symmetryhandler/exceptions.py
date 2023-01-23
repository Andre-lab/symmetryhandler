#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exception classes
@Author: Mads Jeppesen
@Date: 9/21/22
"""

class UndefinedParameters(BaseException):
    """An exception that reports that the user has undefined parameters. Used in the CoordinateFrame class."""
    def __init__(self, message=None):
        self.message = message
        print(message)
