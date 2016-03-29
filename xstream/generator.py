#!/usr/bin/env python
# -*- coding: utf-8 -*-

# #########################################################################
# Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.         #
#                                                                         #
# Copyright 2016. UChicago Argonne, LLC. This software was produced       #
# under U.S. Government contract DE-AC02-06CH11357 for Argonne National   #
# Laboratory (ANL), which is operated by UChicago Argonne, LLC for the    #
# U.S. Department of Energy. The U.S. Government has rights to use,       #
# reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR    #
# UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR        #
# ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is     #
# modified to produce derivative works, such modified software should     #
# be clearly marked, so as not to confuse it with the version available   #
# from ANL.                                                               #
#                                                                         #
# Additionally, redistribution and use in source and binary forms, with   #
# or without modification, are permitted provided that the following      #
# conditions are met:                                                     #
#                                                                         #
#     * Redistributions of source code must retain the above copyright    #
#       notice, this list of conditions and the following disclaimer.     #
#                                                                         #
#     * Redistributions in binary form must reproduce the above copyright #
#       notice, this list of conditions and the following disclaimer in   #
#       the documentation and/or other materials provided with the        #
#       distribution.                                                     #
#                                                                         #
#     * Neither the name of UChicago Argonne, LLC, Argonne National       #
#       Laboratory, ANL, the U.S. Government, nor the names of its        #
#       contributors may be used to endorse or promote products derived   #
#       from this software without specific prior written permission.     #
#                                                                         #
# THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS     #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT       #
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS       #
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago     #
# Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,        #
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,    #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;        #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER        #
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      #
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN       #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE         #
# POSSIBILITY OF SUCH DAMAGE.                                             #
# #########################################################################

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


__author__ = "Doga Gursoy"
__contact__ = "dgursoy@aps.anl.gov"
__copyright__ = "Copyright (c) 2016, UChicago Argonne, LLC."
__license__ = "BSD-3"
__version__ = "0.1.0"
__status__ = "Development"
__docformat__ = "restructuredtext en"
__all__ = ['project',
           'Point',
           'Line',
           'random_point',
           'plot_points',
           'show']


import logging
import numpy as np
import math
import random
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def project(x0, y0, x1, y1, obj):
    """Project single x-ray beam.
    """

    x0, y0, x1, y1 = float(x0), float(y0), float(x1), float(y1)
    sx, sy = obj.shape

    # grid frame (gx, gy)
    gx = np.arange(0, sx + 1)
    gy = np.arange(0, sy + 1)

    # avoid upper-right boundary errors
    if (x1 - x0) == 0:
        x0 += 1e-6
    if (y1 - y0) == 0:
        y0 += 1e-6

    # vector lengths (ax, ay)
    ax = (gx - x0) / (x1 - x0)
    ay = (gy - y0) / (y1 - y0)

    # edges of alpha (a0, a1)
    ax0 = min(ax[0], ax[-1])
    ax1 = max(ax[0], ax[-1])
    ay0 = min(ay[0], ay[-1])
    ay1 = max(ay[0], ay[-1])
    a0 = max(max(ax0, ay0), 0)
    a1 = min(min(ax1, ay1), 1)

    # sorted alpha vector
    cx = (ax >= a0) & (ax <= a1)
    cy = (ay >= a0) & (ay <= a1)
    alpha = np.sort(np.r_[ax[cx], ay[cy]])

    # lengths
    xv = x0 + alpha * (x1 - x0)
    yv = y0 + alpha * (y1 - y0)
    lx = np.ediff1d(xv)
    ly = np.ediff1d(yv)
    dist = np.sqrt(lx**2 + ly**2)
    ind = dist != 0

    # indexing
    mid = alpha[:-1] + np.ediff1d(alpha) / 2.
    xm = x0 + mid * (x1 - x0)
    ym = y0 + mid * (y1 - y0)
    ix = np.floor(xm).astype('int')
    iy = np.floor(ym).astype('int')

    # projection
    return np.dot(dist[ind], obj[ix[ind], iy[ind]])


class Point(object):
    """Definition of a point in 2-D Cartesian space.
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ")"

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __add__(self, other):
        """Addition."""
        return Point(self.x + other.x, self.y + other.y)

    def __rmul__(self, c):
        """Scalar multiplication."""
        return Point(c * self.x, c * self.y)

    def is_close(self, point, epsilon=1e-6):
        """Checks if is close to a point."""
        return self.dist(point) < epsilon

    def dist(self, point):
        """Returns the distance from a point."""
        return math.sqrt(self.dist2(point))

    def dist2(self, point):
        """Returns the square of distance from a point."""
        dx = self.x - point.x
        dy = self.y - point.y
        return dx * dx + dy * dy

    def list(self):
        """Returns the point's list representation."""
        return [self.x, self.y]

    def numpy(self):
        """Returns the Numpy representation."""
        return np.array([self.x, self.y])


class Line(object):
    """Definition of a line in 2-D Cartesian space.
    """

    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

        if p1.x == p2.x:
            self.slope = None
            self.intercept = None
            self.vertical = True
        else:
            self.slope = float(p2.y - p1.y) / (p2.x - p1.x)
            self.intercept = p1.y - self.slope * p1.x
            self.vertical = False

    def __str__(self):
        if self.vertical:
            return "x = " + str(self.p1.x)
        return "y = " + str(self.slope) + "x + " + str(self.intercept)

    def __eq__(self, line):
        if self.vertical != line.vertical:
            return False

        if self.vertical:
            return self.p1.x == line.p1.x

        return self.slope == line.slope and self.intercept == line.intercept

    def at_x(self, x):
        if self.vertical:
            return None

        return Point(x, self.slope * x + self.intercept)

    def at_y(self, y):
        return Point((y - self.intercept) / self.slope, y)

    def dist(self, point):
        """Returns the distance from a point."""
        return sqrt(self.dist2(point))

    def dist2(self, point):
        """Returns the square of distance from a point."""
        numerator = float(self.p2.x - self.p1.x) * (self.p1.y - point.y) - \
            (self.p1.x - point.x) * (self.p2.y - self.p1.y)
        numerator *= numerator
        denominator = float(self.p2.x - self.p1.x) * (self.p2.x - self.p1.x) + \
            (self.p2.y - self.p1.y) * (self.p2.y - self.p1.y)
        return numerator / denominator

    def intersection(self, line):
        """Returns the intersection point with a line."""
        if line.slope == self.slope:
            return None

        if self.vertical:
            return line.at_x(self.p1.x)
        elif line.vertical:
            return self.at_x(line.p1.x)

        x = float(self.intercept - line.intercept) / (line.slope - self.slope)
        return self.at_x(x)

    def midpoint(self):
        """Returns the midpoint of two points describing a line."""
        x = float(self.p1.x + self.p2.x) / 2
        y = float(self.p1.y + self.p2.y) / 2
        return Point(x, y)


def random_point(k=1):
    """ Generates a random point in 2-D Cartesian space. """
    return Point(k * random.random(), k * random.random())


def plot_points(points, style='bo'):
    _plot_points(points, style=style)


def _plot_points(points, style):
    if not type(points) == list:
        points = [points]

    points = _points_to_numpy(points)
    plt.plot(points[:, 0], points[:, 1], style)


def plot_lines(lines, style='b-'):
    _plot_lines(lines, style=style)


def show():
    plt.show()


def _points_to_numpy(points):
    return np.array(map(lambda p: p.numpy(), points), np.float32)
