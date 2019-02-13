"""
Utility functions for geometry.
Created on Feb 14, 2019

@author: Mahyar Osanlouy
"""


import math

from scaffoldmaker.utils.vector import *

def rotationMatrix(angle, direction, point=None):
    '''
    :param angle: Theta angle
    :param direction: Axis
    :param point: If roation is not around the origin, the input the point
    :return: The rotation matrix
    '''
    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = normalise(direction[:3])
    # rotation matrix around unit vector
    R = diagonal([cosa, cosa, cosa])
    R += outerProduct(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += ([[0.0, -direction[2], direction[1]],
         [direction[2], 0.0, -direction[0]],
         [-direction[1], direction[0], 0.0]])
    M = identity(4)
    M[:3, :3] = R
    if point is not None:  # rotation not around origin
        point = (point[:3])
        M[:3, 3] = point - dotproduct(R, point)
    return M