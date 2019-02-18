"""
Utility functions for geometry.
Created on Feb 14, 2019
"""


import math

from scaffoldmaker.utils.vector import *

def rotationMatrix(angle, direction, fourbyfour=None, point=None):
    '''
    :param angle: Theta angle
    :param direction: Axis
    :param point: If roation is not around the origin, the input the point
    :return: The rotation matrix
    '''
    from numpy.core.numeric import asarray

    sina = math.sin(angle)
    cosa = math.cos(angle)
    # direction = unit_vector(direction[:3])
    direction = normalise(direction[:3])
    # rotation matrix around unit vector
    R = makeDiagonal([cosa, cosa, cosa])
    # c = outerProduct(direction, direction).tolist()
    R += outerProduct(direction, direction) * (1.0 - cosa)
    direction = [x * sina for x in direction]
    R += ([[0.0, -direction[2], direction[1]],
         [direction[2], 0.0, -direction[0]],
         [-direction[1], direction[0], 0.0]])

    if fourbyfour:
        M = asarray(identity(4))
        M[:3, :3] = R

    M = asarray(identity(3))
    M[:, :] = R

    if point is not None:  # rotation not around origin
        point = (point[:3])
        M[:3, 3] = point - dotproduct(R, point)

    return M.tolist()