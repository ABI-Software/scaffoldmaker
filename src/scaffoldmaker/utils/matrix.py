'''
Utility functions for matrix.
'''

import math

def getRotationMatrixFromAxisAngle(rotAxis, theta):
    """
    Generate the rotation matrix for rotation about an axis.
    :param rotAxis: axis of rotation
    :param theta: angle of rotation
    :return: rotation matrix
    """
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)
    C = 1 - cosTheta
    rotMatrix = ([[rotAxis[0]*rotAxis[0]*C + cosTheta, rotAxis[0]*rotAxis[1]*C - rotAxis[2]*sinTheta, rotAxis[0]*rotAxis[2]*C + rotAxis[1]*sinTheta],
        [rotAxis[1]*rotAxis[0]*C + rotAxis[2]*sinTheta, rotAxis[1]*rotAxis[1]*C + cosTheta, rotAxis[1]*rotAxis[2]*C - rotAxis[0]*sinTheta],
        [rotAxis[2]*rotAxis[0]*C - rotAxis[1]*sinTheta, rotAxis[2]*rotAxis[1]*C + rotAxis[0]*sinTheta, rotAxis[2]*rotAxis[2]*C + cosTheta]])

    return rotMatrix

def rotateAboutZAxis(x, theta):
    """
    Rotates matrix about z-axis.
    : param x: matrix to be rotated.
    : param theta: angle of rotation.
    :return rotated matrix
    """
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)

    xRot = [x[0]*cosTheta - x[1]*sinTheta,
            x[0]*sinTheta + x[1]*cosTheta,
            x[2]]

    return xRot
