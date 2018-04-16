'''
Interpolation functions shared by mesh generators.
Created on Nov 15, 2017

@author: Richard Christie
'''

import math
import scaffoldmaker.utils.vector as vector

def interpolateCubicHermite(v1, d1, v2, d2, xi):
    """
    Return cubic Hermite interpolated value of tuples v1, d1 (end 1) to v2, d2 (end 2) for xi in [0,1]
    :return: tuple containing result
    """
    xi2 = xi*xi
    xi3 = xi2*xi
    f1 = 1.0 - 3.0*xi2 + 2.0*xi3
    f2 = xi - 2.0*xi2 + xi3
    f3 = 3.0*xi2 - 2.0*xi3
    f4 = -xi2 + xi3
    return tuple([ (f1*v1[i] + f2*d1[i] + f3*v2[i] + f4*d2[i]) for i in range(len(v1)) ])

def interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi):
    """
    Return cubic Hermite interpolated derivatives of tuples v1, d1 (end 1) to v2, d2 (end 2) for xi in [0,1]
    :return: tuple containing result
    """
    xi2 = xi*xi
    f1 = -6.0*xi + 6.0*xi2
    f2 = 1.0 - 4.0*xi + 3.0*xi2
    f3 = 6.0*xi - 6.0*xi2
    f4 = -2.0*xi + 3.0*xi2
    return tuple([ (f1*v1[i] + f2*d1[i] + f3*v2[i] + f4*d2[i]) for i in range(len(v1)) ])

def interpolateCubicHermiteSecondDerivative(v1, d1, v2, d2, xi):
    """
    Return cubic Hermite interpolated second derivatives of tuples v1, d1 (end 1) to v2, d2 (end 2) for xi in [0,1]
    :return: tuple containing result
    """
    f1 = -6.0 + 12.0*xi
    f2 = -4.0 +  6.0*xi
    f3 =  6.0 - 12.0*xi
    f4 = -2.0 +  6.0*xi
    return tuple([ (f1*v1[i] + f2*d1[i] + f3*v2[i] + f4*d2[i]) for i in range(len(v1)) ])

def computeCubicHermiteArcLength(v1, d1, v2, d2):
    """
    Compute arc length between v1 and v2, scaling unit d1 and d2.
    Iterative; not optimised.
    :param d1: Initially unit length derivative at v1.
    :param d2: Initially unit length derivative at v2.
    :return: Arc length.
    """
    lastArcLength = math.sqrt(sum((v2[i] - v1[i])*(v2[i] - v1[i]) for i in range(len(v1))))
    xiTol = 1.0E-6
    iters = 0
    while True:
        iters += 1
        #print('iter',iters,'=',lastArcLength)
        d1s = [lastArcLength*d for d in d1]
        d2s = [lastArcLength*d for d in d2]
        dm = interpolateCubicHermiteDerivative(v1, d1s, v2, d2s, 0.5)
        magdm = math.sqrt(sum(d*d for d in dm))
        arcLength = 0.5*(magdm + lastArcLength)
        if (math.fabs(arcLength - lastArcLength) < xiTol) or (iter == 100):
            #print('iter',iters,'=',arcLength,', closeness', math.fabs(arcLength - lastArcLength)/xiTol)
            return arcLength
        lastArcLength = arcLength

def getCubicHermiteCurvature(v1, d1, v2, d2, radialVector, xi):
    """
    :param radialVector: Radial direction, assumed unit normal to curve tangent at point.
    :return: Scalar curvature (1/R) of the 1-D cubic Hermite curve.
    """
    tangent = interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi)
    dTangent = interpolateCubicHermiteSecondDerivative(v1, d1, v2, d2, xi)
    #tangentVector = vector.normalise(tangent)
    #tangentCurvature = vector.dotproduct(dTangent, tangentVector)
    radialCurvature = vector.dotproduct(dTangent, radialVector)
    magTangent = vector.magnitude(tangent)
    curvature = radialCurvature/(magTangent*magTangent)
    return curvature
