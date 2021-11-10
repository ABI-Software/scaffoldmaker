'''
Utility functions for geometry.
'''

from __future__ import division

import copy
import math

from scaffoldmaker.utils import vector


def getApproximateEllipsePerimeter(a, b):
    '''
    Get perimeter of ellipse using Ramanujan II approximation.
    :param a: Major axis length.
    :param b: Minor axis length.
    :return: Perimeter length.
    '''
    h = ((a-b)/(a+b))**2
    return math.pi*(a + b)*(1.0 + 3.0*h/(10.0 + math.sqrt(4.0 - 3.0*h)))

def getEllipseAngleFromVector(a, b, x, y):
    '''
    Given an ellipse with major and minor axes a and b cented at (0.0, 0.0),
    get angle in radians in direction of vector (dx, dy)
    :param a: Major axis length.
    :param b: Minor axis length.
    :param x, y: Non-zero vector direction
    :return: Angle in radians.
    '''
    # Compute theta such that
    # r*x = a*cos_theta
    # r*y = b*sin_theta
    # hence:
    # y/x = b/a * tan(theta)
    # (a*y)/(b*x) = tan(theta)
    return math.atan2(a*y, b*x)

def getEllipseArcLength(a, b, angle1Radians, angle2Radians):
    '''
    Calculates perimeter distance between two angles by summing line segments at regular angles.
    :param a: Major axis length (On x, 0 / PI).
    :param b: Minor axis length.(On y, PI/2, 3PI/2).
    :param angle1Radians: First angle anticlockwise from major axis.
    :param angle2Radians: Second angle anticlockwise from major axis.
    :return: Perimeter length, positive if anticlockwise, otherwise negative.
    '''
    angle1 = min(angle1Radians, angle2Radians)
    angle2 = max(angle1Radians, angle2Radians)
    # Max 100 segments around ellipse
    segmentCount = int(math.ceil(50*(angle2-angle1)/math.pi))
    length = 0.0
    for i in range(segmentCount + 1):
        r = i/segmentCount
        angle = r*angle1 + (1.0 - r)*angle2
        x = ( a*math.cos(angle), b*math.sin(angle) )
        if i > 0:
            delta = math.sqrt((x[0] - lastX[0])*(x[0] - lastX[0]) + (x[1] - lastX[1])*(x[1] - lastX[1]))
            length += delta
        lastX = x
    if angle1Radians < angle2Radians:
        return length
    else:
        return -length

def updateEllipseAngleByArcLength(a, b, inAngleRadians, arcLength):
    '''
    Update angle around ellipse to subtend arcLength around the perimeter.
    Iterates using Newton's method.
    :param inAngleRadians: Initial angle anticlockwise from major axis.
    :param arcLength: Arc length to traverse. Positive=anticlockwise, negative=clockwise.
    :param a: Major axis length (On x, 0 / PI).
    :param b: Minor axis length.(On y, PI/2, 3PI/2).
    :return: New angle, in radians.
    '''
    angle = inAngleRadians
    lengthMoved = 0.0
    lengthTol = (a + b)*1.0E-4  # broader tolerance due to reliance on inexact getEllipseArcLength()
    counter=0
    #print('updateEllipseAngleByArcLength', a, b, 'inAngleRadians', inAngleRadians, ', arcLength', arcLength)
    while math.fabs(arcLength - lengthMoved) > lengthTol:
        angleOld = angle
        t = ( -a*math.sin(angle), b*math.cos(angle) )
        dlength_dangle = math.sqrt(t[0]*t[0] + t[1]*t[1])
        angle += (arcLength - lengthMoved)/dlength_dangle
        if counter >= 100:
            angle=(angle+angleOld)/2
        lengthMoved = getEllipseArcLength(a, b, inAngleRadians, angle)
        counter+=1
        #print('lengthMoved', lengthMoved)
    #print('updateEllipseAngleByArcLength a', a, 'b', b, ', angle', inAngleRadians, ', arcLength', arcLength, ' -> ', angle)
    return angle

def getEllipseRadiansToX(ax, bx, dx, initialTheta):
    '''
    Iteratively compute theta in ax*cos(theta) + bx*sin(theta) = dx
    from initialTheta. Uses Newton's method.
    :param ax: x component of major axis.
    :param bx: x component of minor axis.
    :param dx: x distance from centre of ellipse.
    :param initialTheta: Initial value of theta needed since multi-valued.
    :return: theta in radians
    '''
    theta = initialTheta
    iters = 0
    fTol = math.sqrt(ax*ax + bx*bx)*1.0E-10
    while True:
        cosAngle = math.cos(theta)
        sinAngle = math.sin(theta)
        f = ax*cosAngle + bx*sinAngle - dx
        if math.fabs(f) < fTol:
            break;
        df = -ax*sinAngle + bx*cosAngle
        #print(iters, '. theta', theta, 'f', f,'df',df,'-->',theta - f/df)
        theta -= f/df
        iters += 1
        if iters == 100:
            print('getEllipseRadiansToX: did not converge!')
            break
    return theta

def createCirclePoints(cx, axis1, axis2, elementsCountAround, startRadians = 0.0):
    '''
    Create circular ring of points centred at cx, from axis1 around through axis2.
    Assumes axis1 and axis2 are orthogonal and equal magnitude.
    Dimension 3 only.
    :param cx: centre
    :param axis1:  Vector from cx to inside at zero angle
    :param axis2:  Vector from cx to inside at 90 degree angle.
    :param elementsCountAround: Number of elements around.
    :return: lists px, pd1
    '''
    px = []
    pd1 = []
    radiansPerElementAround = 2.0*math.pi/elementsCountAround
    radiansAround = startRadians
    for n in range(elementsCountAround):
        cosRadiansAround = math.cos(radiansAround)
        sinRadiansAround = math.sin(radiansAround)
        px.append([ (cx[c] + cosRadiansAround*axis1[c] + sinRadiansAround*axis2[c]) for c in range(3) ])
        pd1.append([ radiansPerElementAround*(-sinRadiansAround*axis1[c] + cosRadiansAround*axis2[c]) for c in range(3) ])
        radiansAround += radiansPerElementAround
    return px, pd1

def createEllipsoidPoints(centre, poleAxis, sideAxis, elementsCountAround, elementsCountUp, height):
    '''
    Generate a set of points and derivatives for circle of revolution of an ellipse
    starting at pole poleAxis from centre.
    :param centre: Centre of full ellipsoid.
    :param poleAxis: Vector in direction of starting pole, magnitude is ellipse axis length.
    :param sideAxis: Vector normal to poleAxis, magnitude is ellipse side axis length.
    :param height: Height of arc of ellipsoid from starting pole along poleAxis.
    :return: Lists nx, nd1, nd2. Ordered fastest around, starting at pole. Suitable for passing to TrackSurface.
    '''
    nx  = []
    nd1 = []
    nd2 = []
    magPoleAxis = vector.magnitude(poleAxis)
    magSideAxis = vector.magnitude(sideAxis)
    unitPoleAxis = vector.normalise(poleAxis)
    unitSideAxis1 = vector.normalise(sideAxis)
    unitSideAxis2 = vector.normalise(vector.crossproduct3(sideAxis, poleAxis))
    useHeight = min(max(0.0, height), 2.0*magPoleAxis)
    totalRadiansUp = getEllipseRadiansToX(magPoleAxis, 0.0, magPoleAxis - useHeight, initialTheta = 0.5*math.pi*useHeight/magPoleAxis)
    radiansUp = 0.0
    lengthUp = getEllipseArcLength(magPoleAxis, magSideAxis, radiansUp, totalRadiansUp)
    elementLengthUp = lengthUp/elementsCountUp
    radiansPerElementAround = 2.0*math.pi/elementsCountAround
    for n2 in range(elementsCountUp + 1):
        cosRadiansUp = math.cos(radiansUp)
        sinRadiansUp = math.sin(radiansUp)
        radius = sinRadiansUp*magSideAxis
        d2r, d2z = vector.setMagnitude([ cosRadiansUp*magSideAxis, sinRadiansUp*magPoleAxis ], elementLengthUp)
        cx = [ (centre[c] + cosRadiansUp*poleAxis[c]) for c in range(3) ]
        elementLengthAround = radius*radiansPerElementAround
        radiansAround = 0.0
        for n in range(elementsCountAround):
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            nx .append([ (cx[c] + radius*(cosRadiansAround*unitSideAxis1[c] + sinRadiansAround*unitSideAxis2[c])) for c in range(3) ])
            nd1.append([ (elementLengthAround*(-sinRadiansAround*unitSideAxis1[c] + cosRadiansAround*unitSideAxis2[c])) for c in range(3) ])
            nd2.append([ (d2r*(cosRadiansAround*unitSideAxis1[c] + sinRadiansAround*unitSideAxis2[c]) - d2z*unitPoleAxis[c]) for c in range(3) ])
            radiansAround += radiansPerElementAround
        radiansUp = updateEllipseAngleByArcLength(magPoleAxis, magSideAxis, radiansUp, elementLengthUp)
    return nx, nd1, nd2

def getCircleProjectionAxes(ax, ad1, ad2, ad3, length, angle1radians, angle2radians, angle3radians = None):
    '''
    Project coordinates and orthogonal unit axes ax, ad1, ad2, ad3 by length in
    direction of ad3, rotated towards ad1 by angle1radians and ad2 by
    angle2radians such that it conforms to a circular arc of the given length,
    tangential to ad3 at the start and the final bd3 at the end.
    All vectors must be length 3.
    Assumes angles 1 and 2 are not large e.g. less than 90 degrees.
    Note: not robust for all inputs.
    :param angle3radians: Optional final rotation of projection axes about bd3.
    :return: Final coordinates and orthogonal unit axes: bx, bd1, bd2, bd3
    '''
    if (math.fabs(angle1radians) < 0.000001) and (math.fabs(angle2radians) < 0.000001):
        bx = [ (ax[c] + length*ad3[c]) for c in range(3) ]
        bd1 = copy.deepcopy(ad1)
        bd2 = copy.deepcopy(ad2)
        bd3 = copy.deepcopy(ad3)
    else:
        cosAngle1 = math.cos(angle1radians)
        sinAngle1 = math.sin(angle1radians)
        cosAngle2 = math.cos(angle2radians)
        sinAngle2 = math.sin(angle2radians)
        f1 = sinAngle1*cosAngle2
        f2 = cosAngle1*sinAngle2
        f3 = cosAngle1*cosAngle2
        angleAroundRadians = math.atan2(f2, f1)
        fh = math.sqrt(f1*f1 + f2*f2)
        arcAngleRadians = 0.5*math.pi - math.atan2(f3, fh)
        arcRadius = length/arcAngleRadians
        br = arcRadius*(1.0 - math.cos(arcAngleRadians))
        w1 = br*math.cos(angleAroundRadians)  # f1/fh
        w2 = br*math.sin(angleAroundRadians)  # f2/fh
        w3 = arcRadius*math.sin(arcAngleRadians)
        bx = [ (ax[c] + w1*ad1[c] + w2*ad2[c] + w3*ad3[c]) for c in range(3) ]
        bd3 = vector.normalise([ (f1*ad1[c] + f2*ad2[c] + f3*ad3[c]) for c in range(3) ])
        bd1 = vector.normalise(vector.crossproduct3(ad2, bd3))
        bd2 = vector.crossproduct3(bd3, bd1)
    if angle3radians:
        cosAngle3 = math.cos(angle3radians)
        sinAngle3 = math.sin(angle3radians)
        bd1, bd2 = [ (cosAngle3*bd1[c] + sinAngle3*bd2[c]) for c in range(3) ], [ (cosAngle3*bd2[c] - sinAngle3*bd1[c]) for c in range(3) ]
    return bx, bd1, bd2, bd3

def getSurfaceProjectionAxes(ax, ad1, ad2, ad3, angle1radians, angle2radians, length):
    '''
    Project coordinates and orthogonal unit axes ax, ad1, ad2, ad3 in direction
    of ad3, rotated towards ad1 by angle1radians and ad2 by angle2radians by
    simple rotation.
    All vectors are/must be length 3.
    Assumes angles are not large e.g. less than 90 degrees.
    Note: not robust for all inputs.
    :return: Final coordinates and orthogonal unit axes: bx, bd1, bd2, bd3
    '''
    cosAngle1 = math.cos(angle1radians)
    sinAngle1 = math.sin(angle1radians)
    cosAngle2 = math.cos(angle2radians)
    sinAngle2 = math.sin(angle2radians)
    f1 = sinAngle1*cosAngle2
    f2 = cosAngle1*sinAngle2
    f3 = cosAngle1*cosAngle2
    bd3 = [ (f1*ad1[c] + f2*ad2[c] + f3*ad3[c]) for c in range (3) ]
    bx = [ (ax[c] + length*bd3[c]) for c in range(3) ]
    bd1 = vector.crossproduct3(ad2, bd3)
    bd2 = vector.crossproduct3(bd3, bd1)
    return bx, bd1, bd2, bd3

def createEllipsePoints(cx, radian, axis1, axis2, elementsCountAround, startRadians = 0.0):
    '''
    Create ellipse points centred at cx, from axis1 around through axis2.
    Assumes axis1 and axis2 are orthogonal.
    :param cx: centre
    :param radian: Part of ellipse to be created based on radian (will be 2*math.pi for a complete ellipse).
    :param axis1: Vector from cx to inside at zero angle
    :param axis2: Vector from cx to inside at 90 degree angle.
    :param elementsCountAround: Number of elements around.
    :param startRadians: Angle from axis1 to start creating the ellipse.
    :return: lists px, pd1
    '''
    px = []
    pd1 = []
    magAxis1 = vector.magnitude(axis1)
    magAxis2 = vector.magnitude(axis2)
    totalEllipsePerimeter = getApproximateEllipsePerimeter(magAxis1, magAxis2)
    partOfEllipsePerimeter = radian * totalEllipsePerimeter / (2 * math.pi)
    elementLength = partOfEllipsePerimeter / elementsCountAround
    if radian != 2 * math.pi:
        elementsCountAround = elementsCountAround + 1

    unitSideAxis1 = vector.normalise(axis1)
    unitSideAxis2 = vector.normalise(axis2)
    for n in range(elementsCountAround):
        angle = startRadians
        arcLength = n * elementLength
        newAngle = updateEllipseAngleByArcLength(magAxis1, magAxis2, angle, arcLength)
        cosRadiansAround = math.cos(newAngle)
        sinRadiansAround = math.sin(newAngle)
        x = [
            cx[0] + cosRadiansAround * axis1[0] - sinRadiansAround * axis2[0],
            cx[1] + cosRadiansAround * axis1[1] + sinRadiansAround * axis2[1],
            cx[2] + cosRadiansAround * axis1[2] + sinRadiansAround * axis2[2]
        ]
        px.append(x)
        rd1 = [magAxis1 * (-sinRadiansAround * unitSideAxis1[c]) + magAxis2 * (cosRadiansAround * unitSideAxis2[c]) for c in range(3)]
        rd1Norm = vector.normalise(rd1)
        pd1.append([elementLength * (rd1Norm[c])for c in range(3)])

    return px, pd1

