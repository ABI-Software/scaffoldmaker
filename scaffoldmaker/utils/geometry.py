'''
Utility functions for geometry.
'''

from __future__ import division
import math

def getApproximateEllipsePerimeter(a, b):
    '''
    Get perimeter of ellipse using Ramanujan II approximation.
    :param a: Major axis length.
    :param b: Minor axis length.
    :return: Perimeter length.
    '''
    h = ((a-b)/(a+b))**2
    return math.pi*(a + b)*(1.0 + 3.0*h/(10.0 + math.sqrt(4.0 - 3.0*h)))

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
    #print('updateEllipseAngleByArcLength', a, b, 'inAngleRadians', inAngleRadians, ', arcLength', arcLength)
    while math.fabs(arcLength - lengthMoved) > lengthTol:
        t = ( -a*math.sin(angle), b*math.cos(angle) )
        dlength_dangle = math.sqrt(t[0]*t[0] + t[1]*t[1])
        angle += (arcLength - lengthMoved)/dlength_dangle
        lengthMoved = getEllipseArcLength(a, b, inAngleRadians, angle)
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
    :param wallThickness: Constant wall thickness around.
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
