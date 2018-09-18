'''
Interpolation functions shared by mesh generators.
Created on Nov 15, 2017

@author: Richard Christie
'''

from __future__ import division
import math
import scaffoldmaker.utils.vector as vector

gaussXi3 = ( (-math.sqrt(0.6)+1.0)/2.0, 0.5, (+math.sqrt(0.6)+1.0)/2.0 )
gaussWt3 = ( 5.0/18.0, 4.0/9.0, 5.0/18.0 )

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

def computeCubicHermiteArcLength(v1, d1, v2, d2, rescaleDerivatives):
    """
    Compute arc length between v1 and v2, scaling unit d1 and d2.
    Iterative; not optimised.
    :param d1: Initial derivative at v1.
    :param d2: Initial derivative at v2.
    :param rescaleDerivatives: If True, rescale initial d1 and d2 to |v2 - v|
    :return: Arc length.
    """
    if rescaleDerivatives:
        lastArcLength = math.sqrt(sum((v2[i] - v1[i])*(v2[i] - v1[i]) for i in range(len(v1))))
    else:
        lastArcLength = getCubicHermiteArcLength(v1, d1, v2, d2)
    d1 = vector.normalise(d1)
    d2 = vector.normalise(d2)
    tol = 1.0E-6
    for iters in range(100):
        #print('iter',iters,'=',lastArcLength)
        d1s = [lastArcLength*d for d in d1]
        d2s = [lastArcLength*d for d in d2]
        arcLength = getCubicHermiteArcLength(v1, d1s, v2, d2s)
        if iters > 9:
            arcLength = 0.8*arcLength + 0.2*lastArcLength
        if math.fabs(arcLength - lastArcLength) < tol*arcLength:
            #print('computeCubicHermiteArcLength converged at iter',iters,'=',arcLength,', closeness', math.fabs(arcLength - lastArcLength))
            return arcLength
        lastArcLength = arcLength
    print('computeCubicHermiteArcLength Max iters reached:',iters,'=',arcLength,', closeness', math.fabs(arcLength - lastArcLength))
    return arcLength

def getCubicHermiteArcLength(v1, d1, v2, d2):
    '''
    Note this is approximate.
    :return: Arc length of cubic curve using 3 point Gaussian quadrature.
    '''
    arcLength = 0.0
    for i in range(3):
        dm = interpolateCubicHermiteDerivative(v1, d1, v2, d2, gaussXi3[i])
        arcLength += gaussWt3[i]*math.sqrt(sum(d*d for d in dm))
    return arcLength

def getCubicHermiteArcLengthToXi(v1, d1, v2, d2, xi):
    '''
    Note this is approximate.
    :return: Arc length of cubic curve up to given xi coordinate.
    '''
    d1m = [ d*xi for d in d1 ]
    v2m = interpolateCubicHermite(v1, d1, v2, d2, xi)
    d2m = interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi)
    d2m = [ d*xi for d in d2m ]
    return getCubicHermiteArcLength(v1, d1m, v2m, d2m)

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

def getLagrangeHermiteStartDerivative(v1, v2, d2):
    """
    Computes the derivative at v2 from quadratic Hermite-Lagrange interpolation
    from v1 to v2, d2.
    :return: d1 (dx/dxi) at v1
    """
    xi = 0.0
    #phi1 = 1 - 2.0*xi + xi*xi
    #phi2 = 2.0*xi - xi*xi
    #phi3 = -xi + xi*xi
    dphi1 = -2.0 + 2.0*xi
    dphi2 = 2.0 - 2.0*xi
    dphi3 = -1.0 + 2.0*xi
    d1 = [ (v1[c]*dphi1 + v2[c]*dphi2 + d2[c]*dphi3) for c in range(3) ]
    return d1

def getHermiteLagrangeEndDerivative(v1, d1, v2):
    """
    Computes the derivative at v2 from quadratic Hermite-Lagrange interpolation
    from v1, d1.
    :return: d2 (dx/dxi) at v2
    """
    xi = 1.0
    #phi1 = 1 - xi*xi
    #phi2 = xi - xi*xi
    #phi3 = xi*xi
    dphi1 = -2.0*xi
    dphi2 = 1 - 2.0*xi
    dphi3 = 2.0*xi
    d2 = [ (v1[c]*dphi1 + d1[c]*dphi2 + v2[c]*dphi3) for c in range(3) ]
    return d2

def sampleCubicHermiteCurves(nx, nd1, nd2, elementsCountOut,
    addLengthStart = 0.0, addLengthEnd = 0.0,
    lengthFractionStart = 1.0, lengthFractionEnd = 1.0,
    elementLengthStartEndRatio = 1.0):
    """
    Get even-spaced points through cubic Hermite nodes nx with derivatives nd1
    in line and nd2 across. Derivatives nd1 are rescaled to give arc length
    scaling across each input element.
    :param nx: Coordinates of nodes along curves.
    :param nd1: Derivatives of nodes along curves.
    :param nd2: List of cross derivatives to interpolate, or None to ignore.
    :param addLengthStart, addLengthEnd: Extra length to add to start and end elements.
    :param lengthFractionStart, lengthFractionEnd: Fraction of mid element length for
        start and end elements. Can use in addition to AddLengths: If LengthFraction
        is 0.5 and AddLength is derivative/2.0 can blend into known derivative at start
        or end.
    :param elementLengthStartEndRatio: Start/end element length ratio, with lengths
        smoothly varying in between. Requires at least 2 elements. Applied in proportion
        to lengthFractionStart, lengthFractionEnd.
    :return: px[], pd1[], pd2[]
    """
    elementsCountIn = len(nx) - 1
    lengths = [ 0.0 ]
    nd1a = []
    nd1b = []
    length = 0.0
    for e in range(elementsCountIn):
        arcLength = computeCubicHermiteArcLength(nx[e], nd1[e], nx[e + 1], nd1[e + 1], rescaleDerivatives = True)
        length += arcLength
        lengths.append(length)
        nd1a.append(vector.setMagnitude(nd1[e], arcLength))
        nd1b.append(vector.setMagnitude(nd1[e + 1], arcLength))
    proportionEnd = 2.0/(elementLengthStartEndRatio + 1)
    proportionStart = elementLengthStartEndRatio*proportionEnd
    elementLengthMid = (length - addLengthStart - addLengthEnd) / \
        (elementsCountOut - 2.0 + proportionStart*lengthFractionStart + proportionEnd*lengthFractionEnd)
    elementLengthProportionStart = proportionStart*lengthFractionStart*elementLengthMid
    elementLengthProportionEnd = proportionEnd*lengthFractionEnd*elementLengthMid
    # get smoothly varying element lengths, not accounting for start and end
    elementLengths = []
    #cumulativeLengths = [ 0.0 ]
    for eOut in range(elementsCountOut):
        xi = eOut/(elementsCountOut - 1)
        elementLengths.append(((1.0 - xi)*proportionStart + xi*proportionEnd)*elementLengthMid)
        #cummulativeLengths.append(cummulativeLengths[eOut] + elementLengths[eOut])
    # get middle derivative magnitudes
    nodeDerivativeMagnitudes = [ None ]*(elementsCountOut + 1)  # start and end determined below
    for n in range(1, elementsCountOut):
        nodeDerivativeMagnitudes[n] = 0.5*(elementLengths[n - 1] + elementLengths[n])
    # fix end lengths:
    elementLengths[ 0] = addLengthStart + elementLengthProportionStart
    elementLengths[-1] = addLengthEnd + elementLengthProportionEnd
    #print('\nsampleCubicHermiteCurves:')
    #print('  elementLengths', elementLengths, 'addLengthStart', addLengthStart, 'addLengthEnd', addLengthEnd)
    #print('  sum lengths', sum(elementLengths), 'vs. length', length, 'diff', sum(elementLengths) - length)
    # set end derivatives:
    if elementsCountOut == 1:
        nodeDerivativeMagnitudes[0] = nodeDerivativeMagnitudes[1] = elementLengths[0]
    else:
        nodeDerivativeMagnitudes[0] = elementLengths[ 0]*2.0 - nodeDerivativeMagnitudes[ 1]
        nodeDerivativeMagnitudes[-1]   = elementLengths[-1]*2.0 - nodeDerivativeMagnitudes[-2]

    px = []
    pd1 = []
    if nd2 is None:
        pd2 = None
    else:
        pd2 = []
    distance = 0.0
    e = 0
    for eOut in range(elementsCountOut):
        while e < elementsCountIn:
            if distance < lengths[e + 1]:
                xi = (distance - lengths[e])/(lengths[e + 1] - lengths[e])
                px.append(list(interpolateCubicHermite(nx[e], nd1a[e], nx[e + 1], nd1b[e], xi)))
                pd1.append(vector.setMagnitude(interpolateCubicHermiteDerivative(nx[e], nd1a[e], nx[e + 1], nd1b[e], xi), nodeDerivativeMagnitudes[eOut]))
                if pd2 is not None:
                    pd2.append([ (nd2[e][c]*(1.0 - xi) + nd2[e + 1][c]*xi) for c in range(3) ])
                break
            e += 1
        distance += elementLengths[eOut]
    px.append(nx[-1])
    pd1.append(vector.setMagnitude(nd1[-1], nodeDerivativeMagnitudes[-1]))
    if pd2 is not None:
        pd2.append(nd2[-1])
    return px, pd1, pd2

def getCubicHermiteCurvesPointAtArcDistance(nx, nd, arcDistance):
    """
    Get the coordinates, derivatives at distance along cubic Hermite curves. Note this is approximate.
    :param nx: Coordinates of nodes along curves.
    :param nd: Derivatives of nodes along curves.
    :param distance: Distance along curves.
    :return: coordinates, derivatives; clamped to first or last nx if distance is beyond curves
    """
    elementsCount = len(nx) - 1
    assert elementsCount > 0, 'getCubicHermiteCurvesPointAtArcDistance.  Invalid number of points'
    if arcDistance < 0.0:
        return nx[0], nd[0]
    length = 0.0
    xiDelta = 1.0E-6
    xiTol = 1.0E-6
    #print('elementsCount',elementsCount,'arcDistance',arcDistance)
    for e in range(elementsCount):
        partDistance = arcDistance - length
        v1 = nx[e]
        d1 = nd[e]
        v2 = nx[e + 1]
        d2 = nd[e + 1]
        arcLength = getCubicHermiteArcLength(v1, d1, v2, d2)
        #print('e',e,'partDistance',partDistance,'arcLength',arcLength)
        if partDistance <= arcLength:
            xiLast = 100.0
            xi = partDistance/arcLength
            for iter in range(100):
                xiLast = xi
                dist = getCubicHermiteArcLengthToXi(v1, d1, v2, d2, xi)
                #print('iter',iter,'xi',xi,'--> dist',dist)
                distp = getCubicHermiteArcLengthToXi(v1, d1, v2, d2, xi + xiDelta)
                distm = getCubicHermiteArcLengthToXi(v1, d1, v2, d2, xi - xiDelta)
                dxi_ddist = (2.0*xiDelta)/(distp - distm)
                xi -= dxi_ddist*(dist - partDistance)
                if math.fabs(xi - xiLast) <= xiTol:
                    #print('converged xi',xi)
                    return list(interpolateCubicHermite(v1, d1, v2, d2, xi)), list(interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi))
            print('getCubicHermiteCurvesPointAtArcDistance Max iters reached:',iter,': e', e, ', xi',xi,', closeness', math.fabs(dist - partDistance))
            return v2
        length += arcLength
    return v2, d2

def getDoubleCubicHermiteCurvesMidDerivative(ax, ad1, mx, bx, bd1):
    """
    Get derivative at centre of two cubic curves.
    :return: Derivative at mx to balance ax, ad1 with bx, bd1.
    """
    md1 = [ (bx[c] - ax[c]) for c in range(3) ]
    arcLengtha = computeCubicHermiteArcLength(ax, ad1, mx, md1, rescaleDerivatives = True)
    arcLengthb = computeCubicHermiteArcLength(mx, md1, bx, bd1, rescaleDerivatives = True)
    maga = vector.magnitude(ad1)
    magb = vector.magnitude(bd1)
    magm = arcLengtha + arcLengthb - 0.5*(maga + magb)
    return vector.setMagnitude(md1, magm)
