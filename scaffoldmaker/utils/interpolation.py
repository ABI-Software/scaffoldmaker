'''
Interpolation functions shared by mesh generators.
Created on Nov 15, 2017

@author: Richard Christie
'''

from __future__ import division
import collections
import copy
from enum import Enum
import math
import scaffoldmaker.utils.vector as vector

gaussXi3 = ( (-math.sqrt(0.6)+1.0)/2.0, 0.5, (+math.sqrt(0.6)+1.0)/2.0 )
gaussWt3 = ( 5.0/18.0, 4.0/9.0, 5.0/18.0 )

def interpolateCubicHermite(v1, d1, v2, d2, xi):
    """
    Get values of cubic Hermite interpolated from v1, d1 to v2, d2.
    :param v1, v2: Values at xi = 0.0 and xi = 1.0, respectively.
    :param d1, d2: Derivatives w.r.t. xi at xi = 0.0 and xi = 1.0, respectively.
    :param xi: Position in curve, nominally in [0.0, 1.0].
    :return: List of interpolated values at xi.
    """
    xi2 = xi*xi
    xi3 = xi2*xi
    f1 = 1.0 - 3.0*xi2 + 2.0*xi3
    f2 = xi - 2.0*xi2 + xi3
    f3 = 3.0*xi2 - 2.0*xi3
    f4 = -xi2 + xi3
    return [ (f1*v1[i] + f2*d1[i] + f3*v2[i] + f4*d2[i]) for i in range(len(v1)) ]

def interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi):
    """
    Get derivatives of cubic Hermite interpolated from v1, d1 to v2, d2.
    :param v1, v2: Values at xi = 0.0 and xi = 1.0, respectively.
    :param d1, d2: Derivatives w.r.t. xi at xi = 0.0 and xi = 1.0, respectively.
    :param xi: Position in curve, nominally in [0.0, 1.0].
    :return: List of interpolated derivatives at xi.
    """
    xi2 = xi*xi
    f1 = -6.0*xi + 6.0*xi2
    f2 = 1.0 - 4.0*xi + 3.0*xi2
    f3 = 6.0*xi - 6.0*xi2
    f4 = -2.0*xi + 3.0*xi2
    return [ (f1*v1[i] + f2*d1[i] + f3*v2[i] + f4*d2[i]) for i in range(len(v1)) ]

def interpolateCubicHermiteSecondDerivative(v1, d1, v2, d2, xi):
    """
    Get second derivatives of cubic Hermite interpolated from v1, d1 to v2, d2.
    :param v1, v2: Values at xi = 0.0 and xi = 1.0, respectively.
    :param d1, d2: Derivatives w.r.t. xi at xi = 0.0 and xi = 1.0, respectively.
    :param xi: Position in curve, nominally in [0.0, 1.0].
    :return: List of interpolated second derivatives at xi.
    """
    f1 = -6.0 + 12.0*xi
    f2 = -4.0 +  6.0*xi
    f3 =  6.0 - 12.0*xi
    f4 = -2.0 +  6.0*xi
    return [ (f1*v1[i] + f2*d1[i] + f3*v2[i] + f4*d2[i]) for i in range(len(v1)) ]

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
    :param v1, v2: Values at xi = 0.0 and xi = 1.0, respectively.
    :param d1, d2: Derivatives w.r.t. xi at xi = 0.0 and xi = 1.0, respectively.
    :param radialVector: Radial direction, assumed unit normal to curve tangent at point.
    :param xi: Position in curve, nominally in [0.0, 1.0].
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

def interpolateHermiteLagrange(v1, d1, v2, xi):
    """
    Get value at xi for quadratic Hermite-Lagrange interpolation from v1, d1 to v2.
    :return: List of values at xi.
    """
    f1 = 1 - xi*xi
    f2 = xi - xi*xi
    f3 = xi*xi
    return [ (v1[c]*f1 + d1[c]*f2 + v2[c]*f3) for c in range(len(v1)) ]

def interpolateHermiteLagrangeDerivative(v1, d1, v2, xi):
    """
    Get derivative at xi for quadratic Hermite-Lagrange interpolation from v1, d1 to v2.
    :return: List of derivatives w.r.t. xi at xi.
    """
    df1 = -2.0*xi
    df2 = 1 - 2.0*xi
    df3 = 2.0*xi
    return [ (v1[c]*df1 + d1[c]*df2 + v2[c]*df3) for c in range(len(v1)) ]

def interpolateLagrangeHermite(v1, v2, d2, xi):
    """
    Get value at xi for quadratic Lagrange-Hermite interpolation from v1 to v2, d2.
    :return: List of values at xi.
    """
    f1 = 1 - 2.0*xi + xi*xi
    f2 = 2.0*xi - xi*xi
    f3 = -xi + xi*xi
    return [ (v1[c]*f1 + d1[c]*f2 + v2[c]*f3) for c in range(len(v1)) ]

def interpolateLagrangeHermiteDerivative(v1, v2, d2, xi):
    """
    Get derivative at xi for quadratic Lagrange-Hermite interpolation to from v1 to v2, d2.
    :return: List of derivatives w.r.t. xi at xi.
    """
    df1 = -2.0 + 2.0*xi
    df2 = 2.0 - 2.0*xi
    df3 = -1.0 + 2.0*xi
    return [ (v1[c]*df1 + v2[c]*df2 + d2[c]*df3) for c in range(len(v1)) ]

def sampleCubicHermiteCurves(nx, nd1, elementsCountOut,
    addLengthStart = 0.0, addLengthEnd = 0.0,
    lengthFractionStart = 1.0, lengthFractionEnd = 1.0,
    elementLengthStartEndRatio = 1.0, arcLengthDerivatives = False):
    """
    Get systematically spaced points and derivatives over cubic Hermite interpolated
    curves with nodes nx and derivatives nd1. The first element uses the first two nodes.
    :param nx: Coordinates of nodes along curves.
    :param nd1: Derivatives of nodes along curves.
    :param addLengthStart, addLengthEnd: Extra length to add to start and end elements.
    :param lengthFractionStart, lengthFractionEnd: Fraction of mid element length for
        start and end elements. Can use in addition to AddLengths: If LengthFraction
        is 0.5 and AddLength is derivative/2.0 can blend into known derivative at start
        or end.
    :param elementLengthStartEndRatio: Start/end element length ratio, with lengths
        smoothly varying in between. Requires at least 2 elements. Applied in proportion
        to lengthFractionStart, lengthFractionEnd.
    :param arcLengthDerivatives: If True each cubic section is rescaled to arc length.
    If False (default), derivatives and distances are used as supplied.
    :return: px[], pd1[], pe[], pxi[], psf[], where pe[] and pxi[] are lists of element indices and
    and xi locations in the 'in' elements to pass to partner interpolateSample functions. psf[] is
    a list of scale factors for converting derivatives from old to new xi coordinates: dxi(old)/dxi(new).
    """
    elementsCountIn = len(nx) - 1
    assert (elementsCountIn > 0) and (len(nd1) == (elementsCountIn + 1)) and \
        (elementsCountOut > 0), 'sampleCubicHermiteCurves.  Invalid arguments'
    lengths = [ 0.0 ]
    nd1a = []
    nd1b = []
    length = 0.0
    for e in range(elementsCountIn):
        if arcLengthDerivatives:
            arcLength = computeCubicHermiteArcLength(nx[e], nd1[e], nx[e + 1], nd1[e + 1], rescaleDerivatives=True)
            nd1a.append(vector.setMagnitude(nd1[e], arcLength))
            nd1b.append(vector.setMagnitude(nd1[e + 1], arcLength))
        else:
            arcLength = getCubicHermiteArcLength(nx[e], nd1[e], nx[e + 1], nd1[e + 1])
        length += arcLength
        lengths.append(length)
    proportionEnd = 2.0/(elementLengthStartEndRatio + 1)
    proportionStart = elementLengthStartEndRatio*proportionEnd
    elementLengthMid = (length - addLengthStart - addLengthEnd) / \
        (elementsCountOut - 2.0 + proportionStart*lengthFractionStart + proportionEnd*lengthFractionEnd)
    elementLengthProportionStart = proportionStart*lengthFractionStart*elementLengthMid
    elementLengthProportionEnd = proportionEnd*lengthFractionEnd*elementLengthMid
    # get smoothly varying element lengths, not accounting for start and end
    elementLengths = []
    for eOut in range(elementsCountOut):
        xi = eOut/(elementsCountOut - 1)
        elementLengths.append(((1.0 - xi)*proportionStart + xi*proportionEnd)*elementLengthMid)
    # get middle derivative magnitudes
    nodeDerivativeMagnitudes = [ None ]*(elementsCountOut + 1)  # start and end determined below
    for n in range(1, elementsCountOut):
        nodeDerivativeMagnitudes[n] = 0.5*(elementLengths[n - 1] + elementLengths[n])
    # fix end lengths:
    elementLengths[ 0] = addLengthStart + elementLengthProportionStart
    elementLengths[-1] = addLengthEnd + elementLengthProportionEnd
    # set end derivatives:
    if elementsCountOut == 1:
        nodeDerivativeMagnitudes[0] = nodeDerivativeMagnitudes[1] = elementLengths[0]
    else:
        nodeDerivativeMagnitudes[0] = elementLengths[0]*2.0 - nodeDerivativeMagnitudes[1]
        nodeDerivativeMagnitudes[-1] = elementLengths[-1]*2.0 - nodeDerivativeMagnitudes[-2]

    px = []
    pd1 = []
    pe = []
    pxi = []
    psf = []
    distance = 0.0
    e = 0
    for eOut in range(elementsCountOut):
        while e < elementsCountIn:
            if distance < lengths[e + 1]:
                partDistance = distance - lengths[e]
                if arcLengthDerivatives:
                    xi = partDistance/(lengths[e + 1] - lengths[e])
                    x = interpolateCubicHermite(nx[e], nd1a[e], nx[e + 1], nd1b[e], xi)
                    d1 = interpolateCubicHermiteDerivative(nx[e], nd1a[e], nx[e + 1], nd1b[e], xi)
                else:
                    x, d1, _eIn, xi = getCubicHermiteCurvesPointAtArcDistance(nx[e:e + 2], nd1[e:e + 2], partDistance)
                sf = nodeDerivativeMagnitudes[eOut]/vector.magnitude(d1)
                px.append(x)
                pd1.append([ sf*d for d in d1 ])
                pe.append(e)
                pxi.append(xi)
                psf.append(sf)
                break
            e += 1
        distance += elementLengths[eOut]
    e = elementsCountIn
    eOut = elementsCountOut
    xi = 1.0
    d1 = nd1[e]
    sf = nodeDerivativeMagnitudes[eOut]/vector.magnitude(d1)
    px.append(nx[e])
    pd1.append([ sf*d for d in d1 ])
    pe.append(e - 1)
    pxi.append(xi)
    psf.append(sf)
    return px, pd1, pe, pxi, psf

def interpolateSampleCubicHermite(v, d, pe, pxi, psf):
    '''
    Partner function to sampleCubicHermiteCurves for interpolating additional variables with
    cubic Hermite basis, at the element indexes, xi coordinates and xi scaling returned from that function.
    Note: this does not work for sampleCubicHermiteCurves with arcLengthDerivatives = False.
    :param v, d: List of values and derivatives to interpolate, either scalar or sequence-of-scalar.
    len(v) == len(d) == number of elements in + 1.
    :param pe, pxi: List if integer element indexes and real xi coordinates giving sample positions into v to
    interpolate linearly. len(pe) == len(pxi) == number of values out.
    Indexes in pe start at 0, and are not checked; sampleCubicHermiteCurves() guarantees these are valid for
    the number of elements passed to it.
    :param psf: List of scale factors dxi(old)/dxi(new). Length same as pe, pxi. Used to convert derivatives
    from old to new xi spacing.
    :return: List of interpolated values, list of interpolated derivatives; scalar or vector as for v, d.
    '''
    assert (len(v) > 1) and (len(d) == len(v)), 'interpolateSampleLinear. Invalid values v, d'
    valuesCountOut = len(pe)
    assert (valuesCountOut > 0) and (len(pxi) == valuesCountOut), 'interpolateSampleLinear. Invalid element, xi'
    vOut = []
    dOut = []
    if isinstance(v[0], collections.Sequence):
        for n in range(valuesCountOut):
            e = pe[n]
            v1 = v[e]
            d1 = d[e]
            v2 = v[e + 1]
            d2 = d[e + 1]
            vOut.append(interpolateCubicHermite(v1, d1, v2, d2, pxi[n]))
            dOut.append([ psf[n]*d for d in interpolateCubicHermiteDerivative(v1, d1, v2, d2, pxi[n]) ])
    else:
        for n in range(valuesCountOut):
            e = pe[n]
            v1 = [ v[e] ]
            d1 = [ d[e] ]
            v2 = [ v[e + 1] ]
            d2 = [ d[e + 1] ]
            vOut.append(interpolateCubicHermite(v1, d1, v2, d2, pxi[n])[0])
            dOut.append(psf[n]*interpolateCubicHermiteDerivative(v1, d1, v2, d2, pxi[n])[0])
    return vOut, dOut

def interpolateSampleLinear(v, pe, pxi):
    '''
    Partner function to sampleCubicHermiteCurves for linearly interpolating additional variables based on the 
    element indexes and element xi coordinates returned from that function.
    :param v: List of scalar values or sequence-of-values to interpolate. len(v) == number of elements in + 1.
    :param pe, pxi: List if integer element indexes and real xi coordinates giving sample positions into v to
    interpolate linearly. len(pe) == len(pxi) == number of values out.
    Indexes in pe start at 0, and are not checked; sampleCubicHermiteCurves() guarantees these are valid for
    the number of elements passed to it.
    :return: List of interpolated values, scalar or vector as for v.
    '''
    assert len(v) > 1, 'interpolateSampleLinear. Invalid values v: not enough data'
    valuesCountOut = len(pe)
    assert (valuesCountOut > 0) and (len(pxi) == valuesCountOut), 'interpolateSampleLinear. Invalid element, xi'
    vOut = []
    if isinstance(v[0], collections.Sequence):
        vLen = len(v[0])
        for n in range(valuesCountOut):
            wp = pxi[n]
            wm = 1.0 - wp
            vp = v[pe[n] + 1]
            vm = v[pe[n]]
            vOut.append([ (wm*vm[c] + wp*vp[c]) for c in range(vLen) ])
    else:
        for n in range(valuesCountOut):
            wp = pxi[n]
            wm = 1.0 - wp
            vOut.append(wm*v[pe[n]] + wp*v[pe[n] + 1])
    return vOut

def getCubicHermiteCurvesPointAtArcDistance(nx, nd, arcDistance):
    """
    Get the coordinates, derivatives at distance along cubic Hermite curves.
    Supplied derivatives are used i.e. not rescaled to arc length.
    Note this is approximate.
    :param nx: Coordinates of nodes along curves.
    :param nd: Derivatives of nodes along curves.
    :param distance: Distance along curves.
    :return: coordinates, derivatives, element index, xi; clamped to first or last nx if distance is beyond curves
    """
    elementsCount = len(nx) - 1
    assert elementsCount > 0, 'getCubicHermiteCurvesPointAtArcDistance.  Invalid number of points'
    if arcDistance < 0.0:
        return nx[0], nd[0], 0, 0.0
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
                    return interpolateCubicHermite(v1, d1, v2, d2, xi), interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi), e, xi
            print('getCubicHermiteCurvesPointAtArcDistance Max iters reached:',iter,': e', e, ', xi',xi,', closeness', math.fabs(dist - partDistance))
            return v2, d2, e, xi
        length += arcLength
    return nx[-1], nd[-1], elementsCount - 1, 1.0

class DerivativeScalingMode(Enum):
    ARITHMETIC_MEAN = 1  # derivative is half of sum of arclengths on either side
    HARMONIC_MEAN   = 2  # derivative is reciprocal of arithmetic mean of reciprocals of arclengths = arc lengths weighted by proportion from other side

def smoothCubicHermiteDerivativesLine(nx, nd1,
        fixAllDirections = False,
        fixStartDerivative = False, fixEndDerivative = False,
        fixStartDirection = False, fixEndDirection = False,
        magnitudeScalingMode = DerivativeScalingMode.ARITHMETIC_MEAN, instrument=False):
    """
    Modifies derivatives nd1 to be smoothly varying and near arc length.
    Values are treated as being in a line.
    Assumes initial derivatives are zero or reasonable.
    Where directions are smoothed the weighted/harmonic mean is used.
    :param nx: List of coordinates of nodes along curves.
    :param nd1: List of derivatives of nodes along curves.
    :param fixAllDirections: Set to True to only smooth magnitudes, otherwise both direction and magnitude are adjusted.
    :param fixStartDerivative, fixEndDerivative: Set to True to fix derivative direction and magnitude at respective end.
    :param fixStartDirection, fixEndDirection: Set to True to fix direction at respective end.
    Redundant if fixAllDirections or respective fixStart/EndDerivative is True.
    :param magnitudeScalingMode: A value from enum DerivativeScalingMode specifying
    expression used to get derivative magnitude from adjacent arc lengths.
    :return: Modified nd1
    """
    nodesCount = len(nx)
    elementsCount = nodesCount - 1
    assert elementsCount > 0, 'smoothCubicHermiteDerivativesLine.  Too few nodes/elements'
    assert len(nd1) == nodesCount, 'smoothCubicHermiteDerivativesLine.  Mismatched number of derivatives'
    arithmeticMeanMagnitude = magnitudeScalingMode is DerivativeScalingMode.ARITHMETIC_MEAN
    assert arithmeticMeanMagnitude or (magnitudeScalingMode is DerivativeScalingMode.HARMONIC_MEAN), \
        'smoothCubicHermiteDerivativesLine. Invalid magnitude scaling mode'
    md1 = copy.copy(nd1)
    componentsCount = len(nx[0])
    componentRange = range(componentsCount)
    tol = 1.0E-6
    # special case where equal derivatives at each end are sought
    equalDerivatives = (elementsCount == 1) and not (fixStartDerivative or fixEndDerivative)
    if instrument:
        print('iter 0', md1)
    for iter in range(100):
        lastmd1 = copy.copy(md1)
        arcLengths = [ getCubicHermiteArcLength(nx[e], md1[e], nx[e + 1], md1[e + 1]) for e in range(elementsCount) ]
        # start
        if not fixStartDerivative:
            if fixAllDirections or fixStartDirection:
                md1[0] = vector.setMagnitude(lastmd1[0], 2.0*arcLengths[0] - vector.magnitude(lastmd1[1]))
            else:
                md1[0] = interpolateLagrangeHermiteDerivative(nx[0], nx[1], lastmd1[1], 0.0)
        # middle
        for n in range(1, nodesCount - 1):
            nm = n - 1
            if not fixAllDirections:
                # get mean of directions from point n to points (n - 1) and (n + 1)
                np = n + 1
                dirm = [ (nx[n ][c] - nx[nm][c]) for c in componentRange ]
                dirp = [ (nx[np][c] - nx[n ][c]) for c in componentRange ]
                # mean weighted by fraction towards that end, equivalent to harmonic mean
                arcLengthmp = arcLengths[nm] + arcLengths[n]
                wm = arcLengths[n ]/arcLengthmp
                wp = arcLengths[nm]/arcLengthmp
                md1[n] = [ (wm*dirm[c] + wp*dirp[c]) for c in componentRange ]
            if arithmeticMeanMagnitude:
                mag = 0.5*(arcLengths[nm] + arcLengths[n])
            else: # harmonicMeanMagnitude
                mag = 2.0/(1.0/arcLengths[nm] + 1.0/arcLengths[n])
            md1[n] = vector.setMagnitude(md1[n], mag)
        # end
        if not fixEndDerivative:
            if fixAllDirections or fixEndDirection:
                md1[-1] = vector.setMagnitude(lastmd1[-1], 2.0*arcLengths[-1] - vector.magnitude(lastmd1[-2]))
            else:
                md1[-1] = interpolateHermiteLagrangeDerivative(nx[-2], lastmd1[-2], nx[-1], 1.0)
        if equalDerivatives:
            mag = getCubicHermiteArcLength(nx[0], md1[0], nx[1], md1[1])
            md1[0] = vector.setMagnitude(md1[0], mag)
            md1[1] = vector.setMagnitude(md1[1], mag)
        if instrument:
            print('iter', iter + 1, md1)
        dtol = tol*sum(arcLengths)/len(arcLengths)
        for n in range(nodesCount):
            for c in componentRange:
                if math.fabs(md1[n][c] - lastmd1[n][c]) > dtol:
                    break
            else:
                continue
            break
        else:
            if instrument:
                print('smoothCubicHermiteDerivativesLine converged after iter:',iter)
            return md1

    print('smoothCubicHermiteDerivativesLine max iters reached:',iter)
    return md1

def smoothCubicHermiteDerivativesLoop(nx, nd1,
        fixAllDirections = False,
        magnitudeScalingMode = DerivativeScalingMode.ARITHMETIC_MEAN, instrument=False):
    """
    Modifies derivatives nd1 to be smoothly varying and near arc length.
    Values are treated as being in a loop, so the first point follows the last.
    Assumes initial derivatives are zero or reasonable.
    Where directions are smoothed the weighted/harmonic mean is used.
    :param nx: List of coordinates of nodes along curves.
    :param nd1: List of derivatives of nodes along curves.
    :param fixAllDirections: Set to True to only smooth magnitudes, otherwise both direction and magnitude are adjusted.
    :param magnitudeScalingMode: A value from enum DerivativeScalingMode specifying
    expression used to get derivative magnitude from adjacent arc lengths.
    :return: Modified nd1
    """
    nodesCount = elementsCount = len(nx)
    assert elementsCount > 1, 'smoothCubicHermiteDerivativesLoop.  Too few nodes/elements'
    assert len(nd1) == elementsCount, 'smoothCubicHermiteDerivativesLoop.  Mismatched number of derivatives'
    arithmeticMeanMagnitude = magnitudeScalingMode is DerivativeScalingMode.ARITHMETIC_MEAN
    assert arithmeticMeanMagnitude or (magnitudeScalingMode is DerivativeScalingMode.HARMONIC_MEAN), \
        'smoothCubicHermiteDerivativesLine. Invalid magnitude scaling mode'
    md1 = copy.copy(nd1)
    componentsCount = len(nx[0])
    componentRange = range(componentsCount)
    tol = 1.0E-6
    if instrument:
        print('iter 0', md1)
    for iter in range(100):
        lastmd1 = copy.copy(md1)
        arcLengths = [ getCubicHermiteArcLength(nx[e], md1[e], nx[(e + 1)%elementsCount], md1[(e + 1)%elementsCount]) for e in range(elementsCount) ]
        for n in range(nodesCount):
            nm = n - 1
            if not fixAllDirections:
                # get mean of directions from point n to points (n - 1) and (n + 1)
                np = (n + 1)%nodesCount
                dirm = [ (nx[n ][c] - nx[nm][c]) for c in componentRange ]
                dirp = [ (nx[np][c] - nx[n ][c]) for c in componentRange ]
                # mean weighted by fraction towards that end, equivalent to harmonic mean
                arcLengthmp = arcLengths[nm] + arcLengths[n]
                wm = arcLengths[n ]/arcLengthmp
                wp = arcLengths[nm]/arcLengthmp
                md1[n] = [ (wm*dirm[c] + wp*dirp[c]) for c in componentRange ]
            if arithmeticMeanMagnitude:
                mag = 0.5*(arcLengths[nm] + arcLengths[n])
            else: # harmonicMeanMagnitude
                mag = 2.0/(1.0/arcLengths[nm] + 1.0/arcLengths[n])
            md1[n] = vector.setMagnitude(md1[n], mag)
        if instrument:
            print('iter', iter + 1, md1)
        dtol = tol*sum(arcLengths)/len(arcLengths)
        for n in range(nodesCount):
            for c in componentRange:
                if math.fabs(md1[n][c] - lastmd1[n][c]) > dtol:
                    break
            else:
                continue
            break
        else:
            if instrument:
                print('smoothCubicHermiteDerivativesLoop converged after iter:',iter)
            return md1

    print('smoothCubicHermiteDerivativesLoop max iters reached:',iter)
    return md1

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
