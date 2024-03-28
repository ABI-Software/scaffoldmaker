"""
Utilities for building bifurcating network meshes.
"""

from cmlibs.maths.vectorops import add, cross, dot, magnitude, mult, normalize, sub
from cmlibs.utils.zinc.general import ChangeManager
# from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils.interpolation import computeCubicHermiteArcLength, computeCubicHermiteDerivativeScaling, \
    DerivativeScalingMode, getCubicHermiteCurvesLength, getNearestLocationBetweenCurves, getNearestLocationOnCurve, \
    interpolateCubicHermite, interpolateCubicHermiteDerivative, interpolateCubicHermiteSecondDerivative, \
    interpolateLagrangeHermiteDerivative, smoothCubicHermiteDerivativesLine, smoothCubicHermiteDerivativesLoop, \
    smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.networkmesh import NetworkMesh, getPathRawTubeCoordinates, resampleTubeCoordinates
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_delta_xi
from scaffoldmaker.utils.zinc_utils import generateCurveMesh, get_nodeset_path_ordered_field_parameters
import copy
import math


def get_curve_circle_points(x1, xd1, x2, xd2, r1, rd1, r2, rd2, xi, dmag, side, elementsCountAround):
    """
    :param dmag: Magnitude of derivative on curve.
    :param side: Vector in side direction of first node around.
    Need not be unit or exactly normal to curve at xi.
    :return: x[], d1[] around, d2[] along
    """
    cx = interpolateCubicHermite(x1, xd1, x2, xd2, xi)
    cxd = interpolateCubicHermiteDerivative(x1, xd1, x2, xd2, xi)
    mag_cxd = magnitude(cxd)
    cxd2 = interpolateCubicHermiteSecondDerivative(x1, xd1, x2, xd2, xi)
    mag_cxd2 = magnitude(cxd2)
    r = interpolateCubicHermite([ r1 ], [ rd1 ], [ r2 ], [ rd2 ], xi)[0]
    rd = interpolateCubicHermiteDerivative([ r1 ], [ rd1 ], [ r2 ], [ rd2 ], xi)[0]
    axis1 = normalize(cxd)
    axis3 = normalize(cross(axis1, side))
    axis2 = cross(axis3, axis1)
    x, d1 = createCirclePoints(cx, mult(axis2, r), mult(axis3, r), elementsCountAround)
    curvatureVector = mult(cxd2, 1.0/(mag_cxd*mag_cxd))
    d2 = []
    radialGrowth = rd/(mag_cxd*r)
    for e in range(elementsCountAround):
        radialVector = sub(x[e], cx)
        dmagFinal = dmag*(1.0 - dot(radialVector, curvatureVector))
        # add curvature and radius change components:
        d2.append(add(mult(cxd, dmagFinal/mag_cxd), mult(radialVector, dmagFinal*radialGrowth)))
    return x, d1, d2

def get_bifurcation_triple_point(p1x, p1d, p2x, p2d, p3x, p3d):
    """
    Get coordinates and derivatives of triple point between p1, p2 and p3 with derivatives.
    :param p1x..p3d: Point coordinates and derivatives, numbered anticlockwise around triple point.
    All derivatives point away from triple point.
    Returned d1 points from triple point to p2, d2 points from triple point to p3.
    :return: x, d1, d2
    """
    scaling12 = computeCubicHermiteDerivativeScaling(p1x, [-d for d in p1d], p2x, p2d)
    scaling23 = computeCubicHermiteDerivativeScaling(p2x, [-d for d in p2d], p3x, p3d)
    scaling31 = computeCubicHermiteDerivativeScaling(p3x, [-d for d in p3d], p1x, p1d)
    trx1 = interpolateCubicHermite(p1x, mult(p1d, -scaling12), p2x, mult(p2d, scaling12), 0.5)
    trx2 = interpolateCubicHermite(p2x, mult(p2d, -scaling23), p3x, mult(p3d, scaling23), 0.5)
    trx3 = interpolateCubicHermite(p3x, mult(p3d, -scaling31), p1x, mult(p1d, scaling31), 0.5)
    trx = [(trx1[c] + trx2[c] + trx3[c]) / 3.0 for c in range(3)]
    td1 = interpolateLagrangeHermiteDerivative(trx, p1x, p1d, 0.0)
    td2 = interpolateLagrangeHermiteDerivative(trx, p2x, p2d, 0.0)
    td3 = interpolateLagrangeHermiteDerivative(trx, p3x, p3d, 0.0)
    n12 = cross(td1, td2)
    n23 = cross(td2, td3)
    n31 = cross(td3, td1)
    norm = normalize([(n12[c] + n23[c] + n31[c]) for c in range(3)])
    sd1 = smoothCubicHermiteDerivativesLine([trx, p1x], [normalize(cross(norm, cross(td1, norm))), p1d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    sd2 = smoothCubicHermiteDerivativesLine([trx, p2x], [normalize(cross(norm, cross(td2, norm))), p2d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    sd3 = smoothCubicHermiteDerivativesLine([trx, p3x], [normalize(cross(norm, cross(td3, norm))), p3d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    trd1 = mult(sub(sd2, add(sd3, sd1)), 0.5)
    trd2 = mult(sub(sd3, add(sd1, sd2)), 0.5)
    return trx, trd1, trd2


def get_tube_bifurcation_connection_elements_counts(tCounts):
    """
    Get number of elements directly connecting tubes 1, 2 and 3 from the supplied number around.
    :param tCounts: Number of elements around tubes in order.
    :return: List of elements connect tube with its next neighbour, looping back to first.
    """
    assert len(tCounts) == 3
    return [(tCounts[i] + tCounts[i - 2] - tCounts[i - 1]) // 2 for i in range(3)]


def make_tube_bifurcation_points(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):
    """
    Gets first ring of coordinates and derivatives between parent pa and
    children c1, c2, and over the crotch between c1 and c2.
    :return rox, rod1, rod2, cox, cod1, cod2
    """
    paCount = len(pax)
    c1Count = len(c1x)
    c2Count = len(c2x)
    pac1Count, c1c2Count, pac2Count = get_tube_bifurcation_connection_elements_counts([paCount, c1Count, c2Count])
    # convert to number of nodes, includes both 6-way points
    pac1NodeCount = pac1Count + 1
    pac2NodeCount = pac2Count + 1
    c1c2NodeCount = c1c2Count + 1
    paStartIndex = 0
    c1StartIndex = 0
    c2StartIndex = 0
    pac1x  = [ None ]*pac1NodeCount
    pac1d1 = [ None ]*pac1NodeCount
    pac1d2 = [ None ]*pac1NodeCount
    for n in range(pac1NodeCount):
        pan = (paStartIndex + n) % paCount
        c1n = (c1StartIndex + n) % c1Count
        x1, d1, x2, d2 = pax[pan], mult(pad2[pan], 2.0), c1x[c1n], mult(c1d2[c1n], 2.0)
        pac1x [n] = interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac1d1[n] = [ 0.0, 0.0, 0.0 ]
        pac1d2[n] = mult(interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    paStartIndex2 = paStartIndex + pac1Count
    c1StartIndex2 = c1StartIndex + pac1Count
    c2StartIndex2 = c2StartIndex + c1c2Count
    pac2x  = [ None ]*pac2NodeCount
    pac2d1 = [ None ]*pac2NodeCount
    pac2d2 = [ None ]*pac2NodeCount
    for n in range(pac2NodeCount):
        pan = (paStartIndex2 + n) % paCount
        c2n = (c2StartIndex2 + n) % c2Count
        x1, d1, x2, d2 = pax[pan], mult(pad2[pan], 2.0), c2x[c2n], mult(c2d2[c2n], 2.0)
        pac2x [n] = interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac2d1[n] = [ 0.0, 0.0, 0.0 ]
        pac2d2[n] = mult(interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    c1c2x  = [ None ]*c1c2NodeCount
    c1c2d1 = [ None ]*c1c2NodeCount
    c1c2d2 = [ None ]*c1c2NodeCount
    for n in range(c1c2NodeCount):
        c1n = (c1StartIndex2 + n) % c1Count
        c2n = (c2StartIndex2 - n) % c2Count  # note: reversed
        x1, d1, x2, d2 = c2x[c2n], mult(c2d2[c2n], -2.0), c1x[c1n], mult(c1d2[c1n], 2.0)
        c1c2x [n] = interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        c1c2d1[n] = [ 0.0, 0.0, 0.0 ]
        c1c2d2[n] = mult(interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    # get hex triple points
    hex1, hex1d1, hex1d2 = get_bifurcation_triple_point(
        pax[paStartIndex], mult(pad2[paStartIndex], -1.0),
        c1x[c1StartIndex], c1d2[c1StartIndex],
        c2x[c1StartIndex], c2d2[c2StartIndex])
    hex2, hex2d1, hex2d2 = get_bifurcation_triple_point(
        pax[paStartIndex2], mult(pad2[paStartIndex2], -1.0),
        c2x[c2StartIndex2], c2d2[c2StartIndex2],
        c1x[c1StartIndex2], c1d2[c1StartIndex2])
    # smooth around loops through hex points to get d1
    loop1x  = [ hex2 ] + pac2x[1:-1] + [ hex1 ]
    loop1d1 = [ [ -d for d in hex2d2 ] ] + pac2d1[1:-1] + [ hex1d1 ]
    loop2x  = [ hex1 ] + pac1x[1:-1] + [ hex2 ]
    loop2d1 = [ [ -d for d in hex1d2 ] ] + pac1d1[1:-1] + [ hex2d1 ]
    loop1d1 = smoothCubicHermiteDerivativesLine(loop1x, loop1d1, fixStartDirection=True, fixEndDirection=True, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
    loop2d1 = smoothCubicHermiteDerivativesLine(loop2x, loop2d1, fixStartDirection=True, fixEndDirection=True, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
    # smooth over "crotch" between c1 and c2
    crotchx = [ hex2 ] + c1c2x[1:-1] + [ hex1 ]
    crotchd1 = [ add(hex2d1, hex2d2) ] + c1c2d1[1:-1] + [ [ (-hex1d1[c] - hex1d2[c]) for c in range(3) ] ]
    crotchd1 = smoothCubicHermiteDerivativesLine(crotchx, crotchd1, fixStartDerivative=True, fixEndDerivative=True, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
    rox  = [ hex1 ] + pac1x[1:-1] + [ hex2 ] + pac2x[1:-1]
    rod1 = [ loop1d1[-1] ] + loop2d1[1:] + loop1d1[1:-1]
    rod2 = [ [ -d for d in loop2d1[ 0] ] ] + pac1d2[1:-1] + [ [ -d for d in loop1d1[0] ] ] + pac2d2[1:-1]
    cox  = crotchx [1:-1]
    cod1 = crotchd1[1:-1]
    cod2 = c1c2d2[1:-1]
    return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex


def make_tube_bifurcation_elements_2d(region, coordinates, elementIdentifier,
        paNodeId, paStartIndex, c1NodeId, c1StartIndex, c2NodeId, c2StartIndex, roNodeId, coNodeId,
        useCrossDerivatives=False):
    """
    Creates elements from parent, ring/row, crotch/column, child1 and child2 nodes.
    Assumes client has active ChangeManager(fieldmodule).
    :param region: Zinc region to create model in.
    :param coordinates: Finite element coordinate field to define.
    :param elementIdentifier: First 2D element identifier to use.
    :param paNodeId, paStartIndex, c1NodeId, c1StartIndex, c2NodeId, c2StartIndex:
    Lists of parent, child1, child2 nodes and their starting index to be at hex2
    :param roNodeId, coNodeId: Lists of ring/row and crotch/column nodes,
    starting at hex2 and between hex1 and hex2, respectively.
    :return next element identifier.
    """
    paCount = len(paNodeId)
    c1Count = len(c1NodeId)
    c2Count = len(c2NodeId)
    pac1Count, c1c2Count, pac2Count = get_tube_bifurcation_connection_elements_counts([paCount, c1Count, c2Count])

    fieldmodule = region.getFieldmodule()
    mesh = fieldmodule.findMeshByDimension(2)

    elementtemplateStd = mesh.createElementtemplate()
    elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
    bicubicHermiteBasis = fieldmodule.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    eftStd = mesh.createElementfieldtemplate(bicubicHermiteBasis)
    if not useCrossDerivatives:
        for n in range(4):
            eftStd.setFunctionNumberOfTerms(n*4 + 4, 0)
    elementtemplateStd.defineField(coordinates, -1, eftStd)
    elementtemplateMod = mesh.createElementtemplate()
    elementtemplateMod.setElementShapeType(Element.SHAPE_TYPE_SQUARE)

    for e1 in range(paCount):
        eft = eftStd
        elementtemplate = elementtemplateStd
        np = e1 + paStartIndex
        nids = [ paNodeId[np % paCount], paNodeId[(np + 1) % paCount], roNodeId[e1], roNodeId[(e1 + 1) % paCount] ]
        scalefactors = None
        meshGroups = [ ]

        if e1 in (0, pac1Count - 1, pac1Count, paCount - 1):
            eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
            if not useCrossDerivatives:
                for n in range(4):
                    eft.setFunctionNumberOfTerms(n*4 + 4, 0)
            if e1 in (0, pac1Count):
                scalefactors = [ -1.0 ]
                setEftScaleFactorIds(eft, [1], [])
                remapEftNodeValueLabel(eft, [ 3 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [ 1 ] ) ])
                remapEftNodeValueLabel(eft, [ 3 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
            elif e1 in (pac1Count - 1, paCount - 1):
                remapEftNodeValueLabel(eft, [ 4 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS2, [] ) ])
            elementtemplateMod.defineField(coordinates, -1, eft)
            elementtemplate = elementtemplateMod

        element = mesh.createElement(elementIdentifier, elementtemplate)
        result2 = element.setNodesByIdentifier(eft, nids)
        if scalefactors:
            result3 = element.setScaleFactors(eft, scalefactors)
        else:
            result3 = '-'
        #print('create element tube bifurcation pa', element.isValid(), elementIdentifier, result2, result3, nids)
        elementIdentifier += 1
        for meshGroup in meshGroups:
            meshGroup.addElement(element)

    for e1 in range(c1Count):
        eft = eftStd
        elementtemplate = elementtemplateStd
        nr = e1
        nc = e1 + c1StartIndex
        nids = [ roNodeId[nr % paCount], roNodeId[(nr + 1) % paCount], c1NodeId[nc % c1Count], c1NodeId[(nc + 1) % c1Count] ]
        if e1 >= pac1Count:
            nids[1] = roNodeId[0] if (e1 == (c1Count - 1)) else coNodeId[e1 - pac1Count]
            if e1 > pac1Count:
                #nids[0] = coNodeId[(e1 - pac1Count - 1)]
                nids[0] = coNodeId[(e1 - pac1Count - 1) % c1Count]
        scalefactors = None
        meshGroups = [ ]

        if e1 in (0, pac1Count, c1Count - 1):
            eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
            if not useCrossDerivatives:
                for n in range(4):
                    eft.setFunctionNumberOfTerms(n*4 + 4, 0)
            if e1 == 0:
                scalefactors = [ -1.0 ]
                setEftScaleFactorIds(eft, [1], [])
                remapEftNodeValueLabel(eft, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [ 1 ] ) ])
                remapEftNodeValueLabel(eft, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
            elif e1 == pac1Count:
                remapEftNodeValueLabel(eft, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [ ] ), ( Node.VALUE_LABEL_D_DS2, [ ] ) ])
            elif e1 == (c1Count - 1):
                scalefactors = [ -1.0 ]
                setEftScaleFactorIds(eft, [1], [])
                remapEftNodeValueLabel(eft, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [ 1 ] ), ( Node.VALUE_LABEL_D_DS2, [ 1 ] ) ])
                remapEftNodeValueLabel(eft, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [] ) ])
            elementtemplateMod.defineField(coordinates, -1, eft)
            elementtemplate = elementtemplateMod

        element = mesh.createElement(elementIdentifier, elementtemplate)
        result2 = element.setNodesByIdentifier(eft, nids)
        if scalefactors:
            result3 = element.setScaleFactors(eft, scalefactors)
        else:
            result3 = '-'
        #print('create element tube bifurcation c1', element.isValid(), elementIdentifier, result2, result3, nids)
        elementIdentifier += 1
        for meshGroup in meshGroups:
            meshGroup.addElement(element)

    for e1 in range(c2Count):
        eft = eftStd
        elementtemplate = elementtemplateStd
        nr = 0 if (e1 == 0) else (paCount - c2Count + e1)
        nc = e1 + c2StartIndex
        nids = [ roNodeId[nr % paCount], roNodeId[(nr + 1) % paCount], c2NodeId[nc % c2Count], c2NodeId[(nc + 1) % c2Count] ]
        if 0 <= e1 < (c1c2Count - 1):
            nids[1] = coNodeId[c1c2Count - e1 - 2]
        if 0 < e1 <= (c1c2Count - 1):
            nids[0] = coNodeId[c1c2Count - e1 - 1]
        scalefactors = None
        meshGroups = [ ]

        if e1 <= c1c2Count:
            eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
            if not useCrossDerivatives:
                for n in range(4):
                    eft.setFunctionNumberOfTerms(n*4 + 4, 0)
            scalefactors = [ -1.0 ]
            setEftScaleFactorIds(eft, [1], [])
            if e1 == 0:
                remapEftNodeValueLabel(eft, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [ ] ), ( Node.VALUE_LABEL_D_DS2, [ ] ) ])
                scaleEftNodeValueLabels(eft, [ 2 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], [ 1 ])
            elif e1 < (c1c2Count - 1):
                scaleEftNodeValueLabels(eft, [ 1, 2 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], [ 1 ])
            elif e1 == (c1c2Count - 1):
                scaleEftNodeValueLabels(eft, [ 1 ], [ Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2 ], [ 1 ])
                remapEftNodeValueLabel(eft, [ 2 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [ 1 ] ), ( Node.VALUE_LABEL_D_DS2, [ 1 ] ) ])
                remapEftNodeValueLabel(eft, [ 2 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [ ] ) ])
            elif e1 == c1c2Count:
                remapEftNodeValueLabel(eft, [ 1 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS2, [ 1 ] ) ])
                remapEftNodeValueLabel(eft, [ 1 ], Node.VALUE_LABEL_D_DS2, [ ( Node.VALUE_LABEL_D_DS1, [ ] ) ])
            elementtemplateMod.defineField(coordinates, -1, eft)
            elementtemplate = elementtemplateMod

        element = mesh.createElement(elementIdentifier, elementtemplate)
        result2 = element.setNodesByIdentifier(eft, nids)
        if scalefactors:
            result3 = element.setScaleFactors(eft, scalefactors)
        else:
            result3 = '-'
        #print('create element tube bifurcation c2', element.isValid(), elementIdentifier, result2, result3, nids)
        for meshGroup in meshGroups:
            meshGroup.addElement(element)
        elementIdentifier += 1

    return elementIdentifier


def getBifurcationCrossPoint(cx, p1x, p1d, p2x, p2d, p3x, p3d):
    """
    Get derivatives of cross point between p1, p2 and p3 with derivatives.
    :param cx: Cross point coordinates.
    :param p1x..p3d: Point coordinates and derivatives, numbered anticlockwise around triple point.
    All derivatives point away from triple point.
    Returned d1 points from triple point to p2, d2 points from triple point to p3.
    :return: d1, d2
    """
    # cx1 = interpolateCubicHermite(p1x, mult(p1d, -2.0), p2x, mult(p2d, 2.0), 0.5)
    # cx2 = interpolateCubicHermite(p2x, mult(p2d, -2.0), p3x, mult(p3d, 2.0), 0.5)
    # cx3 = interpolateCubicHermite(p3x, mult(p3d, -2.0), p1x, mult(p1d, 2.0), 0.5)
    # cx = [ (cx1[c] + cx2[c] + cx3[c])/3.0 for c in range(3) ]
    td1 = interpolateLagrangeHermiteDerivative(cx, p1x, p1d, 0.0)
    td2 = interpolateLagrangeHermiteDerivative(cx, p2x, p2d, 0.0)
    td3 = interpolateLagrangeHermiteDerivative(cx, p3x, p3d, 0.0)
    n12 = cross(td1, td2)
    n23 = cross(td2, td3)
    n31 = cross(td3, td1)
    norm = normalize([(n12[c] + n23[c] + n31[c]) for c in range(3)])
    sd1 = smoothCubicHermiteDerivativesLine([cx, p1x], [normalize(cross(norm, cross(td1, norm))), p1d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    sd2 = smoothCubicHermiteDerivativesLine([cx, p2x], [normalize(cross(norm, cross(td2, norm))), p2d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    sd3 = smoothCubicHermiteDerivativesLine([cx, p3x], [normalize(cross(norm, cross(td3, norm))), p3d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    cd1 = mult(sub(sd2, add(sd3, sd1)), 0.5)
    cd2 = mult(sub(sd3, add(sd1, sd2)), 0.5)
    return cd1, cd2


def getTubeBifurcationCoordinates2D(tCoords, inCount):
    """
    Compute half loops of middle coordinates between tube and its next neighbour.
    :param tCoords: List over 3 tubes (starting with "in" tubes) of [tx, td1, td2] for last/first tube ring in/out.
    :param inCount: Number of tubes which are directed in to the junction.
    :return: mCoords (list over tubes of mid coordinates x, d1, d2 to next tube, around half loop from first to
    second cross point).
    """
    tCount = len(tCoords)
    assert tCount == 3
    taCounts = [len(v[0]) for v in tCoords]
    # get numbers of elements directly connecting t0-t1, t1-t2, t2-t1
    teCounts = [(taCounts[i] + taCounts[i - tCount + 1] - taCounts[i - 1]) // 2 for i in range(tCount)]
    # get midside coordinates for edges directly connecting neighbouring tubes, plus start/end cross point contributions
    mCoords = [[] for _ in range(tCount)]
    for ia in range(tCount):
        ib = (ia + 1) % tCount
        teCount = teCounts[ia]
        tnCount = teCount + 1
        tax, tad2 = tCoords[ia][0], tCoords[ia][2]
        tbx, tbd2 = tCoords[ib][0], tCoords[ib][2]
        for n in range(tnCount):
            aIn = (ia < inCount)
            an = n if aIn else -n
            ax = tax[an]
            ad2 = mult(tad2[an], 2.0 if aIn else -2.0)
            bIn = (ib < inCount)
            bn = -n if bIn else n
            bx = tbx[bn]
            bd2 = mult(tbd2[bn], -2.0 if bIn else 2.0)
            if not aIn:
                ax, ad2, bx, bd2 = bx, bd2, ax, ad2
            mx = interpolateCubicHermite(ax, ad2, bx, bd2, 0.5)
            md2 = mult(interpolateCubicHermiteDerivative(ax, ad2, bx, bd2, 0.5), 0.5)
            mCoords[ia].append([mx, [0.0, 0.0, 0.0], md2])

    # get cross points
    cCoords = []
    for ic in (0, -1):
        x = [sum(mCoords[it][ic][0][c] for it in range(tCount)) / tCount for c in range(3)]
        d1, d2 = getBifurcationCrossPoint(
            x,
            tCoords[0][ic][0], mult(tCoords[0][ic][2], -1.0),
            tCoords[1][ic][0], tCoords[1][ic][2],
            tCoords[2][ic][0], tCoords[2][ic][2])
        cCoords.append([x, d1, d2])

    # smooth around curves between cross points and create interior nodes
    for ia in range(tCount):
        teCount = teCounts[ia]
        tnCount = teCount + 1
        lx = [cCoords[0][0]] + [p[0] for p in mCoords[ia][1:-1]] + [cCoords[1][0]]
        c1d1 = [-d for d in cCoords[0][2]] if (ia == 0) else \
            [cCoords[0][1][c] + cCoords[0][2][c] for c in range(3)] if (ia == 1) else \
                [-d for d in cCoords[0][2]]
        c2d1 = cCoords[1][1] if (ia == 0) else \
            [-cCoords[1][1][c] - cCoords[1][2][c] for c in range(3)] if (ia == 1) else \
                cCoords[1][2]
        ld1 = [c1d1] + [p[1] for p in mCoords[ia][1:-1]] + [c2d1]
        ld1 = smoothCubicHermiteDerivativesLine(lx, ld1, fixStartDirection=True, fixEndDirection=True,
                                                magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
        for il in range(1, teCount):
            mCoords[ia][il][1] = ld1[il]

    return mCoords


def generateTube(outerTubeCoordinates, innerTubeCoordinates, elementsCountThroughWall,
                 region, fieldcache, coordinates: Field, nodeIdentifier, elementIdentifier,
                 startSkipCount: int=0, endSkipCount:int=0, startNodeIds: list=None, endNodeIds: list=None,
                 annotationMeshGroups=[], loop=False, serendipity=False):
    """
    Generate a 2D or 3D thick walled tube from supplied coordinates.
    Assumes client has active ChangeManager(fieldmodule).
    :param outerTubeCoordinates: Coordinates of outside of tube, or only coordinates for 2D.
    [ox, od1, od2, od12], each indexed by e.g. ox[along][around].
    :param innerTubeCoordinates: Coordinates of inside of tube, if 3D elements; like outer coordinates, or None if 2D
    :param elementsCountThroughWall: Number of elements through wall. Must be 1 if 2D.
    :param region: Zinc region to create model in.
    :param fieldcache: Field evaluation cache.
    :param coordinates: Finite element coordinate field to define.
    :param nodeIdentifier: First node identifier to use.
    :param elementIdentifier: First 2D element identifier to use.
    :param startSkipCount: Number of element rows to skip along at start.
    :param endSkipCount: Number of element rows to skip along at end.
    :param startNodeIds: Optional existing node identifiers [wall outward][around] to use at start, or None.
    :param endNodeIds: Optional existing node identifiers [wall outward][around] to use at end, or None.
    :param annotationMeshGroups: Mesh groups to add elements to.
    :param loop: Set to true to loop back to start coordinates at end.
    :param serendipity: True to use Hermite serendipity basis, False for regular Hermite with zero cross derivatives.
    :return: next node identifier, next element identifier, startNodeIds, endNodeIds
    """
    ox, od1, od2, od12 = outerTubeCoordinates
    ix, id1, id2, id12 = innerTubeCoordinates if innerTubeCoordinates else (None, None, None, None)
    dimension = 3 if innerTubeCoordinates else 2
    assert ((dimension == 2) and (elementsCountThroughWall == 1)) or \
           ((dimension == 3) and (elementsCountThroughWall >= 1))
    nodesCountThroughWall = (elementsCountThroughWall + 1) if (dimension == 3) else 1
    elementsCountAlong = len(ox) - 1
    elementsCountAround = len(ox[0])
    assert (not loop) or ((startSkipCount == 0) and (endSkipCount == 0))

    fieldmodule = region.getFieldmodule()

    # create nodes

    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if od12 and not serendipity:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

    tubeNodeIds = []
    oFactor = iFactor = 0.0
    for n2 in range(elementsCountAlong + 1):
        if (n2 < startSkipCount) or (n2 > elementsCountAlong - endSkipCount):
            tubeNodeIds.append(None)
            continue
        if startNodeIds and (n2 == startSkipCount):
            tubeNodeIds.append(startNodeIds)
            continue
        if endNodeIds and (n2 == (elementsCountAlong - endSkipCount)):
            tubeNodeIds.append(endNodeIds)
            continue
        if loop and (n2 == elementsCountAlong):
            tubeNodeIds.append(tubeNodeIds[0])
            continue
        tubeNodeIds.append([])
        for n3 in range(nodesCountThroughWall):
            ringNodeIds = []
            otx, otd1, otd2 = (ox[n2], od1[n2], od2[n2])
            otd12 = od12[n2] if (od12 and not serendipity) else None
            itx, itd1, itd2 = (ix[n2], id1[n2], id2[n2]) if innerTubeCoordinates else (None, None, None)
            itd12 = id12[n2] if (innerTubeCoordinates and otd12) else None
            if innerTubeCoordinates:
                oFactor = n3 / elementsCountThroughWall
                iFactor = 1.0 - oFactor
            for n1 in range(elementsCountAround):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                fieldcache.setNode(node)
                if (not innerTubeCoordinates) or (n3 == elementsCountThroughWall):
                    rx, rd1, rd2 = otx[n1], otd1[n1], otd2[n1]
                elif n3 == 0:
                    rx, rd1, rd2 = itx[n1], itd1[n1], itd2[n1]
                else:
                    rx = add(mult(otx[n1], oFactor), mult(itx[n1], iFactor))
                    rd1 = add(mult(otd1[n1], oFactor), mult(itd1[n1], iFactor))
                    rd2 = add(mult(otd2[n1], oFactor), mult(itd2[n1], iFactor))
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                if otd12:
                    rd12 = otd12[n1] if ((not innerTubeCoordinates) or (n3 == elementsCountThroughWall)) else \
                        itd12[n1] if (n3 == 0) else \
                        add(mult(otd12[n1], oFactor), mult(itd12[n1], iFactor))
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, rd12)
                ringNodeIds.append(nodeIdentifier)
                nodeIdentifier += 1
            tubeNodeIds[-1].append(ringNodeIds)

    # create elements

    mesh = fieldmodule.findMeshByDimension(dimension)
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE if (dimension == 3) else Element.SHAPE_TYPE_SQUARE)
    basis = fieldmodule.createElementbasis(
        dimension, (Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY if serendipity
            else Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE))
    if dimension == 3:
        basis.setFunctionType(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
    eft = mesh.createElementfieldtemplate(basis)
    if (not serendipity) and (not od12):
        # remove cross derivative terms for full bicubic Hermite
        ln = [1, 2, 3, 4] if (dimension == 2) else [1, 2, 3, 4, 5, 6, 7, 8]
        remapEftNodeValueLabel(eft, ln, Node.VALUE_LABEL_D2_DS1DS2, [])
    elementtemplate.defineField(coordinates, -1, eft)

    for e2 in range(startSkipCount, elementsCountAlong - endSkipCount):
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                e2p = e2 + 1
                e1p = (e1 + 1) % elementsCountAround
                nids = []
                for n3 in [e3, e3 + 1] if (dimension == 3) else [0]:
                    nids += [tubeNodeIds[e2][n3][e1], tubeNodeIds[e2][n3][e1p],
                             tubeNodeIds[e2p][n3][e1], tubeNodeIds[e2p][n3][e1p]]
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, nids)
                # print("Tube element", elementIdentifier, "nodes", nids)
                for annotationMeshGroup in annotationMeshGroups:
                    annotationMeshGroup.addElement(element)
                elementIdentifier += 1

    return nodeIdentifier, elementIdentifier, tubeNodeIds[startSkipCount], \
        tubeNodeIds[elementsCountAlong - endSkipCount]


def generateTubeBifurcation(outerTubeCoordinates, innerTubeCoordinates, inward, elementsCountThroughWall,
                            outerMidCoordinates, innerMidCoordinates, crossIndexes,
                            region, fieldcache, coordinates: Field,
                            nodeIdentifier, elementIdentifier, tubeNodeIds,
                            annotationMeshGroups, serendipity=False):
    """
    Generate a 2D tube bifurcation as elements connecting 3 rings of coordinates, optionally using existing nodes.
    Assumes client has active ChangeManager(fieldmodule).
    :param outerTubeCoordinates: List over 3 tubes (starting with "in" tubes) of outer coordinates
    [[ox], [od1], [od2], [od12]] for last/first tube rings in/out.
    :param innerTubeCoordinates: List over 3 tubes (starting with "in" tubes) of inner coordinates
    [[ix], [id1], [id2], [id12]] for last/first tube rings in/out, or None if 2D.
    :param inward: List over 3 tubes of True if inward, False if not. All inward tubes precede outward.
    :param elementsCountThroughWall: Number of elements through wall. Must be 1 if 2D.
    :param outerMidCoordinates: List of 3 middle half rings between tubes 1-2, 2-3 and 3-1 of outer coordinates
    [[omx], [omd1], [omd2]], from first to second cross point. First/last coordinates in each are the same with d1
    pointing towards tube 2 and d2 pointing towards tube 3 at first cross point, flipped on second cross point.
    :param innerMidCoordinates: List of 3 middle half rings between tubes 1-2, 2-3 and 3-1 of inner coordinates
    [[imx], [imd1], [imd2]], from first to second cross point. First/last coordinates in each are the same with d1
    pointing towards tube 2 and d2 pointing towards tube 3 at first cross point, flipped on second cross point.
    :param crossIndexes: List of 3 indexes around tube coordinates which link to first cross point.
    :param region: Zinc region to create model in.
    :param fieldcache: Field evaluation cache.
    :param coordinates: Finite element coordinate field to define.
    :param nodeIdentifier: First node identifier to use.
    :param elementIdentifier: First 2D element identifier to use.
    :param tubeNodeIds: List over 3 tubes of existing node identifiers [wall outward][around] to use at that inlet/outlet,
    any of which may be None. On return, None indexes are filled with new node identifiers.
    :param annotationMeshGroups: List over 3 tubes of lists of meshGroups to add elements to for the part of the
    bifurcation that is part of that segment.
    :param serendipity: True to use Hermite serendipity basis, False for regular Hermite with zero cross derivatives.
    :return: next node identifier, next element identifier
    """
    dimension = 3 if innerTubeCoordinates else 2
    assert ((dimension == 2) and (elementsCountThroughWall == 1)) or \
           ((dimension == 3) and (elementsCountThroughWall >= 1))
    nodesCountThroughWall = (elementsCountThroughWall + 1) if (dimension == 3) else 1

    fieldmodule = region.getFieldmodule()

    # create nodes

    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if (not serendipity) and (len(outerTubeCoordinates[0]) > 3):
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
    nodetemplateCross = nodes.createNodetemplate()
    nodetemplateCross.defineField(coordinates)
    nodetemplateCross.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplateCross.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplateCross.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

    midNodeIds = []
    oFactor = iFactor = 0.0
    for s in range(3):

        # ensure tube nodes are supplied or create here
        if not tubeNodeIds[s]:
            tubeNodeIds[s] = []
            otx, otd1, otd2 = outerTubeCoordinates[s][:3]
            otd12 = None if (serendipity or (len(outerTubeCoordinates[s]) == 3)) \
                else outerTubeCoordinates[s][3]
            itx, itd1, itd2 = innerTubeCoordinates[s][:3] if innerTubeCoordinates else (None, None, None)
            itd12 = None if ((not innerTubeCoordinates) or serendipity or (len(innerTubeCoordinates[s]) == 3)) \
                else innerTubeCoordinates[s][3]
            elementsCountAround = len(otx)
            for n3 in range(nodesCountThroughWall):
                ringNodeIds = []
                if innerTubeCoordinates:
                    oFactor = n3 / elementsCountThroughWall
                    iFactor = 1.0 - oFactor
                for n1 in range(elementsCountAround):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    fieldcache.setNode(node)
                    if (not innerTubeCoordinates) or (n3 == elementsCountThroughWall):
                        rx, rd1, rd2 = otx[n1], otd1[n1], otd2[n1]
                    elif n3 == 0:
                        rx, rd1, rd2 = itx[n1], itd1[n1], itd2[n1]
                    else:
                        rx = add(mult(otx[n1], oFactor), mult(itx[n1], iFactor))
                        rd1 = add(mult(otd1[n1], oFactor), mult(itd1[n1], iFactor))
                        rd2 = add(mult(otd2[n1], oFactor), mult(itd2[n1], iFactor))
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                    if otd12:
                        rd12 = otd12[n1] if ((not innerTubeCoordinates) or (n3 == elementsCountThroughWall)) else \
                            itd12[n1] if (n3 == 0) else add(mult(otd12[n1], oFactor), mult(itd12[n1], iFactor))
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, rd12)
                    ringNodeIds.append(nodeIdentifier)
                    nodeIdentifier += 1
                tubeNodeIds[s].append(ringNodeIds)

        # create nodes around middle half rings
        midNodeIds.append([])
        omx, omd1, omd2 = outerMidCoordinates[s][:3]
        omd12 = None if (serendipity or (len(outerMidCoordinates[s]) == 3)) else outerMidCoordinates[s][3]
        imx, imd1, imd2 = innerMidCoordinates[s][:3] if innerMidCoordinates else (None, None, None)
        imd12 = None if ((not innerMidCoordinates) or serendipity or (len(innerMidCoordinates[s]) == 3)) \
            else innerMidCoordinates[s][3]
        elementsCountAroundHalf = len(omx) - 1
        for n3 in range(nodesCountThroughWall):
            ringNodeIds = []
            if innerTubeCoordinates:
                oFactor = n3 / elementsCountThroughWall
                iFactor = 1.0 - oFactor
            for n1 in range(elementsCountAroundHalf + 1):
                cross1 = n1 == 0
                cross2 = n1 == elementsCountAroundHalf
                if s > 0:
                    if cross1:
                        ringNodeIds.append(midNodeIds[0][n3][0])
                        continue
                    if cross2:
                        ringNodeIds.append(midNodeIds[0][n3][-1])
                        continue
                node = nodes.createNode(nodeIdentifier, nodetemplateCross if (cross1 or cross2) else nodetemplate)
                fieldcache.setNode(node)
                if (not innerMidCoordinates) or (n3 == elementsCountThroughWall):
                    rx, rd1, rd2 = omx[n1], omd1[n1], omd2[n1]
                elif n3 == 0:
                    rx, rd1, rd2 = imx[n1], imd1[n1], imd2[n1]
                else:
                    rx = add(mult(omx[n1], oFactor), mult(imx[n1], iFactor))
                    rd1 = add(mult(omd1[n1], oFactor), mult(imd1[n1], iFactor))
                    rd2 = add(mult(omd2[n1], oFactor), mult(imd2[n1], iFactor))
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                if omd12 and not (cross1 or cross2):
                    rd12 = omd12[n1] if ((not innerMidCoordinates) or (n3 == elementsCountThroughWall)) else \
                        imd12[n1] if (n3 == 0) else add(mult(omd12[n1], oFactor), mult(imd12[n1], iFactor))
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, rd12)
                ringNodeIds.append(nodeIdentifier)
                nodeIdentifier = nodeIdentifier + 1

            midNodeIds[s].append(ringNodeIds)

    # create elements

    mesh = fieldmodule.findMeshByDimension(dimension)
    elementtemplateStd = mesh.createElementtemplate()
    elementtemplateMidInward = mesh.createElementtemplate()
    elementtemplateMidOutward = mesh.createElementtemplate()
    elementtemplateCross = mesh.createElementtemplate()
    for et in [elementtemplateStd, elementtemplateMidInward, elementtemplateMidOutward, elementtemplateCross]:
        et.setElementShapeType(Element.SHAPE_TYPE_CUBE if (dimension == 3) else Element.SHAPE_TYPE_SQUARE)
    basis = fieldmodule.createElementbasis(
        dimension, (Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY if serendipity
            else Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE))
    if dimension == 3:
        basis.setFunctionType(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
    eftStd = mesh.createElementfieldtemplate(basis)
    eftMidInward = mesh.createElementfieldtemplate(basis)
    eftMidOutward = mesh.createElementfieldtemplate(basis)
    if (not serendipity) and (len(outerTubeCoordinates[0]) < 4):
        # remove cross derivative terms for full bicubic Hermite
        ln = [1, 2, 3, 4] if (dimension == 2) else [1, 2, 3, 4, 5, 6, 7, 8]
        for eft in [eftStd, eftMidInward, eftMidOutward]:
            remapEftNodeValueLabel(eft, ln, Node.VALUE_LABEL_D2_DS1DS2, [])
    elementtemplateStd.defineField(coordinates, -1, eftStd)
    setEftScaleFactorIds(eftMidInward, [1], [])
    scaleEftNodeValueLabels(eftMidInward, [3, 4] if (dimension == 2) else [3, 4, 7, 8],
                            [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
    elementtemplateMidInward.defineField(coordinates, -1, eftMidInward)
    setEftScaleFactorIds(eftMidOutward, [1], [])
    scaleEftNodeValueLabels(eftMidOutward, [1, 2] if (dimension == 2) else [1, 2, 5, 6],
                            [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
    elementtemplateMidOutward.defineField(coordinates, -1, eftMidOutward)

    eftCrossForward = []
    eftCrossForwardScalefactors = []
    for s in range(3):
        eftCrossForward.append([])
        eftCrossForwardScalefactors.append([])
        for ce in range(2):  # cross end 0 or 1
            eft = mesh.createElementfieldtemplate(basis)
            if inward[s]:
                lnCross = [3] if (ce == 0) else [4]
                lnOther = [4] if (ce == 0) else [3]
                lnTube = [1, 2]
                lnMid = [3, 4]
            else:
                lnCross = [1] if (ce == 0) else [2]
                lnOther = [2] if (ce == 0) else [1]
                lnTube = [3, 4]
                lnMid = [1, 2]
            if dimension == 3:
                lnCross.append(lnCross[0] + 4)
                lnOther.append(lnOther[0] + 4)
                for i in range(2):
                    lnTube.append(lnTube[i] + 4)
                    lnMid.append(lnMid[i] + 4)
            if not serendipity:
                # remove cross derivative terms for full bicubic Hermite where not supplied
                if len(outerTubeCoordinates[0]) < 4:
                    remapEftNodeValueLabel(eft, lnTube, Node.VALUE_LABEL_D2_DS1DS2, [])
                if len(outerMidCoordinates[0]) < 4:
                    remapEftNodeValueLabel(eft, lnMid, Node.VALUE_LABEL_D2_DS1DS2, [])
                else:
                    remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D2_DS1DS2, [])
            if (not inward[s]) or ((s != 1) and (ce == 0)) or ((s == 1) and (inward[s] or (ce != 0))) or \
                    ((s == 2) and inward[s] and (ce != 0)):
                setEftScaleFactorIds(eft, [1], [])
                scalefactors = [-1.0]
            else:
                scalefactors = None
            eftCrossForwardScalefactors[s].append(scalefactors)
            if s == 0:
                if ce == 0:
                    if inward[s]:
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                    else:
                        scaleEftNodeValueLabels(eft, lnCross, [Node.VALUE_LABEL_D_DS1], [1])
                scaling = [] if inward[s] else [1]
                remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS2,
                                       [(Node.VALUE_LABEL_D_DS1, scaling), (Node.VALUE_LABEL_D_DS2, scaling)])
                if (ce != 0) and not inward[s]:
                    remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS2, [])])
            elif s == 1:
                scaling = [] if (ce == 0) else [1]
                remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                       [(Node.VALUE_LABEL_D_DS1, scaling), (Node.VALUE_LABEL_D_DS2, scaling)])
                if (ce == 0) and inward[s]:
                    remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                if ce != 0:
                    if inward[s]:
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                    else:
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS1, [])])
            else:  # s == 2:
                if ce == 0:
                    if inward[s]:
                        remapEftNodeValueLabel(
                            eft, lnCross, Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(
                            eft, lnCross, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                    else:
                        remapEftNodeValueLabel(
                            eft, lnCross, Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                        remapEftNodeValueLabel(
                            eft, lnCross, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                elif inward[s]:
                    remapEftNodeValueLabel(
                        eft, lnCross, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(
                        eft, lnCross, Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
            if not inward[s]:
                scaleEftNodeValueLabels(eft, lnOther, [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
            eftCrossForward[s].append(eft)

    eftCrossReverse = []
    eftCrossReverseScalefactors = []
    for s in range(3):
        eftCrossReverse.append([])
        eftCrossReverseScalefactors.append([])
        for ce in range(2):  # cross end 0 or 1
            eft = mesh.createElementfieldtemplate(basis)
            if inward[s]:
                lnCross = [3] if (ce == 0) else [4]
                lnOther = [4] if (ce == 0) else [3]
                lnTube = [1, 2]
                lnMid = [3, 4]
            else:
                lnCross = [1] if (ce == 0) else [2]
                lnOther = [2] if (ce == 0) else [1]
                lnTube = [3, 4]
                lnMid = [1, 2]
            if dimension == 3:
                lnCross.append(lnCross[0] + 4)
                lnOther.append(lnOther[0] + 4)
                for i in range(2):
                    lnTube.append(lnTube[i] + 4)
                    lnMid.append(lnMid[i] + 4)
            if not serendipity:
                # remove cross derivative terms for full bicubic Hermite where not supplied
                if len(outerTubeCoordinates[0]) < 4:
                    remapEftNodeValueLabel(eft, lnTube, Node.VALUE_LABEL_D2_DS1DS2, [])
                if len(outerMidCoordinates[0]) < 4:
                    remapEftNodeValueLabel(eft, lnMid, Node.VALUE_LABEL_D2_DS1DS2, [])
                else:
                    remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D2_DS1DS2, [])
            if inward[s] or (s == 0) or ((s != 2) and (ce == 0)) or ((s == 2) and (ce != 0)):
                setEftScaleFactorIds(eft, [1], [])
                scalefactors = [-1.0]
            else:
                scalefactors = None
            eftCrossReverseScalefactors[s].append(scalefactors)
            if s == 0:
                if ce == 0:
                    if inward[s]:
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                    else:
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                scaling = [] if inward[s] else [1]
                remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS2,
                                       [(Node.VALUE_LABEL_D_DS1, scaling), (Node.VALUE_LABEL_D_DS2, scaling)])
                if (ce != 0) and not inward[s]:
                    remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
            elif s == 1:
                if ce == 0:
                    if inward[s]:
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                    else:
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                else:
                    if inward[s]:
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS2, [])])
            else:  # s == 2:
                scaling = [] if (ce == 0) else [1]
                remapEftNodeValueLabel(eft, lnCross, Node.VALUE_LABEL_D_DS1,
                                       [(Node.VALUE_LABEL_D_DS1, scaling), (Node.VALUE_LABEL_D_DS2, scaling)])
                if ce == 0:
                    if inward[s]:
                        remapEftNodeValueLabel(
                            eft, lnCross, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                else:
                    if inward[s]:
                        remapEftNodeValueLabel(
                            eft, lnCross, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                    else:
                        remapEftNodeValueLabel(
                            eft, lnCross, Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
            if inward[s]:
                scaleEftNodeValueLabels(eft, lnOther, [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
            eftCrossReverse[s].append(eft)

    aroundCounts = [len(outerTubeCoordinates[s][0]) for s in range(3)]
    connectionCounts = get_tube_bifurcation_connection_elements_counts(aroundCounts)

    for s in range(3):

        # forward connections

        eftMid = eftStd
        elementtemplateMid = elementtemplateStd
        scalefactorsMid = None
        if not inward[s]:
            eftMid = eftMidOutward
            elementtemplateMid = elementtemplateMidOutward
            scalefactorsMid = [-1.0]

        for e3 in range(elementsCountThroughWall):
            for e1 in range(connectionCounts[s]):
                eft = eftMid
                elementtemplate = elementtemplateMid
                scalefactors = scalefactorsMid

                if e1 == 0:
                    eft = eftCrossForward[s][0]
                    elementtemplateCross.defineField(coordinates, -1, eft)
                    elementtemplate = elementtemplateCross
                    scalefactors = eftCrossForwardScalefactors[s][0]
                elif e1 == (connectionCounts[s] - 1):
                    eft = eftCrossForward[s][1]
                    elementtemplateCross.defineField(coordinates, -1, eft)
                    elementtemplate = elementtemplateCross
                    scalefactors = eftCrossForwardScalefactors[s][1]

                nids = []
                for n3 in ([0] if (dimension == 2) else [e3, e3 + 1]):
                    if inward[s]:
                        nStart = crossIndexes[s] - aroundCounts[s]
                        nids += [tubeNodeIds[s][n3][nStart + e1], tubeNodeIds[s][n3][nStart + e1 + 1],
                                 midNodeIds[s][n3][e1], midNodeIds[s][n3][e1 + 1]]
                    else:
                        nStart = crossIndexes[s] - connectionCounts[s]
                        re1 = connectionCounts[s] - e1
                        nids += [midNodeIds[s][n3][re1], midNodeIds[s][n3][re1 - 1],
                                 tubeNodeIds[s][n3][nStart + e1], tubeNodeIds[s][n3][nStart + e1 + 1]]

                element = mesh.createElement(elementIdentifier, elementtemplate)
                result1 = element.setNodesByIdentifier(eft, nids)
                result2 = element.setScaleFactors(eft, scalefactors) if scalefactors else "-"
                # print('create element tube bifurcation forward s', s, elementIdentifier, element.isValid(), result1, result2, nids)
                for meshGroup in annotationMeshGroups[s]:
                    meshGroup.addElement(element)
                elementIdentifier += 1

        # reverse connections

        eftMid = eftStd
        elementtemplateMid = elementtemplateStd
        scalefactorsMid = None
        if inward[s]:
            eftMid = eftMidInward
            elementtemplateMid = elementtemplateMidInward
            scalefactorsMid = [-1.0]

        for e3 in range(elementsCountThroughWall):
            for e1 in range(connectionCounts[s - 1]):
                eft = eftMid
                elementtemplate = elementtemplateMid
                scalefactors = scalefactorsMid

                if e1 == 0:
                    eft = eftCrossReverse[s][0]
                    elementtemplateCross.defineField(coordinates, -1, eft)
                    elementtemplate = elementtemplateCross
                    scalefactors = eftCrossReverseScalefactors[s][0]
                elif e1 == (connectionCounts[s - 1] - 1):
                    eft = eftCrossReverse[s][1]
                    elementtemplateCross.defineField(coordinates, -1, eft)
                    elementtemplate = elementtemplateCross
                    scalefactors = eftCrossReverseScalefactors[s][1]

                nids = []
                for n3 in ([0] if (dimension == 2) else [e3, e3 + 1]):
                    if inward[s]:
                        nStart = crossIndexes[s] - connectionCounts[s - 1]
                        re1 = connectionCounts[s - 1] - e1
                        nids += [tubeNodeIds[s][n3][nStart + e1], tubeNodeIds[s][n3][nStart + e1 + 1],
                                 midNodeIds[s - 1][n3][re1], midNodeIds[s - 1][n3][re1 - 1]]
                    else:
                        nStart = crossIndexes[s] - aroundCounts[s]
                        nids += [midNodeIds[s - 1][n3][e1], midNodeIds[s - 1][n3][e1 + 1],
                                 tubeNodeIds[s][n3][nStart + e1], tubeNodeIds[s][n3][nStart + e1 + 1]]

                element = mesh.createElement(elementIdentifier, elementtemplate)
                result1 = element.setNodesByIdentifier(eft, nids)
                result2 = element.setScaleFactors(eft, scalefactors) if scalefactors else "-"
                # print('Tube bifurcation element reverse s', s, elementIdentifier, element.isValid(), result1, result2, nids)
                for meshGroup in annotationMeshGroups[s]:
                    meshGroup.addElement(element)
                elementIdentifier += 1

    return nodeIdentifier, elementIdentifier


class SegmentTubeData:

    def __init__(self, pathParameters, elementsCountAround):
        self._pathParameters = pathParameters
        self._elementsCountAround = elementsCountAround
        self._segmentLength = getCubicHermiteCurvesLength(pathParameters[0], pathParameters[1])
        self._rawTubeCoordinates = None
        self._rawTrackSurface = None
        self._sampledTubeCoordinates = None
        self._sampledTrackSurface = None
        self._sampledNodeIds = []  # indexed along sample nodes. Only ever set on rows where junctions occur
        self._annotationMeshGroups = []

    def getPathParameters(self):
        return self._pathParameters

    def getElementsCountAround(self):
        return self._elementsCountAround

    def getSegmentLength(self):
        return self._segmentLength

    def getRawTrackSurface(self):
        """
        Available after calling setRawTubeCoordinates().
        :return: TrackSurface
        """
        return self._rawTrackSurface

    def getRawTubeCoordinates(self):
        return self._rawTubeCoordinates

    def setRawTubeCoordinates(self, rawTubeCoordinates):
        """
        Set raw tube coordinates at network layout spacing.
        Creates TrackSurface internally.
        :param rawTubeCoordinates: px, pd1, pd2, pd12
        """
        assert len(rawTubeCoordinates[0][0]) == self._elementsCountAround
        self._rawTubeCoordinates = rawTubeCoordinates
        px, pd1, pd2, pd12 = rawTubeCoordinates
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(px)):
            nx += px[i]
            nd1 += pd1[i]
            nd2 += pd2[i]
            nd12 += pd12[i]
        self._rawTrackSurface = TrackSurface(len(px[0]), len(px) - 1, nx, nd1, nd2, nd12, loop1=True)

    def getSampledTrackSurface(self):
        """
        Available after calling setSampledTubeCoordinates().
        :return: TrackSurface
        """
        return self._sampledTrackSurface

    def getSampledTubeCoordinates(self):
        return self._sampledTubeCoordinates

    def setSampledTubeCoordinates(self, sampledTubeCoordinates):
        """
        :param sampledTubeCoordinates: sx, sd1, sd2, sd12
        """
        assert len(sampledTubeCoordinates[0][0]) == self._elementsCountAround
        self._sampledTubeCoordinates = sampledTubeCoordinates
        self._sampledNodeIds = [None] * len(self._sampledTubeCoordinates[0])
        px, pd1, pd2, pd12 = sampledTubeCoordinates
        nx = []
        nd1 = []
        nd2 = []
        nd12 = []
        for i in range(len(px)):
            nx += px[i]
            nd1 += pd1[i]
            nd2 += pd2[i]
            nd12 += pd12[i]
        self._sampledTrackSurface = TrackSurface(len(px[0]), len(px) - 1, nx, nd1, nd2, nd12, loop1=True)

    def getSampledElementsCountAlong(self):
        """
        Must have previously called setSampledTubeCoordinates
        :return: Number of elements along sampled tube.
        """
        return len(self._sampledTubeCoordinates[0]) - 1

    def getStartNodeIds(self, startSkipCount):
        """
        Get start nodes for supplying to adjacent tube or bifurcation.
        :param startSkipCount: Row in from start that node ids are for.
        """
        return self._sampledNodeIds[startSkipCount]

    def setStartNodeIds(self, startNodeIds, startSkipCount):
        """
        :param startNodeIds: list of node ids around tube start row.
        :param startSkipCount: Row in from start that ids are for.
        """
        self._sampledNodeIds[startSkipCount] = startNodeIds

    def getEndNodeIds(self, endSkipCount):
        """
        Get end nodes for supplying to adjacent tube or bifurcation.
        :param endSkipCount: Row in from end that node ids are for.
        """
        return self._sampledNodeIds[self.getSampledElementsCountAlong() - endSkipCount]

    def setEndNodeIds(self, endNodeIds, endSkipCount):
        """
        :param endNodeIds: list of node ids around tube end row.
        :param endSkipCount: Row in from end that node ids are for.
        """
        self._sampledNodeIds[self.getSampledElementsCountAlong() - endSkipCount] = endNodeIds

    def addAnnotationMeshGroup(self, annotationMeshGroup):
        """
        Add an annotation mesh group for segment elements to be added to.
        :param annotationMeshGroup: Mesh group to add.
        """
        self._annotationMeshGroups.append(annotationMeshGroup)

    def getAnnotationMeshGroups(self):
        return self._annotationMeshGroups


class TubeBifurcationData:
    """
    Describes junction between three segments.
    Used to get intersection curves and points between them, and trim surfaces for segments.
    """

    def __init__(self, networkSegmentsIn: list, networkSegmentsOut: list, segmentTubeData):
        """
        :param networkSegmentsIn: List of input segments.
        :param networkSegmentsOut: List of output segments.
        :param segmentTubeData: dict NetworkSegment -> SegmentTubeData.
        """
        self._networkSegmentsIn = networkSegmentsIn
        self._networkSegmentsOut = networkSegmentsOut
        self._networkSegments = networkSegmentsIn + networkSegmentsOut
        segmentCount = len(self._networkSegments)
        assert segmentCount == 3
        self._tubeData = [segmentTubeData[networkSegment] for networkSegment in self._networkSegments]
        self._segmentsIn = [self._networkSegments[s] in self._networkSegmentsIn for s in range(3)]
        # following are calculated in determineCrossIndexes()
        self._connectingCoordinateRings = [[]] * 3  # second row of coordinates from end, made into nodes
        self._endCoordinateRings = [[]] * 3  # row of coordinates at end, NOT made into nodes
        self._aroundCounts = [0] * 3
        self._connectionCounts = [0] * 3
        self._aCrossIndexes = None
        self._bCrossIndexes = None
        self._intersectionCurves = []
        self._trimSurfaces = []
        # self._calculateTrimSurfacesFromIntersectionCurves()
        self._calculateTrimSurfacesNew()

    def _calculateTrimSurfacesFromIntersectionCurves(self):
        assert (not self._intersectionCurves) and (not self._trimSurfaces)
        segmentCount = len(self._networkSegments)
        assert segmentCount == 3

        # get intersection curves between pairs of segments
        for s in range(segmentCount):
            tubeData1 = self._tubeData[s]
            tubeTrackSurface1 = tubeData1.getRawTrackSurface()
            tubeData2 = self._tubeData[(s + 1) % segmentCount]
            tubeTrackSurface2 = tubeData2.getRawTrackSurface()
            cx, cd1, cProportions, loop = tubeTrackSurface1.findIntersectionCurve(tubeTrackSurface2)
            self._intersectionCurves.append((cx, cd1, cProportions, loop))

        # get trim surfaces
        for s in range(segmentCount):
            networkSegment = self._networkSegments[s]
            tubeData = self._tubeData[s]
            pathParameters = tubeData.getPathParameters()
            path_d1 = pathParameters[1]
            along = path_d1[-1 if networkSegment in self._networkSegmentsIn else 0]
            ax, ad1, _, aloop = self._intersectionCurves[s - 1]
            bx, bd1, _, bloop = self._intersectionCurves[s]
            trimSurface = None
            if ax and bx:
                d2a = None
                cax = None
                ha = len(ax) // 2  # halfway index in a
                assert ha > 0
                if not aloop:
                    # get coordinates and 2nd derivative halfway along a
                    assert ha * 2 == len(ax) - 1  # must be an even number of elements
                    d2am = interpolateCubicHermiteSecondDerivative(ax[ha - 1], ad1[ha - 1], ax[ha], ad1[ha], 1.0)
                    d2ap = interpolateCubicHermiteSecondDerivative(ax[ha], ad1[ha], ax[ha + 1], ad1[ha + 1], 0.0)
                    d2a = add(d2am, d2ap)
                    cax = ax[4]
                d2b = None
                cbx = None
                hb = len(bx) // 2  # halfway index in b
                assert hb > 0
                if not bloop:
                    # get coordinates and 2nd derivative halfway along b
                    assert hb * 2 == len(bx) - 1  # must be an even number of elements
                    d2bm = interpolateCubicHermiteSecondDerivative(bx[hb - 1], bd1[hb - 1], bx[hb], bd1[hb], 1.0)
                    d2bp = interpolateCubicHermiteSecondDerivative(bx[hb], bd1[hb], bx[hb + 1], bd1[hb + 1], 0.0)
                    d2b = add(d2bm, d2bp)
                    cbx = bx[4]
                aNearestLocation = None  # (element index, xi)
                bNearestLocation = None  # (element index, xi)
                if aloop and bloop:
                    aNearestLocation, bNearestLocation = \
                        getNearestLocationBetweenCurves(ax, ad1, bx, bd1, aloop, bloop)[:2]
                else:
                    if aloop:
                        aNearestLocation = getNearestLocationOnCurve(ax, ad1, cbx, aloop)[0]
                    if bloop:
                        bNearestLocation = getNearestLocationOnCurve(bx, bd1, cax, bloop)[0]
                if aloop and aNearestLocation:
                    # get coordinates and 2nd derivative at opposite to nearest location on a
                    assert ha * 2 == len(ax)  # must be an even number of elements
                    oa = aNearestLocation[0] - ha  # opposite element index in a
                    xi = aNearestLocation[1]
                    d2a = interpolateCubicHermiteSecondDerivative(ax[oa], ad1[oa], ax[oa + 1], ad1[oa + 1], xi)
                    cax = interpolateCubicHermite(ax[oa], ad1[oa], ax[oa + 1], ad1[oa + 1], xi)
                if bloop and bNearestLocation:
                    # get coordinates and 2nd derivative at opposite to nearest location on b
                    assert hb * 2 == len(bx)  # must be an even number of elements
                    ob = bNearestLocation[0] - hb  # opposite element index in b
                    xi = bNearestLocation[1]
                    d2b = interpolateCubicHermiteSecondDerivative(bx[ob], bd1[ob], bx[ob + 1], bd1[ob + 1], xi)
                    cbx = interpolateCubicHermite(bx[ob], bd1[ob], bx[ob + 1], bd1[ob + 1], xi)
                if cax and cbx:
                    d2b = [-s for s in d2b]  # reverse derivative on second point
                    across = sub(cbx, cax)
                    sizeAcross = magnitude(across)
                    normAcross = normalize(across)
                    up = cross(across, along)
                    cad2 = normalize(ad1[4])
                    if dot(up, cad2) < 0.0:
                        cad2 = [-d for d in cad2]
                    cbd2 = normalize(bd1[4])
                    if dot(up, cbd2) < 0.0:
                        cbd2 = [-d for d in cbd2]
                    # move coordinates out slightly to guarantee intersection
                    arcLength = computeCubicHermiteArcLength(cax, d2a, cbx, d2b, rescaleDerivatives=True)
                    d2a = mult(d2a, arcLength / magnitude(d2a))
                    # cax = sub(cax, mult(d2a, 0.05))
                    d2b = mult(d2b, arcLength / magnitude(d2b))
                    # cbx = add(cbx, mult(d2b, 0.05))
                    nx = []
                    nd1 = []
                    nd2 = []
                    elementsCount2 = 2
                    for i in range(elementsCount2 + 1):
                        delta = (i - elementsCount2 // 2) * 1.2 * sizeAcross
                        vax = add(cax, mult(cad2, delta))
                        vbx = add(cbx, mult(cbd2, delta))
                        deltaAcross = sub(vbx, vax)
                        normDeltaAcross = normalize(deltaAcross)
                        use_d2a = d2a
                        use_d2b = d2b
                        if dot(normAcross, normDeltaAcross) <= 0.5:
                            use_d2a = deltaAcross
                            use_d2b = deltaAcross
                        arcLength = 1.5 * computeCubicHermiteArcLength(vax, use_d2a, vbx, use_d2b, rescaleDerivatives=True)
                        vad1 = mult(use_d2a, arcLength / magnitude(use_d2a))
                        vbd1 = mult(use_d2b, arcLength / magnitude(use_d2b))
                        vad2 = mult(cad2, sizeAcross)
                        vbd2 = mult(cbd2, sizeAcross)
                        # add extensions before / after
                        nx.append(sub(vax, vad1))
                        nx.append(vax)
                        nd1.append(vad1)
                        nd1.append(vad1)
                        nd2.append(vad2)
                        nd2.append(vad2)
                        nx.append(vbx)
                        nx.append(add(vbx, vbd1))
                        nd1.append(vbd1)
                        nd1.append(vbd1)
                        nd2.append(vbd2)
                        nd2.append(vbd2)
                    # smooth extension d2
                    for j in [0, 3]:
                        tx = [nx[i * 4 + j] for i in range(3)]
                        td2 = smoothCubicHermiteDerivativesLine(tx, [nd2[i * 4 + j] for i in range(3)])
                        for i in range(3):
                            nd2[i * 4 + j] = td2[i]
                    trimSurface = TrackSurface(3, elementsCount2, nx, nd1, nd2)
            elif ax:
                # use previous tube surface as trim surface
                trimSurface = self._tubeData[s - 1].getRawTrackSurface()
            elif bx:
                # use next tube surface as trim surface
                trimSurface = self._tubeData[s - 2].getRawTrackSurface()
            self._trimSurfaces.append(trimSurface)

    def _calculateTrimSurfacesNew(self):
        assert not self._trimSurfaces
        segmentCount = len(self._networkSegments)
        assert segmentCount == 3

        for s in range(segmentCount):
            self._intersectionCurves.append(None)

        dirEnd = []
        for s in range(segmentCount):
            tubeData = self._tubeData[s]
            pathParameters = tubeData.getPathParameters()
            endIndex = -1 if self._segmentsIn[s] else 0
            dirEnd.append(normalize(pathParameters[1][endIndex]))

        for s in range(segmentCount):
            networkSegment = self._networkSegments[s]
            tubeData = self._tubeData[s]
            pathParameters = tubeData.getPathParameters()
            endIndex = -1 if self._segmentsIn[s] else 0
            # d1End = pathParameters[1][endIndex]
            d2End = pathParameters[2][endIndex]
            d3End = pathParameters[4][endIndex]
            # get index of other direction most normal to this segment direction
            sos = ((s - 1), (s - 2)) if (abs(dot(dirEnd[s], dirEnd[s - 1])) < abs(dot(dirEnd[s], dirEnd[s - 2]))) \
                else ((s - 2), (s - 1))
            xEnd = pathParameters[0][endIndex]
            endEllipseNormal = normalize(cross(d2End, d3End))
            # get components of so1 direction aligned with d2, d3 and compute initial phase
            oDirEnd = dirEnd[sos[0]]
            if self._segmentsIn[sos[0]]:
                oDirEnd = [-d for d in oDirEnd]
            dx = dot(oDirEnd, d2End)
            dy = dot(oDirEnd, d3End)
            phaseAngle = math.atan2(dy, dx)
            # if s == 0:
            #     print("oDirEnd", oDirEnd, "d2", dx, "d3", dy, "phaseAngle", math.degrees(phaseAngle))
            # GRC Future: get exact phase angle with non-linear optimisation!
            elementsCountAround = 6
            tubeCoordinates = getPathRawTubeCoordinates(
                pathParameters, elementsCountAround, radius=1.0, phaseAngle=phaseAngle)
            px, pd1, pd2, pd12 = tubeCoordinates
            nx = []
            nd1 = []
            nd2 = []
            nd12 = []
            for i in range(len(px)):
                nx += px[i]
                nd1 += pd1[i]
                nd2 += pd2[i]
                nd12 += pd12[i]
            tmpTrackSurface = TrackSurface(len(px[0]), len(px) - 1, nx, nd1, nd2, nd12, loop1=True)
            pointsCountAlong = len(pathParameters[0])
            # get coordinates and directions of intersection points of longitudinal lines and other track surfaces
            rx = []
            rd1 = []
            trim = False
            # if s == 0:
            #     print("sos", (sos[0] % 3) + 1, (sos[1] % 3) + 1)
            trimIndex = 0
            for n1 in range(elementsCountAround):
                proportion1 = n1 / elementsCountAround
                cx = [tubeCoordinates[0][n2][n1] for n2 in range(pointsCountAlong)]
                cd2 = [tubeCoordinates[2][n2][n1] for n2 in range(pointsCountAlong)]
                ox = x = tubeCoordinates[0][endIndex][n1]
                d1 = tubeCoordinates[1][endIndex][n1]
                maxProportionFromEnd = 0.0
                for so in sos:
                    otherTrackSurface = self._tubeData[so].getRawTrackSurface()
                    otherSurfacePosition, curveLocation, isIntersection = \
                        otherTrackSurface.findNearestPositionOnCurve(
                            cx, cd2, loop=False, sampleEnds=False, sampleHalf=2 if self._segmentsIn[s] else 1)
                    if isIntersection:
                        proportion2 = (curveLocation[0] + curveLocation[1]) / (pointsCountAlong - 1)
                        proportionFromEnd = abs(proportion2 - (1.0 if self._segmentsIn[s] else 0.0))
                        if proportionFromEnd > maxProportionFromEnd:
                            trim = True
                            trimIndex = (so % 3) + 1
                            surfacePosition = tmpTrackSurface.createPositionProportion(proportion1, proportion2)
                            x, d1, d2 = tmpTrackSurface.evaluateCoordinates(surfacePosition, derivatives=True)
                            n = cross(d1, d2)
                            ox, od1, od2 = otherTrackSurface.evaluateCoordinates(
                                otherSurfacePosition, derivatives=True)
                            on = cross(od1, od2)
                            d1 = cross(n, on)
                            maxProportionFromEnd = proportionFromEnd
                # if s == 0:
                #     print("    ", n1, "trim", trimIndex, "max", maxProportionFromEnd, "x", x, "d1", d1, sub(ox, x))
                # ensure d1 directions go around in same direction as loop
                if dot(endEllipseNormal, cross(sub(x, xEnd), d1)) < 0.0:
                    d1 = [-d for d in d1]
                rx.append(x)
                rd1.append(d1)
            if trim:
                rd1 = smoothCubicHermiteDerivativesLoop(rx, rd1, fixAllDirections=True,
                                                        magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
                rd2 = [sub(rx[n1], xEnd) for n1 in range(elementsCountAround)]
                rd12 = smoothCurveSideCrossDerivatives(rx, rd1, [rd2], loop=True)[0]
                nx = []
                nd1 = []
                nd2 = []
                nd12 = []
                for factor in (0.5, 1.5):
                    for n1 in range(elementsCountAround):
                        d2 = sub(rx[n1], xEnd)
                        x = add(xEnd, mult(d2, factor))
                        d1 = mult(rd1[n1], factor)
                        d12 = mult(rd12[n1], factor)
                        nx.append(x)
                        nd1.append(d1)
                        nd2.append(d2)
                        nd12.append(d12)
                trimSurface = TrackSurface(elementsCountAround, 1, nx, nd1, nd2, nd12, loop1=True)
                self._trimSurfaces.append(trimSurface)
            else:
                self._trimSurfaces.append(None)

    def getIntersectionCurve(self, s):
        """
        :param s: Index from 0 to 2
        :return: (cx, cd1, cProportions on TrackSurface s, loop)
        """
        return self._intersectionCurves[s]

    def getCrossIndexes(self):
        return self._aCrossIndexes

    def getMidCoordinates(self):
        return self._midCoordinates

    def getSegmentsIn(self):
        return self._segmentsIn

    def getSegmentTrimSurface(self, networkSegment):
        """
        :return: Trim TrackSurface
        """
        return self.getTrimSurface(self._networkSegments.index(networkSegment))

    def getTrimSurface(self, s):
        """
        :param s: Index from 0 to 2
        :return: Trim TrackSurface
        """
        return self._trimSurfaces[s]

    def getConnectingTubeCoordinates(self):
        """
        Get ring of coordinates of attached tubes -- 1 row in from end sampled coordinate rings.
        :return: List[4] of coordinates around ring: 4 values are x, d1, d2, d12.
        """
        return self._connectingCoordinateRings

    def getTubeData(self):
        return self._tubeData

    def _getConnectionCounts(self):
        for s in range(3):
            tubeData = self._tubeData[s]
            row = -2 if self._segmentsIn[s] else 1
            sampledTubeCoordinates = tubeData.getSampledTubeCoordinates()
            self._connectingCoordinateRings[s] = [value[row] for value in sampledTubeCoordinates]
            endRow = -1 if self._segmentsIn[s] else 0
            self._endCoordinateRings[s] = [value[endRow] for value in sampledTubeCoordinates]
            self._aroundCounts[s] = len(self._connectingCoordinateRings[s][0])
        self._connectionCounts = get_tube_bifurcation_connection_elements_counts(self._aroundCounts)

    def copyCrossIndexes(self, sourceTubeBifurcationData):
        assert self._aCrossIndexes is None  # should only call once
        self._getConnectionCounts()
        self._aCrossIndexes = sourceTubeBifurcationData._aCrossIndexes
        self._bCrossIndexes = sourceTubeBifurcationData._bCrossIndexes

    def determineCrossIndexes(self):
        """
        Determine the node indexes around the tubes which give the least distorted 3-way crossing point.
        Only call for one of inner or outer bifurcation data (outer is recommendd). For the other (inner)
        use copyCrossIndexes(outerTubeBifurcationData)
        """
        assert self._aCrossIndexes is None  # should only call once
        self._getConnectionCounts()

        # get cross points where there is the shortest sum of distance directly connected points on adjacent tubes
        minDist = None
        aCrossIndexes = [None, None, None]
        for i in range(self._aroundCounts[0]):
            aCrossIndexes[0] = i
            for j in range(self._aroundCounts[1]):
                aCrossIndexes[1] = j
                for k in range(self._aroundCounts[2]):
                    aCrossIndexes[2] = k
                    dist = 0.0
                    for s in range(3):
                        ring1 = self._endCoordinateRings[s]
                        ring2 = self._endCoordinateRings[s - 2]
                        ic1 = aCrossIndexes[s]
                        ic2 = aCrossIndexes[s - 2]
                        for n in range(1, self._connectionCounts[s]):
                            i1 = (ic1 + n) % self._aroundCounts[s] if self._segmentsIn[s] else ic1 - n
                            i2 = (ic2 - n) if self._segmentsIn[s - 2] else (ic2 + n) % self._aroundCounts[s - 2]
                            delta = sub(ring1[0][i1], ring2[0][i2])
                            dist += magnitude(delta)
                    if (minDist is None) or (dist < minDist):
                        self._aCrossIndexes = copy.copy(aCrossIndexes)
                        minDist = dist
        self._bCrossIndexes = [
            (self._aCrossIndexes[s] + (self._connectionCounts[s] if self._segmentsIn[s] else -self._connectionCounts[s])) % self._aroundCounts[s] for s in range(3)]
        # print("aCrossIndexes", self._aCrossIndexes, "bCrossIndexes", self._bCrossIndexes)

    def _interpolateMidPoint(self, s1, i1, s2, i2):
        """
        Calculate position and derivatives of mid-point between two tubes.
        Algorithm goes halfway from connecting node to end of sampled track surface,
        calculates in-plane direction to the same point on other tube, then blends that direction
        with the outward d2 to get a curve on which to sample the coordinates and direction of the mid-point.
        The mid-point derivative is smoothed from the connecting node coordinates and derivatives.
        :param s1: First tube index from 0 to 2.
        :param i1: Node index around first tube.
        :param s2: Second tube index from 0 to 2.
        :param i2: Node index around second tube.
        :return: Mid-point x, d2.
        """
        trackSurface1 = self._tubeData[s1].getSampledTrackSurface()
        position1 = TrackSurfacePosition(
            i1, (trackSurface1.getElementsCount2() - 1) if self._segmentsIn[s1] else 0,0.0, 0.5)
        p1x, p1d1, p1d2 = trackSurface1.evaluateCoordinates(position1, derivatives=True)
        trackSurface2 = self._tubeData[s2].getSampledTrackSurface()
        position2 = TrackSurfacePosition(
            i2, (trackSurface2.getElementsCount2() - 1) if self._segmentsIn[s2] else 0,0.0, 0.5)
        p2x, p2d1, p2d2 = trackSurface2.evaluateCoordinates(position2, derivatives=True)
        sideFactor = 1.0
        outFactor = 0.5
        direction = normalize(sub(p2x, p1x))
        p1dxi1, p1dxi2 = mult(normalize(calculate_surface_delta_xi(p1d1, p1d2, direction)), sideFactor)
        p1dxi2 += outFactor if self._segmentsIn[s1] else -outFactor
        s1d2 = add(mult(p1d1, p1dxi1), mult(p1d2, p1dxi2))
        p2dxi1, p2dxi2 = mult(normalize(calculate_surface_delta_xi(p2d1, p2d2, direction)), sideFactor)
        p2dxi2 += -outFactor if self._segmentsIn[s2] else outFactor
        s2d2 = add(mult(p2d1, p2dxi1), mult(p2d2, p2dxi2))
        scaling = computeCubicHermiteDerivativeScaling(p1x, s1d2, p2x, s2d2)
        f1d2 = mult(s1d2, scaling)
        f2d2 = mult(s2d2, scaling)
        # print("mag", magnitude(f1d2), "-", magnitude(f2d2))
        mx = interpolateCubicHermite(p1x, f1d2, p2x, f2d2, 0.5)
        md2 = interpolateCubicHermiteDerivative(p1x, f1d2, p2x, f2d2, 0.5)
        return mx, md2

    def _interpolateCrossPoint(self, crossIndexes, swap23):
        """
        Calculate position and derivatives of mid-point between three tubes.
        Algorithm gets the mean of the mid-points between each pair of tubes,
        computes and averages mid-point derivatives from Hermite-Lagrange interpolation,
        the re-smooths to fit the tube connecting coordinates and derivatives.
        :param crossIndexes: Index of nodes for each ring which connect to cross point.
        :param swap23: True if on the first cross point, False if on the other to swap the second and third tubes.
        :return: Cross point x, d1, d2.
        """
        s1i = [0, 2, 1] if swap23 else [0, 1, 2]
        s2i = [2, 1, 0] if swap23 else [1, 2, 0]
        px = []
        pd = []
        hx = []
        hd2 = []
        for s in range(3):
            s1 = s1i[s]
            px.append(self._connectingCoordinateRings[s1][0][crossIndexes[s1]])
            d2 = self._connectingCoordinateRings[s1][2][crossIndexes[s1]]
            pd.append([-d for d in d2] if self._segmentsIn[s1] else d2)
            s2 = s2i[s]
            x, d2 = self._interpolateMidPoint(s1, crossIndexes[s1], s2, crossIndexes[s2])
            hx.append(x)
            hd2.append(d2)
        cx = [(hx[0][c] + hx[1][c] + hx[2][c]) / 3.0 for c in range(3)]
        hd = [interpolateLagrangeHermiteDerivative(cx, px[s], pd[s], 0.0) for s in range(3)]
        ns = [cross(hd[s1i[s]], hd[s2i[s]]) for s in range(3)]
        # ns = [cross(hd2[s1i[s]], hd2[s2i[s]]) for s in range(3)]
        normal = normalize([(ns[0][c] + ns[1][c] + ns[2][c]) for c in range(3)])
        sd = [smoothCubicHermiteDerivativesLine([cx, px[s]], [normalize(cross(normal, cross(hd[s], normal))), pd[s]],
                                                fixStartDirection=True, fixEndDerivative=True)[0] for s in range(3)]
        cd1 = mult(sub(sd[1], add(sd[2], sd[0])), 0.5)
        cd2 = mult(sub(sd[2], add(sd[0], sd[1])), 0.5)
        return cx, cd1, cd2

    def determineMidCoordinates(self):
        """
        Get 3 half-rings of coordinates between the tube coordinate rings.
        The order is between the 1-2, 2-3 and 3-1 tubes.
        Each half ring goes from the first cross point to the second.
        Coordinates at the cross points are averaged.
        """
        assert self._aCrossIndexes is not None  # must call copyCrossIndexes() or determineCrossIndexes() first
        # get cross points a and b
        acx, acd1, acd2 = self._interpolateCrossPoint(self._aCrossIndexes, swap23=False)
        bcx, bcd1, bcd2 = self._interpolateCrossPoint(self._bCrossIndexes, swap23=True)

        self._midCoordinates = []
        for s in range(3):
            self._midCoordinates.append([[acx], [acd1], [acd2]])
            s1 = s
            s2 = (s + 1) % 3
            ring1 = self._connectingCoordinateRings[s1]
            ring2 = self._connectingCoordinateRings[s2]
            for n in range(1, self._connectionCounts[s]):
                i1 = (self._aCrossIndexes[s1] + n) % self._aroundCounts[s1] if self._segmentsIn[s1] \
                    else self._aCrossIndexes[s] - n
                i2 = self._aCrossIndexes[s2] - n if self._segmentsIn[s2] \
                    else (self._aCrossIndexes[s2] + n) % self._aroundCounts[s2]
                hx, hd2 = self._interpolateMidPoint(s1, i1, s2, i2)
                # # print("s", s, "n", n, "/", self._connectionCounts[s], "i1", i1, "i2", i2)
                t1x = ring1[0][i1]
                r1d2 = ring1[2][i1] if self._segmentsIn[s] else [-d for d in ring1[2][i1]]
                t2x = ring2[0][i2]
                r2d2 = [-d for d in ring2[2][i2]] if self._segmentsIn[s - 2] else ring2[2][i2]
                hd1 = [0.0, 0.0, 0.0]
                hd2 = smoothCubicHermiteDerivativesLine([t1x, hx, t2x], [r1d2, hd2, r2d2],
                                                        fixStartDerivative=True, fixEndDerivative=True)[1]
                self._midCoordinates[s][0].append(hx)
                self._midCoordinates[s][1].append(hd1)
                self._midCoordinates[s][2].append(hd2)
            self._midCoordinates[s][0].append(bcx)
            self._midCoordinates[s][1].append(bcd1)
            self._midCoordinates[s][2].append(bcd2)

            # smooth around loops through hex points to get d1
            loopx = self._midCoordinates[s][0]
            startd1 = [-d for d in acd2] if (s == 0) else \
                add(acd1, acd2) if (s == 1) else \
                [-d for d in acd1]
            endd1 = acd1 if (s == 0) else \
                [-d for d in add(acd1, acd2)] if (s == 1) else \
                acd2
            loopd1 = [startd1] + self._midCoordinates[s][1][1:-1] + [endd1]
            loopd1 = smoothCubicHermiteDerivativesLine(loopx, loopd1, fixStartDerivative=True, fixEndDerivative=True,
                                                       magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
            self._midCoordinates[s][1][1:-1] = loopd1[1:-1]


def generateTubeBifurcationTree(networkMesh: NetworkMesh, region, coordinates, nodeIdentifier, elementIdentifier,
                                defaultElementsCountAround: int, targetElementDensityAlongLongestSegment: float,
                                elementsCountThroughWall: int, layoutAnnotationGroups: list=[],
                                annotationElementsCountsAround: list=[],
                                serendipity=False, showTrimSurfaces=False):
    """
    Generate a 2D, or 3D (thick walled) tube bifurcation tree mesh.
    :param networkMesh: Specification of network path and lateral sizes.
    :param region: Zinc region to generate mesh in.
    :param coordinates: Finite element coordinate field to define.
    :param nodeIdentifier: First node identifier to use.
    :param elementIdentifier: First 2D element identifier to use.
    :param defaultElementsCountAround: Number of elements around tube, optionally modified per segment.
    :param targetElementDensityAlongLongestSegment: Target number of elements along the longest segment, used to
    calculate target element length along all segments of network.
    :param elementsCountThroughWall: Number of elements through wall if inner coordinates provided to make 3D elements,
    otherwise 1.
    :param layoutAnnotationGroups: Optional list of annotations defined on layout mesh for networkMesh.
    :param annotationElementsCountsAround: Optional list of elements around annotation groups in the same order as
    layoutAnnotationGroups. A value 0 means ignore; a shorter list means assume 0 for all remaining annotation groups.
    :param serendipity: True to use Hermite serendipity basis, False for regular Hermite with zero cross derivatives.
    :param showTrimSurfaces: Set to True to make surfaces from inter-tube trim surfaces. For diagnostic use.
    :return: next node identifier, next element identifier, annotationGroups.
    """
    layoutRegion = networkMesh.getRegion()
    layoutFieldmodule = layoutRegion.getFieldmodule()
    layoutNodes = layoutFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    layoutCoordinates = layoutFieldmodule.findFieldByName("coordinates").castFiniteElement()
    layoutInnerCoordinates = layoutFieldmodule.findFieldByName("inner coordinates").castFiniteElement()
    if not layoutInnerCoordinates.isValid():
        layoutInnerCoordinates = None
    dimension = 3 if layoutInnerCoordinates else 2
    layoutMesh = layoutFieldmodule.findMeshByDimension(1)
    assert (elementsCountThroughWall == 1) or (layoutInnerCoordinates and (elementsCountThroughWall >= 1))

    fieldmodule = region.getFieldmodule()
    mesh = fieldmodule.findMeshByDimension(dimension)
    fieldcache = fieldmodule.createFieldcache()

    # make tube mesh annotation groups from 1D network layout annotation groups
    annotationGroups = []
    layoutAnnotationMeshGroupMap = []  # List of tuples of layout annotation mesh group to final mesh group
    for layoutAnnotationGroup in layoutAnnotationGroups:
        if layoutAnnotationGroup.getDimension() == 1:
            annotationGroup = AnnotationGroup(region, layoutAnnotationGroup.getTerm())
            annotationGroups.append(annotationGroup)
            layoutAnnotationMeshGroupMap.append(
                (layoutAnnotationGroup.getMeshGroup(layoutMesh), annotationGroup.getMeshGroup(mesh)))

    valueLabels = [
        Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
        Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
        Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

    networkSegments = networkMesh.getNetworkSegments()

    # map from NetworkSegment to SegmentTubeData
    outerSegmentTubeData = {}
    innerSegmentTubeData = {} if layoutInnerCoordinates else None
    longestSegmentLength = 0.0
    for networkSegment in networkSegments:
        pathParameters = get_nodeset_path_ordered_field_parameters(
            layoutNodes, layoutCoordinates, valueLabels,
            networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())
        elementsCountAround = defaultElementsCountAround
        i = 0
        for layoutAnnotationGroup in layoutAnnotationGroups:
            if i >= len(annotationElementsCountsAround):
                break
            annotationElementsCountAround = annotationElementsCountsAround[i]
            if annotationElementsCountAround > 0:
                if networkSegment.hasLayoutElementsInMeshGroup(layoutAnnotationGroup.getMeshGroup(layoutMesh)):
                    elementsCountAround = annotationElementsCountAround
                    # print("Segment in group", i, "using", elementsCountAround, "elements around")
                    break
            i += 1
        outerSegmentTubeData[networkSegment] = tubeData = SegmentTubeData(pathParameters, elementsCountAround)
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(pathParameters, elementsCountAround)
        tubeData.setRawTubeCoordinates((px, pd1, pd2, pd12))
        segmentLength = tubeData.getSegmentLength()
        if segmentLength > longestSegmentLength:
            longestSegmentLength = segmentLength
        if layoutInnerCoordinates:
            innerPathParameters = get_nodeset_path_ordered_field_parameters(
                layoutNodes, layoutInnerCoordinates, valueLabels,
                networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())
            innerSegmentTubeData[networkSegment] = innerTubeData = SegmentTubeData(
                innerPathParameters, elementsCountAround)
            px, pd1, pd2, pd12 = getPathRawTubeCoordinates(innerPathParameters, elementsCountAround)
            innerTubeData.setRawTubeCoordinates((px, pd1, pd2, pd12))
        for layoutAnnotationMeshGroup, annotationMeshGroup in layoutAnnotationMeshGroupMap:
            if networkSegment.hasLayoutElementsInMeshGroup(layoutAnnotationMeshGroup):
                tubeData.addAnnotationMeshGroup(annotationMeshGroup)
                if layoutInnerCoordinates:
                    innerTubeData.addAnnotationMeshGroup(annotationMeshGroup)
    if longestSegmentLength > 0.0:
        targetElementLength = longestSegmentLength / targetElementDensityAlongLongestSegment
    else:
        targetElementLength = 1.0

    # map from NetworkNodes to bifurcation data, resample tube coordinates to fit bifurcation
    outerNodeTubeBifurcationData = {}
    innerNodeTubeBifurcationData = {} if layoutInnerCoordinates else None
    allSegmentTubeData = [outerSegmentTubeData]
    if layoutInnerCoordinates:
        allSegmentTubeData.append(innerSegmentTubeData)
    for segmentTubeData in allSegmentTubeData:
        nodeTubeBifurcationData = innerNodeTubeBifurcationData if (segmentTubeData is innerSegmentTubeData) else \
            outerNodeTubeBifurcationData
        for networkSegment in networkSegments:
            tubeData = segmentTubeData[networkSegment]
            rawTubeCoordinates = tubeData.getRawTubeCoordinates()
            segmentNodes = networkSegment.getNetworkNodes()
            startSegmentNode = segmentNodes[0]
            startTubeBifurcationData = nodeTubeBifurcationData.get(startSegmentNode)
            startSurface = None
            newBifurcationData = []
            if not startTubeBifurcationData:
                startInSegments = startSegmentNode.getInSegments()
                startOutSegments = startSegmentNode.getOutSegments()
                if (len(startInSegments) + len(startOutSegments)) == 3:
                    # print("create start", networkSegment, startSegmentNode)
                    startTubeBifurcationData = TubeBifurcationData(startInSegments, startOutSegments, segmentTubeData)
                    nodeTubeBifurcationData[startSegmentNode] = startTubeBifurcationData
                    newBifurcationData.append(startTubeBifurcationData)
            if startTubeBifurcationData:
                startSurface = startTubeBifurcationData.getSegmentTrimSurface(networkSegment)
            endSegmentNode = segmentNodes[-1]
            endTubeBifurcationData = nodeTubeBifurcationData.get(endSegmentNode)
            endSurface = None
            if not endTubeBifurcationData:
                endInSegments = endSegmentNode.getInSegments()
                endOutSegments = endSegmentNode.getOutSegments()
                if (len(endInSegments) + len(endOutSegments)) == 3:
                    # print("create end", networkSegment, endSegmentNode)
                    endTubeBifurcationData = TubeBifurcationData(endInSegments, endOutSegments, segmentTubeData)
                    nodeTubeBifurcationData[endSegmentNode] = endTubeBifurcationData
                    newBifurcationData.append(endTubeBifurcationData)
            if endTubeBifurcationData:
                endSurface = endTubeBifurcationData.getSegmentTrimSurface(networkSegment)
            if segmentTubeData is outerSegmentTubeData:
                segmentLength = tubeData.getSegmentLength()
                # Previous code setting number of elements along to satisfy targetElementAspectRatio
                # ringCount = len(rawTubeCoordinates[0])
                # sumRingLength = 0.0
                # for n in range(ringCount):
                #     ringLength = getCubicHermiteCurvesLength(rawTubeCoordinates[0][n], rawTubeCoordinates[1][n], loop=True)
                #     sumRingLength += ringLength
                # meanElementLengthAround = sumRingLength / (ringCount * tubeData.getElementsCountAround())
                # targetElementLength = targetElementAspectRatio * meanElementLengthAround
                elementsCountAlong = max(1, math.ceil(segmentLength / targetElementLength))
                loop = (len(startSegmentNode.getInSegments()) == 1) and \
                       (startSegmentNode.getInSegments()[0] is networkSegment) and \
                       (networkSegment.getNodeVersions()[0] == networkSegment.getNodeVersions()[-1])
                if (elementsCountAlong == 1) and startTubeBifurcationData and endTubeBifurcationData:
                    # at least 2 segments if bifurcating at both ends
                    elementsCountAlong = 2
                elif (elementsCountAlong < 2) and loop:
                    # at least 2 segments around loop
                    elementsCountAlong = 2
            else:
                # must match count from outer surface!
                outerTubeData = outerSegmentTubeData[networkSegment]
                elementsCountAlong = outerTubeData.getSampledElementsCountAlong()
            # print("Resample startSurface", startSurface is not None, "endSurface", endSurface is not None)
            sx, sd1, sd2, sd12 = resampleTubeCoordinates(
                rawTubeCoordinates, elementsCountAlong, startSurface=startSurface, endSurface=endSurface)
            tubeData.setSampledTubeCoordinates((sx, sd1, sd2, sd12))
    del segmentTubeData

    completedBifurcations = set()  # record so only done once

    with ChangeManager(fieldmodule):
        for networkSegment in networkSegments:
            segmentNodes = networkSegment.getNetworkNodes()
            startSegmentNode = segmentNodes[0]
            startInSegments = startSegmentNode.getInSegments()
            startOutSegments = startSegmentNode.getOutSegments()
            startSkipCount = 1 if ((len(startInSegments) > 1) or (len(startOutSegments) > 1)) else 0
            endSegmentNode = segmentNodes[-1]
            endInSegments = endSegmentNode.getInSegments()
            endOutSegments = endSegmentNode.getOutSegments()
            endSkipCount = 1 if ((len(endInSegments) > 1) or (len(endOutSegments) > 1)) else 0

            for stage in range(3):
                if stage == 1:
                    # tube
                    outerTubeData = outerSegmentTubeData[networkSegment]
                    outerTubeCoordinates = outerTubeData.getSampledTubeCoordinates()
                    innerTubeData = innerSegmentTubeData[networkSegment] if layoutInnerCoordinates else None
                    innerTubeCoordinates = innerTubeData.getSampledTubeCoordinates() if layoutInnerCoordinates else None
                    startNodeIds = outerTubeData.getStartNodeIds(startSkipCount)
                    endNodeIds = outerTubeData.getEndNodeIds(endSkipCount)
                    loop = (len(startInSegments) == 1) and (startInSegments[0] is networkSegment) and \
                           (networkSegment.getNodeVersions()[0] == networkSegment.getNodeVersions()[-1])
                    nodeIdentifier, elementIdentifier, startNodeIds, endNodeIds = generateTube(
                        outerTubeCoordinates, innerTubeCoordinates, elementsCountThroughWall,
                        region, fieldcache, coordinates, nodeIdentifier, elementIdentifier,
                        startSkipCount=startSkipCount, endSkipCount=endSkipCount,
                        startNodeIds=startNodeIds, endNodeIds=endNodeIds,
                        annotationMeshGroups=outerTubeData.getAnnotationMeshGroups(),
                        loop=loop, serendipity=serendipity)
                    outerTubeData.setStartNodeIds(startNodeIds, startSkipCount)
                    outerTubeData.setEndNodeIds(endNodeIds, endSkipCount)

                    if (len(startInSegments) == 1) and (startSkipCount == 0):
                        # copy startNodeIds to end of last segment
                        inTubeData = outerSegmentTubeData[startInSegments[0]]
                        inTubeData.setEndNodeIds(startNodeIds, 0)
                    if (len(endOutSegments) == 1) and (endSkipCount == 0):
                        # copy endNodesIds to start of next segment
                        outTubeData = outerSegmentTubeData[endOutSegments[0]]
                        outTubeData.setStartNodeIds(endNodeIds, 0)
                else:
                    # start, end bifurcation
                    outerTubeBifurcationData = outerNodeTubeBifurcationData.get(
                        startSegmentNode if (stage == 0) else endSegmentNode)
                    if outerTubeBifurcationData and not outerTubeBifurcationData in completedBifurcations:
                        # if showIntersectionCurves:
                        #     lineIdentifier = None
                        #     for s in range(3):
                        #         curve = outerTubeBifurcationData.getIntersectionCurve(s)
                        #         cx, cd1, cProportions, loop = curve
                        #         if cx:
                        #             nodeIdentifier, lineIdentifier = \
                        #                 generateCurveMesh(region, cx, cd1, loop, nodeIdentifier, lineIdentifier)
                        if showTrimSurfaces:
                            faceIdentifier = elementIdentifier if (dimension == 2) else None
                            for s in range(3):
                                trimSurface = outerTubeBifurcationData.getTrimSurface(s)
                                if trimSurface:
                                    nodeIdentifier, faceIdentifier = \
                                        trimSurface.generateMesh(region, nodeIdentifier, faceIdentifier)
                            if dimension == 2:
                                elementIdentifier = faceIdentifier
                        innerTubeBifurcationData = None
                        if innerNodeTubeBifurcationData:
                            innerTubeBifurcationData = innerNodeTubeBifurcationData.get(
                                startSegmentNode if (stage == 0) else endSegmentNode)

                        crossIndexes = outerTubeBifurcationData.getCrossIndexes()  # only get these from outer
                        if not crossIndexes:
                            outerTubeBifurcationData.determineCrossIndexes()
                            outerTubeBifurcationData.determineMidCoordinates()
                            if innerTubeBifurcationData:
                                innerTubeBifurcationData.copyCrossIndexes(outerTubeBifurcationData)
                                innerTubeBifurcationData.determineMidCoordinates()
                            crossIndexes = outerTubeBifurcationData.getCrossIndexes()

                        outerTubeCoordinates = outerTubeBifurcationData.getConnectingTubeCoordinates()
                        outerMidCoordinates = outerTubeBifurcationData.getMidCoordinates()
                        inward = outerTubeBifurcationData.getSegmentsIn()
                        outerTubeData = outerTubeBifurcationData.getTubeData()
                        tubeNodeIds = [outerTubeData[s].getEndNodeIds(1) if inward[s] else \
                                           outerTubeData[s].getStartNodeIds(1) for s in range(3)]
                        innerTubeCoordinates = None
                        innerMidCoordinates = None
                        if innerTubeBifurcationData:
                            innerTubeCoordinates = innerTubeBifurcationData.getConnectingTubeCoordinates()
                            innerMidCoordinates = innerTubeBifurcationData.getMidCoordinates()
                        annotationMeshGroups = [outerTubeData[s].getAnnotationMeshGroups() for s in range(3)]

                        nodeIdentifier, elementIdentifier = generateTubeBifurcation(
                            outerTubeCoordinates, innerTubeCoordinates, inward, elementsCountThroughWall,
                            outerMidCoordinates, innerMidCoordinates, crossIndexes,
                            region, fieldcache, coordinates, nodeIdentifier, elementIdentifier, tubeNodeIds,
                            annotationMeshGroups, serendipity=serendipity)

                        for s in range(3):
                            if inward[s]:
                                if not outerTubeData[s].getEndNodeIds(1):
                                    outerTubeData[s].setEndNodeIds(tubeNodeIds[s], 1)
                            else:
                                if not outerTubeData[s].getStartNodeIds(1):
                                    outerTubeData[s].setStartNodeIds(tubeNodeIds[s], 1)

                        completedBifurcations.add(outerTubeBifurcationData)

    return nodeIdentifier, elementIdentifier, annotationGroups
