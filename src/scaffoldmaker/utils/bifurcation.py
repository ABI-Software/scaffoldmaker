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
from scaffoldmaker.utils import vector, geometry
from scaffoldmaker.utils.cylindermesh import Ellipse2D
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds, \
    determineTricubicHermiteEft
from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils.interpolation import computeCubicHermiteArcLength, computeCubicHermiteDerivativeScaling, \
    DerivativeScalingMode, getCubicHermiteCurvesLength, getNearestLocationBetweenCurves, getNearestLocationOnCurve, \
    interpolateCubicHermite, interpolateCubicHermiteDerivative, interpolateCubicHermiteSecondDerivative, \
    interpolateLagrangeHermiteDerivative, smoothCubicHermiteDerivativesLine, smoothCubicHermiteDerivativesLoop, \
    smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.networkmesh import NetworkMesh, getPathRawTubeCoordinates, resampleTubeCoordinates
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_delta_xi
from scaffoldmaker.utils.vector import crossproduct3, dotproduct, normalise
from scaffoldmaker.utils.zinc_utils import generateCurveMesh, get_nodeset_path_ordered_field_parameters, \
    get_nodeset_path_field_parameters
import copy
import math
import numpy as np


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
                 region, fieldcache, coordinates: Field, nodeIdentifier, elementIdentifier, segmentIdentifier,
                 innerTubeData, ellipseParameters: list = None, nodeParametersList: list = None,
                 startSkipCount: int = 0, endSkipCount: int = 0, startNodeIds: list = None, endNodeIds: list = None,
                 annotationMeshGroups=[], loop=False, serendipity=False, isCore=False):
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
    :param segmentIdentifier: Identifier used to track tube segments.
    :param innerTubeData: Segment tube data containing inner path parameters.
    :param ellipseParameters: A list of parameters needed to generate an ellipse using Ellipse2D.
    [elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossTransition, elementsCountAlong]
    :param nodeParametersList: A list of node parameters used to determine eft [nodeIdentifier, x, d1, d2, d3]
    :param startSkipCount: Number of element rows to skip along at start.
    :param endSkipCount: Number of element rows to skip along at end.
    :param startNodeIds: Optional existing node identifiers [wall outward][around] to use at start, or None.
    :param endNodeIds: Optional existing node identifiers [wall outward][around] to use at end, or None.
    :param annotationMeshGroups: Mesh groups to add elements to.
    :param loop: Set to true to loop back to start coordinates at end.
    :param serendipity: True to use Hermite serendipity basis, False for regular Hermite with zero cross derivatives.
    :param isCore: True for generating a solid core inside the tube, False for regular tube network
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
    if isCore:
        ellipseParameters[3] = elementsCountAlong
        n = 1 if ellipseParameters[2] > 1 else 0
        coreNodeCount = int(
            (ellipseParameters[0] - n - 2 ** (ellipseParameters[2] - 1)) *
            (ellipseParameters[1] - n - 2 ** (ellipseParameters[2] - 1)) +
            elementsCountAround * ellipseParameters[2])
        coreElementsCount = ((ellipseParameters[0] - 2 * ellipseParameters[2]) *
                             (ellipseParameters[1] - 2 * ellipseParameters[2]))

        corex, cored1, cored2, cored3 = findCoreCoordinates(coreNodeCount, ellipseParameters, innerTubeData,
                                                            ox, ix, od1, segmentIdentifier)

        # Add od3 and id3
        od3 = []
        for n1 in range(len(od2)):
            od3.append([])
            for n2 in range(len(od2[n1])):
                d3 = crossproduct3(od1[n1][n2], od2[n1][n2])
                od3[n1].append(d3)

    assert (not loop) or ((startSkipCount == 0) and (endSkipCount == 0))

    fieldmodule = region.getFieldmodule()

    # Create nodes
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
    if od12 and not serendipity:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

    tubeNodeIds = []
    oFactor = iFactor = 0.0
    if isCore:
        coreNodeIds, coreRimNodeIds = [], []
        nodeParametersList = [] if not nodeParametersList else nodeParametersList
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
        if isCore:
            coreNodeIds.append([])

        for n3 in range(nodesCountThroughWall):
            ringNodeIds = []
            if not isCore:
                otx, otd1, otd2 = (ox[n2], od1[n2], od2[n2])
            else:
                otx, otd1, otd2, otd3 = (ox[n2], od1[n2], od2[n2], od3[n2])
            otd12 = od12[n2] if (od12 and not serendipity) else None
            itx, itd1, itd2 = (ix[n2], id1[n2], id2[n2]) if innerTubeCoordinates else (None, None, None)
            ctx, ctd1, ctd2, ctd3 = (corex[n2], cored1[n2], cored2[n2], cored3[n2]) if isCore else (
                None, None, None, None)
            itd12 = id12[n2] if (innerTubeCoordinates and otd12) else None

            # This needs checking
            # ctd12 = cd12[n2] if (innerTubeCoordinates and otd12 and isCore) else None

            if innerTubeCoordinates:
                oFactor = n3 / elementsCountThroughWall
                iFactor = 1.0 - oFactor
            if n3 == 0 and isCore:
                nodeIds = []
                # Create nodes for the core
                for n1 in range(coreNodeCount):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    fieldcache.setNode(node)

                    rx, rd1, rd2, rd3 = ctx[n1], ctd1[n1], ctd2[n1], ctd3[n1]

                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3)
                    # if itd12:
                    #     rd12 = ctd12[n1]
                    #     coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, rd12)
                    nodeIds.append(nodeIdentifier)
                    nodeParametersList.append([nodeIdentifier, rx, rd1, rd2, rd3])

                    nodeIdentifier += 1

                coreNodeIds, coreRimNodeIds, tubeNodeIds = getCoreTubeNodeIds(coreNodeCount, ellipseParameters,
                        elementsCountAround, coreNodeIds, nodeIds, ringNodeIds, tubeNodeIds, coreRimNodeIds)
            else:
                # Create nodes for non-core layers
                for n1 in range(elementsCountAround):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    fieldcache.setNode(node)
                    if (not innerTubeCoordinates) or (n3 == elementsCountThroughWall):
                        rx, rd1, rd2 = otx[n1], otd1[n1], otd2[n1]
                        rd3 = otd3[n1] if isCore else None
                    elif n3 == 0:
                        rx, rd1, rd2 = itx[n1], itd1[n1], itd2[n1]
                    else:
                        rx = add(mult(otx[n1], oFactor), mult(itx[n1], iFactor))
                        rd1 = add(mult(otd1[n1], oFactor), mult(itd1[n1], iFactor))
                        rd2 = add(mult(otd2[n1], oFactor), mult(itd2[n1], iFactor))
                        rd3 = crossproduct3(rd1, rd2) if isCore else None

                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                    if isCore:
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3)
                    if otd12:
                        rd12 = otd12[n1] if ((not innerTubeCoordinates) or (n3 == elementsCountThroughWall + 1)) else \
                            itd12[n1] if (n3 == 0) else \
                                add(mult(otd12[n1], oFactor), mult(itd12[n1], iFactor))
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, rd12)
                    ringNodeIds.append(nodeIdentifier)
                    if isCore:
                        nodeParametersList.append([nodeIdentifier, rx, rd1, rd2, rd3])

                    nodeIdentifier += 1

                tubeNodeIds[-1].append(ringNodeIds)

    # Create elements
    # element template for bicubic hermite
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

    # Generate elements
    for e2 in range(startSkipCount, elementsCountAlong - endSkipCount):
        for e3 in range(elementsCountThroughWall + 1 + ellipseParameters[2] if isCore else elementsCountThroughWall):
            if e3 == 0 and isCore:
                # Generate elements for the core
                for e1 in range(coreElementsCount):
                    nids, nodeParameters = determineCoreElementNids(e1, e2, e3, ellipseParameters, tubeNodeIds,
                                                                    nodeParametersList)
                    elementIdentifier = generateCoreElements(mesh, coordinates, nids, nodeParameters, elementIdentifier)
            elif e3 == 1 and isCore:
                # Generate elements for the core rim
                for e1 in range(elementsCountAround):
                    nids, nodeParameters = determineCoreRimElementsNids(e1, e2, e3, ellipseParameters,
                            elementsCountAround, tubeNodeIds, nodeParametersList)

                    elementIdentifier = generateCoreRimElements(e1, ellipseParameters, elementsCountAround,
                            elementIdentifier, mesh, coordinates, nids, nodeParameters)
            elif e3 > 1 and isCore:
                # Generate elements for the outer shell
                for e1 in range(elementsCountAround):
                    e2p = e2 + 1
                    e1p = (e1 + 1) % elementsCountAround
                    nids = []
                    n3range = [e3, e3 + 1]

                    for n3 in n3range:
                        nids += [tubeNodeIds[e2][n3][e1], tubeNodeIds[e2][n3][e1p],
                                 tubeNodeIds[e2p][n3][e1], tubeNodeIds[e2p][n3][e1p]]

                    nodeParameters = []
                    for id in nids:
                        for n in range(len(nodeParametersList)):
                            if id == nodeParametersList[n][0]:
                                nodeParameters.append(nodeParametersList[n][1:])

                    # Determine eft and scalefactors
                    eftShell, scalefactorsShell = determineTricubicHermiteEft(mesh, nodeParameters, serendipity=True)

                    elementtemplateShell = mesh.createElementtemplate()
                    elementtemplateShell.setElementShapeType(Element.SHAPE_TYPE_CUBE)

                    elementtemplateShell.defineField(coordinates, -1, eftShell)
                    element = mesh.createElement(elementIdentifier, elementtemplateShell)
                    element.setNodesByIdentifier(eftShell, nids)
                    element.setScaleFactors(eftShell, scalefactorsShell) if scalefactorsShell else "-"

                    elementIdentifier += 1
            else:
                # Create elements for non-core layers
                for e1 in range(elementsCountAround):
                    e2p = e2 + 1
                    e1p = (e1 + 1) % elementsCountAround
                    nids = []
                    n3range = [e3, e3 + 1]

                    for n3 in n3range if (dimension == 3) else [0]:
                        nids += [tubeNodeIds[e2][n3][e1], tubeNodeIds[e2][n3][e1p],
                                 tubeNodeIds[e2p][n3][e1], tubeNodeIds[e2p][n3][e1p]]

                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    element.setNodesByIdentifier(eft, nids)

                    for annotationMeshGroup in annotationMeshGroups:
                        annotationMeshGroup.addElement(element)
                    elementIdentifier += 1

    return nodeIdentifier, elementIdentifier, tubeNodeIds[startSkipCount], \
        tubeNodeIds[elementsCountAlong - endSkipCount], nodeParametersList


def findCoreCentre(ox, ix):
    """
    Finds cartesian coordinates for the centre of the solid core.
    A number of direction vectors (straight lines) are projected based on the outer and inner coordinates of the tube.
    Then the function finds the intersection point of these direction vectors using least squares method.
    :param ox: a list of outer tube coordinates.
    :param ix: a list of inner tube coordinates.
    :return: coordinates at the centre of the core.
    """
    P0 = []
    P1 = []
    elementsCountAround = len(ox[0]) if isinstance(ox[0], list) else len(ox)

    for n3 in range(elementsCountAround):
        P0.append(ox[n3])
        P1.append(ix[n3])
    P0 = np.array(P0)
    P1 = np.array(P1)

    # generate all line direction vectors
    n = (P1 - P0) / np.linalg.norm(P1 - P0, axis=1)[:, np.newaxis]

    # generate the array of all projectors
    projs = np.eye(n.shape[1]) - n[:, :, np.newaxis] * n[:, np.newaxis]

    # generate R matrix and q vector
    R = projs.sum(axis=0)
    q = (projs @ P0[:, :, np.newaxis]).sum(axis=0)

    # solve the least squares problem for the intersection point p: Rp = q
    p = np.linalg.lstsq(R, q, rcond=None)[0]

    Px = p[0][0]
    Py = p[1][0]
    Pz = p[2][0]

    return Px, Py, Pz


def findCoreCoordinates(coreNodeCount, ellipseParameters, pathParameters, ox, ix, od1, segmentIdentifier):
    """
    Blackbox function for generating core coordinates and derivatives using Ellipse2D.
    Steps are:
    1. Find direction vector of the 1D network, centre points of the tube, and major/minor axes for the 2D ellipse.
    2. Generate 2D ellipses.
    3. Extract coordinates and derivatives from 2D ellipses generated in the previous step.
    4. Adjust (rotate and/or scale) derivatives to get the correct direction and magnitude.
    :param coreNodeCount: number of nodes that form the solid core, including the rim layer, e.g. for 8 elements around
    and 4 elements across major axis, there are 17 core nodes (9 nodes forming the square mesh, and 8 rim nodes).
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param pathParameters: parameters describing a path formed by 1D network layout.
    :param ox, ix, od1: outer tube coordinates, inner tube coordinates, and outer tube D1 derivatives, respectively.
    :param segmentIdentifier: An identifier tracking tube segments.
    :return: Core coordinates and derivatives - corex, cored1, cored2, cored3
    """

    # Find directionVector of a tube
    directionVector = findDirectionVector(pathParameters)

    # Find centre points
    centres = []
    for n2 in range(ellipseParameters[3] + 1):
        centre = [0.0, 0.0, 0.0]
        centre[0], centre[1], centre[2] = findCoreCentre(ox[n2], ix[n2])
        centres.append(centre)

    # Find major and minor axes for the ellipse
    majorAxes = []
    minorAxes = []

    for n2 in range(ellipseParameters[3] + 1):
        majorAxes.append([])
        minorAxes.append([])
        for lft in [0, 1]:
            majorAxis, minorAxis = findEllipseAxes(ix[n2], centres[n2], lft)
            majorAxes[-1].append(majorAxis)
            minorAxes[-1].append(minorAxis)

    # Generate ellipses
    lftEllipses = []
    rhtEllipses = []

    for n2 in range(ellipseParameters[3] + 1):
        # Find coreRadii
        coreMajorRadius, coreMinorRadius = findEllipseRadii(ix[n2], centres[n2])

        # Generate ellipse using Ellipse2D
        for lft in [0, 1]:
            ellipse = Ellipse2D(centres[n2], majorAxes[n2][lft], minorAxes[n2][lft],
                                ellipseParameters[0], ellipseParameters[1], 0,
                                ellipseParameters[2], 1.0,
                                coreMajorRadius, coreMinorRadius, isCore=True)

            lftEllipses.append(ellipse) if lft == 0 else rhtEllipses.append(ellipse)

    # Get core coordinates and derivatives from ellipses
    corex, cored1, cored2, cored3 = ([] for i in range(4))

    for n2 in range(ellipseParameters[3] + 1):
        cx, cd1, cd2, cd3 = extractCoordinatesFromEllipse(ellipseParameters, lftEllipses[n2], rhtEllipses[n2])
        [x.append(y) for x, y in zip([corex, cored1, cored2, cored3], [cx, cd1, cd2, cd3])]

    # Adjust D1 and D2 derivatives to correct angles and magnitudes - Needs fixing for trifurcation
    for n2 in range(ellipseParameters[3] + 1):
        cored1[n2] = adjustD1Derivatives(cored1[n2], ellipseParameters, coreNodeCount, od1[n2], segmentIdentifier)
        cored2[n2] = adjustD2Derivatives(cored2[n2], ellipseParameters[3], directionVector)

    return corex, cored1, cored2, cored3


def findDirectionVector(innerTubeData):
    """
    Calculates a direction vector for a particular network segment.
    :param innerTubeData: Segment tube data containing inner path parameters.
    :return: direction vector describing the direction of a tube segment.
    """
    pathParameters = innerTubeData.getPathParameters()
    cpx = pathParameters[0]
    directionVector = (vector.addVectors([cpx[1], cpx[0]], [1, -1]))

    return directionVector  # May no longer need alongAxis


def findEllipseRadii(ix, centre, mode=0):
    """
    Calculates major / minor radii for ellipses along the tube.
    :param ix: inner tube coordinates.
    :param centre: centre point of a tube at a specific section along a tube.
    :param mode: 0: for half-ellipse on the left-hand side of a tube; 1: right-hand side; 2: special mode used for
    mid-section of a bifurcation.
    :return: major and minor radii of an ellipse
    """
    # Calculate coreMajorRadius and coreMinorRadius
    elementsCount = len(ix)
    if mode == 0:
        r2 = 0
        r1 = elementsCount // 4
    elif mode == 1:
        r2 = elementsCount // 4 * 2
        r1 = elementsCount // 4 * 3
    elif mode == 2:
        r2 = elementsCount // 2
        r1 = 0

    majorRadius = vector.magnitude(vector.addVectors([ix[r2], centre], [1, -1]))
    minorRadius = vector.magnitude(vector.addVectors([ix[r1], centre], [1, -1]))

    return majorRadius, minorRadius


def findEllipseAxes(ix, centre, mode=0, coreMajorRadius=None, coreMinorRadius=None):
    """
    Calculates major and minor axes of an ellipse.
    :param ix: inner tube coordinates.
    :param centre: centre point of a tube at a specific section along a tube.
    :param mode: for half-ellipse on the left-hand side of a tube; 1: right-hand side; 2: special mode used for
    mid-section of a bifurcation.
    :param coreMajorRadius: major radius of an ellipse that forms the solid core.
    :param coreMinorRadius: minor radius of an ellipse that forms the solid core.
    :return: major and minor axes of an ellipse.
    """
    elementsCount = len(ix)
    assert mode in [0, 1, 2]

    if mode == 0:
        a1 = 0
        a2 = elementsCount // 4
        scaleFactor = 1
    elif mode == 1:
        a1 = elementsCount // 4 * 2
        a2 = elementsCount // 4 * 3
        scaleFactor = -1
    elif mode == 2:
        a1 = elementsCount // 2
        a2 = 0
        scaleFactor = 1

    # Calculate coreMajorRadius and coreMinorRadius if it doesn't exist
    if coreMajorRadius == None and coreMinorRadius == None:
        coreMajorRadius, coreMinorRadius = findEllipseRadii(ix, centre, mode)
    else:
        coreMajorRadius = coreMajorRadius
        coreMinorRadius = coreMinorRadius

    majorAxis = vector.addVectors([ix[a1], centre], [1, -1])
    minorAxis = [0, 0, coreMinorRadius * scaleFactor]

    return majorAxis, minorAxis


def extractCoordinatesFromEllipse(ellipseParameters, lftEllipse, rhtEllipse):
    """
    Extracts coordinates and derivatives from ellipses generated using Ellipse2D.
    Two ellipses are generated (left-hand side and right-hand side) for each slice along the tube.
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param lftEllipse: A 2D ellipse object on the left-hand side of a tube.
    :param rhtEllipse: A 2D ellipse object on the right-hand side of a tube.
    :return: coordinates and derivatives - cx, cd1, cd2, cd3
    """
    # Left-hand side half ellipse
    x, d1, d2, d3 = [], [], [], []

    for n2 in range(ellipseParameters[0] + 1):
        for n1 in range(ellipseParameters[1] + 1):
            if lftEllipse.px[n2][n1]:
                x.append(lftEllipse.px[n2][n1])
                d1.append(lftEllipse.pd1[n2][n1])
                d2.append(lftEllipse.pd2[n2][n1])
                d3.append(lftEllipse.pd3[n2][n1])

    lftx, lftd1, lftd2, lftd3 = x, d1, d2, d3

    # Right-hand side half ellipse
    x, d1, d2, d3 = [], [], [], []
    lst = (list(range(0, ellipseParameters[2])) +
           list(range(ellipseParameters[1] - (ellipseParameters[2] - 1), ellipseParameters[1] + 1)))

    for n2 in range(ellipseParameters[0] + 1):
        for n1 in range(ellipseParameters[1] + 1):
            if rhtEllipse.px[n2][n1]:
                # switch direction of D1 and D3 derivatives to match the D1 and D3 on left-hand side
                if n2 > ellipseParameters[2] - 1 and n1 not in lst:
                    rhtEllipse.pd1[n2][n1] = vector.scaleVector(rhtEllipse.pd1[n2][n1], -1)
                    rhtEllipse.pd3[n2][n1] = vector.scaleVector(rhtEllipse.pd3[n2][n1], -1)

                x.append(rhtEllipse.px[n2][n1])
                d1.append(rhtEllipse.pd1[n2][n1])
                d2.append(rhtEllipse.pd2[n2][n1])
                d3.append(rhtEllipse.pd3[n2][n1])

    # reverse the order of the list and slice only the nodes not already included on the left-hand side
    [lst.reverse() for lst in [x, d1, d2, d3]]

    idx = int(ellipseParameters[1] + 1)
    x, d1, d2, d3 = x[idx:], d1[idx:], d2[idx:], d3[idx:]

    rhtx, rhtd1, rhtd2, rhtd3 = x, d1, d2, d3
    cx, cd1, cd2, cd3 = (lftx + rhtx), (lftd1 + rhtd1), (lftd2 + rhtd2), (lftd3 + rhtd3)

    return cx, cd1, cd2, cd3


def adjustD2Derivatives(cd2, elementsCountAlong, directionVector):
    """
    Rotates and/or scales d2 derivatives of tube core to match d2 derivatives of surrounding inner tube layer.
    :param cd2: D2 derivatives of tube core.
    :param elementsCountAlong: A number of elements along a tube.
    :param directionVector: A vector describing the direction of a 1D network layout that form the tube.
    :return: cd2
    """
    # Adjust d2 derivatives to correct angles and magnitudes
    nodeCount = len(cd2)

    for n2 in range(nodeCount):
        d2 = vector.normalise(cd2[n2])
        angleD2 = vector.angleBetweenVectors(d2, directionVector)
        if angleD2 > 0:
            normD2 = directionVector
            scaleD2 = vector.scaleVector(normD2, 1 / elementsCountAlong)
            cd2[n2] = scaleD2
        else:
            cd2[n2] = vector.scaleVector(cd2[n2], 1 / elementsCountAlong)

    return cd2


def adjustD1Derivatives(cd1, ellipseParameters, coreNodeCount, od1, segmentIdentifier):
    """
    Rotates d1 derivatives of tube core to match d1 derivatives of surrounding inner tube layer.
    :param cd1: D1 derivatives of tube core.
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param coreNodeCount: A number of nodes that form the solid core, including the rim layer.
    :param od1: D1 derivatives of outer tube, used as a reference point.
    :param segmentIdentifier: An identifier tracking tube segments.
    :return: d1 derivatives of tube core - cd1
    """
    start = int(coreNodeCount // 2 - (ellipseParameters[1] / 2 - (ellipseParameters[2] - 1)))
    end = start + (ellipseParameters[1] - (ellipseParameters[2] - 1)) + 1
    pos1 = len(od1) // 4
    pos2 = pos1 * 3

    for n2 in range(start, end):
        d1 = cd1[n2]
        if start < n2 < end:
            if dotproduct(normalise(d1), normalise(od1[pos1])) > 1.0:
                angleD1 = 0
            else:
                angleD1 = vector.angleBetweenVectors(d1, od1[pos1])
            k = [0, 0, 1] if segmentIdentifier > 1 else [0, 0, -1]
            scaleFactor = 1
        elif n2 == end:
            if dotproduct(normalise(d1), normalise(od1[pos1])) > 1.0:
                angleD1 = 0
            else:
                angleD1 = vector.angleBetweenVectors(d1, od1[pos1])
            k = [0, 0, -1] if segmentIdentifier > 1 else [0, 0, 1]
            scaleFactor = -1
        else:
            # if normalise(od1[pos2]) == normalise(d1):
            if dotproduct(normalise(od1[pos2]), normalise(d1)) > 1.0:
                angleD1 = 0
            else:
                angleD1 = vector.angleBetweenVectors(od1[pos2], d1)
            k = [0, 0, -1] if segmentIdentifier > 1 else [0, 0, 1]
            scaleFactor = -1
        if angleD1 > 0:
            rotateD1 = vector.rotateVectorAroundVector(d1, k, angleD1 * scaleFactor)
            cd1[n2] = rotateD1

    return cd1


def findRingNodeIndicesForSorting(coreNodeCount, ellipseParameters, halfEllipse=False):
    """
    Determines indices of nodeIds list, which are used to sort coreNodeIds and ringNodeIds from more general nodeIds list.
    :param coreNodeCount: A number of nodes that form the solid core, including the rim layer.
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param halfEllipse: True if an ellipse is half-shaped, false if an ellipse is full-shaped.
    :return: A list of indices indicating which nodes are core nodes and which are ring nodes.
    """
    indicesList = []
    idx = ellipseParameters[1] - 1 - 2 * (ellipseParameters[2] - 1)

    for k in range(ellipseParameters[2]):
        end1 = coreNodeCount - 1 - (idx * k)
        start1 = end1 - idx + 1 if not halfEllipse else end1 - 1

        for i in range(coreNodeCount):
            if idx * k <= i < idx * (1 + k):
                indicesList.append(i)
            elif start1 <= i <= end1 and not halfEllipse:
                indicesList.append(i)
            else:
                start2 = (ellipseParameters[1] - 1 - 2 * (ellipseParameters[2] - 1)) * (1 + ellipseParameters[2])
                for j in range(ellipseParameters[0] - 3 - 2 * (ellipseParameters[2] - 1)):
                    start2p = start2 + k + (ellipseParameters[1] + 1) * j
                    end2 = start2p + (ellipseParameters[1] - 2 * k)
                    if i in [start2p, end2]:
                        indicesList.append(i)

    return indicesList


def sortCoreNodeIds(indexList, coreNodeCount, coreNodeIds, nodeIds):
    """
    A helper function for sorting between core node ids and rim node ids in a general list of node ids.
    :param indexList: A list of indices indicating which nodes are core nodes and which are ring nodes, determined from
    findRingNodeIndicesForSorting function.
    :param coreNodeCount: A number of nodes that form the solid core, including the rim layer.
    :param coreNodeIds: A list of node ids that form the solid core. Initially empty.
    :nodeIds: A flat list of all node ids excluding the node ids that form the outer layer.
    :return: coreNodeIds contain all ids that collectively form the solid core. tmpNodeIds are appended to ringNodeIds
    at a later stage.
    """
    tmpNodeIds = [None for i in range(len(indexList))]

    for i in range(coreNodeCount):
        if i in indexList:
            tmpNodeIds[indexList.index(i)] = nodeIds[i]
        else:
            coreNodeIds[-1].append(nodeIds[i])

    return tmpNodeIds, coreNodeIds


def rearrangeNodeIds(nodeIds, ellipseParameters):
    """
    A helper function for rearranging node ids in a list to be circular, similar to other tube node ids rotating in a
    clockwise direction.
    :param nodeIds: A list of node ids.
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :return: A list of node ids rearranged in a circular format.
    """
    idx = ellipseParameters[1] - 1 - 2 * (ellipseParameters[2] - 1)
    lst = nodeIds[-idx::]
    lst.reverse()
    nodeIds[-idx::] = lst

    start = ellipseParameters[0] - 4 - 2 * (ellipseParameters[2] - 1)
    for j in range(int(start), -1, -1):
        nodeIds.append(nodeIds.pop(idx + 2 * j))

    nloop = int(ellipseParameters[1] / 2 - ellipseParameters[2])
    for i in range(nloop):
        nodeIds.insert(len(nodeIds), nodeIds.pop(0))

    return nodeIds


def rearrangeMidNodeIds(nodeIds, ellipseParameters):
    """
    A helper function for rearranging node ids in a list to be circular. This works in a similar way to
    rearrangeNodeIds function but only used specifically for mid-section nodes in a bifurcation.
    :param nodeIds: A list of node ids.
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :return: A list of node ids rearranged in a circular format.
    """
    idx = ellipseParameters[1] - 2 * (ellipseParameters[2] - 1)
    nodeIds.reverse()
    end = ellipseParameters[0] // 2 - ellipseParameters[2]
    end = 1 if end <= 0 else end

    for j in range(int(end)):
        nodeIds.append(nodeIds.pop(-idx - 2 * j))

    return nodeIds


def createCoreRimIdsList(coreNodeIds, ellipseParameters, mid=False):
    """
    A helper function for creating a list of node ids that collectively form the rim layer around the core. These are
    the nodes that form the outer boundaris of a square mesh. e.g. for the simplest 2 x 2 solid core, the rim node ids
    will be all core node ids excluding the node id at the centre.
        r - r - r
        r - c - r
        r - r - r
    :param coreNodeIds: A list of core node ids.
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param mid: False if it is a tube, True if it is a mid-section of a bifurcation.
    :return: A list of rim node ids.
    """
    elementsCountAcrossMajor = ellipseParameters[0]
    elementsCountAcrossMinor = ellipseParameters[1]
    rimNodeIds = coreNodeIds[-1].copy()
    boundaryIndices = []
    em = elementsCountAcrossMinor - 1

    if not mid or elementsCountAcrossMajor == 4:
        eM = elementsCountAcrossMajor - 2
        eS = 2 * (ellipseParameters[2] - 1)
    else:
        eM = elementsCountAcrossMajor // 2
        eS = (ellipseParameters[2] - 1)

    for i in range(1, eM - eS):
        r1 = (em - 2 * (ellipseParameters[2] - 1)) * i + 1
        r2 = r1 + (em - 2 - 2 * (ellipseParameters[2] - 1))
        boundaryIndices += list(range(r1, r2))

    for idx in sorted(boundaryIndices, reverse=True):
        if rimNodeIds[idx]:
            del rimNodeIds[idx]

    return rimNodeIds


def getCoreTubeNodeIds(coreNodeCount, ellipseParameters, elementsCountAround, coreNodeIds, nodeIds, ringNodeIds,
                       tubeNodeIds, coreRimNodeIds, s=None, isMid=False):
    """
    Sorts and rearranges node ids and returns three sets of node id lists required for generating elements.
    :param coreNodeCount: Number of nodes that form the solid core, including the rim layer
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param elementsCountAround: Number of elements around the solid core.
    :param coreNodeIds, nodeIds, ringNodeIds, tubeNodeIds, coreRimNodeIds: Lists containing node ids.
    :param s: segment identifier. This is not used for a straight cylinder.
    :param isMid: False if generating a straight cylinder, True if generating bifurcation.
    :return: Three sets of node ids lists - coreNodeIds, coreRimNodeIds, tubeNodeIds
    """
    # Rearrange nodes
    sortedNodeIndices = findRingNodeIndicesForSorting(coreNodeCount, ellipseParameters)
    tmpNodeIds, coreNodeIds = sortCoreNodeIds(sortedNodeIndices, coreNodeCount, coreNodeIds, nodeIds)

    for k in range(ellipseParameters[2]):
        ringNodeIds.append([])
        for i in range(elementsCountAround):
            ringNodeIds[k].insert(i, tmpNodeIds.pop(0))
    del tmpNodeIds

    for k in range(ellipseParameters[2]):
        ringNodeIds[k] = rearrangeNodeIds(ringNodeIds[k], ellipseParameters)

    # Create a new list for core boundary nodes
    tmpBoundaryIds = createCoreRimIdsList(coreNodeIds, ellipseParameters)
    rearrangedBoundaryIds = rearrangeNodeIds(tmpBoundaryIds, ellipseParameters)
    coreRimNodeIds.append(rearrangedBoundaryIds)

    m = -1 if not isMid else s
    tubeNodeIds[m].append(coreNodeIds[-1])
    tubeNodeIds[m].append(coreRimNodeIds[-1])

    for k in range((ellipseParameters[2] - 1), -1, -1):
        tubeNodeIds[m].append(ringNodeIds[k])

    return coreNodeIds, coreRimNodeIds, tubeNodeIds


def determineCoreElementNids(e1, e2, e3, ellipseParameters, tubeNodeIds, nodeParametersList):
    """
    Determines a set of 8 node ids required to form a solid core element, and searches for node parameters matching
    the selected set of node ids from nodeParametersList.
    :param e1, e2, e3: index for elementsCountAround, elementsCountAlong, and elementsCountThrough, respectively.
    :param ellipseParameters: A list of parameters needed to generate an ellipse.
    :param tubeNodeIds: A list of tube node ids.
    :param nodeParametersList: A list of node parameters.
    :return: a set of node ids and their node parameter values
    """
    # Determine 8 node ids that form elements of the core
    nodeParameters = []
    e1p = int(e1 + e1 // (ellipseParameters[1] - 2 * (ellipseParameters[2])))
    skipIndex = ellipseParameters[1] - 2 * ellipseParameters[2] + 1

    n1 = tubeNodeIds[e2][e3][e1p]
    n2 = tubeNodeIds[e2][e3][e1p + skipIndex]
    n3 = tubeNodeIds[e2 + 1][e3][e1p]
    n4 = tubeNodeIds[e2 + 1][e3][e1p + skipIndex]

    n1p = tubeNodeIds[e2][e3][e1p + 1]
    n2p = tubeNodeIds[e2][e3][e1p + skipIndex + 1]
    n3p = tubeNodeIds[e2 + 1][e3][e1p + 1]
    n4p = tubeNodeIds[e2 + 1][e3][e1p + skipIndex + 1]

    nids = [n1, n2, n3, n4, n1p, n2p, n3p, n4p]

    for id in nids:
        for n in range(len(nodeParametersList)):
            if id == nodeParametersList[n][0]:
                nodeParameters.append(nodeParametersList[n][1:])

    return nids, nodeParameters

def generateCoreElements(mesh, coordinates, nids, nodeParameters, elementIdentifier):
    """
    First determines the required eft and scalefactors using determineTricubicHermite function based on nodeParameters,
    and then generates a core element for a given set of node ids.
    :param mesh:  A Zinc mesh of dimension 3.
    :param coordinates: Finite element coordinate field to define.
    :param nids: A set of 8 node ids that forms an element.
    :param nodeParameters: A list of node parameters for given set of nids.
    :param elementIdentifier: First 2D element identifier to use.
    :return: elementIdentifier for the next element to be generated.
    """
    # Determine eft and scalefactors
    eftCore, scalefactorsCore = determineTricubicHermiteEft(mesh, nodeParameters, serendipity=True)

    # Create element template
    elementtemplateCore = mesh.createElementtemplate()
    elementtemplateCore.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplateCore.defineField(coordinates, -1, eftCore)

    # Generate core elements
    element = mesh.createElement(elementIdentifier, elementtemplateCore)
    result1 = element.setNodesByIdentifier(eftCore, nids)
    result2 = element.setScaleFactors(eftCore, scalefactorsCore) if scalefactorsCore else "-"

    elementIdentifier += 1

    return elementIdentifier


def determineCoreRimElementsNids(e1, e2, e3, ellipseParameters, elementsCountAround, tubeNodeIds, nodeParametersList):
    """
    Determines a set of 8 node ids required to form a core rim element, and searches for node parameters matching
    the selected set of node ids from nodeParametersList.
    :param e1, e2, e3: index for elementsCountAround, elementsCountAlong, and elementsCountThrough, respectively.
    :param ellipseParameters: A list of parameters needed to generate an ellipse.
    :param elementsCountAround: Number of elements around solid core.
    :param tubeNodeIds: A list of tube node ids.
    :param nodeParametersList: A list of node parameters.
    :return: a set of node ids and their node parameter values
    """
    em = (ellipseParameters[1] - 2) // 2 - (ellipseParameters[2] - 1)
    eM = (ellipseParameters[0] - 2) // 2 - (ellipseParameters[2] - 1)
    ec = elementsCountAround // 4

    topRowElements = list(range(0, ec - eM)) + list(range(3 * ec + eM, elementsCountAround))
    btmRowElements = list((range(2 * ec - em, 2 * ec + em)))
    lftColumnElements = list(range(3 * ec - eM, 3 * ec + eM))
    rhtColumnElements = list(range(ec - eM, ec + eM))

    nodeParameters = []

    # Find nids for the core rim elements
    n1 = tubeNodeIds[e2][e3 + 1][e1]
    n2 = tubeNodeIds[e2][e3 + 1][(e1 + 1) % elementsCountAround]
    n3 = tubeNodeIds[e2 + 1][e3 + 1][e1]
    n4 = tubeNodeIds[e2 + 1][e3 + 1][(e1 + 1) % elementsCountAround]

    n1p = tubeNodeIds[e2][e3][e1]
    n2p = tubeNodeIds[e2][e3][(e1 + 1) % elementsCountAround]
    n3p = tubeNodeIds[e2 + 1][e3][e1]
    n4p = tubeNodeIds[e2 + 1][e3][(e1 + 1) % elementsCountAround]

    if e1 in topRowElements:
        nids = [n1, n1p, n3, n3p, n2, n2p, n4, n4p]
    elif e1 in btmRowElements:
        nids = [n2p, n2, n4p, n4, n1p, n1, n3p, n3]
    elif e1 in lftColumnElements:
        nids = [n2, n1, n4, n3, n2p, n1p, n4p, n3p]
    elif e1 in rhtColumnElements:
        nids = [n1p, n2p, n3p, n4p, n1, n2, n3, n4]

    for id in nids:
        for n in range(len(nodeParametersList)):
            if id == nodeParametersList[n][0]:
                nodeParameters.append(nodeParametersList[n][1:])

    return nids, nodeParameters


def generateCoreRimElements(e1, ellipseParameters, elementsCountAround, elementIdentifier, mesh, coordinates,
                            nids, nodeParameters):
    """
    First determines the required eft and scalefactors using determineTricubicHermite function based on nodeParameters,
    and then generates a core rim element for a given set of node ids.
    :param e1: index for elementsCountAround
    :param ellipseParameters: A list of parameters needed to generate an ellipse.
    :param elementsCountAround: Number of elements around solid core.
    :param elementIdentifier: First 2D element identifier to use.
    :param mesh:  A Zinc mesh of dimension 3.
    :param coordinates: Finite element coordinate field to define.
    :param nids: A set of 8 node ids that forms an element.
    :param nodeParameters: A list of node parameters for given set of nids.
    :return: elementIdentifier for the next element to be generated.
    """
    em = (ellipseParameters[1] - 2) // 2 - (ellipseParameters[2] - 1)
    eM = (ellipseParameters[0] - 2) // 2 - (ellipseParameters[2] - 1)
    ec = elementsCountAround // 4

    topRowElements = list(range(0, ec - eM)) + list(range(3 * ec + eM, elementsCountAround))
    btmRowElements = list((range(2 * ec - em, 2 * ec + em)))
    lftColumnElements = list(range(3 * ec - eM, 3 * ec + eM))
    rhtColumnElements = list(range(ec - eM, ec + eM))

    # Top row rim elements of the solid core
    if e1 in topRowElements:
        idx = len(topRowElements) // 2
        if e1 == topRowElements[idx]:
            nodeDerivativeFixedWeights = [[], [[1.0, 0.0, 1.0]], [], [[1.0, 0.0, 1.0]],
                                          [], [], [], []]
        elif e1 == topRowElements[idx - 1]:
            nodeDerivativeFixedWeights = [[], [], [], [],
                                          [], [[1.0, 0.0, -1.0]], [], [[1.0, 0.0, -1.0]]]
        else:
            nodeDerivativeFixedWeights = None

        # Determine eft and scalefactors
        eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                      serendipity=True)
    # Bottom row rim elements of the solid core
    elif e1 in btmRowElements:
        if e1 == btmRowElements[0]:
            nodeDerivativeFixedWeights = [[], [], [], [],
                                          [[1.0, 0.0, 1.0]], [], [[1.0, 0.0, 1.0]], []]
        elif e1 == btmRowElements[-1]:
            nodeDerivativeFixedWeights = [[[1.0, 0.0, -1.0]], [], [[1.0, 0.0, -1.0]], [],
                                          [], [], [], []]
        else:
            nodeDerivativeFixedWeights = None

        # Determine eft and scalefactors
        eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                      serendipity=True)
    # Left-hand side column elements of the solid core
    elif e1 in lftColumnElements:
        if e1 == lftColumnElements[0]:
            nodeDerivativeFixedWeights = [[], [], [], [],
                                          [], [None, None, [-1.0, 0.0, 1.0]], [], [None, None, [-1.0, 0.0, 1.0]]]
        elif e1 == lftColumnElements[-1]:
                nodeDerivativeFixedWeights = [[], [], [], [],
                                              [None, None, [1.0, 0.0, 1.0]], [], [None, None, [1.0, 0.0, 1.0]], []]
        else:
            nodeDerivativeFixedWeights = None

        # Determine eft and scalefactors
        eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                      serendipity=True)
    # Right-hand side column elements of the solid core
    elif e1 in rhtColumnElements:
        if e1 == rhtColumnElements[0]:
            nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, 1.0]], [], [None, None, [-1.0, 0.0, 1.0]], [],
                                          [], [], [], []]
        elif e1 == rhtColumnElements[-1]:
            nodeDerivativeFixedWeights = [[], [None, None, [1.0, 0.0, 1.0]], [], [None, None, [1.0, 0.0, 1.0]],
                                          [], [], [], []]
        else:
            nodeDerivativeFixedWeights = None

        # Determine eft and scalefactors
        eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                      serendipity=True)

    elementtemplateCoreRim = mesh.createElementtemplate()
    elementtemplateCoreRim.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    elementtemplateCoreRim.defineField(coordinates, -1, eftCoreRim)
    element = mesh.createElement(elementIdentifier, elementtemplateCoreRim)
    element.setNodesByIdentifier(eftCoreRim, nids)
    element.setScaleFactors(eftCoreRim, scalefactorsCoreRim) if scalefactorsCoreRim else "-"

    elementIdentifier += 1

    return elementIdentifier


def generateTubeBifurcation(outerTubeCoordinates, innerTubeCoordinates, inward, elementsCountThroughWall,
                            outerMidCoordinates, innerMidCoordinates, crossIndexes,
                            region, fieldcache, coordinates: Field,
                            nodeIdentifier, elementIdentifier, tubeNodeIds,
                            segmentIdentifier, outerTubeData,
                            annotationMeshGroups, ellipseParameters: list = None, nodeParametersList : list = None,
                            serendipity=False, isCore=False):
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
    :param segmentIdentifier: Identifier used to track tube segments.
    :param outerTubeData: Segment tube data containing outer tube path parameters.
    :param annotationMeshGroups: List over 3 tubes of lists of meshGroups to add elements to for the part of the
    bifurcation that is part of that segment.
    :param serendipity: True to use Hermite serendipity basis, False for regular Hermite with zero cross derivatives.
    between the first generated tube and the bifurcation. The value is None for regular tubes.
    :param ellipseParameters:A list of parameters needed to generate an ellipse using Ellipse2D.
    [elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossTransition, elementsCountAlong]
    :param nodeParametersList: A list of node parameters used to determine eft [nodeIdentifier, x, d1, d2, d3]
    :param isCore: True for generating a solid core inside the tube, False for regular tube network
    :return: next node identifier, next element identifier
    """
    dimension = 3 if innerTubeCoordinates else 2
    assert ((dimension == 2) and (elementsCountThroughWall == 1)) or \
           ((dimension == 3) and (elementsCountThroughWall >= 1))
    nodesCountThroughWall = (elementsCountThroughWall + 1) if (dimension == 3) else 1

    fieldmodule = region.getFieldmodule()

    segmentCount = len(outerTubeCoordinates)
    aroundCounts = [len(outerTubeCoordinates[s][0]) for s in range(segmentCount)]
    connectionCounts = get_tube_bifurcation_connection_elements_counts(aroundCounts)

    if isCore:
        isMid = True
        coreElementsCount = ((ellipseParameters[0] - 2 * ellipseParameters[2]) *
                             (ellipseParameters[1] - 2 * ellipseParameters[2]))
        coreElementsCountHalf = coreElementsCount // 2
        coreCentre = getCoreMidCentre(outerMidCoordinates, innerMidCoordinates)

        coreNodeCount = []
        n = 1 if ellipseParameters[2] > 1 else 0
        for s in range(segmentCount):
            nodeCount = int((ellipseParameters[0] - n - 2 ** (ellipseParameters[2] - 1)) *
                            (ellipseParameters[1] - n - 2 ** (ellipseParameters[2] - 1)) +
                            aroundCounts[s] * ellipseParameters[2])
            coreNodeCount.append(nodeCount)

        # Find coordinates and derivatives for mid core nodes
        coremx, coremd1, coremd2, coremd3 = findMidCoreCoordinates(coreCentre, innerMidCoordinates, ellipseParameters,
                                                                   segmentCount, inward)

    # Create nodes

    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
    if (not serendipity) and (len(outerTubeCoordinates[0]) > 3):
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
    nodetemplateCross = nodes.createNodetemplate()
    nodetemplateCross.defineField(coordinates)
    nodetemplateCross.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplateCross.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplateCross.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplateCross.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    midNodeIds = []
    oFactor = iFactor = 0.0
    if isCore:
        coreTubeNodeIds, coreTubeRimNodeIds, coreMidNodeIds, coreMidRimNodeIds = [], [], [], []
        midColumnNodeIds = None
    for s in range(segmentCount):
        # ensure tube nodes are supplied or create here
        if not tubeNodeIds[s]:
            tubeNodeIds[s] = []
            otx, otd1, otd2 = outerTubeCoordinates[s][:3]
            otd12 = None if (serendipity or (len(outerTubeCoordinates[s]) == 3)) \
                else outerTubeCoordinates[s][3]
            itx, itd1, itd2 = innerTubeCoordinates[s][:3] if innerTubeCoordinates else (None, None, None)
            itd12 = None if ((not innerTubeCoordinates) or serendipity or (len(innerTubeCoordinates[s]) == 3)) \
                else innerTubeCoordinates[s][3]

            if isCore:
                ctx, ctd1, ctd2, ctd3 = findCoreTubeCoordinates(coreTubeNodeIds, otx, otd1, itx, ellipseParameters,
                        aroundCounts[s], outerTubeData[s], segmentCount, coreNodeCount[s])

            for n3 in range(nodesCountThroughWall):
                ringNodeIds = []
                if innerTubeCoordinates:
                    oFactor = n3 / elementsCountThroughWall
                    iFactor = 1.0 - oFactor

                if n3 == 0 and isCore:
                    nodeIds = []
                    # Create nodes for tube core
                    for n1 in range(coreNodeCount[s]):
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)

                        rx, rd1, rd2, rd3 = ctx[n1], ctd1[n1], ctd2[n1], ctd3[n1]

                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3)

                        nodeIds.append(nodeIdentifier)
                        nodeParametersList.append([nodeIdentifier, rx, rd1, rd2, rd3])

                        nodeIdentifier += 1

                    coreTubeNodeIds, coreTubeRimNodeIds, tubeNodeIds = getCoreTubeNodeIds(coreNodeCount[s],
                            ellipseParameters, aroundCounts[s], coreTubeNodeIds, nodeIds, ringNodeIds, tubeNodeIds,
                            coreTubeRimNodeIds, s, isMid)
                else:
                    # Create nodes for non-core layers
                    for n1 in range(aroundCounts[s]):
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        if (not innerTubeCoordinates) or (n3 == elementsCountThroughWall):
                            rx, rd1, rd2 = otx[n1], otd1[n1], otd2[n1]
                            rd3 = crossproduct3(rd1, rd2) if isCore else None
                        elif n3 == 0:
                            rx, rd1, rd2 = itx[n1], itd1[n1], itd2[n1]
                        else:
                            rx = add(mult(otx[n1], oFactor), mult(itx[n1], iFactor))
                            rd1 = add(mult(otd1[n1], oFactor), mult(itd1[n1], iFactor))
                            rd2 = add(mult(otd2[n1], oFactor), mult(itd2[n1], iFactor))
                            rd3 = crossproduct3(rd1, rd2) if isCore else None

                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                        if isCore:
                            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3)
                        if otd12:
                            rd12 = otd12[n1] if ((not innerTubeCoordinates) or (n3 == elementsCountThroughWall)) else \
                                itd12[n1] if (n3 == 0) else add(mult(otd12[n1], oFactor), mult(itd12[n1], iFactor))
                            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, rd12)

                        ringNodeIds.append(nodeIdentifier)
                        if isCore:
                            nodeParametersList.append([nodeIdentifier, rx, rd1, rd2, rd3])

                        nodeIdentifier += 1

                    tubeNodeIds[s].append(ringNodeIds)

        # Create nodes around middle half rings
        midNodeIds.append([])
        omx, omd1, omd2 = outerMidCoordinates[s][:3]
        omd12 = None if (serendipity or (len(outerMidCoordinates[s]) == 3)) else outerMidCoordinates[s][3]
        imx, imd1, imd2 = innerMidCoordinates[s][:3] if innerMidCoordinates else (None, None, None)
        imd12 = None if ((not innerMidCoordinates) or serendipity or (len(innerMidCoordinates[s]) == 3)) \
            else innerMidCoordinates[s][3]
        elementsCountAroundHalf = len(omx) - 1
        if isCore:
            coreMidNodeIds.append([])
            cmx, cmd1, cmd2, cmd3 = coremx[s], coremd1[s], coremd2[s], coremd3[s]
            # Adjust d1 & d2 derivatives
            cmd1, cmd2 = adjustMidCoreDerivatives(ellipseParameters[1], cmd1, cmd2, imd1, imd2, s, inward)

        for n3 in range(nodesCountThroughWall):
            ringNodeIds = []
            if innerTubeCoordinates:
                oFactor = n3 / elementsCountThroughWall
                iFactor = 1.0 - oFactor
            if n3 == 0 and isCore:
                # Create mid nodes for the core
                nodeIds = []
                coreMidNodeCount = len(cmx)
                for n1 in range(coreMidNodeCount):
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    fieldcache.setNode(node)

                    rx, rd1, rd2, rd3 = cmx[n1], cmd1[n1], cmd2[n1], cmd3[n1]

                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3)

                    nodeIds.append(nodeIdentifier)
                    nodeParametersList.append([nodeIdentifier, rx, rd1, rd2, rd3])

                    nodeIdentifier = nodeIdentifier + 1

                # Rearrange nodes for generating elements
                coreMidNodeIds, coreMidRimNodeIds, midColumnNodeIds, midNodeIds = getCoreMidNodeIds(s, ellipseParameters,
                        coreMidNodeCount, coreMidNodeIds, coreMidRimNodeIds, midColumnNodeIds, nodeIds, midNodeIds,
                        ringNodeIds, connectionCounts)

            # Create nodes for non-core layers
            else:
                for n1 in range(elementsCountAroundHalf + 1):
                    cross1 = n1 == 0
                    cross2 = n1 == elementsCountAroundHalf
                    if s > 0:
                        n3p = n3 + (ellipseParameters[2] - 1) if isCore else n3
                        if cross1:
                            ringNodeIds.append(midNodeIds[0][n3p][0])
                            continue
                        if cross2:
                            ringNodeIds.append(midNodeIds[0][n3p][-1])
                            continue
                    node = nodes.createNode(nodeIdentifier, nodetemplateCross if (cross1 or cross2) else nodetemplate)
                    fieldcache.setNode(node)
                    if (not innerMidCoordinates) or (n3 == elementsCountThroughWall):
                        rx, rd1, rd2 = omx[n1], omd1[n1], omd2[n1]
                        rd3 = crossproduct3(rd1, rd2) if isCore else None
                    elif n3 == 0:
                        rx, rd1, rd2 = imx[n1], imd1[n1], imd2[n1]
                    else:
                        rx = add(mult(omx[n1], oFactor), mult(imx[n1], iFactor))
                        rd1 = add(mult(omd1[n1], oFactor), mult(imd1[n1], iFactor))
                        rd2 = add(mult(omd2[n1], oFactor), mult(imd2[n1], iFactor))
                        rd3 = crossproduct3(rd1, rd2) if isCore else None

                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2)
                    if isCore:
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3)

                    if omd12 and not (cross1 or cross2):
                        rd12 = omd12[n1] if ((not innerMidCoordinates) or (n3 == elementsCountThroughWall)) else \
                            imd12[n1] if (n3 == 0) else add(mult(omd12[n1], oFactor), mult(imd12[n1], iFactor))
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, rd12)

                    ringNodeIds.append(nodeIdentifier)
                    if isCore:
                        nodeParametersList.append([nodeIdentifier, rx, rd1, rd2, rd3])

                    nodeIdentifier = nodeIdentifier + 1

                midNodeIds[s].append(ringNodeIds)

    # Create elements
    # element template for bicubic hermite
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

    # Creating element templates for cross elements
    eftCrossForward = []
    eftCrossForwardScalefactors = []

    for s in range(segmentCount):
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

    for s in range(segmentCount):
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

    for s in range(segmentCount):
        # Create elements for the core
        if isCore:
            # Forward connections
            for e3 in range(2):
                if e3 == 0:
                    # Core elements
                    pass
                    for e1 in range(coreElementsCountHalf):
                        nids, nodeParameters = determineMidCoreElementsNids(e1, e3, s, ellipseParameters, inward,
                                                coreElementsCountHalf, tubeNodeIds, coreMidNodeIds, nodeParametersList)
                        elementIdentifier = generateMidCoreElements(s, mesh, coordinates, nids, nodeParameters,
                                                                    elementIdentifier, inward)
                elif e3 == 1:
                    # Rim elements
                    pass
                    for e1 in range(connectionCounts[s]):
                        nids, nodeParameters = determineMidCoreRimElementsNids(e1, e3, s, inward, crossIndexes,
                                aroundCounts, tubeNodeIds, coreMidRimNodeIds, midNodeIds, nodeParametersList)
                        elementIdentifier = generateMidCoreRimElements(e1, s, connectionCounts, mesh, coordinates, nids,
                                                                       nodeParameters, elementIdentifier)
            # Shell elements
            elementsCountThrough = elementsCountThroughWall + (ellipseParameters[2] - 1)
            for e3 in range(elementsCountThrough):
                for e1 in range(connectionCounts[s]):
                    nids, nodeParameters = determineMidCoreOuterShellNids(s, e1, e3, inward, crossIndexes, aroundCounts,
                            connectionCounts, tubeNodeIds, midNodeIds, nodeParametersList)
                    elementIdentifier = generateMidCoreOuterShellElements(e1, s, connectionCounts, mesh, coordinates,
                                                                          nids, nodeParameters, elementIdentifier)
            # Reverse connections
            for e3 in range(2):
                isForward = False
                if e3 == 0:
                    pass
                    for e1 in range(coreElementsCountHalf):
                        nids, nodeParameters = determineMidCoreElementsNids(e1, e3, s, ellipseParameters, inward,
                                                coreElementsCountHalf, tubeNodeIds, coreMidNodeIds, nodeParametersList,
                                                isForward)
                        elementIdentifier = generateMidCoreElements(s, mesh, coordinates, nids, nodeParameters,
                                                                    elementIdentifier, inward, isForward)
                elif e3 == 1:
                    pass
                    for e1 in range(connectionCounts[s]):
                        nids, nodeParameters = determineMidCoreRimElementsNids(e1, e3, s, inward, crossIndexes,
                                aroundCounts, tubeNodeIds, coreMidRimNodeIds, midNodeIds, nodeParametersList,
                                isForward)
                        elementIdentifier = generateMidCoreRimElements(e1, s, connectionCounts, mesh, coordinates, nids,
                                                                       nodeParameters, elementIdentifier, isForward)
            # Shell elements
            elementsCountThrough = elementsCountThroughWall + (ellipseParameters[2] - 1)
            for e3 in range(elementsCountThrough):
                for e1 in range(connectionCounts[s]):
                    nids, nodeParameters = determineMidCoreOuterShellNids(s, e1, e3, inward, crossIndexes, aroundCounts,
                            connectionCounts, tubeNodeIds, midNodeIds, nodeParametersList, isForward)
                    elementIdentifier = generateMidCoreOuterShellElements(e1, s, connectionCounts, mesh, coordinates,
                                                                          nids, nodeParameters, elementIdentifier, isForward)

        # Non-core elements
        # forward connections
        else:
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
                        n3p = n3 + 2 if isCore else n3
                        if inward[s]:
                            nStart = crossIndexes[s] - aroundCounts[s]
                            nids += [tubeNodeIds[s][n3p][nStart + e1], tubeNodeIds[s][n3p][nStart + e1 + 1],
                                     midNodeIds[s][n3][e1], midNodeIds[s][n3][e1 + 1]]
                        else:
                            nStart = crossIndexes[s] - connectionCounts[s]
                            re1 = connectionCounts[s] - e1
                            nids += [midNodeIds[s][n3][re1], midNodeIds[s][n3][re1 - 1],
                                     tubeNodeIds[s][n3p][nStart + e1], tubeNodeIds[s][n3p][nStart + e1 + 1]]

                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result1 = element.setNodesByIdentifier(eft, nids)
                    result2 = element.setScaleFactors(eft, scalefactors) if scalefactors else "-"
                    # print('create element tube bifurcation forward s', s, elementIdentifier, element.isValid(), result1, result2, nids)
                    for meshGroup in annotationMeshGroups[s]:
                        meshGroup.addElement(element)
                    elementIdentifier += 1

            # Non-core elements
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
                        n3p = n3 + 2 if isCore else n3
                        if inward[s]:
                            nStart = crossIndexes[s] - connectionCounts[s - 1]
                            re1 = connectionCounts[s - 1] - e1
                            nids += [tubeNodeIds[s][n3p][nStart + e1], tubeNodeIds[s][n3p][nStart + e1 + 1],
                                     midNodeIds[s - 1][n3][re1], midNodeIds[s - 1][n3][re1 - 1]]
                        else:
                            nStart = crossIndexes[s] - aroundCounts[s]
                            nids += [midNodeIds[s - 1][n3][e1], midNodeIds[s - 1][n3][e1 + 1],
                                     tubeNodeIds[s][n3p][nStart + e1], tubeNodeIds[s][n3p][nStart + e1 + 1]]

                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result1 = element.setNodesByIdentifier(eft, nids)
                    result2 = element.setScaleFactors(eft, scalefactors) if scalefactors else "-"
                    # print('Tube bifurcation element reverse s', s, elementIdentifier, element.isValid(), result1, result2, nids)
                    for meshGroup in annotationMeshGroups[s]:
                        meshGroup.addElement(element)
                    elementIdentifier += 1

    return nodeIdentifier, elementIdentifier


def getCoreMidCentre(outerMidCoordinates, innerMidCoordinates):
    """
    Get the centre coordinates of the core at the bifurcation using outer and inner bifurcation coordinates.
    :param outerMidCoordinates: A list of outer bifurcation coordinates.
    :param innerMidCoordinates: A list of inner bifurcation coordinates.
    :return: Centre coordinates for the core within the bifurcation.
    """
    P0 = []
    P1 = []

    P0.append(outerMidCoordinates[0][0][0])
    P0.append(outerMidCoordinates[0][0][-1])
    P1.append(innerMidCoordinates[0][0][0])
    P1.append(innerMidCoordinates[0][0][-1])

    P0 = np.array(P0)
    P1 = np.array(P1)

    # generate all line direction vectors
    n = (P1 - P0) / np.linalg.norm(P1 - P0, axis=1)[:, np.newaxis]  # normalized

    # generate the array of all projectors
    projs = np.eye(n.shape[1]) - n[:, :, np.newaxis] * n[:, np.newaxis]  # I - n*n.T

    # generate R matrix and q vector
    R = projs.sum(axis=0)
    q = (projs @ P0[:, :, np.newaxis]).sum(axis=0)

    # solve the least squares problem for the intersection point p: Rp = q
    p = np.linalg.lstsq(R, q, rcond=None)[0]

    Px = p[0][0]
    Py = p[1][0]
    Pz = p[2][0]

    centre = [Px, Py, Pz]

    return centre


def findMidCoreCoordinates(coreCentre, innerMidCoordinates, ellipseParameters, segmentCount, inward):
    """
    Blackbox function for generating core coordinates and derivatives at the mid-section of a bifurcation
    using Ellipse2D.
    Steps are:
    1. Find direction vector of the 1D network, centre points of the tube, and major/minor axes for the 2D ellipse.
    2. Generate 2D ellipses.
    3. Extract coordinates and derivatives from 2D ellipses generated in the previous step.
    4. Adjust (rotate and/or scale) derivatives to get the correct direction and magnitude.
    :param coreCentre: Coordinates for the core centre within the bifurcation.
    :param innerMidCoordinates: List of 3 middle half rings between tubes 1-2, 2-3 and 3-1 of inner coordinates
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param segmentCount: An identifier tracking tube segments.
    :param inward: List over 3 tubes of True if inward, False if not. All inward tubes precede outward.
    :return: Core coordinates and derivatives - coremx, coremd1, coremd2, coremd3
    """
    midMajorAxes = []
    midMinorAxes = []
    coreMidMajorRadii = []
    coreMidMinorRadii = []

    for s in range(segmentCount):
        imx, _, imd2 = innerMidCoordinates[s][:3]

        # Calculate midMajorRadius and midMinorRadius
        midMajorRadius, midMinorRadius = findEllipseRadii(imx, coreCentre, mode=2)

        # Calculate major and minor axes for mid section
        majorAxis, minorAxis = findEllipseAxes(imx, coreCentre, mode=2)
        midMajorAxes.append(majorAxis)
        midMinorAxes.append(minorAxis)

        coreMidMajorRadii.append(midMajorRadius)
        coreMidMinorRadii.append(midMinorRadius)

    # Generate ellipses for mid-section
    ellipses = []
    for s in range(segmentCount):
        ellipse = Ellipse2D(coreCentre, midMajorAxes[s], midMinorAxes[s], ellipseParameters[0], ellipseParameters[1],
                            0, ellipseParameters[2], 1.0, coreMidMajorRadii[s],
                            coreMidMinorRadii[s], isCore=True)
        ellipses.append(ellipse)

    # Extract coordinates and derivatives
    x, d1, d2, d3 = ([] for i in range(4))
    for s in range(segmentCount):
        [lst.append([]) for lst in [x, d1, d2, d3]]
        for n2 in range(ellipseParameters[0] + 1):
            for n1 in range(ellipseParameters[1] + 1):
                if ellipses[s].px[n2][n1]:
                    [lst.append(value) for lst, value in zip([x[s], d1[s], d2[s], d3[s]],
                                                             [ellipses[s].px[n2][n1], ellipses[s].pd1[n2][n1],
                                                              ellipses[s].pd2[n2][n1], ellipses[s].pd3[n2][n1]])]

    coremx, coremd1, coremd2, coremd3 = ([] for i in range(4))
    for s in range(segmentCount):
        endIndex = len(x[s]) if s == 0 else (len(x[s]) - ellipseParameters[1] - 1)
        [x.append(y) for x, y in zip([coremx, coremd1, coremd2, coremd3],
                                     [x[s][0:endIndex], d1[s][0:endIndex], d2[s][0:endIndex], d3[s][0:endIndex]])]

    for s in range(segmentCount):
        for n2 in range(len(coremd2[s])):
            d2 = coremd2[s][n2]
            _, _, imd2 = innerMidCoordinates[s][:3]
            scalefactor = vector.magnitude(imd2[len(imd2) // 2])
            d2 = vector.scaleVector(d2, -scalefactor)
            coremd2[s][n2] = d2

    return coremx, coremd1, coremd2, coremd3


def findCoreTubeCoordinates(coreTubeNodeIds, otx, otd1, itx, ellipseParameters, elementsCountAlong, outerTubeData,
                            segmentCount, coreNodeCount):
    """
    Blackbox function for generating core coordinates and derivatives at the starting ends of child tubes that connect
    with the bifurcation using Ellipse2D.
    Steps are:
    1. Find direction vector of the 1D network, centre points of the tube, and major/minor axes for the 2D ellipse.
    2. Generate 2D ellipses.
    3. Extract coordinates and derivatives from 2D ellipses generated in the previous step.
    4. Adjust (rotate and/or scale) derivatives to get the correct direction and magnitude.
    :param coreTubeNodeIds: A list of core node ids within a straight tube.
    :param otx, otd1, itx: outer tube coordinates, outer tube d1 derivatives, and inner tube coordinates, respectively.
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param elementsCountAlong: Number of elements along the tube.
    :param outerTubeData: Segment tube data containing outer tube path parameters.
    :param segmentCount: An identifier tracking tube segments.
    :param coreNodeCount: Number of nodes that form the solid core, including the rim layer.
    :return: Core coordinates and derivatives - ctx, ctd1, ctd2, ctd3
    """
    coreTubeNodeIds.append([])

    # find tube core coordinates and derivatives
    centre = [0, 0, 0]
    centre[0], centre[1], centre[2] = findCoreCentre(otx, itx)

    # Find directionVector of a tube
    directionVector = findDirectionVector(outerTubeData)

    for lft in [0, 1]:
        majorRadius, minorRadius = findEllipseRadii(itx, centre, lft)
        majorAxis, minorAxis = findEllipseAxes(itx, centre, lft, majorRadius, minorRadius)

        ellipse = Ellipse2D(centre, majorAxis, minorAxis, ellipseParameters[0], ellipseParameters[1],
                            0, ellipseParameters[2], 1.0, majorRadius, minorRadius,
                            isCore=True)
        if lft == 0:
            lftEllipse = ellipse
        elif lft == 1:
            rhtEllipse = ellipse

    ctx, ctd1, ctd2, ctd3 = extractCoordinatesFromEllipse(ellipseParameters, lftEllipse, rhtEllipse)

    ctd1 = adjustD1Derivatives(ctd1, ellipseParameters, coreNodeCount, otd1, segmentCount)
    ctd2 = adjustD2Derivatives(ctd2, ellipseParameters[3], directionVector)

    return ctx, ctd1, ctd2, ctd3


def adjustMidCoreDerivatives(elementsCountAcrossMinor, cmd1, cmd2, imd1, imd2, s, inward):
    """
    Rotates d2 derivatives of the solid core within the bifurcation to match d2 derivatives of surrounding inner layer.
    :param elementsCountAcrossMinor: Number of elements across minor axis of an ellipse.
    :param cmd1, cmd2, imd1, imd2: D1 and D2 derivatives for the core and inner bifurcation, respectively.
    :param s: Segment count.
    :param inward: List over 3 tubes of True if inward, False if not.
    :return: d2 derivatives of solid core within bifurcation - cmd2
    """
    segmentCount = len(inward)
    if segmentCount == 3 and inward[1] == True:
        isConverging = True
    elif segmentCount == 4 and inward[1] == True and inward[2] == True:
        isConverging = True
    else:
        isConverging = False

    coreMidNodeCount = len(cmd1)

    # Adjust d1 & d2 derivatives
    end = coreMidNodeCount - 1
    start = end - elementsCountAcrossMinor

    # Adjust d1 derivatives to fit the right-hand rule convention
    for n2 in range(0, elementsCountAcrossMinor-1):
        cmd1[n2] = vector.scaleVector(cmd1[n2], -1)

    # # d1 derivatives
    # if s == 0:
    #     for n in [start, end]:
    #         d1 = cmd1[n]
    #         if not isConverging:
    #             pos = -1 if n == start else 0
    #         else:
    #             pos = 0 if n == start else -1
    #         cmd1[n] = imd1[pos]

    # d2 derivatives
    for n2 in range(len(cmd2)):
        d2 = cmd2[n2]
        if s == 0:
            if n2 == start:
                pos1 = len(imd2) - 1
            elif n2 == end:
                pos1 = 0
            else:
                pos1 = 1
        else:
            pos1 = 2
        if vector.magnitude(imd2[pos1]) != 0:
            angleD2 = vector.angleBetweenVectors(d2, imd2[pos1])
        else:
            angleD2 = 0

        if s == 0:
            k = [0, 0, 1] if n2 != start else [0, 0, -1]
        else:
            if not isConverging:
                k = [0,0,-1]
            else:
                k = [0, 0, 1] if not inward[s] else [0, 0, -1]
        scaleFactor = 1

        if angleD2 > 0:
            rotateD2 = vector.rotateVectorAroundVector(d2, k, angleD2 * scaleFactor)
            cmd2[n2] = rotateD2

    if s == 0:
        for n in range(start, end + 1):
        # for n in range(start + 1, end):
            v = cmd1[n]
            k = [0,1,0]
            mirror_d2 = [(v[c] - 2 * vector.dotproduct(v, k) * k[c]) for c in range(3)]
            cmd2[n] = mirror_d2

    return cmd1, cmd2


def getCoreMidNodeIds(s, ellipseParameters, coreMidNodeCount, coreMidNodeIds, coreMidRimNodeIds, midColumnNodeIds,
                      nodeIds, midNodeIds, ringNodeIds, connectionCounts):
    """
    Sorts and rearranges node ids and returns four sets of node id lists required for generating elements at the
    mid-section of a bifurcation.
    :param s: Segment count.
    :param ellipseParameters: A list of parameters required to form an ellipse.
    :param coreMidNodeCount: Number of nodes that form the solid core, including the rim layer at the mid-section.
    :param coreMidNodeIds, coreMidRimNodeIds: Lists of core node ids at the mid-section.
    :param midColumnNodeIds: A list of core node ids that form a column through the centre of a bifurcation.
    :param nodeIds: A flat list of all node ids, excluding nodes at the outer layer.
    :param midNodeIds, ringNodeIds: Lists of node ids that are used to form the shell surrounding the core.
    :param connectionCounts: A list of number of connecting elements for each segment pair, 1-2, 2-3, and 3-1.
    :return: Four sets of node ids lists - coreMidNodeIds, coreMidRimNodeIds, midColumnNodeIds, midNodeIds
    """
    # Rearrange nodes for generating elements
    sortedMidNodeIndices = findRingNodeIndicesForSorting(coreMidNodeCount, ellipseParameters,
                                                         halfEllipse=True)
    tmpNodeIds, coreMidNodeIds = sortCoreNodeIds(sortedMidNodeIndices, coreMidNodeCount, coreMidNodeIds,
                                                 nodeIds)
    if s == 0:
        nColumnNodes = ellipseParameters[1] + 1
        midColumnNodeIds = nodeIds[-nColumnNodes::]
    else:
        idx1 = ellipseParameters[2]
        coreMidNodeIds[s] += midColumnNodeIds[idx1:-idx1]
        em = connectionCounts[s] - 1

        for i in range(idx1):
            value = midColumnNodeIds[i]
            # idx2 = (ellipseParameters[0] - 1 - 2 * (ellipseParameters[2] - 1)) * (i + 1) + i
            idx2 = em + (em + 1) * i
            tmpNodeIds.insert(idx2, value)
        for i in range(idx1):
            value = midColumnNodeIds[-(i + 1)]
            # idx2 = (ellipseParameters[0] - 2 * (ellipseParameters[2] - 1)) * (i + 1) + i
            idx2 = (em + 1) + (em + 2) * i
            tmpNodeIds.insert(idx2, value)
    for k in range(ellipseParameters[2]):
        ringNodeIds.append([])
        for i in range(connectionCounts[s] + 1):
            ringNodeIds[k].insert(i, tmpNodeIds.pop(0))
    del tmpNodeIds

    for k in range(ellipseParameters[2]):
        ringNodeIds[k] = rearrangeMidNodeIds(ringNodeIds[k], ellipseParameters)
    for k in range((ellipseParameters[2] - 1), -1, -1):
        midNodeIds[-1].append(ringNodeIds[k])

    # Create a new list for mid core rim nodes
    tmpRimIds = createCoreRimIdsList(coreMidNodeIds, ellipseParameters, mid=True)
    rearrangedRimIds = rearrangeMidNodeIds(tmpRimIds, ellipseParameters)
    coreMidRimNodeIds.append(rearrangedRimIds)

    return coreMidNodeIds, coreMidRimNodeIds, midColumnNodeIds, midNodeIds


def determineMidCoreElementsNids(e1, e3, s, ellipseParameters, inward, coreElementsCountHalf,
                                 tubeNodeIds, coreMidNodeIds, nodeParametersList, isForward=True):
    """
    Determines a set of 8 node ids required to form a solid core element, and searches for node parameters matching
    the selected set of node ids from nodeParametersList.
    :param e1, e2, s: index for elementsCountAround, elementsCountAlong, and segmentCount, respectively.
    :param ellipseParameters: A list of parameters needed to generate an ellipse.
    :param inward: List over 3 tubes of True if inward, False if not.
    :param coreElementCountHalf: Half the number of elements that form the solid core.
    :param tubeNodeIds, coreMidNodeIds: A list of tube node ids and core node ids, respectively.
    :param nodeParametersList: A list of node parameters.
    :param isForward: True if the connection is in the forward direction, False if in reverse.
    :return: a set of node ids and their node parameter values
    """
    segmentCount = len(inward)
    if segmentCount == 3 and inward[1] == True:
        isConverging = True
    elif segmentCount == 4 and inward[1] == True and inward[2] == True:
        isConverging = True
    else:
        isConverging = False

    e1p = int(e1 + e1 // (ellipseParameters[1] - 2 * ellipseParameters[2]))
    skipIndex = ellipseParameters[1] - 1 - 2 * (ellipseParameters[2] - 1)
    eS = 2 * (ellipseParameters[2] - 1)

    # nStart & cStart
    if isForward:
        # Forward connection
        cStart = -(ellipseParameters[1] - 1 - eS) if not isConverging else 0
        if inward[s]:
            if (ellipseParameters[0] - eS) > 4:
                nStart = int((ellipseParameters[1] - 1 - eS) * (ellipseParameters[0] // 2 - ellipseParameters[2]))
                if not isConverging:
                    cStart += int(cStart * (e1 // (ellipseParameters[1] - 2 * ellipseParameters[2])))
                else:
                    cStart += int((ellipseParameters[1] - 1 - eS) * (e1 // (coreElementsCountHalf / 2)))
                e1c = e1 % (coreElementsCountHalf // (ellipseParameters[0] // 2 - ellipseParameters[2]))
            else:
                nStart = int(ellipseParameters[1] - 1 - eS)
                e1c = e1
            if isConverging:
                nStart = 0
        else:
            # if isConverging:
            #     if (ellipseParameters[0] - eS) > 4:
            #         nStart = int(-(ellipseParameters[1] - 1 - eS))
            #     else:
            #         nStart = int(-(ellipseParameters[1] - 1 - eS))
            # else:
            #     nStart = 0
            # cStart = 0
            # skipIndex = -skipIndex if isConverging else skipIndex
            # e1c = (e1 % (coreElementsCountHalf // 2)) if (ellipseParameters[0] - eS) > 4 else e1
            cStart = -(ellipseParameters[1] - 1 - eS)
            skipIndex = skipIndex if isConverging else -skipIndex
            e1c = (e1 % (coreElementsCountHalf // 2)) if (ellipseParameters[0] - eS) > 4 else e1
            if (ellipseParameters[0] - eS) > 4:
                cStart += int(cStart * (e1 // (coreElementsCountHalf / 2)))
            nStart = int((ellipseParameters[1] - 1 - eS) * (ellipseParameters[0] // 2 - ellipseParameters[2]))

    else:
        # Reverse connection
        if inward[s]:
            if not isConverging:
                nStart = 0
                cStart = 0
                if (ellipseParameters[0] - eS) > 4:
                    cStart += int((ellipseParameters[1] - 1 - eS) * (e1 // (coreElementsCountHalf / 2)))
                e1c = (e1 % (coreElementsCountHalf // 2)) if (ellipseParameters[0] - eS) > 4 else e1
            else:
                nStart = int((ellipseParameters[1] - 1 - eS) * (ellipseParameters[0] // 2 - ellipseParameters[2]))
                cStart = -(ellipseParameters[1] - 1 - eS)
                if (ellipseParameters[0] - eS) > 4:
                    cStart += int(cStart * (e1 // (coreElementsCountHalf / 2)))
                e1c = (e1 % (coreElementsCountHalf // 2)) if (ellipseParameters[0] - eS) > 4 else e1
                skipIndex = -skipIndex
        else:
            cStart = -(ellipseParameters[1] - 1 - eS) if not isConverging else 0
            skipIndex = -skipIndex if isConverging else skipIndex
            e1c = (e1 % (coreElementsCountHalf // 2)) if (ellipseParameters[0] - eS) > 4 else e1
            if (ellipseParameters[0] - eS) > 4:
                cStart += int(cStart * (e1 // (coreElementsCountHalf / 2))) if not isConverging \
                    else int((ellipseParameters[1] - 1 - eS) * (e1 // (coreElementsCountHalf / 2)))
                nStart = int((ellipseParameters[1] - 1 - eS) * (ellipseParameters[0] - 3 - ellipseParameters[2]))
            else:
                nStart = int(ellipseParameters[1] - 1 - eS)
            if isConverging:
                nStart = 0

    # nids
    if isForward:
        # Foward connection
        scalefactor = -1 if isConverging else 1
        if inward[s]:
            nids = [tubeNodeIds[s][e3][nStart + e1p], tubeNodeIds[s][e3][nStart + e1p + skipIndex],
                    coreMidNodeIds[s][cStart + e1c], coreMidNodeIds[s][cStart + e1c - skipIndex*scalefactor],
                    tubeNodeIds[s][e3][nStart + e1p + 1], tubeNodeIds[s][e3][nStart + e1p + skipIndex + 1],
                    coreMidNodeIds[s][cStart + e1c + 1], coreMidNodeIds[s][cStart + e1c - skipIndex*scalefactor + 1]]
        else:
            # nids = [coreMidNodeIds[s][cStart + e1p], coreMidNodeIds[s][cStart + e1p + skipIndex*scalefactor],
            #         tubeNodeIds[s][e3][nStart + e1p], tubeNodeIds[s][e3][nStart + e1p + skipIndex],
            #         coreMidNodeIds[s][cStart + e1p + 1], coreMidNodeIds[s][cStart + e1p + skipIndex*scalefactor + 1],
            #         tubeNodeIds[s][e3][nStart + e1p + 1], tubeNodeIds[s][e3][nStart + e1p + skipIndex + 1]]

            nids = [coreMidNodeIds[s][cStart+e1c], coreMidNodeIds[s][cStart+e1c-skipIndex],
                    tubeNodeIds[s][e3][nStart+e1p], tubeNodeIds[s][e3][nStart+e1p+skipIndex],
                    coreMidNodeIds[s][cStart+e1c+1], coreMidNodeIds[s][cStart+e1c-skipIndex+1],
                    tubeNodeIds[s][e3][nStart+e1p+1], tubeNodeIds[s][e3][nStart+e1p+skipIndex+1]]
    else:
        # Reverse connection
        scalefactor = -1 if isConverging else 1
        if inward[s]:
            nids = [tubeNodeIds[s][0][nStart+e1p], tubeNodeIds[s][0][nStart+e1p+skipIndex*scalefactor],
                    coreMidNodeIds[s-1][cStart+e1c], coreMidNodeIds[s-1][cStart+e1c+skipIndex],
                    tubeNodeIds[s][0][nStart+e1p+1], tubeNodeIds[s][0][nStart+e1p+skipIndex*scalefactor+1],
                    coreMidNodeIds[s-1][cStart+e1c+1], coreMidNodeIds[s-1][cStart+e1c+skipIndex+1]]
        else:
            nids = [coreMidNodeIds[s-1][cStart+e1c], coreMidNodeIds[s-1][cStart+e1c-skipIndex],
                    tubeNodeIds[s][0][nStart+e1p], tubeNodeIds[s][0][nStart+e1p+skipIndex*scalefactor],
                    coreMidNodeIds[s-1][cStart+e1c+1], coreMidNodeIds[s-1][cStart+e1c-skipIndex+1],
                    tubeNodeIds[s][0][nStart+e1p+1], tubeNodeIds[s][0][nStart+e1p+skipIndex*scalefactor+1]]

    nodeParameters = []
    for id in nids:
        for n in range(len(nodeParametersList)):
            if id == nodeParametersList[n][0]:
                nodeParameters.append(nodeParametersList[n][1:])

    return nids, nodeParameters


def generateMidCoreElements(s, mesh, coordinates, nids, nodeParameters, elementIdentifier, inward, isForward=True):
    """
    First determines the required eft and scalefactors using determineTricubicHermite function based on nodeParameters,
    and then generates a core element at the mid-section of a bifurcation for a given set of node ids.
    :param s: segment count.
    :param mesh:  A Zinc mesh of dimension 3.
    :param coordinates: Finite element coordinate field to define.
    :param nids: A set of 8 node ids that forms an element.
    :param nodeParameters: A list of node parameters for given set of nids.
    :param elementIdentifier: First 2D element identifier to use.
    :param segmentCount:
    :param inward:
    :param isForward:
    :return: elementIdentifier for the next element to be generated.
    """
    segmentCount = len(inward)
    if segmentCount == 3 and inward[1] == True:
        isConverging = True
    elif segmentCount == 4 and inward[1] == True and inward[2] == True:
        isConverging = True
    else:
        isConverging = False

    # Create element template
    if s == 0:
        if not isConverging:
            if isForward:
                nodeDerivativeFixedWeights = [[], [], [None, [1.0, 1.0, 0.0]], [], [], [], [None, [1.0, 1.0, 0.0]], []]
            else:
                nodeDerivativeFixedWeights = [[], [], [], [None, [1.0, 1.0, 0.0]], [], [], [], [None, [1.0, 1.0, 0.0]]]
        else:
            nodeDerivativeFixedWeights = None
        eftCore, scalefactorsCore = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                      serendipity=True)
    elif s == 1:
        if not isConverging:
            if isForward:
                nodeDerivativeFixedWeights = [[[1.0, 1.0, 0.0]], [], [], [], [[1.0, 1.0, 0.0]], [], [], []]
            else:
                nodeDerivativeFixedWeights = None
        else:
            nodeDerivativeFixedWeights = None
        eftCore, scalefactorsCore = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                      serendipity=True)
    else:
        if not isConverging:
            if isForward:
                nodeDerivativeFixedWeights = None
            else:
                nodeDerivativeFixedWeights = [[[1.0, 1.0, 0.0]], [], [], [], [[1.0, 1.0, 0.0]], [], [], []]
        else:
            nodeDerivativeFixedWeights = None
        eftCore, scalefactorsCore = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                      serendipity=True)

    elementtemplateCore = mesh.createElementtemplate()
    elementtemplateCore.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    elementtemplateCore.defineField(coordinates, -1, eftCore)
    element = mesh.createElement(elementIdentifier, elementtemplateCore)
    element.setNodesByIdentifier(eftCore, nids)
    element.setScaleFactors(eftCore, scalefactorsCore) if scalefactorsCore else "-"

    elementIdentifier += 1

    return elementIdentifier


def determineMidCoreRimElementsNids(e1, e3, s, inward, crossIndexes, aroundCounts, tubeNodeIds, coreMidRimNodeIds,
                                    midNodeIds, nodeParametersList, isForward=True):
    """
    Determines a set of 8 node ids required to form a rim element around the solid core, and searches for node
    parameters matching the selected set of node ids from nodeParametersList.
    :param e1, e2, s: index for connectionCounts, elementsCountAlong, and segmentCount, respectively.
    :param inward: List over 3 tubes of True if inward, False if not.
    :param crossIndexes: A list of indexes that are at cross junction of a bifurcation.
    :param aroundCounts: Number of elements around bifurcation.
    :param tubeNodeIds, coreMidRimNodeIds, midNodeIds: A list of tube node ids, core rim node ids, and
     mid-section node ids, respectively.
    :param nodeParametersList: A list of node parameters.
    :param isForward: True if the connection is in the forward direction, False if in reverse.
    :return: a set of node ids and their node parameter values
    """
    if isForward:
        # Forward connection
        if inward[s]:
            nStart = crossIndexes[s] - aroundCounts[s]
            nids = [tubeNodeIds[s][e3][nStart+e1], tubeNodeIds[s][e3][nStart+e1+1],
                    coreMidRimNodeIds[s][e1], coreMidRimNodeIds[s][e1+1],
                    tubeNodeIds[s][e3+1][nStart+e1], tubeNodeIds[s][e3+1][nStart+e1+1],
                    midNodeIds[s][e3-1][e1], midNodeIds[s][e3-1][e1+1]]

        else:
            nStart = crossIndexes[s]
            nids = [coreMidRimNodeIds[s][e1+1], coreMidRimNodeIds[s][e1],
                    tubeNodeIds[s][e3][nStart-e1-1], tubeNodeIds[s][e3][nStart-e1],
                    midNodeIds[s][e3-1][e1+1], midNodeIds[s][e3-1][e1],
                    tubeNodeIds[s][e3+1][nStart-e1-1], tubeNodeIds[s][e3+1][nStart-e1]]
    else:
        # Reverse connection
        if inward[s]:
            nStart = crossIndexes[s] - aroundCounts[s]
            e1t = (nStart - e1) % aroundCounts[s]
            nids = [tubeNodeIds[s][e3][e1t], tubeNodeIds[s][e3][e1t - 1],
                    coreMidRimNodeIds[s - 1][e1], coreMidRimNodeIds[s - 1][e1 + 1],
                    tubeNodeIds[s][e3 + 1][e1t], tubeNodeIds[s][e3 + 1][e1t - 1],
                    midNodeIds[s - 1][e3 - 1][e1], midNodeIds[s - 1][e3 - 1][e1 + 1]]
        else:
            nStart = crossIndexes[s]
            e1t = (nStart + e1) % aroundCounts[s]
            # print("e1t", e1t)
            nids = [coreMidRimNodeIds[s - 1][e1], coreMidRimNodeIds[s - 1][e1 + 1],
                    tubeNodeIds[s][e3][e1t], tubeNodeIds[s][e3][(e1t + 1)],
                    midNodeIds[s - 1][e3 - 1][e1], midNodeIds[s - 1][e3 - 1][e1 + 1],
                    tubeNodeIds[s][e3 + 1][e1t], tubeNodeIds[s][e3 + 1][(e1t + 1)]]
        # print("s", s, "nids", nids)
    nodeParameters = []
    for id in nids:
        for n in range(len(nodeParametersList)):
            if id == nodeParametersList[n][0]:
                nodeParameters.append(nodeParametersList[n][1:])

    return nids, nodeParameters


def generateMidCoreRimElements(e1, s, connectionCounts, mesh, coordinates, nids, nodeParameters, elementIdentifier,
                               isForward=True):
    """
    First determines the required eft and scalefactors using determineTricubicHermite function based on nodeParameters,
    and then generates a core rim element for a given set of node ids.
    :param e1, s: index for elementsCountAround and segment count, respectively.
    :param connectionCounts: A list of number of connecting elements for each segment pair, 1-2, 2-3, and 3-1.
    :param mesh:  A Zinc mesh of dimension 3.
    :param coordinates: Finite element coordinate field to define.
    :param nids: A set of 8 node ids that forms an element.
    :param nodeParameters: A list of node parameters for given set of nids.
    :param elementIdentifier: First 2D element identifier to use.
    :param isForward: True if the connection is in the forward direction, False if in reverse.
    :return: elementIdentifier for the next element to be generated.
    """

    # Forward elements
    if isForward:
        # Create element template
        if s == 0:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[], [None, None, [1.0, 0.0, 1.0]],
                                              [None, [1.0, 1.0, 0.0]], [None, None, [-1.0, 0.0, 1.0]],
                                              [], [], [None, [1.0, 1.0, 0.0]], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[None, None, [1.0, 0.0, 1.0]], [], [None, None, [-1.0, 0.0, 1.0]], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[], [None, None, [1.0, 0.0, -1.0]], [], [None, None, [-1.0, 0.0, -1.0]],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[None, None, [1.0, 0.0, -1.0]], [],
                                              [None, None, [-1.0, 0.0, -1.0]], [None, [1.0, 1.0, 0.0]],
                                              [], [], [], [None, [-1.0, -1.0, 0.0]]]
            else:
                nodeDerivativeFixedWeights = None

            eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters,
                    nodeDerivativeFixedWeights, serendipity=True)
        elif s == 1:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, 1.0]], [[-1.0, -1.0, 0.0]],
                                              [None, None, [-1.0, 0.0, 1.0]], [],
                                              [], [[-1.0, -1.0, 0.0]], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[], [None, None, [-1.0, 0.0, 1.0]], [], [None, None, [-1.0, 0.0, 1.0]],
                                              [[-1.0, 0.0, 0.0]], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, -1.0]], [], [None, None, [-1.0, 0.0, -1.0]], [],
                                              [], [[-1.0, 0.0, 0.0]], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[[1.0, 1.0, 0.0]], [None, None, [-1.0, 0.0, -1.0]],
                                              [], [None, None, [-1.0, 0.0, -1.0]],
                                              [[-1.0, -1.0, 0.0]], [], [], []]
            else:
                nodeDerivativeFixedWeights = None

            eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                          serendipity=True)
        else:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, 1.0]], [], [None, None, [-1.0, 0.0, 1.0]], [],
                                              [], [], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[], [None, None, [-1.0, 0.0, 1.0]], [], [None, None, [-1.0, 0.0, 1.0]],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, -1.0]], [], [None, None, [-1.0, 0.0, -1.0]], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[], [None, None, [-1.0, 0.0, -1.0]], [], [None, None, [-1.0, 0.0, -1.0]],
                                              [], [], [], []]
            else:
                nodeDerivativeFixedWeights = None
            eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                          serendipity=True)
    # Reverse elements
    else:
        if s == 0:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[], [None, None, [-1.0, 0.0, 1.0]], [], [None, None, [-1.0, 0.0, 1.0]],
                                              [], [], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, 1.0]], [], [None, None, [-1.0, 0.0, 1.0]], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[], [None, None, [-1.0, 0.0, -1.0]], [], [None, None, [-1.0, 0.0, -1.0]],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, -1.0]], [], [None, None, [-1.0, 0.0, -1.0]], [],
                                              [], [], [], [None, [-1.0, -1.0, 0.0]]]
            else:
                nodeDerivativeFixedWeights = None

            eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters,
                    nodeDerivativeFixedWeights, serendipity=True)

            setEftScaleFactorIds(eftCoreRim, [1], [])
            if e1 == 0:
                remapEftNodeValueLabel(eftCoreRim, [3,7], Node.VALUE_LABEL_D_DS1,
                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
            elif e1 == (connectionCounts[s] - 1):
                remapEftNodeValueLabel(eftCoreRim, [4], Node.VALUE_LABEL_D_DS1,
                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
        elif s == 1:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[], [None, None, [-1.0, 0.0, 1.0]], [], [None, None, [1.0, 0.0, 1.0]],
                                              [], [], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, 1.0]], [], [None, None, [1.0, 0.0, 1.0]], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[], [None, None, [-1.0, 0.0, -1.0]], [], [None, None, [1.0, 0.0, -1.0]],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, -1.0]], [], [None, None, [1.0, 0.0, -1.0]], [],
                                              [], [], [], []]
            else:
                nodeDerivativeFixedWeights = None

            eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters,
                    nodeDerivativeFixedWeights, serendipity=True)
        else:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[[1.0, 1.0, 0.0]], [None, None, [-1.0, 0.0, 1.0]],
                                              [], [None, None, [1.0, 0.0, 1.0]],
                                              [[1.0, 1.0, 0.0]], [], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, 1.0]], [], [None, None, [1.0, 0.0, 1.0]], [],
                                              [], [[1.0, 0.0, 0.0]], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[], [None, None, [-1.0, 0.0, -1.0]], [], [None, None, [1.0, 0.0, -1.0]],
                                              [[1.0, 0.0, 0.0]], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[None, None, [-1.0, 0.0, -1.0]], [[-1.0, -1.0, 0.0]],
                                              [None, None, [1.0, 0.0, -1.0]], [],
                                              [], [[1.0, 1.0, 0.0]], [], []]
            else:
                nodeDerivativeFixedWeights = None

            eftCoreRim, scalefactorsCoreRim = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                          serendipity=True)

    elementtemplateCoreRim = mesh.createElementtemplate()
    elementtemplateCoreRim.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    elementtemplateCoreRim.defineField(coordinates, -1, eftCoreRim)
    element = mesh.createElement(elementIdentifier, elementtemplateCoreRim)
    element.setNodesByIdentifier(eftCoreRim, nids)
    element.setScaleFactors(eftCoreRim, scalefactorsCoreRim) if scalefactorsCoreRim else "-"

    elementIdentifier += 1

    return elementIdentifier


def determineMidCoreOuterShellNids(s, e1, e3, inward, crossIndexes, aroundCounts, connectionCounts, tubeNodeIds, midNodeIds,
                              nodeParametersList, isForward=True):
    """
    Determines a set of 8 node ids required to form an outer shell element, and searches for node parameters matching
    the selected set of node ids from nodeParametersList.
    :param s, e1, e3: index for segmentCount, connectionCounts, and elementsCountThrough, respectively.
    :param inward: List over 3 tubes of True if inward, False if not.
    :param crossIndexes: A list of indexes that are at cross junction of a bifurcation.
    :param aroundCounts: Number of elements around bifurcation.
    :param connectionCounts: A list of number of connecting elements for each segment pair, 1-2, 2-3, and 3-1.
    :param tubeNodeIds, midNodeIds: A list of tube node ids, and mid-section node ids, respectively.
    :param nodeParametersList: A list of node parameters.
    :param isForward: True if the connection is in the forward direction, False if in reverse.
    :return: a set of node ids and their node parameter values
    """
    nids = []
    if isForward:
        for n3 in [e3, e3 + 1]:
            n3p = n3 + 2
            if inward[s]:
                nStart = crossIndexes[s] - aroundCounts[s]
                nids += [tubeNodeIds[s][n3p][nStart + e1], tubeNodeIds[s][n3p][nStart + e1 + 1],
                         midNodeIds[s][n3][e1], midNodeIds[s][n3][e1 + 1]]
            else:
                nStart = crossIndexes[s] - connectionCounts[s]
                re1 = connectionCounts[s] - e1
                nids += [midNodeIds[s][n3][re1], midNodeIds[s][n3][re1 - 1],
                         tubeNodeIds[s][n3p][nStart + e1], tubeNodeIds[s][n3p][nStart + e1 + 1]]
    else:
        for n3 in [e3, e3 + 1]:
            n3p = n3 + 2
            if inward[s]:
                nStart = crossIndexes[s] - connectionCounts[s - 1]
                re1 = connectionCounts[s - 1] - e1
                nids += [tubeNodeIds[s][n3p][nStart + e1], tubeNodeIds[s][n3p][nStart + e1 + 1],
                         midNodeIds[s - 1][n3][re1], midNodeIds[s - 1][n3][re1 - 1]]
            else:
                nStart = crossIndexes[s] - aroundCounts[s]
                nids += [midNodeIds[s - 1][n3][e1], midNodeIds[s - 1][n3][e1 + 1],
                         tubeNodeIds[s][n3p][nStart + e1], tubeNodeIds[s][n3p][nStart + e1 + 1]]

    nodeParameters = []
    for id in nids:
        for n in range(len(nodeParametersList)):
            if id == nodeParametersList[n][0]:
                nodeParameters.append(nodeParametersList[n][1:])

    return nids, nodeParameters


def generateMidCoreOuterShellElements(e1, s, connectionCounts, mesh, coordinates, nids, nodeParameters,
                                      elementIdentifier, isForward=True):
    """

    """
    # Forward elements
    if isForward:
        if s == 0:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[], [], [None, [1.0, 1.0, 0.0]], [None, None, [0.0, 0.0, 1.0]],
                                              [], [], [None, [1.0, 1.0, 0.0]], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[], [], [None, None, [0.0, 0.0, 1.0]], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[], [], [], [None, None, [0.0, 0.0, 1.0]],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[], [], [None, None, [0.0, 0.0, 1.0]], [None, [-1.0, -1.0, 0.0]],
                                              [], [], [], [None, [1.0, 1.0, 0.0]]]
            else:
                nodeDerivativeFixedWeights = None

            eftShell, scalefactorsShell = determineTricubicHermiteEft(mesh, nodeParameters,
                    nodeDerivativeFixedWeights, serendipity=True)
        elif s == 1:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[[-1.0, -1.0, 0.0]], [None, None, [0.0, 0.0, 1.0]], [], [],
                                              [[1.0, 1.0, 0.0]], [None, None, [0.0, 0.0, 1.0]], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[None, None, [0.0, 0.0, 1.0]], [[-1.0, 0.0, 0.0]], [], [],
                                              [None, None, [0.0, 0.0, 1.0]], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[[-1.0, 0.0, 0.0]], [None, None, [0.0, 0.0, 1.0]], [], [],
                                              [], [None, None, [0.0, 0.0, 1.0]], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[None, None, [0.0, 0.0, 1.0]], [[-1.0, -1.0, 0.0]], [], [],
                                              [None, None, [0.0, 0.0, 1.0]], [[-1.0, -1.0, 0.0]], [], []]
            else:
                nodeDerivativeFixedWeights = None

            eftShell, scalefactorsShell = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                          serendipity=True)
        else:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[], [None, None, [0.0, 0.0, 1.0]], [], [],
                                              [], [], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[None, None, [0.0, 0.0, 1.0]], [], [], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[], [None, None, [0.0, 0.0, 1.0]], [], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[None, None, [0.0, 0.0, 1.0]], [], [], [],
                                              [], [], [], []]
            else:
                nodeDerivativeFixedWeights = None
            eftShell, scalefactorsShell = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                          serendipity=True)
    # Reverse elements
    else:
        if s == 0:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[], [], [None, [-1.0, -1.0, 0.0]], [None, None, [0.0, 0.0, 1.0]],
                                              [], [], [None, [1.0, 1.0, 0.0]], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[], [], [None, None, [0.0, 0.0, 1.0]], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[], [], [], [None, None, [0.0, 0.0, 1.0]],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[], [], [None, None, [0.0, 0.0, 1.0]], [None, [1.0, 1.0, 0.0]],
                                              [], [], [], [None, [1.0, 1.0, 0.0]]]
            else:
                nodeDerivativeFixedWeights = None

            eftShell, scalefactorsShell = determineTricubicHermiteEft(mesh, nodeParameters,
                    nodeDerivativeFixedWeights, serendipity=True)
        elif s == 2:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[[1.0, 1.0, 0.0]], [None, None, [0.0, 0.0, 1.0]], [], [],
                                              [[1.0, 1.0, 0.0]], [None, None, [0.0, 0.0, 1.0]], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[None, None, [0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]], [], [],
                                              [None, None, [0.0, 0.0, 1.0]], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[[1.0, 0.0, 0.0]], [None, None, [0.0, 0.0, 1.0]], [], [],
                                              [], [None, None, [0.0, 0.0, 1.0]], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[None, None, [0.0, 0.0, 1.0]], [[1.0, 1.0, 0.0]], [], [],
                                              [None, None, [0.0, 0.0, 1.0]], [[-1.0, -1.0, 0.0]], [], []]
            else:
                nodeDerivativeFixedWeights = None

            eftShell, scalefactorsShell = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                          serendipity=True)
        else:
            if e1 == 0:
                nodeDerivativeFixedWeights = [[], [None, None, [0.0, 0.0, 1.0]], [], [],
                                              [], [], [], []]
            elif e1 == 1:
                nodeDerivativeFixedWeights = [[None, None, [0.0, 0.0, 1.0]], [], [], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 2):
                nodeDerivativeFixedWeights = [[], [None, None, [0.0, 0.0, 1.0]], [], [],
                                              [], [], [], []]
            elif e1 == (connectionCounts[s] - 1):
                nodeDerivativeFixedWeights = [[None, None, [0.0, 0.0, 1.0]], [], [], [],
                                              [], [], [], []]
            else:
                nodeDerivativeFixedWeights = None
            eftShell, scalefactorsShell = determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights,
                                                                          serendipity=True)

    elementtemplateShell = mesh.createElementtemplate()
    elementtemplateShell.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    elementtemplateShell.defineField(coordinates, -1, eftShell)
    element = mesh.createElement(elementIdentifier, elementtemplateShell)
    element.setNodesByIdentifier(eftShell, nids)
    element.setScaleFactors(eftShell, scalefactorsShell) if scalefactorsShell else "-"

    elementIdentifier += 1

    return elementIdentifier


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


def blendNetworkNodeCoordinates(networkNode, segmentTubeDataList):
    """
    Blend coordinate d2 between segments connecting at networkNode if sharing the same version.
    Must only call after all tube data has been resampled.
    :param networkNode: The node to blend coordinates at.
    :param segmentTubeDataList: List of dict mapping NetworkSegment to TubeData. Assumes all same structure.
    :return: None
    """
    inSegments = networkNode.getInSegments()
    outSegments = networkNode.getOutSegments()
    nodeVersionSegments = {}  # map from version number to list of segments using it
    for segment in (inSegments + outSegments):
        nodeVersions = segment.getNodeVersions()
        nodeIndex = -1 if (segment in inSegments) else 0
        nodeVersion = nodeVersions[nodeIndex]
        nodeVersionSegment = nodeVersionSegments.get(nodeVersion)
        if not nodeVersionSegment:
            nodeVersionSegments[nodeVersion] = nodeVersionSegment = []
        nodeVersionSegment.append(segment)
    for nodeVersion, segments in nodeVersionSegments.items():
        if len(segments) < 2:
            continue  # no blending required
        for segmentTubeData in segmentTubeDataList:
            d2Rings = []
            nodesCountAround = None
            for segment in segments:
                tubeData = segmentTubeData[segment]
                nodeIndex = -1 if segment in inSegments else 0
                d2Ring = tubeData.getSampledTubeCoordinates()[2][nodeIndex]
                if (not d2Rings) or (len(d2Ring) == nodesCountAround):
                    d2Rings.append(d2Ring)
                    nodesCountAround = len(d2Ring)
                else:
                    print("Cannot blend d1 version " + str(nodeVersion) + " at layout node " +
                          str(networkNode.getNodeIdentifier()) + " due to mismatch in number around")
                    break
            ringCount = len(d2Rings)
            if ringCount < 2:
                break
            for n in range(nodesCountAround):
                # harmonic mean magnitude; directions are the same as same version
                sum = 0.0
                ringMags = []
                for d2Ring in d2Rings:
                    ringMag = magnitude(d2Ring[n])
                    ringMags.append(ringMag)
                    if ringMag == 0.0:
                        sum = 0.0
                        break
                    sum += 1.0 / ringMag
                mag = (ringCount / sum) if (sum != 0.0) else 0.0
                for r in range(ringCount):
                    d2Ring = d2Rings[r]
                    if ringMags[r] > 0.0:
                        d2 = mult(d2Ring[n], mag / ringMags[r])
                        for c in range(3):
                            d2Ring[n][c] = d2[c]


def generateTubeBifurcationTree(networkMesh: NetworkMesh, region, coordinates, nodeIdentifier, elementIdentifier,
                                defaultElementsCountAround: int, targetElementDensityAlongLongestSegment: float,
                                elementsCountThroughWall: int, layoutAnnotationGroups: list = [],
                                annotationElementsCountsAround: list = [], annotationElementsCountsAcrossMajor: list = [],
                                ellipseParameters: list = [],
                                serendipity=False, showTrimSurfaces=False, isCore=False):
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
    :param ellipseParameters: A list of parameters needed to generate an ellipse using Ellipse2D.
    [elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossTransition, elementsCountAlong]
    :param serendipity: True to use Hermite serendipity basis, False for regular Hermite with zero cross derivatives.
    :param showTrimSurfaces: Set to True to make surfaces from inter-tube trim surfaces. For diagnostic use.
    :param isCore: True to generate a solid core inside the tube, False for regular tube network
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

    # make 2D annotation groups from 1D network layout annotation groups
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
    segmentCount = len(networkSegments)
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
            sx, sd1, sd2, sd12 = resampleTubeCoordinates(
                rawTubeCoordinates, elementsCountAlong, startSurface=startSurface, endSurface=endSurface)
            tubeData.setSampledTubeCoordinates((sx, sd1, sd2, sd12))
    del segmentTubeData

    # blend coordinates where versions are shared between segments
    blendedNetworkNodes = set()
    for networkSegment in networkSegments:
        segmentNodes = networkSegment.getNetworkNodes()
        for segmentNode in [segmentNodes[0], segmentNodes[-1]]:
            if segmentNode not in blendedNetworkNodes:
                blendNetworkNodeCoordinates(segmentNode, allSegmentTubeData)
                blendedNetworkNodes.add(segmentNode)
    del blendedNetworkNodes

    completedBifurcations = set()  # record so only done once
    segmentIdentifier = 1
    nodeParametersList = None

    with ChangeManager(fieldmodule):
        for networkSegment in networkSegments:
            segmentNodes = networkSegment.getNetworkNodes()
            startSegmentNode = segmentNodes[0]
            startInSegments = startSegmentNode.getInSegments()
            startOutSegments = startSegmentNode.getOutSegments()
            startSkipCount = 1 if ((len(startInSegments) + len(startOutSegments)) > 2) else 0
            endSegmentNode = segmentNodes[-1]
            endInSegments = endSegmentNode.getInSegments()
            endOutSegments = endSegmentNode.getOutSegments()
            endSkipCount = 1 if ((len(endInSegments) + len(endOutSegments)) > 2) else 0

            for stage in range(3):
                if stage == 1:
                    # tube
                    outerTubeData = outerSegmentTubeData[networkSegment]
                    outerTubeCoordinates = outerTubeData.getSampledTubeCoordinates()
                    innerTubeData = innerSegmentTubeData[networkSegment] if layoutInnerCoordinates else None
                    innerTubeCoordinates = innerTubeData.getSampledTubeCoordinates() if layoutInnerCoordinates else None
                    startNodeIds = outerTubeData.getStartNodeIds(startSkipCount)
                    if (not startNodeIds) and (startSkipCount == 0) and (startInSegments or startOutSegments):
                        # discover start nodes from single adjacent segment
                        if startInSegments:
                            startNodeIds = outerSegmentTubeData[startInSegments[0]].getEndNodeIds(0)
                        else:
                            startNodeIds = outerSegmentTubeData[startOutSegments[0]].getStartNodeIds(0)
                        if startNodeIds:
                            outerTubeData.setStartNodeIds(startNodeIds, startSkipCount)
                    endNodeIds = outerTubeData.getEndNodeIds(endSkipCount)
                    if (not endNodeIds) and (endSkipCount == 0) and (endOutSegments or endInSegments):
                        # discover end nodes from single adjacent segment
                        if endOutSegments:
                            endNodeIds = outerSegmentTubeData[endOutSegments[0]].getStartNodeIds(0)
                        elif endInSegments:
                            endNodeIds = outerSegmentTubeData[endInSegments[0]].getEndNodeIds(0)
                        if endNodeIds:
                            outerTubeData.setEndNodeIds(endNodeIds, endSkipCount)
                    loop = (len(startInSegments) == 1) and (startInSegments[0] is networkSegment) and \
                           (networkSegment.getNodeVersions()[0] == networkSegment.getNodeVersions()[-1])

                    nodeIdentifier, elementIdentifier, startNodeIds, endNodeIds, nodeParametersList = generateTube(
                        outerTubeCoordinates, innerTubeCoordinates, elementsCountThroughWall, region, fieldcache,
                        coordinates, nodeIdentifier, elementIdentifier, segmentIdentifier, innerTubeData,
                        ellipseParameters, nodeParametersList, startSkipCount=startSkipCount, endSkipCount=endSkipCount,
                        startNodeIds=startNodeIds, endNodeIds=endNodeIds,
                        annotationMeshGroups=outerTubeData.getAnnotationMeshGroups(), loop=loop,
                        serendipity=serendipity, isCore=isCore)
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

                    segmentIdentifier += 1
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
                        tubeNodeIds = [outerTubeData[s].getEndNodeIds(1) if inward[s] else
                                       outerTubeData[s].getStartNodeIds(1) for s in range(3)]

                        if isCore:
                            tmpNodeParametersList = []
                            for i in range(len(tubeNodeIds[0])):
                                for j in range(len(nodeParametersList)):
                                    if i == 1:
                                        pass
                                    else:
                                        if nodeParametersList[j][0] in tubeNodeIds[0][i]:
                                            tmpNodeParametersList.append(nodeParametersList[j])
                            nodeParametersList = tmpNodeParametersList

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
                            segmentIdentifier, outerTubeData,
                            annotationMeshGroups, ellipseParameters, nodeParametersList,
                            serendipity=serendipity, isCore=isCore)

                        for s in range(3):
                            if inward[s]:
                                if not outerTubeData[s].getEndNodeIds(1):
                                    outerTubeData[s].setEndNodeIds(tubeNodeIds[s], 1)
                            else:
                                if not outerTubeData[s].getStartNodeIds(1):
                                    outerTubeData[s].setStartNodeIds(tubeNodeIds[s], 1)

                        completedBifurcations.add(outerTubeBifurcationData)

    return nodeIdentifier, elementIdentifier, annotationGroups
