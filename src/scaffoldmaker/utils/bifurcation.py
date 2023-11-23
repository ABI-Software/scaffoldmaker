"""
Utilities for building bifurcating network meshes.
"""

from cmlibs.maths.vectorops import add, cross, dot, magnitude, mult, normalize, sub
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils.interpolation import computeCubicHermiteArcLength, computeCubicHermiteDerivativeScaling, \
    DerivativeScalingMode, getCubicHermiteCurvesLength, interpolateCubicHermite, interpolateCubicHermiteDerivative, \
    interpolateCubicHermiteSecondDerivative, smoothCubicHermiteDerivativesLine, interpolateLagrangeHermiteDerivative
from scaffoldmaker.utils.networkmesh import NetworkMesh, getPathRawTubeCoordinates, resampleTubeCoordinates
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import generateCurveMesh, get_nodeset_path_ordered_field_parameters, \
    print_node_field_parameters
import math
import copy


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

    def __init__(self, pathParameters):
        self._pathParameters = pathParameters
        self._rawTubeCoordinates = None
        self._rawTrackSurface = None
        self._sampledTubeCoordinates = None
        self._sampledNodeIds = []  # indexed along sample nodes. Only ever set on rows where junctions occur
        self._annotationMeshGroups = []

    def getPathParameters(self):
        return self._pathParameters

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

    def getSampledTubeCoordinates(self):
        return self._sampledTubeCoordinates

    def setSampledTubeCoordinates(self, sampledTubeCoordinates):
        """
        :param sampledTubeCoordinates: sx, sd1, sd2, sd12
        """
        self._sampledTubeCoordinates = sampledTubeCoordinates
        self._sampledNodeIds = [None] * len(self._sampledTubeCoordinates[0])

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

    def __init__(self, networkSegmentsIn: list, networkSegmentsOut: list, segmentTubeData,
                 outerTubeBifurcationData=None):
        """
        :param networkSegmentsIn: List of input segments.
        :param networkSegmentsOut: List of output segments.
        :param segmentTubeData: dict NetworkSegment -> SegmentTubeData.
        :param outerTubeData: Optional reference to outer TubeBifurcationData to inherit trim surfaces and
        cross indexes from.
        """
        self._networkSegmentsIn = networkSegmentsIn
        self._networkSegmentsOut = networkSegmentsOut
        self._networkSegments = networkSegmentsIn + networkSegmentsOut
        segmentCount = len(self._networkSegments)
        assert segmentCount == 3
        self._tubeData = [segmentTubeData[networkSegment] for networkSegment in self._networkSegments]
        self._segmentsIn = [self._networkSegments[s] in self._networkSegmentsIn for s in range(3)]
        # following are calculated in determineCrossIndexes()
        self._coordinateRings = [[]] * 3
        self._aroundCounts = [0] * 3
        self._connectionCounts = [0] * 3
        self._aCrossIndexes = None
        self._bCrossIndexes = None
        self._outerTubeBifurcationData = outerTubeBifurcationData

        if outerTubeBifurcationData:
            # inherit trim surfaces
            self._intersectionCurves = [(None, None, None, False)] * segmentCount  # not calculated
            self._trimSurfaces = outerTubeBifurcationData._trimSurfaces
        else:
            # get intersection curves between pairs of segments
            self._intersectionCurves = []
            for s in range(segmentCount):
                tubeData1 = self._tubeData[s]
                tubeTrackSurface1 = tubeData1.getRawTrackSurface()
                tubeData2 = self._tubeData[(s + 1) % segmentCount]
                tubeTrackSurface2 = tubeData2.getRawTrackSurface()
                cx, cd1, cProportions, loop = tubeTrackSurface1.findIntersectionCurve(tubeTrackSurface2)
                self._intersectionCurves.append((cx, cd1, cProportions, loop))

            # get trim surfaces
            self._trimSurfaces = []
            for s in range(segmentCount):
                networkSegment = self._networkSegments[s]
                tubeData = segmentTubeData[networkSegment]
                pathParameters = tubeData.getPathParameters()
                path_d1 = pathParameters[1]
                along = path_d1[-1 if networkSegment in self._networkSegmentsIn else 0]
                ax, ad1, _, aloop = self._intersectionCurves[s - 1]
                bx, bd1, _, bloop = self._intersectionCurves[s]
                trimSurface = None
                if ax and (len(ax) == 9) and not aloop and bx and (len(bx) == 9) and not bloop:
                    d2am = interpolateCubicHermiteSecondDerivative(ax[3], ad1[3], ax[4], ad1[4], 1.0)
                    d2ap = interpolateCubicHermiteSecondDerivative(ax[4], ad1[4], ax[5], ad1[5], 0.0)
                    d2a = add(d2am, d2ap)
                    cax = ax[4]
                    d2bm = interpolateCubicHermiteSecondDerivative(bx[3], bd1[3], bx[4], bd1[4], 1.0)
                    d2bp = interpolateCubicHermiteSecondDerivative(bx[4], bd1[4], bx[5], bd1[5], 0.0)
                    d2b = [-s for s in add(d2bm, d2bp)]
                    cbx = bx[4]
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
                    cax = sub(cax, mult(d2a, 0.05))
                    d2b = mult(d2b, arcLength / magnitude(d2b))
                    cbx = add(cbx, mult(d2b, 0.05))
                    nx = []
                    nd1 = []
                    nd2 = []
                    eCount = 2
                    for i in range(eCount + 1):
                        delta = (i - eCount // 2) * 1.2 * sizeAcross
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
                        nx.append(vax)
                        nd1.append(vad1)
                        nd2.append(mult(cad2, sizeAcross))
                        nx.append(vbx)
                        nd1.append(vbd1)
                        nd2.append(mult(cbd2, sizeAcross))
                    trimSurface = TrackSurface(1, eCount, nx, nd1, nd2)
                self._trimSurfaces.append(trimSurface)


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

    def getTubeCoordinates(self):
        return self._coordinateRings

    def getTubeData(self):
        return self._tubeData

    def determineCrossIndexes(self):
        assert self._aCrossIndexes is None  # should only call once

        for s in range(3):
            tubeData = self._tubeData[s]
            row = -2 if self._segmentsIn[s] else 1
            sampledTubeCoordinates = tubeData.getSampledTubeCoordinates()
            coordinateRing = [value[row] for value in sampledTubeCoordinates]
            self._coordinateRings[s] = coordinateRing
            self._aroundCounts[s] = len(coordinateRing[0])
        self._connectionCounts = get_tube_bifurcation_connection_elements_counts(self._aroundCounts)

        if self._outerTubeBifurcationData:
            self._aCrossIndexes = self._outerTubeBifurcationData.getCrossIndexes()
            assert self._aCrossIndexes, "Outer bifurcation data must be generated first"
        else:
            min_dist = None
            ibCrossIndex = 2 if self._segmentsIn[0] else 1
            jbCrossIndex = 0 if self._segmentsIn[1] else 2
            kbCrossIndex = 1 if self._segmentsIn[2] else 0
            irx = self._coordinateRings[0][0]
            jrx = self._coordinateRings[1][0]
            krx = self._coordinateRings[2][0]
            for ia in range(self._aroundCounts[0]):
                ib = ia - self._connectionCounts[ibCrossIndex]
                for ja in range(self._aroundCounts[1]):
                    jb = ja - self._connectionCounts[jbCrossIndex]
                    a12 = magnitude(sub(jrx[ja], irx[ia]))
                    b12 = magnitude(sub(jrx[jb], irx[ib]))
                    for ka in range(self._aroundCounts[2]):
                        kb = ka - self._connectionCounts[kbCrossIndex]
                        a23 = magnitude(sub(krx[ka], jrx[ja]))
                        b23 = magnitude(sub(krx[kb], jrx[jb]))
                        a31 = magnitude(sub(irx[ia], krx[ka]))
                        b31 = magnitude(sub(irx[ib], krx[kb]))
                        dist_a = a12 + a23 + a31
                        dist_b = b12 + b23 + b31
                        dist = dist_a + dist_b
                        if (min_dist is None) or (dist < min_dist):
                            self._aCrossIndexes = [ia, ja, ka]
                            self._bCrossIndexes = [ib, jb, kb]
                            # switch if wrong sense relative to b_to_a axis
                            aNormal = cross(sub(jrx[ja], irx[ia]), sub(krx[ka], jrx[ja]))
                            aMean = add(add(irx[ia], jrx[ja]), krx[ka])
                            bMean = add(add(irx[ib], jrx[jb]), krx[kb])
                            b_to_a = sub(aMean, bMean)
                            if dot(aNormal, b_to_a) < 0.0:
                                self._aCrossIndexes = [((self._bCrossIndexes[s] + self._aroundCounts[s])
                                                        if (self._bCrossIndexes[s] < 0) else self._bCrossIndexes[s])
                                                       for s in range(3)]
                                self._bCrossIndexes = [
                                    self._aCrossIndexes[0] - self._connectionCounts[ibCrossIndex],
                                    self._aCrossIndexes[1] - self._connectionCounts[jbCrossIndex],
                                    self._aCrossIndexes[2] - self._connectionCounts[kbCrossIndex]]
                            min_dist = dist

            # incomplete version trying to get shortest distance between directly connected points on adjacent tubes
            # minDist = None
            # aCrossIndexes = [None, None, None]
            # for i in range(self._aroundCounts[0]):
            #     aCrossIndexes[0] = i
            #     for j in range(self._aroundCounts[1]):
            #         aCrossIndexes[1] = j
            #         for k in range(self._aroundCounts[2]):
            #             aCrossIndexes[2] = k
            #             dist = 0.0
            #             for s in range(3):
            #                 ring1 = self._coordinateRings[s]
            #                 ring2 = self._coordinateRings[s - 2]
            #                 ic1 = aCrossIndexes[s]
            #                 ic2 = aCrossIndexes[s - 2]
            #                 for n in range(1, self._connectionCounts[s]):
            #                     i1 = (ic1 + n) % self._aroundCounts[s] if self._segmentsIn[s] else ic1 - n
            #                     i2 = (ic2 - n) if self._segmentsIn[s - 2] else (ic2 + n) % self._aroundCounts[s - 2]
            #                     delta = sub(ring1[0][i1], ring2[0][i2])
            #                     dist += magnitude(delta)
            #             if (minDist is None) or (dist < minDist):
            #                 self._aCrossIndexes = copy.copy(aCrossIndexes)
            #                 self._bCrossIndexes = [ib, jb, kb]

            # print("aCrossIndexes", self._aCrossIndexes, "bCrossIndexes", self._bCrossIndexes)

    def determineMidCoordinates(self):
        """
        Get 3 half-rings of coordinates between the tube coordinate rings.
        The order is between the 1-2, 2-3 and 3-1 tubes.
        Each half ring goes from the first cross point to the second.
        Coordinates at the cross points are averaged.
        """
        assert self._aCrossIndexes is not None  # must call determineCrossIndexes() first
        # get cross points a and b
        at1d2 = self._coordinateRings[0][2][self._aCrossIndexes[0]]
        at2d2 = self._coordinateRings[1][2][self._aCrossIndexes[1]]
        at3d2 = self._coordinateRings[2][2][self._aCrossIndexes[2]]
        acx, acd1, acd2 = get_bifurcation_triple_point(
            self._coordinateRings[0][0][self._aCrossIndexes[0]], [-d for d in at1d2] if self._segmentsIn[0] else at1d2,
            self._coordinateRings[1][0][self._aCrossIndexes[1]], [-d for d in at2d2] if self._segmentsIn[1] else at2d2,
            self._coordinateRings[2][0][self._aCrossIndexes[2]], [-d for d in at3d2] if self._segmentsIn[2] else at3d2)
        bt1d2 = self._coordinateRings[0][2][self._bCrossIndexes[0]]
        bt2d2 = self._coordinateRings[1][2][self._bCrossIndexes[1]]
        bt3d2 = self._coordinateRings[2][2][self._bCrossIndexes[2]]
        bcx, bcd1, bcd2 = get_bifurcation_triple_point(
            self._coordinateRings[0][0][self._bCrossIndexes[0]], [-d for d in bt1d2] if self._segmentsIn[0] else bt1d2,
            self._coordinateRings[2][0][self._bCrossIndexes[2]], [-d for d in bt3d2] if self._segmentsIn[2] else bt3d2,
            self._coordinateRings[1][0][self._bCrossIndexes[1]], [-d for d in bt2d2] if self._segmentsIn[1] else bt2d2)

        self._midCoordinates = []
        for s in range(3):
            self._midCoordinates.append([[acx], [acd1], [acd2]])
            ring1 = self._coordinateRings[s]
            ring2 = self._coordinateRings[s - 2]
            for n in range(1, self._connectionCounts[s]):
                i1 = (self._aCrossIndexes[s] + n) % self._aroundCounts[s] if self._segmentsIn[s] \
                    else self._aCrossIndexes[s] - n
                i2 = self._aCrossIndexes[s - 2] - n if self._segmentsIn[s - 2] \
                    else (self._aCrossIndexes[s - 2] + n) % self._aroundCounts[s - 2]
                # print("s", s, "n", n, "/", self._connectionCounts[s], "i1", i1, "i2", i2)
                t1x = ring1[0][i1]
                r1d2 = ring1[2][i1] if self._segmentsIn[s] else [-d for d in ring1[2][i1]]
                t2x = ring2[0][i2]
                r2d2 = [-d for d in ring2[2][i2]] if self._segmentsIn[s - 2] else ring2[2][i2]
                # 2.0 (original) looks better
                scaling = 2.0  # computeCubicHermiteDerivativeScaling(t1x, r1d2, t2x, r2d2)
                t1d2 = mult(r1d2, scaling)
                t2d2 = mult(r2d2, scaling)
                hx = interpolateCubicHermite(t1x, t1d2, t2x, t2d2, 0.5)
                hd1 = [0.0, 0.0, 0.0]
                hd2 = [0.5 * d for d in interpolateCubicHermiteDerivative(t1x, t1d2, t2x, t2d2, 0.5)]
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
                                elementsCountAround: int, targetElementAspectRatio: float,
                                elementsCountThroughWall: int, layoutAnnotationGroups: list=[], serendipity=False):
    """
    Generate a 2D, or 3D (thick walled) tube bifurcation tree mesh.
    :param networkMesh: Specification of network path and lateral sizes.
    :param region: Zinc region to generate mesh in.
    :param coordinates: Finite element coordinate field to define.
    :param nodeIdentifier: First node identifier to use.
    :param elementIdentifier: First 2D element identifier to use.
    :param elementsCountAround: Number of elements around tube.
    :param targetElementAspectRatio: Target ratio of element size along over around. Approximately satisfied.
    :param elementsCountThroughWall: Number of elements through wall if inner coordinates provided to make 3D elements,
    otherwise 1.
    :param layoutAnnotationGroups: Optional list of annotations defined on layout mesh for networkMesh.
    :param serendipity: True to use Hermite serendipity basis, False for regular Hermite with zero cross derivatives.
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
    # map from NetworkSegment to SegmentTubeData
    outerSegmentTubeData = {}
    innerSegmentTubeData = {} if layoutInnerCoordinates else None
    for networkSegment in networkSegments:
        pathParameters = get_nodeset_path_ordered_field_parameters(
            layoutNodes, layoutCoordinates, valueLabels,
            networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())
        outerSegmentTubeData[networkSegment] = tubeData = SegmentTubeData(pathParameters)
        px, pd1, pd2, pd12 = getPathRawTubeCoordinates(pathParameters, elementsCountAround)
        tubeData.setRawTubeCoordinates((px, pd1, pd2, pd12))
        if layoutInnerCoordinates:
            innerPathParameters = get_nodeset_path_ordered_field_parameters(
                layoutNodes, layoutInnerCoordinates, valueLabels,
                networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())
            innerSegmentTubeData[networkSegment] = innerTubeData = SegmentTubeData(innerPathParameters)
            px, pd1, pd2, pd12 = getPathRawTubeCoordinates(innerPathParameters, elementsCountAround)
            innerTubeData.setRawTubeCoordinates((px, pd1, pd2, pd12))
        for layoutAnnotationMeshGroup, annotationMeshGroup in layoutAnnotationMeshGroupMap:
            if networkSegment.hasLayoutElementsInMeshGroup(layoutAnnotationMeshGroup):
                tubeData.addAnnotationMeshGroup(annotationMeshGroup)
                if layoutInnerCoordinates:
                    innerTubeData.addAnnotationMeshGroup(annotationMeshGroup)

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
                if ((len(startInSegments) + len(startOutSegments)) == 3):
                    # print("create start", networkSegment, startSegmentNode)
                    startTubeBifurcationData = TubeBifurcationData(startInSegments, startOutSegments, segmentTubeData)
                    nodeTubeBifurcationData[startSegmentNode] = startTubeBifurcationData
                    newBifurcationData.append(startTubeBifurcationData)
            if startTubeBifurcationData:
                startSurface = startTubeBifurcationData.getSegmentTrimSurface(networkSegment)
            endSegmentNode = segmentNodes[-1]
            endTubeBifurcationData = nodeTubeBifurcationData.get(endSegmentNode)
            endSurface = None
            createEndBifurcationData = not endTubeBifurcationData
            if createEndBifurcationData:
                endInSegments = endSegmentNode.getInSegments()
                endOutSegments = endSegmentNode.getOutSegments()
                if ((len(endInSegments) + len(endOutSegments)) == 3):
                    # print("create end", networkSegment, endSegmentNode)
                    endTubeBifurcationData = TubeBifurcationData(endInSegments, endOutSegments, segmentTubeData)
                    nodeTubeBifurcationData[endSegmentNode] = endTubeBifurcationData
                    newBifurcationData.append(endTubeBifurcationData)
            if endTubeBifurcationData:
                endSurface = endTubeBifurcationData.getSegmentTrimSurface(networkSegment)
            if segmentTubeData is outerSegmentTubeData:
                pathParameters = tubeData.getPathParameters()
                segmentLength = getCubicHermiteCurvesLength(pathParameters[0], pathParameters[1])
                ringCount = len(rawTubeCoordinates[0])
                sumRingLength = 0.0
                for n in range(ringCount):
                    ringLength = getCubicHermiteCurvesLength(rawTubeCoordinates[0][n], rawTubeCoordinates[1][n], loop=True)
                    sumRingLength += ringLength
                meanElementLengthAround = sumRingLength / (ringCount * elementsCountAround)
                targetElementLength = targetElementAspectRatio * meanElementLengthAround
                elementsCountAlong = max(1, math.ceil(segmentLength / targetElementLength))
                loop = (len(startSegmentNode.getInSegments()) == 1) and \
                       (startSegmentNode.getInSegments()[0] is networkSegment) and \
                       (networkSegment.getNodeVersions()[0] == networkSegment.getNodeVersions()[-1])
                if (elementsCountAlong == 1) and (startTubeBifurcationData or endTubeBifurcationData):
                    # at least 2 segments if bifurcating at either end, or loop
                    elementsCountAlong = 2
                elif (elementsCountAlong < 3) and loop:
                    # at least 3 segments around loop; 2 should work, but zinc currently makes incorrect faces
                    elementsCountAlong = 3
            else:
                # must match count from outer surface!
                outerTubeData = outerSegmentTubeData[networkSegment]
                elementsCountAlong = outerTubeData.getSampledElementsCountAlong()
            # print("startSurface", startSurface is not None, "endSurface", endSurface is not None)
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
                        # for s in range(3):
                        #     curve = outerTubeBifurcationData.getIntersectionCurve(s)
                        #     cx, cd1, cProportions, loop = curve
                        #     if cx:
                        #         nodeIdentifier, elementIdentifier = \
                        #             generateCurveMesh(region, cx, cd1, loop, nodeIdentifier, elementIdentifier)
                        # for s in range(3):
                        #     trimSurface = outerTubeBifurcationData.getTrimSurface(s)
                        #     if trimSurface:
                        #         nodeIdentifier = trimSurface.generateMesh(region, nodeIdentifier)[0]
                        innerTubeBifurcationData = None
                        if innerNodeTubeBifurcationData:
                            innerTubeBifurcationData = innerNodeTubeBifurcationData.get(
                                startSegmentNode if (stage == 0) else endSegmentNode)

                        crossIndexes = outerTubeBifurcationData.getCrossIndexes()  # only get these from outer
                        if not crossIndexes:
                            outerTubeBifurcationData.determineCrossIndexes()
                            outerTubeBifurcationData.determineMidCoordinates()
                            if innerTubeBifurcationData:
                                innerTubeBifurcationData.determineCrossIndexes()
                                innerTubeBifurcationData.determineMidCoordinates()
                            crossIndexes = outerTubeBifurcationData.getCrossIndexes()

                        outerTubeCoordinates = outerTubeBifurcationData.getTubeCoordinates()
                        outerMidCoordinates = outerTubeBifurcationData.getMidCoordinates()
                        inward = outerTubeBifurcationData.getSegmentsIn()
                        outerTubeData = outerTubeBifurcationData.getTubeData()
                        tubeNodeIds = [outerTubeData[s].getEndNodeIds(1) if inward[s] else \
                                           outerTubeData[s].getStartNodeIds(1) for s in range(3)]
                        innerTubeCoordinates = None
                        innerMidCoordinates = None
                        if innerTubeBifurcationData:
                            innerTubeCoordinates = innerTubeBifurcationData.getTubeCoordinates()
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
