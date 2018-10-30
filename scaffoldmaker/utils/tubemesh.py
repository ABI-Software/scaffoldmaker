'''
Utility functions for generating tubular mesh.
Created on Oct 9, 2018

@author: Mabelle Lin
'''
from __future__ import division
import math
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.vector import *
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

def generatetubemesh(region,
    elementsCountAlong,
    elementsCountAround,
    elementsCountThroughWall,
    cx, cd1, cd2, cd3, t2, t2d, t3, t3d,
    useCrossDerivatives,
    useCubicHermiteThroughWall, # or Zinc Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE etc.
    nextNodeIdentifier = 1, nextElementIdentifier = 1
    ):
    '''
    Generates a 3-D tubular mesh with variable numbers of elements
    around, along the central axis, and radially through wall.
    The tubular mesh follows the profile of a central line and 
    has a elliptical cross-section.
    :param elementsCountAlong: number of elements along tube
    :param elementsCountAround: number of elements around tube
    :param elementsCountThroughWall: number of elements through wall thickness
    :param cx: coordinates on central line
    :param cd1: derivative 1 along central line
    :param cd2: derivative 2 to the side
    :param cd3: derivative 3 to the side, perpendicular to direction of cd2
    :param t2: wall thickness in cd2 direction
    :param t2d: wall thickness derivative in cd2 direction (rate of change)
    :param t3: wall thickness in cd3 direction
    :param t3d: wall thickness derivative in cd3 direction (rate of change)
    :param useCubicHermiteThroughWall: use linear when false
    :param nextNodeIdentifier, nextElementIdentifier: Next identifiers to use and increment.
    :return: Final values of nextNodeIdentifier, nextElementIdentifier
    '''

    zero = [ 0.0, 0.0, 0.0 ]

    # Sample central line to get same number of elements as elementsCountAlong
    sx, sd1, se, sxi, _ = sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
    sd2 = interpolateSampleLinear(cd2, se, sxi)
    sd3 = interpolateSampleLinear(cd3, se, sxi)
    st2 = interpolateSampleLinear(t2, se, sxi)
    st2d = interpolateSampleLinear(t2d, se, sxi)
    st3 = interpolateSampleLinear(t3, se, sxi)
    st3d = interpolateSampleLinear(t3d, se, sxi)

    # Find unit normals and binormals at each sample points
    sNormal = []
    sBinormal = []

    # Set up normal and binormal for first frame
    firstUnitTangent = normalise(sd1[0])
    if magnitude(crossproduct3(firstUnitTangent,[0.0, 0.0, 1.0])) > 0.0:
        firstBinormal = crossproduct3(firstUnitTangent,[0.0, 0.0, 1.0])
    else:
        firstBinormal = crossproduct3(firstUnitTangent,[0.0, -1.0, 0.0])
    firstUnitBinormal = [c/magnitude(firstBinormal) for c in firstBinormal]
    firstUnitNormal = crossproduct3(firstUnitBinormal, firstUnitTangent)
    sNormal.append(firstUnitNormal)
    sBinormal.append(firstUnitBinormal)

    # Step through central line and rotate axes to align tangent to tangent from first frame
    for n in range(1, elementsCountAlong+1):
        unitTangent = normalise(sd1[n])
        cp = crossproduct3(firstUnitTangent, unitTangent)
        if magnitude(cp)> 0.0:
            axisRot = [c/magnitude(cp) for c in cp]
            thetaRot = math.acos(dotproduct(firstUnitTangent, unitTangent))
            rotFrame = rotationMatrixAboutAxis(axisRot, thetaRot)
            rotNormal = [rotFrame[j][0]*firstUnitNormal[0] + rotFrame[j][1]*firstUnitNormal[1] + rotFrame[j][2]*firstUnitNormal[2] for j in range(3)]
            unitNormal = [c/magnitude(rotNormal) for c in rotNormal]
            unitBinormal = crossproduct3(unitTangent, unitNormal)
        else:
            unitBinormal = firstUnitBinormal
            unitNormal = firstUnitNormal
        sNormal.append(unitNormal)
        sBinormal.append(unitBinormal)

    fm = region.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()
    coordinates = getOrCreateCoordinateField(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
    if useCubicHermiteThroughWall:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

    mesh = fm.findMeshByDimension(3)

    if useCubicHermiteThroughWall:
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    else:
        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    result = elementtemplate.defineField(coordinates, -1, eft)

    # create nodes
    nodeIdentifier = nextNodeIdentifier
    x = [ 0.0, 0.0, 0.0 ]
    dx_ds1 = [ 0.0, 0.0, 0.0 ]
    dx_ds2 = [ 0.0, 0.0, 0.0 ]
    dx_ds3 = [ 0.0, 0.0, 0.0 ]

    for n3 in range(elementsCountThroughWall + 1):

        for n2 in range(elementsCountAlong+1):

            aThroughWallElement = sd2[n2][1] + st2[n2]*(n3/elementsCountThroughWall) # Check sd2, sd3
            bThroughWallElement = sd3[n2][2] + st3[n2]*(n3/elementsCountThroughWall)
            perimeterAroundWallElement = getApproximateEllipsePerimeter(aThroughWallElement, bThroughWallElement)
            arcLengthPerElementAround = perimeterAroundWallElement / elementsCountAround
            prevRadiansAround = updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement, 0.0, arcLengthPerElementAround)
            st2PerWallElement = st2[n2]/elementsCountThroughWall
            st3PerWallElement = st3[n2]/elementsCountThroughWall

            # Pre-calculate next node downstream for arclength calculation
            if n2 < elementsCountAlong:
                aThroughWallElementNext = sd2[n2][1] + st2[n2]*(n3/elementsCountThroughWall) # Check sd2, sd3 
                bThroughWallElementNext = sd3[n2][2] + st3[n2]*(n3/elementsCountThroughWall)
                perimeterAroundWallElementNext = getApproximateEllipsePerimeter(aThroughWallElementNext, bThroughWallElementNext)
                arcLengthPerElementAroundNext = perimeterAroundWallElementNext / elementsCountAround

            for n1 in range(elementsCountAround):
                arcLengthAround = n1*arcLengthPerElementAround
                radiansAround = -1*updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement, 0.0, arcLengthAround)
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)

                x = [sx[n2][j] + aThroughWallElement*cosRadiansAround*sBinormal[n2][j] + bThroughWallElement*sinRadiansAround*sNormal[n2][j] for j in range(3)]
                dx_ds1 = [(radiansAround - prevRadiansAround)*(aThroughWallElement*-sinRadiansAround*sBinormal[n2][j] + bThroughWallElement*cosRadiansAround*sNormal[n2][j]) for j in range(3)]

                if n2 < elementsCountAlong:
                    arcLengthAroundNext = n1*arcLengthPerElementAroundNext
                    radiansAroundNext = -1*updateEllipseAngleByArcLength(aThroughWallElementNext, bThroughWallElementNext, 0.0, arcLengthAroundNext)
                    cosRadiansAroundNext = math.cos(radiansAroundNext)
                    sinRadiansAroundNext = math.sin(radiansAroundNext)
                    xNext = [sx[n2+1][j] + aThroughWallElementNext*cosRadiansAroundNext*sBinormal[n2+1][j] + bThroughWallElementNext*sinRadiansAroundNext*sNormal[n2+1][j] for j in range(3)]
                    cubicArcLength = computeCubicHermiteArcLength(x, sd1[n2], xNext, sd1[n2+1], True)
                dx_ds2 = [cubicArcLength*sd1[n2][j]/magnitude(sd1[n2]) for j in range(3)]

                dx_ds3 = [ st2PerWallElement*cosRadiansAround*sBinormal[n2][j] + st3PerWallElement*sinRadiansAround*sNormal[n2][j] for j in range(3)] # Double check if correct

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if useCubicHermiteThroughWall:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    if useCubicHermiteThroughWall:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                nodeIdentifier = nodeIdentifier + 1
                prevRadiansAround = radiansAround

    # # For debugging - Nodes along central line
    # for pt in range(len(sx)):
        # node = nodes.createNode(nodeIdentifier, nodetemplate)
        # cache.setNode(node)
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sx[pt])
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1[pt])
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sNormal[pt])
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, sBinormal[pt])
        # nodeIdentifier = nodeIdentifier + 1

    # create elements
    elementIdentifier = nextElementIdentifier
    now = (elementsCountAlong + 1)*elementsCountAround
    for e3 in range(elementsCountThroughWall):
        for e2 in range(elementsCountAlong):
            for e1 in range(elementsCountAround):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni11 = e3*now + e2*elementsCountAround + e1 + 1
                bni12 = e3*now + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
                bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

    fm.endChange()

    return nodeIdentifier, elementIdentifier

def rotationMatrixAboutAxis(rotAxis, theta):
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)
    C = 1 - cosTheta
    rotMatrix = ([[rotAxis[0]*rotAxis[0]*C + cosTheta, rotAxis[0]*rotAxis[1]*C - rotAxis[2]*sinTheta, rotAxis[0]*rotAxis[2]*C + rotAxis[1]*sinTheta],
        [rotAxis[1]*rotAxis[0]*C + rotAxis[2]*sinTheta, rotAxis[1]*rotAxis[1]*C + cosTheta, rotAxis[1]*rotAxis[2]*C - rotAxis[0]*sinTheta],
        [rotAxis[2]*rotAxis[0]*C - rotAxis[1]*sinTheta, rotAxis[2]*rotAxis[1]*C + rotAxis[0]*sinTheta, rotAxis[2]*rotAxis[2]*C + cosTheta]])
    return rotMatrix
