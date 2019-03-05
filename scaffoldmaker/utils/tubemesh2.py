'''
Utility function for generating tubular mesh from a central line
using a unit profile.
Created on 18 Feb, 2019

@author: Mabelle Lin
'''
from __future__ import division
import math
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.matrix import *
from scaffoldmaker.utils.vector import *
from scaffoldmaker.utils.zinc_utils import *
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

def generatetubemesh(region,
    elementsCountAround,
    elementsCountAlongUnit,
    elementsCountThroughWall,
    unitsCountAlong,
    cx, cd1,
    xOuter, d1Outer, d2Outer, wallThickness,
    unitAxis, unitLength,
    useCrossDerivatives,
    useCubicHermiteThroughWall, # or Zinc Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE etc.
    nextNodeIdentifier = 1, nextElementIdentifier = 1
    ):
    '''
    Generates a 3-D tubular mesh with variable numbers of elements
    around, along the central axis, and radially through wall. The
    tubular mesh is created from a unit profile which is mapped onto
    the central line and lateral axes data
    :param elementsCountAround: number of elements around tube
    :param elementsCountAlongUnit: number of elements along unit profile
    :param elementsCountThroughWall: number of elements through wall thickness
    :param unitsCountAlong: number of units along the tube
    :param cx: coordinates on central line
    :param cd1: derivative along central line
    :param xOuter: coordinates on outer surface of unit profile
    :param d1Outer: derivatives around outer surface of unit profile
    :param d2Outer: derivatives along outer surface of unit profile
    :param wallThickness: thickness of wall
    :param unitAxis: axis of unit profile
    :param unitLength: length of unit profile
    :param useCubicHermiteThroughWall: use linear when false
    '''
    zero  = [0.0, 0.0, 0.0]
    elementsCountAlong = elementsCountAlongUnit*unitsCountAlong

    # Sample central line to get same number of elements as elementsCountAlong
    sx, sd1, se, sxi, _ = sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)

    # Find unit normals and binormals at each sample points
    sNormal = []
    sBinormal = []

    # Set up normal and binormal for first frame
    prevUnitTangent = normalise(sd1[0])
    if magnitude(crossproduct3(prevUnitTangent,[0.0, 0.0, 1.0])) > 0.0:
        prevBinormal = crossproduct3(prevUnitTangent,[0.0, 0.0, 1.0])
    else:
        prevBinormal = crossproduct3(prevUnitTangent,[0.0, -1.0, 0.0])
    prevUnitBinormal = normalise(prevBinormal)
    prevUnitNormal = crossproduct3(prevUnitBinormal, prevUnitTangent)
    sNormal.append(prevUnitNormal)
    sBinormal.append(prevUnitBinormal)

    # Step through central line and rotate central line axes to align tangent 
    # to tangent from previous frame
    for n in range(1, elementsCountAlong+1):
        unitTangent = normalise(sd1[n])
        cp = crossproduct3(prevUnitTangent, unitTangent)
        if magnitude(cp)> 0.0:
            axisRot = normalise(cp)
            thetaRot = math.acos(dotproduct(prevUnitTangent, unitTangent))
            rotFrame = getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            rotNormal = [rotFrame[j][0]*prevUnitNormal[0] + rotFrame[j][1]*prevUnitNormal[1] + rotFrame[j][2]*prevUnitNormal[2] for j in range(3)]
            unitNormal = normalise(rotNormal)
            unitBinormal = crossproduct3(unitTangent, unitNormal)
            prevUnitTangent = unitTangent
            prevUnitNormal = unitNormal
        else:
            unitBinormal = prevUnitBinormal
            unitNormal = prevUnitNormal
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
    xOuterList = []
    d1OuterList = []
    d2OuterList = []
    d3OuterUnitList = []
    xList = []
    dx_ds1List = []
    dx_ds2List = []
    dx_ds3List = []

    # Map each face along unit profile to central line
    for nUnit in range(unitsCountAlong):
        for nAlongUnit in range(elementsCountAlongUnit+1):
            n2 = nUnit*elementsCountAlongUnit + nAlongUnit
            if nUnit == 0 or (nUnit > 0 and nAlongUnit > 0):
                # Rotate to align unit axis with tangent of central line
                unitMid = [0.0, 0.0, unitLength/elementsCountAlongUnit* nAlongUnit]
                unitTangent = normalise(sd1[n2])
                cp = crossproduct3(unitAxis, unitTangent)
                if magnitude(cp)> 0.0:
                    axisRot = normalise(cp)
                    thetaRot = math.acos(dotproduct(unitAxis, unitTangent))
                    rotFrame = getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                    midRot = [rotFrame[j][0]*unitMid[0] + rotFrame[j][1]*unitMid[1] + rotFrame[j][2]*unitMid[2] for j in range(3)]
                    translateMatrix = [sx[n2][j] - midRot[j] for j in range(3)]
                else:
                    midRot = unitMid

                for n1 in range(elementsCountAround):
                    n = nAlongUnit*elementsCountAround + n1
                    x = xOuter[n]
                    d1 = d1Outer[n]
                    d2 = d2Outer[n]
                    if magnitude(cp)> 0.0:
                        xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                        d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                        d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]
                        # Rotate to align first vector on face with binormal axis
                        if n1 == 0:
                            firstVector = normalise([xRot1[j] - midRot[j] for j in range(3)])
                            thetaRot2 = math.acos(dotproduct(normalise(sBinormal[n2]), firstVector))
                            cp2 = normalise(crossproduct3(normalise(sBinormal[n2]), firstVector))
                            signThetaRot2 = dotproduct(unitTangent, cp2)
                            axisRot2 = unitTangent
                            rotFrame2 = getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                        xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
                        d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
                        d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
                    else:
                        xRot2 = x
                        d1Rot2 = d1
                        d2Rot2 = d2

                    xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

                    xOuterList.append(xTranslate)
                    d1OuterList.append(d1Rot2)
                    d2OuterList.append(d2Rot2)
                    d3Unit = normalise(crossproduct3(normalise(d2Rot2), normalise(d1Rot2)))
                    d3OuterUnitList.append(d3Unit)

    # Pre-calculate node locations and derivatives on outer boundary
    xInnerList, d1Inner, curvatureOuter = getInnerCoordinatesAndDerivativesFromOuter(xOuterList, d1OuterList, d3OuterUnitList, wallThickness, elementsCountAlong, elementsCountAround)

    # Interpolate to get nodes through wall
    for n3 in range(elementsCountThroughWall + 1):
        xi3 = 1/elementsCountThroughWall * n3
        x, dx_ds1, dx_ds2, dx_ds3 = interpolatefromInnerAndOuter( xInnerList, xOuterList, wallThickness, xi3, curvatureOuter, d1OuterList, d2OuterList, d3OuterUnitList, elementsCountAround, elementsCountAlong, elementsCountThroughWall)
        xList = xList + x
        dx_ds1List = dx_ds1List + dx_ds1
        dx_ds2List = dx_ds2List + dx_ds2
        dx_ds3List = dx_ds3List + dx_ds3

    for n in range(len(xList)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1List[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2List[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3List[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        # print('NodeIdentifier = ', nodeIdentifier, xList[n])
        nodeIdentifier = nodeIdentifier + 1

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

def getInnerCoordinatesAndDerivativesFromOuter(xOuter, d1Outer, d3Outer, wallThickness, elementsCountAlong, elementsCountAround):
    """
    Generates coordinates and derivatives of inner surface from 
    coordinates and derivatives of outer surface using wall thickness
    and normals.
    param xOuter: Coordinates on outer surface
    param d1Outer: Derivatives on outer surface around tube
    param d3Outer: Derivatives on outer surface through wall
    param wallThickness: Thickness of wall
    param elementsCountAlong: Number of elements along tube
    param elementsCountAround: Number of elements around tube
    return xInner: Coordinates on inner surface
    return nd1: Derivatives on inner surface around tube
    return curvatureOuter: Curvature of coordinates on outer surface
    """
    xInner = []
    nd1 = []
    dWall = []
    curvatureOuter = []
    for n2 in range(elementsCountAlong + 1):
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            x = [xOuter[n][i] + d3Outer[n][i]*wallThickness for i in range(3)]
            prevIdx = n-1 if (n1 != 0) else (n2+1)*elementsCountAround - 1
            nextIdx = n+1 if (n1 < elementsCountAround-1) else n2*elementsCountAround
            norm = [-1.0*c for c in d3Outer[n]]
            curvatureAround = 0.5*(
                getCubicHermiteCurvature(xOuter[prevIdx], d1Outer[prevIdx], xOuter[n], d1Outer[n], norm, 1.0) +
                getCubicHermiteCurvature(xOuter[n], d1Outer[n], xOuter[nextIdx], d1Outer[nextIdx], norm, 0.0))
            factor = 1.0 - curvatureAround*wallThickness
            nd1Inner = [ factor*c for c in d1Outer[n]]
            xInner.append(x)
            nd1.append(nd1Inner)
            curvatureOuter.append(curvatureAround)

    return xInner, nd1, curvatureOuter

def interpolatefromInnerAndOuter( xInner, xOuter, thickness, xi3, curvatureOuter, d1Outer, d2Outer, d3OuterUnit,
    elementsCountAround, elementsCountAlong, elementsCountThroughWall):
    """
    Generate coordinates and derivatives at xi3 by interpolating with 
    inner and outer coordinates and derivatives.
    param xInner: Coordinates on inner surface
    param xOuter: Coordinates on outer surface
    param thickness: Thickness of wall
    param curvatureOuter: Curvature of coordinates on outer surface
    param d1Outer: Derivatives on outer surface around tube
    param d2outer: Derivatives on outer surface along tube
    param d3OuterUnit: Unit derivatives on inner surface through wall
    param elementsCountAround: Number of elements around tube
    param elementsCountAlong: Number of elements along tube
    param elementsCountThroughWall: Number of elements through wall
    return xList, dx_ds1List, dx_ds2List, dx_ds3List: Coordinates and derivatives on xi3
    """
    xList = []
    dx_ds1List = []
    dx_ds2List = []
    dx_ds3List =[]

    for n2 in range(elementsCountAlong+1):
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = [-1* c for c in d3OuterUnit[n]]
            # x
            innerx = xInner[n]
            outerx = xOuter[n]
            dWall = [thickness*c for c in norm]
            x = interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
            xList.append(x)
            # dx_ds1
            factor = 1.0 - curvatureOuter[n]*thickness*xi3
            dx_ds1 = [ factor*c for c in d1Outer[n]]
            dx_ds1List.append(dx_ds1)
            # dx_ds2
            if n2 > 0 and n2 < elementsCountAlong:
                prevIdx = (n2-1)*elementsCountAround + n1
                nextIdx = (n2+1)*elementsCountAround + n1
                curvatureAround = 0.5*(
                    getCubicHermiteCurvature(xOuter[prevIdx], d2Outer[prevIdx], xOuter[n], d2Outer[n], norm, 1.0)+
                    getCubicHermiteCurvature(xOuter[n], d2Outer[n], xOuter[nextIdx], d2Outer[nextIdx], norm, 0.0))
            elif n2 == 0:
                nextIdx = (n2+1)*elementsCountAround + n1
                curvatureAround = getCubicHermiteCurvature(xOuter[n], d2Outer[n], xOuter[nextIdx], d2Outer[nextIdx], norm, 0.0)
            else:
                prevIdx = (n2-1)*elementsCountAround + n1
                curvatureAround = getCubicHermiteCurvature(xOuter[prevIdx], d2Outer[prevIdx], xOuter[n], d2Outer[n], norm, 1.0)

            factor = 1.0 - curvatureAround*thickness*xi3
            dx_ds2 = [ factor*c for c in d2Outer[n]]
            dx_ds2List.append(dx_ds2)
            #dx_ds3
            dx_ds3 = [c * thickness/elementsCountThroughWall for c in norm]
            dx_ds3List.append(dx_ds3)

    return xList, dx_ds1List, dx_ds2List, dx_ds3List
