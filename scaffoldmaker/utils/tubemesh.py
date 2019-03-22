'''
Utility function for generating tubular mesh from a central line
using a segment profile.
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
    elementsCountAlongSegment,
    elementsCountThroughWall,
    segmentCountAlong,
    cx, cd1,
    xInner, d1Inner, d2Inner, wallThickness,
    segmentAxis, segmentLength,
    useCrossDerivatives,
    useCubicHermiteThroughWall, # or Zinc Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE etc.
    firstNodeIdentifier = 1, firstElementIdentifier = 1
    ):
    '''
    Generates a 3-D tubular mesh with variable numbers of elements
    around, along the central axis, and radially through wall. The
    tubular mesh is created from a segment profile which is mapped onto
    the central line and lateral axes data
    :param elementsCountAround: number of elements around tube
    :param elementsCountAlongSegment: number of elements along segment profile
    :param elementsCountThroughWall: number of elements through wall thickness
    :param segmentCountAlong: number of segments along the tube
    :param cx: coordinates on central line
    :param cd1: derivative along central line
    :param xInner: coordinates on inner surface of segment profile
    :param d1Inner: derivatives around inner surface of segment profile
    :param d2Inner: derivatives along inner surface of segment profile
    :param wallThickness: thickness of wall
    :param segmentAxis: axis of segment profile
    :param segmentLength: length of segment profile
    :param useCubicHermiteThroughWall: use linear when false
    :return: annotationGroups, nodeIdentifier, elementIdentifier
    '''
    zero  = [0.0, 0.0, 0.0]
    annotationGroups = []
    elementsCountAlong = elementsCountAlongSegment*segmentCountAlong

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
    nodeIdentifier = firstNodeIdentifier
    x = [ 0.0, 0.0, 0.0 ]
    dx_ds1 = [ 0.0, 0.0, 0.0 ]
    dx_ds2 = [ 0.0, 0.0, 0.0 ]
    dx_ds3 = [ 0.0, 0.0, 0.0 ]
    xInnerList = []
    d1InnerList = []
    d2InnerList = []
    d3InnerUnitList = []
    xList = []
    dx_ds1List = []
    dx_ds2List = []
    dx_ds3List = []

    # Map each face along segment profile to central line
    for nSegment in range(segmentCountAlong):
        for nAlongSegment in range(elementsCountAlongSegment+1):
            n2 = nSegment*elementsCountAlongSegment + nAlongSegment
            if nSegment == 0 or (nSegment > 0 and nAlongSegment > 0):
                # Rotate to align segment axis with tangent of central line
                segmentMid = [0.0, 0.0, segmentLength/elementsCountAlongSegment* nAlongSegment]
                unitTangent = normalise(sd1[n2])
                cp = crossproduct3(segmentAxis, unitTangent)
                if magnitude(cp)> 0.0:
                    axisRot = normalise(cp)
                    thetaRot = math.acos(dotproduct(segmentAxis, unitTangent))
                    rotFrame = getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                    midRot = [rotFrame[j][0]*segmentMid[0] + rotFrame[j][1]*segmentMid[1] + rotFrame[j][2]*segmentMid[2] for j in range(3)]
                    translateMatrix = [sx[n2][j] - midRot[j] for j in range(3)]
                else:
                    midRot = segmentMid

                for n1 in range(elementsCountAround):
                    n = nAlongSegment*elementsCountAround + n1
                    x = xInner[n]
                    d1 = d1Inner[n]
                    d2 = d2Inner[n]
                    if magnitude(cp)> 0.0:
                        xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                        d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                        d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]
                        # Rotate to align first vector on face with binormal axis
                        if n1 == 0:
                            firstVector = normalise([xRot1[j] - midRot[j] for j in range(3)])
                            thetaRot2 = math.acos(dotproduct(normalise(sBinormal[n2]), firstVector))
                            cp2 = crossproduct3(normalise(sBinormal[n2]), firstVector)
                            if magnitude(cp2) > 0.0:
                                cp2 = normalise(cp2)
                                signThetaRot2 = dotproduct(unitTangent, cp2)
                                axisRot2 = unitTangent
                                rotFrame2 = getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                            else:
                                rotFrame2 = [ [1, 0, 0], [0, 1, 0], [0, 0, 1]]
                        xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
                        d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
                        d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
                    else:
                        xRot2 = x
                        d1Rot2 = d1
                        d2Rot2 = d2

                    xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

                    xInnerList.append(xTranslate)
                    d1InnerList.append(d1Rot2)
                    d2InnerList.append(d2Rot2)
                    d3Unit = normalise(crossproduct3(normalise(d1Rot2), normalise(d2Rot2)))
                    d3InnerUnitList.append(d3Unit)

    # Pre-calculate node locations and derivatives on outer boundary
    xOuterList, curvatureInner = getOuterCoordinatesAndCurvatureFromInner(xInnerList, d1InnerList, d3InnerUnitList, wallThickness, elementsCountAlong, elementsCountAround)

    # Interpolate to get nodes through wall
    for n3 in range(elementsCountThroughWall + 1):
        xi3 = 1/elementsCountThroughWall * n3
        x, dx_ds1, dx_ds2, dx_ds3 = interpolatefromInnerAndOuter( xInnerList, xOuterList, wallThickness, xi3, curvatureInner, d1InnerList, d2InnerList, d3InnerUnitList, elementsCountAround, elementsCountAlong, elementsCountThroughWall)
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
    elementIdentifier = firstElementIdentifier
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

    # Create texture coordinates
    textureCoordinates = getOrCreateTextureCoordinateField(fm)
    textureNodetemplate1 = nodes.createNodetemplate()
    textureNodetemplate1.defineField(textureCoordinates)
    textureNodetemplate1.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    textureNodetemplate1.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    textureNodetemplate1.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        textureNodetemplate1.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

    textureNodetemplate2 = nodes.createNodetemplate()
    textureNodetemplate2.defineField(textureCoordinates)
    textureNodetemplate2.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_VALUE, 2)
    textureNodetemplate2.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
    textureNodetemplate2.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D_DS2, 2)
    if useCrossDerivatives:
        textureNodetemplate2.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 2)

    bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eftTexture = bicubichermitelinear.createEftBasic()

    elementtemplate1 = mesh.createElementtemplate()
    elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate1.defineField(textureCoordinates, -1, eftTexture)

    eftTexture2 = bicubichermitelinear.createEftOpenTube()
    elementtemplate2 = mesh.createElementtemplate()
    elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate2.defineField(textureCoordinates, -1, eftTexture2)

    # Calculate texture coordinates and derivatives
    uTexture = []
    d1Texture = []
    d2Texture = []
    for n3 in range(elementsCountThroughWall + 1):
        for n2 in range(elementsCountAlong + 1):
            for n1 in range(elementsCountAround + 1):
                u = [ 1.0 / elementsCountAround * n1,
                    1.0 / elementsCountAlong * n2,
                    1.0 / elementsCountThroughWall * n3]
                d1 = [1.0 / elementsCountAround, 0.0, 0.0]
                d2 = [0.0, 1.0 / elementsCountAlong, 0.0]
                uTexture.append(u)
                d1Texture.append(d1)
                d2Texture.append(d2)

    nodeIdentifier = firstNodeIdentifier
    for n in range(len(uTexture)):
        if n%(elementsCountAround+1) == 0.0:
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            node.merge(textureNodetemplate2)
            cache.setNode(node)
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, uTexture[n])
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Texture[n])
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Texture[n])
            endIdx = n + elementsCountAround
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2, uTexture[endIdx])
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1Texture[endIdx])
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2, d2Texture[endIdx])
            if useCrossDerivatives:
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
            nodeIdentifier = nodeIdentifier + 1
        elif (n+1)%(elementsCountAround+1) > 0:
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            node.merge(textureNodetemplate1)
            cache.setNode(node)
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, uTexture[n])
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Texture[n])
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Texture[n])
            if useCrossDerivatives:
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

    # define texture coordinates field over elements
    elementIdentifier = firstElementIdentifier
    now = (elementsCountAlong + 1)*elementsCountAround
    for e3 in range(elementsCountThroughWall):
        for e2 in range(elementsCountAlong):
            for e1 in range(elementsCountAround):
                if e1 < elementsCountAround - 1:
                    element = mesh.findElementByIdentifier(elementIdentifier)
                    element.merge(elementtemplate1)
                    bni11 = e3*now + e2*elementsCountAround + e1 + 1
                    bni12 = e3*now + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                    bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
                    bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                    nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                    result = element.setNodesByIdentifier(eftTexture, nodeIdentifiers)
                else:
                    element = mesh.findElementByIdentifier(elementIdentifier)
                    element.merge(elementtemplate2)
                    # element = mesh.createElement(elementIdentifier, elementtemplate2)
                    bni11 = e3*now + e2*elementsCountAround + e1 + 1
                    bni12 = e3*now + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                    bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
                    bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                    nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now]
                    result2 = element.setNodesByIdentifier(eftTexture2, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

    fm.endChange()

    return annotationGroups, nodeIdentifier, elementIdentifier

def getOuterCoordinatesAndCurvatureFromInner(xInner, d1Inner, d3Inner, wallThickness, elementsCountAlong, elementsCountAround):
    """
    Generates coordinates on outer surface and curvature of inner
    surface from coordinates and derivatives of inner surface using
    wall thickness and normals.
    param xInner: Coordinates on inner surface
    param d1Inner: Derivatives on inner surface around tube
    param d3Inner: Derivatives on inner surface through wall
    param wallThickness: Thickness of wall
    param elementsCountAlong: Number of elements along tube
    param elementsCountAround: Number of elements around tube
    return xOuter: Coordinates on outer surface
    return curvatureInner: Curvature of coordinates on inner surface
    """
    xOuter = []
    curvatureInner = []
    for n2 in range(elementsCountAlong + 1):
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            x = [xInner[n][i] + d3Inner[n][i]*wallThickness for i in range(3)]
            prevIdx = n-1 if (n1 != 0) else (n2+1)*elementsCountAround - 1
            nextIdx = n+1 if (n1 < elementsCountAround-1) else n2*elementsCountAround
            norm = d3Inner[n]
            curvatureAround = 0.5*(
                getCubicHermiteCurvature(xInner[prevIdx], d1Inner[prevIdx], xInner[n], d1Inner[n], norm, 1.0) +
                getCubicHermiteCurvature(xInner[n], d1Inner[n], xInner[nextIdx], d1Inner[nextIdx], norm, 0.0))
            xOuter.append(x)
            curvatureInner.append(curvatureAround)

    return xOuter, curvatureInner

def interpolatefromInnerAndOuter( xInner, xOuter, thickness, xi3, curvatureInner, d1Inner, d2Inner, d3InnerUnit,
    elementsCountAround, elementsCountAlong, elementsCountThroughWall):
    """
    Generate coordinates and derivatives at xi3 by interpolating with 
    inner and outer coordinates and derivatives.
    param xInner: Coordinates on inner surface
    param xOuter: Coordinates on outer surface
    param thickness: Thickness of wall
    param curvatureInner: Curvature of coordinates on inner surface
    param d1Inner: Derivatives on inner surface around tube
    param d2Inner: Derivatives on inner surface along tube
    param d3InnerUnit: Unit derivatives on inner surface through wall
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
            norm = d3InnerUnit[n]
            # x
            innerx = xInner[n]
            outerx = xOuter[n]
            dWall = [thickness*c for c in norm]
            x = interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
            xList.append(x)
            # dx_ds1
            factor = 1.0 - curvatureInner[n]*thickness*xi3
            dx_ds1 = [ factor*c for c in d1Inner[n]]
            dx_ds1List.append(dx_ds1)
            # dx_ds2
            if n2 > 0 and n2 < elementsCountAlong:
                prevIdx = (n2-1)*elementsCountAround + n1
                nextIdx = (n2+1)*elementsCountAround + n1
                curvatureAround = 0.5*(
                    getCubicHermiteCurvature(xInner[prevIdx], d2Inner[prevIdx], xInner[n], d2Inner[n], norm, 1.0)+
                    getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[nextIdx], d2Inner[nextIdx], norm, 0.0))
            elif n2 == 0:
                nextIdx = (n2+1)*elementsCountAround + n1
                curvatureAround = getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[nextIdx], d2Inner[nextIdx], norm, 0.0)
            else:
                prevIdx = (n2-1)*elementsCountAround + n1
                curvatureAround = getCubicHermiteCurvature(xInner[prevIdx], d2Inner[prevIdx], xInner[n], d2Inner[n], norm, 1.0)

            factor = 1.0 - curvatureAround*thickness*xi3
            dx_ds2 = [ factor*c for c in d2Inner[n]]
            dx_ds2List.append(dx_ds2)
            #dx_ds3
            dx_ds3 = [c * thickness/elementsCountThroughWall for c in norm]
            dx_ds3List.append(dx_ds3)

    return xList, dx_ds1List, dx_ds2List, dx_ds3List
