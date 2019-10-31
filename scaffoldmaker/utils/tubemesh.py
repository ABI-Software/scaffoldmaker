'''
Utility function for generating tubular mesh from a central line
using a segment profile.
'''
from __future__ import division
import math
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector
from scaffoldmaker.utils import zinc_utils
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


def warpSegmentPoints(xList, d1List, d2List, segmentAxis, segmentLength,
    sx, sd1, sd2, elementsCountAround, elementsCountAlongSegment, nSegment):
    """
    Warps points in segment to account for bending and twisting
    along central path defined by nodes sx and derivatives sd1 and sd2.
    :param xList: coordinates of segment points.
    :param d1List: derivatives around axis of segment.
    :param d2List: derivatives along axis of segment.
    :param segmentAxis: axis perpendicular to segment plane.
    :param sx: coordinates of points on central path.
    :param sd1: derivatives of points along central path.
    :param sd2: derivatives representing cross axes.
    :param elementsCountAround: Number of elements around segment.
    :param elementsCountAlong: Number of elements along segment.
    :param nSegment: Segment index along central path.
    :return coordinates and derivatives of warped points.
    """
    xWarpedList = []
    d1WarpedList = []
    d2WarpedList = []
    smoothd2WarpedList = []
    d3WarpedUnitList = []

    for nAlongSegment in range(elementsCountAlongSegment + 1):
        n2 = elementsCountAlongSegment * nSegment + nAlongSegment
        xElementAlongSegment = xList[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d1ElementAlongSegment = d1List[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d2ElementAlongSegment = d2List[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        xMid = [0.0, 0.0, segmentLength/elementsCountAlongSegment* nAlongSegment]

        # Rotate to align segment axis with tangent of central line
        unitTangent = vector.normalise(sd1[n2])
        cp = vector.crossproduct3(segmentAxis, unitTangent)
        dp = vector.dotproduct(segmentAxis, unitTangent)
        if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
            axisRot = vector.normalise(cp)
            thetaRot = math.acos(vector.dotproduct(segmentAxis, unitTangent))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
        else: # path tangent parallel to segment axis (z-axis)
            if dp == -1.0: # path tangent opposite direction to segment axis
                thetaRot = math.pi
                axisRot = [1.0, 0, 0]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
            else: # segment axis in same direction as unit tangent
                midRot = xMid
        translateMatrix = [sx[n2][j] - midRot[j] for j in range(3)]

        for n1 in range(elementsCountAround):
            x = xElementAlongSegment[n1]
            d1 = d1ElementAlongSegment[n1]
            d2 = d2ElementAlongSegment[n1]

            if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]

                if n1 == 0:
                    # Project sd2 onto plane normal to sd1
                    v = sd2[n2]
                    pt = [midRot[j] + sd2[n2][j] for j in range(3)]
                    dist = vector.dotproduct(v, unitTangent)
                    ptOnPlane = [pt[j] - dist*unitTangent[j] for j in range(3)]
                    newVector = [ptOnPlane[j] - midRot[j] for j in range(3)]
                    # Rotate first point to align with planar projection of sd2
                    firstVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    thetaRot2 = math.acos(vector.dotproduct(vector.normalise(newVector), firstVector))
                    cp2 = vector.crossproduct3(vector.normalise(newVector), firstVector)
                    if vector.magnitude(cp2) > 0.0:
                        cp2 = vector.normalise(cp2)
                        signThetaRot2 = vector.dotproduct(unitTangent, cp2)
                        axisRot2 = unitTangent
                        rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                    else:
                        rotFrame2 = [ [1, 0, 0], [0, 1, 0], [0, 0, 1]]

            else: # path tangent parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)] if dp == -1.0 else x
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)] if dp == -1.0 else d1
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)] if dp == -1.0 else d2

                # Rotate to align start of elementsAround with sd2
                if n1 == 0:
                    v = vector.normalise(sd2[n2])
                    startVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    axisRot2 = unitTangent
                    thetaRot2 = dp*-math.acos(vector.dotproduct(v, startVector))
                    rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, thetaRot2)

            xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
            d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
            d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
            xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

            xWarpedList.append(xTranslate)
            d1WarpedList.append(d1Rot2)
            d2WarpedList.append(d2Rot2)

    # Smooth d2 for segment
    smoothd2Raw = []
    for n1 in range(elementsCountAround):
        nx = []
        nd2 = []
        for n2 in range(elementsCountAlongSegment + 1):
            n = n2*elementsCountAround + n1
            nx.append(xWarpedList[n])
            nd2.append(d2WarpedList[n])
        smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
        smoothd2Raw.append(smoothd2)

    # Re-arrange smoothd2
    for n2 in range(elementsCountAlongSegment + 1):
        for n1 in range(elementsCountAround):
            smoothd2WarpedList.append(smoothd2Raw[n1][n2])

    # Calculate unit d3
    for n in range(len(xWarpedList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1WarpedList[n]), vector.normalise(smoothd2WarpedList[n])))
        d3WarpedUnitList.append(d3Unit)

    return xWarpedList, d1WarpedList, smoothd2WarpedList, d3WarpedUnitList


def getCoordinatesFromInner(xInner, d1Inner, d2Inner, d3Inner,
    sx, wallThicknessList, elementsCountAround,
    elementsCountAlong, elementsCountThroughWall, transitElementList):
    """
    Generates coordinates from inner to outer surface using coordinates
    and derivatives of inner surface.
    :param xInner: Coordinates on inner surface
    :param d1Inner: Derivatives on inner surface around tube
    :param d2Inner: Derivatives on inner surface along tube
    :param d3Inner: Derivatives on inner surface through wall
    :param sx: coordinates of points on central path.
    :param wallThicknessList: Wall thickness for each element along tube
    :param elementsCountAround: Number of elements around tube
    :param elementsCountAlong: Number of elements along tube
    :param elementsCountThroughWall: Number of elements through tube wall
    :param transitElementList: stores true if element around is a transition
    element that is between a big and a small element.
    return nodes and derivatives for mesh, and curvature along inner surface.
    """

    xOuter = []
    curvatureAroundInner = []
    curvatureAlong = []
    curvatureList = []
    xList = []
    d1List = []
    d2List = []
    d3List = []

    for n2 in range(elementsCountAlong + 1):
        wallThickness = wallThicknessList[n2]
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3Inner[n]
            # Calculate outer coordinates
            x = [xInner[n][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)
            # Calculate curvature along elements around
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            kappam = interp.getCubicHermiteCurvatureSimple(xInner[prevIdx], d1Inner[prevIdx], xInner[n], d1Inner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xInner[n], d1Inner[n], xInner[nextIdx], d1Inner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            curvatureAroundInner.append(curvatureAround)

            # Calculate curvature along
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[n + elementsCountAround], d2Inner[n + elementsCountAround], vector.normalise(d3Inner[n]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(xInner[n - elementsCountAround], d2Inner[n - elementsCountAround], xInner[n], d2Inner[n], vector.normalise(d3Inner[n]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(xInner[n - elementsCountAround], d2Inner[n - elementsCountAround], xInner[n], d2Inner[n], vector.normalise(d3Inner[n]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[n + elementsCountAround], d2Inner[n + elementsCountAround], vector.normalise(d3Inner[n]), 0.0)))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for n1 in range(elementsCountAround):
                n = n2*elementsCountAround + n1
                norm = d3Inner[n]
                innerx = xInner[n]
                outerx = xOuter[n]
                dWall = [wallThickness*c for c in norm]
                # x
                x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
                xList.append(x)

                # dx_ds1
                factor = 1.0 + wallThickness*xi3 * curvatureAroundInner[n]
                d1 = [ factor*c for c in d1Inner[n]]
                d1List.append(d1)

                # dx_ds2
                curvature = curvatureAlong[n]
                distance = vector.magnitude([x[i] - xInner[n][i] for i in range(3)])
                factor = 1.0 - curvature*distance
                d2 = [ factor*c for c in d2Inner[n]]
                d2List.append(d2)
                curvatureList.append(curvature)

                #dx_ds3
                d3 = [c * wallThickness/elementsCountThroughWall for c in norm]
                d3List.append(d3)

    return xList, d1List, d2List, d3List, curvatureList

def createFlatAndTextureCoordinates(xiList, lengthAroundList,
    totalLengthAlong, wallThickness, elementsCountAround,
    elementsCountAlong, elementsCountThroughWall, transitElementList):
    """
    Calculates flat coordinates and texture coordinates
    for a tube when it is opened into a flat preparation.
    :param xiList: List containing xi for each point around
    the outer surface of the tube.
    :param lengthAroundList: List of total arclength around the
    outer surface for each element along.
    :param totalLengthAlong: Total length along tube.
    :param wallThickness: Thickness of wall.
    :param elementsCountAround: Number of elements around tube.
    :param elementsCountAlong: Number of elements along tube.
    :param elementsCountThroughWall: Number of elements through wall.
    :param transitElementList: stores true if element around is a
    transition element between a big and small element.
    :return: coordinates and derivatives of flat and texture coordinates fields.
    """

    # Calculate flat coordinates and derivatives
    xFlatList = []
    d1FlatList = []
    d2FlatList = []
    for n2 in range(elementsCountAlong + 1):
        xiFace = xiList[n2]
        lengthAround = lengthAroundList[n2]
        d1List = []
        for n1 in range(len(xiFace)):
            d1 = (xiFace[n1] - xiFace[n1-1]) if n1 > 0 else (xiFace[n1+1] - xiFace[n1])
            d1List.append(d1)

        # To modify derivative along transition elements
        for i in range(len(transitElementList)):
            if transitElementList[i]:
                d1List[i+1] = d1List[i+2]

        xPad = (lengthAroundList[0] - lengthAround)*0.5
        for n3 in range(elementsCountThroughWall + 1):
            z = wallThickness / elementsCountThroughWall * n3
            for n1 in range(elementsCountAround + 1):
                xFlat = [xPad + xiFace[n1] * lengthAround,
                        totalLengthAlong / elementsCountAlong * n2,
                        z]
                d1Flat = [ d1List[n1]*lengthAround, 0.0, 0.0 ]
                xFlatList.append(xFlat)
                d1FlatList.append(d1Flat)

    for n2 in range(elementsCountAlong):
        for n3 in range(elementsCountThroughWall + 1):
            for n1 in range(elementsCountAround + 1 ):
                nodeIdx = n2*(elementsCountAround + 1)*(elementsCountThroughWall + 1) + n3*(elementsCountAround + 1) + n1
                nodeNextElementAlong = nodeIdx + (elementsCountAround+1)*(elementsCountThroughWall + 1)
                # print(nodeIdx + 1, nodeNextElementAlong + 1)
                v1 = xFlatList[nodeNextElementAlong]
                v2 = xFlatList[nodeIdx]
                d1 = d2 = [v1[i] - v2[i] for i in range(3)]
                arclength = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                d2Flat = vector.setMagnitude(d1, arclength)
                d2FlatList.append(d2Flat)
    d2FlatList = d2FlatList + d2FlatList[-(elementsCountAround+1)*(elementsCountThroughWall+1):]

    # Calculate texture coordinates and derivatives
    xTextureList = []
    d1Texture = []
    d1TextureList = []
    d2TextureList = []

    d2 = [0.0, 1.0 / elementsCountAlong, 0.0]
    xiTexture = xiList[0]

    for n1 in range(len(xiTexture)):
        d1 = [xiTexture[n1] - xiTexture[n1-1] if n1 > 0 else xiTexture[n1+1] - xiTexture[n1],
              0.0,
              0.0]
        d1Texture.append(d1)

    # To modify derivative along transition elements
    for i in range(len(transitElementList)):
        if transitElementList[i]:
            d1Texture[i+1] = d1Texture[i+2]

    for n2 in range(elementsCountAlong + 1):
        for n3 in range(elementsCountThroughWall + 1):
            for n1 in range(elementsCountAround + 1):
                u = [ xiTexture[n1],
                      1.0 / elementsCountAlong * n2,
                      1.0 / elementsCountThroughWall * n3]
                d1 = d1Texture[n1]
                xTextureList.append(u)
                d1TextureList.append(d1)
                d2TextureList.append(d2)

    return xFlatList, d1FlatList, d2FlatList, xTextureList, d1TextureList, d2TextureList

def createNodesAndElements(region,
    x, d1, d2, d3,
    xFlat, d1Flat, d2Flat,
    xTexture, d1Texture, d2Texture,
    elementsCountAround, elementsCountAlong, elementsCountThroughWall,
    annotationGroups, annotationArray,
    firstNodeIdentifier, firstElementIdentifier,
    useCubicHermiteThroughWall, useCrossDerivatives):
    """
    Create nodes and elements for the coordinates, flat coordinates,
    and texture coordinates field.
    :param x, d1, d2, d3: coordinates and derivatives of coordinates field.
    :param xFlat, d1Flat, d2Flat, d3Flat: coordinates and derivatives of
    flat coordinates field.
    :param xTexture, d1Texture, d2Texture, d3Texture: coordinates and derivatives
    of texture coordinates field.
    :param elementsCountAround: Number of elements around tube.
    :param elementsCountAlong: Number of elements along tube.
    :param elementsCountThroughWall: Number of elements through wall.
    :param annotationGroups: stores information about annotation groups.
    :param annotationArray: stores annotation names of elements around.
    :param firstNodeIdentifier, firstElementIdentifier: first node and
    element identifier to use.
    :param useCubicHermiteThroughWall: use linear when false
    :param useCrossDerivatives: use cross derivatives when true
    :return nodeIdentifier, elementIdentifier, annotationGroups
    """

    nodeIdentifier = firstNodeIdentifier
    elementIdentifier = firstElementIdentifier
    zero = [ 0.0, 0.0, 0.0 ]

    fm = region.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()

    # Coordinates field
    coordinates = zinc_utils.getOrCreateCoordinateField(fm)
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

    # Flat coordinates field
    bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eftTexture1 = bicubichermitelinear.createEftBasic()
    eftTexture2 = bicubichermitelinear.createEftOpenTube()

    flatCoordinates = zinc_utils.getOrCreateFlatCoordinateField(fm)
    flatNodetemplate1 = nodes.createNodetemplate()
    flatNodetemplate1.defineField(flatCoordinates)
    flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

    flatNodetemplate2 = nodes.createNodetemplate()
    flatNodetemplate2.defineField(flatCoordinates)
    flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_VALUE, 2)
    flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
    flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS2, 2)
    if useCrossDerivatives:
        flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 2)

    flatElementtemplate1 = mesh.createElementtemplate()
    flatElementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    flatElementtemplate1.defineField(flatCoordinates, -1, eftTexture1)

    flatElementtemplate2 = mesh.createElementtemplate()
    flatElementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    flatElementtemplate2.defineField(flatCoordinates, -1, eftTexture2)

    # Texture coordinates field
    textureCoordinates = zinc_utils.getOrCreateTextureCoordinateField(fm)
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

    elementtemplate1 = mesh.createElementtemplate()
    elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate1.defineField(textureCoordinates, -1, eftTexture1)

    elementtemplate2 = mesh.createElementtemplate()
    elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate2.defineField(textureCoordinates, -1, eftTexture2)

    # Create nodes
    # Coordinates field
    for n in range(len(x)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        # print('NodeIdentifier = ', nodeIdentifier, xList[n])
        nodeIdentifier = nodeIdentifier + 1

    # Flat and texture coordinates fields
    nodeIdentifier = firstNodeIdentifier
    for n2 in range(elementsCountAlong + 1):
        for n3 in range(elementsCountThroughWall + 1):
            for n1 in range(elementsCountAround):
                i = n2*(elementsCountAround + 1)*(elementsCountThroughWall + 1) + (elementsCountAround + 1)*n3 + n1
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                node.merge(flatNodetemplate2 if n1 == 0 else flatNodetemplate1)
                node.merge(textureNodetemplate2 if n1 == 0 else textureNodetemplate1)
                cache.setNode(node)
                # print('NodeIdentifier', nodeIdentifier, 'version 1, xList Index =', i+1)
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xFlat[i])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Flat[i])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Flat[i])
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTexture[i])
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Texture[i])
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Texture[i])
                if useCrossDerivatives:
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                if n1 == 0:
                    # print('NodeIdentifier', nodeIdentifier, 'version 2, xList Index =', i+elementsCountAround+1)
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2, xFlat[i+elementsCountAround])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1Flat[i+elementsCountAround])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2, d2Flat[i+elementsCountAround])
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2, xTexture[i+elementsCountAround])
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1Texture[i+elementsCountAround])
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2, d2Texture[i+elementsCountAround])
                    if useCrossDerivatives:
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
                        textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
                nodeIdentifier = nodeIdentifier + 1

    # create elements
    now = elementsCountAround*(elementsCountThroughWall+1)
    for e2 in range(elementsCountAlong):
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                bni11 = e2*now + e3*elementsCountAround + e1 + 1
                bni12 = e2*now + e3*elementsCountAround + (e1 + 1) % elementsCountAround + 1
                bni21 = e2*now + (e3 + 1)*elementsCountAround + e1 + 1
                bni22 = e2*now + (e3 + 1)*elementsCountAround + (e1 + 1) % elementsCountAround + 1
                nodeIdentifiers = [ bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now ]
                onOpening = e1 > elementsCountAround - 2
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                element.merge(flatElementtemplate2 if onOpening else flatElementtemplate1)
                element.merge(elementtemplate2 if onOpening else elementtemplate1)
                element.setNodesByIdentifier(eftTexture2 if onOpening else eftTexture1, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1
                if annotationGroups:
                    for annotationGroup in annotationGroups:
                        if annotationArray[e1] == annotationGroup._name:
                            meshGroup = annotationGroup.getMeshGroup(mesh)
                            meshGroup.addElement(element)

    fm.endChange()

    return nodeIdentifier, elementIdentifier, annotationGroups
