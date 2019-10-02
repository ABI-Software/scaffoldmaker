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

def generatetubemesh(region,
    elementsCountAround,
    elementsCountAlongSegment,
    elementsCountThroughWall,
    segmentCountAlong,
    cx, cd1, cd2, cd12,
    tubeMeshSegmentInnerPoints,
    wallThickness,
    segmentLength,
    useCrossDerivatives,
    useCubicHermiteThroughWall, # or Zinc Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE etc.
    firstNodeIdentifier = 1, firstElementIdentifier = 1):
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
    :param cd2: derivative representing cross axis
    :param cd12: rate of change of cd2 along cd1
    :param tubeMeshSegmentInnerPoints: class function for generating
    coordinates and derivatives on inner surface of segment profile
    :param wallThickness: thickness of wall
    :param segmentLength: length of segment profile
    :param useCubicHermiteThroughWall: use linear when false
    :return nodeIdentifier, elementIdentifier
    :return xList, dx_ds1List, dx_ds2List, dx_ds3List: list of coordinates and derivatives
    on tube
    :return sx: list of coordinates sampled from central line
    :return curvatureAlong, factorList: list of curvature and scale factor along mesh
    for each node on inner surface of mesh.
    :return uList: list of xi for each node around mid-length haustra.
    :return relaxedLengthList: list of lengths around elements along tube length.
    '''

    zero  = [0.0, 0.0, 0.0]
    elementsCountAlong = elementsCountAlongSegment*segmentCountAlong

    # Sample central line to get same number of elements as elementsCountAlong
    sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
    sd2, _ = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

    fm = region.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()
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
    curvatureAlong = []
    smoothd2Raw = []
    smoothd2InnerList = []
    d1List = []
    flatWidthListOuter = []
    uList = []
    relaxedLengthList = []

    for nSegment in range(segmentCountAlong):
        # Create inner points
        annotationGroups, annotationArray, transitElementList, uSegment, relaxedLengthSegment, xInner, d1Inner, d2Inner, segmentAxis, sRadius = tubeMeshSegmentInnerPoints.getTubeMeshSegmentInnerPoints(
            nSegment)

        startIdx = 0 if nSegment == 0 else 1
        u = uSegment[startIdx:len(sRadius)]
        uList.append(u)

        relaxedLength = relaxedLengthSegment[startIdx:len(sRadius)]
        relaxedLengthList.append(relaxedLength)

        # Map each face along inner points to central line
        for nAlongSegment in range(elementsCountAlongSegment + 1):
            n2 = elementsCountAlongSegment*nSegment + nAlongSegment
            xFace = xInner[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
            d1Face = d1Inner[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
            d2Face = d2Inner[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]

            if nSegment == 0 or (nSegment > 0 and nAlongSegment > 0): # take whole of 1st segment and then all along except last for rest
                # Rotate to align segment axis with tangent of central line
                segmentMid = [0.0, 0.0, segmentLength/elementsCountAlongSegment* nAlongSegment]
                unitTangent = vector.normalise(sd1[n2])
                cp = vector.crossproduct3(segmentAxis, unitTangent)
                dp = vector.dotproduct(segmentAxis, unitTangent)
                if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                    axisRot = vector.normalise(cp)
                    thetaRot = math.acos(vector.dotproduct(segmentAxis, unitTangent))
                    rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                    midRot = [rotFrame[j][0]*segmentMid[0] + rotFrame[j][1]*segmentMid[1] + rotFrame[j][2]*segmentMid[2] for j in range(3)]
                else: # path tangent parallel to segment axis (z-axis)
                    if dp == -1.0: # path tangent opposite direction to segment axis
                        thetaRot = math.pi
                        axisRot = [1.0, 0, 0]
                        rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                        midRot = [rotFrame[j][0]*segmentMid[0] + rotFrame[j][1]*segmentMid[1] + rotFrame[j][2]*segmentMid[2] for j in range(3)]
                    else: # segment axis in same direction as unit tangent
                        midRot = segmentMid
                translateMatrix = [sx[n2][j] - midRot[j] for j in range(3)]

                for n1 in range(elementsCountAround):
                    n = nAlongSegment*elementsCountAround + n1
                    x = xFace[n1]
                    d1 = d1Face[n1]
                    d2 = d2Face[n1]
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

                    xInnerList.append(xTranslate)
                    d1InnerList.append(d1Rot2)
                    d2InnerList.append(d2Rot2)
                    d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1Rot2), vector.normalise(d2Rot2)))
                    d3InnerUnitList.append(d3Unit)

        # Smooth d2 for segment
        smoothd2Raw = []
        for n1 in range(elementsCountAround):
            nx = []
            nd2 = []
            nIdx = []
            for n2 in range(elementsCountAlongSegment + 1):
                n = nSegment*(elementsCountAlongSegment+1 if nSegment == 0 else elementsCountAlongSegment)*elementsCountAround + n2*elementsCountAround + n1
                nx.append(xInnerList[n])
                nd2.append(d2InnerList[n])
            smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
            smoothd2Raw.append(smoothd2)

        # Re-arrange smoothd2
        for n2 in range(startIdx, elementsCountAlongSegment + 1):
            for n1 in range(elementsCountAround):
                smoothd2InnerList.append(smoothd2Raw[n1][n2])

    for n2 in range(elementsCountAlong + 1):
        for n1 in range(elementsCountAround):
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2+1], sd1[n2+1], vector.normalise(sd2[n2]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(sx[n2-1], sd1[n2-1], sx[n2], sd1[n2], vector.normalise(sd2[n2]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(sx[n2-1], sd1[n2-1], sx[n2], sd1[n2], vector.normalise(sd2[n2]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2+1], sd1[n2+1], vector.normalise(sd2[n2]), 0.0)))
            curvatureAlong.append(curvature)

    # Pre-calculate node locations and derivatives on outer boundary
    xOuterList, curvatureInner = getOuterCoordinatesAndCurvatureFromInner(xInnerList, d1InnerList, d3InnerUnitList, wallThickness, elementsCountAlong, elementsCountAround, transitElementList)

    # Interpolate to get nodes through wall
    for n3 in range(elementsCountThroughWall + 1):
        xi3 = 1/elementsCountThroughWall * n3
        x, dx_ds1, dx_ds2, dx_ds3, factorList = interpolatefromInnerAndOuter(xInnerList, xOuterList,
            wallThickness, xi3, sx, curvatureInner, curvatureAlong, d1InnerList, smoothd2InnerList, d3InnerUnitList,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall)
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
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, vector.normalise(sd2[pt]))
        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, vector.normalise(testRot[pt]))
        # nodeIdentifier = nodeIdentifier + 1

    # create elements
    elementIdentifier = firstElementIdentifier
    now = (elementsCountAlong + 1)*elementsCountAround
    for e3 in range(elementsCountThroughWall):
        for e2 in range(elementsCountAlong):
            for e1 in range(elementsCountAround):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni11 = e3*now + e2*elementsCountAround + e1 + 1
                bni12 = e3*now + e2*elementsCountAround + (e1 + 1) % elementsCountAround + 1
                bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
                bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1) % elementsCountAround + 1
                nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1
                if annotationGroups:
                    for annotationGroup in annotationGroups:
                        if annotationArray[e1] == annotationGroup._name:
                            meshGroup = annotationGroup.getMeshGroup(mesh)
                            meshGroup.addElement(element)

    # Define texture coordinates field
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

    bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eftTexture1 = bicubichermitelinear.createEftBasic()

    elementtemplate1 = mesh.createElementtemplate()
    elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate1.defineField(textureCoordinates, -1, eftTexture1)

    eftTexture2 = bicubichermitelinear.createEftOpenTube()
    elementtemplate2 = mesh.createElementtemplate()
    elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate2.defineField(textureCoordinates, -1, eftTexture2)

    # Calculate texture coordinates and derivatives
    d2 = [0.0, 1.0 / elementsCountAlong, 0.0]
    uTexture = uList[0][0]

    for n1 in range(len(uTexture)):
        d1 = [uTexture[n1] - uTexture[n1-1] if n1 > 0 else uTexture[n1+1] - uTexture[n1],
              0.0,
              0.0]
        d1List.append(d1)

    # To modify derivative along transition elements
    for i in range(len(transitElementList)):
        if transitElementList[i]:
            d1List[i+1] = d1List[i+2]

    nodeIdentifier = firstNodeIdentifier
    for n3 in range(elementsCountThroughWall + 1):
        for n2 in range(elementsCountAlong + 1):
            for n1 in range(elementsCountAround):
                u = [ uTexture[n1],
                      1.0 / elementsCountAlong * n2,
                      1.0 / elementsCountThroughWall * n3]
                d1 = d1List[n1]
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                node.merge(textureNodetemplate2 if n1 == 0 else textureNodetemplate1)
                cache.setNode(node)
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, u)
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                if useCrossDerivatives:
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                if n1 == 0:
                    u = [ 1.0, 1.0 / elementsCountAlong * n2, 1.0 / elementsCountThroughWall * n3]
                    d1 = d1List[-1]
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2, u)
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1)
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2, d2)
                    if useCrossDerivatives:
                        textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
                nodeIdentifier = nodeIdentifier + 1

    # Define flat coordinates field
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

    totalLengthAlong = segmentLength*segmentCountAlong

    v = []
    xFlatList = []
    d1FlatList = []
    d2FlatList = []
    nodeIdentifier = firstNodeIdentifier

    for n3 in range(elementsCountThroughWall + 1):
        z = wallThickness / elementsCountThroughWall * n3
        faceCount = 0
        for nSegment in range(segmentCountAlong):
            uSegment = uList[nSegment]
            relaxedLengthSegment = relaxedLengthList[nSegment]
            for n2 in range(elementsCountAlongSegment + 1 if nSegment == 0 else elementsCountAlongSegment):
                uFace = uSegment[n2]
                relaxedLength = relaxedLengthSegment[n2]
                d1List = []
                for n1 in range(len(uFace)):
                    d1 = (uFace[n1] - uFace[n1-1]) if n1 > 0 else (uFace[n1+1] - uFace[n1])
                    d1List.append(d1)

                # To modify derivative along transition elements
                for i in range(len(transitElementList)):
                    if transitElementList[i]:
                        d1List[i+1] = d1List[i+2]

                xPad = (relaxedLengthList[0][0] - relaxedLength)*0.5
                for n1 in range(elementsCountAround + 1):
                    xFlat = [xPad + uFace[n1] * relaxedLength,
                            totalLengthAlong / elementsCountAlong * faceCount,
                            z]
                    d1Flat = [ d1List[n1]*relaxedLength, 0.0, 0.0 ]
                    xFlatList.append(xFlat)
                    d1FlatList.append(d1Flat)
                faceCount += 1

    for n3 in range(elementsCountThroughWall + 1):
        for n2 in range(elementsCountAlong):
            for n1 in range(elementsCountAround + 1 ):
                nodeIdx = n3*(elementsCountAround + 1)*(elementsCountAlong + 1) + n2*(elementsCountAround + 1) + n1
                nodeNextElementAlong = nodeIdx + (elementsCountAround+1)
                # print(nodeIdx + 1, nodeNextElementAlong + 1)
                v1 = xFlatList[nodeNextElementAlong]
                v2 = xFlatList[nodeIdx]
                d1 = d2 = [v1[i] - v2[i] for i in range(3)]
                arclength = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                d2Flat = vector.setMagnitude(d1, arclength)
                d2FlatList.append(d2Flat)
        d2FlatList = d2FlatList + d2FlatList[-elementsCountAround-1:]

    # Create nodes for flat coordinate field
    for n3 in range(elementsCountThroughWall + 1):
        for n2 in range(elementsCountAlong + 1):
            for n1 in range(elementsCountAround):
                i = (elementsCountAround + 1)*(elementsCountAlong + 1)*n3 + (elementsCountAround + 1)*n2 + n1
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                node.merge(flatNodetemplate2 if n1 == 0 else flatNodetemplate1)
                cache.setNode(node)
                # print('NodeIdentifier', nodeIdentifier, 'version 1, xList Index =', i+1)
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xFlatList[i])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1FlatList[i])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2FlatList[i])
                if useCrossDerivatives:
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                if n1 == 0:
                    # print('NodeIdentifier', nodeIdentifier, 'version 2, xList Index =', i+elementsCountAround+1)
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2, xFlatList[i+elementsCountAround])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1FlatList[i+elementsCountAround])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2, d2FlatList[i+elementsCountAround])
                    if useCrossDerivatives:
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
                nodeIdentifier = nodeIdentifier + 1

    # Define flat coordinates field & texture coordinates field over elements
    elementIdentifier = firstElementIdentifier
    now = (elementsCountAlong + 1)*elementsCountAround

    for e3 in range(elementsCountThroughWall):
        for e2 in range(elementsCountAlong):
            for e1 in range(elementsCountAround):
                onOpening = e1 > elementsCountAround - 2
                element = mesh.findElementByIdentifier(elementIdentifier)
                element.merge(elementtemplate2 if onOpening else elementtemplate1)
                element.merge(flatElementtemplate2 if onOpening else flatElementtemplate1)
                bni11 = e3*now + e2*elementsCountAround + e1 + 1
                bni12 = e3*now + e2*elementsCountAround + (e1 + 1) % elementsCountAround + 1
                bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
                bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1) % elementsCountAround + 1
                nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                element.setNodesByIdentifier(eftTexture2 if onOpening else eftTexture1, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

    fm.endChange()

    return annotationGroups, nodeIdentifier, elementIdentifier, xList, dx_ds1List, dx_ds2List, dx_ds3List, sx, curvatureAlong, factorList, uList, relaxedLengthList

def getOuterCoordinatesAndCurvatureFromInner(xInner, d1Inner, d3Inner, wallThickness, elementsCountAlong, elementsCountAround, transitElementList):
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
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            norm = d3Inner[n]
            kappam = interp.getCubicHermiteCurvatureSimple(xInner[prevIdx], d1Inner[prevIdx], xInner[n], d1Inner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xInner[n], d1Inner[n], xInner[nextIdx], d1Inner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            xOuter.append(x)
            curvatureInner.append(curvatureAround)

    return xOuter, curvatureInner

def interpolatefromInnerAndOuter( xInner, xOuter, thickness, xi3, sx, curvatureInner, curvatureAlong,
    d1Inner, d2Inner, d3InnerUnit, elementsCountAround, elementsCountAlong, elementsCountThroughWall):
    """
    Generate coordinates and derivatives at xi3 by interpolating with 
    inner and outer coordinates and derivatives.
    :param xInner: Coordinates on inner surface
    :param xOuter: Coordinates on outer surface
    :param thickness: Thickness of wall
    :param sx: List of coordinates sampled from central line
    :param curvatureInner: Curvature of coordinates on inner surface
    :param curvatureAlong: Curvature of coordinates on inner surface along mesh
    :param d1Inner: Derivatives on inner surface around tube
    :param d2Inner: Derivatives on inner surface along tube
    :param d3InnerUnit: Unit derivatives on inner surface through wall
    :param elementsCountAround: Number of elements around tube
    :param elementsCountAlong: Number of elements along tube
    :param elementsCountThroughWall: Number of elements through wall
    :return xList, dx_ds1List, dx_ds2List, dx_ds3List: Coordinates and derivatives on xi3
    :return factorList: List of factors used for scaling d2
    """
    xList = []
    dx_ds1List = []
    dx_ds2List = []
    dx_ds3List =[]
    factorList = []

    for n2 in range(elementsCountAlong + 1):
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3InnerUnit[n]
            # x
            innerx = xInner[n]
            outerx = xOuter[n]
            dWall = [thickness*c for c in norm]
            x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
            xList.append(x)
            # dx_ds1
            factor = 1.0 + thickness*xi3 * curvatureInner[n]
            dx_ds1 = [ factor*c for c in d1Inner[n]]
            dx_ds1List.append(dx_ds1)
            # dx_ds2
            curvature = curvatureAlong[n]
            distance = vector.magnitude([x[i] - sx[n2][i] for i in range(3)])
            factor = 1.0 + curvature*distance
            dx_ds2 = [ factor*c for c in d2Inner[n]]
            dx_ds2List.append(dx_ds2)
            factorList.append(factor)

            #dx_ds3
            dx_ds3 = [c * thickness/elementsCountThroughWall for c in norm]
            dx_ds3List.append(dx_ds3)

    return xList, dx_ds1List, dx_ds2List, dx_ds3List, factorList
