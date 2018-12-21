'''
Utility function for generating tubular mesh from a central line.
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
    around, along the central axis, and radially through wall. The
    ellipsoidal tubular mesh is created from central line and lateral
    axes data
    :param elementsCountAlong: number of elements along tube
    :param elementsCountAround: number of elements around tube
    :param elementsCountThroughWall: number of elements through wall thickness
    :param cx: coordinates on central line
    :param cd1: derivative 1 along central line
    :param cd2: derivative 2 to the lateral axis
    :param cd3: derivative 3 to the lateral axis perpendicular to direction of cd2
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
    # sd1Smoothed = smoothCubicHermiteDerivativesLine(sx, sd1)
    # sd1 = []
    # sd1 = sd1Smoothed
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
            rotFrame = rotationMatrixAboutAxis(axisRot, thetaRot)
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
    innerx = []
    innerdx_ds1 = []
    outerx = []
    dWall = []
    unitNormalInner = []
    curvatureAlong = []
    ellipsex = []
    ellipsedx_ds1 = []
    tubex = []
    tubedx_ds1 = []
    tubedx_ds2 = []
    tubedx_ds3 = []

    # Calculate node location and derivatives for inner and outer surface
    for n2 in range(elementsCountAlong+1):
        aInner = sd2[n2][1]
        bInner = sd3[n2][2]
        perimeterInner = getApproximateEllipsePerimeter(aInner, bInner)
        arcLengthPerElementAroundInner = perimeterInner / elementsCountAround
        prevRadiansAroundInner = updateEllipseAngleByArcLength(aInner, bInner, 0.0, arcLengthPerElementAroundInner)

        aOuter = sd2[n2][1] + st2[n2]
        bOuter = sd3[n2][2] + st3[n2]
        perimeterOuter = getApproximateEllipsePerimeter(aOuter, bOuter)
        arcLengthPerElementAroundOuter = perimeterOuter / elementsCountAround
        prevRadiansAroundOuter = updateEllipseAngleByArcLength(aOuter, bOuter, 0.0, arcLengthPerElementAroundOuter)

        for n1 in range(elementsCountAround):
            arcLengthAroundInner = n1*arcLengthPerElementAroundInner
            radiansAroundInner = -1*updateEllipseAngleByArcLength(aInner, bInner, 0.0, arcLengthAroundInner)
            cosRadiansAroundInner = math.cos(radiansAroundInner)
            sinRadiansAroundInner = math.sin(radiansAroundInner)
            xInner = [sx[n2][j] + aInner*cosRadiansAroundInner*sBinormal[n2][j] + bInner*sinRadiansAroundInner*sNormal[n2][j] for j in range(3)]
            innerx.append(xInner)
            dx_ds1Inner = [(radiansAroundInner - prevRadiansAroundInner)*(aInner*-sinRadiansAroundInner*sBinormal[n2][j] + bInner*cosRadiansAroundInner*sNormal[n2][j]) for j in range(3)]
            innerdx_ds1.append(dx_ds1Inner)
            prevRadiansAroundInner = radiansAroundInner
            unitNormalInnerNode = normalise([aInner*cosRadiansAroundInner*sBinormal[n2][j] + bInner*sinRadiansAroundInner*sNormal[n2][j] for j in range(3)])
            unitNormalInner.append(unitNormalInnerNode)

            arcLengthAroundOuter = n1*arcLengthPerElementAroundOuter
            radiansAroundOuter = -1*updateEllipseAngleByArcLength(aOuter, bOuter, 0.0, arcLengthAroundOuter)
            cosRadiansAroundOuter = math.cos(radiansAroundOuter)
            sinRadiansAroundOuter = math.sin(radiansAroundOuter)
            xOuter = [sx[n2][j] + aOuter*cosRadiansAroundOuter*sBinormal[n2][j] + bOuter*sinRadiansAroundOuter*sNormal[n2][j] for j in range(3)]
            outerx.append(xOuter)
            prevRadiansAroundOuter = radiansAroundOuter

            d = [ xOuter[i] - xInner[i] for i in range(3)]
            dWall.append(d)

            if n2 == 0:
                curvature = getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2+1], sd1[n2+1], unitNormalInnerNode, 0.0)
            elif n2 == elementsCountAlong:
                curvature = getCubicHermiteCurvature(sx[n2-1], sd1[n2-1], sx[n2], sd1[n2], unitNormalInnerNode, 1.0)
            else:
                curvature = 0.5*(
                    getCubicHermiteCurvature(sx[n2-1], sd1[n2-1], sx[n2], sd1[n2], unitNormalInnerNode, 1.0) +
                    getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2+1], sd1[n2+1], unitNormalInnerNode, 0.0))
            curvatureAlong.append(curvature)

    # Interpolate to get nodes through wall
    for n3 in range(elementsCountThroughWall + 1):
        xi = 1/elementsCountThroughWall*n3

        for n2 in range(elementsCountAlong+1):
            for n1 in range(elementsCountAround):
                n = n2*elementsCountAround + n1

                # x
                x = interpolateCubicHermite(innerx[n], dWall[n], outerx[n], dWall[n], xi)
                tubex.append(x)

                # dx_ds1
                if n1 == 0:
                    prevIdx = (((n+1)//elementsCountAround + 1)*elementsCountAround)-1
                    ellipsex = []
                    ellipsedx_ds1 = []
                else:
                    prevIdx = n - 1
                nextIdx = ((n+1)//(elementsCountAround +1) * elementsCountAround + ((n+1)%elementsCountAround) + 1)-1
                curvatureAround = 0.5*(
                    getCubicHermiteCurvature(innerx[prevIdx], innerdx_ds1[prevIdx], innerx[n], innerdx_ds1[n], unitNormalInner[n], 1.0) +
                    getCubicHermiteCurvature(innerx[n], innerdx_ds1[n], innerx[nextIdx], innerdx_ds1[nextIdx], unitNormalInner[n], 0.0))
                wallDistance = magnitude(dWall[n])*xi
                factor = 1.0 - curvatureAround*wallDistance
                dx_ds1 = [ factor*c for c in innerdx_ds1[n]]
                ellipsedx_ds1.append(dx_ds1)
                ellipsex.append(x)
                if n1 == elementsCountAround-1:
                    smoothdx_ds1 = smoothCubicHermiteDerivativesLoop(ellipsex, ellipsedx_ds1)
                    tubedx_ds1 = tubedx_ds1 + smoothdx_ds1

                # dx_ds2
                factor = 1.0 - curvatureAlong[n]*wallDistance 
                d1AlongTube = [ factor*c for c in sd1[n2]]
                tubedx_ds2.append(d1AlongTube)

                #dx_ds3
                dx_ds3 = interpolateCubicHermiteDerivative(innerx[n], dWall[n], outerx[n], dWall[n], xi)
                dx_ds3 = [c/elementsCountThroughWall for c in dx_ds3]
                tubedx_ds3.append(dx_ds3)

    for n in range((elementsCountAround)*(elementsCountAlong+1)*(elementsCountThroughWall+1)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, tubex[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, tubedx_ds1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, tubedx_ds2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, tubedx_ds3[n])
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

def rotationMatrixAboutAxis(rotAxis, theta):
    """
    Generate the rotation matrix for rotation about an axis.
    :param rotAxis: axis of rotation
    :param theta: angle of rotation
    :return: rotation matrix
    """
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)
    C = 1 - cosTheta
    rotMatrix = ([[rotAxis[0]*rotAxis[0]*C + cosTheta, rotAxis[0]*rotAxis[1]*C - rotAxis[2]*sinTheta, rotAxis[0]*rotAxis[2]*C + rotAxis[1]*sinTheta],
        [rotAxis[1]*rotAxis[0]*C + rotAxis[2]*sinTheta, rotAxis[1]*rotAxis[1]*C + cosTheta, rotAxis[1]*rotAxis[2]*C - rotAxis[0]*sinTheta],
        [rotAxis[2]*rotAxis[0]*C - rotAxis[1]*sinTheta, rotAxis[2]*rotAxis[1]*C + rotAxis[0]*sinTheta, rotAxis[2]*rotAxis[2]*C + cosTheta]])
    return rotMatrix
