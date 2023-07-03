'''
Utility function for generating tubular mesh from a central line
using a segment profile.
'''
from __future__ import division

import math

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import mergeAnnotationGroups
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabelsVersion
from scaffoldmaker.utils.geometry import createCirclePoints


def getPlaneProjectionOnCentralPath(x, elementsCountAround, elementsCountAlong,
                                    segmentLength, sx, sd1, sd2, sd12):
    """
    Projects reference point used for warping onto the central path and find coordinates
    and derivatives at projected location.
    :param x: coordinates of nodes.
    :param elementsCountAround: number of elements around.
    :param elementsCountAlong: number of elements along.
    :param segmentLength: Length of segment.
    :param sx: coordinates of equally spaced points on central path.
    :param sd1: tangent of equally spaced points on central path.
    :param sd2: derivative representing cross axis at equally spaced points on central path.
    :param sd12: rate of change of cross axis at equally spaced points on central path.
    :return: coordinates and derivatives on project points and z-coordinates of reference points.
    """

    # Use first node in each group of elements along as reference for warping later
    zRefList = []
    for n2 in range(elementsCountAlong + 1):
        zFirstNodeAlong = x[n2*elementsCountAround][2]
        zRefList.append(zFirstNodeAlong)

    # Find sx, sd1, sd2 at projection of reference points on central path
    lengthElementAlong = segmentLength / elementsCountAlong

    # Append values from first node on central path
    sxRefList = []
    sd1RefList = []
    sd2RefList = []
    sxRefList.append(sx[0])
    sd1RefList.append(sd1[0])
    sd2RefList.append(sd2[0])

    # Interpolate the ones in between
    for n2 in range(1, elementsCountAlong):
        ei = int(zRefList[n2]//lengthElementAlong + 1)
        xi = (zRefList[n2] - lengthElementAlong*(ei-1))/lengthElementAlong
        sxRef = interp.interpolateCubicHermite(sx[ei - 1], sd1[ei - 1], sx[ei], sd1[ei], xi)
        sd1Ref = interp.interpolateCubicHermiteDerivative(sx[ei - 1], sd1[ei - 1], sx[ei], sd1[ei], xi)
        sd2Ref = interp.interpolateCubicHermite(sd2[ei - 1], sd12[ei - 1], sd2[ei], sd12[ei], xi)
        sxRefList.append(sxRef)
        sd1RefList.append(sd1Ref)
        sd2RefList.append(sd2Ref)

    # Append values from last node on central path
    sxRefList.append(sx[-1])
    sd1RefList.append(sd1[-1])
    sd2RefList.append(sd2[-1])

    # Smooth sd1
    sd1RefList = interp.smoothCubicHermiteDerivativesLine(sxRefList, sd1RefList, fixAllDirections=True,
                                                          fixStartDerivative=True, fixEndDerivative=True)

    # Project sd2 to plane orthogonal to sd1
    sd2ProjectedListRef = []

    for n in range(len(sd2RefList)):
        sd1Normalised = vector.normalise(sd1RefList[n])
        dp = vector.dotproduct(sd2RefList[n], sd1Normalised)
        dpScaled = [dp * c for c in sd1Normalised]
        sd2Projected = vector.normalise([sd2RefList[n][c] - dpScaled[c] for c in range(3)])
        sd2ProjectedListRef.append(sd2Projected)

    return sxRefList, sd1RefList, sd2ProjectedListRef, zRefList

def warpSegmentPoints(xList, d1List, d2List, segmentAxis, sx, sd1, sd2, elementsCountAround, elementsCountAlongSegment,
                      refPointZ):
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
    :param elementsCountAlongSegment: Number of elements along segment.
    :param refPointZ: z-coordinate of reference point for each element
    groups along the segment to be used for transformation.
    :return coordinates and derivatives of warped points.
    """

    xWarpedList = []
    d1WarpedList = []
    d2WarpedList = []
    d2WarpedListFinal = []
    d3WarpedUnitList = []

    for nAlongSegment in range(elementsCountAlongSegment + 1):
        xElementAlongSegment = xList[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d1ElementAlongSegment = d1List[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d2ElementAlongSegment = d2List[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]

        centroid = [0.0, 0.0, refPointZ[nAlongSegment]]

        # Rotate to align segment axis with tangent of central line
        unitTangent = vector.normalise(sd1[nAlongSegment])
        cp = vector.crossproduct3(segmentAxis, unitTangent)
        dp = vector.dotproduct(segmentAxis, unitTangent)
        if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
            axisRot = vector.normalise(cp)
            thetaRot = math.acos(vector.dotproduct(segmentAxis, unitTangent))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            centroidRot = [rotFrame[j][0]*centroid[0] + rotFrame[j][1]*centroid[1] + rotFrame[j][2]*centroid[2] for j in range(3)]

        else: # path tangent parallel to segment axis (z-axis)
            if dp == -1.0: # path tangent opposite direction to segment axis
                thetaRot = math.pi
                axisRot = [1.0, 0, 0]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                centroidRot = [rotFrame[j][0] * centroid[0] + rotFrame[j][1] * centroid[1] + rotFrame[j][2] * centroid[2] for j in range(3)]

            else: # segment axis in same direction as unit tangent
                rotFrame = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                centroidRot = centroid

        translateMatrix = [sx[nAlongSegment][j] - centroidRot[j] for j in range(3)]

        for n1 in range(elementsCountAround):
            x = xElementAlongSegment[n1]
            d1 = d1ElementAlongSegment[n1]
            d2 = d2ElementAlongSegment[n1]

            if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]
                # xTranslate = [xRot1[j] + translateMatrix[j] for j in range(3)]

            else: # path tangent parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)] if dp == -1.0 else x
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)] if dp == -1.0 else d1
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)] if dp == -1.0 else d2
                # xTranslate = [xRot1[j] + translateMatrix[j] for j in range(3)]

            if n1 == 0:  # Find angle between xCentroidRot and first node in the face
                vectorToFirstNode = [xRot1[c] - centroidRot[c] for c in range(3)]
                if vector.magnitude(vectorToFirstNode) > 0.0:
                    cp = vector.crossproduct3(vector.normalise(vectorToFirstNode), vector.normalise(sd2[nAlongSegment]))
                    if vector.magnitude(cp) > 1e-7:
                        cp = vector.normalise(cp)
                        signThetaRot2 = vector.dotproduct(unitTangent, cp)
                        thetaRot2 = math.acos(
                            vector.dotproduct(vector.normalise(vectorToFirstNode), vector.normalise(sd2[nAlongSegment])))
                        axisRot2 = unitTangent
                        rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, signThetaRot2*thetaRot2)
                    else:
                        rotFrame2 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                else:
                    rotFrame2 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

            xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
            d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
            d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
            xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

            xWarpedList.append(xTranslate)
            d1WarpedList.append(d1Rot2)
            d2WarpedList.append(d2Rot2)

    # Scale d2 with curvature of central path
    d2WarpedListScaled = []
    for nAlongSegment in range(elementsCountAlongSegment + 1):
        for n1 in range(elementsCountAround):
            n = nAlongSegment * elementsCountAround + n1
            # Calculate norm
            sd1Normalised = vector.normalise(sd1[nAlongSegment])
            v = [xWarpedList[n][c] - sx[nAlongSegment][c] for c in range(3)]
            dp = vector.dotproduct(v, sd1Normalised)
            dpScaled = [dp * c for c in sd1Normalised]
            vProjected = [v[c] - dpScaled[c] for c in range(3)]
            if vector.magnitude(vProjected) > 0.0:
                vProjectedNormalised = vector.normalise(vProjected)
            else:
                vProjectedNormalised = [0.0, 0.0, 0.0]

            # Calculate curvature along at each node
            if nAlongSegment == 0:
                curvature = interp.getCubicHermiteCurvature(sx[0], sd1[0], sx[1], sd1[1], vProjectedNormalised, 0.0)
            elif nAlongSegment == elementsCountAlongSegment:
                curvature = interp.getCubicHermiteCurvature(sx[-2], sd1[-2], sx[-1], sd1[-1], vProjectedNormalised, 1.0)
            else:
                curvature = 0.5 * (interp.getCubicHermiteCurvature(sx[nAlongSegment - 1], sd1[nAlongSegment - 1],
                                                                   sx[nAlongSegment], sd1[nAlongSegment],
                                                                   vProjectedNormalised, 1.0) +
                                   interp.getCubicHermiteCurvature(sx[nAlongSegment], sd1[nAlongSegment],
                                                                   sx[nAlongSegment + 1], sd1[nAlongSegment + 1],
                                                                   vProjectedNormalised, 0.0))
            # Scale
            if nAlongSegment < elementsCountAlongSegment:
                factor = 1.0 - curvature * vector.magnitude(v)
                d2 = [factor * c for c in d2WarpedList[n]]
                d2WarpedListScaled.append(d2)
            else:
                d2WarpedListScaled.append(d2WarpedList[n])

    # Smooth d2 for segment
    smoothd2Raw = []
    for n1 in range(elementsCountAround):
        nx = []
        nd2 = []
        for n2 in range(elementsCountAlongSegment + 1):
            n = n2*elementsCountAround + n1
            nx.append(xWarpedList[n])
            nd2.append(d2WarpedListScaled[n])
        smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative=True, fixEndDerivative=True)
        smoothd2Raw.append(smoothd2)

    # Re-arrange smoothd2
    for n2 in range(elementsCountAlongSegment + 1):
        for n1 in range(elementsCountAround):
            d2WarpedListFinal.append(smoothd2Raw[n1][n2])

    # Calculate unit d3
    for n in range(len(xWarpedList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1WarpedList[n]),
                                                       vector.normalise(d2WarpedListFinal[n])))
        d3WarpedUnitList.append(d3Unit)

    return xWarpedList, d1WarpedList, d2WarpedListFinal, d3WarpedUnitList

def getCoordinatesFromInner(xInner, d1Inner, d2Inner, d3Inner,
    wallThicknessList, relativeThicknessList, elementsCountAround,
    elementsCountAlong, elementsCountThroughWall, transitElementList):
    """
    Generates coordinates from inner to outer surface using coordinates
    and derivatives of inner surface.
    :param xInner: Coordinates on inner surface
    :param d1Inner: Derivatives on inner surface around tube
    :param d2Inner: Derivatives on inner surface along tube
    :param d3Inner: Derivatives on inner surface through wall
    :param wallThicknessList: Wall thickness for each element along tube
    :param relativeThicknessList: Relative wall thickness for each element through wall
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

    if relativeThicknessList:
        xi3 = 0.0
        xi3List = [0.0]
        for n3 in range(elementsCountThroughWall):
            xi3 += relativeThicknessList[n3]
            xi3List.append(xi3)
        relativeThicknessList.append(relativeThicknessList[-1])

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
                curvature = interp.getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[n + elementsCountAround],
                                                            d2Inner[n + elementsCountAround],
                                                            vector.normalise(d3Inner[n]), 0.0)
            elif n2 == elementsCountAlong:
                curvature = interp.getCubicHermiteCurvature(xInner[n - elementsCountAround],
                                                            d2Inner[n - elementsCountAround],
                                                            xInner[n], d2Inner[n], vector.normalise(d3Inner[n]), 1.0)
            else:
                curvature = 0.5*(
                    interp.getCubicHermiteCurvature(xInner[n - elementsCountAround], d2Inner[n - elementsCountAround],
                                                    xInner[n], d2Inner[n], vector.normalise(d3Inner[n]), 1.0) +
                    interp.getCubicHermiteCurvature(xInner[n], d2Inner[n],
                                                    xInner[n + elementsCountAround], d2Inner[n + elementsCountAround],
                                                    vector.normalise(d3Inner[n]), 0.0))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = xi3List[n3] if relativeThicknessList else 1.0/elementsCountThroughWall * n3
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
                d3 = [c * wallThickness * (relativeThicknessList[n3] if relativeThicknessList else 1.0/elementsCountThroughWall) for c in norm]
                d3List.append(d3)

    return xList, d1List, d2List, d3List, curvatureList

def createFlatCoordinates(xiList, lengthAroundList, totalLengthAlong, wallThickness, relativeThicknessList,
                          elementsCountAround, elementsCountAlong, elementsCountThroughWall, transitElementList):
    """
    Calculates flat coordinates for a tube when it is opened into a flat preparation.
    :param xiList: List containing xi for each point around the outer surface of the tube.
    :param lengthAroundList: List of total arclength around the outer surface for each element along.
    :param totalLengthAlong: Total length along tube.
    :param wallThickness: Thickness of wall.
    :param elementsCountAround: Number of elements around tube.
    :param elementsCountAlong: Number of elements along tube.
    :param elementsCountThroughWall: Number of elements through wall.
    :param transitElementList: stores true if element around is a
    transition element between a big and small element.
    :return: coordinates and derivatives of flat coordinates field.
    """
    if relativeThicknessList:
        xi3 = 0.0
        xi3List = [0.0]
        for n3 in range(elementsCountThroughWall):
            xi3 += relativeThicknessList[n3]
            xi3List.append(xi3)
        relativeThicknessList.append(relativeThicknessList[-1])

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
            z = wallThickness * (xi3List[n3] if relativeThicknessList else 1.0 / elementsCountThroughWall * n3)
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

    return xFlatList, d1FlatList, d2FlatList

def createOrganCoordinates(xiList, relativeThicknessList, lengthToDiameterRatio, wallThicknessToDiameterRatio,
                           elementsCountAround, elementsCountAlong, elementsCountThroughWall, transitElementList):
    """
    Calculates organ coordinates and derivatives represented by a cylindrical tube with a unit inner diameter,
    length equivalent to lengthToDiameterRatio and wall thickness of wallThicknessToDiameterRatio.
    :param xiList: List containing xi for each point around the outer surface of the tube.
    :param relativeThicknessList: Relative thickness of each element through wall for organ coordinates.
    :param lengthToDiameterRatio: Ratio of total length along organ to inner diameter of organ
    :param wallThicknessToDiameterRatio: Ratio of wall thickness to inner diameter of organ.
    :param elementsCountAround: Number of elements around tube.
    :param elementsCountAlong: Number of elements along tube.
    :param elementsCountThroughWall: Number of elements through wall.
    :param transitElementList: stores true if element around is a transition element between a big and small element.
    :return: coordinates and derivatives of organ coordinates field.
    """
    if relativeThicknessList:
        xi3 = 0.0
        xi3List = [0.0]
        for n3 in range(elementsCountThroughWall):
            xi3 += relativeThicknessList[n3]
            xi3List.append(xi3)
        relativeThicknessList.append(relativeThicknessList[-1])

    # Calculate organ coordinates and derivatives
    xOrganList = []
    d1OrganList = []
    d2OrganList = []
    d2 = [0.0, lengthToDiameterRatio / elementsCountAlong, 0.0]

    for n2 in range(elementsCountAlong + 1):
        cx = [0.0, lengthToDiameterRatio / elementsCountAlong * n2, 0.0]
        xiFace = xiList[n2]
        for n3 in range(elementsCountThroughWall + 1):
            xi3 = xi3List[n3] if relativeThicknessList else 1.0/elementsCountThroughWall * n3
            radius = 0.5 + wallThicknessToDiameterRatio * xi3
            axis1 = [0.0, 0.0, radius]
            axis2 = [radius, 0.0, 0.0]
            d1List = []
            for n1 in range(len(xiFace) - 1):
                dTheta = (xiFace[n1 + 1 if n1 < len(xiFace) - 1 else 0] - xiFace[n1]) * 2.0 * math.pi
                radiansAround = 2.0 * math.pi * xiFace[n1]
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                xOrganList.append([(cx[c] + cosRadiansAround * axis1[c] + sinRadiansAround * axis2[c]) for c in range(3)])
                d1List.append([ dTheta*(-sinRadiansAround*axis1[c] + cosRadiansAround*axis2[c]) for c in range(3) ])
                d2OrganList.append(d2)

            # To modify derivative along transition elements
            for i in range(len(transitElementList)):
                if transitElementList[i]:
                    d1List[i] = vector.setMagnitude(d1List[i], vector.magnitude(d1List[i - 1]))
                    d1List[i + 1] = vector.setMagnitude(d1List[i+ 1], vector.magnitude(d1List[(i + 2) % elementsCountAround]))

            d1OrganList += d1List

    return xOrganList, d1OrganList, d2OrganList

def createNodesAndElements(region,
    x, d1, d2, d3,
    xFlat, d1Flat, d2Flat,
    xOrgan, d1Organ, d2Organ, organCoordinateFieldName,
    elementsCountAround, elementsCountAlong, elementsCountThroughWall,
    annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
    firstNodeIdentifier, firstElementIdentifier,
    useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd):
    """
    Create nodes and elements for the coordinates and flat coordinates fields.
    :param x, d1, d2, d3: coordinates and derivatives of coordinates field.
    :param xFlat, d1Flat, d2Flat, d3Flat: coordinates and derivatives of
    flat coordinates field.
    :param xOrgan, d1Organ, d2Organ, d3Organ, organCoordinateFieldName: coordinates, derivatives and name of organ
    coordinates field.
    :param elementsCountAround: Number of elements around tube.
    :param elementsCountAlong: Number of elements along tube.
    :param elementsCountThroughWall: Number of elements through wall.
    :param annotationGroupsAround: Annotation groups of elements around.
    :param annotationGroupsAlong: Annotation groups of elements along.
    :param annotationGroupsThroughWall: Annotation groups of elements through wall.
    :param firstNodeIdentifier, firstElementIdentifier: first node and
    element identifier to use.
    :param useCubicHermiteThroughWall: use linear when false
    :param useCrossDerivatives: use cross derivatives when true
    :return nodeIdentifier, elementIdentifier, allAnnotationGroups
    """

    nodeIdentifier = firstNodeIdentifier
    elementIdentifier = firstElementIdentifier
    zero = [ 0.0, 0.0, 0.0 ]

    fm = region.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()

    # Coordinates field
    coordinates = findOrCreateFieldCoordinates(fm)
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

    elementtemplateApex = mesh.createElementtemplate()
    elementtemplateApex.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    if xFlat:
        # Flat coordinates field
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eftFlat1 = bicubichermitelinear.createEftBasic()
        eftFlat2 = bicubichermitelinear.createEftOpenTube()

        flatCoordinates = findOrCreateFieldCoordinates(fm, name="flat coordinates")
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

        flatNodetemplate3 = nodes.createNodetemplate()
        flatNodetemplate3.defineField(flatCoordinates)
        flatNodetemplate3.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        flatNodetemplate3.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
        flatNodetemplate3.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            flatNodetemplate3.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

        flatElementtemplate1 = mesh.createElementtemplate()
        flatElementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        flatElementtemplate1.defineField(flatCoordinates, -1, eftFlat1)

        flatElementtemplate2 = mesh.createElementtemplate()
        flatElementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        flatElementtemplate2.defineField(flatCoordinates, -1, eftFlat2)

        if closedProximalEnd:
            flatElementtemplateApex1 = mesh.createElementtemplate()
            flatElementtemplateApex1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

            flatElementtemplateApex2 = mesh.createElementtemplate()
            flatElementtemplateApex2.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    if xOrgan:
        # Organ coordinates field
        bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eftOrgan = bicubichermitelinear.createEftBasic()

        organCoordinates = findOrCreateFieldCoordinates(fm, name=organCoordinateFieldName)
        organNodetemplate = nodes.createNodetemplate()
        organNodetemplate.defineField(organCoordinates)
        organNodetemplate.setValueNumberOfVersions(organCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        organNodetemplate.setValueNumberOfVersions(organCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        organNodetemplate.setValueNumberOfVersions(organCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            organNodetemplate.setValueNumberOfVersions(organCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

        organElementtemplate = mesh.createElementtemplate()
        organElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        organElementtemplate.defineField(organCoordinates, -1, eftOrgan)

        organElementtemplateApex = mesh.createElementtemplate()
        organElementtemplateApex.setElementShapeType(Element.SHAPE_TYPE_CUBE)

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
        # print('NodeIdentifier = ', nodeIdentifier, x[n], d1[n], d2[n])
        nodeIdentifier = nodeIdentifier + 1

    # Flat coordinates field
    if xFlat:
        nodeIdentifier = firstNodeIdentifier
        if closedProximalEnd:
            # Apexes
            for n in range(elementsCountThroughWall + 1):
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                node.merge(flatNodetemplate3)
                cache.setNode(node)
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xFlat[n])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Flat[n])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1Flat[n])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Flat[n])
                if useCrossDerivatives:
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                nodeIdentifier = nodeIdentifier + 1
        for n2 in range(elementsCountAlong if closedProximalEnd else elementsCountAlong + 1):
            for n3 in range(elementsCountThroughWall + 1):
                for n1 in range(elementsCountAround):
                    if closedProximalEnd:
                        i = n2*(elementsCountAround + 1)*(elementsCountThroughWall + 1) + (elementsCountAround + 1)*n3 + n1 + elementsCountThroughWall + 1
                    else:
                        i = n2*(elementsCountAround + 1)*(elementsCountThroughWall + 1) + (elementsCountAround + 1)*n3 + n1
                    node = nodes.findNodeByIdentifier(nodeIdentifier)
                    node.merge(flatNodetemplate2 if n1 == 0 else flatNodetemplate1)
                    cache.setNode(node)
                    # print('NodeIdentifier', nodeIdentifier, 'version 1, xList Index =', i+1)
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xFlat[i])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Flat[i])
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Flat[i])
                    if useCrossDerivatives:
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    if n1 == 0:
                        # print('NodeIdentifier', nodeIdentifier, 'version 2, xList Index =', i+elementsCountAround+1)
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2, xFlat[i+elementsCountAround])
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, d1Flat[i+elementsCountAround])
                        flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2, d2Flat[i+elementsCountAround])
                        if useCrossDerivatives:
                            flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
                    nodeIdentifier = nodeIdentifier + 1

    # Organ coordinates field
    if xOrgan:
        nodeIdentifier = firstNodeIdentifier
        for n in range(len(xOrgan)):
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            node.merge(organNodetemplate)
            cache.setNode(node)
            organCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOrgan[n])
            organCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Organ[n])
            organCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Organ[n])
            if useCrossDerivatives:
                organCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

    # create elements
    radiansPerElementAround = math.pi*2.0 / elementsCountAround
    radiansPerElementAroundFlat = math.pi / elementsCountAround
    allValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                      Node.VALUE_LABEL_D2_DS1DS2]

    allAnnotationGroups = []

    if closedProximalEnd:
        # Create apex elements
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                va = e1
                vb = (e1 + 1) % elementsCountAround
                eftApex = eftfactory.createEftShellPoleBottom(va * 100, vb * 100)
                elementtemplateApex.defineField(coordinates, -1, eftApex)
                element = mesh.createElement(elementIdentifier, elementtemplateApex)
                bni1 = e3 + 1
                bni2 = elementsCountThroughWall + 1 + elementsCountAround*e3 + e1 + 1
                bni3 = elementsCountThroughWall + 1 + elementsCountAround*e3 + (e1 + 1) % elementsCountAround + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni1 + 1, bni2 + elementsCountAround, bni3 + elementsCountAround]
                element.setNodesByIdentifier(eftApex, nodeIdentifiers)
                onOpening = e1 > elementsCountAround - 2
                if xOrgan:
                    vao = e1
                    vbo = (e1 + 1) % elementsCountAround
                    eftApexOrgan = eftfactory.createEftShellPoleBottom(vao * 100, vbo * 100)
                    organElementtemplateApex.defineField(organCoordinates, -1, eftApexOrgan)
                    element.merge(organElementtemplateApex)
                    element.setNodesByIdentifier(eftApexOrgan, nodeIdentifiers)
                if xFlat:
                    vaf = e1 + elementsCountAround
                    vbf = vaf + 1
                    eftApexFlat = eftfactory.createEftShellPoleBottom(vaf * 100, vbf * 100)
                    if e1 >= (elementsCountAround // 2):
                        remapEftNodeValueLabelsVersion(eftApexFlat, [1, 4], [Node.VALUE_LABEL_D_DS1], 2)
                    flatElementtemplateApex1.defineField(flatCoordinates, -1, eftApexFlat)
                    element.merge(flatElementtemplateApex1)
                    element.setNodesByIdentifier(eftApexFlat, nodeIdentifiers)
                    # set general linear map coefficients
                    e1b = e1 - (elementsCountAround // 2)
                    radiansAroundApexFlat = e1b * radiansPerElementAroundFlat
                    radiansAroundApexFlatNext = (e1b + 1) * radiansPerElementAroundFlat
                    scalefactorsFlat = [
                        -1.0,
                        math.sin(radiansAroundApexFlat), math.cos(radiansAroundApexFlat), radiansPerElementAroundFlat,
                        math.sin(radiansAroundApexFlatNext), math.cos(radiansAroundApexFlatNext), radiansPerElementAroundFlat,
                        math.sin(radiansAroundApexFlat), math.cos(radiansAroundApexFlat), radiansPerElementAroundFlat,
                        math.sin(radiansAroundApexFlatNext), math.cos(radiansAroundApexFlatNext), radiansPerElementAroundFlat
                    ]
                    result = element.setScaleFactors(eftApexFlat, scalefactorsFlat)
                    if onOpening:
                        lnRemapV2 = [3, 6]
                        eftApexOpen = eftApexFlat
                        remapEftNodeValueLabelsVersion(eftApexOpen, lnRemapV2, allValueLabels, 2)
                        flatElementtemplateApex2.defineField(flatCoordinates, -1, eftApexOpen)
                        element.merge(flatElementtemplateApex2)
                        element.setNodesByIdentifier(eftApexOpen, nodeIdentifiers)
                        if eftApexOpen.getNumberOfLocalScaleFactors() > 0:
                            element.setScaleFactor(eftApexOpen, 1, -1.0)
                # set general linear map coefficients
                radiansAround = e1 * radiansPerElementAround
                radiansAroundNext = ((e1 + 1) % elementsCountAround) * radiansPerElementAround
                scalefactors = [
                    -1.0,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                    math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                    math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
                ]
                result = element.setScaleFactors(eftApex, scalefactors)
                elementIdentifier = elementIdentifier + 1
                annotationGroups = annotationGroupsAround[e1] + annotationGroupsAlong[0] + \
                                   annotationGroupsThroughWall[e3]
                if annotationGroups:
                    allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                    for annotationGroup in annotationGroups:
                        meshGroup = annotationGroup.getMeshGroup(mesh)
                        meshGroup.addElement(element)

    # Create regular elements
    now = elementsCountAround * (elementsCountThroughWall + 1)
    for e2 in range(1 if closedProximalEnd else 0, elementsCountAlong):
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                if closedProximalEnd:
                    bni11 = (e2-1) * now + e3 * elementsCountAround + e1 + 1 + (elementsCountThroughWall + 1)
                    bni12 = (e2-1) * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + 1 + \
                            (elementsCountThroughWall + 1)
                    bni21 = (e2-1) * now + (e3 + 1) * elementsCountAround + e1 + 1 + (elementsCountThroughWall + 1)
                    bni22 = (e2-1) * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + 1 + \
                            (elementsCountThroughWall + 1)
                else:
                    bni11 = e2 * now + e3 * elementsCountAround + e1 + 1
                    bni12 = e2 * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                    bni21 = e2 * now + (e3 + 1) * elementsCountAround + e1 + 1
                    bni22 = e2 * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + 1
                nodeIdentifiers = [bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now]
                onOpening = e1 > elementsCountAround - 2
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, nodeIdentifiers)
                if xFlat:
                    element.merge(flatElementtemplate2 if onOpening else flatElementtemplate1)
                    element.setNodesByIdentifier(eftFlat2 if onOpening else eftFlat1, nodeIdentifiers)
                if xOrgan:
                    element.merge(organElementtemplate)
                    element.setNodesByIdentifier(eftOrgan, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

                annotationGroups = annotationGroupsAround[e1] + annotationGroupsAlong[e2] + \
                                   annotationGroupsThroughWall[e3]
                if annotationGroups:
                    allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                    for annotationGroup in annotationGroups:
                        meshGroup = annotationGroup.getMeshGroup(mesh)
                        meshGroup.addElement(element)

    fm.endChange()

    return nodeIdentifier, elementIdentifier, allAnnotationGroups

class CylindricalSegmentTubeMeshInnerPoints:
    """
    Generates inner profile of a cylindrical segment for use by tubemesh.
    """

    def __init__(self, elementsCountAround, elementsCountAlongSegment,
                 segmentLength, wallThickness, innerRadiusList, startPhase):

        self._elementsCountAround = elementsCountAround
        self._elementsCountAlongSegment = elementsCountAlongSegment
        self._segmentLength = segmentLength
        self._wallThickness = wallThickness
        self._innerRadiusList = innerRadiusList
        self._xiList = []
        self._flatWidthList = []
        self._startPhase = startPhase

    def getCylindricalSegmentTubeMeshInnerPoints(self, nSegment):

        elementAlongStartIdx = nSegment * self._elementsCountAlongSegment
        elementAlongEndIdx = (nSegment + 1) * self._elementsCountAlongSegment

        xInner, d1Inner, d2Inner, transitElementList, xiSegment, flatWidthSegment, segmentAxis, radiusAlongSegmentList \
            = getCylindricalSegmentInnerPoints(self._elementsCountAround, self._elementsCountAlongSegment,
                                               self._segmentLength, self._wallThickness,
                                               self._innerRadiusList[elementAlongStartIdx: elementAlongEndIdx + 1],
                                               self._startPhase)

        startIdx = 0 if nSegment == 0 else 1
        xi = xiSegment[startIdx:self._elementsCountAlongSegment + 1]
        self._xiList += xi

        flatWidth = flatWidthSegment[startIdx:self._elementsCountAlongSegment + 1]
        self._flatWidthList += flatWidth

        return xInner, d1Inner, d2Inner, transitElementList, segmentAxis, radiusAlongSegmentList

    def getFlatWidthAndXiList(self):
        return self._flatWidthList, self._xiList

def getCylindricalSegmentInnerPoints(elementsCountAround, elementsCountAlongSegment, segmentLength,
                                     wallThickness, radiusList, startPhase):
    """
    Generates a 3-D cylindrical segment mesh with variable numbers of elements
    around, along the central path, and through wall.
    :param elementsCountAround: Number of elements around.
    :param elementsCountAlongSegment: Number of elements along cylindrical segment.
    :param segmentLength: Length of a cylindrical segment.
    :param wallThickness: Thickness of wall.
    :param radiusList: Inner radius at elements along tube length.
    :param startPhase: Phase at start.
    :return coordinates, derivatives on inner surface of a cylindrical segment.
    :return transitElementList: stores true if element around is an element that
    transits between a big and small element.
    :return xiList: List of xi for each node around. xi refers to node position
    along the width when cylindrical segment is opened into a flat preparation,
    nominally in [0.0, 1.0].
    :return flatWidthList: List of width around elements for each element
    along cylindrical segment when the segment is opened into a flat preparation.
    :return segmentAxis: Axis of segment.
    :return sRadiusAlongSegment: radius of each element along segment.
    """

    transitElementList = [0] * elementsCountAround

    # create nodes
    segmentAxis = [0.0, 0.0, 1.0]

    xFinal = []
    d1Final = []
    d2Final = []
    xiList = []
    flatWidthList = []
    sRadiusAlongSegment = []

    for n2 in range(elementsCountAlongSegment + 1):
        radius = radiusList[n2]
        sRadiusAlongSegment.append(radius)
        z = segmentLength / elementsCountAlongSegment * n2 + startPhase / 360.0 * segmentLength

        xLoop, d1Loop = createCirclePoints([0.0, 0.0, z], [radius, 0.0, 0.0], [0.0, radius, 0.0],
                                           elementsCountAround, startRadians=0.0)
        xFinal = xFinal + xLoop
        d1Final = d1Final + d1Loop

    # Smooth d2 for segment
    smoothd2Raw = []
    for n1 in range(elementsCountAround):
        nx = []
        nd2 = []
        for n2 in range(elementsCountAlongSegment + 1):
            n = n2 * elementsCountAround + n1
            nx.append(xFinal[n])
            nd2.append(segmentAxis)
        smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2)
        smoothd2Raw.append(smoothd2)

    # Re-arrange smoothd2
    for n2 in range(elementsCountAlongSegment + 1):
        radius = sRadiusAlongSegment[n2]
        flatWidth = 2.0*math.pi*(radius + wallThickness)
        flatWidthList.append(flatWidth)
        xiFace = []
        for n1 in range(elementsCountAround):
            d2Final.append(smoothd2Raw[n1][n2])
        for n1 in range(elementsCountAround + 1):
            xi = 1.0/elementsCountAround * n1
            xiFace.append(xi)
        xiList.append(xiFace)

    return xFinal, d1Final, d2Final, transitElementList, xiList, flatWidthList, segmentAxis, sRadiusAlongSegment
