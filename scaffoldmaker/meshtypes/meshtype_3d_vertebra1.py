"""
Generates a 3-D vertebra model including body, spinous process, transverse processes,
articular processes, vertebral arch and pedicle.
"""

from __future__ import division

import math

from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.zinc_utils import *
from scaffoldmaker.utils.geometry import *
from scaffoldmaker.utils.interpolation import *
from scaffoldmaker.utils.vector import *
from scaffoldmaker.utils.transformation import *

from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node


class MeshType_3d_vertebra1(Scaffold_base):
    """
    Generates a 3-D vertebral mesh
    """

    @staticmethod
    def getName():
        return '3D Vertebra 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {'Number of elements around': 7, 'Number of elements up': 3,  # 'Number of elements through wall': 2,
                # 'Major diameter': 1.1,
                # 'Minor diameter': 0.9,
                'Height': 0.3, 'Body thickness ratio': 0.85, 'Body posterior surface curvature factor': 0.5,
                'Body posterior surface depth factor': 1.1, 'Pedicle arch factor': 0.3, 'Pedicle length factor': 1.2,
                'Arch peak length factor': 2., 'Arch length factor': 1.3, 'Use cross derivatives': False,
                'Refine': False, 'Refine number of elements around': 1, 'Refine number of elements through wall': 1,
                'Refine number of elements up': 1}

    @staticmethod
    def getOrderedOptionNames():
        return ['Number of elements around', 'Number of elements up',  # 'Number of elements through wall',
                # 'Major diameter',
                # 'Minor diameter',
                'Height', 'Body thickness ratio', 'Body posterior surface curvature factor',
                'Body posterior surface depth factor', 'Pedicle arch factor', 'Pedicle length factor',
                'Arch peak length factor', 'Arch length factor', 'Refine', 'Use cross derivatives',
                'Refine number of elements around', 'Refine number of elements up',
                'Refine number of elements through wall']

    @staticmethod
    def checkOptions(options):
        for key in ['Refine number of elements around', 'Refine number of elements up',
                    'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements up'] < 2:
            options['Number of elements up'] = 2
        elif options['Number of elements up'] > 9:
            options['Number of elements up'] = 9
        if options['Number of elements around'] < 4:
            options['Number of elements around'] = 4
        # if options['Number of elements through wall'] < 2:
        #     options['Number of elements through wall'] = 2
        # elif options['Number of elements through wall'] > 4:
        #     options['Number of elements through wall'] = 4
        # if options['Major diameter'] < 0.0:
        #     options['Major diameter'] = 0.0
        # if options['Minor diameter'] < 0.0:
        #     options['Minor diameter'] = 0.0
        if options['Height'] < 0.0:
            options['Height'] = 0.0
        if options['Body thickness ratio'] <= 0.5:
            options['Body thickness ratio'] = 0.55
        elif options['Body thickness ratio'] >= 1.0:
            options['Body thickness ratio'] = 0.9
        if options['Body posterior surface curvature factor'] < 0.5:
            options['Body posterior surface curvature factor'] = 0.5
        elif options['Body posterior surface curvature factor'] > 1.5:
            options['Body posterior surface curvature factor'] = 1.5
        if options['Body posterior surface depth factor'] < 1.0:
            options['Body posterior surface depth factor'] = 1.1
        elif options['Body posterior surface depth factor'] > 1.8:
            options['Body posterior surface depth factor'] = 1.8
            options['Body posterior surface curvature factor'] = 0.7
        if options['Pedicle arch factor'] < 0.0:
            options['Pedicle arch factor'] = 0.0
        elif options['Pedicle arch factor'] > 1.0:
            options['Pedicle arch factor'] = 1.0
        if options['Pedicle length factor'] < 1.0:
            options['Pedicle length factor'] = 1.0
        elif options['Pedicle length factor'] > 1.3:
            options['Pedicle length factor'] = 1.3
        if options['Arch peak length factor'] < 1.0:
            options['Arch peak length factor'] = 1.1
        if options['Arch length factor'] < 1.0:
            options['Arch length factor'] = 1.1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generates the base tricubic Hermite mesh for this MeshType. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAround = options['Number of elements around']
        elementsCountUp = options['Number of elements up']
        # elementsCountThroughWall = options['Number of elements through wall']
        elementsCountThroughWall = 2
        useCrossDerivatives = options['Use cross derivatives']
        # axisA = options['Major diameter']
        # axisB = options['Minor diameter']
        axisA = 1.3
        axisB = 0.9
        height = options['Height']
        useCubicHermiteThroughWall = True
        foldFactor = options['Body posterior surface curvature factor']
        depthFactor = options['Body posterior surface depth factor']
        pedicleArchFactor = options['Pedicle arch factor']
        pedicleLenghtFactor = options['Pedicle length factor']
        vertebralArchLengthFactor = options['Arch length factor']
        vertebralArchPeakLengthFactor = options['Arch peak length factor']

        zero = [0.0, 0.0, 0.0]

        """ Get nodes along the centre line """
        # Straight tube
        cx = [zero, [0.0, 0.0, height]]
        cd1 = [[0.0, 0.0, 0.1], [0.0, 0.0, 0.1]]
        cd2 = [[0.0, 0.3, 0.0], [0.0, 0.3, 0.0]]
        cd3 = [[0.35, 0.0, 0.0], [0.35, 0.0, 0.0]]
        # thickness in cd2 and cd3 directions and derivatives
        t2 = [1.0 - options['Body thickness ratio']] * 2
        t3 = [1.0 - options['Body thickness ratio']] * 2

        sx, sd1, se, sxi, _ = sampleCubicHermiteCurves(cx, cd1, elementsCountUp)
        sd2 = interpolateSampleLinear(cd2, se, sxi)
        sd3 = interpolateSampleLinear(cd3, se, sxi)
        st2 = interpolateSampleLinear(t2, se, sxi)
        st3 = interpolateSampleLinear(t3, se, sxi)

        # Find unit normals and binormals at each sample points
        sNormal = list()
        sBinormal = list()

        # Normal and binormal
        prevUnitTangent = normalise(sd1[0])
        if magnitude(crossproduct3(prevUnitTangent, [0.0, 0.0, 1.0])) > 0.0:
            prevBinormal = crossproduct3(prevUnitTangent, [0.0, 0.0, 1.0])
        else:
            prevBinormal = crossproduct3(prevUnitTangent, [0.0, -1.0, 0.0])
        prevUnitBinormal = normalise(prevBinormal)
        prevUnitNormal = crossproduct3(prevUnitBinormal, prevUnitTangent)
        sNormal.append(prevUnitNormal)
        sBinormal.append(prevUnitBinormal)
        for n in range(1, elementsCountUp + 1):
            unitTangent = normalise(sd1[n])
            cp = crossproduct3(prevUnitTangent, unitTangent)
            if magnitude(cp) > 0.0:
                axisRot = normalise(cp)
                thetaRot = math.acos(dotproduct(prevUnitTangent, unitTangent))
                rotFrame = rotationMatrix(thetaRot, axisRot)
                rotNormal = [rotFrame[j][0] * prevUnitNormal[0] + rotFrame[j][1] * prevUnitNormal[1] + rotFrame[j][2] *
                             prevUnitNormal[2] for j in range(3)]
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
        coordinates = getOrCreateCoordinateField(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateApex = nodes.createNodetemplate()
        nodetemplateApex.defineField(coordinates)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
        else:
            nodetemplate = nodetemplateApex

        cache = fm.createFieldcache()

        """ Create nodes """
        nodeIdentifier = 1

        nodesForLeftPedicleElement = list()
        nodesForRightPedicleElement = list()

        nodesForVertebralArchElement = list()

        nodesOfTheVertebralArchCetreAxis = list()

        # Create bottom node
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, zero)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1,
                                      [0.0, axisB / elementsCountThroughWall, 0.0])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, 0.0, height / elementsCountUp])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                      [axisA / elementsCountThroughWall, 0.0, 0.0])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

        # Create node along the centre line between bottom and top nodes
        for n2 in range(1, elementsCountUp):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [0, 0, n2 * height / elementsCountUp])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1,
                                          [0.0, axisB / elementsCountThroughWall, 0.0])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, 0.0, height / elementsCountUp])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                          [axisA / elementsCountThroughWall, 0.0, 0.0])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        # Create top node
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [0, 0, height])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, 0.0, height / elementsCountUp])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1,
                                      [0.0, axisB / elementsCountThroughWall, 0.0])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                      [axisA / elementsCountThroughWall, 0.0, 0.0])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

        for n3 in range(elementsCountThroughWall):
            for n2 in range(elementsCountUp + 1):
                aThroughWallElement = sd2[n2][1] + st2[n2] * (n3 / elementsCountThroughWall)
                bThroughWallElement = sd3[n2][0] + st3[n2] * (n3 / elementsCountThroughWall)
                perimeterAroundWallElement = getApproximateEllipsePerimeter(aThroughWallElement, bThroughWallElement)
                arcLengthPerElementAround = perimeterAroundWallElement / elementsCountAround
                prevRadiansAround = updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement, 0.0,
                                                                  arcLengthPerElementAround)
                st2PerWallElement = st2[n2] / elementsCountThroughWall
                st3PerWallElement = st3[n2] / elementsCountThroughWall

                if n2 < elementsCountUp:
                    aThroughWallElementNext = sd2[n2 + 1][1] + st2[n2 + 1] * (n3 / elementsCountThroughWall)
                    bThroughWallElementNext = sd3[n2 + 1][0] + st3[n2 + 1] * (n3 / elementsCountThroughWall)
                    perimeterAroundWallElementNext = getApproximateEllipsePerimeter(aThroughWallElementNext,
                                                                                    bThroughWallElementNext)
                    arcLengthPerElementAroundNext = perimeterAroundWallElementNext / elementsCountAround

                for n1 in range(elementsCountAround):
                    arcLengthAround = n1 * arcLengthPerElementAround
                    radiansAround = -1 * updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement, 0.0,
                                                                       arcLengthAround)
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [sx[n2][j] + aThroughWallElement * cosRadiansAround * sBinormal[n2][
                        j] + bThroughWallElement * sinRadiansAround * sNormal[n2][j] for j in range(3)]
                    dx_ds1 = [(radiansAround - prevRadiansAround) * (
                            aThroughWallElement * -sinRadiansAround * sBinormal[n2][
                        j] + bThroughWallElement * cosRadiansAround * sNormal[n2][j]) for j in range(3)]

                    # Calculate curvature to find d1 for node
                    unitNormal = normalise([aThroughWallElement * cosRadiansAround * sBinormal[n2][
                        j] + bThroughWallElement * sinRadiansAround * sNormal[n2][j] for j in range(3)])

                    if n2 == 0:
                        curvature = getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2 + 1], sd1[n2 + 1], unitNormal, 0.0)
                    elif n2 == elementsCountUp:
                        curvature = getCubicHermiteCurvature(sx[n2 - 1], sd1[n2 - 1], sx[n2], sd1[n2], unitNormal, 1.0)
                    else:
                        curvature = 0.5 * (
                                getCubicHermiteCurvature(sx[n2 - 1], sd1[n2 - 1], sx[n2], sd1[n2], unitNormal,
                                                         1.0) + getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2 + 1],
                                                                                         sd1[n2 + 1], unitNormal, 0.0))
                    wallDistance = magnitude([aThroughWallElement * cosRadiansAround * sBinormal[n2][
                        j] + bThroughWallElement * sinRadiansAround * sNormal[n2][j] for j in range(3)])
                    factor = 1.0 - curvature * wallDistance
                    d1Wall = [factor * c for c in sd1[n2]]

                    # Calculate curvature to find d1 for downstream node
                    if n2 < elementsCountUp:
                        arcLengthAroundNext = n1 * arcLengthPerElementAroundNext
                        radiansAroundNext = -1 * updateEllipseAngleByArcLength(aThroughWallElementNext,
                                                                               bThroughWallElementNext, 0.0,
                                                                               arcLengthAroundNext)
                        cosRadiansAroundNext = math.cos(radiansAroundNext)
                        sinRadiansAroundNext = math.sin(radiansAroundNext)
                        xNext = [sx[n2 + 1][j] + aThroughWallElementNext * cosRadiansAroundNext * sBinormal[n2 + 1][
                            j] + bThroughWallElementNext * sinRadiansAroundNext * sNormal[n2 + 1][j] for j in range(3)]
                        unitNormalNext = normalise([aThroughWallElementNext * cosRadiansAroundNext * sBinormal[n2 + 1][
                            j] + bThroughWallElementNext * sinRadiansAroundNext * sNormal[n2 + 1][j] for j in range(3)])
                        if n2 + 1 == elementsCountUp:
                            curvatureNext = getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2 + 1], sd1[n2 + 1],
                                                                     unitNormalNext, 1.0)
                        else:
                            curvatureNext = 0.5 * (
                                    getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2 + 1], sd1[n2 + 1], unitNormalNext,
                                                             1.0) + getCubicHermiteCurvature(sx[n2 + 1], sd1[n2 + 1],
                                                                                             sx[n2 + 2], sd1[n2 + 2],
                                                                                             unitNormalNext, 0.0))
                        wallDistanceNext = magnitude([aThroughWallElementNext * cosRadiansAroundNext *
                                                      sBinormal[n2 + 1][
                                                          j] + bThroughWallElementNext * sinRadiansAroundNext *
                                                      sNormal[n2 + 1][j] for j in range(3)])
                        factorNext = 1.0 - curvatureNext * wallDistanceNext
                        d1WallNext = [factorNext * c for c in sd1[n2 + 1]]
                        arcLength = computeCubicHermiteArcLength(x, d1Wall, xNext, d1WallNext, True)
                        dx_ds2 = [arcLength * c for c in normalise(d1Wall)]
                        if n2 == elementsCountUp - 1:
                            secondLastX = x
                            secondLastd1Wall = d1Wall
                            lastX = xNext
                            lastd1Wall = d1WallNext
                    else:
                        arcLength = computeCubicHermiteArcLength(secondLastX, secondLastd1Wall, lastX, lastd1Wall, True)
                        dx_ds2 = [arcLength * c for c in normalise(lastd1Wall)]

                    dx_ds3 = [
                        st2PerWallElement * cosRadiansAround * sBinormal[n2][j] + st3PerWallElement * sinRadiansAround *
                        sNormal[n2][j] for j in range(3)]  # Modify later to calculate with interpolation

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

        generalNodeIndentifier = nodeIdentifier

        # Collecting the posterior curve nodes
        iNodeList = list()
        jNodeList = list()
        kNodeList = list()

        diff = elementsCountAround - elementsCountUp

        for e3 in range(elementsCountThroughWall):
            if e3 == 0:
                for e2 in range(elementsCountUp + 1):
                    if e2 == 0:
                        for e1 in range(elementsCountAround):
                            i = elementsCountUp + 2
                            j = i + 1
                            k = (i + elementsCountAround) - 1
                        iNodeList.append(i)
                        jNodeList.append(j)
                        kNodeList.append(k)
                    else:
                        for e1 in range(elementsCountAround):
                            i = (e2 * elementsCountAround) + 2 + elementsCountUp
                            j = i + 1
                            k = (i + elementsCountAround) - 1
                        iNodeList.append(i)
                        jNodeList.append(j)
                        kNodeList.append(k)
            else:
                nx = e3 * ((elementsCountUp * elementsCountAround) + elementsCountAround)
                for e2 in range(elementsCountUp):
                    for e1 in range(elementsCountAround - diff + 1):
                        if e2 == 1:
                            i = (e1 * e2 * elementsCountAround + 2 + elementsCountUp) + nx
                            j = i + 1
                            k = (i + elementsCountAround) - 1
                            iNodeList.append(i)
                            jNodeList.append(j)
                            kNodeList.append(k)

        # Sampling arclength on the posterior curves
        samplePoints = 4

        for n1 in range(len(iNodeList)):
            dx3List = list()

            # first node
            cache.setNode(nodes.findNodeByIdentifier(kNodeList[n1]))
            result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, dx1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, dx2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            result, dx3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)

            v1 = x1
            d1 = dx1
            dx3List.append(dx3)

            # second node
            cache.setNode(nodes.findNodeByIdentifier(iNodeList[n1]))
            result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, dx1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, dx2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            result, dx3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)

            v2 = x1
            d2 = dx1
            dx3List.append(dx3)

            # third node
            cache.setNode(nodes.findNodeByIdentifier(jNodeList[n1]))
            result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
            result, dx1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
            result, dx2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
            result, dx3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)

            v3 = x1
            d3 = dx1
            dx3List.append(dx3)

            nx = [v1, v2, v3]
            nd1 = [d1, d2, d3]
            sx, sd, se, sxi, _ = sampleCubicHermiteCurves(nx, nd1, samplePoints)

            sxn, sdn = [sx[1], sx[3]], [sd[1], sd[3]]

            # move the i nodes slightly anteriorly
            v2[0] = v2[0] / depthFactor
            cache.setNode(nodes.findNodeByIdentifier(iNodeList[n1]))
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, v2)

            # estimate the dx3 for the newly generated nodes in sxn.
            dx31 = [sum(x) for x in zip(dx3List[0], dx3List[1])]
            dx311 = [x / 2. for x in dx31]
            dx32 = [sum(x) for x in zip(dx3List[1], dx3List[2])]
            dx321 = [x / 2. for x in dx32]
            dx3All = [dx311, dx321]

            # fold derivative dx1
            for i in range(len(sdn)):
                for j in range(len(sdn[i])):
                    sdn[i][j] = sdn[i][j] * foldFactor

                # sdnSmooth = smoothCubicHermiteDerivativesLoop(sxn[i], sdn[i])

                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sxn[i])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sdn[i])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx2)
                if useCubicHermiteThroughWall:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx3All[i])
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                    if useCubicHermiteThroughWall:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                nodeIdentifier = nodeIdentifier + 1

        """ Create elements """
        mesh = fm.findMeshByDimension(3)

        if useCubicHermiteThroughWall:
            eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        else:
            eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        # Regular element template
        elementtemplateRegular = mesh.createElementtemplate()
        elementtemplateRegular.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # Tetrahedron element template
        elementtemplateWedge = mesh.createElementtemplate()
        elementtemplateWedge.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        elementIdentifier = 1

        no2 = elementsCountAround
        no3 = elementsCountAround * (elementsCountUp - 1)
        rni = (1 + elementsCountUp) - no3 - no2 + 1  # regular node identifier

        for e3 in range(elementsCountThroughWall):
            # Create elements on bottom pole
            radiansIncline = math.pi * 0.5 * e3 / elementsCountThroughWall
            radiansInclineNext = math.pi * 0.5 * (e3 + 1) / elementsCountThroughWall

            if e3 == 0:  # Create wedge elements on the bottom pole
                for e2 in range(elementsCountUp):
                    for e1 in range(elementsCountAround + 2):
                        va = e1
                        vb = (e1 + 1) % elementsCountAround
                        eft2 = eftfactory.createEftWedgeRadial(va * 100, vb * 100)
                        elementtemplateWedge.defineField(coordinates, -1, eft2)
                        element = mesh.createElement(elementIdentifier, elementtemplateWedge)
                        if e2 == 0:
                            bni1 = elementsCountUp + 2
                            bni2 = elementsCountUp + 3
                            bni3 = elementsCountUp + elementsCountAround + 2
                            bni4 = elementsCountUp + elementsCountAround + 3
                            if e1 == 0:
                                bni21 = generalNodeIndentifier + 1
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1 + e1, bni21, bni3 + e1, bni21 + 2 + e1]
                            elif e1 == 1:
                                bni21 = generalNodeIndentifier + 1
                                nodeIdentifiers = [e2 + 1, 2, bni21, bni1 + e1, bni21 + 1 + e1, bni3 + e1]
                            elif e1 <= elementsCountAround - 1 and e1 != 0 and e1 != 1:
                                nodeIdentifiers = [e2 + 1, 2, bni1 + e1 - 1, bni2 + e1 - 1, bni3 + e1 - 1,
                                                   bni4 + e1 - 1]
                            elif e1 == elementsCountAround:
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1 + e1 - 1, generalNodeIndentifier, bni3 + e1 - 1,
                                                   generalNodeIndentifier + 2]
                            else:
                                nodeIdentifiers = [e2 + 1, e2 + 2, generalNodeIndentifier, bni1,
                                                   generalNodeIndentifier + 2, bni3]
                        else:
                            bni1 = (e2 * elementsCountAround) + 2 + elementsCountUp
                            bni2 = e2 * 2 + generalNodeIndentifier + 1
                            bni3 = bni1 + elementsCountAround
                            bni4 = e2 * 2 + generalNodeIndentifier + 3
                            if e1 == 0:
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1, bni2, bni3, bni4]
                            elif e1 == 1:
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni2, bni1 + 1, bni4, bni3 + 1]
                            elif e1 < elementsCountAround and e1 != 0 and e1 != 1:
                                bni1 = (e2 * elementsCountAround) + 2 + elementsCountUp + 1
                                bni2 = bni1 + 1
                                bni3 = bni1 + elementsCountAround
                                bni4 = bni3 + 1
                                nx = e1 - 2
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1 + nx, bni2 + nx, bni3 + nx, bni4 + nx]
                            elif e1 == elementsCountAround:
                                nx = e1 - 2
                                bni1 = (e2 * elementsCountAround) + 2 + elementsCountUp + 1
                                bni2 = generalNodeIndentifier + e2 * 2
                                bni3 = bni1 + elementsCountAround
                                bni4 = bni2 + 2
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1 + nx, bni2, bni3 + nx, bni4]
                            else:
                                # pass
                                bni1 = generalNodeIndentifier + e2 * 2
                                bni2 = elementsCountUp + 2 + elementsCountAround * e2
                                bni3 = bni1 + 2
                                bni4 = elementsCountUp + 2 + elementsCountAround * e2 + elementsCountAround
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1, bni2, bni3, bni4]
                        result1 = element.setNodesByIdentifier(eft2, nodeIdentifiers)
                        scalefactors = [-1.0, 1.0, 1.0, 1.0, 1.0]
                        result2 = element.setScaleFactors(eft2, scalefactors)
                        elementIdentifier = elementIdentifier + 1

            # Create regular elements
            elif e3 == 1:
                for e2 in range(elementsCountUp):
                    for e1 in range(elementsCountAround + 2):

                        elementtemplateRegular.defineField(coordinates, -1, eft)
                        element = mesh.createElement(elementIdentifier, elementtemplateRegular)

                        if e1 == 0:
                            bni1 = (e2 * elementsCountAround) + 2 + elementsCountUp
                            bni2 = e2 * 2 + generalNodeIndentifier + 1
                            bni3 = elementsCountUp + 2 + elementsCountAround * e2 + elementsCountAround
                            bni4 = bni2 + 2
                            bni5 = bni1 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni6 = (e2 * 2) + (generalNodeIndentifier + (elementsCountUp * 2 + 3))
                            bni7 = bni3 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni8 = bni6 + 2
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                            if e2 >= elementsCountUp // 2:
                                nodesOfTheVertebralArchCetreAxis.append(bni5)
                                nodesOfTheVertebralArchCetreAxis.append(bni7)
                        elif e1 == 1:
                            bni1 = e2 * 2 + generalNodeIndentifier + 1
                            bni2 = (e2 * elementsCountAround) + 2 + elementsCountUp + 1
                            bni3 = bni1 + 2
                            bni4 = elementsCountUp + 2 + elementsCountAround * e2 + elementsCountAround + 1
                            bni5 = (e2 * 2) + (generalNodeIndentifier + (elementsCountUp * 2 + 3))
                            bni6 = bni2 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni7 = bni5 + 2
                            bni8 = bni6 + elementsCountAround
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                            if e2 >= elementsCountUp // 2:
                                nodesForRightPedicleElement.append(bni5)
                                nodesForRightPedicleElement.append(bni6)
                            if e2 == elementsCountUp - 1:
                                nodesForRightPedicleElement.append(bni7)
                                nodesForRightPedicleElement.append(bni8)

                        elif e1 <= elementsCountAround - 1 and e1 != 0 and e1 != 1:
                            bni1 = e1 + (e2 * elementsCountAround) + 2 + elementsCountUp - 1
                            bni2 = bni1 + 1
                            bni3 = bni2 + (elementsCountAround - 1)
                            bni4 = bni3 + 1
                            bni5 = bni1 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni6 = bni5 + 1
                            bni7 = bni3 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni8 = bni7 + 1
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                        elif e1 == elementsCountAround:
                            bni1 = e1 + (e2 * elementsCountAround) + 2 + elementsCountUp - 1
                            bni2 = generalNodeIndentifier + e2 * 2
                            bni3 = bni1 + elementsCountAround
                            bni4 = bni2 + 2
                            bni5 = bni1 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni6 = (e2 * 2) + (generalNodeIndentifier + (elementsCountUp * 2 + 2))
                            bni7 = bni5 + elementsCountAround
                            bni8 = bni6 + 2
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                            if e2 >= elementsCountUp // 2:
                                nodesForLeftPedicleElement.append(bni5)
                                nodesForLeftPedicleElement.append(bni6)
                            if e2 == elementsCountUp - 1:
                                nodesForLeftPedicleElement.append(bni7)
                                nodesForLeftPedicleElement.append(bni8)
                        else:
                            bni1 = generalNodeIndentifier + e2 * 2
                            bni2 = (e2 * elementsCountAround) + 2 + elementsCountUp
                            bni3 = bni1 + 2
                            bni4 = elementsCountUp + 2 + elementsCountAround * e2 + elementsCountAround
                            bni5 = (e2 * 2) + (generalNodeIndentifier + (elementsCountUp * 2 + 2))
                            bni6 = bni2 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni7 = bni5 + 2
                            bni8 = bni4 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                            if e2 >= elementsCountUp // 2:
                                nodesForVertebralArchElement.append(bni6)
                            if e2 == elementsCountUp - 1:
                                nodesForVertebralArchElement.append(bni8)
                        result1 = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1
            else:
                pass
                # for e2 in range(elementsCountUp):  #     for e1 in range(elementsCountAround + 2):  #  #
                # elementtemplateRegular.defineField(coordinates, -1, eft)  #         element = mesh.createElement(
                # elementIdentifier, elementtemplateRegular)  #  #         nx = (e3-1)*((elementsCountUp *
                # elementsCountAround) + elementsCountAround)  #         if e1 == 0:  #             bni1 = (e1 + (e2
                # * elementsCountAround) + 2 + elementsCountUp) + nx  #             bni2 = (e2 * 2) + (
                # generalNodeIndentifier + (elementsCountUp * 2 + 3))  #             bni3 = bni1 +
                # elementsCountAround  #             bni4 = bni2 + 2  #             bni5 = bni1 + ((elementsCountUp *
                #  elementsCountAround) + elementsCountAround)  #             bni6 = bni2 + elementsCountAround +
                # elementsCountUp  #             bni7 = bni3 + ((elementsCountUp * elementsCountAround) +
                # elementsCountAround)  #             bni8 = bni6 + 2  #             nodeIdentifiers = [bni1, bni2,
                # bni3, bni4, bni5, bni6, bni7, bni8]  #         elif e1 == 1:  # pass  # bni1 = (e1 + (e2 *
                # elementsCountAround) + 2 + elementsCountUp) + nx  # bni2 = e2 * 2 + (generalNodeIndentifier +
                # elementsCountAround + 4)  # bni3 = bni1 + elementsCountAround  # bni4 = bni2 + 2  # bni5 = bni1 + (
                # (elementsCountUp * elementsCountAround) + elementsCountAround)  # bni6 = (e2 * 2) + (
                # generalNodeIndentifier + (elementsCountUp * 2 + 3) + elementsCountAround  #                    +
                # elementsCountUp)  # bni7 = bni3 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                #  bni8 = bni6 + 2  # nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]  # else:  #
                #     pass  # bni1 = (e1 + (e2 * elementsCountAround) + 2 + elementsCountUp) + nx  # bni2 = bni1 - (
                # elementsCountAround - 1)  # bni3 = bni1 + elementsCountAround  # bni4 = bni3 - (elementsCountAround
                #  - 1)  # bni5 = bni1 + ((elementsCountUp * elementsCountAround) + elementsCountAround)  # bni6 =
                # bni5 - (elementsCountAround - 1)  # bni7 = bni3 + ((elementsCountUp * elementsCountAround) +
                # elementsCountAround)  # bni8 = bni7 - (elementsCountAround - 1)  # nodeIdentifiers = [bni1, bni2,
                # bni3, bni4, bni5, bni6, bni7, bni8]  # result1 = element.setNodesByIdentifier(eft, nodeIdentifiers)
                #  elementIdentifier = elementIdentifier + 1

        """ Create vertebral arch ring """
        _createVertebralArchCentreLineNodes(nodesOfTheVertebralArchCetreAxis, cache, nodes, coordinates,
                                            options['Body thickness ratio'], elementsCountAround, elementsCountUp,
                                            nodeIdentifier, nodetemplate, elementtemplateRegular, eft, mesh,
                                        elementIdentifier, nodesForLeftPedicleElement, nodesForRightPedicleElement)


        """ Create left and right pedicles """
        # # Left
        # leftPedicleNodes, nID, eID = _createPedicleElements(nodesForLeftPedicleElement, cache, nodes, coordinates,
        #                                                     nodeIdentifier, nodetemplate, useCubicHermiteThroughWall,
        #                                                     useCrossDerivatives, pedicleArchFactor, pedicleLenghtFactor,
        #                                                     elementsCountUp, elementtemplateRegular, eft, mesh,
        #                                                     elementIdentifier, zero)
        #
        # nodeIdentifier, elementIdentifier = nID, eID
        #
        # # Right
        # rightPedicleNodes, nID, eID = _createPedicleElements(nodesForRightPedicleElement, cache, nodes, coordinates,
        #                                                      nodeIdentifier, nodetemplate, useCubicHermiteThroughWall,
        #                                                      useCrossDerivatives, -pedicleArchFactor,
        #                                                      pedicleLenghtFactor, elementsCountUp,
        #                                                      elementtemplateRegular, eft, mesh, elementIdentifier, zero)
        #
        # nodeIdentifier, elementIdentifier = nID, eID
        #
        # """ Create the vertebral arch """
        # # Nodes of the peak of the arch
        # peakArchNodes, nID = _createVertebralArchNodes(nodesForVertebralArchElement, cache, nodes, coordinates,
        #                                                nodeIdentifier, nodetemplate, useCubicHermiteThroughWall,
        #                                                useCrossDerivatives, vertebralArchPeakLengthFactor, elongated=False,
        #                                                elongationNodes=None, elongationFactor=0)
        #
        # nodeIdentifier = nID
        #
        # # Nodes of the left side of the arch
        # leftArchNodes, nID = _createVertebralArchNodes(leftPedicleNodes, cache, nodes, coordinates, nodeIdentifier,
        #                                                nodetemplate, useCubicHermiteThroughWall, useCrossDerivatives,
        #                                                vertebralArchLengthFactor, elongated=True,
        #                                                elongationNodes=[leftPedicleNodes[0], leftPedicleNodes[2]],
        #                                                elongationFactor=1.4)
        #
        # # leftArchNodes, nID, eID = _createPedicleElements(leftPedicleNodes, cache, nodes, coordinates,
        # #                                                     nodeIdentifier, nodetemplate, useCubicHermiteThroughWall,
        # #                                                     useCrossDerivatives, pedicleArchFactor, pedicleLenghtFactor,
        # #                                                     elementsCountUp, elementtemplateRegular, eft, mesh,
        # #                                                     elementIdentifier, zero)
        # # nodeIdentifier, elementIdentifier = nID, eID
        #
        # nodeIdentifier = nID
        #
        # # Elements of the left side of the arch
        # eID = _createVertebralArchElements(leftPedicleNodes, leftArchNodes, coordinates, elementsCountUp,
        #                              elementtemplateRegular, eft, mesh, elementIdentifier)
        #
        # elementIdentifier = eID
        #
        # # Nodes of the right side of the arch
        # rightArchNodes, nID = _createVertebralArchNodes(rightPedicleNodes, cache, nodes, coordinates, nodeIdentifier,
        #                                                 nodetemplate, useCubicHermiteThroughWall, useCrossDerivatives,
        #                                                 vertebralArchLengthFactor, elongated=True,
        #                                                 elongationNodes=[rightPedicleNodes[1], rightPedicleNodes[3]],
        #                                                 elongationFactor=1.4)
        #
        # nodeIdentifier = nID

        fm.endChange()

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountUp,
                                                       refineElementsCountThroughWall)


def _createPedicleElements(nodeList, cache, nodes, coordinates, nodeIdentifier, nodetemplate,
                           useCubicHermiteThroughWall, useCrossDerivatives, archFactor, lengthFactor, elementsCountUp,
                           elementtemplateRegular, eft, mesh, elementIdentifier, dx3Values):
    newNodeList = list()
    for node in nodeList:
        cache.setNode(nodes.findNodeByIdentifier(node))
        result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, dx1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        result, dx2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        result, dx3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)

        dx4 = [(dx3[x] + (archFactor * dx1[x])) / 2 for x in range(len(dx1))]
        x4 = [x1[x] + dx4[x] for x in range(len(dx4))]
        x4[0] = x4[0] * lengthFactor

        dx3[1] = -.5*dx3[1]
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x4)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx2)
        if useCubicHermiteThroughWall:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx3)
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, [0.0, 0.0, 0.0])
            if useCubicHermiteThroughWall:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, [0.0, 0.0, 0.0])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, [0.0, 0.0, 0.0])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, [0.0, 0.0, 0.0])

        newNodeList.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1

    nodeIdStarter = 0

    for e2 in range(int(math.ceil(elementsCountUp / 2))):
        elementtemplateRegular.defineField(coordinates, -1, eft)
        element = mesh.createElement(elementIdentifier, elementtemplateRegular)

        nodeIdentifiers_1 = [nodeList[nodeIdStarter], nodeList[nodeIdStarter + 1], nodeList[nodeIdStarter + 2],
                             nodeList[nodeIdStarter + 3]]
        nodeIdentifiers_2 = [newNodeList[nodeIdStarter], newNodeList[nodeIdStarter + 1], newNodeList[nodeIdStarter + 2],
                             newNodeList[nodeIdStarter + 3]]
        # nodeIdentifiers_1 = [nodeList[nodeIdStarter + 1], newNodeList[nodeIdStarter + 1], nodeList[nodeIdStarter + 3],
        #                      newNodeList[nodeIdStarter + 3]]
        # nodeIdentifiers_2 = [nodeList[nodeIdStarter], newNodeList[nodeIdStarter], nodeList[nodeIdStarter + 2],
        #                      newNodeList[nodeIdStarter + 2]]

        nodeIdentifiers = nodeIdentifiers_1 + nodeIdentifiers_2
        result1 = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        nodeIdStarter += 2

    return newNodeList, nodeIdentifier, elementIdentifier


def _createVertebralArchNodes(nodeList, cache, nodes, coordinates, nodeIdentifier, nodetemplate, useCubicHermiteThroughWall,
                              useCrossDerivatives, lengthFactor, elongated=False, elongationNodes=None,
                              elongationFactor=0.0):
    newNodeList = list()
    for node in nodeList:
        cache.setNode(nodes.findNodeByIdentifier(node))
        result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, dx1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        result, dx2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        result, dx3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)

        if elongated:
            if node in elongationNodes:
                x1[0] = lengthFactor * x1[0] * elongationFactor
                # dx1[1] = dx1[1]*-0.1
                # dx1[0] = dx1[0]*-0.1
                # dx1 = [-1. * x for x in dx1]
            else:
                x1[0] = lengthFactor * x1[0]
                # dx1[1] = dx1[1]*-0.1
                # dx1 = [-1. * x for x in dx1]
        else:
            x1[0] = lengthFactor * x1[0]

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx1)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx2)
        if useCubicHermiteThroughWall:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx3)
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, [0.0, 0.0, 0.0])
            if useCubicHermiteThroughWall:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, [0.0, 0.0, 0.0])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, [0.0, 0.0, 0.0])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, [0.0, 0.0, 0.0])

        newNodeList.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1

    return newNodeList, nodeIdentifier


def _createVertebralArchElements(nodeList1, nodeList2, coordinates, elementsCountUp, elementtemplateRegular, eft, mesh,
                                 elementIdentifier):
    nodeIdStarter = 0
    for e2 in range(int(math.ceil(elementsCountUp / 2))):
        elementtemplateRegular.defineField(coordinates, -1, eft)
        element = mesh.createElement(elementIdentifier, elementtemplateRegular)

        nodeIdentifiers_1 = [nodeList1[nodeIdStarter], nodeList1[nodeIdStarter + 1], nodeList1[nodeIdStarter + 2],
                             nodeList1[nodeIdStarter + 3]]
        nodeIdentifiers_2 = [nodeList2[nodeIdStarter], nodeList2[nodeIdStarter + 1], nodeList2[nodeIdStarter + 2],
                             nodeList2[nodeIdStarter + 3]]

        nodeIdentifiers = nodeIdentifiers_1 + nodeIdentifiers_2
        result1 = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        nodeIdStarter += 2

    return elementIdentifier


def _createVertebralArchCentreLineNodes(nodeList, cache, nodes, coordinates, thickness, elementsCountAround,
                                        elementsCountUp, nodeIdentifier, nodetemplate, elementtemplateRegular, eft, mesh,
                                        elementIdentifier, nodesForLeftPedicleElement, nodesForRightPedicleElement):
    X = list()
    DX1 = list()
    DX2 = list()
    DX3 = list()

    finalBodyNodeIdentifier = nodeIdentifier - 1

    for node in nodeList:
        cache.setNode(nodes.findNodeByIdentifier(node))
        result, x1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
        result, dx1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
        result, dx2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
        result, dx3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
        X.append(x1), DX1.append(dx1), DX2.append(dx2), DX3.append(dx3)

    cx = [X[0], X[-1]]
    cd1 = [DX2[0], DX2[-1]]
    cd2 = [DX1[0], DX1[-1]]
    cd3 = [DX3[0], DX3[-1]]

    t2 = [1.0 - thickness] * 2
    t3 = [1.0 - thickness] * 2

    if elementsCountUp > 2:
        sx, sd1, se, sxi, _ = sampleCubicHermiteCurves(cx, cd1, int(math.ceil(elementsCountUp / 2)))
        sd2 = interpolateSampleLinear(cd2, se, sxi)
        sd3 = interpolateSampleLinear(cd3, se, sxi)
        st2 = interpolateSampleLinear(t2, se, sxi)
        st3 = interpolateSampleLinear(t3, se, sxi)
    else:
        sx, sd1, se = cx, cd1, int(math.ceil(elementsCountUp / 2))
        sd2, sd3 = cd2, cd3
        st2, st3 = t2, t3

    elementsCountUpNew = int(math.ceil(elementsCountUp / 2))

    # Find unit normals and binormals at each sample points
    sNormal = list()
    sBinormal = list()

    # Normal and binormal
    prevUnitTangent = normalise(sd1[0])
    if magnitude(crossproduct3(prevUnitTangent, [0.0, 0.0, 1.0])) > 0.0:
        prevBinormal = crossproduct3(prevUnitTangent, [0.0, 0.0, 1.0])
    else:
        prevBinormal = crossproduct3(prevUnitTangent, [0.0, -1.0, 0.0])
    prevUnitBinormal = normalise(prevBinormal)
    prevUnitNormal = crossproduct3(prevUnitBinormal, prevUnitTangent)
    sNormal.append(prevUnitNormal)
    sBinormal.append(prevUnitBinormal)
    for n in range(1, elementsCountUpNew + 1):
        unitTangent = normalise(sd1[n])
        cp = crossproduct3(prevUnitTangent, unitTangent)
        if magnitude(cp) > 0.0:
            axisRot = normalise(cp)
            thetaRot = math.acos(dotproduct(prevUnitTangent, unitTangent))
            rotFrame = rotationMatrix(thetaRot, axisRot)
            rotNormal = [rotFrame[j][0] * prevUnitNormal[0] + rotFrame[j][1] * prevUnitNormal[1] + rotFrame[j][2] *
                         prevUnitNormal[2] for j in range(3)]
            unitNormal = normalise(rotNormal)
            unitBinormal = crossproduct3(unitTangent, unitNormal)
            prevUnitTangent = unitTangent
            prevUnitNormal = unitNormal
        else:
            unitBinormal = prevUnitBinormal
            unitNormal = prevUnitNormal
        sNormal.append(unitNormal)
        sBinormal.append(unitBinormal)

    elementsCountThroughWall = 2
    useCubicHermiteThroughWall = True
    useCrossDerivatives = False

    if elementsCountAround > 7:
        pedicleCoefficient = 2.5 * (elementsCountAround / (7. * 3.))
    elif elementsCountAround < 7:
        pedicleCoefficient = 0.3 * (elementsCountAround * (7. / 3.))
    else:
        pedicleCoefficient = elementsCountAround / (7. / 3.)
    print(elementsCountAround)
    print(pedicleCoefficient)

    """ Nodes """
    for n3 in range(elementsCountThroughWall):
        for n2 in range(elementsCountUpNew + 1):
            aThroughWallElement = sd2[n2][1] + st2[n2] * (n3 / elementsCountThroughWall)
            bThroughWallElement = sd3[n2][0] + st3[n2] * (n3 / elementsCountThroughWall)
            perimeterAroundWallElement = getApproximateEllipsePerimeter(aThroughWallElement, bThroughWallElement)
            arcLengthPerElementAround = perimeterAroundWallElement / elementsCountAround
            prevRadiansAround = updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement, 0.0,
                                                              arcLengthPerElementAround)
            st2PerWallElement = st2[n2] / elementsCountThroughWall
            st3PerWallElement = st3[n2] / elementsCountThroughWall

            if n2 < elementsCountUpNew:
                aThroughWallElementNext = sd2[n2 + 1][1] + st2[n2 + 1] * (n3 / elementsCountThroughWall)
                bThroughWallElementNext = sd3[n2 + 1][0] + st3[n2 + 1] * (n3 / elementsCountThroughWall)
                perimeterAroundWallElementNext = getApproximateEllipsePerimeter(aThroughWallElementNext,
                                                                                bThroughWallElementNext)
                arcLengthPerElementAroundNext = perimeterAroundWallElementNext / elementsCountAround

            for n1 in range(elementsCountAround):
                if n1 == 0 or n1 == 1 or n1 == elementsCountAround - 1:
                    arcLengthAround = n1 * arcLengthPerElementAround
                    radiansAround = -1 * updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement, 0.0,
                                                                       arcLengthAround)
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x = [sx[n2][j] + aThroughWallElement * cosRadiansAround * sBinormal[n2][
                        j] + bThroughWallElement * sinRadiansAround * sNormal[n2][j] for j in range(3)]
                    dx_ds1 = [(radiansAround - prevRadiansAround) * (
                            aThroughWallElement * -sinRadiansAround * sBinormal[n2][
                        j] + bThroughWallElement * cosRadiansAround * sNormal[n2][j]) for j in range(3)]

                    # Calculate curvature to find d1 for node
                    unitNormal = normalise([aThroughWallElement * cosRadiansAround * sBinormal[n2][
                        j] + bThroughWallElement * sinRadiansAround * sNormal[n2][j] for j in range(3)])

                    if n2 == 0:
                        curvature = getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2 + 1], sd1[n2 + 1], unitNormal, 0.0)
                    elif n2 == elementsCountUpNew:
                        curvature = getCubicHermiteCurvature(sx[n2 - 1], sd1[n2 - 1], sx[n2], sd1[n2], unitNormal, 1.0)
                    else:
                        curvature = 0.5 * (getCubicHermiteCurvature(sx[n2 - 1], sd1[n2 - 1], sx[n2], sd1[n2], unitNormal,
                                                                    1.0) + getCubicHermiteCurvature(sx[n2], sd1[n2],
                                                                                                    sx[n2 + 1], sd1[n2 + 1],
                                                                                                    unitNormal, 0.0))
                    wallDistance = magnitude([aThroughWallElement * cosRadiansAround * sBinormal[n2][
                        j] + bThroughWallElement * sinRadiansAround * sNormal[n2][j] for j in range(3)])
                    factor = 1.0 - curvature * wallDistance
                    d1Wall = [factor * c for c in sd1[n2]]

                    # Calculate curvature to find d1 for downstream node
                    if n2 < elementsCountUpNew:
                        arcLengthAroundNext = n1 * arcLengthPerElementAroundNext
                        radiansAroundNext = -1 * updateEllipseAngleByArcLength(aThroughWallElementNext,
                                                                               bThroughWallElementNext, 0.0,
                                                                               arcLengthAroundNext)
                        cosRadiansAroundNext = math.cos(radiansAroundNext)
                        sinRadiansAroundNext = math.sin(radiansAroundNext)
                        xNext = [sx[n2 + 1][j] + aThroughWallElementNext * cosRadiansAroundNext * sBinormal[n2 + 1][
                            j] + bThroughWallElementNext * sinRadiansAroundNext * sNormal[n2 + 1][j] for j in range(3)]
                        unitNormalNext = normalise([aThroughWallElementNext * cosRadiansAroundNext * sBinormal[n2 + 1][
                            j] + bThroughWallElementNext * sinRadiansAroundNext * sNormal[n2 + 1][j] for j in range(3)])
                        if n2 + 1 == elementsCountUpNew:
                            curvatureNext = getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2 + 1], sd1[n2 + 1],
                                                                     unitNormalNext, 1.0)
                        else:
                            curvatureNext = 0.5 * (
                                    getCubicHermiteCurvature(sx[n2], sd1[n2], sx[n2 + 1], sd1[n2 + 1], unitNormalNext,
                                                             1.0) + getCubicHermiteCurvature(sx[n2 + 1], sd1[n2 + 1],
                                                                                             sx[n2 + 2], sd1[n2 + 2],
                                                                                             unitNormalNext, 0.0))
                        wallDistanceNext = magnitude([aThroughWallElementNext * cosRadiansAroundNext * sBinormal[n2 + 1][
                            j] + bThroughWallElementNext * sinRadiansAroundNext * sNormal[n2 + 1][j] for j in range(3)])
                        factorNext = 1.0 - curvatureNext * wallDistanceNext
                        d1WallNext = [factorNext * c for c in sd1[n2 + 1]]
                        arcLength = computeCubicHermiteArcLength(x, d1Wall, xNext, d1WallNext, True)
                        dx_ds2 = [arcLength * c for c in normalise(d1Wall)]
                        if n2 == elementsCountUpNew - 1:
                            secondLastX = x
                            secondLastd1Wall = d1Wall
                            lastX = xNext
                            lastd1Wall = d1WallNext
                    else:
                        arcLength = computeCubicHermiteArcLength(secondLastX, secondLastd1Wall, lastX, lastd1Wall, True)
                        dx_ds2 = [arcLength * c for c in normalise(lastd1Wall)]

                    dx_ds3 = [
                        st2PerWallElement * cosRadiansAround * sBinormal[n2][j] + st3PerWallElement * sinRadiansAround *
                        sNormal[n2][j] for j in range(3)]  # Modify later to calculate with interpolation

                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    if n1 == 1:
                        x[1] = x[1]*pedicleCoefficient
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [x*0.5 for x in dx_ds1])
                    elif n1 == elementsCountAround - 1:
                        x[1] = x[1]*pedicleCoefficient
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [x*0.1 for x in dx_ds1])
                    else:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    if useCubicHermiteThroughWall:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, [0., 0., 0.])
                        if useCubicHermiteThroughWall:
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, [0., 0., 0.])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, [0., 0., 0.])
                            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, [0., 0., 0.])
                    nodeIdentifier = nodeIdentifier + 1
                    prevRadiansAround = radiansAround
                else:
                    pass

    """ Left pedicle elements """
    nodeIdStarter = 0
    for e2 in range(elementsCountUpNew):
        elementtemplateRegular.defineField(coordinates, -1, eft)
        element = mesh.createElement(elementIdentifier, elementtemplateRegular)

        if elementsCountUp == 3:
            bni1 = 3*e2+(finalBodyNodeIdentifier + 3) + (2 * elementsCountUp + (elementsCountUp%3) + 3)
        else:
            bni1 = 3*e2+(finalBodyNodeIdentifier + 3) + (2 * elementsCountUp + (elementsCountUp%3))
        bni2 = 3*e2+(finalBodyNodeIdentifier + 3)
        if elementsCountUp == 3:
            bni3 = (bni2 + 3) + (2 * elementsCountUp + (elementsCountUp%3) + 3)
        else:
            bni3 = (bni2 + 3) + (2 * elementsCountUp + (elementsCountUp%3))
        bni4 = bni2 + 3

        nodeIdentifiers1 = [nodesForLeftPedicleElement[nodeIdStarter], nodesForLeftPedicleElement[nodeIdStarter + 1],
                            nodesForLeftPedicleElement[nodeIdStarter + 2], nodesForLeftPedicleElement[nodeIdStarter + 3]]
        nodeIdentifiers2 = [bni1, bni2, bni3, bni4]
        nodeIdentifiers = nodeIdentifiers1 + nodeIdentifiers2
        result1 = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1
        nodeIdStarter += 2

    """ Right pedicle elements """
    nodeIdStarter = 0
    for e2 in range(elementsCountUpNew):
        elementtemplateRegular.defineField(coordinates, -1, eft)
        element = mesh.createElement(elementIdentifier, elementtemplateRegular)

        bni1 = 3*e2+(finalBodyNodeIdentifier + 2)
        if elementsCountUp == 3:
            bni2 = 3*e2+(finalBodyNodeIdentifier + 2) + (2 * elementsCountUp + (elementsCountUp%3) + 3)
        else:
            bni2 = 3*e2+(finalBodyNodeIdentifier + 2) + (2 * elementsCountUp + (elementsCountUp%3))
        if elementsCountUp == 3:
            bni3 = bni2 - (2 * elementsCountUp + (elementsCountUp%3))
        else:
            bni3 = bni2 - (2 * elementsCountUp + (elementsCountUp%3) - 3)
        if elementsCountUp == 3:
            bni4 = bni3 + (2 * elementsCountUp + (elementsCountUp%3) + 3)
        else:
            bni4 = bni3 + (2 * elementsCountUp + (elementsCountUp%3))

        nodeIdentifiers1 = [nodesForRightPedicleElement[nodeIdStarter], nodesForRightPedicleElement[nodeIdStarter + 1],
                            nodesForRightPedicleElement[nodeIdStarter + 2], nodesForRightPedicleElement[nodeIdStarter + 3]]
        nodeIdentifiers2 = [bni1, bni2, bni3, bni4]
        nodeIdentifiers = nodeIdentifiers1 + nodeIdentifiers2
        result1 = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1
        nodeIdStarter += 2

    return None