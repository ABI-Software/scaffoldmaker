"""
Generates a 3-D elliptic cylinder model.
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
        return {'Number of elements around': 6, 'Number of elements up': 3, 'Number of elements through wall': 2,
                'Major diameter': 1.1, 'Minor diameter': 0.9, 'Height': 0.3, 'Body thickness ratio': 0.7,
                'Use cross derivatives': False, 'Refine': False, 'Refine number of elements around': 1,
                'Refine number of elements through wall': 1, 'Refine number of elements up': 1}

    @staticmethod
    def getOrderedOptionNames():
        return ['Number of elements around', 'Number of elements up',
                'Number of elements through wall', 'Major diameter', 'Minor diameter', 'Height', 'Body thickness ratio',
                'Use cross derivatives', 'Refine', 'Refine number of elements around', 'Refine number of elements up',
                'Refine number of elements through wall']

    @staticmethod
    def checkOptions(options):
        for key in ['Refine number of elements around', 'Refine number of elements up',
                    'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements up'] < 2:
            options['Number of elements up'] = 2
        if options['Number of elements around'] < 4:
            options['Number of elements around'] = 4
        if options['Number of elements through wall'] < 2:
            options['Number of elements through wall'] = 2
        elif options['Number of elements through wall'] > 4:
            options['Number of elements through wall'] = 4
        if options['Major diameter'] < 0.0:
            options['Diameter'] = 0.0
        if options['Minor diameter'] < 0.0:
            options['Diameter'] = 0.0
        if options['Height'] < 0.0:
            options['Height'] = 0.0
        if options['Body thickness ratio'] <= 0.0:
            options['Body thickness ratio'] = 0.55
        elif options['Body thickness ratio'] > 1.0:
            options['Body thickness ratio'] = 0.9

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
        elementsCountThroughWall = options['Number of elements through wall']
        useCrossDerivatives = options['Use cross derivatives']
        axisA = options['Major diameter']
        axisB = options['Minor diameter']
        height = options['Height']
        useCubicHermiteThroughWall = True

        zero = [0.0, 0.0, 0.0]

        """ Get nodes along the centre line """
        cx = [zero, [0.0, 0.0, height]]
        cd1 = [[0.0, 0.0, height / elementsCountUp], [0.0, 0.0, height / elementsCountUp]]
        cd2 = [[0.0, axisB / elementsCountThroughWall, 0.0], [0.0, axisB / elementsCountThroughWall, 0.0]]
        cd3 = [[axisA / elementsCountThroughWall, 0.0, 0.0], [axisA / elementsCountThroughWall, 0.0, 0.0]]

        t = options['Body thickness ratio']
        t2 = [1.0 - options['Body thickness ratio']] * 2
        t2d = [0.0, 0.0]
        t3 = [1.0 - options['Body thickness ratio']] * 2
        t3d = [0.0, 0.0]

        sx, sd1, se, sxi, _ = sampleCubicHermiteCurves(cx, cd1, elementsCountUp)
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
        if magnitude(crossproduct3(prevUnitTangent, [0.0, 0.0, 1.0])) > 0.0:
            prevBinormal = crossproduct3(prevUnitTangent, [0.0, 0.0, 1.0])
        else:
            prevBinormal = crossproduct3(prevUnitTangent, [0.0, -1.0, 0.0])
        prevUnitBinormal = normalise(prevBinormal)
        prevUnitNormal = crossproduct3(prevUnitBinormal, prevUnitTangent)
        sNormal.append(prevUnitNormal)
        sBinormal.append(prevUnitBinormal)

        # Step through central line and rotate central line axes to align tangent
        # to tangent from previous level
        for n in range(1, elementsCountUp + 1):
            unitTangent = normalise(sd1[n])
            cp = crossproduct3(prevUnitTangent, unitTangent)
            if magnitude(cp) > 0.0:
                axisRot = normalise(cp)
                thetaRot = math.acos(dotproduct(prevUnitTangent, unitTangent))
                rotFrame = rotationMatrixAboutAxis(axisRot, thetaRot)
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
            # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1,
            #                               [0.0, 0.0, 0.0])
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

        # Create nodes on the surface
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
                bni1 = elementsCountUp + 2
                for e2 in range(elementsCountUp):
                    aThroughWallElement = sd2[e2][1] + st2[e2] * (e3 / elementsCountThroughWall)
                    bThroughWallElement = sd3[e2][0] + st3[e2] * (e3 / elementsCountThroughWall)
                    perimeterAroundWallElement = getApproximateEllipsePerimeter(aThroughWallElement,
                                                                                bThroughWallElement)
                    arcLengthPerElementAround = perimeterAroundWallElement / elementsCountAround
                    for e1 in range(elementsCountAround):
                        # arcLengthAround = e1 * arcLengthPerElementAround
                        # radiansAround = -1 * updateEllipseAngleByArcLength(aThroughWallElement, bThroughWallElement,
                        #                                                    0.0, arcLengthAround)
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
                            if e1 <= elementsCountAround - 2:
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1 + e1, bni2 + e1, bni3 + e1, bni4 + e1]
                            else:
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1 + e1, bni2 + e1 - elementsCountAround,
                                                   bni3 + e1, bni4 + e1 - elementsCountAround]
                        else:
                            if e1 <= elementsCountAround - 2:
                                bni1 = (e2 * elementsCountAround) + 2 + elementsCountUp
                                bni2 = bni1 + 1
                                bni3 = bni2 + (elementsCountAround - 1)
                                bni4 = bni3 + 1
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1 + e1, bni2 + e1, bni3 + e1, bni4 + e1]
                            else:
                                bni1 = (e2 * elementsCountAround) + 2 + elementsCountUp
                                bni2 = bni1 - (elementsCountAround - 1)
                                bni3 = bni1 + 1 + (elementsCountAround - 1)
                                bni4 = bni3 - (elementsCountAround - 1)
                                nodeIdentifiers = [e2 + 1, e2 + 2, bni1 + e1, bni2 + e1, bni3 + e1, bni4 + e1]
                        result1 = element.setNodesByIdentifier(eft2, nodeIdentifiers)
                        scalefactors = [-1.0, 1.0, 1.0, 1.0, 1.0]
                        result2 = element.setScaleFactors(eft2, scalefactors)
                        elementIdentifier = elementIdentifier + 1

            # Create regular elements
            elif e3 == 1:
                for e2 in range(elementsCountUp):
                    for e1 in range(elementsCountAround):

                        elementtemplateRegular.defineField(coordinates, -1, eft)
                        element = mesh.createElement(elementIdentifier, elementtemplateRegular)

                        if e1 <= elementsCountAround - 2:
                            bni1 = e1 + (e2 * elementsCountAround) + 2 + elementsCountUp
                            bni2 = bni1 + 1
                            bni3 = bni2 + (elementsCountAround - 1)
                            bni4 = bni3 + 1
                            bni5 = bni1 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni6 = bni5 + 1
                            bni7 = bni3 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni8 = bni7 + 1
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                        else:
                            bni1 = e1 + (e2 * elementsCountAround) + 2 + elementsCountUp
                            bni2 = bni1 - (elementsCountAround - 1)
                            bni3 = bni1 + elementsCountAround
                            bni4 = bni3 - (elementsCountAround - 1)
                            bni5 = bni1 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni6 = bni5 - (elementsCountAround - 1)
                            bni7 = bni3 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni8 = bni7 - (elementsCountAround - 1)
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                        result1 = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1
            else:
                for e2 in range(elementsCountUp):
                    for e1 in range(elementsCountAround):

                        elementtemplateRegular.defineField(coordinates, -1, eft)
                        element = mesh.createElement(elementIdentifier, elementtemplateRegular)

                        nx = (e3-1)*((elementsCountUp * elementsCountAround) + elementsCountAround)
                        if e1 <= elementsCountAround - 2:
                            bni1 = (e1 + (e2 * elementsCountAround) + 2 + elementsCountUp) + nx
                            bni2 = bni1 + 1
                            bni3 = bni2 + (elementsCountAround - 1)
                            bni4 = bni3 + 1
                            bni5 = bni1 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni6 = bni5 + 1
                            bni7 = bni3 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni8 = bni7 + 1
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                        else:
                            bni1 = (e1 + (e2 * elementsCountAround) + 2 + elementsCountUp) + nx
                            bni2 = bni1 - (elementsCountAround - 1)
                            bni3 = bni1 + elementsCountAround
                            bni4 = bni3 - (elementsCountAround - 1)
                            bni5 = bni1 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni6 = bni5 - (elementsCountAround - 1)
                            bni7 = bni3 + ((elementsCountUp * elementsCountAround) + elementsCountAround)
                            bni8 = bni7 - (elementsCountAround - 1)
                            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                        result1 = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1

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
    rotMatrix = ([[rotAxis[0] * rotAxis[0] * C + cosTheta, rotAxis[0] * rotAxis[1] * C - rotAxis[2] * sinTheta,
                   rotAxis[0] * rotAxis[2] * C + rotAxis[1] * sinTheta],
                  [rotAxis[1] * rotAxis[0] * C + rotAxis[2] * sinTheta, rotAxis[1] * rotAxis[1] * C + cosTheta,
                   rotAxis[1] * rotAxis[2] * C - rotAxis[0] * sinTheta],
                  [rotAxis[2] * rotAxis[0] * C - rotAxis[1] * sinTheta,
                   rotAxis[2] * rotAxis[1] * C + rotAxis[0] * sinTheta, rotAxis[2] * rotAxis[2] * C + cosTheta]])
    return rotMatrix
