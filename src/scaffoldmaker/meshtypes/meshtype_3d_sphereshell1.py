"""
Generates a 3-D unit sphere shell mesh with variable numbers of elements
around, up and through the thickness.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_sphereshell1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Sphere Shell 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around' : 4,
            'Number of elements up' : 4,
            'Number of elements through wall' : 1,
            'Wall thickness' : 0.25,
            'Wall thickness ratio apex' : 1.0,
            'Exclude bottom rows' : 0,
            'Exclude top rows' : 0,
            'Length ratio' : 1.0,
            'Element length ratio equator/apex' : 1.0,
            'Use cross derivatives' : False,
            'Use linear through wall' : False,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements up' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around',
            'Number of elements up',
            'Number of elements through wall',
            'Exclude bottom rows',
            'Exclude top rows',
            'Wall thickness',
            'Wall thickness ratio apex',
            'Length ratio',
            'Element length ratio equator/apex',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements up',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements up',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements up'] < 2:
            options['Number of elements up'] = 2
        if options['Number of elements around'] < 2:
            options['Number of elements around'] = 2
        if options['Exclude bottom rows'] < 0:
            options['Exclude bottom rows'] = 0
        if options['Exclude top rows'] < 0:
            options['Exclude top rows'] = 0
        if options['Wall thickness'] < 0.0:
            options['Wall thickness'] = 0.0
        elif options['Wall thickness'] > 0.5:
            options['Wall thickness'] = 0.5
        if options['Wall thickness ratio apex'] < 0.0:
            options['Wall thickness ratio apex'] = 0.0
        if options['Length ratio'] < 1.0E-6:
            options['Length ratio'] = 1.0E-6
        if options['Element length ratio equator/apex'] < 1.0E-6:
            options['Element length ratio equator/apex'] = 1.0E-6

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite or bicubic Hermite linear mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        elementsCountAround = options['Number of elements around']
        elementsCountUp = options['Number of elements up']
        elementsCountThroughWall = options['Number of elements through wall']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        excludeBottomRows = options['Exclude bottom rows']
        excludeTopRows = options['Exclude top rows']
        wallThickness = options['Wall thickness']
        wallThicknessRatioApex = options['Wall thickness ratio apex']
        lengthRatio = options['Length ratio']
        elementLengthRatioEquatorApex = options['Element length ratio equator/apex']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateApex = nodes.createNodetemplate()
        nodetemplateApex.defineField(coordinates)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCubicHermiteThroughWall:
            nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
            if useCubicHermiteThroughWall:
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
        else:
            nodetemplate = nodetemplateApex

        mesh = fm.findMeshByDimension(3)

        if useCubicHermiteThroughWall:
            eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        else:
            eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 0.0 ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]

        # pre-calculate positions and tangent/normal vectors up (elementsCountUp + 1) node layers
        outerWidth = 0.5
        outerLength = outerWidth*lengthRatio
        bOuter = 2.0 / (1.0 + elementLengthRatioEquatorApex / lengthRatio)
        aOuter = 1.0 - bOuter

        innerWidth = outerWidth - wallThickness
        innerLength = outerLength - wallThickness*wallThicknessRatioApex
        lengthRatioInner = innerLength/innerWidth if (innerWidth > 0.0) else lengthRatio
        bInner = 2.0 / (1.0 + elementLengthRatioEquatorApex / lengthRatio)
        aInner = 1.0 - bInner

        positionOuterArray = [(0,0)]*(elementsCountUp + 1)
        positionInnerArray = [(0,0)]*(elementsCountUp + 1)
        radiansUpOuterArray = [0]*(elementsCountUp + 1)
        radiansUpInnerArray = [0]*(elementsCountUp + 1)
        vector2OuterArray = [(0,0)]*(elementsCountUp + 1)
        vector2InnerArray = [(0,0)]*(elementsCountUp + 1)
        for n2 in range(elementsCountUp + 1):
            if n2*2 <= elementsCountUp:
                xi = n2*2 / elementsCountUp
            else:
                xi = 2.0 - (n2*2 / elementsCountUp)

            nxiOuter = aOuter*xi*xi + bOuter*xi
            dnxiOuter = 2.0*aOuter*xi + bOuter
            radiansUpOuterArray[n2] = radiansUpOuter = nxiOuter*math.pi*0.5 if (n2*2 <= elementsCountUp) else (math.pi - nxiOuter*math.pi*0.5)
            dRadiansUpOuter = dnxiOuter*math.pi/elementsCountUp
            cosRadiansUpOuter = math.cos(radiansUpOuter);
            sinRadiansUpOuter = math.sin(radiansUpOuter);
            positionOuterArray[n2] = positionOuter = ( outerWidth*sinRadiansUpOuter, -outerLength*cosRadiansUpOuter )
            vector2OuterArray[n2] = ( outerWidth*cosRadiansUpOuter*dRadiansUpOuter, outerLength*sinRadiansUpOuter*dRadiansUpOuter )

            nxiInner = aInner*xi*xi + bInner*xi
            dnxiInner = 2.0*aInner*xi + bInner
            radiansUpInnerArray[n2] = radiansUpInner = nxiInner*math.pi*0.5 if (n2*2 <= elementsCountUp) else (math.pi - nxiInner*math.pi*0.5)
            dRadiansUpInner = dnxiInner*math.pi/elementsCountUp
            cosRadiansUpInner = math.cos(radiansUpInner);
            sinRadiansUpInner = math.sin(radiansUpInner);
            positionInnerArray[n2] = positionInner = ( innerWidth*sinRadiansUpInner, -innerLength*cosRadiansUpInner )
            vector2InnerArray[n2] = ( innerWidth*cosRadiansUpInner*dRadiansUpInner, innerLength*sinRadiansUpInner*dRadiansUpInner )

        # now create the nodes
        for n3 in range(elementsCountThroughWall + 1):

            n3_fraction = n3 / elementsCountThroughWall

            for n2 in range(excludeBottomRows, elementsCountUp + 1 - excludeTopRows):

                positionOuter = positionOuterArray[n2]
                positionInner = positionInnerArray[n2]
                position = (
                    positionOuter[0]*n3_fraction + positionInner[0]*(1.0 - n3_fraction),
                    positionOuter[1]*n3_fraction + positionInner[1]*(1.0 - n3_fraction)
                )

                radiansUpOuter = radiansUpOuterArray[n2]
                sinRadiansUpOuter = math.sin(radiansUpOuter)
                radiansUpInner = radiansUpInnerArray[n2]
                sinRadiansUpInner = math.sin(radiansUpInner)

                vector2Outer = vector2OuterArray[n2]
                vector2Inner = vector2InnerArray[n2]
                vector2 = (
                    vector2Outer[0]*n3_fraction + vector2Inner[0]*(1.0 - n3_fraction),
                    vector2Outer[1]*n3_fraction + vector2Inner[1]*(1.0 - n3_fraction)
                )

                vector3 = (
                    (positionOuter[0] - positionInner[0])/elementsCountThroughWall,
                    (positionOuter[1] - positionInner[1])/elementsCountThroughWall
                )

                if n2 == 0:
                    # create bottom apex node
                    node = nodes.createNode(nodeIdentifier, nodetemplateApex)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, position[1] ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, vector2[0], 0.0 ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ vector2[0], 0.0, 0.0 ])
                    if useCubicHermiteThroughWall:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, vector3[1] ])
                    nodeIdentifier = nodeIdentifier + 1

                elif n2 < elementsCountUp:
                    # create regular rows between apexes
                    for n1 in range(elementsCountAround):
                        radiansAround = n1*radiansPerElementAround
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)
                        x = [
                            position[0]*cosRadiansAround,
                            position[0]*sinRadiansAround,
                            position[1]
                        ]
                        dx_ds1 = [
                            position[0]*-sinRadiansAround*radiansPerElementAround,
                            position[0]*cosRadiansAround*radiansPerElementAround,
                            0.0
                        ]
                        dx_ds2 = [
                            vector2[0]*cosRadiansAround,
                            vector2[0]*sinRadiansAround,
                            vector2[1]
                        ]
                        dx_ds3 = [
                            vector3[0]*cosRadiansAround,
                            vector3[0]*sinRadiansAround,
                            vector3[1]
                        ]
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

                else:
                    # create top apex node
                    node = nodes.createNode(nodeIdentifier, nodetemplateApex)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, position[1] ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, vector2[0], 0.0 ])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ -vector2[0], 0.0, 0.0 ])
                    if useCubicHermiteThroughWall:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [ 0.0, 0.0, vector3[1] ])
                    nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        # now (node offset through wall) varies with number of excluded rows
        now = 0;
        if excludeBottomRows == 0:
            now += 1  # bottom apex node
        if excludeTopRows == 0:
            now += 1  # top apex node
        fullNodeRows = elementsCountUp - 1
        if excludeBottomRows > 1:
            fullNodeRows -= (excludeBottomRows - 1)
        if excludeTopRows > 1:
            fullNodeRows -= (excludeTopRows - 1)
        if fullNodeRows > 0:
            now += fullNodeRows*elementsCountAround
        row2NodeOffset = 2 if (excludeBottomRows == 0) else 1
        for e3 in range(elementsCountThroughWall):
            no = e3*now

            if excludeBottomRows == 0:
                # create bottom apex elements, editing eft scale factor identifiers around apex
                # scale factor identifiers follow convention of offsetting by 100 for each 'version'
                elementtemplate1 = mesh.createElementtemplate()
                elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                for e1 in range(elementsCountAround):
                    va = e1
                    vb = (e1 + 1)%elementsCountAround
                    eft1 = eftfactory.createEftShellPoleBottom(va*100, vb*100)
                    elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    bni1 = no + 1
                    bni2 = no + e1 + 2
                    bni3 = no + (e1 + 1)%elementsCountAround + 2
                    nodeIdentifiers = [ bni1, bni2, bni3, bni1 + now, bni2 + now, bni3 + now ]
                    element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    # set general linear map coefficients
                    radiansAround = e1*radiansPerElementAround
                    radiansAroundNext = ((e1 + 1)%elementsCountAround)*radiansPerElementAround
                    scalefactors = [
                        -1.0,
                        math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                        math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                        math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                        math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
                    ]
                    result = element.setScaleFactors(eft1, scalefactors)
                    elementIdentifier = elementIdentifier + 1

            # create regular rows between apexes
            rowLimit = (elementsCountUp - 2)
            if excludeBottomRows > 1:
                rowLimit -= (excludeBottomRows - 1)
            if excludeTopRows > 1:
                rowLimit -= (excludeTopRows - 1)
            for e2 in range(0, rowLimit):
                for e1 in range(elementsCountAround):
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    bni11 = no + e2*elementsCountAround + e1 + row2NodeOffset
                    bni12 = no + e2*elementsCountAround + (e1 + 1)%elementsCountAround + row2NodeOffset
                    bni21 = no + (e2 + 1)*elementsCountAround + e1 + row2NodeOffset
                    bni22 = no + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + row2NodeOffset
                    nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                    elementIdentifier = elementIdentifier + 1

            if excludeTopRows == 0:
                # create top apex elements, editing eft scale factor identifiers around apex
                # scale factor identifiers follow convention of offsetting by 100 for each 'version'
                elementtemplate1 = mesh.createElementtemplate()
                elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                for e1 in range(elementsCountAround):
                    va = e1
                    vb = (e1 + 1)%elementsCountAround
                    eft1 = eftfactory.createEftShellPoleTop(va*100, vb*100)
                    elementtemplate1.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    bni3 = no + now
                    bni1 = bni3 - elementsCountAround + e1
                    bni2 = bni3 - elementsCountAround + (e1 + 1)%elementsCountAround
                    nodeIdentifiers = [ bni1, bni2, bni3, bni1 + now, bni2 + now, bni3 + now ]
                    element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    # set general linear map coefficients
                    radiansAround = math.pi + e1*radiansPerElementAround
                    radiansAroundNext = math.pi + ((e1 + 1)%elementsCountAround)*radiansPerElementAround
                    scalefactors = [
                        -1.0,
                        -math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                        -math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                        -math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                        -math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
                    ]
                    result = element.setScaleFactors(eft1, scalefactors)
                    elementIdentifier = elementIdentifier + 1

        if False:  # (wallThickness < 0.5):  # and (wallThicknessRatioApex != 1.0):
            r = fm.createFieldMagnitude(coordinates)
            const05 = fm.createFieldConstant([0.5])
            d = fm.createFieldSubtract(const05, r)
            d_r = fm.createFieldDivide(d, r)
            rRatio = 1.0 - wallThicknessRatioApex
            rScale = fm.createFieldConstant([rRatio])
            zScale = fm.createFieldMultiply(d_r, rScale)
            one = fm.createFieldConstant([1.0])
            one_plus_zScale = fm.createFieldAdd(one, zScale)
            scale = fm.createFieldConcatenate([one, one, one_plus_zScale])
            newCoordinates = fm.createFieldMultiply(coordinates, scale)
            fieldassignment = coordinates.createFieldassignment(newCoordinates)
            fieldassignment.assign()

        fm.endChange()
        return []

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshrefinement, MeshRefinement)
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountThroughWall = options['Refine number of elements through wall']
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountUp, refineElementsCountThroughWall)
