"""
Generates 3-D Left and Right ventricles mesh starting from modified sphere shell mesh.
"""

import math
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_3d_sphereshell1 import MeshType_3d_sphereshell1
from mapclientplugins.meshgeneratorstep.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_heartventricles1:
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Heart Ventricles 1'

    @staticmethod
    def getDefaultOptions():
        return {
            'Number of elements up' : 4,
            'Number of elements around' : 10,
            'Number of elements across septum' : 4,
            'Number of elements below septum' : 2,
            'Number of elements through LV wall' : 1,
            'LV wall thickness' : 0.15,
            'LV wall thickness ratio apex' : 0.5,
            'RV free wall thickness' : 0.05,
            'RV width' : 0.2,
            'Length ratio' : 2.0,
            'Element length ratio equator/apex' : 1.0,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up',
            'Number of elements around',
            'Number of elements across septum',
            'Number of elements below septum',
            'Number of elements through LV wall',
            'LV wall thickness',
            'LV wall thickness ratio apex',
            'RV free wall thickness',
            'RV width',
            'Length ratio',
            'Element length ratio equator/apex'
        ]

    @staticmethod
    def checkOptions(options):
        if options['Number of elements up'] < 4:
            options['Number of elements up'] = 4
        if options['Number of elements around'] < 4:
            options['Number of elements around'] = 4
        if options['Number of elements across septum'] < 2:
            options['Number of elements across septum'] = 2
        elif options['Number of elements across septum'] > (options['Number of elements around'] - 2):
            options['Number of elements across septum'] = options['Number of elements around'] - 2
        if options['Number of elements below septum'] < 2:
            options['Number of elements below septum'] = 2
        elif options['Number of elements below septum'] > (options['Number of elements up'] - 1):
            options['Number of elements below septum'] = options['Number of elements up'] - 1
        if (options['Number of elements through LV wall'] < 1) :
            options['Number of elements through LV wall'] = 1
        if options['LV wall thickness'] < 0.0:
            options['LV wall thickness'] = 0.0
        elif options['LV wall thickness'] > 0.5:
            options['LV wall thickness'] = 0.5
        if options['LV wall thickness ratio apex'] < 0.0:
            options['LV wall thickness ratio apex'] = 0.0
        if options['RV free wall thickness'] < 0.0:
            options['RV free wall thickness'] = 0.0
        if options['RV width'] < 0.0:
            options['RV width'] = 0.0
        if options['Length ratio'] < 1.0E-6:
            options['Length ratio'] = 1.0E-6
        if options['Element length ratio equator/apex'] < 1.0E-6:
            options['Element length ratio equator/apex'] = 1.0E-6

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountUp = options['Number of elements up']
        elementsCountAround = options['Number of elements around']
        elementsCountAcrossSeptum = options['Number of elements across septum']
        elementsCountBelowSeptum = options['Number of elements below septum']
        elementsCountThroughLVWall = options['Number of elements through LV wall']
        LVWallThickness = options['LV wall thickness']
        LVWallThicknessRatioApex = options['LV wall thickness ratio apex']
        RVWidth = options['RV width']
        RVFreeWallThickness = options['RV free wall thickness']
        useCrossDerivatives = options['Use cross derivatives']

        # generate a half sphere shell which will be edited
        sphereShellOptions = MeshType_3d_sphereshell1.getDefaultOptions()
        sphereShellOptions['Number of elements up'] = elementsCountUp*2
        sphereShellOptions['Number of elements around'] = elementsCountAround
        sphereShellOptions['Number of elements through wall'] = elementsCountThroughLVWall
        sphereShellOptions['Exclude top rows'] = elementsCountUp
        sphereShellOptions['Wall thickness'] = LVWallThickness
        sphereShellOptions['Wall thickness ratio apex'] = LVWallThicknessRatioApex
        sphereShellOptions['Length ratio'] = options['Length ratio']
        sphereShellOptions['Element length ratio equator/apex'] = options['Element length ratio equator/apex']
        MeshType_3d_sphereshell1.generateMesh(region, sphereShellOptions)

        fm = region.getFieldmodule()
        fm.beginChange()
        # find the coordinates field created for the sphere shell
        coordinates = fm.findFieldByName('coordinates').castFiniteElement()

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        mesh = fm.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()

        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        # RV free wall inner bottom elements
        eftRVFreeWallInnerBottom = tricubichermite.createEftBasic()
        eftRVFreeWallInnerBottom.setNumberOfLocalScaleFactors(4)
        for s in range(4):
            eftRVFreeWallInnerBottom.setScaleFactorType(s + 1, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eftRVFreeWallInnerBottom.setScaleFactorIdentifier(s + 1, (s % 2) + 3)  # allow for node sf ids 1,2 on inside elements
        # map d/dxi2 = d/ds3, d/dxi3 = a.d/ds2 + b.d/ds3 on RV apex
        for n in range(2):
            ln = n + 1
            eftRVFreeWallInnerBottom.setTermNodeParameter(n*8 + 3, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            # 2 terms for d/dxi3 via general linear map
            eftRVFreeWallInnerBottom.setFunctionNumberOfTerms(n*8 + 5, 2)
            eftRVFreeWallInnerBottom.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS2, 1)
            eftRVFreeWallInnerBottom.setTermScaling(n*8 + 5, 1, [n*2 + 1])
            eftRVFreeWallInnerBottom.setTermNodeParameter(n*8 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftRVFreeWallInnerBottom.setTermScaling(n*8 + 5, 2, [n*2 + 2])

        # RV free wall inner side 2 elements
        eftRVFreeWallInnerSide1 = tricubichermite.createEftBasic()
        eftRVFreeWallInnerSide1.setNumberOfLocalScaleFactors(5)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eftRVFreeWallInnerSide1.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eftRVFreeWallInnerSide1.setScaleFactorIdentifier(1, 1)
        for s in range(4):
            eftRVFreeWallInnerSide1.setScaleFactorType(s + 2, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eftRVFreeWallInnerSide1.setScaleFactorIdentifier(s + 2, (s % 2) + 3)  # allow for node sf ids 1,2 on inside elements
        # map d/dxi1 = d/ds3, d/dxi3 = -a.d/ds1 + b.d/ds3 on RV apex
        for s in range(2):
            n = s*2
            ln = n + 1
            eftRVFreeWallInnerSide1.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            # 2 terms for d/dxi3 via general linear map
            eftRVFreeWallInnerSide1.setFunctionNumberOfTerms(n*8 + 5, 2)
            eftRVFreeWallInnerSide1.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftRVFreeWallInnerSide1.setTermScaling(n*8 + 5, 1, [1, s*2 + 2])
            eftRVFreeWallInnerSide1.setTermNodeParameter(n*8 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftRVFreeWallInnerSide1.setTermScaling(n*8 + 5, 2, [s*2 + 3])
        print('eftRVFreeWallInnerSide1.isValid()',eftRVFreeWallInnerSide1.isValid())

        # RV free wall inner side 2 elements
        eftRVFreeWallInnerSide2 = tricubichermite.createEftBasic()
        eftRVFreeWallInnerSide2.setNumberOfLocalScaleFactors(5)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eftRVFreeWallInnerSide2.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eftRVFreeWallInnerSide2.setScaleFactorIdentifier(1, 1)
        for s in range(4):
            eftRVFreeWallInnerSide2.setScaleFactorType(s + 2, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eftRVFreeWallInnerSide2.setScaleFactorIdentifier(s + 2, (s % 2) + 3)
        # map d/dxi1 = -d/ds3, d/dxi3 = a.d/ds1 + b.d/ds3 on RV apex
        for s in range(2):
            n = s*2 + 1
            ln = n + 1
            eftRVFreeWallInnerSide2.setTermNodeParameter(n*8 + 2, 1, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftRVFreeWallInnerSide2.setTermScaling(n*8 + 2, 1, [1])
            # 2 terms for d/dxi3 via general linear map
            eftRVFreeWallInnerSide2.setFunctionNumberOfTerms(n*8 + 5, 2)
            eftRVFreeWallInnerSide2.setTermNodeParameter(n*8 + 5, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
            eftRVFreeWallInnerSide2.setTermScaling(n*8 + 5, 1, [s*2 + 2])
            eftRVFreeWallInnerSide2.setTermNodeParameter(n*8 + 5, 2, ln, Node.VALUE_LABEL_D_DS3, 1)
            eftRVFreeWallInnerSide2.setTermScaling(n*8 + 5, 2, [s*2 + 3])
        print('eftRVFreeWallInnerSide2.isValid()',eftRVFreeWallInnerSide2.isValid())

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)
        elementtemplateRVFreeWallInnerBottom = mesh.createElementtemplate()
        elementtemplateRVFreeWallInnerBottom.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVFreeWallInnerBottom.defineField(coordinates, -1, eftRVFreeWallInnerBottom)
        elementtemplateRVFreeWallInnerSide1 = mesh.createElementtemplate()
        elementtemplateRVFreeWallInnerSide1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVFreeWallInnerSide1.defineField(coordinates, -1, eftRVFreeWallInnerSide1)
        elementtemplateRVFreeWallInnerSide2 = mesh.createElementtemplate()
        elementtemplateRVFreeWallInnerSide2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVFreeWallInnerSide2.defineField(coordinates, -1, eftRVFreeWallInnerSide2)

        cache = fm.createFieldcache()

        edgeOffsetInner = RVWidth*0.25
        edgeOffsetOuter = edgeOffsetInner*2.0 + RVFreeWallThickness

        # create RV nodes and modify adjoining LV nodes
        elementsCountUpRV = elementsCountUp - elementsCountBelowSeptum
        now = 1 + elementsCountUp*elementsCountAround
        nodeIdentifier = 1 + 2*elementsCountThroughLVWall*now
        rv_nidsInner = []
        rv_nidsOuter = []
        for n3 in range(1, -1, -1):  # only 1 element through RV free wall

            onInside = n3 == 0

            for n2 in range(elementsCountUpRV + 1):

                onBottomEdge = (n2 == 0)

                for n1 in range(elementsCountAcrossSeptum + 1):

                    onSideEdge = (n1 == 0) or (n1 == elementsCountAcrossSeptum)

                    nid = 3 + elementsCountThroughLVWall*now + (elementsCountBelowSeptum + n2 - 1)*elementsCountAround + n1

                    baseNode = nodes.findNodeByIdentifier(nid)
                    cache.setNode(baseNode)
                    result, base_x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                    result, base_dx_ds1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                    result, base_dx_ds2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                    result, base_dx_ds3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
                    #print('node ', nid, 'dx_ds3', result, base_dx_ds3)
                    mag = math.sqrt(base_dx_ds3[0]*base_dx_ds3[0] + base_dx_ds3[1]*base_dx_ds3[1] + base_dx_ds3[2]*base_dx_ds3[2])
                    unitOutward = [ base_dx_ds3[0]/mag, base_dx_ds3[1]/mag, base_dx_ds3[2]/mag ]
                    baseRadiusAround = unitOutward[0]*base_x[0] + unitOutward[1]*base_x[1] + unitOutward[2]*base_x[2]
                    if onInside:
                        if onBottomEdge or onSideEdge:
                            offset = edgeOffsetInner
                        else:
                            offset = RVWidth
                    else:
                        if onBottomEdge or onSideEdge:
                            offset = edgeOffsetOuter
                        else:
                            offset = RVWidth + RVFreeWallThickness
                    x = [
                        base_x[0] + unitOutward[0]*offset,
                        base_x[1] + unitOutward[1]*offset,
                        base_x[2] + unitOutward[2]*offset,
                    ]
                    scale1 = (baseRadiusAround + offset)/baseRadiusAround
                    dx_ds1 = [
                        base_dx_ds1[0]*scale1,
                        base_dx_ds1[1]*scale1,
                        base_dx_ds1[2]*scale1
                    ]
                    # GRC not sure if appropriate to use scale1 here:
                    dx_ds2 = [
                        base_dx_ds2[0]*scale1,
                        base_dx_ds2[1]*scale1,
                        base_dx_ds2[2]*scale1
                    ]
                    dx_ds3 = [
                        unitOutward[0]*RVFreeWallThickness,
                        unitOutward[1]*RVFreeWallThickness,
                        unitOutward[2]*RVFreeWallThickness
                    ]

                    if onInside and (onBottomEdge or onSideEdge):
                        node = baseNode
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        nodeIdentifier += 1
                        cache.setNode(node)
                    if onInside:
                        rv_nidsInner.append(node.getIdentifier())
                    else:
                        rv_nidsOuter.append(node.getIdentifier())
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)

        # create RV elements and modify adjoining LV element fields
        elementIdentifier = elementsCountThroughLVWall*elementsCountUp*elementsCountAround + 1
        crossAngle = math.pi/12
        sinCrossAngle = math.sin(crossAngle)
        cosCrossAngle = math.cos(crossAngle)
        scalefactors4 = [ sinCrossAngle, cosCrossAngle, sinCrossAngle, cosCrossAngle ]
        scalefactors4r2 = [ -sinCrossAngle, cosCrossAngle, -sinCrossAngle, cosCrossAngle ]
        scalefactors5 = [ -1.0, sinCrossAngle, cosCrossAngle, sinCrossAngle, cosCrossAngle ]

        for n2 in range(elementsCountUpRV):

            for n1 in range(elementsCountAcrossSeptum):

                bni = n2*(elementsCountAcrossSeptum + 1) + n1
                nodeIdentifiers = [
                    rv_nidsInner[bni], rv_nidsInner[bni + 1], rv_nidsInner[bni + elementsCountAcrossSeptum + 1], rv_nidsInner[bni + elementsCountAcrossSeptum + 2],
                    rv_nidsOuter[bni], rv_nidsOuter[bni + 1], rv_nidsOuter[bni + elementsCountAcrossSeptum + 1], rv_nidsOuter[bni + elementsCountAcrossSeptum + 2]
                ]
                if n2 == 0:
                    element = mesh.createElement(elementIdentifier, elementtemplateRVFreeWallInnerBottom)
                    result = element.setNodesByIdentifier(eftRVFreeWallInnerBottom, nodeIdentifiers)
                    result = element.setScaleFactors(eftRVFreeWallInnerBottom, scalefactors4r2)
                elif n1 == 0:
                    element = mesh.createElement(elementIdentifier, elementtemplateRVFreeWallInnerSide1)
                    result = element.setNodesByIdentifier(eftRVFreeWallInnerSide1, nodeIdentifiers)
                    result = element.setScaleFactors(eftRVFreeWallInnerSide1, scalefactors5)
                elif n1 == (elementsCountAcrossSeptum - 1):
                    element = mesh.createElement(elementIdentifier, elementtemplateRVFreeWallInnerSide2)
                    result = element.setNodesByIdentifier(eftRVFreeWallInnerSide2, nodeIdentifiers)
                    result = element.setScaleFactors(eftRVFreeWallInnerSide2, scalefactors5)
                else:
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                print('element ', elementIdentifier, result)
                elementIdentifier += 1

        fm.endChange()
