"""
Generates 3-D Left and Right ventricles mesh starting from modified sphere shell mesh.
"""

import math
from mapclientplugins.meshgeneratorstep.meshtypes.meshtype_3d_sphereshell1 import MeshType_3d_sphereshell1
from mapclientplugins.meshgeneratorstep.utils.eft_utils import *
from mapclientplugins.meshgeneratorstep.utils.zinc_utils import *
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
        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        # RV septal wall inner bottom elements
        eftRVSeptalWallInnerBottom = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVSeptalWallInnerBottom, [1], [3, 4, 3, 4])
        for s in range(2):
            n = 4 + s
            ln = n + 1
            # d/dxi2 = -d/ds3, d/dxi3 = sa.d/ds2 + sb.d/ds3
            mapEftFunction1Node1Term(eftRVSeptalWallInnerBottom, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
            mapEftFunction1Node2Terms(eftRVSeptalWallInnerBottom, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [s*2 + 2], Node.VALUE_LABEL_D_DS3, 1, [s*2 + 3])
        elementtemplateRVSeptalWallInnerBottom = mesh.createElementtemplate()
        elementtemplateRVSeptalWallInnerBottom.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVSeptalWallInnerBottom.defineField(coordinates, -1, eftRVSeptalWallInnerBottom)
        print('eftRVSeptalWallInnerBottom', result)

        # RV septal wall inner side 1 elements
        eftRVSeptalWallInnerSide1 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVSeptalWallInnerSide1, [1], [3, 4, 3, 4])
        for s in range(2):
            n = 4 + s*2
            ln = n + 1
            # d/dxi1 = -d/ds3, d/dxi3 = sa.d/ds1 + sb.d/ds3
            mapEftFunction1Node1Term(eftRVSeptalWallInnerSide1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
            mapEftFunction1Node2Terms(eftRVSeptalWallInnerSide1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [s*2 + 2], Node.VALUE_LABEL_D_DS3, 1, [s*2 + 3])
        elementtemplateRVSeptalWallInnerSide1 = mesh.createElementtemplate()
        elementtemplateRVSeptalWallInnerSide1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVSeptalWallInnerSide1.defineField(coordinates, -1, eftRVSeptalWallInnerSide1)
        print('eftRVSeptalWallInnerSide1', result)

        # RV septal wall inner side 2 elements
        eftRVSeptalWallInnerSide2 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVSeptalWallInnerSide2, [1], [3, 4, 3, 4])
        for s in range(2):
            n = 5 + s*2
            ln = n + 1
            # d/dxi1 = d/ds3, d/dxi3 = -sa.d/ds1 + sb.d/ds3
            mapEftFunction1Node1Term(eftRVSeptalWallInnerSide2, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [])
            mapEftFunction1Node2Terms(eftRVSeptalWallInnerSide2, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [1, s*2 + 2], Node.VALUE_LABEL_D_DS3, 1, [s*2 + 3])
        elementtemplateRVSeptalWallInnerSide2 = mesh.createElementtemplate()
        elementtemplateRVSeptalWallInnerSide2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVSeptalWallInnerSide2.defineField(coordinates, -1, eftRVSeptalWallInnerSide2)
        print('eftRVSeptalWallInnerSide2', result)

        # RV free wall inner bottom elements
        eftRVFreeWallInnerBottom = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVFreeWallInnerBottom, [1], [3, 4, 3, 4])
        for s in range(2):
            n = s
            ln = n + 1
            # d/dxi2 = d/ds3, d/dxi3 = -sa.d/ds2 + sb.d/ds3
            mapEftFunction1Node1Term(eftRVFreeWallInnerBottom, n*8 + 3, ln, Node.VALUE_LABEL_D_DS3, 1, [])
            mapEftFunction1Node2Terms(eftRVFreeWallInnerBottom, n*8 + 5, ln, Node.VALUE_LABEL_D_DS2, 1, [1, s*2 + 2], Node.VALUE_LABEL_D_DS3, 1, [s*2 + 3])
        elementtemplateRVFreeWallInnerBottom = mesh.createElementtemplate()
        elementtemplateRVFreeWallInnerBottom.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVFreeWallInnerBottom.defineField(coordinates, -1, eftRVFreeWallInnerBottom)
        print('eftRVFreeWallInnerBottom', result)

        # RV free wall inner side 1 elements
        eftRVFreeWallInnerSide1 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVFreeWallInnerSide1, [1], [3, 4, 3, 4])
        for s in range(2):
            n = s*2
            ln = n + 1
            # d/dxi1 = d/ds3, d/dxi3 = -sa.d/ds1 + sb.d/ds3
            mapEftFunction1Node1Term(eftRVFreeWallInnerSide1, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [])
            mapEftFunction1Node2Terms(eftRVFreeWallInnerSide1, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [1, s*2 + 2], Node.VALUE_LABEL_D_DS3, 1, [s*2 + 3])
        elementtemplateRVFreeWallInnerSide1 = mesh.createElementtemplate()
        elementtemplateRVFreeWallInnerSide1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVFreeWallInnerSide1.defineField(coordinates, -1, eftRVFreeWallInnerSide1)
        print('eftRVFreeWallInnerSide1', result)

        # RV free wall inner side 2 elements
        eftRVFreeWallInnerSide2 = tricubichermite.createEftBasic()
        setEftScaleFactorIds(eftRVFreeWallInnerSide2, [1], [3, 4, 3, 4])
        for s in range(2):
            n = s*2 + 1
            ln = n + 1
            # d/dxi1 = -d/ds3, d/dxi3 = sa.d/ds1 + sb.d/ds3
            mapEftFunction1Node1Term(eftRVFreeWallInnerSide2, n*8 + 2, ln, Node.VALUE_LABEL_D_DS3, 1, [1])
            mapEftFunction1Node2Terms(eftRVFreeWallInnerSide2, n*8 + 5, ln, Node.VALUE_LABEL_D_DS1, 1, [s*2 + 2], Node.VALUE_LABEL_D_DS3, 1, [s*2 + 3])
        elementtemplateRVFreeWallInnerSide2 = mesh.createElementtemplate()
        elementtemplateRVFreeWallInnerSide2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplateRVFreeWallInnerSide2.defineField(coordinates, -1, eftRVFreeWallInnerSide2)
        print('eftRVFreeWallInnerSide2', result)

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
                    scale3 = RVFreeWallThickness
                    if onInside and (onBottomEdge or onSideEdge):
                        scale3 += 2.0*edgeOffsetInner
                    dx_ds3 = [
                        unitOutward[0]*scale3,
                        unitOutward[1]*scale3,
                        unitOutward[2]*scale3
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
        scalefactors5 = [ -1.0, sinCrossAngle, cosCrossAngle, sinCrossAngle, cosCrossAngle ]

        # Tweak RV septal wall
        eow = elementsCountUp*elementsCountAround
        RVSeptumElementIdBase = eow*(elementsCountThroughLVWall - 1) + elementsCountBelowSeptum*elementsCountAround + 2
        for n2 in range(elementsCountUpRV):

            for n1 in range(elementsCountAcrossSeptum):
                existingElementIdentifier = RVSeptumElementIdBase + n2*elementsCountAround + n1

                element = mesh.findElementByIdentifier(existingElementIdentifier)
                eftTemp = element.getElementfieldtemplate(coordinates, -1)
                nodeIdentifiers = getElementNodeIdentifiers(element, eft)
                #print('existing element', existingElementIdentifier, eftTemp.isValid(), nodeIdentifiers)
                if n2 == 0:
                    if n1 == 0:
                        # RV free wall inner bottom side 1 elements: doubly curved but can't avoid a sharp corner at node 1
                        eftTemp = tricubichermite.createEftBasic()
                        setEftScaleFactorIds(eftTemp, [1], [3, 4, 3, 4])
                        # node 6 d/dxi2 = -d/ds3, d/dxi3 = sa.d/ds2 + sb.d/ds3
                        mapEftFunction1Node1Term(eftTemp, 5*8 + 3, 6, Node.VALUE_LABEL_D_DS3, 1, [1])
                        mapEftFunction1Node2Terms(eftTemp, 5*8 + 5, 6, Node.VALUE_LABEL_D_DS2, 1, [2], Node.VALUE_LABEL_D_DS3, 1, [3])
                        # node 7 d/dxi1 = -d/ds3, d/dxi3 = sa.d/ds1 + sb.d/ds3
                        mapEftFunction1Node1Term(eftTemp, 6*8 + 2, 7, Node.VALUE_LABEL_D_DS3, 1, [1])
                        mapEftFunction1Node2Terms(eftTemp, 6*8 + 5, 7, Node.VALUE_LABEL_D_DS1, 1, [1, 4], Node.VALUE_LABEL_D_DS3, 1, [5])
                        #print('inner bottom side 1 eftTemp.isValid()',eftTemp.isValid())
                        elementtemplateTemp = mesh.createElementtemplate()
                        elementtemplateTemp.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                        result = elementtemplateTemp.defineField(coordinates, -1, eftTemp)
                        result = element.merge(elementtemplateTemp)
                        result = element.setNodesByIdentifier(eftTemp, nodeIdentifiers)
                        result = element.setScaleFactors(eftTemp, scalefactors5)
                    elif n1 == (elementsCountAcrossSeptum - 1):
                        # RV free wall inner bottom side 1 elements: doubly curved but can't avoid a sharp corner at node 1
                        eftTemp = tricubichermite.createEftBasic()
                        setEftScaleFactorIds(eftTemp, [1], [3, 4, 3, 4])
                        # node 5 d/dxi2 = -d/ds3, d/dxi3 = sa.d/ds2 + sb.d/ds3
                        mapEftFunction1Node1Term(eftTemp, 4*8 + 3, 5, Node.VALUE_LABEL_D_DS3, 1, [1])
                        mapEftFunction1Node2Terms(eftTemp, 4*8 + 5, 5, Node.VALUE_LABEL_D_DS2, 1, [2], Node.VALUE_LABEL_D_DS3, 1, [3])
                        # node 8 d/dxi1 = d/ds3, d/dxi3 = s2.d/ds1 + s3.d/ds3
                        mapEftFunction1Node1Term(eftTemp, 7*8 + 2, 8, Node.VALUE_LABEL_D_DS3, 1, [])
                        mapEftFunction1Node2Terms(eftTemp, 7*8 + 5, 8, Node.VALUE_LABEL_D_DS1, 1, [1, 4], Node.VALUE_LABEL_D_DS3, 1, [5])
                        #print('inner bottom side 2 eftTemp.isValid()',eftTemp.isValid())
                        elementtemplateTemp = mesh.createElementtemplate()
                        elementtemplateTemp.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                        result = elementtemplateTemp.defineField(coordinates, -1, eftTemp)
                        result = element.merge(elementtemplateTemp)
                        result = element.setNodesByIdentifier(eftTemp, nodeIdentifiers)
                        result = element.setScaleFactors(eftTemp, scalefactors5)
                    else:
                        result = element.merge(elementtemplateRVSeptalWallInnerBottom)
                        result = element.setNodesByIdentifier(eftRVSeptalWallInnerBottom, nodeIdentifiers)
                        result = element.setScaleFactors(eftRVSeptalWallInnerBottom, scalefactors5)
                else:
                    if n1 == 0:
                        result = element.merge(elementtemplateRVSeptalWallInnerSide1)
                        result = element.setNodesByIdentifier(eftRVSeptalWallInnerSide1, nodeIdentifiers)
                        result = element.setScaleFactors(eftRVSeptalWallInnerSide1, scalefactors5)
                    elif n1 == (elementsCountAcrossSeptum - 1):
                        result = element.merge(elementtemplateRVSeptalWallInnerSide2)
                        result = element.setNodesByIdentifier(eftRVSeptalWallInnerSide2, nodeIdentifiers)
                        result = element.setScaleFactors(eftRVSeptalWallInnerSide2, scalefactors5)
                    else:
                        pass
                print('RV septal wall element merge', existingElementIdentifier, result)

        # RV free wall
        for n2 in range(elementsCountUpRV):

            for n1 in range(elementsCountAcrossSeptum):

                bni = n2*(elementsCountAcrossSeptum + 1) + n1
                nodeIdentifiers = [
                    rv_nidsInner[bni], rv_nidsInner[bni + 1], rv_nidsInner[bni + elementsCountAcrossSeptum + 1], rv_nidsInner[bni + elementsCountAcrossSeptum + 2],
                    rv_nidsOuter[bni], rv_nidsOuter[bni + 1], rv_nidsOuter[bni + elementsCountAcrossSeptum + 1], rv_nidsOuter[bni + elementsCountAcrossSeptum + 2]
                ]
                if n2 == 0:
                    if n1 == 0:
                        # RV free wall inner bottom side 1 elements: doubly curved but can't avoid a sharp corner at node 1
                        eftTemp = tricubichermite.createEftBasic()
                        setEftScaleFactorIds(eftTemp, [1], [3, 4, 3, 4])
                        # node 2 d/dxi2 = d/ds3, d/dxi3 = -sa.d/ds2 + sb.d/ds3
                        mapEftFunction1Node1Term(eftTemp, 1*8 + 3, 2, Node.VALUE_LABEL_D_DS3, 1, [])
                        mapEftFunction1Node2Terms(eftTemp, 1*8 + 5, 2, Node.VALUE_LABEL_D_DS2, 1, [1, 2], Node.VALUE_LABEL_D_DS3, 1, [3])
                        # node 3 d/dxi1 = d/ds3, d/dxi3 = -sa.d/ds1 + sb.d/ds3
                        mapEftFunction1Node1Term(eftTemp, 2*8 + 2, 3, Node.VALUE_LABEL_D_DS3, 1, [])
                        mapEftFunction1Node2Terms(eftTemp, 2*8 + 5, 3, Node.VALUE_LABEL_D_DS1, 1, [1, 4], Node.VALUE_LABEL_D_DS3, 1, [5])
                        #print('inner bottom side 1 eftTemp.isValid()',eftTemp.isValid())
                        elementtemplateTemp = mesh.createElementtemplate()
                        elementtemplateTemp.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                        result = elementtemplateTemp.defineField(coordinates, -1, eftTemp)
                        element = mesh.createElement(elementIdentifier, elementtemplateTemp)
                        result = element.setNodesByIdentifier(eftTemp, nodeIdentifiers)
                        result = element.setScaleFactors(eftTemp, scalefactors5)
                    elif n1 == (elementsCountAcrossSeptum - 1):
                        # RV free wall inner bottom side 1 elements: doubly curved but can't avoid a sharp corner at node 1
                        eftTemp = tricubichermite.createEftBasic()
                        setEftScaleFactorIds(eftTemp, [1], [3, 4, 3, 4])
                        # node 1 d/dxi2 = d/ds3, d/dxi3 = -sa.d/ds2 + sb.d/ds3
                        mapEftFunction1Node1Term(eftTemp, 0*8 + 3, 1, Node.VALUE_LABEL_D_DS3, 1, [])
                        mapEftFunction1Node2Terms(eftTemp, 0*8 + 5, 1, Node.VALUE_LABEL_D_DS2, 1, [1, 2], Node.VALUE_LABEL_D_DS3, 1, [3])
                        # node 4 d/dxi1 = -d/ds3, d/dxi3 = sa.d/ds1 + sb.d/ds3
                        mapEftFunction1Node1Term(eftTemp, 3*8 + 2, 4, Node.VALUE_LABEL_D_DS3, 1, [1])
                        mapEftFunction1Node2Terms(eftTemp, 3*8 + 5, 4, Node.VALUE_LABEL_D_DS1, 1, [4], Node.VALUE_LABEL_D_DS3, 1, [5])
                        #print('inner bottom side 2 eftTemp.isValid()',eftTemp.isValid())
                        elementtemplateTemp = mesh.createElementtemplate()
                        elementtemplateTemp.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                        result = elementtemplateTemp.defineField(coordinates, -1, eftTemp)
                        element = mesh.createElement(elementIdentifier, elementtemplateTemp)
                        result = element.setNodesByIdentifier(eftTemp, nodeIdentifiers)
                        result = element.setScaleFactors(eftTemp, scalefactors5)
                    else:
                        element = mesh.createElement(elementIdentifier, elementtemplateRVFreeWallInnerBottom)
                        result = element.setNodesByIdentifier(eftRVFreeWallInnerBottom, nodeIdentifiers)
                        result = element.setScaleFactors(eftRVFreeWallInnerBottom, scalefactors5)
                else:
                    if n1 == 0:
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
                print('RV free wall element create', elementIdentifier, result)
                elementIdentifier += 1

        fm.endChange()
