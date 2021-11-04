"""
Generates a 2-D unit sphere mesh with variable numbers of elements around and up.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base


class MeshType_2d_sphere1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '2D Sphere 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements up' : 4,
            'Number of elements around' : 4,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up',
            'Number of elements around',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Number of elements up'] < 2) :
            options['Number of elements up'] = 2
        if (options['Number of elements around'] < 2) :
            options['Number of elements around'] = 2

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        elementsCountUp = options['Number of elements up']
        elementsCountAround = options['Number of elements around']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplateApex = nodes.createNodetemplate()
        nodetemplateApex.defineField(coordinates)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        else:
            nodetemplate = nodetemplateApex


        mesh = fm.findMeshByDimension(2)
        bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
        if not useCrossDerivatives:
            for n in range(4):
                eft.setFunctionNumberOfTerms(n*4 + 4, 0)

        # Apex1: collapsed on xi2 = 0
        eftApex1 = mesh.createElementfieldtemplate(bicubicHermiteBasis)
        eftApex1.setNumberOfLocalNodes(3)
        eftApex1.setNumberOfLocalScaleFactors(7)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eftApex1.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eftApex1.setScaleFactorIdentifier(1, 1)
        for s in range(6):
            si = s + 2
            # 3 scale factors per node: cos(theta), sin(theta), arc angle radians
            sid = (s // 3)*100 + (s % 3) + 1  # add 100 for different 'version'
            eftApex1.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eftApex1.setScaleFactorIdentifier(si, sid)
        # basis nodes 1, 2 -> local node 1
        ln = 1
        eftApex1.setTermNodeParameter(1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
        # 0 terms = zero parameter for d/dxi1 basis
        eftApex1.setFunctionNumberOfTerms(2, 0)
        # 2 terms for d/dxi2 via general linear map:
        eftApex1.setFunctionNumberOfTerms(3, 2)
        eftApex1.setTermNodeParameter(3, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
        eftApex1.setTermScaling(3, 1, [2])
        eftApex1.setTermNodeParameter(3, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
        eftApex1.setTermScaling(3, 2, [3])
        # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
        eftApex1.setFunctionNumberOfTerms(4, 2)
        eftApex1.setTermNodeParameter(4, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
        eftApex1.setTermScaling(4, 1, [3, 4])
        eftApex1.setTermNodeParameter(4, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
        eftApex1.setTermScaling(4, 2, [1, 2, 4])
        # basis node 2 -> local node 1
        eftApex1.setTermNodeParameter(5, 1, ln, Node.VALUE_LABEL_VALUE, 1)
        # 0 terms = zero parameter for d/dxi1 basis
        eftApex1.setFunctionNumberOfTerms(6, 0)
        # 2 terms for d/dxi2 via general linear map:
        eftApex1.setFunctionNumberOfTerms(7, 2)
        eftApex1.setTermNodeParameter(7, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
        eftApex1.setTermScaling(7, 1, [5])
        eftApex1.setTermNodeParameter(7, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
        eftApex1.setTermScaling(7, 2, [6])
        # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
        eftApex1.setFunctionNumberOfTerms(8, 2)
        eftApex1.setTermNodeParameter(8, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
        eftApex1.setTermScaling(8, 1, [6, 7])
        eftApex1.setTermNodeParameter(8, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
        eftApex1.setTermScaling(8, 2, [1, 5, 7])

        # basis nodes 3, 4 -> regular local nodes 2, 3
        for bn in range(2,4):
            fo = bn*4
            ni = bn
            eftApex1.setTermNodeParameter(fo + 1, 1, ni, Node.VALUE_LABEL_VALUE, 1)
            eftApex1.setTermNodeParameter(fo + 2, 1, ni, Node.VALUE_LABEL_D_DS1, 1)
            eftApex1.setTermNodeParameter(fo + 3, 1, ni, Node.VALUE_LABEL_D_DS2, 1)
            if useCrossDerivatives:
                eftApex1.setTermNodeParameter(fo + 4, 1, ni, Node.VALUE_LABEL_D2_DS1DS2, 1)
            else:
                eftApex1.setFunctionNumberOfTerms(fo + 4, 0)

        # Apex2: collapsed on xi2 = 1
        eftApex2 = mesh.createElementfieldtemplate(bicubicHermiteBasis)
        eftApex2.setNumberOfLocalNodes(3)
        eftApex2.setNumberOfLocalScaleFactors(7)
        # GRC: allow scale factor identifier for global -1.0 to be prescribed
        eftApex2.setScaleFactorType(1, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eftApex2.setScaleFactorIdentifier(1, 1)
        for s in range(6):
            si = s + 2
            # 3 scale factors per node: cos(theta), sin(theta), arc angle radians
            sid = (s // 3)*100 + (s % 3) + 1  # add 100 for different 'version'
            eftApex2.setScaleFactorType(si, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
            eftApex2.setScaleFactorIdentifier(si, sid)
        # basis nodes 1, 2 -> regular local nodes 1, 2 (for each layer)
        for bn in range(2):
            fo = bn*4
            ni = bn + 1
            eftApex2.setTermNodeParameter(fo + 1, 1, ni, Node.VALUE_LABEL_VALUE, 1)
            eftApex2.setTermNodeParameter(fo + 2, 1, ni, Node.VALUE_LABEL_D_DS1, 1)
            eftApex2.setTermNodeParameter(fo + 3, 1, ni, Node.VALUE_LABEL_D_DS2, 1)
            if useCrossDerivatives:
                eftApex2.setTermNodeParameter(fo + 4, 1, ni, Node.VALUE_LABEL_D2_DS1DS2, 1)
            else:
                eftApex2.setFunctionNumberOfTerms(fo + 4, 0)

        # basis nodes 3, 4 -> local node 3
        ln = 3
        fo3 = 8
        eftApex2.setTermNodeParameter(fo3 + 1, 1, ln, Node.VALUE_LABEL_VALUE, 1)
        # 0 terms = zero parameter for d/dxi1 basis
        eftApex2.setFunctionNumberOfTerms(fo3 + 2, 0)
        # 2 terms for d/dxi2 via general linear map:
        eftApex2.setFunctionNumberOfTerms(fo3 + 3, 2)
        eftApex2.setTermNodeParameter(fo3 + 3, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
        eftApex2.setTermScaling(fo3 + 3, 1, [2])
        eftApex2.setTermNodeParameter(fo3 + 3, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
        eftApex2.setTermScaling(fo3 + 3, 2, [3])
        # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
        eftApex2.setFunctionNumberOfTerms(fo3 + 4, 2)
        eftApex2.setTermNodeParameter(fo3 + 4, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
        eftApex2.setTermScaling(fo3 + 4, 1, [1, 3, 4])
        eftApex2.setTermNodeParameter(fo3 + 4, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
        eftApex2.setTermScaling(fo3 + 4, 2, [2, 4])
        # basis node 4 -> local node 3
        eftApex2.setTermNodeParameter(fo3 + 5, 1, ln, Node.VALUE_LABEL_VALUE, 1)
        # 0 terms = zero parameter for d/dxi1 basis
        eftApex2.setFunctionNumberOfTerms(fo3 + 6, 0)
        # 2 terms for d/dxi2 via general linear map:
        eftApex2.setFunctionNumberOfTerms(fo3 + 7, 2)
        eftApex2.setTermNodeParameter(fo3 + 7, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
        eftApex2.setTermScaling(fo3 + 7, 1, [5])
        eftApex2.setTermNodeParameter(fo3 + 7, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
        eftApex2.setTermScaling(fo3 + 7, 2, [6])
        # 2 terms for cross derivative 1 2 to correct circular apex: -sin(theta).phi, cos(theta).phi
        eftApex2.setFunctionNumberOfTerms(fo3 + 8, 2)
        eftApex2.setTermNodeParameter(fo3 + 8, 1, ln, Node.VALUE_LABEL_D_DS1, 1)
        eftApex2.setTermScaling(fo3 + 8, 1, [1, 6, 7])
        eftApex2.setTermNodeParameter(fo3 + 8, 2, ln, Node.VALUE_LABEL_D_DS2, 1)
        eftApex2.setTermScaling(fo3 + 8, 2, [5, 7])

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        elementtemplate.defineField(coordinates, -1, eft)
        elementtemplateApex1 = mesh.createElementtemplate()
        elementtemplateApex1.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        elementtemplateApex1.defineField(coordinates, -1, eftApex1)
        elementtemplateApex2 = mesh.createElementtemplate()
        elementtemplateApex2.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        elementtemplateApex2.defineField(coordinates, -1, eftApex2)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementUp = math.pi/elementsCountUp
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        radius = 0.5

        # create apex1 node
        node = nodes.createNode(nodeIdentifier, nodetemplateApex)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, -radius ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, radius*radiansPerElementUp, 0.0 ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
        nodeIdentifier = nodeIdentifier + 1

        # create regular rows between apexes
        for n2 in range(1, elementsCountUp):
            radiansUp = n2*radiansPerElementUp
            cosRadiansUp = math.cos(radiansUp);
            sinRadiansUp = math.sin(radiansUp);
            for n1 in range(elementsCountAround):
                radiansAround = n1*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x = [
                    radius*cosRadiansAround*sinRadiansUp,
                    radius*sinRadiansAround*sinRadiansUp,
                    -radius*cosRadiansUp
                ]
                dx_ds1 = [
                    radius*-sinRadiansAround*sinRadiansUp*radiansPerElementAround,
                    radius*cosRadiansAround*sinRadiansUp*radiansPerElementAround,
                    0.0
                ]
                dx_ds2 = [
                    radius*cosRadiansAround*cosRadiansUp*radiansPerElementUp,
                    radius*sinRadiansAround*cosRadiansUp*radiansPerElementUp,
                    radius*sinRadiansUp*radiansPerElementUp
                ]
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                nodeIdentifier = nodeIdentifier + 1

        # create apex2 node
        node = nodes.createNode(nodeIdentifier, nodetemplateApex)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, [ 0.0, 0.0, radius ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [ 0.0, -radius*radiansPerElementUp, 0.0 ])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [ radius*radiansPerElementUp, 0.0, 0.0 ])
        nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        # create Apex1 elements, editing eft scale factor identifiers around apex
        # scale factor identifiers follow convention of offsetting by 100 for each 'version'
        bni1 = 1
        for e1 in range(elementsCountAround):
            va = e1
            vb = (e1 + 1)%elementsCountAround
            eftApex1.setScaleFactorIdentifier(2, va*100 + 1)
            eftApex1.setScaleFactorIdentifier(3, va*100 + 2)
            eftApex1.setScaleFactorIdentifier(4, va*100 + 3)
            eftApex1.setScaleFactorIdentifier(5, vb*100 + 1)
            eftApex1.setScaleFactorIdentifier(6, vb*100 + 2)
            eftApex1.setScaleFactorIdentifier(7, vb*100 + 3)
            # redefine field in template for changes to eftApex1:
            elementtemplateApex1.defineField(coordinates, -1, eftApex1)
            element = mesh.createElement(elementIdentifier, elementtemplateApex1)
            bni2 = e1 + 2
            bni3 = (e1 + 1)%elementsCountAround + 2
            nodeIdentifiers = [ bni1, bni2, bni3 ]
            element.setNodesByIdentifier(eftApex1, nodeIdentifiers)
            # set general linear map coefficients
            radiansAround = e1*radiansPerElementAround
            radiansAroundNext = ((e1 + 1)%elementsCountAround)*radiansPerElementAround
            scalefactors = [
                -1.0,
                math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
            ]
            result = element.setScaleFactors(eftApex1, scalefactors)
            elementIdentifier = elementIdentifier + 1

        # create regular rows between apexes
        for e2 in range(elementsCountUp - 2):
            for e1 in range(elementsCountAround):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni1 = e2*elementsCountAround + e1 + 2
                bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + 2
                nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

        # create Apex2 elements, editing eft scale factor identifiers around apex
        # scale factor identifiers follow convention of offsetting by 100 for each 'version'
        bni3 = 2 + (elementsCountUp - 1)*elementsCountAround
        for e1 in range(elementsCountAround):
            va = e1
            vb = (e1 + 1)%elementsCountAround
            eftApex2.setScaleFactorIdentifier(2, va*100 + 1)
            eftApex2.setScaleFactorIdentifier(3, va*100 + 2)
            eftApex2.setScaleFactorIdentifier(4, va*100 + 3)
            eftApex2.setScaleFactorIdentifier(5, vb*100 + 1)
            eftApex2.setScaleFactorIdentifier(6, vb*100 + 2)
            eftApex2.setScaleFactorIdentifier(7, vb*100 + 3)
            # redefine field in template for changes to eftApex2:
            elementtemplateApex1.defineField(coordinates, -1, eftApex2)
            element = mesh.createElement(elementIdentifier, elementtemplateApex1)
            bni1 = bni3 - elementsCountAround + e1
            bni2 = bni3 - elementsCountAround + (e1 + 1)%elementsCountAround
            nodeIdentifiers = [ bni1, bni2, bni3 ]
            result = element.setNodesByIdentifier(eftApex2, nodeIdentifiers)
            # set general linear map coefficients
            radiansAround = math.pi + e1*radiansPerElementAround
            radiansAroundNext = math.pi + ((e1 + 1)%elementsCountAround)*radiansPerElementAround
            scalefactors = [
                -1.0,
                -math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                -math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
            ]
            result = element.setScaleFactors(eftApex2, scalefactors)
            elementIdentifier = elementIdentifier + 1

        fm.endChange()
        return []
