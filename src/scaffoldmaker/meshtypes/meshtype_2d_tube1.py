"""
Generates a 2-D unit tube mesh with variable numbers of elements around, along.
"""

from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base


class MeshType_2d_tube1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '2D Tube 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements along' : 1,
            'Number of elements around' : 4,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along',
            'Number of elements around',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Number of elements along'] < 1) :
            options['Number of elements along'] = 1
        if (options['Number of elements around'] < 2) :
            options['Number of elements around'] = 2

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        elementsCountAlong = options['Number of elements along']
        elementsCountAround = options['Number of elements around']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

        mesh = fm.findMeshByDimension(2)
        bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
        if not useCrossDerivatives:
            for n in range(4):
                eft.setFunctionNumberOfTerms(n*4 + 4, 0)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 1.0 / elementsCountAlong ]
        zero = [ 0.0, 0.0, 0.0 ]
        radius = 0.5
        for n2 in range(elementsCountAlong + 1):
            x[2] = n2 / elementsCountAlong
            for n1 in range(elementsCountAround):
                radiansAround = n1*radiansPerElementAround
                cosRadiansAround = math.cos(radiansAround)
                sinRadiansAround = math.sin(radiansAround)
                x[0] = radius*cosRadiansAround
                x[1] = radius*sinRadiansAround
                dx_ds1[0] = radiansPerElementAround*radius*-sinRadiansAround
                dx_ds1[1] = radiansPerElementAround*radius*cosRadiansAround
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                if useCrossDerivatives:
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        for e2 in range(elementsCountAlong):
            for e1 in range(elementsCountAround):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni1 = e2*elementsCountAround + e1 + 1
                bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

        fm.endChange()
        return []
