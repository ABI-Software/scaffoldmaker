"""
Generates a 2-D unit plate mesh with variable numbers of elements in 2 directions.
"""

from __future__ import division

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base


class MeshType_2d_plate1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '2D Plate 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Coordinate dimensions' : 3,
            'Number of elements 1' : 1,
            'Number of elements 2' : 1,
            'Use cross derivatives' : False
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Coordinate dimensions',
            'Number of elements 1',
            'Number of elements 2',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Coordinate dimensions'] < 2) :
            options['Coordinate dimensions'] = 2
        elif (options['Coordinate dimensions'] > 3) :
            options['Coordinate dimensions'] = 3
        if (options['Number of elements 1'] < 1) :
            options['Number of elements 1'] = 1
        if (options['Number of elements 2'] < 1) :
            options['Number of elements 2'] = 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        coordinateDimensions = options['Coordinate dimensions']
        elementsCount1 = options['Number of elements 1']
        elementsCount2 = options['Number of elements 2']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm, components_count=coordinateDimensions)

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
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 1.0 / elementsCount1, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 1.0 / elementsCount2, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        for n2 in range(elementsCount2 + 1):
            x[1] = n2 / elementsCount2
            for n1 in range(elementsCount1 + 1):
                x[0] = n1 / elementsCount1
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
        no2 = (elementsCount1 + 1)
        for e2 in range(elementsCount2):
            for e1 in range(elementsCount1):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni = e2*no2 + e1 + 1
                nodeIdentifiers = [ bni, bni + 1, bni + no2, bni + no2 + 1 ]
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

        fm.endChange()
        return []

