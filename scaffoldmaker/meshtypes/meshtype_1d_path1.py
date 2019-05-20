"""
Generates a 1-D path mesh.
"""

from __future__ import division
import math
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import zinc_utils
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_1d_path1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '1D Path 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Coordinate dimensions' : 3,
            'Length' : 1.0,
            'Number of elements' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Coordinate dimensions',
            'Length',
            'Number of elements'
        ]

    @staticmethod
    def checkOptions(options):
        if (options['Coordinate dimensions'] < 1) :
            options['Coordinate dimensions'] = 1
        elif (options['Coordinate dimensions'] > 3) :
            options['Coordinate dimensions'] = 3
        if (options['Number of elements'] < 1) :
            options['Number of elements'] = 1

    @staticmethod
    def generateMesh(region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        coordinateDimensions = options['Coordinate dimensions']
        length = options['Length']
        elementsCount = options['Number of elements']

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = zinc_utils.getOrCreateCoordinateField(fm, componentsCount=coordinateDimensions)
        cache = fm.createFieldcache()

        #################
        # Create nodes
        #################

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

        nodeIdentifier = 1
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ length/elementsCount, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 1.0, 0.0 ]
        d2x_ds1ds2 = [ 0.0, 0.0, 0.0 ]
        for n in range(elementsCount + 1):
            x[0] = length*n/elementsCount
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            nodeIdentifier = nodeIdentifier + 1

        #################
        # Create elements
        #################

        mesh = fm.findMeshByDimension(1)
        cubicHermiteBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(cubicHermiteBasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        for e in range(elementsCount):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1, e + 2 ])
            elementIdentifier = elementIdentifier + 1

        fm.endChange()


def extractPathParametersFromRegion(region):
    '''
    Returns parameters of all nodes in region in identifier order.
    Assumes nodes in region have field coordinates (1 to 3 components).
    Currently limited to nodes with exactly value, d_ds1, d_ds2, d2_ds12,
    same as path 1 scaffold.
    :return: cx, cd1, cd2, cd12 (all padded with zeroes to 3 components)
    '''
    fm = region.getFieldmodule()
    coordinates = fm.findFieldByName('coordinates').castFiniteElement()
    componentsCount = coordinates.getNumberOfComponents()
    assert componentsCount in [ 1, 2, 3 ], 'extractPathParametersFromRegion.  Invalid coordinates number of components'
    cache = fm.createFieldcache()
    cx = []
    cd1 = []
    cd2 = []
    cd12 = []
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodeIter = nodes.createNodeiterator()
    node = nodeIter.next()
    while node.isValid():
        cache.setNode(node)
        result, x  = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, componentsCount)
        result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, componentsCount)
        result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, componentsCount)
        result, d12 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, componentsCount)
        for c in range(componentsCount, 3):
            x.append(0.0)
            d1.append(0.0)
            d2.append(0.0)
            d12.append(0.0)
        cx.append(x)
        cd1.append(d1)
        cd2.append(d2)
        cd12.append(d12)
        node = nodeIter.next()
    return cx, cd1, cd2, cd12
