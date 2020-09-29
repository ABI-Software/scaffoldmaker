"""
Generates a 1-D path mesh.
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.interpolation import DerivativeScalingMode, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils import vector
from opencmiss.zinc.result import RESULT_OK


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
            'D2 derivatives': False,
            'D3 derivatives': False,
            'Length' : 1.0,
            'Number of elements' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Coordinate dimensions',
            'D2 derivatives',
            'D3 derivatives',
            'Length',
            'Number of elements'
        ]

    @staticmethod
    def checkOptions(options):
        dependentChanges = False
        if (options['Coordinate dimensions'] < 1) :
            options['Coordinate dimensions'] = 1
        elif (options['Coordinate dimensions'] > 3) :
            options['Coordinate dimensions'] = 3
        if (options['Number of elements'] < 1) :
            options['Number of elements'] = 1

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """
        coordinateDimensions = options['Coordinate dimensions']
        d2derivatives = options['D2 derivatives']
        d3derivatives = options['D3 derivatives']
        length = options['Length']
        elementsCount = options['Number of elements']

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule, components_count=coordinateDimensions)
        cache = fieldmodule.createFieldcache()

        #################
        # Create nodes
        #################

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        if d2derivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        if d3derivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)

        nodeIdentifier = 1
        x = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ length/elementsCount, 0.0, 0.0 ]
        if d2derivatives:
            dx_ds2 = [ 0.0, 1.0, 0.0 ]
            d2x_ds1ds2 = [ 0.0, 0.0, 0.0 ]
        if d3derivatives:
            dx_ds3 = [ 0.0, 0.0, 1.0 ]
            d2x_ds1ds3 = [ 0.0, 0.0, 0.0 ]
        for n in range(elementsCount + 1):
            x[0] = length*n/elementsCount
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            if d2derivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            if d3derivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, d2x_ds1ds3)
            nodeIdentifier = nodeIdentifier + 1

        #################
        # Create elements
        #################

        mesh = fieldmodule.findMeshByDimension(1)
        cubicHermiteBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(cubicHermiteBasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        for e in range(elementsCount):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1, e + 2 ])
            elementIdentifier = elementIdentifier + 1

        return []

    @classmethod
    def smoothPath(cls, region, options, mode : DerivativeScalingMode):
        x, d1 = extractPathParametersFromRegion(region, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1])
        d1 = smoothCubicHermiteDerivativesLine(x, d1, magnitudeScalingMode=mode)
        setPathParameters(region, [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1 ], [ x, d1 ])

    @classmethod
    def makeD2Normal(cls, region, options):
        if not options['D2 derivatives']:
            return
        d1, d2 = extractPathParametersFromRegion(region, [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2])
        for c in range(len(d1)):
            td2 = vector.vectorRejection(d2[c], d1[c])
            d2[c] = vector.setMagnitude(td2, vector.magnitude(d2[c]))
        setPathParameters(region, [Node.VALUE_LABEL_D_DS2], [d2])

    @classmethod
    def makeD3Normal(cls, region, options):
        if not options['D3 derivatives']:
            return
        if options['D2 derivatives']:
            d1, d2, d3 = extractPathParametersFromRegion(region, [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                                                                  Node.VALUE_LABEL_D_DS3])
            for c in range(len(d1)):
                d3[c] = vector.setMagnitude(vector.crossproduct3(d1[c], d2[c]), vector.magnitude(d3[c]))
            setPathParameters(region, [Node.VALUE_LABEL_D_DS3], [d3])
        else:
            d1, d3 = extractPathParametersFromRegion(region, [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS3])
            for c in range(len(d1)):
                td3 = vector.vectorRejection(d3[c], d1[c])
                d3[c] = vector.setMagnitude(td3, vector.magnitude(d3[c]))
            setPathParameters(region, [Node.VALUE_LABEL_D_DS3], [d3])


    @classmethod
    def getInteractiveFunctions(cls):
        """
        Supply client with functions for smoothing path parameters.
        """
        return [
            ("Smooth D1 arithmetic", lambda region, options: cls.smoothPath(region, options, DerivativeScalingMode.ARITHMETIC_MEAN)),
            ("Smooth D1 harmonic", lambda region, options: cls.smoothPath(region, options, DerivativeScalingMode.HARMONIC_MEAN)),
            ("Make D2 normal", lambda region, options: cls.makeD2Normal(region, options)),
            ("Make D3 normal", lambda region, options: cls.makeD3Normal(region, options))
        ]


def extractPathParametersFromRegion(region, valueLabels, groupName=None):
    '''
    Returns parameters of all nodes in region in identifier order.
    Assumes nodes in region have field coordinates (1 to 3 components).
    Currently limited to nodes with exactly value, d_ds1, d_ds2, d2_ds12,
    same as path 1 scaffold.
    :param valueLabels: List of parameters required as list of node value labels. e.g. [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1].
    :param groupName: Optional name of Zinc group to get parameters from.
    :return: cx, cd1, cd2, cd12 (all padded with zeroes to 3 components)
    '''
    fieldmodule = region.getFieldmodule()
    coordinates = fieldmodule.findFieldByName('coordinates').castFiniteElement()
    componentsCount = coordinates.getNumberOfComponents()
    assert componentsCount in [ 1, 2, 3 ], 'extractPathParametersFromRegion.  Invalid coordinates number of components'
    cache = fieldmodule.createFieldcache()

    valueLabelsCount = len(valueLabels)
    returnValues = [[] for i in range(valueLabelsCount)]
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    if groupName:
        group = fieldmodule.findFieldByName(groupName).castGroup()
        nodeGroup = group.getFieldNodeGroup(nodes)
        if nodeGroup.isValid():
            nodes = nodeGroup.getNodesetGroup()
        else:
            print('extractPathParametersFromRegion: missing group "' + groupName + '"')
    nodeIter = nodes.createNodeiterator()
    node = nodeIter.next()
    while node.isValid():
        cache.setNode(node)
        for i in range(valueLabelsCount):
            result, values = coordinates.getNodeParameters(cache, -1, valueLabels[i], 1, componentsCount)
            for c in range(componentsCount, 3):
                values.append(0.0)
            returnValues[i].append(values)
        node = nodeIter.next()

    return returnValues


def setPathParameters(region, nodeValueLabels, nodeValues):
    '''
    Set node parameters for coordinates field in path from listed values.
    :param nodeValueLabels: List of nodeValueLabels to set e.g. [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1 ]
    :param nodeValues: List of values for each type e.g. [ xlist, d1list ]
    '''
    fieldmodule = region.getFieldmodule()
    coordinates = fieldmodule.findFieldByName('coordinates').castFiniteElement()
    componentsCount = coordinates.getNumberOfComponents()
    # following requires at least one value label and node, assumes consistent values and components counts
    nodeValueLabelsCount = len(nodeValueLabels)
    nodesCount = len(nodeValues[0])
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    assert nodesCount == nodes.getSize()
    with ChangeManager(fieldmodule):
        cache = fieldmodule.createFieldcache()
        nodeIter = nodes.createNodeiterator()
        node = nodeIter.next()
        n = 0
        while node.isValid():
            cache.setNode(node)
            for v in range(nodeValueLabelsCount):
                coordinates.setNodeParameters(cache, -1, nodeValueLabels[v], 1, nodeValues[v][n])
            node = nodeIter.next()
            n += 1
