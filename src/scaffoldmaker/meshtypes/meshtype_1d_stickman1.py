"""
Generates a 1-D stickman that is used for body posture.
"""

from __future__ import division

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.interpolation import smoothCubicHermiteCrossDerivativesLine


class MeshType_1d_stickman1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '1D Stickman 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Coordinate dimensions': 3,
            'Arm length': 1.0,
            'Leg length': 1.0,
            'Left arm angle': 0.0,
            'Left leg angle': 1.2,
            'Right arm angle': 0.0,
            'Right leg angle': 1.2,
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Coordinate dimensions',
            'Arm length',
            'Leg length',
            'Left arm angle',
            'Left leg angle',
            'Right arm angle',
            'Right leg angle',
        ]

    @staticmethod
    def checkOptions(options):
        dependentChanges = False

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup
        """

        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule, components_count=3)
        cache = fieldmodule.createFieldcache()

        arm_length = options['Arm length']
        leg_length = options['Leg length']
        left_arm_angle = options['Left arm angle']
        left_leg_angle = options['Left leg angle']
        Right_arm_angle = options['Right arm angle']
        Right_leg_angle = options['Right leg angle']

        six_stick = True
        #################
        # Create nodes
        #################

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)

        nodeIdentifier = 1
        height = arm_length + arm_length + leg_length
        longitudinal = [0.0, 0.0, 1.0]
        transverse = [1.0, 0.0, 0.0]
        sagittal = vector.crossproduct3(longitudinal, transverse)
        x_head = vector.scaleVector(longitudinal, height)
        x_shoulder = vector.scaleVector(longitudinal, 0.75*height)
        x_hip = vector.scaleVector(longitudinal, 0.5*height)
        x_centre = vector.scaleVector(vector.addVectors([x_shoulder, x_hip], [1, 1]), 0.5)
        x_left_hand = vector.scaleVector(vector.rotateVectorAroundVector(transverse,
                                                                         sagittal, left_arm_angle), arm_length)
        x_left_hand = vector.addVectors([x_shoulder, x_left_hand], [1, 1])
        x_right_hand = vector.scaleVector(vector.rotateVectorAroundVector(vector.scaleVector(transverse, -1),
                                                                          sagittal, -Right_arm_angle), arm_length)
        x_right_hand = vector.addVectors([x_shoulder, x_right_hand], [1, 1])
        x_left_leg = vector.scaleVector(vector.rotateVectorAroundVector(transverse, sagittal,
                                                                        left_leg_angle), leg_length)
        x_left_leg = vector.addVectors([x_hip, x_left_leg], [1, 1])
        x_right_leg = vector.scaleVector(vector.rotateVectorAroundVector(vector.scaleVector(transverse, -1),
                                                                          sagittal, -Right_leg_angle), leg_length)
        x_right_leg = vector.addVectors([x_hip, x_right_leg], [1, 1])
        if six_stick:
            x_list = [x_head, x_shoulder, x_hip, x_right_hand, x_right_leg, x_left_hand, x_left_leg]
        else:
            x_list = [[0.0, 0.0, 1.0], [0.0, 0.0, 0.75], [0.0, 0.0, 0.5],
                      [-0.25, 0.0, 0.81], [-0.25, 0.0, 0.25], [0.25, 0.0, 0.81], [0.25, 0.0, 0.25], [-0.5, 0.0, 0.87] , [-0.5, 0.0, 0.0],
                      [0.5, 0.0, 0.87], [0.5, 0.0, 0.0]]

        x_list = [vector.addVectors([x, x_centre], [1, -1]) for x in x_list]
        for x in x_list:
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            nodeIdentifier += 1

        #################
        # Create elements
        #################

        mesh = fieldmodule.findMeshByDimension(1)
        linearLagrangeBasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        eft = mesh.createElementfieldtemplate(linearLagrangeBasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        if six_stick:
            lines = [[1, 2], [2, 3],
                     [2, 4], [3, 5],
                     [2, 6], [3, 7]]
        else:
            lines = [[1, 2], [2, 3],
                     [2, 4], [3, 5], [2, 6], [3, 7],
                     [4, 8], [5, 9], [6, 10], [7, 11]]
        for e in lines:
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, e)
            elementIdentifier += 1

        return []


def extractPathParametersFromRegion(region, valueLabels, groupName=None):
    '''
    Returns parameters of all nodes in region in identifier order.
    Assumes nodes in region have field coordinates (1 to 3 components).
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


def setPathParameters(region, nodeValueLabels, nodeValues, editGroupName=None):
    '''
    Set node parameters for coordinates field in path from listed values.
    :param nodeValueLabels: List of nodeValueLabels to set e.g. [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1 ]
    :param nodeValues: List of values for each type e.g. [ xlist, d1list ]
    :param editGroupName: Optional name of existing or new Zinc group to record modified nodes in.
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
        if editGroupName:
            editGroup = findOrCreateFieldGroup(fieldmodule, editGroupName, managed=True)
            editNodeGroup = editGroup.getFieldNodeGroup(nodes)
            if not editNodeGroup.isValid():
                editNodeGroup = editGroup.createFieldNodeGroup(nodes)
            editNodesetGroup = editNodeGroup.getNodesetGroup()
        cache = fieldmodule.createFieldcache()
        nodeIter = nodes.createNodeiterator()
        node = nodeIter.next()
        n = 0
        while node.isValid():
            cache.setNode(node)
            for v in range(nodeValueLabelsCount):
                coordinates.setNodeParameters(cache, -1, nodeValueLabels[v], 1, nodeValues[v][n])
            if editGroupName:
                editNodesetGroup.addNode(node)
            node = nodeIter.next()
            n += 1
