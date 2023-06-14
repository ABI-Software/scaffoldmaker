"""
Generates a 2-D diverging bifurcation network mesh from a 1-D network layout, with variable
numbers of elements around, and through wall.
"""

import copy
import math

from cmlibs.maths.vectorops import add, mult
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.bifurcation import get_tube_bifurcation_connection_elements_counts, \
    get_bifurcation_triple_point
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters


class MeshType_3d_divergingbifurcation1(Scaffold_base):
    """
    Generates a 2-D diverging bifurcation network mesh from a 1-D network layout, with variable
    numbers of elements around, and through wall.
    Magnitude of D2 and D3 are the radii of the tube in the respective directions.
    """
    parameterSetStructureStrings = {
        'Mouse 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5-6, 6.2-7-8-9-10, 6.3-11-12-13-14"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    # (1, [[13.84, 2.07, 20.50], [-0.53, -0.68, -6.63], [1.00, 0.07, -0.09], [-0.02, -0.05, -0.24], [0.08, -0.99, 0.10], [-0.02, -0.01, -0.01]]),
                    # (2, [[12.56, 1.40, 14.10], [-2.03, -0.65, -6.11], [0.95, 0.03, -0.32], [-0.08, -0.03, -0.22], [0.06, -1.00, 0.09], [-0.01, -0.00, -0.01]]),
                    # (3, [[9.85, 0.78, 8.39], [-3.33, -0.59, -5.33], [0.84, 0.01, -0.53], [-0.15, -0.01, -0.22], [0.06, -1.00, 0.07], [-0.00, -0.00, -0.03]]),
                    # (4, [[5.97, 0.24, 3.53], [-4.99, -0.41, -4.33], [0.66, 0.01, -0.76], [-0.24, 0.08, -0.20], [0.05, -1.00, 0.03], [0.01, 0.01, -0.11]]),
                    # (5, [[0.00, 0.00, 0.00], [[-6.82, -0.07, -2.67], [6.82, -0.07, -2.67], [0.00, 0.00, -2.00]],
                    #      [[0.66, 0.01, -0.76], [0.66, -0.01, 0.76], [1.00, 0.00, 0.00]], [-0.10, 0.09, 0.10],
                    #      [[0.05, -1.00, 0.03], [0.05, -1.00, 0.03], [0.05, -1.00, 0.03]], [-0.01, 0.01, -0.10]]),
                    # (6, [[-13.84, 2.07, 20.50], [0.53, -0.68, -6.63], [1.00, -0.04, 0.08], [0.10, -0.04, 0.40], [-0.05, -0.99, 0.10], [-0.03, -0.00, 0.04]]),
                    # (7, [[-12.56, 1.40, 14.10], [2.03, -0.65, -6.11], [0.96, -0.03, 0.32], [-0.08, 0.01, 0.22], [-0.06, -0.99, 0.09], [-0.01, -0.01, -0.01]]),
                    # (8, [[-9.85, 0.78, 8.39], [3.33, -0.59, -5.33], [0.85, -0.02, 0.53], [-0.15, 0.01, 0.22], [-0.07, -1.00, 0.07], [0.00, -0.01, -0.03]]),
                    # (9, [[-5.97, 0.24, 3.53], [4.99, -0.41, -4.33], [0.66, -0.01, 0.76], [0.04, 0.01, -0.20], [-0.05, -1.00, 0.03], [0.03, 0.00, -0.03]]),
                    # (10, [[0.00, 0.00, -2.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.07, 0.00, -0.15], [0.00, -1.00, 0.00], [0.01, -0.00, -0.01]]),
                    # (11, [[0.00, 0.00, -4.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    # (12, [[0.00, 0.00, -6.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    # (13, [[0.00, 0.00, -8.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    # (14, [[0.00, 0.00, -10.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]])]),

                    # (1, [[0.00, 0.00, -10.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    # (2, [[0.00, 0.00, -8.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    # (3, [[0.00, 0.00, -6.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    # (4, [[0.00, 0.00, -4.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    # (5, [[0.00, 0.00, -2.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.07, 0.00, -0.15], [0.00, -1.00, 0.00], [0.01, -0.00, -0.01]]),
                    # (6, [[0.00, 0.00, 0.00], [[0.00, 0.00, 2.00], [-6.82, -0.07, 2.67], [6.82, -0.07, -2.67]],
                    #      [[1.00, 0.00, 0.00], [0.66, -0.01, 0.76], [0.66, 0.01, -0.76]], [-0.10, 0.09, 0.10],
                    #      [[0.05, -1.00, 0.03], [0.05, -1.00, 0.03], [0.05, -1.00, 0.03]], [-0.01, 0.01, -0.10]]),
                    # (7, [[-5.97, 0.24, 3.53], [-4.99, -0.41, 4.33], [0.66, -0.01, 0.76], [0.04, 0.01, -0.20], [-0.05, -1.00, 0.03], [0.03, 0.00, -0.03]]),
                    # (8, [[-9.85, 0.78, 8.39], [-3.33, -0.59, 5.33], [0.85, -0.02, 0.53], [-0.15, 0.01, 0.22], [-0.07, -1.00, 0.07], [0.00, -0.01, -0.03]]),
                    # (9, [[-12.56, 1.40, 14.10], [-2.03, -0.65, 6.11], [0.96, -0.03, 0.32], [-0.08, 0.01, 0.22], [-0.06, -0.99, 0.09], [-0.01, -0.01, -0.01]]),
                    # (10, [[-13.84, 2.07, 20.50], [-0.53, -0.68, 6.63], [1.00, -0.04, 0.08], [0.10, -0.04, 0.40], [-0.05, -0.99, 0.10], [-0.03, -0.00, 0.04]]),
                    # (11, [[5.97, 0.24, 3.53], [4.99, -0.41, 4.33], [0.66, 0.01, -0.76], [-0.24, 0.08, -0.20], [0.05, -1.00, 0.03], [0.01, 0.01, -0.11]]),
                    # (12, [[9.85, 0.78, 8.39], [3.33, -0.59, 5.33], [0.84, 0.01, -0.53], [-0.15, -0.01, -0.22], [0.06, -1.00, 0.07], [-0.00, -0.00, -0.03]]),
                    # (13, [[12.56, 1.40, 14.10], [2.03, -0.65, 6.11], [0.95, 0.03, -0.32], [-0.08, -0.03, -0.22], [0.06, -1.00, 0.09], [-0.01, -0.00, -0.01]]),
                    # (14, [[13.84, 2.07, 20.50], [0.53, -0.68, 6.63], [1.00, 0.07, -0.09], [-0.02, -0.05, -0.24], [0.08, -0.99, 0.10], [-0.02, -0.01, -0.01]])]),

                    (1, [[0.00, 0.00, -10.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00],
                         [0.00, 1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (2, [[0.00, 0.00, -8.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00],
                         [0.00, 1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (3, [[0.00, 0.00, -6.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00],
                         [0.00, 1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (4, [[0.00, 0.00, -4.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00],
                         [0.00, 1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (5, [[0.00, 0.00, -2.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.07, 0.00, -0.15],
                         [0.00, 1.00, 0.00], [0.01, -0.00, -0.01]]),
                    (6, [[0.00, 0.00, 0.00], [[0.00, 0.00, 2.00], [-6.82, 0.07, 2.67], [6.82, 0.07, 2.67]],
                         [[1.00, 0.00, 0.00], [0.37, -0.01, 0.94], [0.37, 0.01, -0.94]], [-0.10, 0.09, 0.10],
                         [[0.00, 1.00, 0.00], [0.01, 1.00, 0.00], [-0.01, 1.00, 0.00]], [-0.01, 0.01, -0.10]]),
                    (7, [[-5.97, 0.24, 3.53], [-4.99, 0.41, 4.33], [0.66, -0.01, 0.76], [0.04, 0.01, -0.20],
                         [0.05, 1.00, -0.03], [0.03, 0.00, -0.03]]),
                    (8, [[-9.85, 0.78, 8.39], [-3.33, 0.59, 5.33], [0.85, -0.02, 0.53], [-0.15, 0.01, 0.22],
                         [0.07, 1.00, -0.07], [0.00, -0.01, -0.03]]),
                    (9, [[-12.56, 1.40, 14.10], [-2.03, 0.65, 6.11], [0.96, -0.03, 0.32], [-0.08, 0.01, 0.22],
                         [0.06, 0.99, -0.09], [-0.01, -0.01, -0.01]]),
                    (10, [[-13.84, 2.07, 20.50], [-0.53, 0.68, 6.63], [1.00, -0.04, 0.08], [0.10, -0.04, 0.40],
                          [0.05, 0.99, -0.10], [-0.03, -0.00, 0.04]]),
                    (11, [[5.97, 0.24, 3.53], [4.99, 0.41, 4.33], [0.66, 0.01, -0.76], [-0.24, 0.08, -0.20],
                          [-0.05, 1.00, -0.03], [0.01, 0.01, -0.11]]),
                    (12, [[9.85, 0.78, 8.39], [3.33, 0.59, 5.33], [0.84, 0.01, -0.53], [-0.15, -0.01, -0.22],
                          [-0.06, 1.00, -0.07], [-0.00, -0.00, -0.03]]),
                    (13, [[12.56, 1.40, 14.10], [2.03, 0.65, 6.11], [0.95, 0.03, -0.32], [-0.08, -0.03, -0.22],
                          [-0.06, 1.00, -0.09], [-0.01, -0.00, -0.01]]),
                    (14, [[13.84, 2.07, 20.50], [0.53, 0.68, 6.63], [1.00, 0.07, -0.09], [-0.02, -0.05, -0.24],
                          [-0.08, 0.99, -0.10], [-0.02, -0.01, -0.01]])]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-5',
                    'name': 'parent',
                    'ontId': 'None'
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-9',
                    'name': 'child 1',
                    'ontId': 'None'
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '10-13',
                    'name': 'child 2',
                    'ontId': 'None'
                }]
        })
    }

    @staticmethod
    def getName():
        return '3D Diverged Bifurcation 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        networkLayoutOption = cls.parameterSetStructureStrings['Mouse 1']
        options = {
            'Network layout': copy.deepcopy(networkLayoutOption),
            'Target element length': 4.0,
            'Number of elements around': 8,
            'Wall thickness': 1.0,
            'Number of elements through wall': 1,
            'Use linear through wall': True,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Network layout',
            'Target element length',
            'Number of elements around',
            'Wall thickness',
            'Number of elements through wall',
            'Use linear through wall',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall']
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Network layout':
            return [MeshType_1d_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Network layout':
            return list(cls.parameterSetStructureStrings.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        """
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        """
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Network layout':
            if not parameterSetName:
                parameterSetName = list(cls.parameterSetStructureStrings.keys())[1]
            return copy.deepcopy(cls.parameterSetStructureStrings[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Network layout'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Network layout'):
            options['Network layout'] = cls.getOptionScaffoldPackage('Network layout', MeshType_1d_network_layout1)
        for key in [
            'Number of elements around',
            'Number of elements through wall',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall'
        ]:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements around'] % 2 != 0:
            options['Number of elements around'] += 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        networkLayout = options['Network layout']
        elementsCountAround = options['Number of elements around']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        targetElementLength = options['Target element length']
        useCrossDerivatives = options['Use cross derivatives']

        fm = region.getFieldmodule()
        fm.beginChange()
        mesh = fm.findMeshByDimension(3)
        coordinates = findOrCreateFieldCoordinates(fm)

        firstNodeIdentifier = 1
        elementIdentifier = 1

        # Create annotation groups
        child2Group = AnnotationGroup(region, ("child 2", "None"))
        child1Group = AnnotationGroup(region, ("child 1", "None"))
        parentGroup = AnnotationGroup(region, ("parent", "None"))
        bifurcationGroup = AnnotationGroup(region, ("diverging bifurcation", "None"))
        annotationGroups = [parentGroup, child1Group, child2Group, bifurcationGroup]

        child2MeshGroup = child2Group.getMeshGroup(mesh)
        child1MeshGroup = child1Group.getMeshGroup(mesh)
        parentMeshGroup = parentGroup.getMeshGroup(mesh)
        bifurcationMeshGroup = bifurcationGroup.getMeshGroup(mesh)

        # Geometric coordinates
        geometricNetworkLayout = BifurcationNetworkLayout(region, networkLayout, targetElementLength)

        child2Length = geometricNetworkLayout.arcLengthOfGroupsAlong[2]
        child1Length = geometricNetworkLayout.arcLengthOfGroupsAlong[1]
        parentLength = geometricNetworkLayout.arcLengthOfGroupsAlong[0]

        elementsCountInChild2 = math.ceil(child2Length / targetElementLength)
        elementsCountInChild1 = math.ceil(child1Length / targetElementLength)
        elementsCountInParent = math.ceil(parentLength / targetElementLength)

        cx_child2_group = geometricNetworkLayout.cxGroups[2]
        cx_child1_group = geometricNetworkLayout.cxGroups[1]
        cx_parent_group = geometricNetworkLayout.cxGroups[0]

        sx_child2_group = geometricNetworkLayout.sxGroups[2]
        sx_child1_group = geometricNetworkLayout.sxGroups[1]
        sx_parent_group = geometricNetworkLayout.sxGroups[0]

        # Get parent nodes
        parentCoordinates = getCoordinatesAlongTube3D(cx_parent_group, elementsCountAround, elementsCountInParent,
                                                      elementsCountThroughWall, wallThickness, startRadian=-math.pi/2)

        parentLastRingNodeCoordinates = getTargetedRingNodesCoordinates(parentCoordinates, elementsCountAround,
                                                                    elementsCountInParent, elementsCountThroughWall,
                                                                    omitStartRows=0, omitEndRows=1)

        # Get child1 nodes
        child1Coordinates = getCoordinatesAlongTube3D(cx_child1_group, elementsCountAround,
                                                        elementsCountInChild1, elementsCountThroughWall,
                                                        wallThickness, startRadian=-math.pi/2)

        child1FirstRingNodeCoordinates = getTargetedRingNodesCoordinates(child1Coordinates, elementsCountAround,
                                                                    elementsCountInChild1, elementsCountThroughWall,
                                                                    omitStartRows=1, omitEndRows=0)

        # Get child2 nodes
        child2Coordinates = getCoordinatesAlongTube3D(cx_child2_group, elementsCountAround,
                                                         elementsCountInChild2, elementsCountThroughWall,
                                                         wallThickness, startRadian=-math.pi/2)

        child2FirstRingNodeCoordinates = getTargetedRingNodesCoordinates(child2Coordinates, elementsCountAround,
                                                                    elementsCountInChild2, elementsCountThroughWall,
                                                                    omitStartRows=1, omitEndRows=0)

        # Create nodes
        # Create parent nodes
        nodeIdentifier = generateTubeNodes(fm, firstNodeIdentifier, parentCoordinates, elementsCountInParent,
                                           elementsCountAround, elementsCountThroughWall, omitStartRows=0,
                                           omitEndRows=1, startNodes=None)

        # Create bifurcation nodes
        paCentre = sx_parent_group[0][1]
        c1Centre = sx_child2_group[0][-2]
        c2Centre = sx_child1_group[0][-2]
        paxList = parentLastRingNodeCoordinates[0]
        # pad1List = cFirstRingNodeCoordinates[1]
        # pad2List = cFirstRingNodeCoordinates[2]
        pad2 = parentLastRingNodeCoordinates[2]
        c1xList = child2FirstRingNodeCoordinates[0]
        c1d2 = child2FirstRingNodeCoordinates[2]
        c2xList = child1FirstRingNodeCoordinates[0]
        c2d2 = child1FirstRingNodeCoordinates[2]
        nodeIdentifier, roNodeId, coNodeId, nextNodeId = create3dBifurcationNodes(fm, nodeIdentifier, paCentre, paxList, pad2,
                                                                      c1Centre, c1xList, c1d2, c2Centre, c2xList, c2d2,
                                                                      elementsCountThroughWall)


        # Create child1 nodes
        nodeIdentifier = generateTubeNodes(fm, nodeIdentifier, child1Coordinates, elementsCountInChild1,
                                           elementsCountAround, elementsCountThroughWall, omitStartRows=1,
                                           omitEndRows=0, startNodes=None)

        # Create child2 nodes
        nodeIdentifier = generateTubeNodes(fm, nodeIdentifier,child2Coordinates, elementsCountInChild2,
                                           elementsCountAround, elementsCountThroughWall, omitStartRows=1,
                                           omitEndRows=0, startNodes=None)

        # Create elements
        # Create parent elements
        startNodeId = firstNodeIdentifier
        elementIdentifier = generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountInParent,
                                                 elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                                 omitStartRows=1, omitEndRows=0,
                                                 meshGroups=[parentMeshGroup, bifurcationMeshGroup])

        # Create bifurcation elements
        parentLastRingNodeId, nodeCount = getTargetedRingNodesId(firstNodeIdentifier, elementsCountAround,
                                                             elementsCountInParent,
                                                             elementsCountThroughWall, omitStartRows=0, omitEndRows=1)
        child1FirstRingNodeId, nodeCount = getTargetedRingNodesId(nextNodeId, elementsCountAround, elementsCountInChild1,
                                                             elementsCountThroughWall, omitStartRows=1, omitEndRows=0)
        child2FirstRingNodeId, nodeCount = getTargetedRingNodesId(nodeCount, elementsCountAround,
                                                             elementsCountInChild2, elementsCountThroughWall,
                                                             omitStartRows=1, omitEndRows=0)

        paNodeId = parentLastRingNodeId
        c1NodeId = child2FirstRingNodeId
        c2NodeId = child1FirstRingNodeId
        elementIdentifier = make_tube_bifurcation_elements_3d_diverging(region, coordinates, elementIdentifier,
                                                              elementsCountAround, elementsCountThroughWall, paNodeId,
                                                              c1NodeId, c2NodeId, roNodeId, coNodeId,
                                                              meshGroups=[parentMeshGroup, child2MeshGroup,
                                                                          child1MeshGroup, bifurcationMeshGroup])

        # Create child1 elements
        startNodeId = child1FirstRingNodeId[0][0]
        elementIdentifier = generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountInChild1,
                                                 elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                                 omitStartRows=1, omitEndRows=0,
                                                 meshGroups=[child1MeshGroup, bifurcationMeshGroup])

        # Create child2 elements
        startNodeId = child2FirstRingNodeId[0][0]
        elementIdentifier = generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountInChild2,
                                                 elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                                 omitStartRows=0, omitEndRows=1,
                                                 meshGroups=[child2MeshGroup, bifurcationMeshGroup])


        fm.endChange()
        return annotationGroups, None

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return


def findDerivativeBetweenPoints(v1, v2):
    """
    Find vector difference between two points and rescale vector difference using cubic hermite arclength
    between the points to derive the derivative between the points.
    :param v1: start vector
    :param v2: end vector
    :return: derivative of between v1 and v2
    """
    d = [v2[c] - v1[c] for c in range(3)]
    arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
    d = [c * arcLengthAround for c in vector.normalise(d)]

    return d


class BifurcationNetworkLayout:
    """
    Generates network layout for bifurcation scaffold.
    """
    def __init__(self, region, networkLayout, targetElementLength):
        """
        :param region: Zinc region to define model in.
        :param networkLayout: Network layout subscaffold from MeshType_1d_network_layout1
        """

        layoutRegion = region.createRegion()
        layoutFieldmodule = layoutRegion.getFieldmodule()
        layoutNodes = layoutFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # layoutMesh = layoutFieldmodule.findMeshByDimension(1)
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        # layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        layoutCoordinates = findOrCreateFieldCoordinates(layoutFieldmodule)
        layoutFieldcache = layoutFieldmodule.createFieldcache()

        networkMesh = networkLayout.getConstructionObject()

        networkSegments = networkMesh.getNetworkSegments()

        arcLengthOfGroupsAlong = []
        elementsCountAlongList = []
        cxGroups = []
        sxGroups = []
        for n1 in range(len(networkSegments)):
            networkSegment = networkSegments[n1]
            segmentNodes = networkSegment.getNetworkNodes()
            segmentNodeCount = len(segmentNodes)
            for n in range(segmentNodeCount):
                segmentNode = segmentNodes[n]
                layoutNodeIdentifier = segmentNode.getNodeIdentifier()
                layoutNode = layoutNodes.findNodeByIdentifier(layoutNodeIdentifier)
                layoutFieldcache.setNode(layoutNode)
                cx, cd1, cd2, cd3 = get_nodeset_path_ordered_field_parameters(
                    layoutNodes, layoutCoordinates, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                     Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3],
                    networkSegments[n1].getNodeIdentifiers(), networkSegments[n1].getNodeVersions())
            cxGroup = [cx, cd1, cd2, [], cd3, []]
            arcLength = 0.0
            for e in range(len(cx) - 1):
                arcLength += interp.getCubicHermiteArcLength(cx[e], cd1[e],
                                                             cx[e + 1], cd1[e + 1])
            arcLengthOfGroupsAlong.append(arcLength)
            elementsAlong = math.ceil(arcLength / targetElementLength)
            elementsCountAlongList.append(elementsAlong)
            cxGroups.append(cxGroup)
            # Sample
            sx, sd1, pe, pxi, psf = interp.sampleCubicHermiteCurves(cx, cd1, elementsAlong)
            # sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, pe, pxi, psf)
            # sd3, sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, pe, pxi, psf)
            sxGroup = [sx, sd1]
            sxGroups.append(sxGroup)

        del layoutCoordinates
        del layoutNodes
        del layoutFieldmodule
        del layoutRegion

        self.cxGroups = cxGroups
        self.sxGroups = sxGroups
        self.arcLengthOfGroupsAlong = arcLengthOfGroupsAlong
        self.elementsCountAlongList = elementsCountAlongList


def getCoordinatesAlongTube3D(cx_group, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                              wallThickness, startRadian):

    # Create ellipses along tube around the central path
    xEllipsesAlong = []
    d1EllipsesAlong = []
    for n in range(len(cx_group[0])):
        px, pd1 = createEllipsePoints(cx_group[0][n], 2 * math.pi, cx_group[2][n], cx_group[4][n], elementsCountAround,
                                      startRadians=startRadian)
        xEllipsesAlong.append(px)
        d1EllipsesAlong.append(pd1)

    # Find d2
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipsesAlong) - 1):
            v1 = xEllipsesAlong[n2][n1]
            v2 = xEllipsesAlong[n2 + 1][n1]
            d2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v1)
            d2Along.append(d2)
        xAlong.append(xEllipsesAlong[-1][n1])
        d2Along.append(d2)
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
        d2Raw.append(d2Smoothed)

    # Rearrange d2
    d2EllipsesAlong = []
    for n2 in range(len(xEllipsesAlong)):
        d2Around = []
        for n1 in range(elementsCountAround):
            d2 = d2Raw[n1][n2]
            d2Around.append(d2)
        d2EllipsesAlong.append(d2Around)

    # Spread out elements along tube
    xRaw = []
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipsesAlong)):
            xAlong.append(xEllipsesAlong[n2][n1])
            d2Along.append(d2EllipsesAlong[n2][n1])
        xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along, elementsCountAlongTube,
                                                                        arcLengthDerivatives=True)[0:2]
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlong, d2SampledAlong)
        xRaw.append(xSampledAlong)
        d2Raw.append(d2SampledAlong)

    # Rearrange x and d2
    xSampledTube = []
    d1SampledTube = []
    d2SampledTube = []
    for n2 in range(elementsCountAlongTube + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAround):
            x = xRaw[n1][n2]
            d2 = d2Raw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)
            # Calculate d1
            v1 = xRaw[n1][n2]
            v2 = xRaw[n1 + 1 if n1 < elementsCountAround - 1 else 0][n2]
            d1 = findDerivativeBetweenPoints(v1, v2)
            d1Around.append(d1)
            # if n2 > 0:
            #     v1 = xRaw[n1][n2]
            #     v2 = xRaw[n1 + 1 if n1 < elementsCountAround - 2 else 0][n2]
            #     d1 = findDerivativeBetweenPoints(v1, v2)
            #     d1Around.append(d1)
            # else:
            #     d1Around.append(d2Raw[n2][0])
        # if n2 > 0:
        d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        # else:
        #     d1Smoothed = d1Around
        xSampledTube.append(xAround)
        d1SampledTube.append(d1Smoothed)
        d2SampledTube.append(d2Around)

    d3Tube = []
    for n2 in range(elementsCountAlongTube + 1):
        d3Around = []
        for n1 in range(elementsCountAround):
            d3Around.append(vector.normalise(
                vector.crossproduct3(vector.normalise(d1SampledTube[n2][n1]), vector.normalise(d2SampledTube[n2][n1]))))
        d3Tube.append(d3Around)

    xInner = []
    d1Inner = []
    d2Inner = []
    d3Inner = []
    for n2 in range(elementsCountAlongTube + 1):
        for n1 in range(elementsCountAround):
            n = n2 * elementsCountAround + n1
            xInner.append(xSampledTube[n2][n1])
            d1Inner.append(d1SampledTube[n2][n1])
            d2Inner.append(d2SampledTube[n2][n1])
            d3Inner.append(d3Tube[n2][n1])

    transitElementList = [0] * elementsCountAround
    # if elementsCountThroughWall == 2:
    #     relativeThicknessList = [0.3, 0.7]
    # else:
    #     relativeThicknessList = []
    relativeThicknessList = []
    xList, d1List, d2List, d3List, curvatureList = \
        tubemesh.getCoordinatesFromInner(xInner, d1Inner, d2Inner, d3Inner, [wallThickness]*(elementsCountAlongTube+1),
                                         relativeThicknessList, elementsCountAround, elementsCountAlongTube,
                                         elementsCountThroughWall, transitElementList)

    coordinatesList = [xList, d1List, d2List, d3List]

    return coordinatesList


def generateTubeNodes(fm, nodeIdentifier, tubeCoordinates, elementsCountAlongTube, elementsCountAround,
                      elementsCountThroughWall, omitStartRows, omitEndRows, startNodes=None):

    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    # Create tube nodes
    lastRingsNodeId = []
    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    d3LastRing = []
    firstRingsNodeId = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    d3FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            xLastRingThroughWall = []
            d1LastRingThroughWall = []
            d2LastRingThroughWall = []
            d3LastRingThroughWall = []
            firstRingNodeIdThroughWall = []
            xFirstRingThroughWall = []
            d1FirstRingThroughWall = []
            d2FirstRingThroughWall = []
            d3FirstRingThroughWall = []
            for n1 in range(elementsCountAround):
                n = n2 * elementsCountAround * (elementsCountThroughWall + 1) + n3 * elementsCountAround + n1
                x = tubeCoordinates[0][n]
                d1 = tubeCoordinates[1][n]
                d2 = tubeCoordinates[2][n]
                d3 = tubeCoordinates[3][n]
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        if n2 == elementsCountAlongTube - 1:
                            lastRingNodeIdThroughWall.append(nodeIdentifier)
                            xLastRingThroughWall.append(x)
                            d1LastRingThroughWall.append(d1)
                            d2LastRingThroughWall.append(d2)
                            d3LastRingThroughWall.append(d3)
                        nodeIdentifier += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        if n2 == 1:
                            firstRingNodeIdThroughWall.append(nodeIdentifier)
                            xFirstRingThroughWall.append(x)
                            d1FirstRingThroughWall.append(d1)
                            d2FirstRingThroughWall.append(d2)
                            d3FirstRingThroughWall.append(d3)
                        nodeIdentifier += 1
                else:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                    nodeIdentifier += 1
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
                    xLastRing.append(xLastRingThroughWall)
                    d1LastRing.append(d1LastRingThroughWall)
                    d2LastRing.append(d2LastRingThroughWall)
                    d3LastRing.append(d3LastRingThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)
                    xFirstRing.append(xFirstRingThroughWall)
                    d1FirstRing.append(d1FirstRingThroughWall)
                    d2FirstRing.append(d2FirstRingThroughWall)
                    d3FirstRing.append(d3FirstRingThroughWall)

    return nodeIdentifier


def generateTubeElements(fm, startNodeId, elementIdentifier, elementsCountAlongTube, elementsCountAround,
                         elementsCountThroughWall, useCrossDerivatives, omitStartRows, omitEndRows, startNodes=None,
                         meshGroups=None):

    mesh = fm.findMeshByDimension(3)
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    # Create tube elements
    now = elementsCountAround * (elementsCountThroughWall + 1)
    for e3 in range(elementsCountThroughWall):
        for e2 in range(elementsCountAlongTube - 1 if omitStartRows or omitEndRows == 1 else elementsCountAlongTube):
            for e1 in range(elementsCountAround):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni11 = e2 * now + e3 * elementsCountAround + e1 + startNodeId
                bni12 = e2 * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + startNodeId
                bni21 = e2 * now + (e3 + 1) * elementsCountAround + e1 + startNodeId
                bni22 = e2 * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + startNodeId
                nodeIdentifiers = [bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now]
                # print('nodeIdentifiers', nodeIdentifiers)
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1

    # if omitStartRows == 1 or omitEndRows == 1:
    #     lastNodeId = elementsCountAlongTube * elementsCountAround * (elementsCountThroughWall + 1) + startNodeId
    # else:
    #     lastNodeId = (elementsCountAlongTube + 1) * elementsCountAround * (elementsCountThroughWall + 1) + startNodeId
    return elementIdentifier


def make_tube_bifurcation_elements_3d(region, coordinates, elementIdentifier, elementsCountAround,
                                      elementsCountThroughWall, paNodeId, c1NodeId, c2NodeId, roNodeId, coNodeId,
                                      meshGroups=None):

    paCount = len(paNodeId[0])
    c1Count = len(c1NodeId[0])
    c2Count = len(c2NodeId[0])
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)

    fm = region.getFieldmodule()
    mesh = fm.findMeshByDimension(3)
    eftfactory = eftfactory_bicubichermitelinear(mesh, None)
    eftStd = eftfactory.createEftBasic()

    elementtemplateStd = mesh.createElementtemplate()
    elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplateStd.defineField(coordinates, -1, eftStd)

    elementtemplateMod = mesh.createElementtemplate()
    elementtemplateMod.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    # Child2 part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(c1Count):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = roNodeId[e3][(e1 + 1) % c1Count]
                bni5 = bni1 + elementsCountAround
                bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            elif e1 == elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = coNodeId[e3][e1 - elementsCountAround // 2]
                bni5 = bni1 + elementsCountAround
                bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = coNodeId[e3][e1 - elementsCountAround // 2 - 1]
                bni5 = bni1 + elementsCountAround
                if e1 == c1Count - 1:
                    bni4 = roNodeId[e3][0]
                    bni6 = bni5 - elementsCountAround + 1
                else:
                    bni4 = bni3 + 1
                    bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 in (0, pac1Count - 1, pac1Count, c1Count - 1):
                eft = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                elif e1 == pac1Count - 1:
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == pac1Count:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif e1 == c1Count - 1:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 1:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # Child1 part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(c2Count):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround//2:
                bni1 = c2NodeId[e3][e1]
                bni2 = c2NodeId[e3][(e1 + 1) % c1Count]
                if e1 == 0:
                    bni3 = roNodeId[e3][e1]
                    bni4 = coNodeId[e3][-1 - e1]
                else:
                    bni3 = coNodeId[e3][-e1]
                    if e1 == elementsCountAround//2 - 1:
                        bni4 = roNodeId[e3][0] + pac1Count
                    else:
                        bni4 = bni3 - 1
                bni5 = bni1 + elementsCountAround
                bni6 = bni2 + elementsCountAround
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni1 = c2NodeId[e3][e1]
                bni2 = c2NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = roNodeId[e3][(e1 + 1) % c2Count]
                bni5 = bni1 + elementsCountAround
                bni6 = bni2 + elementsCountAround
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 <= c1c2Count:
                eft = eftfactory.createEftBasic()
                scalefactors = [-1.0]
                setEftScaleFactorIds(eft, [1], [])
                if e1 == 0:
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    scaleEftNodeValueLabels(eft, [4, 8], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                elif e1 < (c1c2Count - 1):
                    scaleEftNodeValueLabels(eft, [3, 4, 7, 8], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                elif e1 == (c1c2Count - 1):
                    scaleEftNodeValueLabels(eft, [3, 7], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == c1c2Count:
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            elif e1 == c2Count - 1:
                eft = eftfactory.createEftBasic()
                remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2,
                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 2:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # parent part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(paCount):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = roNodeId[e3][e1]
            bni2 = roNodeId[e3][(e1 + 1) % elementsCountAround]
            bni3 = paNodeId[e3][e1]
            bni4 = paNodeId[e3][(e1 + 1) % elementsCountAround]
            bni5 = bni1 + len(roNodeId[0]) + len(coNodeId[0])
            bni6 = bni2 + len(roNodeId[0]) + len(coNodeId[0])
            bni7 = bni3 + elementsCountAround
            bni8 = bni4 + elementsCountAround
            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            # print('nodeIdentifiers', nodeIdentifiers)
            if e1 in (0, pac1Count):
                eft = eftfactory.createEftBasic()
                if e1 in (0, pac1Count):
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 0:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)
    # lastNodeId = startNodeIndex + (elementsCountThroughWall + 1) * (len(roNodeId) + len(coNodeId)) + 1
    # lastNodeId = paNodeId[0][0]
    return elementIdentifier


def getTargetedRingNodesCoordinates(tubeCoordinates, elementsCountAround, elementsCountAlongTube,
                                    elementsCountThroughWall, omitStartRows, omitEndRows):

    # Create tube nodes
    lastRingsNodeId = []
    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    d3LastRing = []
    firstRingsNodeId = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    d3FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        # if n2 == 0:
        #     startNodeId = nodeIdentifier - 1
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            xLastRingThroughWall = []
            d1LastRingThroughWall = []
            d2LastRingThroughWall = []
            d3LastRingThroughWall = []
            firstRingNodeIdThroughWall = []
            xFirstRingThroughWall = []
            d1FirstRingThroughWall = []
            d2FirstRingThroughWall = []
            d3FirstRingThroughWall = []
            for n1 in range(elementsCountAround):
                n = n2 * elementsCountAround * (elementsCountThroughWall + 1) + n3 * elementsCountAround + n1
                x = tubeCoordinates[0][n]
                d1 = tubeCoordinates[1][n]
                d2 = tubeCoordinates[2][n]
                d3 = tubeCoordinates[3][n]
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        if n2 == elementsCountAlongTube - 1:
                            # lastRingNodeIdThroughWall.append(nodeCount)
                            xLastRingThroughWall.append(x)
                            d1LastRingThroughWall.append(d1)
                            d2LastRingThroughWall.append(d2)
                            d3LastRingThroughWall.append(d3)
                        # nodeCount += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        if n2 == 1:
                            # firstRingNodeIdThroughWall.append(nodeCount)
                            xFirstRingThroughWall.append(x)
                            d1FirstRingThroughWall.append(d1)
                            d2FirstRingThroughWall.append(d2)
                            d3FirstRingThroughWall.append(d3)
                        # nodeCount += 1
                # else:
                #     nodeCount += 1
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
                    xLastRing.append(xLastRingThroughWall)
                    d1LastRing.append(d1LastRingThroughWall)
                    d2LastRing.append(d2LastRingThroughWall)
                    d3LastRing.append(d3LastRingThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)
                    xFirstRing.append(xFirstRingThroughWall)
                    d1FirstRing.append(d1FirstRingThroughWall)
                    d2FirstRing.append(d2FirstRingThroughWall)
                    d3FirstRing.append(d3FirstRingThroughWall)

    if omitStartRows == 1:
        targetedRingCoordinates = [xFirstRing, d1FirstRing, d2FirstRing, d3FirstRing]
    elif omitEndRows == 1:
        targetedRingCoordinates = [xLastRing, d1LastRing, d2LastRing, d3LastRing]
    else:
        targetedRingCoordinates = []

    return targetedRingCoordinates


def getTargetedRingNodesId(nodeCount, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                           omitStartRows, omitEndRows):

    # Create tube nodes
    lastRingsNodeId = []
    firstRingsNodeId = []
    for n2 in range(elementsCountAlongTube + 1):
        # if n2 == 0:
        #     startNodeId = nodeIdentifier - 1
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            firstRingNodeIdThroughWall = []
            for n1 in range(elementsCountAround):
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        if n2 == elementsCountAlongTube - 1:
                            lastRingNodeIdThroughWall.append(nodeCount)
                        nodeCount += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        if n2 == 1:
                            firstRingNodeIdThroughWall.append(nodeCount)
                        nodeCount += 1
                else:
                    nodeCount += 1
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)

    if omitStartRows == 1:
        targetedRingNodeId = firstRingsNodeId
    elif omitEndRows == 1:
        targetedRingNodeId = lastRingsNodeId
    else:
        targetedRingNodeId = []

    return targetedRingNodeId, nodeCount


def create3dBifurcationNodes(fm, nodeIdentifier, paCentre, paxList, pad2, c1Centre, c1xList, c1d2, c2Centre, c2xList,
                             c2d2, elementsCountThroughWall):

    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    roxList = []
    rod1List = []
    rod2List = []
    coxList = []
    cod1List = []
    cod2List = []
    for n3 in range(elementsCountThroughWall + 1):
        roxOuter, rod1Outer, rod2Outer, coxOuter, cod1Outer, cod2Outer, paStartIndex, c1StartIndex, c2StartIndex = \
            make_tube_bifurcation_points_diverging(paCentre, paxList[n3], pad2[n3], c1Centre, c1xList[n3], c1d2[n3],
                                                    c2Centre, c2xList[n3], c2d2[n3])
        roxList.append(roxOuter)
        rod1List.append(rod1Outer)
        rod2List.append(rod2Outer)
        coxList.append(coxOuter)
        cod1List.append(cod1Outer)
        cod2List.append(cod2Outer)

    # Create bifurcation nodes
    roNodeId = []
    coNodeId = []
    for n3 in range(elementsCountThroughWall + 1):
        coNodeIdThroughWall = []
        for n in range(len(coxOuter)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, coxList[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1List[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2List[n3][n])
            coNodeIdThroughWall.append(nodeIdentifier)
            nodeIdentifier += 1
        coNodeId.append(coNodeIdThroughWall)
        roNodeIdThroughWall = []
        for n in range(len(roxOuter)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, roxList[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1List[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2List[n3][n])
            roNodeIdThroughWall.append(nodeIdentifier)
            nodeIdentifier += 1
        roNodeId.append(roNodeIdThroughWall)
    nextNodeId = nodeIdentifier
    return nodeIdentifier, roNodeId, coNodeId, nextNodeId


def make_tube_bifurcation_points_converging(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):
    """
    Gets first ring of coordinates and derivatives between parent pa and
    children c1, c2, and over the crotch between c1 and c2.
    :return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex
    """
    paCount = len(pax)
    c1Count = len(c1x)
    c2Count = len(c2x)
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)
    # convert to number of nodes, includes both 6-way points
    pac1NodeCount = pac1Count + 1
    pac2NodeCount = pac2Count + 1
    c1c2NodeCount = c1c2Count + 1
    paStartIndex = 0
    c1StartIndex = 0
    c2StartIndex = 0
    pac1x = [None] * pac1NodeCount
    pac1d1 = [None] * pac1NodeCount
    pac1d2 = [None] * pac1NodeCount
    for n in range(pac1NodeCount):
        pan = (paStartIndex + n) % paCount
        c1n = (c1StartIndex + n) % c1Count
        x1, d1, x2, d2 = c1x[c1n], mult(c1d2[c1n], 2.0), pax[pan], mult(pad2[pan], 2.0)
        pac1x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac1d1[n] = [0.0, 0.0, 0.0]
        pac1d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    paStartIndex2 = paStartIndex + pac1Count
    c1StartIndex2 = c1StartIndex + pac1Count
    c2StartIndex2 = c2StartIndex + c1c2Count
    pac2x = [None] * pac2NodeCount
    pac2d1 = [None] * pac2NodeCount
    pac2d2 = [None] * pac2NodeCount
    for n in range(pac2NodeCount):
        pan = (paStartIndex2 + n) % paCount
        c2n = (c2StartIndex2 + n) % c2Count
        x1, d1, x2, d2 = c2x[c2n], mult(c2d2[c2n], 2.0), pax[pan], mult(pad2[pan], 2.0)
        pac2x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac2d1[n] = [0.0, 0.0, 0.0]
        pac2d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    c1c2x = [None] * c1c2NodeCount
    c1c2d1 = [None] * c1c2NodeCount
    c1c2d2 = [None] * c1c2NodeCount
    for n in range(c1c2NodeCount):
        c1n = (c1StartIndex2 + n) % c1Count
        c2n = (c2StartIndex2 - n) % c2Count  # note: reversed
        x1, d1, x2, d2 = c1x[c1n], mult(c1d2[c1n], 2.0), c2x[c2n], mult(c2d2[c2n], -2.0)
        c1c2x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        c1c2d1[n] = [0.0, 0.0, 0.0]
        c1c2d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    # get hex triple points
    hex1, hex1d1, hex1d2 = get_bifurcation_triple_point(
        c2x[c1StartIndex], mult(c2d2[c2StartIndex], -1.0),
        c1x[c1StartIndex], mult(c1d2[c1StartIndex], -1.0),
        pax[paStartIndex], pad2[paStartIndex])
    hex2, hex2d1, hex2d2 = get_bifurcation_triple_point(
        c1x[c1StartIndex2], mult(c1d2[c1StartIndex2], -1.0),
        c2x[c2StartIndex2], mult(c2d2[c2StartIndex2], -1.0),
        pax[paStartIndex2], pad2[paStartIndex2])
    # smooth around loops through hex points to get d1
    loop1x = [hex2] + pac2x[1:-1] + [hex1]
    loop1d1 = [[-d for d in hex2d2]] + pac2d1[1:-1] + [hex1d1]
    loop2x = [hex1] + pac1x[1:-1] + [hex2]
    loop2d1 = [[-d for d in hex1d2]] + pac1d1[1:-1] + [hex2d1]
    loop1d1 = interp.smoothCubicHermiteDerivativesLine(loop1x, loop1d1, fixStartDirection=True, fixEndDirection=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    loop2d1 = interp.smoothCubicHermiteDerivativesLine(loop2x, loop2d1, fixStartDirection=True, fixEndDirection=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    # smooth over "crotch" between c1 and c2
    crotchx = [hex2] + c1c2x[1:-1] + [hex1]
    crotchd1 = [add(hex2d1, hex2d2)] + c1c2d1[1:-1] + [[(-hex1d1[c] - hex1d2[c]) for c in range(3)]]
    crotchd1 = interp.smoothCubicHermiteDerivativesLine(crotchx, crotchd1, fixStartDerivative=True,
                                                        fixEndDerivative=True,
                                                        magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    rox = [hex1] + pac1x[1:-1] + [hex2] + pac2x[1:-1]
    rod1 = [hex1d1] + loop2d1[1:-1] + [hex2d1] + loop1d1[1:-1]
    rod2 = [hex1d2] + pac1d2[1:-1] + [hex2d2] + pac2d2[1:-1]
    cox = crotchx[1:-1]
    cod1 = crotchd1[1:-1]
    cod2 = c1c2d2[1:-1]
    return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex


def make_tube_bifurcation_points_diverging(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):
    """
    Gets first ring of coordinates and derivatives between parent pa and
    children c1, c2, and over the crotch between c1 and c2.
    :return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex
    """
    paCount = len(pax)
    c1Count = len(c1x)
    c2Count = len(c2x)
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)
    # convert to number of nodes, includes both 6-way points
    pac1NodeCount = pac1Count + 1
    pac2NodeCount = pac2Count + 1
    c1c2NodeCount = c1c2Count + 1
    paStartIndex = 0
    c1StartIndex = 0
    c2StartIndex = 0
    pac1x = [None] * pac1NodeCount
    pac1d1 = [None] * pac1NodeCount
    pac1d2 = [None] * pac1NodeCount
    for n in range(pac1NodeCount):
        pan = (paStartIndex + n) % paCount
        c1n = (c1StartIndex + n) % c1Count
        # x1, d1, x2, d2 = c1x[c1n], mult(c1d2[c1n], 2.0), pax[pan], mult(pad2[pan], 2.0)
        x1, d1, x2, d2 = pax[pan], mult(pad2[pan], 2.0), c1x[c1n], mult(c1d2[c1n], 2.0)
        pac1x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac1d1[n] = [0.0, 0.0, 0.0]
        pac1d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    paStartIndex2 = paStartIndex + pac1Count
    c1StartIndex2 = c1StartIndex + pac1Count
    c2StartIndex2 = c2StartIndex + c1c2Count
    pac2x = [None] * pac2NodeCount
    pac2d1 = [None] * pac2NodeCount
    pac2d2 = [None] * pac2NodeCount
    for n in range(pac2NodeCount):
        pan = (paStartIndex2 + n) % paCount
        c2n = (c2StartIndex2 + n) % c2Count
        # x1, d1, x2, d2 = c2x[c2n], mult(c2d2[c2n], 2.0), pax[pan], mult(pad2[pan], 2.0)
        x1, d1, x2, d2 = pax[pan], mult(pad2[pan], 2.0), c2x[c2n], mult(c2d2[c2n], 2.0)
        pac2x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac2d1[n] = [0.0, 0.0, 0.0]
        pac2d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    c1c2x = [None] * c1c2NodeCount
    c1c2d1 = [None] * c1c2NodeCount
    c1c2d2 = [None] * c1c2NodeCount
    for n in range(c1c2NodeCount):
        c1n = (c1StartIndex2 + n) % c1Count
        c2n = (c2StartIndex2 - n) % c2Count  # note: reversed
        # x1, d1, x2, d2 = c1x[c1n], mult(c1d2[c1n], 2.0), c2x[c2n], mult(c2d2[c2n], -2.0)
        x1, d1, x2, d2 = c2x[c2n], mult(c2d2[c2n], -2.0), c1x[c1n], mult(c1d2[c1n], 2.0)
        c1c2x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        c1c2d1[n] = [0.0, 0.0, 0.0]
        c1c2d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    # # get hex triple points
    # hex1, hex1d1, hex1d2 = get_bifurcation_triple_point(
    #     c2x[c1StartIndex], mult(c2d2[c2StartIndex], -1.0),
    #     c1x[c1StartIndex], mult(c1d2[c1StartIndex], -1.0),
    #     pax[paStartIndex], pad2[paStartIndex])
    # hex2, hex2d1, hex2d2 = get_bifurcation_triple_point(
    #     c1x[c1StartIndex2], mult(c1d2[c1StartIndex2], -1.0),
    #     c2x[c2StartIndex2], mult(c2d2[c2StartIndex2], -1.0),
    #     pax[paStartIndex2], pad2[paStartIndex2])
    # # smooth around loops through hex points to get d1
    # loop1x = [hex2] + pac2x[1:-1] + [hex1]
    # loop1d1 = [[-d for d in hex2d2]] + pac2d1[1:-1] + [hex1d1]
    # loop2x = [hex1] + pac1x[1:-1] + [hex2]
    # loop2d1 = [[-d for d in hex1d2]] + pac1d1[1:-1] + [hex2d1]
    # loop1d1 = interp.smoothCubicHermiteDerivativesLine(loop1x, loop1d1, fixStartDirection=True, fixEndDirection=True,
    #                                                    magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    # loop2d1 = interp.smoothCubicHermiteDerivativesLine(loop2x, loop2d1, fixStartDirection=True, fixEndDirection=True,
    #                                                    magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    # # smooth over "crotch" between c1 and c2
    # crotchx = [hex2] + c1c2x[1:-1] + [hex1]
    # crotchd1 = [add(hex2d1, hex2d2)] + c1c2d1[1:-1] + [[(-hex1d1[c] - hex1d2[c]) for c in range(3)]]
    # crotchd1 = interp.smoothCubicHermiteDerivativesLine(crotchx, crotchd1, fixStartDerivative=True,
    #                                                     fixEndDerivative=True,
    #                                                     magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    # rox = [hex1] + pac1x[1:-1] + [hex2] + pac2x[1:-1]
    # rod1 = [hex1d1] + loop2d1[1:-1] + [hex2d1] + loop1d1[1:-1]
    # rod2 = [hex1d2] + pac1d2[1:-1] + [hex2d2] + pac2d2[1:-1]
    # cox = crotchx[1:-1]
    # cod1 = crotchd1[1:-1]
    # cod2 = c1c2d2[1:-1]
    # return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex

    # get hex triple points
    hex1, hex1d1, hex1d2 = get_bifurcation_triple_point(
        pax[paStartIndex], mult(pad2[paStartIndex], -1.0),
        c1x[c1StartIndex], c1d2[c1StartIndex],
        c2x[c1StartIndex], c2d2[c2StartIndex])
    hex2, hex2d1, hex2d2 = get_bifurcation_triple_point(
        pax[paStartIndex2], mult(pad2[paStartIndex2], -1.0),
        c2x[c2StartIndex2], c2d2[c2StartIndex2],
        c1x[c1StartIndex2], c1d2[c1StartIndex2])
    # smooth around loops through hex points to get d1
    loop1x  = [ hex2 ] + pac2x[1:-1] + [ hex1 ]
    loop1d1 = [ [ -d for d in hex2d2 ] ] + pac2d1[1:-1] + [ hex1d1 ]
    loop2x  = [ hex1 ] + pac1x[1:-1] + [ hex2 ]
    loop2d1 = [ [ -d for d in hex1d2 ] ] + pac1d1[1:-1] + [ hex2d1 ]
    loop1d1 = interp.smoothCubicHermiteDerivativesLine(loop1x, loop1d1, fixStartDirection=True, fixEndDirection=True, magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    loop2d1 = interp.smoothCubicHermiteDerivativesLine(loop2x, loop2d1, fixStartDirection=True, fixEndDirection=True, magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    # smooth over "crotch" between c1 and c2
    crotchx = [ hex2 ] + c1c2x[1:-1] + [ hex1 ]
    crotchd1 = [ add(hex2d1, hex2d2) ] + c1c2d1[1:-1] + [ [ (-hex1d1[c] - hex1d2[c]) for c in range(3) ] ]
    crotchd1 = interp.smoothCubicHermiteDerivativesLine(crotchx, crotchd1, fixStartDerivative=True, fixEndDerivative=True, magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    rox  = [ hex1 ] + pac1x[1:-1] + [ hex2 ] + pac2x[1:-1]
    rod1 = [ loop1d1[-1] ] + loop2d1[1:] + loop1d1[1:-1]
    rod2 = [ [ -d for d in loop2d1[ 0] ] ] + pac1d2[1:-1] + [ [ -d for d in loop1d1[0] ] ] + pac2d2[1:-1]
    cox  = crotchx [1:-1]
    cod1 = crotchd1[1:-1]
    cod2 = c1c2d2[1:-1]
    return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex


def make_tube_bifurcation_elements_3d_diverging(region, coordinates, elementIdentifier, elementsCountAround,
                                      elementsCountThroughWall, paNodeId, c1NodeId, c2NodeId, roNodeId, coNodeId,
                                      meshGroups=None):

    paCount = len(paNodeId[0])
    c1Count = len(c1NodeId[0])
    c2Count = len(c2NodeId[0])
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)

    fm = region.getFieldmodule()
    mesh = fm.findMeshByDimension(3)
    eftfactory = eftfactory_bicubichermitelinear(mesh, None)
    eftStd = eftfactory.createEftBasic()

    elementtemplateStd = mesh.createElementtemplate()
    elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplateStd.defineField(coordinates, -1, eftStd)

    elementtemplateMod = mesh.createElementtemplate()
    elementtemplateMod.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    # parent part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(paCount):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = paNodeId[e3][e1]
            bni2 = paNodeId[e3][(e1 + 1) % elementsCountAround]
            bni3 = roNodeId[e3][e1]
            bni4 = roNodeId[e3][(e1 + 1) % elementsCountAround]
            bni5 = bni1 + elementsCountAround
            bni6 = bni2 + elementsCountAround
            bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
            bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            # print('nodeIdentifiers', nodeIdentifiers)
            if e1 in (0, pac1Count - 1, pac1Count, paCount - 1):
                eft = eftfactory.createEftBasic()
                if e1 in (0, pac1Count):
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elif e1 in (pac1Count - 1, paCount - 1):
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 0:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # Child1 part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(c2Count):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround // 2:
                bni3 = c2NodeId[e3][e1]
                bni4 = c2NodeId[e3][(e1 + 1) % c1Count]
                if e1 == 0:
                    bni1 = roNodeId[e3][e1]
                    bni2 = coNodeId[e3][-1 - e1]
                else:
                    bni1 = coNodeId[e3][-e1]
                    if e1 == elementsCountAround // 2 - 1:
                        bni2 = roNodeId[e3][0] + pac1Count
                    else:
                        bni2 = bni1 - 1
                bni5 = bni1 + len(roNodeId[0]) + len(coNodeId[0])
                bni6 = bni2 + len(roNodeId[0]) + len(coNodeId[0])
                bni7 = bni3 + elementsCountAround
                bni8 = bni4 + elementsCountAround
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            else:
                bni3 = c2NodeId[e3][e1]
                bni4 = c2NodeId[e3][(e1 + 1) % c1Count]
                bni1 = roNodeId[e3][e1]
                bni2 = roNodeId[e3][(e1 + 1) % c2Count]
                bni5 = bni1 + len(roNodeId[0]) + len(coNodeId[0])
                bni6 = bni2 + len(roNodeId[0]) + len(coNodeId[0])
                bni7 = bni3 + elementsCountAround
                bni8 = bni4 + elementsCountAround
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 <= c1c2Count:
                eft = eftfactory.createEftBasic()
                scalefactors = [-1.0]
                setEftScaleFactorIds(eft, [1], [])
                if e1 == 0:
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    scaleEftNodeValueLabels(eft, [2, 6], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                elif e1 < (c1c2Count - 1):
                    scaleEftNodeValueLabels(eft, [1, 2, 5, 6], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2],
                                            [1])
                elif e1 == (c1c2Count - 1):
                    scaleEftNodeValueLabels(eft, [1, 5], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                    remapEftNodeValueLabel(eft, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                elif e1 == c1c2Count:
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 2:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # Child2 part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(c1Count):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround // 2:
                bni3 = c1NodeId[e3][e1]
                bni4 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni1 = roNodeId[e3][e1]
                bni2 = roNodeId[e3][(e1 + 1) % c1Count]
                bni5 = bni1 + len(roNodeId[0]) + len(coNodeId[0])
                bni6 = bni5 + 1
                bni7 = bni3 + elementsCountAround
                bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            elif e1 == elementsCountAround // 2:
                bni3 = c1NodeId[e3][e1]
                bni4 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni1 = roNodeId[e3][e1]
                bni2 = coNodeId[e3][e1 - elementsCountAround // 2]
                bni5 = bni1 + len(roNodeId[0]) + len(coNodeId[0])
                bni6 = bni2 + len(roNodeId[0]) + len(coNodeId[0])
                bni7 = bni3 + elementsCountAround
                bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni3 = c1NodeId[e3][e1]
                bni4 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni1 = coNodeId[e3][e1 - elementsCountAround // 2 - 1]
                bni2 = bni1 + 1
                bni5 = bni1 + len(roNodeId[0]) + len(coNodeId[0])
                bni6 = bni2 + len(roNodeId[0]) + len(coNodeId[0])
                bni7 = bni3 + elementsCountAround
                bni8 = bni4 + elementsCountAround
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 in (0, pac1Count, c1Count - 1):
                eft = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                elif e1 == pac1Count:
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == c1Count - 1:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 1:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # lastNodeId = startNodeIndex + (elementsCountThroughWall + 1) * (len(roNodeId) + len(coNodeId)) + 1
    # lastNodeId = paNodeId[0][0]
    return elementIdentifier