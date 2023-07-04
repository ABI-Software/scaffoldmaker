"""
Generates a 2-D diverging bifurcation network mesh from a 1-D network layout.
"""

import copy
import math

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.bifurcation import make_tube_bifurcation_points, make_tube_bifurcation_elements_2d
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters


class MeshType_2d_divergingbifurcation1(Scaffold_base):
    """
    Generates a 2-D diverging bifurcation network mesh from a 1-D network layout.
    Magnitude of D2 and D3 are the radii of the tube in the respective directions.
    """
    parameterSetStructureStrings = {
        'Default': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5-6, 6.2-7-8-9-10, 6.3-11-12-13-14"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[0.00, 0.00, -10.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, 1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (2, [[0.00, 0.00, -8.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, 1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (3, [[0.00, 0.00, -6.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, 1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (4, [[0.00, 0.00, -4.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, 1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (5, [[0.00, 0.00, -2.00], [0.00, 0.00, 2.00], [1.00, 0.00, 0.00], [0.07, 0.00, -0.15], [0.00, 1.00, 0.00], [0.01, -0.00, -0.01]]),
                    (6, [[0.00, 0.00, 0.00], [[0.00, 0.00, 2.00], [-6.82, 0.07, 2.67], [6.82, 0.07, 2.67]],
                         [[1.00, 0.00, 0.00], [0.37, -0.01, 0.94], [0.37, 0.01, -0.94]], [-0.10, 0.09, 0.10],
                         [[0.00, 1.00, 0.00], [0.01, 1.00, 0.00], [-0.01, 1.00, 0.00]], [-0.01, 0.01, -0.10]]),
                    (7, [[-5.97, 0.24, 3.53], [-4.99, 0.41, 4.33], [0.66, -0.01, 0.76], [0.04, 0.01, -0.20], [0.05, 1.00, -0.03], [0.03, 0.00, -0.03]]),
                    (8, [[-9.85, 0.78, 8.39], [-3.33, 0.59, 5.33], [0.85, -0.02, 0.53], [-0.15, 0.01, 0.22], [0.07, 1.00, -0.07], [0.00, -0.01, -0.03]]),
                    (9, [[-12.56, 1.40, 14.10], [-2.03, 0.65, 6.11], [0.96, -0.03, 0.32], [-0.08, 0.01, 0.22], [0.06, 0.99, -0.09], [-0.01, -0.01, -0.01]]),
                    (10, [[-13.84, 2.07, 20.50], [-0.53, 0.68, 6.63], [1.00, -0.04, 0.08], [0.10, -0.04, 0.40], [0.05, 0.99, -0.10], [-0.03, -0.00, 0.04]]),
                    (11, [[5.97, 0.24, 3.53], [4.99, 0.41, 4.33], [0.66, 0.01, -0.76], [-0.24, 0.08, -0.20], [-0.05, 1.00, -0.03], [0.01, 0.01, -0.11]]),
                    (12, [[9.85, 0.78, 8.39], [3.33, 0.59, 5.33], [0.84, 0.01, -0.53], [-0.15, -0.01, -0.22], [-0.06, 1.00, -0.07], [-0.00, -0.00, -0.03]]),
                    (13, [[12.56, 1.40, 14.10], [2.03, 0.65, 6.11], [0.95, 0.03, -0.32], [-0.08, -0.03, -0.22], [-0.06, 1.00, -0.09], [-0.01, -0.00, -0.01]]),
                    (14, [[13.84, 2.07, 20.50], [0.53, 0.68, 6.63], [1.00, 0.07, -0.09], [-0.02, -0.05, -0.24], [-0.08, 0.99, -0.10], [-0.02, -0.01, -0.01]])]),
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
        return '2D Diverging Bifurcation 1'

    @staticmethod
    def getParameterSetNames():
        return ['Default']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        networkLayoutOption = cls.parameterSetStructureStrings['Default']
        options = {
            'Network layout': copy.deepcopy(networkLayoutOption),
            'Target element length': 4.0,
            'Number of elements around': 8,
            'Use cross derivatives': False
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Network layout',
            'Target element length',
            'Number of elements around',
            'Use cross derivatives']
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
        if options['Number of elements around'] < 4:
            options['Number of elements around'] = 4
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
        targetElementLength = options['Target element length']

        fm = region.getFieldmodule()
        fm.beginChange()
        mesh = fm.findMeshByDimension(2)
        coordinates = findOrCreateFieldCoordinates(fm)

        firstNodeIdentifier = 1
        elementIdentifier = 1

        # Create annotation groups
        parentGroup = AnnotationGroup(region, ("parent", "None"))
        child1Group = AnnotationGroup(region, ("child 1", "None"))
        child2Group = AnnotationGroup(region, ("child 2", "None"))
        bifurcationGroup = AnnotationGroup(region, ("diverging bifurcation", "None"))
        annotationGroups = [parentGroup, child1Group, child2Group, bifurcationGroup]

        parentMeshGroup = parentGroup.getMeshGroup(mesh)
        child1MeshGroup = child1Group.getMeshGroup(mesh)
        child2MeshGroup = child2Group.getMeshGroup(mesh)
        bifurcationMeshGroup = bifurcationGroup.getMeshGroup(mesh)

        # Geometric coordinates
        geometricNetworkLayout = BifurcationNetworkLayout(region, networkLayout, targetElementLength)

        parentLength = geometricNetworkLayout.arcLengthOfGroupsAlong[0]
        child1Length = geometricNetworkLayout.arcLengthOfGroupsAlong[1]
        child2Length = geometricNetworkLayout.arcLengthOfGroupsAlong[2]

        elementsCountInParent = math.ceil(parentLength / targetElementLength)
        elementsCountInChild1 = math.ceil(child1Length / targetElementLength)
        elementsCountInChild2 = math.ceil(child2Length / targetElementLength)

        cx_parent_group = geometricNetworkLayout.cxGroups[0]
        cx_child1_group = geometricNetworkLayout.cxGroups[1]
        cx_child2_group = geometricNetworkLayout.cxGroups[2]

        sx_parent_group = geometricNetworkLayout.sxGroups[0]
        sx_child1_group = geometricNetworkLayout.sxGroups[1]
        sx_child2_group = geometricNetworkLayout.sxGroups[2]

        # Get parent nodes
        parentCoordinates = getCoordinatesAlongTube2D(cx_parent_group, elementsCountAround, elementsCountInParent,
                                                      startRadian=-math.pi/2)

        parentLastRingNodeCoordinates = getTargetedRingNodesCoordinates2D(parentCoordinates, elementsCountAround,
                                                                          elementsCountInParent, omitStartRows=0,
                                                                          omitEndRows=1)

        # Get child1 nodes
        child1Coordinates = getCoordinatesAlongTube2D(cx_child1_group, elementsCountAround, elementsCountInChild1,
                                                      startRadian=-math.pi/2)

        child1FirstRingNodeCoordinates = getTargetedRingNodesCoordinates2D(child1Coordinates, elementsCountAround,
                                                                           elementsCountInChild1, omitStartRows=1,
                                                                           omitEndRows=0)

        # Get child2 nodes
        child2Coordinates = getCoordinatesAlongTube2D(cx_child2_group, elementsCountAround,
                                                      elementsCountInChild2, startRadian=-math.pi/2)

        child2FirstRingNodeCoordinates = getTargetedRingNodesCoordinates2D(child2Coordinates, elementsCountAround,
                                                                           elementsCountInChild2, omitStartRows=1,
                                                                           omitEndRows=0)

        # Create nodes
        # Create parent nodes
        nodeIdentifier = generateTubeNodes2D(fm, firstNodeIdentifier, parentCoordinates, elementsCountInParent,
                                             elementsCountAround, omitStartRows=0, omitEndRows=1, startNodes=None)

        # Create bifurcation nodes
        paCentre = sx_parent_group[0][1]
        c1Centre = sx_child2_group[0][-2]
        c2Centre = sx_child1_group[0][-2]
        paxList = parentLastRingNodeCoordinates[0]
        pad2 = parentLastRingNodeCoordinates[2]
        c1xList = child2FirstRingNodeCoordinates[0]
        c1d2 = child2FirstRingNodeCoordinates[2]
        c2xList = child1FirstRingNodeCoordinates[0]
        c2d2 = child1FirstRingNodeCoordinates[2]
        nodeIdentifier, roNodeId, coNodeId, nextNodeId, paStartIndex, c1StartIndex, c2StartIndex = \
            create2DBifurcationNodes(fm, nodeIdentifier, paCentre, paxList, pad2, c1Centre, c1xList, c1d2, c2Centre,
                                     c2xList, c2d2)

        # Create child1 nodes
        nodeIdentifier = generateTubeNodes2D(fm, nodeIdentifier, child1Coordinates, elementsCountInChild1,
                                             elementsCountAround, omitStartRows=1, omitEndRows=0, startNodes=None)

        # Create child2 nodes
        nodeIdentifier = generateTubeNodes2D(fm, nodeIdentifier, child2Coordinates, elementsCountInChild2,
                                             elementsCountAround, omitStartRows=1, omitEndRows=0, startNodes=None)

        # Create elements
        # Create parent elements
        startNodeId = firstNodeIdentifier
        elementIdentifier = generate2DTubeElements(fm, startNodeId, elementIdentifier, elementsCountInParent,
                                                   elementsCountAround, omitStartRows=1, omitEndRows=0,
                                                   meshGroups=[parentMeshGroup, bifurcationMeshGroup])

        # Create bifurcation elements
        parentLastRingNodeId, nodeCount = getTargetedRingNodesId2D(firstNodeIdentifier, elementsCountAround,
                                                                   elementsCountInParent, omitStartRows=0,
                                                                   omitEndRows=1)
        child1FirstRingNodeId, nodeCount = getTargetedRingNodesId2D(nextNodeId, elementsCountAround,
                                                                    elementsCountInChild1, omitStartRows=1,
                                                                    omitEndRows=0)
        child2FirstRingNodeId, nodeCount = getTargetedRingNodesId2D(nodeCount, elementsCountAround,
                                                                    elementsCountInChild2, omitStartRows=1,
                                                                    omitEndRows=0)

        paNodeId = parentLastRingNodeId
        c1NodeId = child2FirstRingNodeId
        c2NodeId = child1FirstRingNodeId
        elementIdentifier = make_tube_bifurcation_elements_2d(region, coordinates, elementIdentifier, paNodeId,
                                                              paStartIndex, c1NodeId, c1StartIndex, c2NodeId,
                                                              c2StartIndex, roNodeId, coNodeId,
                                                              meshGroups=[parentMeshGroup, child2MeshGroup,
                                                                          child1MeshGroup, bifurcationMeshGroup])

        # Create child1 elements
        startNodeId = child1FirstRingNodeId[0]
        elementIdentifier = generate2DTubeElements(fm, startNodeId, elementIdentifier, elementsCountInChild1,
                                                   elementsCountAround, omitStartRows=1, omitEndRows=0,
                                                   meshGroups=[child1MeshGroup, bifurcationMeshGroup])

        # Create child2 elements
        startNodeId = child2FirstRingNodeId[0]
        elementIdentifier = generate2DTubeElements(fm, startNodeId, elementIdentifier, elementsCountInChild2,
                                                   elementsCountAround, omitStartRows=0, omitEndRows=1,
                                                   meshGroups=[child2MeshGroup, bifurcationMeshGroup])

        fm.endChange()
        return annotationGroups, None


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
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
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


def getCoordinatesAlongTube2D(cx_group, elementsCountAround, elementsCountAlongTube, startRadian):

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
        d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        xSampledTube.append(xAround)
        d1SampledTube.append(d1Smoothed)
        d2SampledTube.append(d2Around)

    xList = []
    d1List = []
    d2List = []
    for n2 in range(elementsCountAlongTube + 1):
        for n1 in range(elementsCountAround):
            xList.append(xSampledTube[n2][n1])
            d1List.append(d1SampledTube[n2][n1])
            d2List.append(d2SampledTube[n2][n1])

    coordinatesList = [xList, d1List, d2List]

    return coordinatesList


def getTargetedRingNodesCoordinates2D(tubeCoordinates, elementsCountAround, elementsCountAlongTube, omitStartRows,
                                      omitEndRows):

    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        for n1 in range(elementsCountAround):
            n = n2 * elementsCountAround + n1
            x = tubeCoordinates[0][n]
            d1 = tubeCoordinates[1][n]
            d2 = tubeCoordinates[2][n]
            if omitEndRows == 1:  # merging to the bifurcation
                if n2 == elementsCountAlongTube:
                    pass
                else:
                    if n2 == elementsCountAlongTube - 1:
                        xLastRing.append(x)
                        d1LastRing.append(d1)
                        d2LastRing.append(d2)
            elif omitStartRows == 1:  # diverging from bifurcation
                if n2 == 0:
                    pass
                else:
                    if n2 == 1:
                        xFirstRing.append(x)
                        d1FirstRing.append(d1)
                        d2FirstRing.append(d2)

    if omitStartRows == 1:
        targetedRingCoordinates = [xFirstRing, d1FirstRing, d2FirstRing]
    elif omitEndRows == 1:
        targetedRingCoordinates = [xLastRing, d1LastRing, d2LastRing]
    else:
        targetedRingCoordinates = []

    return targetedRingCoordinates


def generateTubeNodes2D(fm, nodeIdentifier, tubeCoordinates, elementsCountAlongTube, elementsCountAround,
                        omitStartRows, omitEndRows, startNodes=None):

    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

    # Create tube nodes
    lastRingsNodeId = []
    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    firstRingsNodeId = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        lastRingNodeIdThroughWall = []
        xLastRingThroughWall = []
        d1LastRingThroughWall = []
        d2LastRingThroughWall = []
        firstRingNodeIdThroughWall = []
        xFirstRingThroughWall = []
        d1FirstRingThroughWall = []
        d2FirstRingThroughWall = []
        for n1 in range(elementsCountAround):
            n = n2 * elementsCountAround + n1
            x = tubeCoordinates[0][n]
            d1 = tubeCoordinates[1][n]
            d2 = tubeCoordinates[2][n]
            if omitEndRows == 1:  # merging to the bifurcation
                if n2 == elementsCountAlongTube:
                    pass
                else:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    if n2 == elementsCountAlongTube - 1:
                        lastRingNodeIdThroughWall.append(nodeIdentifier)
                        xLastRingThroughWall.append(x)
                        d1LastRingThroughWall.append(d1)
                        d2LastRingThroughWall.append(d2)
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
                    if n2 == 1:
                        firstRingNodeIdThroughWall.append(nodeIdentifier)
                        xFirstRingThroughWall.append(x)
                        d1FirstRingThroughWall.append(d1)
                        d2FirstRingThroughWall.append(d2)
                    nodeIdentifier += 1
            else:
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                nodeIdentifier += 1
        if omitEndRows == 1:
            if n2 == elementsCountAlongTube - 1:
                lastRingsNodeId.append(lastRingNodeIdThroughWall)
                xLastRing.append(xLastRingThroughWall)
                d1LastRing.append(d1LastRingThroughWall)
                d2LastRing.append(d2LastRingThroughWall)
        elif omitStartRows == 1:
            if n2 == 1:
                firstRingsNodeId.append(firstRingNodeIdThroughWall)
                xFirstRing.append(xFirstRingThroughWall)
                d1FirstRing.append(d1FirstRingThroughWall)
                d2FirstRing.append(d2FirstRingThroughWall)

    return nodeIdentifier


def getTargetedRingNodesId2D(nodeCount, elementsCountAround, elementsCountAlongTube, omitStartRows, omitEndRows):

    lastRingsNodeId = []
    firstRingsNodeId = []
    for n2 in range(elementsCountAlongTube + 1):
        for n1 in range(elementsCountAround):
            if omitEndRows == 1:  # merging to the bifurcation
                if n2 == elementsCountAlongTube:
                    pass
                else:
                    if n2 == elementsCountAlongTube - 1:
                        lastRingsNodeId.append(nodeCount)
                    nodeCount += 1
            elif omitStartRows == 1:  # diverging from bifurcation
                if n2 == 0:
                    pass
                else:
                    if n2 == 1:
                        firstRingsNodeId.append(nodeCount)
                    nodeCount += 1
            else:
                nodeCount += 1

    if omitStartRows == 1:
        targetedRingNodeId = firstRingsNodeId
    elif omitEndRows == 1:
        targetedRingNodeId = lastRingsNodeId
    else:
        targetedRingNodeId = []

    return targetedRingNodeId, nodeCount


def create2DBifurcationNodes(fm, nodeIdentifier, paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):

    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

    rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex = \
        make_tube_bifurcation_points(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2)

    # Create bifurcation nodes
    roNodeId = []
    for n in range(len(rox)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rox[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2[n])
        roNodeId.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1
    coNodeId = []
    for n in range(len(cox)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cox[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2[n])
        coNodeId.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1

    nextNodeId = nodeIdentifier
    return nodeIdentifier, roNodeId, coNodeId, nextNodeId, paStartIndex, c1StartIndex, c2StartIndex


def generate2DTubeElements(fm, startNodeId, elementIdentifier, elementsCountAlongTube, elementsCountAround,
                           omitStartRows, omitEndRows, startNodes=None, meshGroups=None):

    mesh = fm.findMeshByDimension(2)
    coordinates = findOrCreateFieldCoordinates(fm)
    bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
    for n in range(4):
        eft.setFunctionNumberOfTerms(n * 4 + 4, 0)
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.defineField(coordinates, -1, eft)
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
    # Create tube elements
    for e2 in range(elementsCountAlongTube - 1 if omitStartRows or omitEndRows == 1 else elementsCountAlongTube):
        for e1 in range(elementsCountAround):
            bni1 = e2 * elementsCountAround + e1 + startNodeId
            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + startNodeId
            bni3 = bni1 + elementsCountAround
            bni4 = bni2 + elementsCountAround
            nodeIdentifiers = [bni1, bni2, bni3, bni4]
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

    return elementIdentifier
