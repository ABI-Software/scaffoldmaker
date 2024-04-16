"""
Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
"""
from cmlibs.maths.vectorops import add, sub
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabelVersion, setEftScaleFactorIds
from scaffoldmaker.utils.interpolation import getCubicHermiteCurvesLength, interpolateSampleCubicHermite, \
    sampleCubicHermiteCurvesSmooth
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters
import math


class BoxSegmentData:
    """
    Definition of box mesh for network segment.
    """

    def __init__(self, pathParameters):
        """
        :param pathParameters: (sx[], sd1[], sd2[], sd12[], sd3[], sd13[])
        """
        self._pathParameters = pathParameters  # original parameters from network layout
        self._segmentLength = getCubicHermiteCurvesLength(pathParameters[0], pathParameters[1])
        self._sampledBoxCoordinates = None
        self._sampledNodeIdentifiers = []
        self._annotationMeshGroups = []

    def getPathParameters(self):
        return self._pathParameters

    def getSegmentLength(self):
        return self._segmentLength

    def getSampledBoxCoordinates(self):
        return self._sampledBoxCoordinates

    def setSampledBoxCoordinates(self, sampledBoxCoordinates):
        """
        :param sampledBoxCoordinates: (sx[], sd1[], sd2[], sd12[], sd3[], sd13[])
        """
        self._sampledBoxCoordinates = sampledBoxCoordinates
        self._sampledNodeIdentifiers = [None] * len(self._sampledBoxCoordinates[0])

    def getSampledElementsCountAlong(self):
        """
        Must have previously called setSampledBoxCoordinates
        :return: Number of elements along sampled box.
        """
        return len(self._sampledBoxCoordinates[0]) - 1

    def getSampledNodeIdentifier(self, nodeIndex):
        """
        Get start node identifier for supplying to adjacent straight or junction.
        :param nodeIndex: Index of sampled node from 0 to sampledElementsCountAlong, or -1 for last.
        """
        return self._sampledNodeIdentifiers[nodeIndex]

    def setSampledNodeIdentifier(self, nodeIndex, nodeIdentifier):
        """
        :param nodeIndex: Index of sampled node from 0 to sampledElementsCountAlong, or -1 for last.
        :param nodeIdentifier: Identifier of node to set at nodeIndex.
        """
        self._sampledNodeIdentifiers[nodeIndex] = nodeIdentifier

    def addAnnotationMeshGroup(self, annotationMeshGroup):
        """
        Add an annotation mesh group for segment elements to be added to.
        :param annotationMeshGroup: Mesh group to add.
        """
        self._annotationMeshGroups.append(annotationMeshGroup)

    def getAnnotationMeshGroups(self):
        return self._annotationMeshGroups


class BoxJunctionData:
    """
    Describes junction between multiple segments, some in, some out.
    """

    def __init__(self, networkSegmentsIn: list, networkSegmentsOut: list, boxSegmentData):
        """
        :param networkSegmentsIn: List of input segments.
        :param networkSegmentsOut: List of output segments.
        :param boxSegmentData: dict NetworkSegment -> BoxSegmentData.
        """
        self._networkSegmentsIn = networkSegmentsIn
        self._networkSegmentsOut = networkSegmentsOut
        self._networkSegments = networkSegmentsIn + networkSegmentsOut
        segmentsCount = len(self._networkSegments)
        assert segmentsCount > 0
        self._boxData = [boxSegmentData[networkSegment] for networkSegment in self._networkSegments]
        self._segmentsIn = [self._networkSegments[s] in self._networkSegmentsIn for s in range(segmentsCount)]
        self._nodeIdentifier = None

    def getSegmentsIn(self):
        return self._segmentsIn

    def getBoxData(self):
        return self._boxData

    def getNodeIdentifier(self):
        return self._nodeIdentifier

    def setNodeIdentifier(self, nodeIdentifier):
        self._nodeIdentifier = nodeIdentifier
        segmentsCount = len(self._networkSegments)
        for s in range(segmentsCount):
            nodeIndex = -1 if self._segmentsIn[s] else 0
            self._boxData[s].setSampledNodeIdentifier(nodeIndex, nodeIdentifier)


class MeshType_3d_boxnetwork1(Scaffold_base):
    """
    Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
    """

    @classmethod
    def getName(cls):
        return "3D Box Network 1"

    @classmethod
    def getParameterSetNames(cls):
        return MeshType_1d_network_layout1.getParameterSetNames()

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            "Network layout": ScaffoldPackage(MeshType_1d_network_layout1, defaultParameterSetName=parameterSetName),
            "Target element density along longest segment": 4.0
        }
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Network layout",
            "Target element density along longest segment"
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Network layout":
            return [MeshType_1d_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + ".getOptionScaffoldTypeParameterSetNames.  " + \
            "Invalid option \"" + optionName + "\" scaffold type " + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()  # use the defaults from the network layout scaffold

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        """
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        """
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                "Invalid parameter set " + str(parameterSetName) + " for scaffold " + str(scaffoldType.getName()) + \
                " in option " + str(optionName) + " of scaffold " + cls.getName()
        if optionName == "Network layout":
            if not parameterSetName:
                parameterSetName = "Default"
            return ScaffoldPackage(scaffoldType, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        if not options["Network layout"].getScaffoldType() in cls.getOptionValidScaffoldTypes("Network layout"):
            options["Network layout"] = cls.getOptionScaffoldPackage("Network layout", MeshType_1d_network_layout1)
        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0
        dependentChanges = False
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        networkLayout = options["Network layout"]
        targetElementDensityAlongLongestSegment = options["Target element density along longest segment"]

        layoutRegion = region.createRegion()
        layoutFieldmodule = layoutRegion.getFieldmodule()
        layoutNodes = layoutFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        layoutMesh = layoutFieldmodule.findMeshByDimension(1)
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        layoutCoordinates = findOrCreateFieldCoordinates(layoutFieldmodule)

        networkMesh = networkLayout.getConstructionObject()

        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeIdentifier = 1
        nodetemplates = {}  # map from derivativeVersionsCount to nodetemplate

        mesh = fieldmodule.findMeshByDimension(3)
        elementIdentifier = 1
        hermiteBilinearBasis = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        hermiteBilinearBasis.setFunctionType(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        elementtemplates = {}  # map from (startVersion, endVersion) to (elementtemplate, eft)

        # make box annotation groups from network layout annotations
        annotationGroups = []
        layoutAnnotationMeshGroupMap = []  # List of tuples of (layout annotation mesh group, box mesh group)
        for layoutAnnotationGroup in layoutAnnotationGroups:
            if layoutAnnotationGroup.getDimension() == 1:
                annotationGroup = AnnotationGroup(region, layoutAnnotationGroup.getTerm())
                annotationGroups.append(annotationGroup)
                layoutAnnotationMeshGroupMap.append(
                    (layoutAnnotationGroup.getMeshGroup(layoutMesh), annotationGroup.getMeshGroup(mesh)))

        valueLabels = [
            Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
            Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

        networkSegments = networkMesh.getNetworkSegments()
        boxSegmentData = {}  # map from NetworkSegment to BoxSegmentData
        longestSegmentLength = 0.0
        for networkSegment in networkSegments:
            pathParameters = get_nodeset_path_ordered_field_parameters(
                layoutNodes, layoutCoordinates, valueLabels,
                networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())
            boxSegmentData[networkSegment] = boxData = BoxSegmentData(pathParameters)
            segmentLength = boxData.getSegmentLength()
            if segmentLength > longestSegmentLength:
                longestSegmentLength = segmentLength
            for layoutAnnotationMeshGroup, annotationMeshGroup in layoutAnnotationMeshGroupMap:
                if networkSegment.hasLayoutElementsInMeshGroup(layoutAnnotationMeshGroup):
                    boxData.addAnnotationMeshGroup(annotationMeshGroup)
        if longestSegmentLength > 0.0:
            targetElementLength = longestSegmentLength / targetElementDensityAlongLongestSegment
        else:
            targetElementLength = 1.0

        # map from NetworkNodes to BoxJunctionData
        boxNodeJunctionDataMap = {}
        for networkSegment in networkSegments:
            boxData = boxSegmentData[networkSegment]
            segmentNodes = networkSegment.getNetworkNodes()
            startSegmentNode = segmentNodes[0]
            startBoxJunctionData = boxNodeJunctionDataMap.get(startSegmentNode)
            if not startBoxJunctionData:
                startInSegments = startSegmentNode.getInSegments()
                startOutSegments = startSegmentNode.getOutSegments()
                startBoxJunctionData = BoxJunctionData(startInSegments, startOutSegments, boxSegmentData)
                boxNodeJunctionDataMap[startSegmentNode] = startBoxJunctionData
            endSegmentNode = segmentNodes[-1]
            endBoxJunctionData = boxNodeJunctionDataMap.get(endSegmentNode)
            if not endBoxJunctionData:
                endInSegments = endSegmentNode.getInSegments()
                endOutSegments = endSegmentNode.getOutSegments()
                endBoxJunctionData = BoxJunctionData(endInSegments, endOutSegments, boxSegmentData)
                boxNodeJunctionDataMap[endSegmentNode] = endBoxJunctionData
            segmentLength = boxData.getSegmentLength()
            elementsCountAlong = max(1, math.ceil(segmentLength / targetElementLength))
            loop = (len(startSegmentNode.getInSegments()) == 1) and \
                   (startSegmentNode.getInSegments()[0] is networkSegment) and \
                   (networkSegment.getNodeVersions()[0] == networkSegment.getNodeVersions()[-1])
            if (elementsCountAlong < 2) and loop:
                elementsCountAlong = 2  # at least 2 segments around loop
            pathParameters = boxData.getPathParameters()
            sx, sd1, pe, pxi, psf = sampleCubicHermiteCurvesSmooth(
                pathParameters[0], pathParameters[1], elementsCountAlong)
            sd2, sd12 = interpolateSampleCubicHermite(pathParameters[2], pathParameters[3], pe, pxi, psf)
            sd3, sd13 = interpolateSampleCubicHermite(pathParameters[4], pathParameters[5], pe, pxi, psf)
            boxData.setSampledBoxCoordinates((sx, sd1, sd2, sd12, sd3, sd13))

        lastNodeIdentifier = None
        for networkSegment in networkSegments:
            boxData = boxSegmentData[networkSegment]
            segmentNodes = networkSegment.getNetworkNodes()
            layoutNodeVersions = networkSegment.getNodeVersions()
            sx, sd1, sd2, sd12, sd3, sd13 = boxData.getSampledBoxCoordinates()
            elementsCountAlong = boxData.getSampledElementsCountAlong()
            annotationMeshGroups = boxData.getAnnotationMeshGroups()
            for n in range(elementsCountAlong + 1):
                thisNodeIdentifier = None
                versionsCount = 1
                version = 1
                boxJunctionData = None
                if n in [0, elementsCountAlong]:
                    nLayout = 0 if (n == 0) else -1
                    segmentNode = segmentNodes[nLayout]
                    boxJunctionData = boxNodeJunctionDataMap[segmentNode]
                    thisNodeIdentifier = boxJunctionData.getNodeIdentifier()
                    version = layoutNodeVersions[nLayout]
                    versionsCount = segmentNode.getVersionsCount()
                if thisNodeIdentifier is None:
                    nodetemplate = nodetemplates.get(versionsCount)
                    if not nodetemplate:
                        nodetemplate = nodes.createNodetemplate()
                        nodetemplate.defineField(coordinates)
                        for valueLabel in valueLabels[1:]:
                            nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, versionsCount)
                        nodetemplates[versionsCount] = nodetemplate
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    if boxJunctionData:
                        boxJunctionData.setNodeIdentifier(nodeIdentifier)
                    thisNodeIdentifier = nodeIdentifier
                    nodeIdentifier += 1
                else:
                    node = nodes.findNodeByIdentifier(thisNodeIdentifier)
                # note following will set shared versions of coordinates or derivatives multiple times
                # future: average derivatives from all adjoining segments
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, version, sd1[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, version, sd2[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, sd12[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, version, sd3[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, sd13[n])

                if n > 0:
                    startVersion = layoutNodeVersions[0] if (n == 1) else 1
                    templateKey = (startVersion, version)
                    elementtemplate_and_eft = elementtemplates.get(templateKey)
                    if elementtemplate_and_eft:
                        elementtemplate, eft = elementtemplate_and_eft
                    else:
                        eft = mesh.createElementfieldtemplate(hermiteBilinearBasis)
                        setEftScaleFactorIds(eft, [1], [])
                        ln = 1
                        for n3 in range(2):
                            s3 = [1] if (n3 == 0) else []
                            for n2 in range(2):
                                s2 = [1] if (n2 == 0) else []
                                for n1 in range(2):
                                    v = startVersion if (n1 == 0) else version
                                    remapEftNodeValueLabelVersion(
                                        eft, [ln], Node.VALUE_LABEL_VALUE,
                                        [(Node.VALUE_LABEL_VALUE, 1, []),
                                         (Node.VALUE_LABEL_D_DS2, v, s2),
                                         (Node.VALUE_LABEL_D_DS3, v, s3)])
                                    remapEftNodeValueLabelVersion(
                                        eft, [ln], Node.VALUE_LABEL_D_DS1,
                                        [(Node.VALUE_LABEL_D_DS1, v, []),
                                         (Node.VALUE_LABEL_D2_DS1DS2, v, s2),
                                         (Node.VALUE_LABEL_D2_DS1DS3, v, s3)])
                                    ln += 1
                        ln_map = [1, 2, 1, 2, 1, 2, 1, 2]
                        remapEftLocalNodes(eft, 2, ln_map)
                        elementtemplate = mesh.createElementtemplate()
                        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                        elementtemplate.defineField(coordinates, -1, eft)
                        elementtemplates[templateKey] = (elementtemplate, eft)

                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    nids = [lastNodeIdentifier, thisNodeIdentifier]
                    element.setNodesByIdentifier(eft, nids)
                    element.setScaleFactors(eft, [-1.0])
                    for annotationMeshGroup in annotationMeshGroups:
                        annotationMeshGroup.addElement(element)
                    elementIdentifier += 1

                lastNodeIdentifier = thisNodeIdentifier

        return annotationGroups, None
