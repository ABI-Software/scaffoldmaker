"""
Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
"""

import copy

from cmlibs.maths.vectorops import add, sub
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.networkmesh import NetworkMesh


class MeshType_3d_boxnetwork1(Scaffold_base):
    """
    Generates a hermite x bilinear 3-D box network mesh from a 1-D network layout.
    """

    @staticmethod
    def getName():
        return "3D Box Network 1"

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            "Network layout": ScaffoldPackage(MeshType_1d_network_layout1)
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            "Network layout"
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
            options["Network layout"] = cls.getOptionScaffoldPackage("Network layout")
        dependentChanges = False
        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        networkLayout = options["Network layout"]

        layoutRegion = region.createRegion()
        layoutFieldmodule = layoutRegion.getFieldmodule()
        layoutNodes = layoutFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        layoutMesh = layoutFieldmodule.findMeshByDimension(1)
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        layoutCoordinates = findOrCreateFieldCoordinates(layoutFieldmodule)
        layoutFieldcache = layoutFieldmodule.createFieldcache()

        networkMesh = networkLayout.getConstructionObject()

        fieldmodule = region.getFieldmodule()
        fieldcache = fieldmodule.createFieldcache()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeIdentifier = 1
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)

        mesh = fieldmodule.findMeshByDimension(3)
        elementIdentifier = 1
        hermiteBilinearBasis = fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        hermiteBilinearBasis.setFunctionType(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(hermiteBilinearBasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        # make box annotation groups from network layout annotations
        annotationGroups = []
        layoutBoxMeshGroups = {}  # map from group name
        for layoutAnnotationGroup in layoutAnnotationGroups:
            if layoutAnnotationGroup.getDimension() == 1:
                annotationGroup = AnnotationGroup(region, layoutAnnotationGroup.getTerm())
                annotationGroups.append(annotationGroup)
                layoutBoxMeshGroups[layoutAnnotationGroup.getName()] = \
                    (layoutAnnotationGroup.getMeshGroup(layoutMesh), annotationGroup.getMeshGroup(mesh))

        networkSegments = networkMesh.getNetworkSegments()
        slices = {}  # map from network layout node identifier to list of 4 box node identifiers on slice
        for networkSegment in networkSegments:
            segmentNodes = networkSegment.getNetworkNodes()
            # segmentVersions = networkSegment.getNodeVersions()
            segmentElementIdentifiers = networkSegment.getElementIdentifiers()
            segmentSlices = []
            segmentNodeCount = len(segmentNodes)
            lastSlice = None
            for n in range(segmentNodeCount):
                segmentNode = segmentNodes[n]
                layoutNodeIdentifier = segmentNode.getNodeIdentifier()
                slice = slices.get(layoutNodeIdentifier)
                if slice:
                    segmentSlices.append(slice)
                    lastSlice = slice
                    continue
                layoutNode = layoutNodes.findNodeByIdentifier(layoutNodeIdentifier)
                layoutFieldcache.setNode(layoutNode)
                # currently only supports node version 1
                _, lx = layoutCoordinates.getNodeParameters(layoutFieldcache, -1, Node.VALUE_LABEL_VALUE, 1, 3)
                _, ld1 = layoutCoordinates.getNodeParameters(layoutFieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)
                _, ld2 = layoutCoordinates.getNodeParameters(layoutFieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)
                _, ld12 = layoutCoordinates.getNodeParameters(layoutFieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)
                _, ld3 = layoutCoordinates.getNodeParameters(layoutFieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)
                _, ld13 = layoutCoordinates.getNodeParameters(layoutFieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)

                slice = []
                for n3 in range(2):
                    for n2 in range(2):
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        x = lx
                        d1 = ld1
                        if n2 == 0:
                            x = sub(x, ld2)
                            d1 = sub(d1, ld12)
                        else:
                            x = add(x, ld2)
                            d1 = add(d1, ld12)
                        if n3 == 0:
                            x = sub(x, ld3)
                            d1 = sub(d1, ld13)
                        else:
                            x = add(x, ld3)
                            d1 = add(d1, ld13)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        slice.append(nodeIdentifier)
                        nodeIdentifier += 1
                slices[layoutNodeIdentifier] = slice

                if lastSlice:
                    element = mesh.createElement(elementIdentifier, elementtemplate);
                    nids = [lastSlice[0], slice[0],
                            lastSlice[1], slice[1],
                            lastSlice[2], slice[2],
                            lastSlice[3], slice[3]]
                    element.setNodesByIdentifier(eft, nids)
                    layoutElementIdentifier = segmentElementIdentifiers[n - 1]
                    layoutElement = layoutMesh.findElementByIdentifier(layoutElementIdentifier)
                    for layoutBoxMeshGroup in layoutBoxMeshGroups.values():
                        if layoutBoxMeshGroup[0].containsElement(layoutElement):
                            layoutBoxMeshGroup[1].addElement(element)
                    elementIdentifier += 1

                lastSlice = slice

        return annotationGroups, None
