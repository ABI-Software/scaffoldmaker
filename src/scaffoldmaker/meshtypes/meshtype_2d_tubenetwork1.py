"""
Generates a 2-D Hermite bifurcating tube network.
"""

import copy

from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.utils.zinc.finiteelement import get_maximum_element_identifier, get_maximum_node_identifier
from cmlibs.zinc.field import Field
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.bifurcation import generateTubeBifurcationTree2D


class MeshType_2d_tubenetwork1(Scaffold_base):
    """
    Generates a 2-D hermite bifurcating tube network.
    """

    @staticmethod
    def getName():
        return "2D Tube Network 1"

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            "Network layout": ScaffoldPackage(MeshType_1d_network_layout1),
            "Elements count around": 8,
            "Target element aspect ratio": 2.0
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            "Network layout",
            "Elements count around",
            "Target element aspect ratio"
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
        elementsCountAround = options["Elements count around"]
        options["Elements count around"] = max(4, elementsCountAround + (elementsCountAround % 2))
        if options["Target element aspect ratio"] < 0.01:
            options["Target element aspect ratio"] = 0.01
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
        elementsCountAround = options["Elements count around"]
        networkLayout = options["Network layout"]
        targetElementAspectRatio = options["Target element aspect ratio"]

        layoutRegion = region.createRegion()
        layoutFieldmodule = layoutRegion.getFieldmodule()
        # layoutNodes = layoutFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # layoutMesh = layoutFieldmodule.findMeshByDimension(1)
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        # layoutCoordinates = find_or_create_field_coordinates(layoutFieldmodule)
        # layoutFieldcache = layoutFieldmodule.createFieldcache()

        fieldmodule = region.getFieldmodule()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = find_or_create_field_coordinates(fieldmodule)
        nodeIdentifier = max(get_maximum_node_identifier(nodes), 0) + 1
        mesh = fieldmodule.findMeshByDimension(2)
        elementIdentifier = max(get_maximum_element_identifier(mesh), 0) + 1

        # make box annotation groups from network layout annotations
        annotationGroups = []
        layoutMeshGroups = {}  # map from group name
        for layoutAnnotationGroup in layoutAnnotationGroups:
            if layoutAnnotationGroup.getDimension() == 1:
                annotationGroup = AnnotationGroup(region, layoutAnnotationGroup.getTerm())
                annotationGroups.append(annotationGroup)
                layoutMeshGroups[layoutAnnotationGroup.getName()] = \
                    (layoutAnnotationGroup.getMeshGroup(layoutMesh), annotationGroup.getMeshGroup(mesh))

        networkMesh = networkLayout.getConstructionObject()

        nodeIdentifier, elementIdentifier = generateTubeBifurcationTree2D(
            networkMesh, region, coordinates, nodeIdentifier, elementIdentifier,
            elementsCountAround, targetElementAspectRatio, serendipity=False)

        return annotationGroups, None
