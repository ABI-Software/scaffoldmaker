"""
Generates a 3-D Hermite bifurcating tube network.
"""

import copy

from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.utils.zinc.finiteelement import get_maximum_element_identifier, get_maximum_node_identifier
from cmlibs.zinc.field import Field
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.bifurcation import generateTubeBifurcationTree


class MeshType_3d_tubenetwork1(Scaffold_base):
    """
    Generates a 3-D hermite bifurcating tube network, with linear basis through wall.
    """

    @staticmethod
    def getName():
        return "3D Tube Network 1"

    @classmethod
    def getParameterSetNames(cls):
        return MeshType_1d_network_layout1.getParameterSetNames()

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            "Network layout": ScaffoldPackage(MeshType_1d_network_layout1,
                                              {"scaffoldSettings": {"Define inner coordinates": True}},
                                              defaultParameterSetName=parameterSetName),
            "Elements count around": 8,
            "Elements count through wall": 1,
            "Target element aspect ratio": 2.0,
            "Serendipity": True
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            "Network layout",
            "Elements count around",
            "Elements count through wall",
            "Target element aspect ratio",
            "Serendipity"
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
            return ScaffoldPackage(scaffoldType,
                                   {"scaffoldSettings": {"Define inner coordinates": True}},
                                   defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        if not options["Network layout"].getScaffoldType() in cls.getOptionValidScaffoldTypes("Network layout"):
            options["Network layout"] = cls.getOptionScaffoldPackage("Network layout")
        elementsCountAround = options["Elements count around"]
        options["Elements count around"] = max(4, elementsCountAround + (elementsCountAround % 2))
        if options["Elements count through wall"] < 1:
            options["Elements count through wall"] = 1
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
        elementsCountThroughWall = options["Elements count through wall"]
        networkLayout = options["Network layout"]
        targetElementAspectRatio = options["Target element aspect ratio"]
        serendipity = options["Serendipity"]

        layoutRegion = region.createRegion()
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()

        fieldmodule = region.getFieldmodule()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = find_or_create_field_coordinates(fieldmodule)
        nodeIdentifier = max(get_maximum_node_identifier(nodes), 0) + 1
        mesh = fieldmodule.findMeshByDimension(2)
        elementIdentifier = max(get_maximum_element_identifier(mesh), 0) + 1

        networkMesh = networkLayout.getConstructionObject()

        try:
            nodeIdentifier, elementIdentifier, annotationGroups = generateTubeBifurcationTree(
                networkMesh, region, coordinates, nodeIdentifier, elementIdentifier,
                elementsCountAround, targetElementAspectRatio, elementsCountThroughWall,
                layoutAnnotationGroups, serendipity=serendipity)
        except Exception as e:
            print(e, "\nException occurred while generating tube network: Please edit network layout")
            return [], None

        return annotationGroups, None
