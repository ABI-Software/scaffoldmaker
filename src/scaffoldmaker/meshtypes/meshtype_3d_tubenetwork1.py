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
            "Annotation elements counts around": [0],
            "Target element density along longest segment": 4.0,
            "Serendipity": True,
            "Show trim surfaces": False
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            "Network layout",
            "Elements count around",
            "Elements count through wall",
            "Annotation elements counts around",
            "Target element density along longest segment",
            "Serendipity",
            "Show trim surfaces"
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
            options["Network layout"] = cls.getOptionScaffoldPackage("Network layout", MeshType_1d_network_layout1)
        elementsCountAround = options["Elements count around"]
        if options["Elements count around"] < 4:
            options["Elements count around"] = 4
        if options["Elements count through wall"] < 1:
            options["Elements count through wall"] = 1
        annotationElementsCountsAround = options["Annotation elements counts around"]
        if len(annotationElementsCountsAround) == 0:
            options["Annotation elements count around"] = [0]
        else:
            for i in range(len(annotationElementsCountsAround)):
                if annotationElementsCountsAround[i] <= 0:
                    annotationElementsCountsAround[i] = 0
                elif annotationElementsCountsAround[i] < 4:
                    annotationElementsCountsAround[i] = 4
        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0
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
        elementsCountAround = options["Elements count around"]
        elementsCountThroughWall = options["Elements count through wall"]
        annotationElementsCountsAround = options["Annotation elements counts around"]
        targetElementDensityAlongLongestSegment = options["Target element density along longest segment"]
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

        nodeIdentifier, elementIdentifier, annotationGroups = generateTubeBifurcationTree(
            networkMesh, region, coordinates, nodeIdentifier, elementIdentifier,
            elementsCountAround, targetElementDensityAlongLongestSegment, elementsCountThroughWall,
            layoutAnnotationGroups, annotationElementsCountsAround,
            serendipity=serendipity, showTrimSurfaces=options["Show trim surfaces"])

        return annotationGroups, None
