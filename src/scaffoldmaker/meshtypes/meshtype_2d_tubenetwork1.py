"""
Generates a 2-D Hermite bifurcating tube network.
"""
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData


class MeshType_2d_tubenetwork1(Scaffold_base):
    """
    Generates a 2-D hermite tube network.
    """

    @classmethod
    def getName(cls):
        return "2D Tube Network 1"

    @classmethod
    def getParameterSetNames(cls):
        return MeshType_1d_network_layout1.getParameterSetNames()

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            "Network layout": ScaffoldPackage(MeshType_1d_network_layout1, defaultParameterSetName=parameterSetName),
            "Number of elements around": 8,
            "Annotation numbers of elements around": [0],
            "Target element density along longest segment": 4.0,
            "Annotation numbers of elements along": [0],
            "Show trim surfaces": False
        }
        if parameterSetName in ["Loop", "Snake", "Vase"]:
            options["Target element density along longest segment"] = 12.0
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Network layout",
            "Number of elements around",
            "Annotation numbers of elements around",
            "Target element density along longest segment",
            "Annotation numbers of elements along",
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
            return ScaffoldPackage(scaffoldType, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        if not options["Network layout"].getScaffoldType() in cls.getOptionValidScaffoldTypes("Network layout"):
            options["Network layout"] = cls.getOptionScaffoldPackage("Network layout", MeshType_1d_network_layout1)
        dependentChanges = False
        if options["Number of elements around"] < 4:
            options["Number of elements around"] = 4
        annotationElementsCountsAround = options["Annotation numbers of elements around"]
        if len(annotationElementsCountsAround) == 0:
            options["Annotation numbers of elements around"] = [0]
        else:
            for i in range(len(annotationElementsCountsAround)):
                if annotationElementsCountsAround[i] <= 0:
                    annotationElementsCountsAround[i] = 0
                elif annotationElementsCountsAround[i] < 4:
                    annotationElementsCountsAround[i] = 4
        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0
        annotationAlongCounts = options["Annotation numbers of elements along"]
        if len(annotationAlongCounts) == 0:
            options["Annotation numbers of elements along"] = [0]
        else:
            for i in range(len(annotationAlongCounts)):
                if annotationAlongCounts[i] <= 0:
                    annotationAlongCounts[i] = 0
                elif annotationAlongCounts[i] < 1:
                    annotationAlongCounts[i] = 1
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base bicubic hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        layoutRegion = region.createRegion()
        networkLayout = options["Network layout"]
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        networkMesh = networkLayout.getConstructionObject()

        tubeNetworkMeshBuilder = TubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            layoutAnnotationGroups=networkLayout.getAnnotationGroups(),
            annotationElementsCountsAlong=options["Annotation numbers of elements along"],
            defaultElementsCountAround=options["Number of elements around"],
            annotationElementsCountsAround=options["Annotation numbers of elements around"],
            elementsCountThroughShell=1)
        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 2,
            isLinearThroughShell=True,
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        return annotationGroups, None
