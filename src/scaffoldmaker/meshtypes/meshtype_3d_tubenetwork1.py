"""
Generates a 3-D Hermite bifurcating tube network (with optional solid core).
"""
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData


def calculateElementsCountAcrossMinor(cls, options):
    # Calculate elementsCountAcrosMinor
    elementsCountAcrossMinor = int(((options["Elements count around"] - 4) / 4 -
                                    (options['Number of elements across core major'] / 2)) * 2 + 6)

    return elementsCountAcrossMinor

def setParametersToDefault(cls, options):
    options["Elements count around"] = 8
    options["Number of elements across core major"] = 4
    options["Number of elements across core transition"] = 1

    return options

def checkCoreParameters(cls, options, dependentChanges=False):
    # Check elements count around are ok
    if options["Elements count around"] < 8:
        dependentChanges = True
        options = setParametersToDefault(cls, options)
    if options["Elements count around"] % 4:
        options["Elements count around"] += 4 - options["Elements count around"] % 4

    # Check elements count across major are ok
    if options['Number of elements across core major'] < 4:
        options['Number of elements across core major'] = 4
    if options['Number of elements across core major'] % 2:
        options['Number of elements across core major'] += 1

    # Check elements count across transition
    if options['Number of elements across core transition'] < 1:
        options['Number of elements across core transition'] = 1

    # Calculate elementsCountAcrossMinor based on the current set of elementsCountAround and elementsCountAcrossMajor
    elementsCountAcrossMinor = calculateElementsCountAcrossMinor(cls, options)

    # Rcrit check
    Rcrit = max(
        min(options['Number of elements across core major'] - 4, elementsCountAcrossMinor - 4) // 2 + 1, 0)
    if Rcrit < options['Number of elements across core transition']:
        dependentChanges = True
        options['Number of elements across core transition'] = Rcrit

    # Number of elements around sanity check
    eM = options["Number of elements across core major"]
    em = elementsCountAcrossMinor
    eC = (eM - 1) * (em - 1) - ((eM - 3) * (em - 3))
    if options["Elements count around"] != eC:
        dependentChanges = True
        # Reset parameter values
        setParametersToDefault(cls, options)

    annotationElementsCountsAround = options["Annotation elements counts around"]
    if len(annotationElementsCountsAround) == 0:
        options["Annotation elements count around"] = [0]
    else:
        for i in range(len(annotationElementsCountsAround)):
            if annotationElementsCountsAround[i] <= 0:
                annotationElementsCountsAround[i] = 0
            elif annotationElementsCountsAround[i] < 4:
                annotationElementsCountsAround[i] = 4

    if options["Use linear through wall"]:
        dependentChanges = True
        options["Use linear through wall"] = False

    return options, dependentChanges


class MeshType_3d_tubenetwork1(Scaffold_base):
    """
    Generates a 3-D hermite tube network, with linear or cubic basis through wall.
    Number of elements generated are primarily controlled by "elements count around" and "elements count across major".
    The elements count around parameter determines the number of elements around a circular rim.
    The elements count across major determines the number of elements across the major axis of an ellipse (y-axis in the scaffold).
    The number of elements across minor axis (z-axis) is calculated internally based on the two elements count parameters.
    """

    @classmethod
    def getName(cls):
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
            "Use linear through wall": False,
            "Show trim surfaces": False,
            "Core": False,
            'Number of elements across core major': 4,
            'Number of elements across core transition': 1,
            "Annotation elements counts across major": [0]
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
            "Use linear through wall",
            "Show trim surfaces",
            "Core",
            'Number of elements across core major',
            'Number of elements across core transition',
            "Annotation elements counts across major"
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
        dependentChanges = False

        # Parameters specific to the core
        if options["Core"] == True:
            options, dependentChanges = checkCoreParameters(cls, options, dependentChanges)

        # When core option is false
        else:
            if options["Elements count around"] < 4:
                options["Elements count around"] = 4
            annotationElementsCountsAround = options["Annotation elements counts around"]
            if len(annotationElementsCountsAround) == 0:
                options["Annotation elements count around"] = [0]
            else:
                for i in range(len(annotationElementsCountsAround)):
                    if annotationElementsCountsAround[i] <= 0:
                        annotationElementsCountsAround[i] = 0
                    elif annotationElementsCountsAround[i] < 4:
                        annotationElementsCountsAround[i] = 4

        # Common parameters
        if options["Elements count through wall"] < 1:
            options["Elements count through wall"] = 1

        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic hermite or bicubic hermite-linear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        layoutRegion = region.createRegion()
        networkLayout = options["Network layout"]
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        tubeNetworkMeshBuilder = TubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            defaultElementsCountAround=options["Elements count around"],
            elementsCountThroughWall=options["Elements count through wall"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAround=options["Annotation elements counts around"],
            defaultElementsCountAcrossMajor=options['Number of elements across core major'],
            elementsCountTransition=options['Number of elements across core transition'],
            annotationElementsCountsAcrossMajor=options["Annotation elements counts across major"],
            isCore=options["Core"])

        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughWall=options["Use linear through wall"],
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        return annotationGroups, None
