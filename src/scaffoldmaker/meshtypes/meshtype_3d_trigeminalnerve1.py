"""
Generates a 3-D Hermite trigeminal nerve.
"""

from cmlibs.maths.vectorops import magnitude, set_magnitude # KM
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.trigeminal_nerve_terms import get_trigeminal_nerve_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters

def getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName):
    assert parameterSetName in cls.getParameterSetNames()
    if parameterSetName in ("Default", "Human 1"):
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            "scaffoldSettings": {
                "Structure": "1-2-3,3-4.1,4.2-5,5-6,4.3-7,7-8,4.4-9,9-10",
                "Define inner coordinates": True
            },
            "meshEdits": exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[0.000,0.000,0.000], [0.357,0.000,0.000], [0.000,0.170,0.000], [0.000,0.005,0.000], [0.000,0.000,0.072], [0.000,0.000,0.000]]), 
                (2, [[0.353,0.000,0.000], [0.349,0.000,0.000], [0.000,0.170,0.000], [0.000,-0.005,0.000], [0.000,0.000,0.072], [0.000,0.000,-0.000]]), 
                (3, [[0.698,0.000,0.000], [0.324,0.000,0.000], [0.000,0.160,0.000], [0.000,0.165,0.000], [0.000,0.000,0.072], [0.000,0.000,-0.000]]), 
                (4, [[1.000,0.000,0.000], [[0.281,0.000,0.000],[0.322,-0.321,0.000],[0.230,0.000,0.000],[0.110,0.122,0.000]], [[0.000,0.500,0.000],[0.353,0.354,0.000],[0.000,0.500,0.000],[-0.371,0.335,0.000]], [[-0.000,0.515,-0.000],[-0.301,-0.300,0.000],[0.048,-0.499,-0.017],[0.411,-0.293,0.000]], [[0.000,0.000,0.072],[-0.000,0.000,0.072],[0.000,0.000,0.072],[0.000,-0.000,0.072]], [[-0.000,-0.000,-0.000],[0.000,0.000,0.015],[-0.010,0.009,0.026],[0.000,0.000,0.011]]]), 
                (5, [[1.307,-0.305,0.000], [0.292,-0.289,0.000], [0.141,0.142,0.000], [-0.123,-0.124,0.000], [-0.000,0.000,0.072], [0.000,0.000,-0.015]]), 
                (6, [[1.584,-0.579,0.000], [0.262,-0.259,0.000], [0.105,0.107,0.000], [0.053,0.053,-0.000], [-0.000,0.000,0.042], [-0.000,-0.000,-0.045]]), 
                (7, [[1.411,-0.005,0.005], [0.500,0.000,0.000], [0.000,0.138,-0.009], [-0.006,-0.225,0.000], [-0.000,0.005,0.072], [0.003,0.000,-0.026]]), 
                (8, [[2.000,0.000,0.000], [0.770,0.000,0.000], [0.000,0.050,0.000], [0.007,0.049,0.018], [0.000,0.000,0.020], [-0.004,-0.009,-0.078]]), 
                (9, [[1.304,0.312,0.000], [0.498,0.502,0.000], [-0.106,0.106,0.000], [0.146,-0.136,0.000], [0.000,-0.000,0.072], [0.000,0.000,-0.011]]), 
                (10, [[2.000,1.000,0.000], [0.894,0.874,0.000], [-0.070,0.072,0.000], [-0.074,0.068,-0.000], [0.000,-0.000,0.050], [-0.000,-0.000,-0.033]])], [
                (1, [[0.000,0.000,0.000], [0.357,0.000,0.000], [0.000,0.128,0.000], [0.000,0.004,0.000], [0.000,0.000,0.054], [0.000,0.000,0.000]]), 
                (2, [[0.353,0.000,0.000], [0.349,0.000,0.000], [0.000,0.128,0.000], [0.000,-0.004,0.000], [0.000,0.000,0.054], [0.000,0.000,-0.000]]), 
                (3, [[0.698,0.000,0.000], [0.324,0.000,0.000], [0.000,0.120,0.000], [0.000,0.124,0.000], [0.000,0.000,0.054], [0.000,0.000,-0.000]]), 
                (4, [[1.000,0.000,0.000], [[0.281,0.000,0.000],[0.322,-0.321,0.000],[0.230,0.000,0.000],[0.110,0.122,0.000]], [[0.000,0.375,0.000],[0.265,0.266,0.000],[0.000,0.375,0.000],[-0.278,0.251,0.000]], [[-0.000,0.386,-0.000],[-0.225,-0.225,0.000],[0.036,-0.374,-0.013],[0.308,-0.219,0.000]], [[0.000,0.000,0.054],[-0.000,0.000,0.054],[0.000,0.000,0.054],[0.000,-0.000,0.054]], [[-0.000,-0.000,-0.000],[0.000,0.000,0.011],[-0.007,0.007,0.019],[0.000,0.000,0.008]]]), 
                (5, [[1.307,-0.305,0.000], [0.292,-0.289,0.000], [0.106,0.107,0.000], [-0.093,-0.093,0.000], [-0.000,0.000,0.054], [0.000,0.000,-0.011]]), 
                (6, [[1.584,-0.579,0.000], [0.262,-0.259,0.000], [0.079,0.080,0.000], [0.040,0.040,-0.000], [-0.000,0.000,0.032], [-0.000,-0.000,-0.034]]), 
                (7, [[1.411,-0.005,0.005], [0.500,0.000,0.000], [0.000,0.103,-0.007], [-0.005,-0.169,0.000], [-0.000,0.003,0.054], [0.002,0.000,-0.020]]), 
                (8, [[2.000,0.000,0.000], [0.770,0.000,0.000], [0.000,0.038,0.000], [0.005,0.037,0.013], [0.000,0.000,0.015], [-0.003,-0.007,-0.058]]), 
                (9, [[1.304,0.312,0.000], [0.498,0.502,0.000], [-0.080,0.079,0.000], [0.110,-0.102,0.000], [0.000,-0.000,0.054], [0.000,0.000,-0.008]]), 
                (10, [[2.000,1.000,0.000], [0.894,0.874,0.000], [-0.053,0.054,0.000], [-0.056,0.051,-0.000], [0.000,-0.000,0.038], [-0.000,-0.000,-0.025]])
                ]]),

            "userAnnotationGroups": [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-9',
                    'name': get_trigeminal_nerve_term('trigeminal nerve')[0],
                    'ontId': get_trigeminal_nerve_term('trigeminal nerve')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-3',
                    'name': get_trigeminal_nerve_term('trigeminal nerve root')[0],
                    'ontId': get_trigeminal_nerve_term('trigeminal nerve root')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4, 6, 8',
                    'name': get_trigeminal_nerve_term('trigeminal ganglion')[0],
                    'ontId': get_trigeminal_nerve_term('trigeminal ganglion')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5',
                    'name': get_trigeminal_nerve_term('mandibular nerve')[0],
                    'ontId': get_trigeminal_nerve_term('mandibular nerve')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_trigeminal_nerve_term('maxillary nerve')[0],
                    'ontId': get_trigeminal_nerve_term('maxillary nerve')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9',
                    'name': get_trigeminal_nerve_term('ophthalmic nerve')[0],
                    'ontId': get_trigeminal_nerve_term('ophthalmic nerve')[1]
                }
            ]
        })

class MeshType_3d_trigeminalnerve1(Scaffold_base):
    """
    Generates a 3-D hermite trigeminal nerve scaffold from a network layout description, with a solid core.
    The number of elements generated are primarily controlled by the "Number of elements around" and
    "Number of elements across core box minor" settings.
    The Number of elements around parameter determines the number of elements around the elliptical shell.
    The core is built out of a regular box with specified number of elements across core box minor.
    The minor direction is in the direction of d3 in the network layout it is defined from.
    The number of elements across the core box major axis is calculated internally from the above.
    The reason for preferring the minor number is that this must currently match when using the core,
    plus it is planned to eventually use an odd number to specify a half phase difference in starting node
    even without a core, as that is required to work with odd minor numbers in the core, once supported.
    There are any number of transition elements between the box and the shell forming further concentric rim elements.
    """

    @classmethod
    def getName(cls):
        return "3D Trigeminal Nerve 1"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {
            "Network layout": getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName),
            "Number of elements around": 8,
            "Number of elements through shell": 1,
            "Annotation numbers of elements around": [0],
            "Target element density along longest segment": 6.0,
            "Use linear through shell": False,
            "Show trim surfaces": False,
            "Number of elements across core box minor": 2,
            "Number of elements across core transition": 1,
            "Annotation numbers of elements across core box minor": [0]
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            "Network layout",
            "Number of elements around",
            "Number of elements through shell",
            "Annotation numbers of elements around",
            "Target element density along longest segment",
            "Use linear through shell",
            "Show trim surfaces",
            "Number of elements across core box minor",
            "Number of elements across core transition",
            "Annotation numbers of elements across core box minor"
        ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Network layout":
            return [MeshType_1d_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Network layout':
            return cls.getParameterSetNames()
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
            return getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        if not options["Network layout"].getScaffoldType() in cls.getOptionValidScaffoldTypes("Network layout"):
            options["Network layout"] = cls.getOptionScaffoldPackage("Network layout", MeshType_1d_network_layout1)
        dependentChanges = False

        if options["Number of elements around"] < 8:
            options["Number of elements around"] = 8
        elif options["Number of elements around"] % 4:
            options["Number of elements around"] += 4 - options["Number of elements around"] % 4

            annotationAroundCounts = options["Annotation numbers of elements around"]
            minAroundCount = options["Number of elements around"]
            if len(annotationAroundCounts) == 0:
                annotationAroundCounts = options["Annotation numbers of elements around"] = [0]
            else:
                for i in range(len(annotationAroundCounts)):
                    if annotationAroundCounts[i] <= 0:
                        annotationAroundCounts[i] = 0
                    else:
                        if annotationAroundCounts[i] < 8:
                            annotationAroundCounts[i] = 8
                        elif annotationAroundCounts[i] % 4:
                            annotationAroundCounts[i] += 4 - annotationAroundCounts[i]
                        if annotationAroundCounts[i] < minAroundCount:
                            minAroundCount = annotationAroundCounts[i]

            if options["Number of elements across core transition"] < 1:
                options["Number of elements across core transition"] = 1

            maxCoreBoxMinorCount = options["Number of elements around"] // 2 - 2
            if options["Number of elements across core box minor"] < 2:
                options["Number of elements across core box minor"] = 2
            elif options["Number of elements across core box minor"] > maxCoreBoxMinorCount:
                options["Number of elements across core box minor"] = maxCoreBoxMinorCount
                dependentChanges = True
            elif options["Number of elements across core box minor"] % 2:
                options["Number of elements across core box minor"] += 1

            annotationCoreBoxMinorCounts = options["Annotation numbers of elements across core box minor"]
            if len(annotationCoreBoxMinorCounts) == 0:
                annotationCoreBoxMinorCounts = options["Annotation numbers of elements across core box minor"] = [0]
            if len(annotationCoreBoxMinorCounts) > len(annotationAroundCounts):
                annotationCoreBoxMinorCounts = options["Annotation numbers of elements across core box minor"] = \
                    annotationCoreBoxMinorCounts[:len(annotationAroundCounts)]
                dependentChanges = True
            for i in range(len(annotationCoreBoxMinorCounts)):
                aroundCount = annotationAroundCounts[i] if annotationAroundCounts[i] \
                    else options["Number of elements around"]
                maxCoreBoxMinorCount = aroundCount // 2 - 2
                if annotationCoreBoxMinorCounts[i] <= 0:
                    annotationCoreBoxMinorCounts[i] = 0
                    # this may reduce the default
                    if maxCoreBoxMinorCount < options["Number of elements across core box minor"]:
                        options["Number of elements across core box minor"] = maxCoreBoxMinorCount
                        dependentChanges = True
                elif annotationCoreBoxMinorCounts[i] < 2:
                    annotationCoreBoxMinorCounts[i] = 2
                elif annotationCoreBoxMinorCounts[i] > maxCoreBoxMinorCount:
                    annotationCoreBoxMinorCounts[i] = maxCoreBoxMinorCount
                    dependentChanges = True
                elif annotationCoreBoxMinorCounts[i] % 2:
                    annotationCoreBoxMinorCounts[i] += 1

            if options["Use linear through shell"]:
                options["Use linear through shell"] = False
                dependentChanges = True

        else:
            if options["Number of elements around"] < 4:
                options["Number of elements around"] = 4
            annotationAroundCounts = options["Annotation numbers of elements around"]
            if len(annotationAroundCounts) == 0:
                options["Annotation numbers of elements around"] = [0]
            else:
                for i in range(len(annotationAroundCounts)):
                    if annotationAroundCounts[i] <= 0:
                        annotationAroundCounts[i] = 0
                    elif annotationAroundCounts[i] < 4:
                        annotationAroundCounts[i] = 4

        if options["Number of elements through shell"] < 1:
            options["Number of elements through shell"] = 1

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

        defaultAroundCount = options["Number of elements around"]
        coreTransitionCount = options["Number of elements across core transition"]
        defaultCoreBoxMinorCount = options["Number of elements across core box minor"]
        # implementation currently uses major count including transition
        defaultCoreMajorCount = defaultAroundCount // 2 - defaultCoreBoxMinorCount + 2 * coreTransitionCount
        annotationAroundCounts = options["Annotation numbers of elements around"]
        annotationCoreBoxMinorCounts = options["Annotation numbers of elements across core box minor"]
        annotationCoreMajorCounts = []
        for i in range(len(annotationCoreBoxMinorCounts)):
            aroundCount = annotationAroundCounts[i] if annotationAroundCounts[i] \
                else defaultAroundCount
            coreBoxMinorCount = annotationCoreBoxMinorCounts[i] if annotationCoreBoxMinorCounts[i] \
                else defaultCoreBoxMinorCount
            annotationCoreMajorCounts.append(aroundCount // 2 - coreBoxMinorCount + 2 * coreTransitionCount)

        tubeNetworkMeshBuilder = TubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            defaultElementsCountAround=defaultAroundCount,
            elementsCountThroughWall=options["Number of elements through shell"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAround=annotationAroundCounts,
            defaultElementsCountAcrossMajor=defaultCoreMajorCount,
            elementsCountTransition=coreTransitionCount,
            annotationElementsCountsAcrossMajor=annotationCoreMajorCounts,
            isCore=True)

        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughWall=options["Use linear through shell"],
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        return annotationGroups, None
