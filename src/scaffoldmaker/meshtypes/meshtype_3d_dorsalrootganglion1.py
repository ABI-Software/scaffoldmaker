"""
Generates a 3-D Hermite dorsal root ganglion.
"""
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.spinal_nerve_terms import get_spinal_nerve_term
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
                "Structure": "1-2.1,2.2-3-4-5, 2.3-6-7, 7-8-9, 9-10-11",
                "Define inner coordinates": True
            },
            "meshEdits": exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[-50.000,-72.000,940.000], [9.659,0.000,2.588], [-0.194,1.299,0.724], [0.000,0.000,0.000], [-0.336,-0.750,1.255], [0.000,0.000,0.000]]), 
                (2, [[-40.341,-72.000,942.588], [[9.659,0.000,2.588],[14.668,1.747,4.898],[3.866,3.938,3.419]], [[-0.194,1.299,0.724],[-0.404,1.248,0.764],[-1.155,0.929,0.236]], [[-0.000,0.000,-0.000],[0.002,0.030,-0.087],[-0.028,0.137,-0.267]], [[-0.336,-0.750,1.255],[-0.303,-0.837,1.207],[-0.346,-0.748,1.253]], [[-0.000,0.000,-0.000],[0.016,0.076,0.057],[-0.126,0.297,-0.280]]]), 
                (3, [[-27.303,-70.235,946.925], [11.405,1.782,3.774], [-0.428,1.262,0.699], [-0.003,0.005,-0.026], [-0.289,-0.787,1.244], [-0.022,0.019,0.007]]), 
                (4, [[-17.530,-68.503,950.142], [11.716,2.577,3.697], [-0.490,1.247,0.684], [-0.073,0.004,-0.076], [-0.226,-0.781,1.261], [0.090,0.070,0.060]]), 
                (5, [[-3.913,-64.921,954.252], [15.505,4.583,4.520], [-0.539,1.275,0.556], [0.058,0.073,-0.158], [-0.193,-0.663,1.334], [-0.120,0.137,0.055]]), 
                (6, [[-36.160,-68.538,945.729], [4.473,2.962,2.843], [-0.957,1.184,0.273], [0.234,0.180,0.173], [-0.344,-0.531,1.095], [0.083,0.089,-0.067]]), 
                (7, [[-31.481,-66.087,948.251], [5.027,2.166,2.464], [-0.784,1.227,0.522], [-0.298,0.768,0.549], [-0.258,-0.620,1.071], [-0.140,-0.514,0.642]]), 
                (8, [[-26.131,-64.238,950.634], [6.321,2.107,2.738], [-1.521,2.750,1.394], [-0.271,0.046,0.055], [-0.540,-1.527,2.422], [-0.008,-0.099,0.035]]), 
                (9, [[-18.830,-61.892,953.709], [7.186,2.087,2.614], [-0.758,1.503,0.882], [0.676,-0.698,-0.289], [-0.198,-0.789,1.175], [0.265,0.399,-0.548]]), 
                (10, [[-11.796,-60.057,955.883], [6.640,1.713,2.132], [-0.544,1.252,0.689], [-0.010,-0.150,-0.143], [-0.203,-0.783,1.262], [-0.106,-0.009,0.021]]), 
                (11, [[-5.550,-58.463,957.964], [5.851,1.475,2.030], [-0.548,1.265,0.662], [0.007,0.178,0.092], [-0.245,-0.768,1.265], [0.050,0.044,-0.005]])], [
                (1, [[-50.000,-72.000,940.000], [9.659,0.000,2.588], [-0.145,0.974,0.543], [0.000,0.000,0.000], [-0.252,-0.562,0.941], [0.000,0.000,0.000]]), 
                (2, [[-40.341,-72.000,942.588], [[9.659,0.000,2.588],[14.668,1.747,4.898],[3.866,3.938,3.419]], [[-0.145,0.974,0.543],[-0.303,0.936,0.573],[-0.866,0.697,0.177]], [[-0.000,0.000,-0.000],[0.002,0.022,-0.065],[-0.021,0.103,-0.200]], [[-0.252,-0.562,0.941],[-0.227,-0.628,0.905],[-0.259,-0.561,0.940]], [[-0.000,0.000,-0.000],[0.012,0.057,0.043],[-0.095,0.223,-0.210]]]), 
                (3, [[-27.303,-70.235,946.925], [11.405,1.782,3.774], [-0.321,0.946,0.524], [-0.003,0.004,-0.020], [-0.217,-0.590,0.933], [-0.017,0.014,0.005]]), 
                (4, [[-17.530,-68.503,950.142], [11.716,2.577,3.697], [-0.368,0.935,0.513], [-0.055,0.003,-0.057], [-0.170,-0.586,0.946], [0.067,0.053,0.045]]), 
                (5, [[-3.913,-64.921,954.252], [15.505,4.583,4.520], [-0.404,0.956,0.417], [0.044,0.055,-0.119], [-0.145,-0.497,1.000], [-0.090,0.102,0.041]]), 
                (6, [[-36.160,-68.538,945.729], [4.473,2.962,2.843], [-0.718,0.888,0.205], [0.175,0.135,0.130], [-0.258,-0.398,0.821], [0.062,0.067,-0.050]]), 
                (7, [[-31.481,-66.087,948.251], [5.027,2.166,2.464], [-0.588,0.920,0.391], [-0.223,0.576,0.412], [-0.193,-0.465,0.803], [-0.105,-0.386,0.481]]), 
                (8, [[-26.131,-64.238,950.634], [6.321,2.107,2.738], [-1.140,2.063,1.046], [-0.203,0.034,0.042], [-0.405,-1.145,1.817], [-0.006,-0.074,0.026]]), 
                (9, [[-18.830,-61.892,953.709], [7.186,2.087,2.614], [-0.568,1.127,0.662], [0.507,-0.523,-0.217], [-0.149,-0.592,0.881], [0.199,0.299,-0.411]]), 
                (10, [[-11.796,-60.057,955.883], [6.640,1.713,2.132], [-0.408,0.939,0.517], [-0.007,-0.112,-0.108], [-0.152,-0.587,0.947], [-0.080,-0.006,0.016]]), 
                (11, [[-5.550,-58.463,957.964], [5.851,1.475,2.030], [-0.411,0.949,0.496], [0.005,0.133,0.069], [-0.184,-0.576,0.949], [0.038,0.033,-0.003]])
                ]]),

            "userAnnotationGroups": [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-10',
                    'name': get_spinal_nerve_term('spinal nerve')[0],
                    'ontId': get_spinal_nerve_term('spinal nerve')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7-8',
                    'name': get_spinal_nerve_term('dorsal root ganglion')[0],
                    'ontId': get_spinal_nerve_term('dorsal root ganglion')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-10',
                    'name': get_spinal_nerve_term('dorsal root of spinal cord')[0],
                    'ontId': get_spinal_nerve_term('dorsal root of spinal cord')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '2-4',
                    'name': get_spinal_nerve_term('ventral root of spinal cord')[0],
                    'ontId': get_spinal_nerve_term('ventral root of spinal cord')[1]
                }
            ]
        })

class MeshType_3d_dorsalrootganglion1(Scaffold_base):
    """
    Generates a 3-D hermite dorsal root ganglion scaffold from a network layout description, with a solid core.
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
        return "3D Dorsal Root Ganglion 1"

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
