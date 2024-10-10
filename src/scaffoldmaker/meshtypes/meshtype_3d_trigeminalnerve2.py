"""
Generates a 3D trigeminal nerve scaffold using tube network mesh.
"""

from cmlibs.zinc.element import Element
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.trigeminal_nerve_terms import get_trigeminal_nerve_term
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters
from cmlibs.zinc.node import Node


def getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName):
    assert parameterSetName in cls.getParameterSetNames()
    if parameterSetName in ("Default", "Human 1"):
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            "scaffoldSettings": {
                "Structure": "1-2-3,3-4.1,4.2-5,5-6,4.3-7,7-8,4.4-9,9-10-11",
                "Define inner coordinates": True
            },
            "meshEdits": exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[0.0000,0.0000,0.0000], [0.3570,0.0000,0.0000], [0.0000,0.1700,0.0000], [0.0000,-0.0132,0.0000], [0.0000,0.0000,0.0720], [0.0000,0.0000,0.0000]]), 
                (2, [[0.3530,0.0000,0.0000], [0.3490,0.0000,0.0000], [0.0000,0.1510,0.0000], [0.0000,-0.0248,0.0000], [0.0000,0.0000,0.0720], [0.0000,0.0000,-0.0000]]), 
                (3, [[0.6980,0.0000,0.0000], [0.3235,0.0000,0.0000], [0.0000,0.1204,0.0000], [0.0000,0.0715,0.0000], [0.0000,0.0000,0.0720], [0.0000,0.0000,-0.0000]]), 
                (4, [[1.0000,0.0000,0.0000], [[0.2805,0.0000,0.0000],[-0.0054,-0.3050,0.0000],[0.2812,-0.0004,0.0000],[0.0007,0.2753,0.0000]], [[0.0000,0.2940,0.0000],[0.3080,-0.0054,0.0000],[0.0003,0.2116,0.0000],[-0.1123,0.0003,0.0000]], [[-0.0000,0.2757,-0.0000],[-0.2114,0.0365,0.0000],[-0.0031,-0.0663,0.0000],[0.0235,0.0039,0.0000]], [[0.0000,0.0000,0.0720],[0.0000,0.0000,0.0720],[-0.0000,0.0000,0.0720],[0.0000,-0.0000,0.0720]], [[-0.0000,-0.0000,-0.0000],[0.0000,0.0000,0.0150],[0.0000,0.0000,0.0260],[0.0000,0.0000,0.0110]]]), 
                (5, [[1.0360,-0.2930,0.0000], [0.0774,-0.2768,0.0000], [0.1570,0.0439,0.0000], [-0.0893,0.0271,0.0000], [-0.0000,0.0000,0.0720], [0.0000,0.0000,-0.0150]]), 
                (6, [[1.1500,-0.5460,0.0000], [0.1481,-0.2255,0.0000], [0.1261,0.0828,0.0000], [0.0471,0.0282,-0.0000], [-0.0000,0.0000,0.0420], [-0.0000,-0.0000,-0.0450]]), 
                (7, [[1.2946,0.0044,0.0000], [0.3080,0.0092,0.0000], [-0.0041,0.1380,0.0000], [-0.0018,-0.0808,0.0000], [0.0000,-0.0000,0.0720], [0.0000,0.0000,-0.0260]]), 
                (8, [[1.6158,0.0189,0.0000], [0.3344,0.0197,0.0000], [-0.0030,0.0500,0.0000], [0.0084,-0.0949,-0.0000], [0.0000,-0.0000,0.0200], [-0.0000,-0.0000,-0.0780]]), 
                (9, [[1.0802,0.2794,0.0000], [0.1596,0.2691,0.0000], [-0.0747,0.0443,0.0000], [0.0533,0.0371,0.0000], [0.0000,-0.0000,0.0720], [0.0000,0.0000,-0.0110]]), 
                (10, [[1.3145,0.5118,0.0000], [0.2840,0.2385,0.0000], [-0.0397,0.0473,0.0000], [0.0163,-0.0171,0.0000], [0.0000,-0.0000,0.0500], [0.0000,0.0000,-0.0260]]), 
                (11, [[1.6488,0.7485,0.0000], [0.3833,0.2341,0.0000], [-0.0194,0.0317,0.0000], [0.0159,-0.0202,-0.0000], [0.0000,-0.0000,0.0200], [-0.0000,-0.0000,-0.0340]])], [
                (1, [[0.0000,0.0000,0.0000], [0.3570,0.0000,0.0000], [0.0000,0.1275,0.0000], [0.0000,-0.0099,0.0000], [0.0000,0.0000,0.0540], [0.0000,0.0000,0.0000]]), 
                (2, [[0.3530,0.0000,0.0000], [0.3490,0.0000,0.0000], [0.0000,0.1133,0.0000], [0.0000,-0.0186,0.0000], [0.0000,0.0000,0.0540], [0.0000,0.0000,-0.0000]]), 
                (3, [[0.6980,0.0000,0.0000], [0.3235,0.0000,0.0000], [0.0000,0.0903,0.0000], [0.0000,0.0536,0.0000], [0.0000,0.0000,0.0540], [0.0000,0.0000,-0.0000]]), 
                (4, [[1.0000,0.0000,0.0000], [[0.2805,0.0000,0.0000],[-0.0054,-0.3050,0.0000],[0.2812,-0.0004,0.0000],[0.0007,0.2753,0.0000]], [[0.0000,0.2205,0.0000],[0.2310,-0.0040,0.0000],[0.0002,0.1587,0.0000],[-0.0842,0.0002,0.0000]], [[-0.0000,0.2068,-0.0000],[-0.1586,0.0274,0.0000],[-0.0023,-0.0497,0.0000],[0.0176,0.0029,0.0000]], [[0.0000,0.0000,0.0540],[0.0000,0.0000,0.0540],[-0.0000,0.0000,0.0540],[0.0000,-0.0000,0.0540]], [[-0.0000,-0.0000,-0.0000],[0.0000,0.0000,0.0113],[0.0000,0.0000,0.0195],[0.0000,0.0000,0.0083]]]), 
                (5, [[1.0360,-0.2930,0.0000], [0.0774,-0.2768,0.0000], [0.1177,0.0329,0.0000], [-0.0670,0.0204,0.0000], [-0.0000,0.0000,0.0540], [0.0000,0.0000,-0.0113]]), 
                (6, [[1.1500,-0.5460,0.0000], [0.1481,-0.2255,0.0000], [0.0946,0.0621,0.0000], [0.0354,0.0211,-0.0000], [-0.0000,0.0000,0.0315], [-0.0000,-0.0000,-0.0338]]), 
                (7, [[1.2946,0.0044,0.0000], [0.3080,0.0092,0.0000], [-0.0031,0.1035,0.0000], [-0.0014,-0.0606,0.0000], [0.0000,-0.0000,0.0540], [0.0000,0.0000,-0.0195]]), 
                (8, [[1.6158,0.0189,0.0000], [0.3344,0.0197,0.0000], [-0.0022,0.0375,0.0000], [0.0063,-0.0711,-0.0000], [0.0000,-0.0000,0.0150], [-0.0000,-0.0000,-0.0585]]), 
                (9, [[1.0802,0.2794,0.0000], [0.1596,0.2691,0.0000], [-0.0560,0.0332,0.0000], [0.0400,0.0278,0.0000], [0.0000,-0.0000,0.0540], [0.0000,0.0000,-0.0083]]), 
                (10, [[1.3145,0.5118,0.0000], [0.2840,0.2385,0.0000], [-0.0298,0.0355,0.0000], [0.0122,-0.0128,0.0000], [0.0000,-0.0000,0.0375], [0.0000,0.0000,-0.0195]]), 
                (11, [[1.6488,0.7485,0.0000], [0.3833,0.2341,0.0000], [-0.0146,0.0238,0.0000], [0.0119,-0.0152,-0.0000], [0.0000,-0.0000,0.0150], [-0.0000,-0.0000,-0.0255]])
                ]]),

            "userAnnotationGroups": [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-10',
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
                    'identifierRanges': '4-5',
                    'name': get_trigeminal_nerve_term('mandibular nerve')[0],
                    'ontId': get_trigeminal_nerve_term('mandibular nerve')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-7',
                    'name': get_trigeminal_nerve_term('maxillary nerve')[0],
                    'ontId': get_trigeminal_nerve_term('maxillary nerve')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '8-10',
                    'name': get_trigeminal_nerve_term('ophthalmic nerve')[0],
                    'ontId': get_trigeminal_nerve_term('ophthalmic nerve')[1]
                }
            ]
        })
    

class MeshType_3d_trigeminalnerve2(Scaffold_base):
    """
    Generates a 3-D hermite bifurcating tube network with core representing the trigeminal nerve.
    """

    @classmethod
    def getName(cls):
        return "3D Trigeminal Nerve 2"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1"            
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):

        options = {
            'Base parameter set': parameterSetName,
            "Network layout": getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName),  
            "Number of elements around trigeminal nerve root": 12,
            "Number of elements around mandibular nerve": 8,
            "Number of elements around maxillary nerve": 12,
            "Number of elements around ophthalmic nerve": 8,
            "Number of elements through shell": 1,
            "Target element density along longest segment": 5.0,
            "Show trim surfaces": False,
            "Use Core": True,
            "Number of elements across core box minor": 2,
            "Number of elements across core transition": 1
        }
        
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            "Network layout",
            "Number of elements around trigeminal nerve root",
            "Number of elements around mandibular nerve",
            "Number of elements around maxillary nerve",
            "Number of elements around ophthalmic nerve",
            "Number of elements through shell",
            "Target element density along longest segment",
            "Show trim surfaces",
            "Use Core",
            "Number of elements across core box minor",
            "Number of elements across core transition"]
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Network layout":
            return [MeshType_1d_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == "Network layout":
            return cls.getParameterSetNames()
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + ".getOptionScaffoldTypeParameterSetNames.  " + \
            "Invalid option \"" + optionName + "\" scaffold type " + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

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
        dependentChanges = False
        if not options["Network layout"].getScaffoldType() in cls.getOptionValidScaffoldTypes("Network layout"):
            options["Network layout"] = cls.getOptionScaffoldPackage('Network layout', MeshType_1d_network_layout1)
        minElementsCountAround = None
        
        for keys in ["Number of elements around trigeminal nerve root",
                     "Number of elements around mandibular nerve",
                     "Number of elements around maxillary nerve",
                     "Number of elements around ophthalmic nerve"]:
            if options[keys] % 4:
                options[keys] += 4 - (options[keys] % 4)
            if (minElementsCountAround is None) or (options[keys] < minElementsCountAround):
                minElementsCountAround = options[keys]

        for keys in ["Number of elements around trigeminal nerve root",
                     "Number of elements around maxillary nerve"]:
            if options[keys] < 12:
                options[keys] = 12

        for keys in ["Number of elements around mandibular nerve",
                     "Number of elements around ophthalmic nerve"]:
            if options[keys] < 8:
                options[keys] = 8
            if options[keys] == options["Number of elements around trigeminal nerve root"]:
                options[keys] = options["Number of elements around trigeminal nerve root"] - 4

        if options["Number of elements through shell"] < 0:
            options["Number of elements through shell"] = 1

        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0

        if options["Number of elements across core transition"] < 1:
            options["Number of elements across core transition"] = 1

        maxElementsCountCoreBoxMinor = minElementsCountAround // 2 - 2
        for key in [
            "Number of elements across core box minor"
        ]:
            if options[key] < 2:
                options[key] = 2
            elif options[key] > maxElementsCountCoreBoxMinor:
                options[key] = maxElementsCountCoreBoxMinor
                dependentChanges = True
            elif options[key] % 2:
                options[key] += options[key] % 2

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        # parameterSetName = options['Base parameter set']
        # isHuman = parameterSetName in ["Default", "Human 1", "Human 2 Coarse", "Human 2 Medium", "Human 2 Fine"]

        layoutRegion = region.createRegion()
        networkLayout = options["Network layout"]
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        coreBoxMinorCount = options["Number of elements across core box minor"]
        coreTransitionCount = options['Number of elements across core transition']
        annotationAroundCounts = []

        # implementation currently uses major count including transition
        annotationCoreMajorCounts = []
        for layoutAnnotationGroup in layoutAnnotationGroups:
            aroundCount = 0
            coreMajorCount = 0
            name = layoutAnnotationGroup.getName()
            if name in ["trigeminal nerve root", "mandibular nerve", "maxillary nerve", "ophthalmic nerve"]:
                aroundCount = options["Number of elements around " + name]
                coreMajorCount = aroundCount // 2 - coreBoxMinorCount + 2 * coreTransitionCount
            annotationAroundCounts.append(aroundCount)
            annotationCoreMajorCounts.append(coreMajorCount)

        isCore = options["Use Core"]

        tubeNetworkMeshBuilder = TubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            defaultElementsCountAround=options["Number of elements around trigeminal nerve root"],
            elementsCountThroughWall=options["Number of elements through shell"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            annotationElementsCountsAround=annotationAroundCounts,
            defaultElementsCountAcrossMajor=annotationCoreMajorCounts[-1],
            elementsCountTransition=coreTransitionCount,
            annotationElementsCountsAcrossMajor=annotationCoreMajorCounts,
            isCore=isCore)

        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughWall=False,
            isShowTrimSurfaces=options["Show trim surfaces"])
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        return annotationGroups, None

