"""
Generates a 3-D Hermite bifurcating solid cylinder network.
"""

from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.utils.zinc.finiteelement import get_maximum_element_identifier, get_maximum_node_identifier
from cmlibs.zinc.field import Field
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
# from scaffoldmaker.utils.bifurcation_practice2 import generateTubeBifurcationTree
from scaffoldmaker.utils.bifurcation import generateTubeBifurcationTree


def findElementsCountAcross(cls, options):
    # Calculate elementsCountAcrosMinor
    options['Number of elements across minor'] = int(((options["Elements count around"] - 4) / 4 -
                                                      (options['Number of elements across major'] / 2)) * 2 + 6)

    return options['Number of elements across minor']

def findInitialElementsCountAround(cls, options, Rcrit):
    if options['Number of elements across transition'] <= Rcrit:
        options["Elements count around"] = (options["Elements count around"] + 8 *
                                            (options['Number of elements across transition'] - 2))
    else:
        options["Elements count around"] = options["Elements count around"] + 8 * (Rcrit - 1)

    return options["Elements count around"]

def setParametersToDefault(cls, options):
    options["Elements count around"] = 8
    options["Number of elements across major"] = 4
    options["Number of elements across minor"] = 4
    options["Number of elements across transition"] = 1

    return options

def checkCoreParameters(cls, options, dependentChanges=False):
    # Store initial parameter values
    Rcrit = max(min(options['Number of elements across major'] - 4,
                options['Number of elements across minor'] - 4) // 2 + 1, 0)
    elementsCountAcrossMajor = options['Number of elements across major']
    elementsCountAcrossMinor = options['Number of elements across minor']
    elementsCountAcrossTransition = options["Number of elements across transition"]
    elementsCountAround = options["Elements count around"] if options["Number of elements across transition"] == 1 \
        or options["Number of elements across transition"] > Rcrit \
        else findInitialElementsCountAround(cls, options, Rcrit)

    # Check elements count around are ok
    if options["Elements count around"] < 8:
        dependentChanges = True
        options = setParametersToDefault(cls, options)
    if options["Elements count around"] % 4:
        options["Elements count around"] += 4 - options["Elements count around"] % 4
        options['Number of elements across minor'] = findElementsCountAcross(cls, options)

    # Check elements count across major are ok
    if options['Number of elements across major'] < 4:
        options['Number of elements across major'] = 4
        options['Number of elements across minor'] = findElementsCountAcross(cls, options)
    if options['Number of elements across major'] % 2:
        options['Number of elements across major'] += 1
        options['Number of elements across minor'] = findElementsCountAcross(cls, options)

    # Calculate elementsCountAcrosMinor
    if elementsCountAcrossTransition > 1:
        options["Number of elements across minor"] = options["Number of elements across major"]
    else:
        options['Number of elements across minor'] = findElementsCountAcross(cls, options)
    if options['Number of elements across minor'] < 4:
        dependentChanges = True
        options = setParametersToDefault(cls, options)
    if options["Number of elements across minor"] % 2:
        options["Number of elements across minor"] += 1

    # Check elements count across transition
    if options['Number of elements across transition'] < 1:
        options['Number of elements across transition'] = 1
    if Rcrit >= options['Number of elements across transition'] > 1:
        dependentChanges = True
        if elementsCountAround > elementsCountAcrossMajor:
            options["Elements count around"] = (elementsCountAround - 8 *
                                                (options['Number of elements across transition'] - 1))
        else:
            options["Elements count around"] = (elementsCountAround + 8 *
                                                (options['Number of elements across transition'] - 1))

    # Rcrit check
    if elementsCountAcrossTransition > Rcrit:
        print("uh oh", Rcrit, options['Number of elements across transition'] - 1)
        dependentChanges = True
        # Set parameters to initial values
        options['Number of elements across transition'] = Rcrit
        options["Elements count around"] = elementsCountAround
        options['Number of elements across major'] = elementsCountAcrossMajor
        options['Number of elements across minor'] = elementsCountAcrossMinor

    # Number of elements around sanity check
    if options["Number of elements across major"] >= options["Number of elements across minor"]:
        if options["Number of elements across minor"] == 4:
            ec = (options["Number of elements across major"] / 2 - 1) * 4 + 4
        else:
            ec = (options["Number of elements across major"] / 2 + (
                    (options["Number of elements across minor"] - 6) / 2) -
                    (options["Number of elements across transition"] * 2 - 2)) * 4 + 4
    elif options["Number of elements across major"] < options["Number of elements across minor"]:
        if options["Number of elements across major"] == 4:
            ec = (options["Number of elements across minor"] / 2 - 1) * 4 + 4
        else:
            ec = (options["Number of elements across minor"] / 2 +
                  ((options["Number of elements across major"] - 6) / 2) -
                  (options["Number of elements across transition"] * 2 - 2)) * 4 + 4

    if options["Elements count around"] != ec:
        dependentChanges = True
        # Reset parameter values
        options['Number of elements across transition'] = 1
        options["Elements count around"] = elementsCountAround
        options['Number of elements across major'] = elementsCountAcrossMajor
        options['Number of elements across minor'] = findElementsCountAcross(cls, options)


    return options, dependentChanges


class MeshType_3d_solidcylindernetwork1(Scaffold_base):
    """
    Generates a 3-D hermite bifurcating tube network, with linear basis through wall.
    Number of elements generated are primarily controlled by "elements count around" and "elements count across major".
    The elements count around parameter determines the number of elements around a circular rim.
    The elements count across major determines the number of elements across the major axis of an ellipse (y-axis in the scaffold).
    The number of elements across minor axis (z-axis) is calculated internally based on the two elements count parameters.
    """

    @staticmethod
    def getName():
        return "3D Solid Cylinder Network 1"

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
            'Number of elements across major': 4,
            'Number of elements across minor': 4,
            "Elements count through wall": 1,
            'Number of elements across transition': 1,
            "Annotation elements counts around": [0],
            "Annotation elements counts across major": [0],
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
            'Number of elements across major',
            "Elements count through wall",
            'Number of elements across transition',
            "Annotation elements counts around",
            "Annotation elements counts across major",
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
            options["Network layout"] = cls.getOptionScaffoldPackage("Network layout")
        dependentChanges = False

        # Parameters specific to the core
        options, dependentChanges = checkCoreParameters(cls, options, dependentChanges)

        ## Add annotationElementsCountAround

        # When core option is false
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
        elementsCountAcrossMajor = options['Number of elements across major']
        elementsCountAcrossMinor = options['Number of elements across minor']  # This should be hidden
        elementsCountThroughWall = options["Elements count through wall"]
        elementsCountAcrossTransition = options['Number of elements across transition']  # Make this global parameter
        annotationElementsCountsAround = options["Annotation elements counts around"]
        annotationElementsCountsAcrossMajor = options["Annotation elements counts across major"]
        targetElementDensityAlongLongestSegment = options["Target element density along longest segment"]
        elementsCountAlong = int(targetElementDensityAlongLongestSegment)
        serendipity = options["Serendipity"]
        isCore = True

        # Create a list of parameters required for generating an ellipse
        ellipseParameters = [elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossTransition,
                             elementsCountAlong]

        layoutRegion = region.createRegion()
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()

        networkMesh = networkLayout.getConstructionObject()

        fieldmodule = region.getFieldmodule()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = find_or_create_field_coordinates(fieldmodule)
        nodeIdentifier = max(get_maximum_node_identifier(nodes), 0) + 1
        mesh = fieldmodule.findMeshByDimension(2)
        elementIdentifier = max(get_maximum_element_identifier(mesh), 0) + 1

        nodeIdentifier, elementIdentifier, annotationGroups = generateTubeBifurcationTree(
            networkMesh, region, coordinates, nodeIdentifier, elementIdentifier,
            elementsCountAround, targetElementDensityAlongLongestSegment, elementsCountThroughWall,
            layoutAnnotationGroups, annotationElementsCountsAround, annotationElementsCountsAcrossMajor,
            ellipseParameters,
            serendipity=serendipity, showTrimSurfaces=options["Show trim surfaces"], isCore=isCore)

        return annotationGroups, None
