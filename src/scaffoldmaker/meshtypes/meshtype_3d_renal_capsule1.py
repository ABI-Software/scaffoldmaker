"""
Generates a 3D renal capsule using tube network mesh.
"""
import math

from cmlibs.maths.vectorops import mult, set_magnitude, cross
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.zinc.field import Field

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm, findAnnotationGroupByName
from scaffoldmaker.annotation.kidney_terms import get_kidney_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine, sampleCubicHermiteCurves, \
    smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.tubenetworkmesh import TubeNetworkMeshBuilder, TubeNetworkMeshGenerateData, \
    KidneyTubeNetworkMeshBuilder
from cmlibs.zinc.node import Node


class MeshType_1d_renal_capsule_network_layout1(MeshType_1d_network_layout1):
    """
    Defines renal capsule network layout.
    """

    @classmethod
    def getName(cls):
        return "1D Renal Capsule Network Layout 1"

    @classmethod
    def getParameterSetNames(cls):
        return ["Default"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Base parameter set"] = "Human 1" if (parameterSetName == "Default") else parameterSetName
        options["Define inner coordinates"] = True
        options["Elements count along"] = 2
        options["Renal capsule length"] = 1.0
        options["Renal capsule diameter"] = 1.5
        options["Renal capsule bend angle degrees"] = 10
        options["Inner proportion default"] = 0.6
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Elements count along",
            "Renal capsule length",
            "Renal capsule diameter",
            "Renal capsule bend angle degrees",
            "Inner proportion default"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        for key in [
            "Renal capsule length",
            "Renal capsule diameter"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1

        if options["Elements count along"] < 2:
            options["Elements count along"] = 2

        if options["Renal capsule bend angle degrees"] < 0.0:
            options["Renal capsule bend angle degrees"] = 0.0
        elif options["Renal capsule bend angle degrees"] > 30.0:
            options["Renal capsule bend angle degrees"] = 30.0

        if options["Inner proportion default"] < 0.1:
            options["Inner proportion default"] = 0.1
        elif options["Inner proportion default"] > 0.9:
            options["Inner proportion default"] = 0.9

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the unrefined mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return [] empty list of AnnotationGroup, NetworkMesh
        """
        # parameters
        structure = options["Structure"] = cls.getLayoutStructure(options)
        capsuleElementsCount = options["Elements count along"]
        capsuleLength = options["Renal capsule length"]
        capsuleRadius = 0.5 * options["Renal capsule diameter"]
        capsuleBendAngle = options["Renal capsule bend angle degrees"]
        innerProportionDefault = options["Inner proportion default"]

        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(1)

        # set up element annotations
        kidneyGroup = AnnotationGroup(region, get_kidney_term("kidney"))
        annotationGroups = [kidneyGroup]

        kidneyMeshGroup = kidneyGroup.getMeshGroup(mesh)
        elementIdentifier = 1
        meshGroups = [kidneyMeshGroup]
        for e in range(capsuleElementsCount):
            element = mesh.findElementByIdentifier(elementIdentifier)
            for meshGroup in meshGroups:
                meshGroup.addElement(element)
            elementIdentifier += 1

        # set coordinates (outer)
        fieldcache = fieldmodule.createFieldcache()
        coordinates = find_or_create_field_coordinates(fieldmodule)
        # need to ensure inner coordinates are at least defined:
        cls.defineInnerCoordinates(region, coordinates, options, networkMesh, innerProportion=innerProportionDefault)
        innerCoordinates = find_or_create_field_coordinates(fieldmodule, "inner coordinates")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        # renal capsule
        nodeIdentifier = 1
        halfCapsuleLength = 0.5 * capsuleLength
        capsuleScale = capsuleLength / capsuleElementsCount
        bendAngleRadians = math.radians(capsuleBendAngle)
        sinBendAngle = math.sin(bendAngleRadians)
        cosBendAngle = math.cos(bendAngleRadians)
        sinCurveAngle = math.sin(3 * bendAngleRadians)
        mx = [0.0, 0.0, 0.0]
        d1 = [capsuleScale, 0.0, 0.0]
        d3 = [0.0, 0.0, capsuleRadius]
        id3 = mult(d3, innerProportionDefault)

        tx = halfCapsuleLength * -cosBendAngle
        ty = halfCapsuleLength * -sinBendAngle
        sx = [tx, ty, 0.0]
        ex = [-tx, ty, 0.0]
        sd1 = mult([1.0, sinCurveAngle, 0.0], capsuleScale)
        ed1 = [sd1[0], -sd1[1], sd1[2]]
        nx, nd1 = sampleCubicHermiteCurves([sx, mx, ex], [sd1, d1, ed1], capsuleElementsCount)[0:2]
        nd1 = smoothCubicHermiteDerivativesLine(nx, nd1)

        sd2_list = []
        sd3_list = []
        sNodeIdentifiers = []
        for e in range(capsuleElementsCount + 1):
            sNodeIdentifiers.append(nodeIdentifier)
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            fieldcache.setNode(node)
            sd2 = set_magnitude(cross(d3, nd1[e]), capsuleRadius)
            sid2 = mult(sd2, innerProportionDefault)
            sd2_list.append(sd2)
            sd3_list.append(d3)
            for field, derivatives in ((coordinates, (nd1[e], sd2, d3)), (innerCoordinates, (nd1[e], sid2, id3))):
                setNodeFieldParameters(field, fieldcache, nx[e], *derivatives)
            nodeIdentifier += 1

        sd12 = smoothCurveSideCrossDerivatives(nx, nd1, [sd2_list])[0]
        sd13 = smoothCurveSideCrossDerivatives(nx, nd1, [sd3_list])[0]
        for e in range(capsuleElementsCount + 1):
            node = nodes.findNodeByIdentifier(sNodeIdentifiers[e])
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sd12[e])
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sd13[e])
            sid12 = mult(sd12[e], innerProportionDefault)
            sid13 = mult(sd13[e], innerProportionDefault)
            innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sid12)
            innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sid13)

        return annotationGroups, networkMesh

    @classmethod
    def getLayoutStructure(cls, options):
        """
        Generate 1D layout structure based on the number of elements count along.
        :param options: Dict containing options. See getDefaultOptions().
        :return string version of the 1D layout structure
        """
        nodesCountAlong = options["Elements count along"] + 1
        assert nodesCountAlong > 1

        return f"({'-'.join(str(i) for i in range(1, nodesCountAlong + 1))})"


class MeshType_3d_renal_capsule1(Scaffold_base):
    """
    Generates a 3-D renal capsule.
    """

    @classmethod
    def getName(cls):
        return "3D Renal Capsule 1"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {}
        useParameterSetName = "Human 1" if (parameterSetName == "Default") else parameterSetName
        options["Base parameter set"] = useParameterSetName
        options["Renal capsule network layout"] = ScaffoldPackage(MeshType_1d_renal_capsule_network_layout1)
        options["Elements count around"] = 8
        options["Elements count through shell"] = 1
        options["Annotation elements counts around"] = [0]
        options["Target element density along longest segment"] = 2.0
        options["Number of elements across core box minor"] = 2
        options["Number of elements across core transition"] = 1
        options["Annotation numbers of elements across core box minor"] = [0]
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            "Renal capsule network layout",
            "Elements count around",
            "Elements count through shell",
            "Annotation elements counts around",
            "Target element density along longest segment",
            "Number of elements across core box minor",
            "Number of elements across core transition",
            "Annotation numbers of elements across core box minor"]
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Renal capsule network layout":
            return [MeshType_1d_renal_capsule_network_layout1]
        return []


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
        if optionName == "Renal capsule network layout":
            if not parameterSetName:
                parameterSetName = "Default"
            return ScaffoldPackage(MeshType_1d_renal_capsule_network_layout1, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if (options["Renal capsule network layout"].getScaffoldType() not in
                cls.getOptionValidScaffoldTypes("Renal capsule network layout")):
            options["Renal capsule network layout"] = ScaffoldPackage(MeshType_1d_renal_capsule_network_layout1)

        if options["Elements count around"] < 8:
            options["Elements count around"] = 8
        elif options["Elements count around"] % 4:
            options["Elements count around"] += 4 - (options["Elements count around"] % 4)

        if options["Elements count through shell"] < 1:
            options["Elements count through shell"] = 1

        if options["Number of elements across core transition"] < 1:
            options["Number of elements across core transition"] = 1

        minElementsCountAround = options["Elements count around"]
        maxElementsCountCoreBoxMinor = minElementsCountAround // 2 - 2
        if options["Number of elements across core box minor"] < 2:
            options["Number of elements across core box minor"] = 2
        elif options["Number of elements across core box minor"] > maxElementsCountCoreBoxMinor:
            options["Number of elements across core box minor"] = maxElementsCountCoreBoxMinor
            dependentChanges = True
        elif options["Number of elements across core box minor"] % 2:
            options["Number of elements across core box minor"] += options["Number of elements across core box minor"] % 2

        annotationElementsCountsAround = options["Annotation elements counts around"]
        if len(annotationElementsCountsAround) == 0:
            options["Annotation elements count around"] = [0]
        else:
            for i in range(len(annotationElementsCountsAround)):
                if annotationElementsCountsAround[i] <= 0:
                    annotationElementsCountsAround[i] = 0
                else:
                    if annotationElementsCountsAround[i] < 8:
                        annotationElementsCountsAround[i] = 8
                    elif annotationElementsCountsAround[i] % 4:
                        annotationElementsCountsAround[i] += 4 - (annotationElementsCountsAround[i] % 4)
                    if annotationElementsCountsAround[i] < minElementsCountAround:
                        minElementsCountAround = annotationElementsCountsAround[i]

        annotationCoreBoxMinorCounts = options["Annotation numbers of elements across core box minor"]
        if len(annotationCoreBoxMinorCounts) == 0:
            annotationCoreBoxMinorCounts = options["Annotation numbers of elements across core box minor"] = [0]
        if len(annotationCoreBoxMinorCounts) > len(annotationElementsCountsAround):
            annotationCoreBoxMinorCounts = options["Annotation numbers of elements across core box minor"] = \
                annotationCoreBoxMinorCounts[:len(annotationElementsCountsAround)]
            dependentChanges = True
        for i in range(len(annotationCoreBoxMinorCounts)):
            aroundCount = annotationElementsCountsAround[i] if annotationElementsCountsAround[i] \
                else options["Elements count around"]
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

        if options["Target element density along longest segment"] < 2.0:
            options["Target element density along longest segment"] = 2.0
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        networkLayout = options["Renal capsule network layout"]
        layoutRegion = region.createRegion()
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()

        tubeNetworkMeshBuilder = KidneyTubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            defaultElementsCountAround=options["Elements count around"],
            elementsCountThroughShell=options["Elements count through shell"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            isCore=True,
            elementsCountTransition=options["Number of elements across core transition"],
            defaultElementsCountCoreBoxMinor=options["Number of elements across core box minor"],
            annotationElementsCountsCoreBoxMinor=options["Annotation numbers of elements across core box minor"]
        )

        tubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughShell=False)
        tubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        # add kidney-specific annotation groups
        fm = region.getFieldmodule()
        mesh = generateData.getMesh()

        coreGroup = getAnnotationGroupForTerm(annotationGroups, ("core", "")).getGroup()
        shellGroup = getAnnotationGroupForTerm(annotationGroups, ("shell", "")).getGroup()
        openingGroup = getAnnotationGroupForTerm(annotationGroups, ("opening", "")).getGroup()

        hilumGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("hilum of kidney"))
        hilumGroup.getMeshGroup(mesh).addElementsConditional(openingGroup)

        cortexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("cortex of kidney"))
        tempGroup = fm.createFieldSubtract(shellGroup, openingGroup)
        cortexGroup.getMeshGroup(mesh).addElementsConditional(tempGroup)

        medullaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("renal medulla"))
        medullaGroup.getMeshGroup(mesh).addElementsConditional(coreGroup)

        for term in ["core", "shell", "opening"]:
            annotationGroups.remove(findAnnotationGroupByName(annotationGroups, term))

        return annotationGroups, None


    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        fm = region.getFieldmodule()
        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)

        is_exterior = fm.createFieldIsExterior()
        kidneyGroup = getAnnotationGroupForTerm(annotationGroups, get_kidney_term("kidney")).getGroup()

        # surface groups
        kidneyCapsuleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("kidney capsule"))
        kidneyCapsuleGroup.getMeshGroup(mesh2d).addElementsConditional(fm.createFieldAnd(kidneyGroup, is_exterior))

        surfaceTypes = [
            "anterior", "posterior",
            "lateral", "medial",
            "dorsal", "ventral"
        ]

        arb_group = {}
        kidney_exterior = {}
        for surfaceType in surfaceTypes:
            group = getAnnotationGroupForTerm(annotationGroups, (surfaceType, ""))
            group2d = group.getGroup()
            group2d_exterior = fm.createFieldAnd(group2d, is_exterior)

            term = f"{surfaceType} surface of kidney"
            surfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(term))
            surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(group2d_exterior)

            arb_group.update({surfaceType: group2d})
            kidney_exterior.update({term: group2d_exterior})

        # edge groups
        dorsalVentralBorderGroup = fm.createFieldAnd(fm.createFieldAnd(arb_group["dorsal"], arb_group["ventral"]), is_exterior)

        for surfaceType in ["lateral", "medial"]:
            term = f"{surfaceType} edge of kidney"
            edgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(term))
            edgeGroup.getMeshGroup(mesh1d).addElementsConditional(fm.createFieldAnd(arb_group[surfaceType], dorsalVentralBorderGroup))

        # for term in ["lateral", "medial", "dorsal", "ventral"]:
        #     annotationGroups.remove(findAnnotationGroupByName(annotationGroups, term))


def setNodeFieldParameters(field, fieldcache, x, d1, d2, d3, d12=None, d13=None):
    """
    Assign node field parameters x, d1, d2, d3 of field.
    :param field: Field parameters to assign.
    :param fieldcache: Fieldcache with node set.
    :param x: Parameters to set for Node.VALUE_LABEL_VALUE.
    :param d1: Parameters to set for Node.VALUE_LABEL_D_DS1.
    :param d2: Parameters to set for Node.VALUE_LABEL_D_DS2.
    :param d3: Parameters to set for Node.VALUE_LABEL_D_DS3.
    :param d12: Optional parameters to set for Node.VALUE_LABEL_D2_DS1DS2.
    :param d13: Optional parameters to set for Node.VALUE_LABEL_D2_DS1DS3.
    :return:
    """
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
    if d12:
        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d12)
    if d13:
        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, d13)


def setNodeFieldVersionDerivatives(field, fieldcache, version, d1, d2, d3, d12=None, d13=None):
    """
    Assign node field parameters d1, d2, d3 of field.
    :param field: Field to assign parameters of.
    :param fieldcache: Fieldcache with node set.
    :param version: Version of d1, d2, d3 >= 1.
    :param d1: Parameters to set for Node.VALUE_LABEL_D_DS1.
    :param d2: Parameters to set for Node.VALUE_LABEL_D_DS2.
    :param d3: Parameters to set for Node.VALUE_LABEL_D_DS3.
    :param d12: Optional parameters to set for Node.VALUE_LABEL_D2_DS1DS2.
    :param d13: Optional parameters to set for Node.VALUE_LABEL_D2_DS1DS3.
    :return:
    """
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, version, d1)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, version, d2)
    field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, version, d3)
    if d12:
        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, version, d12)
    if d13:
        field.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, version, d13)


def setHilumGroupThreshold(fm, coordinates, halfLength):
    """
    Creates a field to identify lung base elements based on y-coordinate threshold.
    Elements with y-coordinates below 45% of the rotated half-breadth are considered part of the lung base region for
    annotation purposes.
    :param fm: Field module used for creating and managing fields.
    :param coordinates: The coordinate field.
    :param halfLength: Half-length of tube.
    :return is_within_threshold: True for elements between the positive and negative x-threshold.
    """
    x_component = fm.createFieldComponent(coordinates, [1])
    x_threshold = 0.25 * halfLength
    minus_one = fm.createFieldConstant(-1)

    x_threshold_field = fm.createFieldConstant(x_threshold)
    is_less_than_threshold = fm.createFieldLessThan(x_component, x_threshold_field)
    is_greater_than_threshold = fm.createFieldGreaterThan(x_component, x_threshold_field * minus_one)
    is_within_threshold = fm.createFieldAnd(is_less_than_threshold, is_greater_than_threshold)

    return is_within_threshold
