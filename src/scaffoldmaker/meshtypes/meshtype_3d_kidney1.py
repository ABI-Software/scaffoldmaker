"""
Generates a 3D kidney using tube network mesh.
"""
import math

from cmlibs.maths.vectorops import mult, set_magnitude, cross, rotate_about_z_axis
from cmlibs.utils.zinc.field import find_or_create_field_coordinates, findOrCreateFieldCoordinates
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.field import Field

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm, findAnnotationGroupByName
from scaffoldmaker.annotation.kidney_terms import get_kidney_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine, sampleCubicHermiteCurves, \
    smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.networkmesh import NetworkMesh
from scaffoldmaker.utils.tubenetworkmesh import KidneyTubeNetworkMeshBuilder, TubeNetworkMeshGenerateData
from cmlibs.zinc.node import Node

from scaffoldmaker.utils.zinc_utils import translate_nodeset_coordinates


class MeshType_1d_kidney_network_layout1(MeshType_1d_network_layout1):
    """
    Defines kidney network layout.
    """

    showKidneys = [False, False]

    @classmethod
    def getName(cls):
        return "1D Kidney Network Layout 1"

    @classmethod
    def getParameterSetNames(cls):
        return ["Default"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        options["Base parameter set"] = "Human 1" if (parameterSetName == "Default") else parameterSetName
        options["Define inner coordinates"] = True
        options["Left kidney"] = True
        options["Right kidney"] = True
        options["Elements count along"] = 2
        options["Kidney length"] = 1.0
        options["Kidney width"] = 0.5
        options["Kidney thickness"] = 0.4
        options["Inner proportion default"] = 0.6
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Left kidney",
            "Right kidney",
            "Elements count along",
            "Kidney length",
            "Kidney width",
            "Kidney thickness",
            "Inner proportion default"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        for key in [
            "Kidney length",
            "Kidney width",
            "Kidney thickness"
        ]:
            if options[key] < 0.1:
                options[key] = 0.1

        if options["Elements count along"] < 2:
            options["Elements count along"] = 2

        if options["Inner proportion default"] < 0.1:
            options["Inner proportion default"] = 0.1
        elif options["Inner proportion default"] > 0.9:
            options["Inner proportion default"] = 0.9

        if not options["Left kidney"] and not options["Right kidney"]:
            dependentChanges = True
            options["Left kidney"] = True

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the unrefined mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: [] empty list of AnnotationGroup, NetworkMesh
        """
        # parameters
        structure = options["Structure"] = cls.getLayoutStructure(options)
        kidneyElementsCount = options["Elements count along"]
        isLeftKidney = options["Left kidney"]
        isRightKidney = options["Right kidney"]
        kidneyLength = options["Kidney length"]
        halfKidneyLength = 0.5 * kidneyLength
        halfKidneyWidth = 0.5 * options["Kidney width"]
        halfKidneyThickness = 0.5 * options["Kidney thickness"]
        innerProportionDefault = options["Inner proportion default"]
        cls.setShowKidneys(options)

        networkMesh = NetworkMesh(structure)
        networkMesh.create1DLayoutMesh(region)

        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(1)

        # set up element annotations
        kidneyGroup = AnnotationGroup(region, get_kidney_term("kidney"))
        kidneyMeshGroup = kidneyGroup.getMeshGroup(mesh)

        leftKidneyGroup = AnnotationGroup(region, get_kidney_term("left kidney"))
        leftKidneyMeshGroup = leftKidneyGroup.getMeshGroup(mesh)

        rightKidneyGroup = AnnotationGroup(region, get_kidney_term("right kidney"))
        rightKidneyMeshGroup = rightKidneyGroup.getMeshGroup(mesh)

        annotationGroups = [kidneyGroup, leftKidneyGroup, rightKidneyGroup]
        meshGroups = [kidneyMeshGroup, leftKidneyMeshGroup, rightKidneyMeshGroup]

        # set coordinates (outer)
        fieldcache = fieldmodule.createFieldcache()
        coordinates = find_or_create_field_coordinates(fieldmodule)
        # need to ensure inner coordinates are at least defined:
        cls.defineInnerCoordinates(region, coordinates, options, networkMesh, innerProportion=innerProportionDefault)
        innerCoordinates = find_or_create_field_coordinates(fieldmodule, "inner coordinates")
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        # Kidney
        nodeIdentifier = 1
        elementIdentifier = 1
        capRadius = cls.getCapRadius(halfKidneyWidth, halfKidneyThickness) * (
                halfKidneyWidth * 0.45 + halfKidneyThickness * 0.55)
        extensionLength = 0.5 * (halfKidneyWidth * 0.45 + halfKidneyThickness * 0.55)
        halfLayoutLength = (halfKidneyLength - capRadius - extensionLength)
        kidneyScale = 2 * halfLayoutLength / kidneyElementsCount

        leftKidney, rightKidney = 0, 1
        kidneys = [kidney for show, kidney in [(isLeftKidney, leftKidney), (isRightKidney, rightKidney)] if show]
        for kidney in kidneys:
            mx = [0.0, 0.0, 0.0]
            d1 = [kidneyScale, 0.0, 0.0]
            d3 = [0.0, 0.0, halfKidneyThickness]
            id3 = mult(d3, innerProportionDefault)

            tx = halfLayoutLength
            sx = [tx, 0.0, 0.0] if kidney is leftKidney else [tx, 0.0, 0.0]
            ex = [-tx, 0.0, 0.0] if kidney is leftKidney else [-tx, 0.0, 0.0]
            sd1 = mult([1.0, 0.0, 0.0], kidneyScale)
            ed1 = [sd1[0], -sd1[1], sd1[2]]
            nx, nd1 = sampleCubicHermiteCurves([sx, mx, ex], [sd1, d1, ed1], kidneyElementsCount)[0:2]
            nd1 = smoothCubicHermiteDerivativesLine(nx, nd1)

            sd2_list = []
            sd3_list = []
            sNodeIdentifiers = []
            for e in range(kidneyElementsCount + 1):
                sNodeIdentifiers.append(nodeIdentifier)
                node = nodes.findNodeByIdentifier(nodeIdentifier)
                fieldcache.setNode(node)
                sd2 = set_magnitude(cross(d3, nd1[e]), halfKidneyWidth)
                sid2 = mult(sd2, innerProportionDefault)
                sd2_list.append(sd2)
                sd3_list.append(d3)
                for field, derivatives in ((coordinates, (nd1[e], sd2, d3)), (innerCoordinates, (nd1[e], sid2, id3))):
                    setNodeFieldParameters(field, fieldcache, nx[e], *derivatives)
                nodeIdentifier += 1

            sd12 = smoothCurveSideCrossDerivatives(nx, nd1, [sd2_list])[0]
            sd13 = smoothCurveSideCrossDerivatives(nx, nd1, [sd3_list])[0]
            for e in range(kidneyElementsCount + 1):
                node = nodes.findNodeByIdentifier(sNodeIdentifiers[e])
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sd12[e])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sd13[e])
                sid12 = mult(sd12[e], innerProportionDefault)
                sid13 = mult(sd13[e], innerProportionDefault)
                innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, sid12)
                innerCoordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, sid13)

            # add annotations
            for e in range(kidneyElementsCount):
                element = mesh.findElementByIdentifier(elementIdentifier)
                meshGroups[0].addElement(element)
                if kidney is leftKidney and isLeftKidney:
                    meshGroups[1].addElement(element)
                if kidney is rightKidney and isRightKidney:
                    meshGroups[2].addElement(element)
                elementIdentifier += 1

        return annotationGroups, networkMesh

    @classmethod
    def getLayoutStructure(cls, options):
        """
        Generate 1D layout structure based on the number of elements count along.
        :param options: Dict containing options. See getDefaultOptions().
        :return: string version of the 1D layout structure
        """
        nodes_count = options["Elements count along"] + 1
        assert nodes_count > 1

        left = f"({'-'.join(map(str, range(1, nodes_count + 1)))})"
        if options["Left kidney"] and options["Right kidney"]:
            right = f"({'-'.join(map(str, range(nodes_count + 1, 2 * nodes_count + 1)))})"
            return f"{left},{right}"
        return left

    @classmethod
    def getCapRadius(cls, majorRadius, minorRadius):
        """
        Calculate the radius of the cap mesh based on the major radius and the minor radius of tube cross-section.
        :param majorRadius: The radius of a tube in the major-axis.
        :param minorRadius: The radius of a tube in the minor-axis.
        :return: Cap radius
        """
        if majorRadius > minorRadius:
            return math.pow((majorRadius / minorRadius), 1 / 3)
        elif majorRadius < minorRadius:
            return math.pow((minorRadius / majorRadius), 1 / 3)
        else:
            return majorRadius

    @classmethod
    def getShowKidneys(cls):
        return cls.showKidneys

    @classmethod
    def setShowKidneys(cls, options):
        cls.showKidneys[0] = True if options["Left kidney"] else False
        cls.showKidneys[1] = True if options["Right kidney"] else False


class MeshType_3d_kidney1(Scaffold_base):
    """
    Generates a 3-D Kidney.
    """

    @classmethod
    def getName(cls):
        return "3D Kidney 1"

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
        options["Kidney network layout"] = ScaffoldPackage(MeshType_1d_kidney_network_layout1)
        options["Number of elements around"] = 8
        options["Number of elements through shell"] = 1
        options["Annotation numbers of elements around"] = [0]
        options["Target element density along longest segment"] = 2.0
        options["Number of elements across core box minor"] = 2
        options["Number of elements across core transition"] = 1
        options["Annotation numbers of elements across core box minor"] = [0]
        options["Kidney spacing"] = 1.0
        options["Kidney curvature"] = 1.0
        options["Refine"] = False
        options["Refine number of elements"] = 4
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        optionNames = [
            "Kidney network layout",
            "Number of elements around",
            "Number of elements through shell",
            "Annotation numbers of elements around",
            "Target element density along longest segment",
            "Number of elements across core box minor",
            "Number of elements across core transition",
            "Annotation numbers of elements across core box minor",
            "Kidney spacing",
            "Kidney curvature",
            "Refine",
            "Refine number of elements"
        ]
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == "Kidney network layout":
            return [MeshType_1d_kidney_network_layout1]
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
        if optionName == "Kidney network layout":
            if not parameterSetName:
                parameterSetName = "Default"
            return ScaffoldPackage(MeshType_1d_kidney_network_layout1, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + ".getOptionScaffoldPackage:  Option " + optionName + " is not a scaffold"

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if (options["Kidney network layout"].getScaffoldType() not in
                cls.getOptionValidScaffoldTypes("Kidney network layout")):
            options["Kidney network layout"] = ScaffoldPackage(MeshType_1d_kidney_network_layout1)

        if options["Number of elements around"] < 8:
            options["Number of elements around"] = 8
        elif options["Number of elements around"] % 4:
            options["Number of elements around"] += 4 - (options["Number of elements around"] % 4)

        if options["Number of elements through shell"] < 1:
            options["Number of elements through shell"] = 1

        if options["Number of elements across core transition"] < 1:
            options["Number of elements across core transition"] = 1

        minElementsCountAround = options["Number of elements around"]
        maxElementsCountCoreBoxMinor = minElementsCountAround // 2 - 2
        if options["Number of elements across core box minor"] < 2:
            options["Number of elements across core box minor"] = 2
        elif options["Number of elements across core box minor"] > maxElementsCountCoreBoxMinor:
            options["Number of elements across core box minor"] = maxElementsCountCoreBoxMinor
            dependentChanges = True
        elif options["Number of elements across core box minor"] % 2:
            options["Number of elements across core box minor"] += options["Number of elements across core box minor"] % 2

        annotationElementsCountsAround = options["Annotation numbers of elements around"]
        if len(annotationElementsCountsAround) == 0:
            options["Annotation numbers of elements around"] = [0]
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

        if options["Target element density along longest segment"] < 2.0:
            options["Target element density along longest segment"] = 2.0

        if options['Refine number of elements'] < 1:
            options['Refine number of elements'] = 1
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base hermite-bilinear mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        networkLayout = options["Kidney network layout"]
        layoutRegion = region.createRegion()
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        networkMesh = networkLayout.getConstructionObject()
        showKidneys = getShowKidneysSettings()
        isLeftKidney = showKidneys[0]
        isRightKidney = showKidneys[1]

        kidneyTubeNetworkMeshBuilder = KidneyTubeNetworkMeshBuilder(
            networkMesh,
            targetElementDensityAlongLongestSegment=options["Target element density along longest segment"],
            defaultElementsCountAround=options["Number of elements around"],
            elementsCountThroughShell=options["Number of elements through shell"],
            layoutAnnotationGroups=layoutAnnotationGroups,
            isCore=True,
            elementsCountTransition=options["Number of elements across core transition"],
            defaultElementsCountCoreBoxMinor=options["Number of elements across core box minor"],
            annotationElementsCountsCoreBoxMinor=options["Annotation numbers of elements across core box minor"],
            showKidneys=showKidneys
        )

        kidneyTubeNetworkMeshBuilder.build()
        generateData = TubeNetworkMeshGenerateData(
            region, 3,
            isLinearThroughShell=False)
        kidneyTubeNetworkMeshBuilder.generateMesh(generateData)
        annotationGroups = generateData.getAnnotationGroups()

        # add kidney-specific annotation groups
        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)
        mesh = generateData.getMesh()
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        coreGroup = getAnnotationGroupForTerm(annotationGroups, ("core", "")).getGroup()
        shellGroup = getAnnotationGroupForTerm(annotationGroups, ("shell", "")).getGroup()
        openingGroup = getAnnotationGroupForTerm(annotationGroups, ("opening", "")).getGroup()

        kidneyGroup = AnnotationGroup(region, get_kidney_term("kidney"))
        kidneyNodesetGroup = kidneyGroup.getNodesetGroup(nodes)

        leftKidneyGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("left kidney"))
        leftKidneyNodesetGroup = leftKidneyGroup.getNodesetGroup(nodes)

        rightKidneyGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("right kidney"))
        rightKidneyNodesetGroup = rightKidneyGroup.getNodesetGroup(nodes)

        hilumGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("hilum of kidney"))
        hilumGroup.getMeshGroup(mesh).addElementsConditional(openingGroup)

        cortexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("cortex of kidney"))
        tempGroup = fm.createFieldSubtract(shellGroup, openingGroup)
        cortexGroup.getMeshGroup(mesh).addElementsConditional(tempGroup)

        medullaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("renal medulla"))
        medullaGroup.getMeshGroup(mesh).addElementsConditional(coreGroup)

        for term in ["core", "shell", "opening"]:
            annotationGroups.remove(findAnnotationGroupByName(annotationGroups, term))

        # marker points
        leftSuperiorPoleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("superior pole of left kidney"))
        leftInferiorPoleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("inferior pole of left kidney"))

        rightSuperiorPoleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("superior pole of right kidney"))
        rightInferiorPoleGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("inferior pole of right kidney"))

        markerList = []
        elementsCountAround = options["Number of elements around"]
        elementsCountAlong = int(options["Target element density along longest segment"] + 2) # extra sections from the cap mesh
        elementsCountThroughShell = options["Number of elements through shell"]
        elementsCountCoreBoxMinor = options["Number of elements across core box minor"]
        elementsCountCoreBoxMajor = (elementsCountAround // 2) - elementsCountCoreBoxMinor
        elementsCountTransition = options["Number of elements across core transition"]

        box_count = elementsCountCoreBoxMinor * elementsCountCoreBoxMajor
        depth = elementsCountThroughShell + elementsCountTransition
        offset = elementsCountCoreBoxMajor // 2 * elementsCountCoreBoxMinor + elementsCountCoreBoxMinor // 2 + 1
        cap_count = box_count * (depth + 1) + elementsCountAround * depth
        tube_section_count = box_count + elementsCountAround * depth
        tube_count = tube_section_count * elementsCountAlong
        kidney_elements_count = cap_count * 2 + tube_count if showKidneys[0] else 0

        if isLeftKidney:
            idx = box_count * depth + offset
            markerList.append({"group": leftSuperiorPoleGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})

            idx = cap_count + tube_count + (box_count * depth + offset)
            markerList.append({"group": leftInferiorPoleGroup, "elementId": idx, "xi": [0.0, 1.0, 1.0]})

        if isRightKidney:
            idx = kidney_elements_count + box_count * depth + offset
            markerList.append({"group": rightSuperiorPoleGroup, "elementId": idx, "xi": [0.0, 0.0, 1.0]})

            idx = kidney_elements_count + cap_count + tube_count + (box_count * depth + offset)
            markerList.append({"group": rightInferiorPoleGroup, "elementId": idx, "xi": [0.0, 1.0, 1.0]})

        nodeIdentifier = generateData.nextNodeIdentifier()
        for marker in markerList:
            annotationGroup = marker["group"]
            markerNode = annotationGroup.createMarkerNode(
                nodeIdentifier, element=mesh.findElementByIdentifier(marker["elementId"]), xi=marker["xi"])
            annotationGroup.setMarkerMaterialCoordinates(coordinates)
            kidneyNodesetGroup.addNode(markerNode)
            nodeIdentifier += 1

        # transformation
        leftKidney, rightKidney = 0, 1
        kidneys = [lung for show, lung in [(isLeftKidney, leftKidney), (isRightKidney, rightKidney)] if show]
        for kidney in kidneys:
            isLeft = True if kidney == leftKidney else False
            isRight = True if kidney == rightKidney else False
            spacing = -options["Kidney spacing"] / 2 if isLeft else options["Kidney spacing"] / 2
            curvature = options["Kidney curvature"] if isLeft else -options["Kidney curvature"]

            kidneyNodeset = leftKidneyNodesetGroup if isLeft else rightKidneyNodesetGroup

            if curvature != 0.0:
                if isLeft:
                    bendKidneyMeshAroundZAxis(curvature, fm, coordinates, kidneyNodeset, stationaryPointXY=[-0.05, 0.0])
                if isRight:
                    bendKidneyMeshAroundZAxis(curvature, fm, coordinates, kidneyNodeset, stationaryPointXY=[0.05, 0.0])

            translate_nodeset_coordinates(kidneyNodeset, coordinates, [0.0, spacing, 0.0])

        return annotationGroups, None


    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCount = options['Refine number of elements']
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCount, refineElementsCount, refineElementsCount)


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
        show_kidneys = getShowKidneysSettings()

        # Initialize field module and meshes
        fm = region.getFieldmodule()
        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)
        mesh3d = fm.findMeshByDimension(3)
        is_exterior = fm.createFieldIsExterior()

        # Get base kidney group
        kidney_group = getAnnotationGroupForTerm(annotationGroups, get_kidney_term("kidney")).getGroup()

        # Create side groups and tracking dictionaries
        side_groups = {}
        side_kidney_groups = {"left": {}, "right": {}}

        for side in ["left", "right"]:
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(f"{side} kidney"))
            side_groups[side] = group.getGroup()
            side_kidney_groups[side][f"{side} kidney"] = group

        # Create kidney part groups (cortex, hilum, medulla)
        create_kidney_part_groups(fm, mesh3d, annotationGroups, region, side_groups, side_kidney_groups)

        # Create capsule groups
        create_capsule_groups(fm, mesh2d, annotationGroups, region, kidney_group, side_groups, side_kidney_groups, is_exterior)

        # Create surface and edge annotation groups
        create_surface_and_edge_groups(fm, mesh2d, mesh1d, annotationGroups, region, side_groups, side_kidney_groups, is_exterior)

        # Remove groups based on kidney visibility settings
        remove_hidden_kidney_groups(show_kidneys, side_kidney_groups, annotationGroups)


def getShowKidneysSettings():
    return MeshType_1d_kidney_network_layout1.getShowKidneys()


def create_kidney_part_groups(fm, mesh3d, annotationGroups, region, side_groups, side_kidney_groups):
    """
    Create cortex, hilum, and medulla groups for each kidney side.
    """
    kidney_parts = ["cortex", "hilum", "medulla"]

    for part in kidney_parts:
        # Get the anatomical term for the part
        arb_term = f"renal {part}" if part == "medulla" else f"{part} of kidney"
        arb_group = getAnnotationGroupForTerm(annotationGroups, get_kidney_term(arb_term)).getGroup()

        # Create side-specific part groups
        for side in ["left", "right"]:
            part_term = f"{part} of {side} kidney"
            part_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(part_term))
            part_group.getMeshGroup(mesh3d).addElementsConditional(fm.createFieldAnd(arb_group, side_groups[side]))
            side_kidney_groups[side][part_term] = part_group


def create_capsule_groups(fm, mesh2d, annotationGroups, region, kidney_group, side_groups, side_kidney_groups, is_exterior):
    """
    Create kidney capsule surface groups.
    """
    # General kidney capsule group
    kidney_capsule_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term("kidney capsule"))
    kidney_exterior = fm.createFieldAnd(kidney_group, is_exterior)
    kidney_capsule_group.getMeshGroup(mesh2d).addElementsConditional(kidney_exterior)

    # Side-specific capsule groups
    for side in ["left", "right"]:
        capsule_term = f"{side} kidney capsule"
        capsule_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(capsule_term))
        capsule_exterior = fm.createFieldAnd(side_groups[side], is_exterior)
        capsule_group.getMeshGroup(mesh2d).addElementsConditional(capsule_exterior)
        side_kidney_groups[side][capsule_term] = capsule_group


def create_surface_and_edge_groups(fm, mesh2d, mesh1d, annotationGroups, region, side_groups, side_kidney_groups, is_exterior):
    """
    Create surface and edge annotation groups.
    """
    surface_types = ["anterior", "posterior", "lateral", "medial", "dorsal", "ventral", "juxtamedullary cortex"]
    surface_fields = {}

    # Create surface groups
    for surface_type in surface_types:
        if surface_type == "juxtamedullary cortex":
            cortex_group = getAnnotationGroupForTerm(annotationGroups, get_kidney_term("cortex of kidney")).getGroup()
            medulla_group = getAnnotationGroupForTerm(annotationGroups, get_kidney_term("renal medulla")).getGroup()
            surface_exterior = fm.createFieldAnd(medulla_group, cortex_group)
        else:
            base_group = getAnnotationGroupForTerm(annotationGroups, (surface_type, ""))
            surface_field = base_group.getGroup()
            surface_exterior = fm.createFieldAnd(surface_field, is_exterior)
            surface_fields[surface_type] = surface_field

        # General kidney surface group
        general_term = f"{surface_type} surface of kidney"
        general_surface_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(general_term))
        general_surface_group.getMeshGroup(mesh2d).addElementsConditional(surface_exterior)

        # Side-specific surface groups
        for side in ["left", "right"]:
            side_term = f"{surface_type} surface of {side} kidney"
            side_surface_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(side_term))
            side_surface_field = fm.createFieldAnd(surface_exterior, side_groups[side])
            side_surface_group.getMeshGroup(mesh2d).addElementsConditional(side_surface_field)
            side_kidney_groups[side][side_term] = side_surface_group

    # Create edge groups at dorsal-ventral intersection
    dorsal_ventral_border = fm.createFieldAnd(
        fm.createFieldAnd(surface_fields["dorsal"], surface_fields["ventral"]), is_exterior)

    edge_types = ["lateral", "medial"]

    for edge_type in edge_types:
        # General edge group
        edge_term = f"{edge_type} edge of kidney"
        edge_field = fm.createFieldAnd(surface_fields[edge_type], dorsal_ventral_border)
        general_edge_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(edge_term))
        general_edge_group.getMeshGroup(mesh1d).addElementsConditional(edge_field)

        # Side-specific edge groups
        for side in ["left", "right"]:
            side_edge_term = f"{edge_type} edge of {side} kidney"
            side_edge_group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_kidney_term(side_edge_term))
            side_edge_field = fm.createFieldAnd(edge_field, side_groups[side])
            side_edge_group.getMeshGroup(mesh1d).addElementsConditional(side_edge_field)
            side_kidney_groups[side][side_edge_term] = side_edge_group


def remove_hidden_kidney_groups(show_kidneys, side_kidney_groups, annotationGroups):
    """
    Remove annotation groups for kidneys that should not be shown.
    """
    if not show_kidneys[0]:  # Left kidney
        for group in side_kidney_groups["left"].values():
            annotationGroups.remove(group)

    if not show_kidneys[1]:  # Right kidney
        for group in side_kidney_groups["right"].values():
            annotationGroups.remove(group)


def bendKidneyMeshAroundZAxis(curvature, fm, coordinates, kidneyNodeset, stationaryPointXY):
    """
    Transform coordinates by bending with curvature about a centre point the radius in
    x direction from stationaryPointXY.
    :param curvature: 1/radius. Must be non-zero.
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param kidneyNodeset: Zinc NodesetGroup containing nodes to transform.
    :param stationaryPointXY: Coordinates x, y which are not displaced by bending.
    """
    rotateKidneyMeshAboutZAxis(90, fm, coordinates, kidneyNodeset)

    radius = 1.0 / curvature
    scale = fm.createFieldConstant([-1.0, -curvature, -1.0])
    centreOffset = [stationaryPointXY[0] - radius, stationaryPointXY[1], 0.0]
    centreOfCurvature = fm.createFieldConstant(centreOffset)
    polarCoordinates = (centreOfCurvature - coordinates) * scale
    polarCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
    rcCoordinates = fm.createFieldCoordinateTransformation(polarCoordinates)
    rcCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
    newCoordinates = rcCoordinates + centreOfCurvature

    fieldassignment = coordinates.createFieldassignment(newCoordinates)
    fieldassignment.setNodeset(kidneyNodeset)
    fieldassignment.assign()

    rotateKidneyMeshAboutZAxis(-90, fm, coordinates, kidneyNodeset)


def rotateKidneyMeshAboutZAxis(rotateAngle, fm, coordinates, kidneyNodeset):
    """
    Rotates the lung mesh coordinates about a specified axis using the right-hand rule.
    :param rotateAngle: Angle of rotation in degrees.
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param kidneyNodeset: Zinc NodesetGroup containing nodes to transform.
    :return: None
    """
    rotateAngle = -math.radians(rotateAngle)  # negative value due to right handed rule
    rotateMatrix = fm.createFieldConstant([math.cos(rotateAngle), math.sin(rotateAngle), 0.0,
                                           -math.sin(rotateAngle), math.cos(rotateAngle), 0.0,
                                           0.0, 0.0, 1.0])

    rotated_coordinates = fm.createFieldMatrixMultiply(3, rotateMatrix, coordinates)

    fieldassignment = coordinates.createFieldassignment(rotated_coordinates)
    fieldassignment.setNodeset(kidneyNodeset)
    fieldassignment.assign()


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
