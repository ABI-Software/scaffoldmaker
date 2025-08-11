"""
Generates a lung scaffold by deforming a hemisphere.
"""
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.ellipsoidmesh import EllipsoidMesh
from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters, disconnectFieldMeshGroupBoundaryNodes

import math


class MeshType_3d_lung3(Scaffold_base):
    """
    Generates a lung scaffold by deforming a hemisphere.
    """

    @classmethod
    def getName(cls):
        return "3D Lung 3"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1 Coarse",
            "Human 1 Medium",
            "Human 1 Fine",
            "Ellipsoid Coarse",
            "Ellipsoid Medium",
            "Ellipsoid Fine"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        useParameterSetName = "Human 1 Coarse" if (parameterSetName == "Default") else parameterSetName
        options["Left lung"] = True
        options["Right lung"] = True
        # options["Open fissures"] = False
        options["Number of left lung lobes"] = 2
        options["Ellipsoid breadth"] = 0.8
        options["Ellipsoid height"] = 1.0
        options["Ellipsoid depth"] = 0.5
        options["Left-right lung spacing"] = 0.6
        options["Refine"] = False
        options["Refine number of elements"] = 4

        if "Coarse" in useParameterSetName:
            options["Number of elements lateral"] = 4
            options["Number of elements normal"] = 6
            options["Number of elements oblique"] = 6
            options["Number of transition elements"] = 1
        elif "Medium" in useParameterSetName:
            options["Number of elements lateral"] = 4
            options["Number of elements normal"] = 10
            options["Number of elements oblique"] = 10
            options["Number of transition elements"] = 1
        elif "Fine" in useParameterSetName:
            options["Number of elements lateral"] = 6
            options["Number of elements normal"] = 14
            options["Number of elements oblique"] = 14
            options["Number of transition elements"] = 1

        if "Human" in useParameterSetName:
            options["Base lateral edge sharpness factor"] = 0.5
            options["Ventral edge sharpness factor"] = 0.8
            options["Left oblique slope degrees"] = 60.0
            options["Right oblique slope degrees"] = 60.0
            options["Medial curvature"] = 3.0
            options["Medial curvature bias"] = 1.0
            options["Dorsal-ventral rotation degrees"] = 20.0
            options["Ventral-medial rotation degrees"] = 0.0
        else:
            options["Base lateral edge sharpness factor"] = 0.0
            options["Ventral edge sharpness factor"] = 0.0
            options["Left oblique slope degrees"] = 0.0
            options["Right oblique slope degrees"] = 0.0
            options["Medial curvature"] = 0.0
            options["Medial curvature bias"] = 0.0
            options["Dorsal-ventral rotation degrees"] = 0.0
            options["Ventral-medial rotation degrees"] = 0.0

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Left lung",
            "Right lung",
            # "Open fissures",
            "Number of left lung lobes",
            "Number of elements lateral",
            "Number of elements normal",
            "Number of elements oblique",
            "Number of transition elements",
            "Ellipsoid breadth",
            "Ellipsoid height",
            "Ellipsoid depth",
            "Left-right lung spacing",
            "Base lateral edge sharpness factor",
            "Ventral edge sharpness factor",
            "Medial curvature",
            "Medial curvature bias",
            "Dorsal-ventral rotation degrees",
            "Ventral-medial rotation degrees",
            "Left oblique slope degrees",
            "Right oblique slope degrees",
            "Refine",
            "Refine number of elements"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if options["Number of left lung lobes"] > 2:
            options["Number of left lung lobes"] = 2
        elif options["Number of left lung lobes"] < 1:
            options["Number of left lung lobes"] = 0

        max_transition_count = None
        for key in [
            "Number of elements lateral",
            "Number of elements normal",
            "Number of elements oblique"
        ]:
            min_elements_count = 4 if key == "Number of elements lateral" else 6
            if options[key] < min_elements_count:
                options[key] = min_elements_count
            elif options[key] % 2:
                options[key] += 1
            transition_count = (options[key] // 2) - 1
            if (max_transition_count is None) or (transition_count < max_transition_count):
                max_transition_count = transition_count

        if options["Number of transition elements"] < 1:
            options["Number of transition elements"] = 1
        elif options["Number of transition elements"] > max_transition_count:
            options["Number of transition elements"] = max_transition_count
            dependentChanges = True

        for dimension in [
            "Ellipsoid breadth",
            "Ellipsoid height",
            "Ellipsoid depth"
        ]:
            if options[dimension] <= 0.0:
                options[dimension] = 1.0

        if options["Left-right lung spacing"] < 0.0:
            options["Left-right lung spacing"] = 0.0

        for dimension in [
            "Base lateral edge sharpness factor",
            "Ventral edge sharpness factor",
            "Medial curvature bias"
        ]:
            if options[dimension] < 0.0:
                options[dimension] = 0.0
            elif options[dimension] > 1.0:
                options[dimension] = 1.0

        for angle in [
            "Dorsal-ventral rotation degrees",
            "Ventral-medial rotation degrees"
        ]:
            if options[angle] < -90.0:
                options[angle] = -90.0
            elif options[angle] > 90.0:
                options[angle] = 90.0

        if options['Refine number of elements'] < 1:
            options['Refine number of elements'] = 1

        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        isLeftLung = options["Left lung"]
        isRightLung = options["Right lung"]
        # isOpenFissure = options["Open fissures"]
        isOpenFissure = False
        hasAccessoryLobe = False
        numberOfLeftLung = options["Number of left lung lobes"]

        elementsCountLateral = options["Number of elements lateral"]
        elementsCountNormal = options["Number of elements normal"]
        elementsCountOblique = options["Number of elements oblique"]
        elementsCountTransition = options["Number of transition elements"]
        lungSpacing = options["Left-right lung spacing"] * 0.5
        baseSharpFactor = options["Base lateral edge sharpness factor"]
        edgeSharpFactor = options["Ventral edge sharpness factor"]
        ellipsoid_breadth = options["Ellipsoid breadth"]
        ellipsoid_height = options["Ellipsoid height"]
        ellipsoid_depth = options["Ellipsoid depth"]
        left_oblique_slope_radians = math.radians(options["Left oblique slope degrees"])
        right_oblique_slope_radians = math.radians(options["Right oblique slope degrees"])
        leftLungMedialCurvature = options["Medial curvature"]
        lungMedialCurvatureBias = options["Medial curvature bias"]
        rotateLeftLungY = options["Dorsal-ventral rotation degrees"]
        rotateLeftLungZ = options["Ventral-medial rotation degrees"]

        fieldmodule = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fieldmodule)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        mesh = fieldmodule.findMeshByDimension(3)
        meshGroups = []
        annotationTerms = []

        boxGroup = AnnotationGroup(region, ("box group", "None"))
        boxMeshGroup = boxGroup.getMeshGroup(mesh)
        boxNodesetGroup = boxGroup.getNodesetGroup(nodes)

        transitionGroup = AnnotationGroup(region, ("transition group", "None"))
        transitionMeshGroup = transitionGroup.getMeshGroup(mesh)
        transitionNodesetGroup = transitionGroup.getNodesetGroup(nodes)

        # annotation groups & nodeset groups
        lungGroup = AnnotationGroup(region, get_lung_term("lung"))
        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        annotationGroups = [lungGroup]

        leftLungGroup = AnnotationGroup(region, get_lung_term("left lung"))
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)
        annotationGroups.append(leftLungGroup)
        annotationTerms.append("left lung")
        meshGroups.append(leftLungMeshGroup)

        rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightLungGroup)
        annotationTerms.append("right lung")
        meshGroups.append(rightLungMeshGroup)

        leftLateralLungGroup = AnnotationGroup(region, ["lateral left lung", "None"])
        leftLateralLungGroupMeshGroup = leftLateralLungGroup.getMeshGroup(mesh)
        annotationGroups.append(leftLateralLungGroup)

        rightLateralLungGroup = AnnotationGroup(region, ["lateral right lung", "None"])
        rightLateralLungGroupMeshGroup = rightLateralLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightLateralLungGroup)

        leftMedialLungGroup = AnnotationGroup(region, ["medial left lung", "None"])
        leftMedialLungGroupMeshGroup = leftMedialLungGroup.getMeshGroup(mesh)
        annotationGroups.append(leftMedialLungGroup)

        rightMedialLungGroup = AnnotationGroup(region, ["medial right lung", "None"])
        rightMedialLungGroupMeshGroup = rightMedialLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightMedialLungGroup)

        leftAnteriorLungGroup = AnnotationGroup(region, ("anterior left lung", "None"))
        rightAnteriorLungGroup = AnnotationGroup(region, ("anterior right lung", "None"))
        leftPosteriorLungGroup = AnnotationGroup(region, ("posterior left lung", "None"))
        rightPosteriorLungGroup = AnnotationGroup(region, ("posterior right lung", "None"))

        lowerRightLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of right lung"))
        lowerRightLungMeshGroup = lowerRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(lowerRightLungGroup)
        annotationTerms.append("lower lobe of right lung")
        meshGroups.append(lowerRightLungMeshGroup)

        upperRightLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of right lung"))
        upperRightLungMeshGroup = upperRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(upperRightLungGroup)
        annotationTerms.append("upper lobe of right lung")
        meshGroups.append(upperRightLungMeshGroup)

        middleRightLungGroup = AnnotationGroup(region, get_lung_term("middle lobe of right lung"))
        middleRightLungMeshGroup = middleRightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(middleRightLungGroup)
        annotationTerms.append("middle lobe of right lung")
        meshGroups.append(middleRightLungMeshGroup)

        leftBaseLungGroup = AnnotationGroup(region, ["base left lung", "None"])
        leftBaseLungGroupMeshGroup = leftBaseLungGroup.getMeshGroup(mesh)
        annotationGroups.append(leftBaseLungGroup)

        rightBaseLungGroup = AnnotationGroup(region, ["base right lung", "None"])
        rightBaseLungGroupMeshGroup = rightBaseLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightBaseLungGroup)

        if numberOfLeftLung == 2:
            lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
            lowerLeftLungMeshGroup = lowerLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(lowerLeftLungGroup)
            annotationTerms.append("lower lobe of left lung")
            meshGroups.append(lowerLeftLungMeshGroup)

            upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
            upperLeftLungMeshGroup = upperLeftLungGroup.getMeshGroup(mesh)
            annotationGroups.append(upperLeftLungGroup)
            annotationTerms.append("upper lobe of left lung")
            meshGroups.append(upperLeftLungMeshGroup)
        #
        # if hasAccessoryLobe:  # currently not used
        #     # Annotation groups
        #     rightLungAccessoryLobeGroup = AnnotationGroup(region, get_lung_term("right lung accessory lobe"))
        #     rightLungAccessoryLobeMeshGroup = rightLungAccessoryLobeGroup.getMeshGroup(mesh)
        #     annotationGroups.append(rightLungAccessoryLobeGroup)
        #     rightLungAccessoryLobeNodesetGroup = rightLungAccessoryLobeGroup.getNodesetGroup(nodes)
        #
        #     # Marker points
        #     accessoryDorsalApexGroup = findOrCreateAnnotationGroupForTerm(
        #         annotationGroups, region, get_lung_term("dorsal apex of right lung accessory lobe"))
        #     accessoryVentralApexGroup = findOrCreateAnnotationGroupForTerm(
        #         annotationGroups, region, get_lung_term("ventral apex of right lung accessory lobe"))
        #     accessoryVentralLeftGroup = findOrCreateAnnotationGroupForTerm(
        #         annotationGroups, region, get_lung_term("left ventral base of right lung accessory lobe"))
        #     accessoryVentralRightGroup = findOrCreateAnnotationGroupForTerm(
        #         annotationGroups, region, get_lung_term("right ventral base of right lung accessory lobe"))
        #     accessoryDorsalLeftGroup = findOrCreateAnnotationGroupForTerm(
        #         annotationGroups, region, get_lung_term("left dorsal base of right lung accessory lobe"))
        #     accessoryDorsalRightGroup = findOrCreateAnnotationGroupForTerm(
        #         annotationGroups, region, get_lung_term("right dorsal base of right lung accessory lobe"))

        # # Marker points/groups # Currently not used
        # leftApexGroup = findOrCreateAnnotationGroupForTerm(
        #     annotationGroups, region, get_lung_term("apex of left lung"))
        # rightApexGroup = findOrCreateAnnotationGroupForTerm(
        #     annotationGroups, region, get_lung_term("apex of right lung"))
        # leftVentralGroup = findOrCreateAnnotationGroupForTerm(
        #     annotationGroups, region, get_lung_term("ventral base of left lung"))
        # rightVentralGroup = findOrCreateAnnotationGroupForTerm(
        #     annotationGroups, region, get_lung_term("ventral base of right lung"))
        # rightLateralGroup = findOrCreateAnnotationGroupForTerm(
        #     annotationGroups, region, get_lung_term("laterodorsal tip of middle lobe of right lung"))
        # leftMedialGroup = findOrCreateAnnotationGroupForTerm(
        #     annotationGroups, region, get_lung_term("medial base of left lung"))
        # rightMedialGroup = findOrCreateAnnotationGroupForTerm(
        #     annotationGroups, region, get_lung_term("medial base of right lung"))

        # Nodeset group
        leftLungNodesetGroup = leftLungGroup.getNodesetGroup(nodes)
        rightLungNodesetGroup = rightLungGroup.getNodesetGroup(nodes)

        elementCounts = [elementsCountLateral, elementsCountOblique, elementsCountNormal]
        halfDepth = ellipsoid_depth * 0.5
        halfBreadth = ellipsoid_breadth * 0.5
        halfHeight = ellipsoid_height * 0.5
        surface_only = False

        leftLung, rightLung = 0, 1
        lungs = [lung for show, lung in [(isLeftLung, leftLung), (isRightLung, rightLung)] if show]
        startNodeIdentifier, startElementIdentifier = 1, 1
        for lung in lungs:
            oblique_slope_radians = left_oblique_slope_radians if lung == leftLung else right_oblique_slope_radians
            axis2_x_rotation_radians = -oblique_slope_radians
            axis3_x_rotation_radians = math.radians(90) - oblique_slope_radians

            ellipsoid = EllipsoidMesh(halfDepth, halfBreadth, halfHeight, elementCounts, elementsCountTransition,
                                      axis2_x_rotation_radians, axis3_x_rotation_radians, surface_only)

            annotationGroups = [lungGroup, leftLungGroup, rightLungGroup,
                                leftLateralLungGroup, leftMedialLungGroup, leftBaseLungGroup, leftPosteriorLungGroup,
                                rightLateralLungGroup, rightMedialLungGroup, rightBaseLungGroup, rightPosteriorLungGroup,
                                lowerLeftLungGroup, upperLeftLungGroup,
                                lowerRightLungGroup, middleRightLungGroup, upperRightLungGroup]

            if lung == leftLung:
                if numberOfLeftLung == 2:
                    octant_group_lists = []
                    for octant in range(8):
                        octant_group_list = []
                        octant_group_list.append(lungGroup.getGroup())
                        octant_group_list.append(leftLungGroup.getGroup())
                        octant_group_list.append((leftMedialLungGroup if (octant & 1) else leftLateralLungGroup).getGroup())
                        octant_group_list.append((leftBaseLungGroup if (octant & 2) else leftPosteriorLungGroup).getGroup())
                        octant_group_list.append((upperLeftLungGroup if (octant & 4) else lowerLeftLungGroup).getGroup())
                        octant_group_lists.append(octant_group_list)
                else:
                    octant_group_lists = []
                    for octant in range(8):
                        octant_group_list = []
                        octant_group_list.append(lungGroup.getGroup())
                        octant_group_list.append(leftLungGroup.getGroup())
            else:
                octant_group_lists = []
                for octant in range(8):
                    octant_group_list = []
                    octant_group_list.append(lungGroup.getGroup())
                    octant_group_list.append(rightLungGroup.getGroup())
                    octant_group_list.append((rightMedialLungGroup if (octant & 1) else rightLateralLungGroup).getGroup())
                    octant_group_list.append((rightBaseLungGroup if (octant & 2) else rightPosteriorLungGroup).getGroup())
                    if octant & 4:
                        octant_group_list.append((middleRightLungGroup if (octant & 2) else upperRightLungGroup).getGroup())
                    else:
                        octant_group_list.append(lowerRightLungGroup.getGroup())
                    octant_group_lists.append(octant_group_list)

            ellipsoid.set_octant_group_lists(octant_group_lists)
            annotationGroups += [boxGroup, transitionGroup]
            ellipsoid.set_box_transition_groups(boxGroup.getGroup(), transitionGroup.getGroup())

            ellipsoid.build()
            nodeIdentifier, elementIdentifier = ellipsoid.generate_mesh(fieldmodule, coordinates, startNodeIdentifier, startElementIdentifier)

            startNodeIdentifier, startElementIdentifier = nodeIdentifier, elementIdentifier

        for lung in lungs:
            isLeft = True if lung == leftLung else False
            lungNodeset = leftLungNodesetGroup if isLeft else rightLungNodesetGroup
            spacing = -lungSpacing if lung == leftLung else lungSpacing
            zOffset = -0.5 * ellipsoid_height
            lungMedialCurvature = -leftLungMedialCurvature if isLeft else leftLungMedialCurvature
            rotateLungAngleY = rotateLeftLungY if isLeft else -rotateLeftLungY
            rotateLungAngleZ = rotateLeftLungZ if isLeft else -rotateLeftLungZ

            if not isLeft:
                mirrorLeftLungToRight(fieldmodule, coordinates, lungNodeset)

            if edgeSharpFactor != 0.0:
                taperLungEdge(edgeSharpFactor, fieldmodule, coordinates, lungNodeset, halfBreadth)

            if baseSharpFactor != 0.0:
                taperLungEdge(baseSharpFactor, fieldmodule, coordinates, lungNodeset, halfHeight, isBase=True)

            dorsalVentralXi = getDorsalVentralXiField(fieldmodule, coordinates, halfBreadth)
            if lungMedialCurvature != 0:
                bendLungMeshAroundZAxis(lungMedialCurvature, fieldmodule, coordinates, lungNodeset,
                                        stationaryPointXY=[0.0, 0.0],
                                        bias=lungMedialCurvatureBias,
                                        dorsalVentralXi=dorsalVentralXi)

            if rotateLungAngleY != 0.0:
                rotateLungMeshAboutAxis(rotateLungAngleY, fieldmodule, coordinates, lungNodeset, axis=2)

            if rotateLungAngleZ != 0.0:
                rotateLungMeshAboutAxis(rotateLungAngleZ, fieldmodule, coordinates, lungNodeset, axis=3)

            translateLungLocation(fieldmodule, coordinates, lungNodeset, spacing, zOffset)

            # if isOpenFissure:
            #     setValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
            #                       Node.VALUE_LABEL_D_DS3]
            #     nodeParameters = get_nodeset_field_parameters(lungNodeset, coordinates, setValueLabels)[1]
            #     nodeIdentifier = nodeParameters[-1][0] + 1
            #     # if numberOfLeftLung > 1:
            #     nodeIdentifier, copyIdentifiersLLU = disconnectFieldMeshGroupBoundaryNodes(
            #         [coordinates], lowerLeftLungMeshGroup, upperLeftLungMeshGroup,
            #         nodeIdentifier)
            #     nodeIdentifier, copyIdentifiersRLM = disconnectFieldMeshGroupBoundaryNodes(
            #         [coordinates], lowerRightLungMeshGroup, middleRightLungMeshGroup,
            #         nodeIdentifier)
            #     nodeIdentifier, copyIdentifiersRLU = disconnectFieldMeshGroupBoundaryNodes(
            #         [coordinates], lowerRightLungMeshGroup, upperRightLungMeshGroup, nodeIdentifier)
            #     nodeIdentifier, copyIdentifiersRMU = disconnectFieldMeshGroupBoundaryNodes(
            #         [coordinates], middleRightLungMeshGroup, upperRightLungMeshGroup,
            #         nodeIdentifier)

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
        numberOfLeftLung = options['Number of left lung lobes']
        # hasAccessoryLobe = options['Accessory lobe']
        hasAccessoryLobe = False
        # openFissures = options['Open fissures']
        openFissures = False

        # create fissure groups
        fm = region.getFieldmodule()
        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)

        # 1D Annotation
        is_exterior = fm.createFieldIsExterior()
        is_xi1_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0)
        is_xi1_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1)
        is_xi2_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0)
        is_xi2_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1)
        is_xi3_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)
        is_xi3_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)

        # Arbitrary terms - are removed from the annotation groups later
        arbLobe_group = {}
        arbLobe_exterior = {}
        arbLobe_2dgroup = {}
        arbTerms = ["upper lobe of left lung", "upper lobe of right lung"]
        for arbTerm in arbTerms:
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(arbTerm))
            group2d = group.getGroup()
            group2d_exterior = fm.createFieldAnd(group2d, is_exterior)
            arbLobe_group.update({arbTerm: group})
            arbLobe_2dgroup.update({arbTerm: group2d})
            arbLobe_exterior.update({arbTerm: group2d_exterior})

        side_group = {}
        side_exterior = {}
        arbSideTerms = ["lateral left lung", "lateral right lung",
                        "medial left lung", "medial right lung",
                        "base left lung", "base right lung"]
        for term in arbSideTerms:
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [term, "None"])
            group2d = group.getGroup()
            group2d_exterior = fm.createFieldAnd(group2d, is_exterior)
            side_group[term] = group
            side_exterior[term] = group2d_exterior

        boxTransition_group = {}
        boxTransition_exterior = {}
        arbBoxTransitionTerms = ["box group", "transition group"]
        for term in arbBoxTransitionTerms:
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [term, "None"])
            group2d = group.getGroup()
            group2d_exterior = fm.createFieldAnd(group2d, is_exterior)
            boxTransition_group[term] = group
            boxTransition_exterior[term] = group2d_exterior

        base_posterior_group = {}
        base_posterior_group_exterior = {}
        base_posterior_group_terms = ["base left lung", "base right lung",
                                      "posterior left lung", "posterior right lung"]
        for term in base_posterior_group_terms:
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [term, "None"])
            group2d = group.getGroup()
            group2d_exterior = fm.createFieldAnd(group2d, is_exterior)
            base_posterior_group[term] = group
            base_posterior_group_exterior[term] = group2d_exterior

        # Exterior surfaces of lungs
        surfaceTerms = [
            "left lung",
            "lower lobe of left lung",
            "upper lobe of left lung",
            "right lung",
            "lower lobe of right lung",
            "middle lobe of right lung",
            "upper lobe of right lung"
            # "right lung accessory lobe"
        ]
        subLeftLungTerms = ["lower lobe of left lung", "upper lobe of left lung"]

        lobe = {}
        lobe_exterior = {}
        for term in surfaceTerms:
            if (numberOfLeftLung == 1) and (term in subLeftLungTerms):
                continue
            if (hasAccessoryLobe is False) and (term == "right lung accessory lobe"):
                continue

            group = getAnnotationGroupForTerm(annotationGroups, get_lung_term(term))
            group2d = group.getGroup()
            group2d_exterior = fm.createFieldAnd(group2d, is_exterior)

            surfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(term + " surface"))
            if (not openFissures) or (term == "right lung accessory lobe"):
                surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(group2d_exterior)

            lobe_exterior.update({term + " surface": group2d_exterior})

            if "lobe of" in term:
                lobe.update({term: group2d})

            # lateral in the subgroup
            for sideTerm in ['lateral surface of ', 'medial surface of ']:
                surfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(sideTerm + term))
                if ('lateral' in sideTerm) and ('left' in term):
                    surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
                        fm.createFieldAnd(group2d_exterior, side_exterior["lateral left lung"]))
                elif ('lateral' in sideTerm) and ('right' in term):
                    surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
                        fm.createFieldAnd(group2d_exterior, side_exterior["lateral right lung"]))
                elif ('medial' in sideTerm) and ('left' in term):
                    surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
                        fm.createFieldAnd(group2d_exterior, side_exterior["medial left lung"]))
                elif ('medial' in sideTerm) and ('right' in term):
                    surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
                        fm.createFieldAnd(group2d_exterior, side_exterior["medial right lung"]))

        # Base surface of lungs (incl. lobes)
        baseGroup = []
        baseTerms = [
            'left lung surface',
            'lower lobe of right lung surface',
            'middle lobe of right lung surface',
            'right lung surface'
        ]

        # Base of left lung
        if numberOfLeftLung > 1:
            baseTerms = ['lower lobe of left lung surface', 'upper lobe of left lung surface'] + baseTerms

            tempGroup = fm.createFieldAnd(lobe_exterior[baseTerms[0]], side_exterior["base left lung"])
            baseLeftLowerLung = fm.createFieldAnd(tempGroup, side_exterior["medial left lung"])
            baseGroup.append(baseLeftLowerLung)

            tempGroup = fm.createFieldAnd(lobe_exterior[baseTerms[1]], side_exterior["base left lung"])
            baseLeftUpperLung = fm.createFieldAnd(tempGroup, side_exterior["medial left lung"])
            baseGroup.append(baseLeftUpperLung)

            baseLeftLung = fm.createFieldOr(baseLeftLowerLung, baseLeftUpperLung)
            baseGroup.append(baseLeftLung)
        else:
            baseLeftLung = side_exterior["base left lung"]
            baseGroup.append(baseLeftLung)

        # Base of right lung
        tempGroup = fm.createFieldAnd(lobe_exterior['lower lobe of right lung surface'], side_exterior["base right lung"])
        baseRightLowerLung = fm.createFieldAnd(tempGroup, side_exterior["medial right lung"])
        baseGroup.append(baseRightLowerLung)

        tempGroup = fm.createFieldAnd(lobe_exterior['middle lobe of right lung surface'], side_exterior["base right lung"])
        baseRightMiddleLung = fm.createFieldAnd(tempGroup, side_exterior["medial right lung"])
        baseGroup.append(baseRightMiddleLung)

        baseRightLung = fm.createFieldOr(baseRightLowerLung, baseRightMiddleLung)
        baseGroup.append(baseRightLung)

        for term in baseTerms:
            baseSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term("base of " + term))
            baseSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(baseGroup[baseTerms.index(term)])

        # Fissures
        if not openFissures:
            fissureTerms = ["oblique fissure of right lung", "horizontal fissure of right lung"]
            if numberOfLeftLung > 1:
                fissureTerms.append("oblique fissure of left lung")
            lobeFissureTerms = [
                "oblique fissure of lower lobe of left lung",
                "oblique fissure of upper lobe of left lung",
                "oblique fissure of lower lobe of right lung",
                "oblique fissure of middle lobe of right lung",
                "oblique fissure of upper lobe of right lung",
                "horizontal fissure of middle lobe of right lung",
                "horizontal fissure of upper lobe of right lung"
            ]
            for fissureTerm in fissureTerms:
                if (fissureTerm == "oblique fissure of left lung") and (numberOfLeftLung > 1):
                    fissureGroup = fm.createFieldAnd(lobe["upper lobe of left lung"], lobe["lower lobe of left lung"])
                elif fissureTerm == "oblique fissure of right lung":
                    fissureGroup = fm.createFieldAnd(
                        fm.createFieldOr(lobe["middle lobe of right lung"], lobe["upper lobe of right lung"]),
                        lobe["lower lobe of right lung"])
                elif fissureTerm == "horizontal fissure of right lung":
                    fissureGroup = fm.createFieldAnd(
                        lobe["upper lobe of right lung"], lobe["middle lobe of right lung"])

                fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(fissureTerm))
                fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(fissureGroup)
                fissureGroup_temp = fissureGroup

                for lobeFissureTerm in lobeFissureTerms:
                    temp_splitTerm = fissureTerm.split("of")
                    if (temp_splitTerm[0] in lobeFissureTerm) and (temp_splitTerm[1] in lobeFissureTerm):
                        if "oblique fissure of upper lobe of right lung" in lobeFissureTerm:
                            fissureGroup = fm.createFieldAnd(fissureGroup_temp,
                                                             arbLobe_2dgroup['upper lobe of right lung'])
                        elif "oblique fissure of middle lobe of right lung" in lobeFissureTerm:
                            fissureGroup = fm.createFieldAnd(fissureGroup_temp, fm.createFieldNot(
                                arbLobe_2dgroup['upper lobe of right lung']))
                        fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(lobeFissureTerm))
                        fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(fissureGroup)

            # add fissures to lobe surface groups
            if numberOfLeftLung > 1:
                obliqueFissureOfLeftLungGroup = getAnnotationGroupForTerm(
                    annotationGroups, get_lung_term("oblique fissure of left lung")).getGroup()
                for lobeSurfaceTerm in ("lower lobe of left lung surface", "upper lobe of left lung surface"):
                    lobeSurfaceGroup = getAnnotationGroupForTerm(
                        annotationGroups, get_lung_term(lobeSurfaceTerm))
                    lobeSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueFissureOfLeftLungGroup)
            horizontalFissureOfRightLungGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("horizontal fissure of right lung")).getGroup()
            obliqueFissureOfRightLungGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("oblique fissure of right lung")).getGroup()
            obliqueFissureOfMiddleLobeOfRightLungGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("oblique fissure of middle lobe of right lung")).getGroup()
            obliqueFissureOfUpperLobeOfRightLungGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("oblique fissure of upper lobe of right lung")).getGroup()
            lobeSurfaceGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("lower lobe of right lung surface"))
            lobeSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueFissureOfRightLungGroup)
            lobeSurfaceGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("middle lobe of right lung surface"))
            lobeSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
                fm.createFieldOr(obliqueFissureOfMiddleLobeOfRightLungGroup, horizontalFissureOfRightLungGroup))
            lobeSurfaceGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("upper lobe of right lung surface"))
            lobeSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
                fm.createFieldOr(obliqueFissureOfUpperLobeOfRightLungGroup, horizontalFissureOfRightLungGroup))
    #
    #     if openFissures:
    #         if numberOfLeftLung > 1:
    #             obliqueLowerLeft = fm.createFieldOr(
    #                 fm.createFieldAnd(lobe_exterior['lower lobe of left lung surface'], is_xi2_1),
    #                 fm.createFieldAnd(lobe_exterior['lower lobe of left lung surface'], is_xi3_1))
    #
    #             tempUpperLeft = fm.createFieldXor(
    #                 fm.createFieldAnd(lobe_exterior['upper lobe of left lung surface'], is_exterior),
    #                 arbLobe_exterior['upper lobe of left lung dorsal'])
    #
    #             obliqueUpperLeft = fm.createFieldOr(
    #                 fm.createFieldAnd(lobe_exterior['upper lobe of left lung surface'], is_xi2_0),
    #                 tempUpperLeft)
    #
    #             fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #                 'oblique fissure of ' + 'lower lobe of left lung'))
    #             fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueLowerLeft)
    #
    #             fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #                 'oblique fissure of ' + 'upper lobe of left lung'))
    #             fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueUpperLeft)
    #
    #             obliqueLeft = fm.createFieldOr(obliqueUpperLeft, obliqueLowerLeft)
    #             fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #                 'oblique fissure of ' + 'left lung'))
    #             fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueLeft)
    #
    #             leftLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                                       get_lung_term('left lung surface'))
    #             leftLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
    #                 fm.createFieldAnd(lobe_exterior['left lung surface'], fm.createFieldNot(obliqueLeft)))
    #
    #             for term in ("lower lobe of left lung surface", "upper lobe of left lung surface"):
    #                 leftLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(
    #                     annotationGroups, region, get_lung_term(term))
    #                 leftLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(lobe_exterior[term])
    #         else:
    #             leftLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                                       get_lung_term('left lung surface'))
    #             leftLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(lobe_exterior['left lung surface'])
    #
    #         obliqueLowerRight = fm.createFieldOr(
    #             fm.createFieldAnd(lobe_exterior['lower lobe of right lung surface'], is_xi2_1),
    #             fm.createFieldAnd(lobe_exterior['lower lobe of right lung surface'], is_xi3_1))
    #
    #         tempObliqueMiddleRight1 = fm.createFieldAnd(
    #             fm.createFieldAnd(lobe_exterior['middle lobe of right lung surface'], is_xi2_0),
    #             fm.createFieldNot(boxTransition_exterior["box group"]))
    #
    #         tempObliqueMiddleRight2 = fm.createFieldAnd(
    #             fm.createFieldAnd(lobe_exterior['middle lobe of right lung surface'], is_xi3_0),
    #             boxTransition_exterior["box group"])
    #
    #         obliqueMiddleRight = fm.createFieldOr(tempObliqueMiddleRight1, tempObliqueMiddleRight2)
    #
    #         tempObliqueUpperRight1 = fm.createFieldAnd(
    #             fm.createFieldAnd(lobe_exterior['upper lobe of right lung surface'], is_xi3_0),
    #             boxTransition_exterior["box group"])
    #
    #         obliqueUpperRight = fm.createFieldOr(
    #             tempObliqueUpperRight1,
    #             fm.createFieldAnd(arbLobe_exterior['upper lobe of right lung dorsal'], is_xi2_0))
    #
    #         obliqueRight = fm.createFieldOr(fm.createFieldOr(obliqueLowerRight, obliqueUpperRight), obliqueMiddleRight)
    #
    #         tempHorizontalMiddleRight1 = fm.createFieldAnd(
    #             lobe_exterior['middle lobe of right lung surface'],
    #             fm.createFieldNot(obliqueMiddleRight))
    #
    #         tempHorizontalMiddleRight2 = fm.createFieldAnd(lobe_exterior['middle lobe of right lung surface'], is_xi3_0)
    #
    #         horizontalMiddleRight = fm.createFieldAnd(tempHorizontalMiddleRight1,
    #                                                   fm.createFieldNot(tempHorizontalMiddleRight2))
    #
    #         tempHorizontalUpperRight1 = fm.createFieldOr(
    #             fm.createFieldAnd(lobe_exterior['upper lobe of right lung surface'], is_xi1_0),
    #             fm.createFieldAnd(lobe_exterior['upper lobe of right lung surface'], is_xi1_1))
    #
    #         horizontalUpperRight = fm.createFieldOr(
    #             tempHorizontalUpperRight1,
    #             fm.createFieldAnd(lobe_exterior['upper lobe of right lung surface'], is_xi2_1))
    #
    #         horizontalRight = fm.createFieldOr(horizontalMiddleRight, horizontalUpperRight)
    #
    #         fissureRight = fm.createFieldOr(obliqueRight, horizontalRight)
    #
    #         rightLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                                    get_lung_term('right lung surface'))
    #         rightLungSurface = fm.createFieldAnd(lobe_exterior['right lung surface'], fm.createFieldNot(fissureRight))
    #         rightLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(rightLungSurface)
    #
    #         for term in ("lower lobe of right lung surface", "middle lobe of right lung surface",
    #                      "upper lobe of right lung surface"):
    #             rightLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(
    #                 annotationGroups, region, get_lung_term(term))
    #             rightLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(lobe_exterior[term])
    #
    #         fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #             'oblique fissure of ' + 'lower lobe of right lung'))
    #         fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueLowerRight)
    #
    #         fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #             'oblique fissure of ' + 'middle lobe of right lung'))
    #         fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueMiddleRight)
    #
    #         fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #             'oblique fissure of ' + 'upper lobe of right lung'))
    #         fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueUpperRight)
    #
    #         fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #             'horizontal fissure of ' + 'middle lobe of right lung'))
    #         fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(horizontalMiddleRight)
    #
    #         fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #             'horizontal fissure of ' + 'upper lobe of right lung'))
    #         fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(horizontalUpperRight)
    #
    #         fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #             'oblique fissure of ' + 'right lung'))
    #         fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueRight)
    #
    #         fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(
    #             'horizontal fissure of ' + 'right lung'))
    #         fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(horizontalRight)

        # 1D edges
        edgeTerms = [
            "anterior edge",
            "antero-posterior edge",
            "base edge",
            "lateral edge",
            "medial edge",
            "posterior edge"
        ]

        fissureTerms = [
            "horizontal fissure",
            "oblique fissure",
        ]

        lobeTerms = [
            "lower lobe",
            "middle lobe",
            "upper lobe"
        ]

        # Define mappings
        edge_lobe_map = {
            "anterior edge": ["middle"],
            "antero-posterior edge": ["upper"],
            "posterior edge": ["lower"]
        }
        surface_edge_terms = ["lateral edge", "medial edge", "base edge"]

        for lung in ["left lung", "right lung"]:
            surfaces = {}
            for surface_type in ["lateral", "medial", "base"]:
                term = f"{surface_type} of {lung} surface" if surface_type == "base" else f"{surface_type} surface of {lung}"
                group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(term))
                surfaces[surface_type] = group.getGroup()

            for edgeTerm in edgeTerms:
                if edgeTerm in surface_edge_terms:
                    surface_type = edgeTerm.split()[0]  # Extract "lateral", "medial", or "base"

                    for fissure in fissureTerms:
                        if "horizontal" in fissure and ("left" in lung or "base" in edgeTerm):
                            continue

                        fissureTerm = f"{fissure} of {lung}"
                        fissureGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                          get_lung_term(fissureTerm))
                        is_fissureGroup = fissureGroup.getGroup()
                        is_fissureGroup_exterior = fm.createFieldAnd(is_fissureGroup, is_exterior)
                        is_surfaceFissureGroup = fm.createFieldAnd(is_fissureGroup_exterior, surfaces[surface_type])

                        edgeTermFull = f"{edgeTerm} of oblique fissure" if "base" in edgeTerm else f"{edgeTerm} of {fissureTerm}"
                        edgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [edgeTermFull, "None"])
                        edgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_surfaceFissureGroup)

                elif edgeTerm in edge_lobe_map:
                    for lobe in lobeTerms:
                        if (("middle" in lobe and "left" in lung) or
                                not any(valid_lobe in lobe for valid_lobe in edge_lobe_map[edgeTerm])):
                            continue

                        lobeTerm = f"{lobe} of {lung}"
                        lobeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                       get_lung_term(lobeTerm))
                        is_lobeGroup = lobeGroup.getGroup()
                        is_lobeGroup_exterior = fm.createFieldAnd(is_lobeGroup, is_exterior)

                        is_edge = fm.createFieldAnd(is_lobeGroup_exterior,
                                                    fm.createFieldAnd(surfaces["medial"], surfaces["lateral"]))

                        edgeTermFull = f"{edgeTerm} of {lobe} of {lung}"
                        edgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [edgeTermFull, "None"])
                        edgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_edge)

            # Process lateral edges of the base (separate from main edge loop)
            for lobe in lobeTerms:
                if ("middle" in lobe and "left" in lung) or ("upper" in lobe and "right" in lung):
                    continue

                lobeTerm = f"base of {lobe} of {lung} surface"
                lobeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(lobeTerm))
                is_lobeGroup = lobeGroup.getGroup()
                is_lobeGroup_exterior = fm.createFieldAnd(is_lobeGroup, is_exterior)

                sideTerm = f"lateral {lung}"
                is_edge = fm.createFieldAnd(is_lobeGroup_exterior, side_exterior[sideTerm])

                edgeTermFull = f"lateral edge of {lobeTerm.replace(' surface', '')}"
                edgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [edgeTermFull, "None"])
                edgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_edge)

        # Medial edge of the base
        halfBreadth = options["Ellipsoid breadth"] * -0.5
        coordinates = find_or_create_field_coordinates(fm)

        for lung in ["left lung", "right lung"]:
            is_base = base_posterior_group_exterior[f"base {lung}"]
            is_posterior = base_posterior_group_exterior[f"posterior {lung}"]
            is_medial = side_exterior[f"medial {lung}"]

            slope_key = "Left oblique slope degrees" if "left" in lung else "Right oblique slope degrees"
            rotationAngle = math.radians(90 - options[slope_key])

            is_horizontal_edge = fm.createFieldAnd(fm.createFieldAnd(is_base, is_posterior), is_medial)
            is_threshold = setBaseGroupThreshold(fm, coordinates, halfBreadth, rotationAngle)
            is_edge = fm.createFieldAnd(is_horizontal_edge, is_threshold)

            edgeTerm = f"medial edge of base of lower lobe of {lung}"
            edgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [edgeTerm, "None"])
            edgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_edge)

        # Remove unnecessary annotations
        for group in [*side_group.values(), *boxTransition_group.values()]:
            annotationGroups.remove(group)

        for key, group in base_posterior_group.items():
            if "posterior" in key:
                annotationGroups.remove(group)


def mirrorLeftLungToRight(fieldmodule, coordinates, lungNodesetGroup):
    """
    The initial right lung is identical to the left lung. This function flips the right lung to be the mirror image
    of the left lung.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :return: None.
    """
    x = fieldmodule.createFieldComponent(coordinates, 1)
    y = fieldmodule.createFieldComponent(coordinates, 2)
    z = fieldmodule.createFieldComponent(coordinates, 3)

    minus_one = fieldmodule.createFieldConstant(-1.0)
    mirrored_x = x * minus_one

    new_coordinates = fieldmodule.createFieldConcatenate([mirrored_x, y, z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def translateLungLocation(fm, coordinates, lungNodesetGroup, spaceFromCentre, zOffset):
    """
    Translates the lung mesh along the x- and z-axes to adjust its position within the anatomical space.
    :param fm: Field module used for creating and managing fields.
    :param coordinates: The coordinate field representing the initial position of the lung, with a circular profile in the y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing the nodes of the lung to be transformed.
    :param spaceFromCentre: The distance to translate the lung along the x-axis, allowing lateral positioning relative
    to the centerline (i.e. x = 0) of the body.
    :param zOffset: The amount to translate the lung along the z-axis, so that the base lateral edge is at z = 0.
    :return: None
    """
    x = fm.createFieldComponent(coordinates, 1)
    y = fm.createFieldComponent(coordinates, 2)
    z = fm.createFieldComponent(coordinates, 3)
    spaceFromCentre = fm.createFieldConstant(spaceFromCentre)

    new_x = x + spaceFromCentre
    new_z = z - fm.createFieldConstant(zOffset)
    new_coordinates = fm.createFieldConcatenate([new_x, y, new_z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def rotateLungMeshAboutAxis(rotateAngle, fm, coordinates, lungNodesetGroup, axis):
    """
    Rotates the lung mesh coordinates about a specified axis using the right-hand rule.
    :param rotateAngle: Angle of rotation in degrees.
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param axis: Axis of rotation.
    :return: None
    """
    if axis not in (2, 3):
        raise ValueError("Axis must be 2 (y), or 3 (z).")

    rotateAngle = -math.radians(rotateAngle)  # negative value due to right handed rule

    if axis == 2:
        rotateMatrix = fm.createFieldConstant([math.cos(rotateAngle), 0.0, -math.sin(rotateAngle),
                                               0.0, 1.0, 0.0,
                                               math.sin(rotateAngle), 0.0, math.cos(rotateAngle)])
    elif axis == 3:
        rotateMatrix = fm.createFieldConstant([math.cos(rotateAngle), math.sin(rotateAngle), 0.0,
                                               -math.sin(rotateAngle), math.cos(rotateAngle), 0.0,
                                               0.0, 0.0, 1.0])

    rotated_coordinates = fm.createFieldMatrixMultiply(3, rotateMatrix, coordinates)

    fieldassignment = coordinates.createFieldassignment(rotated_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def getDorsalVentralXiField(fm, coordinates, halfBreadth):
    """
    Get a field varying from 0.0 on dorsal tip to 1.0 on ventral tip on [-axisLength, axisLength]
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param halfBreadth: Half breadth of lung.
    :return: Scalar Xi field.
    """
    hl = fm.createFieldConstant(halfBreadth)
    fl = fm.createFieldConstant(2.0 * halfBreadth)
    y = fm.createFieldComponent(coordinates, 2)
    return (y + hl) / fl


def bendLungMeshAroundZAxis(curvature, fm, coordinates, lungNodesetGroup, stationaryPointXY, bias=0.0,
                            dorsalVentralXi=None):
    """
    Transform coordinates by bending with curvature about a centre point the radius in
    x direction from stationaryPointXY.
    :param curvature: 1/radius. Must be non-zero.
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param stationaryPointXY: Coordinates x, y which are not displaced by bending.
    :param bias: 0.0 for a simple bend through the whole length, up to 1.0 for no bend at dorsal end.
    :param dorsalVentralXi: Field returned by getDorsalVentralXiField if bias > 0.0:
    """
    radius = 1.0 / curvature
    scale = fm.createFieldConstant([-1.0, -curvature, -1.0])
    centreOffset = [stationaryPointXY[0] - radius, stationaryPointXY[1], 0.0]
    centreOfCurvature = fm.createFieldConstant(centreOffset)
    polarCoordinates = (centreOfCurvature - coordinates) * scale
    polarCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_CYLINDRICAL_POLAR)
    rcCoordinates = fm.createFieldCoordinateTransformation(polarCoordinates)
    rcCoordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
    newCoordinates = rcCoordinates + centreOfCurvature

    if bias > 0.0:
        one = fm.createFieldConstant(1.0)
        xiS = (one - dorsalVentralXi) * fm.createFieldConstant(bias)
        xiC = one - xiS
        newCoordinates = (coordinates * xiS) + (newCoordinates * xiC)

    fieldassignment = coordinates.createFieldassignment(newCoordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def taperLungEdge(sharpeningFactor, fm, coordinates, lungNodesetGroup, halfValue, isBase=False):
    """
    Applies a tapering transformation to the lung geometry to sharpen the anterior edge or the base.
    If isBase is False, it sharpens the anterior edge (along the y-axis).
    If isBase is True, it sharpens the base (along the z-axis), but only for nodes below a certain height.
    :param sharpeningFactor: A value between 0 and 1, where 1 represents the maximum sharpness.
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param halfValue: Half value of lung breadth/height depending on isBase.
    :param isBase: False if transforming the anterior edge, True if transforming the base of the lung.
    """
    # Transformation matrix = [ -k1y + 1, | [x,
    #                                 1, |  y,
    #                                 1] |  z]
    if isBase:
        # Taper the lower portion of the lung near base (z = 0)
        threshold = 0.3 * halfValue
        offset_vector = [0.0, 0.0, -threshold]
        scale_vector = [0.0, 0.0, sharpeningFactor / (halfValue * 1.75)]
        constant_vector = [1.0, 1.0, 1.0]
        component_indices = [3, 2, 2]
    else:
        # Taper the anterior edge of the lung along y-axis
        threshold = halfValue
        offset_vector = [0.0, 0.4 * threshold, 0.0]
        scale_vector = [0.0, -sharpeningFactor / (halfValue * 1.75), 0.0]
        constant_vector = [0.0, 1.0, 1.0]
        component_indices = [2, 3, 3]

    offset = fm.createFieldConstant(offset_vector)
    origin = fm.createFieldAdd(coordinates, offset)

    scale = fm.createFieldConstant(scale_vector)
    scaleFunction = fm.createFieldMultiply(origin, scale)

    constant = fm.createFieldConstant(constant_vector)
    constantFunction = fm.createFieldAdd(scaleFunction, constant)

    transformation_matrix = fm.createFieldComponent(constantFunction, component_indices)
    taper_coordinates = fm.createFieldMultiply(origin, transformation_matrix)
    translated_coordinates = fm.createFieldSubtract(taper_coordinates, offset)

    if isBase:
        z_component = fm.createFieldComponent(coordinates, [3])
        threshold_field = fm.createFieldConstant(threshold)
        is_below_threshold = fm.createFieldLessThan(z_component, threshold_field)
        result = fm.createFieldIf(is_below_threshold, translated_coordinates, coordinates)
    else:
        y_component = fm.createFieldComponent(coordinates, [2])
        threshold_field = fm.createFieldConstant(threshold)
        is_below_threshold = fm.createFieldLessThan(y_component, threshold_field)
        result = fm.createFieldIf(is_below_threshold, translated_coordinates, coordinates)

    fieldassignment = coordinates.createFieldassignment(result)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def setBaseGroupThreshold(fm, coordinates, halfBreadth, rotateAngle):
    """
    Creates a field to identify lung base elements based on y-coordinate threshold.
    Elements with y-coordinates below 45% of the rotated half-breadth are considered part of the lung base region for
    annotation purposes.
    :param fm: Field module used for creating and managing fields.
    :param coordinates: The coordinate field.
    :param halfBreadth: Half breadth of lung.
    :param rotateAngle: The angle of rotation of horizontal line in radians (90 - oblique fissure angle).
    :return is_above_threshold: True for elements below the y-threshold (base region).
    """
    y_component = fm.createFieldComponent(coordinates, [2])
    y_threshold = 0.45 * halfBreadth * math.cos(rotateAngle)

    y_threshold_field = fm.createFieldConstant(y_threshold)
    is_above_threshold = fm.createFieldLessThan(y_component, y_threshold_field)

    return is_above_threshold
