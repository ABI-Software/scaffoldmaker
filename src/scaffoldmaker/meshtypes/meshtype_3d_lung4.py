"""
Generates a lung scaffold with common hilum but open fissures with lobes built from 1/6 ellipsoid segments.
"""
from cmlibs.maths.vectorops import magnitude
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field

from scaffoldmaker.annotation.annotationgroup import (
    AnnotationGroup, findAnnotationGroupByName, findOrCreateAnnotationGroupForTerm, getAnnotationGroupForTerm)
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.geometry import getEllipsePointAtTrueAngle
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.ellipsoidmesh import EllipsoidMesh, EllipsoidSurfaceD3Mode
from scaffoldmaker.utils.zinc_utils import translate_nodeset_coordinates

import copy
import math


class MeshType_3d_lung4(Scaffold_base):
    """
    Generates a lung scaffold with common hilum but open fissures with lobes built from 1/6 ellipsoid segments.
    """

    @classmethod
    def getName(cls):
        return "3D Lung 4"

    @classmethod
    def getParameterSetNames(cls):
        return [
            "Default",
            "Human 1 Coarse",
            "Human 1 Medium",
            "Human 1 Fine",
            "Ellipsoid"
        ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        options = {}
        useParameterSetName = "Human 1 Coarse" if (parameterSetName == "Default") else parameterSetName
        options["Left lung"] = True
        options["Right lung"] = True
        height = options["Ellipsoid height"] = 1.0
        depth = options["Ellipsoid dorsal-ventral size"] = 0.7
        options["Ellipsoid medial-lateral size"] = 0.4
        options["Left-right lung spacing"] = 0.5
        options["Lower lobe extension"] = 0.3
        options["Refine"] = False
        options["Refine number of elements"] = 4

        if "Medium" in useParameterSetName:
            options["Number of elements lateral"] = 4
            options["Number of elements lower extension"] = 3
            options["Number of elements oblique"] = 8
            options["Number of transition elements"] = 1
        elif "Fine" in useParameterSetName:
            options["Number of elements lateral"] = 6
            options["Number of elements lower extension"] = 4
            options["Number of elements oblique"] = 12
            options["Number of transition elements"] = 1
        else:
            options["Number of elements lateral"] = 4
            options["Number of elements lower extension"] = 2
            options["Number of elements oblique"] = 6
            options["Number of transition elements"] = 1

        if "Human" in useParameterSetName:
            options["Ventral edge sharpness factor"] = 0.8
            options["Medial curvature"] = 3.0
            options["Medial curvature bias"] = 1.0
        else:
            options["Ventral edge sharpness factor"] = 0.0
            options["Medial curvature"] = 0.0
            options["Medial curvature bias"] = 0.0

        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Left lung",
            "Right lung",
            # "Number of left lung lobes",
            "Number of elements lateral",
            "Number of elements lower extension",
            "Number of elements oblique",
            "Number of transition elements",
            "Ellipsoid height",
            "Ellipsoid dorsal-ventral size",
            "Ellipsoid medial-lateral size",
            "Left-right lung spacing",
            "Lower lobe extension",
            "Ventral edge sharpness factor",
            "Medial curvature",
            "Medial curvature bias",
            "Refine",
            "Refine number of elements"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependent_changes = False

        max_transition_count = None
        for key in [
            "Number of elements lower extension"
        ]:
            if options[key] < 1:
                options[key] = 1
        for key in [
            "Number of elements lateral",
            "Number of elements oblique"
        ]:
            min_elements_count = 4 if (key == "Number of elements lateral") else 6
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
            dependent_changes = True

        for key in [
            "Ellipsoid height",
            "Ellipsoid dorsal-ventral size",
            "Ellipsoid medial-lateral size",
            "Lower lobe extension"
        ]:
            if options[key] <= 0.0:
                options[key] = 1.0
        depth = options["Ellipsoid dorsal-ventral size"]
        height = options["Ellipsoid height"]
        max_extension = 0.99 * magnitude(getEllipsePointAtTrueAngle(depth / 2.0, height / 2.0, math.pi / 3.0))
        if options["Lower lobe extension"] > max_extension:
            options["Lower lobe extension"] = max_extension

        if options["Left-right lung spacing"] < 0.0:
            options["Left-right lung spacing"] = 0.0

        for dimension in [
            "Ventral edge sharpness factor",
            "Medial curvature bias"
        ]:
            if options[dimension] < 0.0:
                options[dimension] = 0.0
            elif options[dimension] > 1.0:
                options[dimension] = 1.0

        if options['Refine number of elements'] < 1:
            options['Refine number of elements'] = 1

        return dependent_changes

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        is_left_lung = options["Left lung"]
        is_right_lung = options["Right lung"]

        elements_count_lateral = options["Number of elements lateral"]
        elements_count_lower_extension = options["Number of elements lower extension"]
        elements_count_oblique = options["Number of elements oblique"]
        elements_count_transition = options["Number of transition elements"]
        lung_spacing = options["Left-right lung spacing"] * 0.5
        lower_lobe_extension_height = options["Lower lobe extension"]
        ventral_sharpness_factor = options["Ventral edge sharpness factor"]
        ellipsoid_height = options["Ellipsoid height"]
        ellipsoid_breadth = options["Ellipsoid dorsal-ventral size"]
        ellipsoid_depth = options["Ellipsoid medial-lateral size"]
        medial_curvature = options["Medial curvature"]
        medial_curvature_bias = options["Medial curvature bias"]

        fieldmodule = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fieldmodule)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        # annotation groups & nodeset groups
        lungGroup = AnnotationGroup(region, get_lung_term("lung"))

        leftLungGroup = AnnotationGroup(region, get_lung_term("left lung"))
        rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))

        leftLateralLungGroup = AnnotationGroup(region, ("lateral left lung", ""))
        rightLateralLungGroup = AnnotationGroup(region, ("lateral right lung", ""))

        leftMedialLungGroup = AnnotationGroup(region, ("medial left lung", ""))
        rightMedialLungGroup = AnnotationGroup(region, ("medial right lung", ""))

        leftAnteriorLungGroup = AnnotationGroup(region, ("anterior left lung", ""))
        rightAnteriorLungGroup = AnnotationGroup(region, ("anterior right lung", ""))

        lowerRightLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of right lung"))
        upperRightLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of right lung"))
        middleRightLungGroup = AnnotationGroup(region, get_lung_term("middle lobe of right lung"))

        annotation_groups = [lungGroup, leftLungGroup, rightLungGroup,
                            leftAnteriorLungGroup, leftLateralLungGroup, leftMedialLungGroup,
                            rightAnteriorLungGroup, rightLateralLungGroup, rightMedialLungGroup,
                            lowerRightLungGroup, middleRightLungGroup, upperRightLungGroup]

        lowerLeftLungGroup = AnnotationGroup(region, get_lung_term("lower lobe of left lung"))
        upperLeftLungGroup = AnnotationGroup(region, get_lung_term("upper lobe of left lung"))
        annotation_groups += [lowerLeftLungGroup, upperLeftLungGroup]

        box_group = AnnotationGroup(region, ("box", ""))
        transition_group = AnnotationGroup(region, ("transition", ""))
        annotation_groups += [box_group, transition_group]

        # Nodeset group
        leftLungNodesetGroup = leftLungGroup.getNodesetGroup(nodes)
        rightLungNodesetGroup = rightLungGroup.getNodesetGroup(nodes)

        halfDepth = ellipsoid_depth * 0.5
        halfBreadth = ellipsoid_breadth * 0.5
        halfHeight = ellipsoid_height * 0.5

        left_lung, right_lung = 0, 1
        lungs = [lung for show, lung in [(is_left_lung, left_lung), (is_right_lung, right_lung)] if show]
        nodeIdentifier, elementIdentifier = 1, 1

        for lung in lungs:

            if lung == left_lung:
                lower_octant_group_lists = []
                middle_octant_group_lists = None
                upper_octant_group_lists = []
                for octant in range(8):
                    octant_group_list = [group.getGroup() for group in [lungGroup, leftLungGroup, lowerLeftLungGroup] +
                                         [leftMedialLungGroup if (octant & 1) else leftLateralLungGroup]]
                    lower_octant_group_lists.append(octant_group_list)
                    octant_group_list = [group.getGroup() for group in [lungGroup, leftLungGroup, upperLeftLungGroup] +
                                         [leftMedialLungGroup if (octant & 1) else leftLateralLungGroup]]
                    if octant & 2:
                        octant_group_list.append(leftAnteriorLungGroup.getGroup())
                    upper_octant_group_lists.append(octant_group_list)
            else:
                lower_octant_group_lists = []
                middle_octant_group_lists = []
                upper_octant_group_lists = []
                for octant in range(8):
                    octant_group_list = [group.getGroup() for group in [lungGroup, rightLungGroup, lowerRightLungGroup] +
                                         [rightLateralLungGroup if (octant & 1) else rightMedialLungGroup]]
                    lower_octant_group_lists.append(octant_group_list)
                    octant_group_list = [group.getGroup() for group in [lungGroup, rightLungGroup, middleRightLungGroup] +
                                         [rightLateralLungGroup if (octant & 1) else rightMedialLungGroup]]
                    if octant & 2:
                        octant_group_list.append(rightAnteriorLungGroup.getGroup())
                    middle_octant_group_lists.append(octant_group_list)
                    octant_group_list = [group.getGroup() for group in [lungGroup, rightLungGroup, upperRightLungGroup] +
                                         [rightLateralLungGroup if (octant & 1) else rightMedialLungGroup]]
                    if octant & 2:
                        octant_group_list.append(rightAnteriorLungGroup.getGroup())
                    upper_octant_group_lists.append(octant_group_list)

            elementCounts = [elements_count_lateral, elements_count_oblique, elements_count_oblique]
            lower_ellipsoid = EllipsoidMesh(halfDepth, halfBreadth, halfHeight, elementCounts, elements_count_transition)
            lower_ellipsoid.set_surface_d3_mode(EllipsoidSurfaceD3Mode.SURFACE_NORMAL_PLANE_PROJECTION)
            lower_ellipsoid.set_box_transition_groups(box_group.getGroup(), transition_group.getGroup())
            upper_ellipsoid = EllipsoidMesh(halfDepth, halfBreadth, halfHeight, elementCounts, elements_count_transition)
            upper_ellipsoid.set_surface_d3_mode(EllipsoidSurfaceD3Mode.SURFACE_NORMAL_PLANE_PROJECTION)
            upper_ellipsoid.set_box_transition_groups(box_group.getGroup(), transition_group.getGroup())
            if lung == right_lung:
                middle_ellipsoid = EllipsoidMesh(halfDepth, halfBreadth, halfHeight, elementCounts, elements_count_transition)
                middle_ellipsoid.set_surface_d3_mode(EllipsoidSurfaceD3Mode.SURFACE_NORMAL_PLANE_PROJECTION)
                middle_ellipsoid.set_box_transition_groups(box_group.getGroup(), transition_group.getGroup())
            else:
                middle_ellipsoid = upper_ellipsoid
            pi__3 = math.pi / 3.0
            normal_face_factor = 1.0
            half_counts = [count // 2 for count in elementCounts]
            octant1 = middle_ellipsoid.build_octant(half_counts, -pi__3, 0.0, normal_face_factor=normal_face_factor)
            middle_ellipsoid.merge_octant(octant1, quadrant=3)
            if lung == right_lung:
                middle_ellipsoid.copy_to_negative_axis1()

            # save hilum coordinates for all other lobes
            hilum_x = []
            ox = octant1.get_parameters()
            box_count1 = octant1.get_box_counts()[0]
            for n1 in range(elementCounts[0] + 1):
                mirror_x = n1 < half_counts[0]
                o1 = abs(n1 - half_counts[0])
                parameters = ox[0][0][o1]
                obox = o1 <= box_count1
                parameters = [
                    copy.copy(parameters[0]),
                    copy.copy(parameters[3 if obox else 2]),
                    [-d for d in parameters[2 if obox else 1]],
                    copy.copy(parameters[1 if obox else 3])
                ]
                if mirror_x:
                    for i in range(3):
                        parameters[i][0] = -parameters[i][0]
                hilum_x.append(parameters)

            octant2 = upper_ellipsoid.build_octant(half_counts, 0.0, pi__3, normal_face_factor=normal_face_factor)
            upper_ellipsoid.merge_octant(octant2, quadrant=0)
            octant3 = upper_ellipsoid.build_octant(half_counts, pi__3, 2.0 * pi__3, normal_face_factor=normal_face_factor)
            upper_ellipsoid.merge_octant(octant3, quadrant=1)
            upper_ellipsoid.copy_to_negative_axis1()
    
            lower_lobe_extension = 0.6 * halfHeight / math.cos(math.pi / 6.0)
            octant4 = lower_ellipsoid.build_octant(half_counts, 2.0 * pi__3, math.pi,
                                                   lower_lobe_extension_height, elements_count_lower_extension,
                                                   normal_face_factor=normal_face_factor)
            # merge into separate lower ellipsoid to have space for extension elements
            lower_ellipsoid_mesh = EllipsoidMesh(
                halfDepth, halfBreadth, halfHeight,
                [elementCounts[0], elementCounts[1], elementCounts[2] + 2 * elements_count_lower_extension],
                elements_count_transition)
            lower_ellipsoid_mesh.set_surface_d3_mode(EllipsoidSurfaceD3Mode.SURFACE_NORMAL_PLANE_PROJECTION)
            lower_ellipsoid_mesh.set_box_transition_groups(box_group.getGroup(), transition_group.getGroup())
            lower_ellipsoid_mesh.merge_octant(octant4, quadrant=1)
            lower_ellipsoid_mesh.copy_to_negative_axis1()
    
            node_layout_manager = lower_ellipsoid.get_node_layout_manager()
            node_layout_permuted = node_layout_manager.getNodeLayoutRegularPermuted(d3Defined=True)
            for n1 in range(elementCounts[0] + 1):
                lower_ellipsoid_mesh.set_node_parameters(
                    n1, half_counts[1], elementCounts[2] + 2 * elements_count_lower_extension - half_counts[2], hilum_x[n1],
                    node_layout=node_layout_permuted)
            lower_ellipsoid_mesh.set_octant_group_lists(lower_octant_group_lists)
            nodeIdentifier, elementIdentifier = lower_ellipsoid_mesh.generate_mesh(
                fieldmodule, coordinates, nodeIdentifier, elementIdentifier)

            for ellipsoid in [middle_ellipsoid, upper_ellipsoid] if (lung == right_lung) else [upper_ellipsoid]:
                node_layout_manager = ellipsoid.get_node_layout_manager()
                node_layout_6way = node_layout_manager.getNodeLayout6Way12(d3Defined=True)
                for n1 in range(elementCounts[0] + 1):
                    nid = (lower_ellipsoid_mesh.get_node_identifier(
                        n1, half_counts[1], elementCounts[2] + 2 * elements_count_lower_extension - half_counts[2])
                           if (((lung == left_lung) and (n1 >= half_counts[0])) or
                               ((lung == right_lung) and (n1 <= half_counts[0]))) else None)
                    ellipsoid.set_node_parameters(n1, half_counts[1], half_counts[2], hilum_x[n1],
                                                        nid, node_layout=node_layout_6way)
                ellipsoid.set_octant_group_lists(
                    middle_octant_group_lists if ((ellipsoid == middle_ellipsoid) and
                                                  (ellipsoid != upper_ellipsoid)) else upper_octant_group_lists)
                nodeIdentifier, elementIdentifier = ellipsoid.generate_mesh(
                    fieldmodule, coordinates, nodeIdentifier, elementIdentifier)

        for lung in lungs:
            is_left = lung == left_lung
            lungNodeset = leftLungNodesetGroup if is_left else rightLungNodesetGroup
            spacing = -lung_spacing if is_left else lung_spacing
            zOffset = -0.5 * ellipsoid_height
            lungMedialCurvature = -medial_curvature if is_left else medial_curvature

            if ventral_sharpness_factor != 0.0:
                taperLungEdge(ventral_sharpness_factor, fieldmodule, coordinates, lungNodeset, halfBreadth)

            dorsalVentralXi = getDorsalVentralXiField(fieldmodule, coordinates, halfBreadth)
            if lungMedialCurvature != 0:
                bendLungMeshAroundZAxis(lungMedialCurvature, fieldmodule, coordinates, lungNodeset,
                                        stationaryPointXY=[0.0, 0.0],
                                        bias=medial_curvature_bias,
                                        dorsalVentralXi=dorsalVentralXi)

            translate_nodeset_coordinates(lungNodeset, coordinates, [spacing, 0, -zOffset])

        return annotation_groups, None


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
    def defineFaceAnnotations(cls, region, options, annotation_groups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotation_groups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        fm = region.getFieldmodule()
        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)

        is_exterior = fm.createFieldIsExterior()
        is_face_xi1_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0)
        is_face_xi1_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1)
        is_face_xi2_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0)
        is_face_xi2_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1)
        is_face_xi3_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)
        is_face_xi3_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, is_face_xi3_1)

        box_group = findAnnotationGroupByName(annotation_groups, "box")
        is_box = box_group.getGroup()
        transition_group = findAnnotationGroupByName(annotation_groups, "transition")
        is_trans = transition_group.getGroup()
        is_on_ellipsoid = fm.createFieldAnd(fm.createFieldAnd(is_exterior_face_xi3_1, is_trans),
                                            fm.createFieldNot(is_box))

        is_face_box20_trans10 = fm.createFieldOr(
            fm.createFieldAnd(is_box, is_face_xi2_0), fm.createFieldAnd(is_trans, is_face_xi1_0))
        is_face_xi1_0_or_xi1_1_or_xi2_0 = fm.createFieldOr(
            fm.createFieldOr(is_face_xi1_0, is_face_xi1_1), is_face_xi2_0)
        is_face_xi1_0_or_xi1_1_or_xi2_1 = fm.createFieldOr(
            fm.createFieldOr(is_face_xi1_0, is_face_xi1_1), is_face_xi2_1)
        is_face_box21_trans11 = fm.createFieldOr(
            fm.createFieldAnd(is_box, is_face_xi2_1), fm.createFieldAnd(is_trans, is_face_xi1_1))
        is_face_box30_trans20 = fm.createFieldOr(
            fm.createFieldAnd(is_box, is_face_xi3_0), fm.createFieldAnd(is_trans, is_face_xi2_0))
        is_face_box31_trans21 = fm.createFieldOr(
            fm.createFieldAnd(is_box, is_face_xi3_1), fm.createFieldAnd(is_trans, is_face_xi2_1))

        is_lower_left = getAnnotationGroupForTerm(annotation_groups, get_lung_term("lower lobe of left lung")).getGroup()
        is_upper_left = getAnnotationGroupForTerm(annotation_groups, get_lung_term("upper lobe of left lung")).getGroup()
        is_anterior_left = findAnnotationGroupByName(annotation_groups, "anterior left lung").getGroup()
        is_lateral_left = findAnnotationGroupByName(annotation_groups, "lateral left lung").getGroup()
        is_medial_left = findAnnotationGroupByName(annotation_groups, "medial left lung").getGroup()

        is_lower_right = getAnnotationGroupForTerm(annotation_groups, get_lung_term("lower lobe of right lung")).getGroup()
        is_middle_right = getAnnotationGroupForTerm(annotation_groups, get_lung_term("middle lobe of right lung")).getGroup()
        is_upper_right = getAnnotationGroupForTerm(annotation_groups, get_lung_term("upper lobe of right lung")).getGroup()
        is_anterior_right = findAnnotationGroupByName(annotation_groups, "anterior right lung").getGroup()
        is_lateral_right = findAnnotationGroupByName(annotation_groups, "lateral right lung").getGroup()
        is_medial_right = findAnnotationGroupByName(annotation_groups, "medial right lung").getGroup()

        face_term_conditionals_map = {
            "base of lower lobe of left lung surface": (is_lower_left, is_exterior, is_face_box30_trans20),
            "base of lower lobe of right lung surface": (is_lower_right, is_exterior, is_face_box30_trans20),
            "base of middle lobe of right lung surface": (is_middle_right, is_exterior, is_face_xi1_0_or_xi1_1_or_xi2_0),
            "base of upper lobe of left lung surface": (is_upper_left, is_exterior, is_face_xi1_0_or_xi1_1_or_xi2_0, is_anterior_left),
            "horizontal fissure of middle lobe of right lung": (is_middle_right, is_exterior, is_face_box31_trans21),
            "horizontal fissure of upper lobe of right lung": (is_upper_right, is_exterior, is_face_box30_trans20, is_anterior_right),
            "lateral surface of left lung": (is_lateral_left, is_on_ellipsoid),
            "lateral surface of lower lobe of left lung": (is_lower_left, is_on_ellipsoid, is_lateral_left),
            "lateral surface of lower lobe of right lung": (is_lower_right, is_on_ellipsoid, is_lateral_right),
            "lateral surface of middle lobe of right lung": (is_middle_right, is_on_ellipsoid, is_lateral_right),
            "lateral surface of right lung": (is_lateral_right, is_on_ellipsoid),
            "lateral surface of upper lobe of left lung": (is_upper_left, is_on_ellipsoid, is_lateral_left),
            "lateral surface of upper lobe of right lung": (is_upper_right, is_on_ellipsoid, is_lateral_right),
            "lower lobe of left lung surface": (is_lower_left, is_exterior),
            "lower lobe of right lung surface": (is_lower_right, is_exterior),
            "middle lobe of right lung surface": (is_middle_right, is_exterior),
            "upper lobe of left lung surface": (is_upper_left, is_exterior),
            "upper lobe of right lung surface": (is_upper_right, is_exterior),
            "medial surface of left lung": (is_medial_left, is_on_ellipsoid),
            "medial surface of lower lobe of left lung": (is_lower_left, is_on_ellipsoid, is_medial_left),
            "medial surface of lower lobe of right lung": (is_lower_right, is_on_ellipsoid, is_medial_right),
            "medial surface of middle lobe of right lung": (is_middle_right, is_on_ellipsoid, is_medial_right),
            "medial surface of right lung": (is_medial_right, is_on_ellipsoid),
            "medial surface of upper lobe of left lung": (is_upper_left, is_on_ellipsoid, is_medial_left),
            "medial surface of upper lobe of right lung": (is_upper_right, is_on_ellipsoid, is_medial_right),
            "oblique fissure of lower lobe of left lung": (is_lower_left, is_exterior, is_face_xi1_0_or_xi1_1_or_xi2_1),
            "oblique fissure of lower lobe of right lung": (is_lower_right, is_exterior, is_face_xi1_0_or_xi1_1_or_xi2_1),
            "oblique fissure of middle lobe of right lung": (is_middle_right, is_exterior, is_face_xi1_0_or_xi1_1_or_xi2_0),
            "oblique fissure of upper lobe of left lung": (is_upper_left, is_exterior, fm.createFieldOr(is_face_xi1_0_or_xi1_1_or_xi2_0, is_face_box30_trans20)),
            "oblique fissure of upper lobe of right lung": (is_upper_right, is_exterior, fm.createFieldAnd(is_face_box30_trans20, fm.createFieldNot(is_anterior_right))),
        }
        is_face_conditional = {}
        for face_term, conditionals in face_term_conditionals_map.items():
            annotation_group = findOrCreateAnnotationGroupForTerm(annotation_groups, region, get_lung_term(face_term))
            group = annotation_group.getGroup()
            conditional = conditionals[0]
            for add_conditional in conditionals[1:]:
                conditional = fm.createFieldAnd(conditional, add_conditional)
            annotation_group.getMeshGroup(mesh2d).addElementsConditional(conditional)
            print(conditional.isValid(), "Term:", face_term, "sizes", mesh2d.getSize(), group.getMeshGroup(mesh2d).isValid(), group.getMeshGroup(mesh2d).getSize(), group.getMeshGroup(mesh1d).getSize())
            annotation_group.addSubelements()
            is_face_conditional[face_term] = group

        line_term_conditionals_map = {
            "anterior edge of middle lobe of right lung":
                (is_middle_right, is_on_ellipsoid, is_medial_right, is_lateral_right),
            "antero-posterior edge of upper lobe of left lung":
                (is_upper_left, is_on_ellipsoid, is_medial_left, is_lateral_left),
            "antero-posterior edge of upper lobe of right lung":
                (is_upper_right, is_on_ellipsoid, is_medial_right, is_lateral_right),
            "base edge of oblique fissure of lower lobe of left lung": (
                is_face_conditional["base of lower lobe of left lung surface"], is_face_conditional["oblique fissure of lower lobe of left lung"]),
            "base edge of oblique fissure of lower lobe of right lung": (
                is_face_conditional["base of lower lobe of right lung surface"], is_face_conditional["oblique fissure of lower lobe of right lung"]),
            "lateral edge of base of lower lobe of left lung": (
                is_face_conditional["base of lower lobe of left lung surface"], is_lateral_left, is_on_ellipsoid),
            "lateral edge of base of lower lobe of right lung": (
                is_face_conditional["base of lower lobe of right lung surface"], is_lateral_right, is_on_ellipsoid),
            "lateral edge of base of middle lobe of right lung": (
                is_face_conditional["base of middle lobe of right lung surface"], is_lateral_right, is_on_ellipsoid),
            "lateral edge of base of upper lobe of left lung": (
                is_face_conditional["base of upper lobe of left lung surface"], is_lateral_left, is_on_ellipsoid),
            "lateral edge of horizontal fissure of middle lobe of right lung": (
                is_face_conditional["horizontal fissure of middle lobe of right lung"], is_lateral_right, is_on_ellipsoid),
            "lateral edge of horizontal fissure of upper lobe of right lung": (
                is_face_conditional["horizontal fissure of upper lobe of right lung"], is_lateral_right, is_on_ellipsoid),
            "lateral edge of oblique fissure of lower lobe of left lung": (
                is_face_conditional["oblique fissure of lower lobe of left lung"], is_lateral_left, is_on_ellipsoid),
            "lateral edge of oblique fissure of lower lobe of right lung": (
                is_face_conditional["oblique fissure of lower lobe of right lung"], is_lateral_right, is_on_ellipsoid),
            "lateral edge of oblique fissure of middle lobe of right lung": (
                is_face_conditional["oblique fissure of middle lobe of right lung"], is_lateral_right, is_on_ellipsoid),
            "lateral edge of oblique fissure of upper lobe of left lung": (
                is_face_conditional["oblique fissure of upper lobe of left lung"], is_lateral_left, is_on_ellipsoid),
            "lateral edge of oblique fissure of upper lobe of right lung": (
                is_face_conditional["oblique fissure of upper lobe of right lung"], is_lateral_right, is_on_ellipsoid),
            "medial edge of base of lower lobe of left lung": (
                is_face_conditional["base of lower lobe of left lung surface"], is_medial_left, is_on_ellipsoid),
            "medial edge of base of lower lobe of right lung": (
                is_face_conditional["base of lower lobe of right lung surface"], is_medial_right, is_on_ellipsoid),
            "medial edge of base of middle lobe of right lung": (
                is_face_conditional["base of middle lobe of right lung surface"], is_medial_right, is_on_ellipsoid),
            "medial edge of base of upper lobe of left lung": (
                is_face_conditional["base of upper lobe of left lung surface"], is_medial_left, is_on_ellipsoid),
            "medial edge of horizontal fissure of middle lobe of right lung": (
                is_face_conditional["horizontal fissure of middle lobe of right lung"], is_medial_right, is_on_ellipsoid),
            "medial edge of horizontal fissure of upper lobe of right lung": (
                is_face_conditional["horizontal fissure of upper lobe of right lung"], is_medial_right, is_on_ellipsoid),
            "medial edge of oblique fissure of lower lobe of left lung": (
                is_face_conditional["oblique fissure of lower lobe of left lung"], is_medial_left, is_on_ellipsoid),
            "medial edge of oblique fissure of lower lobe of right lung": (
                is_face_conditional["oblique fissure of lower lobe of right lung"], is_medial_right, is_on_ellipsoid),
            "medial edge of oblique fissure of middle lobe of right lung": (
                is_face_conditional["oblique fissure of middle lobe of right lung"], is_medial_right, is_on_ellipsoid),
            "medial edge of oblique fissure of upper lobe of left lung": (
                is_face_conditional["oblique fissure of upper lobe of left lung"], is_medial_left, is_on_ellipsoid),
            "medial edge of oblique fissure of upper lobe of right lung": (
                is_face_conditional["oblique fissure of upper lobe of right lung"], is_medial_right, is_on_ellipsoid),
            "posterior edge of lower lobe of left lung":
                (is_lower_left, is_lateral_left, is_medial_left, is_on_ellipsoid),
            "posterior edge of lower lobe of right lung":
                (is_lower_right, is_lateral_right, is_medial_right, is_on_ellipsoid),

        }
        for line_term, conditionals in line_term_conditionals_map.items():
            annotation_group = findOrCreateAnnotationGroupForTerm(annotation_groups, region, (line_term, ""))
            conditional = conditionals[0]
            for add_conditional in conditionals[1:]:
                conditional = fm.createFieldAnd(conditional, add_conditional)
            annotation_group.getMeshGroup(mesh1d).addElementsConditional(conditional)

        # remove temporary annotation groups
        for group_name in [
            "anterior left lung",
            "anterior right lung",
            "box",
            "lateral left lung",
            "lateral right lung",
            "medial left lung",
            "medial right lung",
            "transition"
        ]:
            annotation_group = findAnnotationGroupByName(annotation_groups, group_name)
            annotation_groups.remove(annotation_group)


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
    :param coordinates: The coordinate field. The anterior edge is towards the +y-axis, and the base is towards the
    -z-axis.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param halfValue: Half value of lung breadth/height depending on isBase.
    :param isBase: False if transforming the anterior edge, True if transforming the base of the lung.
    """
    x = fm.createFieldComponent(coordinates, 1)
    y = fm.createFieldComponent(coordinates, 2)
    z = fm.createFieldComponent(coordinates, 3)

    coord_value = z if isBase else y
    start_value = 0.5 * halfValue if isBase else -0.5 * halfValue
    end_value = -1.1 * halfValue if isBase else 1.1 * halfValue

    start_value_field = fm.createFieldConstant(start_value)
    end_value_field = fm.createFieldConstant(end_value)

    xi = (coord_value - start_value_field) / fm.createFieldConstant(end_value - start_value)
    xi__2 = xi * xi
    one = fm.createFieldConstant(1.0)
    x_scale = one - fm.createFieldConstant(sharpeningFactor) * xi__2
    if isBase:
        new_x = fm.createFieldIf(fm.createFieldLessThan(coord_value, start_value_field), x * x_scale, x)
    else:
        new_x = fm.createFieldIf(fm.createFieldGreaterThan(coord_value, start_value_field), x * x_scale, x)
    new_coordinates = fm.createFieldConcatenate([new_x, y, z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
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
