"""
Generates a lung scaffold by deforming a hemisphere.
"""
from cmlibs.maths.vectorops import mult, sub
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine, \
    sampleCubicHermiteCurves, DerivativeScalingMode
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
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
            "Human 1",
            "Ellipsoid",
            "Teardrop"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        useParameterSetName = "Human 1" if (parameterSetName == "Default") else parameterSetName
        isHuman = "Human" in useParameterSetName
        if parameterSetName in ["Human 1", "Default"]:
            options = {
                "Show left lung": True,
                "Show right lung": True,
                "Open fissures": False,
                "Number of left lung lobes": 2,
                "Number of elements lateral": 6,
                "Number of elements normal": 8,
                "Number of elements oblique": 8,
                "Number of elements shell": 0,
                "Left-right lung spacing": 0.4,
                "Apex edge sharpness factor": 0.4,
                "Ventral edge sharpness factor": 0.95,
                "Left-right apex medial shear displacement": 0.2,
                "Left-right apex ventral shear displacement": -0.2,
                # "Left-right base shear displacement": 0.0,
                "Diaphragm angle degrees": 50.0,
                "Diaphragm proportion": 0.3,
                "Disc breadth": 0.8,
                "Disc height": 1.3,
                "Disc depth": 0.3,
                "Impression breadth proportion": 0.9,
                "Impression height proportion": 0.8,
                "Impression depth proportion": 1.0,
                "Lateral shear rate": -1.5,
                "Left oblique slope degrees": 45.0,
                "Right oblique slope degrees": 45.0,
                "Medial curvature": 1.0,
                "Medial curvature bias": 0.5,
                "Medial rotation about x-axis degrees": 5.0,
                "Ventral-medial rotation degrees": 5.0,
                "Use sizing function": True,
                "Scale factor": 0.7,
                "Refine": False,
                "Refine number of elements": 4,
            }
        elif parameterSetName == "Ellipsoid":
            options = {
                "Show left lung": False,
                "Show right lung": False,
                "Open fissures": False,
                "Number of left lung lobes": 2,
                "Number of elements lateral": 4,
                "Number of elements normal": 8,
                "Number of elements oblique": 8,
                "Number of elements shell": 0,
                "Left-right lung spacing": 1.0,
                "Apex edge sharpness factor": 0.0,
                "Ventral edge sharpness factor": 0.0,
                "Left-right apex medial shear displacement": 0.0,
                "Left-right apex ventral shear displacement": 0.0,
                # "Left-right base shear displacement": 0.0,
                "Diaphragm angle degrees": 0.0,
                "Diaphragm proportion": 0.0,
                "Disc breadth": 0.8,
                "Disc height": 1.1,
                "Disc depth": 0.25,
                "Impression breadth proportion": 0.0,
                "Impression height proportion": 0.0,
                "Impression depth proportion": 0.0,
                "Lateral shear rate": 0.0,
                "Left oblique slope degrees": 45.0,
                "Right oblique slope degrees": 45.0,
                "Medial curvature": 0.0,
                "Medial curvature bias": 0.0,
                "Medial rotation about x-axis degrees": 0.0,
                "Ventral-medial rotation degrees": 0.0,
                "Use sizing function": False,
                "Scale factor": 1.0,
                "Refine": False,
                "Refine number of elements": 4,
            }
        elif parameterSetName == "Teardrop":
            options = {
                "Show left lung": False,
                "Show right lung": False,
                "Open fissures": False,
                "Number of left lung lobes": 2,
                "Number of elements lateral": 4,
                "Number of elements normal": 4,
                "Number of elements oblique": 4,
                "Number of elements shell": 0,
                "Left-right lung spacing": 1.0,
                "Apex edge sharpness factor": 0.0,
                "Ventral edge sharpness factor": 0.95,
                "Left-right apex medial shear displacement": 0.0,
                "Left-right apex ventral shear displacement": 0.0,
                # "Left-right base shear displacement": 0.0,
                "Diaphragm angle degrees": 0.0,
                "Diaphragm proportion": 0.0,
                "Disc breadth": 0.8,
                "Disc height": 1.1,
                "Disc depth": 0.25,
                "Impression breadth proportion": 0.0,
                "Impression height proportion": 0.0,
                "Impression depth proportion": 0.0,
                "Lateral shear rate": 0.0,
                "Left oblique slope degrees": 45.0,
                "Right oblique slope degrees": 45.0,
                "Medial curvature": 0.0,
                "Medial curvature bias": 0.0,
                "Medial rotation about x-axis degrees": 0.0,
                "Ventral-medial rotation degrees": 0.0,
                "Use sizing function": False,
                "Scale factor": 1.0,
                "Refine": False,
                "Refine number of elements": 4,
            }
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Show left lung",
            "Show right lung",
            "Open fissures",
            "Number of left lung lobes",
            "Number of elements lateral",
            "Number of elements normal",
            "Number of elements oblique",
            "Number of elements shell",
            "Left-right lung spacing",
            "Apex edge sharpness factor",
            "Ventral edge sharpness factor",
            "Left-right apex medial shear displacement",
            "Left-right apex ventral shear displacement",
            # "Left-right base shear displacement",
            "Diaphragm angle degrees",
            "Diaphragm proportion",
            "Disc breadth",
            "Disc height",
            "Disc depth",
            "Impression breadth proportion",
            "Impression height proportion",
            "Impression depth proportion",
            "Lateral shear rate",
            "Left oblique slope degrees",
            "Right oblique slope degrees",
            "Medial curvature",
            "Medial curvature bias",
            "Medial rotation about x-axis degrees",
            "Ventral-medial rotation degrees",
            "Use sizing function",
            "Scale factor",
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

        if options["Number of elements lateral"] < 2:
            options["Number of elements lateral"] = 2
        for key in [
            "Number of elements normal",
            "Number of elements oblique"]:
            if options[key] < 4:
                options[key] = 4
            elif options[key] % 2:
                options[key] += 1

        maxShellElements = min(options["Number of elements lateral"] // 2,
                               options["Number of elements normal"] // 2,
                               options["Number of elements oblique"] // 2) - 2
        if options["Number of elements shell"] < 0:
            options["Number of elements shell"] = 0
        elif options["Number of elements shell"] > maxShellElements:
            options["Number of elements shell"] = maxShellElements
            dependentChanges = True

        for dimension in [
            "Disc breadth",
            "Disc height",
            "Disc depth"]:
            if options[dimension] <= 0.0:
                options[dimension] = 1.0
        for dimension in [
            "Diaphragm proportion",
            "Impression depth proportion",
            "Impression height proportion",
            "Impression breadth proportion"]:
            if options[dimension] <= -1.0:
                options[dimension] = 0.0
            elif options[dimension] > 1.0:
                options[dimension] = 1.0

        if options["Left-right lung spacing"] < 0.0:
            options["Left-right lung spacing"] = 0.0

        for dimension in [
            "Apex edge sharpness factor",
            "Ventral edge sharpness factor",
            "Left-right apex medial shear displacement",
            "Medial curvature bias"
        ]:
            if options[dimension] < 0.0:
                options[dimension] = 0.0
            elif options[dimension] > 1.0:
                options[dimension] = 1.0

        for dimension in [
            "Left-right apex ventral shear displacement",
        ]:
            if options[dimension] < -1.0:
                options[dimension] = -1.0
            elif options[dimension] > 1.0:
                options[dimension] = 1.0

        for angle in [
            "Medial rotation about x-axis degrees",
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
        showLeftLung = options["Show left lung"]
        showRightLung = options["Show right lung"]
        isOpenfissure = options["Open fissures"]
        hasAccessoryLobe = False
        numberOfLeftLung = options["Number of left lung lobes"]
        useSerendipity = True
        useSizingFunction = options["Use sizing function"]

        elementsCountLateral = options["Number of elements lateral"]
        elementsCountNormal = options["Number of elements normal"]
        elementsCountOblique = options["Number of elements oblique"]
        elementsCountShell = options["Number of elements shell"]
        elementsCountTransition = 1
        lungSpacing = options["Left-right lung spacing"] * 0.5
        apexSharpFactor = options["Apex edge sharpness factor"]
        edgeSharpFactor = options["Ventral edge sharpness factor"]
        leftApexMedialDisplacement = options["Left-right apex medial shear displacement"]
        forwardLeftRightApex = options["Left-right apex ventral shear displacement"]
        # forwardLeftRightBase = options["Left-right base shear displacement"]
        diaphragm_angle_radians = math.radians(options["Diaphragm angle degrees"])
        diaphragm_proportion = options["Diaphragm proportion"]
        disc_breadth = options["Disc breadth"]
        disc_height = options["Disc height"]
        disc_depth = options["Disc depth"]
        impression_breadth_proportion = options["Impression breadth proportion"]
        impression_height_proportion = options["Impression height proportion"]
        impression_depth_proportion = options["Impression depth proportion"]
        lateral_shear_rate = options["Lateral shear rate"]
        left_oblique_slope_radians = math.radians(options["Left oblique slope degrees"])
        right_oblique_slope_radians = math.radians(options["Right oblique slope degrees"])
        leftLungMedialCurvature = options["Medial curvature"]
        lungMedialCurvatureBias = options["Medial curvature bias"]
        rotateLeftLungX = options["Medial rotation about x-axis degrees"]
        rotateLeftLungZ = options["Ventral-medial rotation degrees"]

        rangeOfRequiredElements = [
            [0, elementsCountLateral],
            [0, elementsCountOblique],
            [0, elementsCountNormal]
        ]

        boxDerivatives = [1, 2, 3]
        sphereBoxDerivatives = [-boxDerivatives[0], boxDerivatives[1], boxDerivatives[2]]
        sphere_shape = SphereShape.SPHERE_SHAPE_FULL

        fieldmodule = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fieldmodule)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        mesh = fieldmodule.findMeshByDimension(3)

        boxGroup = AnnotationGroup(region, ("box group", "None"))
        boxMeshGroup = boxGroup.getMeshGroup(mesh)
        boxNodesetGroup = boxGroup.getNodesetGroup(nodes)

        transitionGroup = AnnotationGroup(region, ("transition group", "None"))
        transitionMeshGroup = transitionGroup.getMeshGroup(mesh)
        transitionNodesetGroup = transitionGroup.getNodesetGroup(nodes)

        annotationGroups = [boxGroup, transitionGroup]
        annotationTerms = ["box group", "transition group"]
        meshGroups = [boxMeshGroup, transitionMeshGroup]

        # annotation groups & nodeset groups
        lungGroup = AnnotationGroup(region, get_lung_term("lung"))
        leftLungGroup = AnnotationGroup(region, get_lung_term("left lung"))
        annotationGroups.append(lungGroup)
        annotationGroups.append(leftLungGroup)
        annotationTerms.append("left lung")
        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)
        meshGroups.append(leftLungMeshGroup)

        rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightLungGroup)
        annotationTerms.append("right lung")
        meshGroups.append(rightLungMeshGroup)

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

        mediastinumLeftGroupX = AnnotationGroup(region,["anterior mediastinum of left lung X", "None"])
        mediastinumLeftGroupMeshGroupX = mediastinumLeftGroupX.getMeshGroup(mesh)
        annotationGroups.append(mediastinumLeftGroupX)

        mediastinumLeftGroupY = AnnotationGroup(region,["anterior mediastinum of left lung Y", "None"])
        mediastinumLeftGroupMeshGroupY = mediastinumLeftGroupY.getMeshGroup(mesh)
        annotationGroups.append(mediastinumLeftGroupY)

        mediastinumLeftGroupZ = AnnotationGroup(region,["anterior mediastinum of left lung Z", "None"])
        mediastinumLeftGroupMeshGroupZ = mediastinumLeftGroupZ.getMeshGroup(mesh)
        annotationGroups.append(mediastinumLeftGroupZ)

        mediastinumRightGroupX = AnnotationGroup(region,["anterior mediastinum of right lung X", "None"])
        mediastinumRightGroupMeshGroupX = mediastinumRightGroupX.getMeshGroup(mesh)
        annotationGroups.append(mediastinumRightGroupX)

        mediastinumRightGroupY = AnnotationGroup(region,["anterior mediastinum of right lung Y", "None"])
        mediastinumRightGroupMeshGroupY = mediastinumRightGroupY.getMeshGroup(mesh)
        annotationGroups.append(mediastinumRightGroupY)

        mediastinumRightGroupZ = AnnotationGroup(region,["anterior mediastinum of right lung Z", "None"])
        mediastinumRightGroupMeshGroupZ = mediastinumRightGroupZ.getMeshGroup(mesh)
        annotationGroups.append(mediastinumRightGroupZ)

        leftLateralLungGroup = AnnotationGroup(region, ["lateral left lung", "None"])
        leftLateralLungGroupMeshGroup = leftLateralLungGroup.getMeshGroup(mesh)
        annotationGroups.append(leftLateralLungGroup)

        rightLateralLungGroup = AnnotationGroup(region, ["lateral right lung", "None"])
        rightLateralLungGroupMeshGroup = rightLateralLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightLateralLungGroup)

        # leftMedialLungGroup = AnnotationGroup(region, ["medial left lung", "None"])
        # leftMedialLungGroupMeshGroup = leftMedialLungGroup.getMeshGroup(mesh)
        # annotationGroups.append(leftMedialLungGroup)
        #
        # rightMedialLungGroup = AnnotationGroup(region, ["medial right lung", "None"])
        # rightMedialLungGroupMeshGroup = rightMedialLungGroup.getMeshGroup(mesh)
        # annotationGroups.append(rightMedialLungGroup)

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

        if hasAccessoryLobe: # currently not used
            # Annotation groups
            rightLungAccessoryLobeGroup = AnnotationGroup(region, get_lung_term("right lung accessory lobe"))
            rightLungAccessoryLobeMeshGroup = rightLungAccessoryLobeGroup.getMeshGroup(mesh)
            annotationGroups.append(rightLungAccessoryLobeGroup)
            rightLungAccessoryLobeNodesetGroup = rightLungAccessoryLobeGroup.getNodesetGroup(nodes)

            # Marker points
            accessoryDorsalApexGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_lung_term("dorsal apex of right lung accessory lobe"))
            accessoryVentralApexGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_lung_term("ventral apex of right lung accessory lobe"))
            accessoryVentralLeftGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_lung_term("left ventral base of right lung accessory lobe"))
            accessoryVentralRightGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_lung_term("right ventral base of right lung accessory lobe"))
            accessoryDorsalLeftGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_lung_term("left dorsal base of right lung accessory lobe"))
            accessoryDorsalRightGroup = findOrCreateAnnotationGroupForTerm(
                annotationGroups, region, get_lung_term("right dorsal base of right lung accessory lobe"))

        # Nodeset group
        leftLungNodesetGroup = leftLungGroup.getNodesetGroup(nodes)
        rightLungNodesetGroup = rightLungGroup.getNodesetGroup(nodes)
        lungNodesetGroup = lungGroup.getNodesetGroup(nodes)

        if numberOfLeftLung == 2:
            lowerLeftLungNodesetGroup = lowerLeftLungGroup.getNodesetGroup(nodes)
            upperLeftLungNodesetGroup = upperLeftLungGroup.getNodesetGroup(nodes)

            lowerRightLungNodesetGroup = lowerRightLungGroup.getNodesetGroup(nodes)
            middleRightLungNodesetGroup = middleRightLungGroup.getNodesetGroup(nodes)
            upperRightLungNodesetGroup = upperRightLungGroup.getNodesetGroup(nodes)

        # Arbitrary anatomical groups and nodesets
        upperLeftDorsalLungGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, ["upper lobe of left lung dorsal", "None"])
        upperLeftDorsalLungMeshGroup = upperLeftDorsalLungGroup.getMeshGroup(mesh)

        middleRightDorsalLungGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, ["middle lobe of right lung dorsal", "None"])
        middleRightDorsalLungMeshGroup = middleRightDorsalLungGroup.getMeshGroup(mesh)

        upperRightDorsalLungGroup = findOrCreateAnnotationGroupForTerm(
            annotationGroups, region, ["upper lobe of right lung dorsal", "None"])
        upperRightDorsalLungMeshGroup = upperRightDorsalLungGroup.getMeshGroup(mesh)


        # # Marker points/groups
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


        centre = [0.0, 0.0, 0.0]
        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]
        axes = [[disc_depth, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]

        elementsCountAcross = [elementsCountLateral, elementsCountOblique, elementsCountNormal]
        shellProportion = 1.0

        leftLung, rightLung = 0, 1
        lungs = [lung for show, lung in [(showLeftLung, leftLung), (showRightLung, rightLung)] if show]
        for i in lungs:
            if i == leftLung:
                if numberOfLeftLung == 2:
                    meshGroups = [boxMeshGroup, transitionMeshGroup,
                                  leftLungMeshGroup, lowerLeftLungMeshGroup, upperLeftLungMeshGroup,
                                  upperLeftDorsalLungMeshGroup,
                                  mediastinumLeftGroupMeshGroupX, mediastinumLeftGroupMeshGroupY, mediastinumLeftGroupMeshGroupZ,
                                  leftLateralLungGroupMeshGroup,] #leftMedialLungGroupMeshGroup]
                    annotationTerms = ["box group", "transition group", "left lung",
                                       "lower lobe of left lung", "upper lobe of left lung", "upper lobe of left lung dorsal",
                                       "anterior mediastinum of left lung X",
                                       "anterior mediastinum of left lung Y",
                                       "anterior mediastinum of left lung Z",
                                       "lateral left lung",] #"medial left lung"]
                else:
                    meshGroups = [boxMeshGroup, transitionMeshGroup,
                                  leftLungMeshGroup,
                                  mediastinumLeftGroupMeshGroupX, mediastinumLeftGroupMeshGroupY, mediastinumLeftGroupMeshGroupZ]
                    annotationTerms = ["box group", "transition group", "left lung",
                                       "anterior mediastinum of left lung X",
                                       "anterior mediastinum of left lung Y",
                                       "anterior mediastinum of left lung Z"]
            else:
                meshGroups = [boxMeshGroup, transitionMeshGroup,
                              rightLungMeshGroup, lowerRightLungMeshGroup, middleRightLungMeshGroup,
                              upperRightLungMeshGroup, upperRightDorsalLungMeshGroup,
                              mediastinumRightGroupMeshGroupX, mediastinumRightGroupMeshGroupY, mediastinumRightGroupMeshGroupZ,
                              rightLateralLungGroupMeshGroup,] #rightMedialLungGroupMeshGroup]
                annotationTerms = ["box group", "transition group", "right lung",
                                   "lower lobe of right lung", "middle lobe of right lung",
                                   "upper lobe of right lung", "upper lobe of right lung dorsal",
                                   "anterior mediastinum of right lung X",
                                   "anterior mediastinum of right lung Y",
                                   "anterior mediastinum of right lung Z",
                                   "lateral right lung",] #"medial right lung"]

            sphere = SphereMesh(fieldmodule, coordinates, centre, axes, elementsCountAcross,
                                elementsCountShell, elementsCountTransition, shellProportion,
                                sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                                boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False,
                                useSerendipity=useSerendipity, meshGroups=meshGroups, annotationTerms=annotationTerms)

        if useSizingFunction:
            scale_factor = options["Scale factor"]
            for i in lungs:
                isLeft = True if i == leftLung else False
                lungNodeset = leftLungNodesetGroup if isLeft else rightLungNodesetGroup
                scaleNodes(fieldmodule, coordinates, lungNodeset, disc_depth, scale_factor)
                if elementsCountLateral > 4:
                    resampleShellNodes(fieldmodule, coordinates, lungNodeset, elementsCountLateral,
                                       elementsCountNormal, elementsCountOblique, elementsCountShell)

        for i in lungs:
            isLeft = True if i == leftLung else False
            lungNodeset = leftLungNodesetGroup if isLeft else rightLungNodesetGroup
            halfBreadth = disc_breadth * 0.5
            height = disc_height
            spacing = -lungSpacing if i == leftLung else lungSpacing
            oblique_slope_radians = left_oblique_slope_radians if i == leftLung else right_oblique_slope_radians
            bendZ = (diaphragm_proportion - 0.5) * disc_height
            lungMedialCurvature = -leftLungMedialCurvature if isLeft else leftLungMedialCurvature
            apexMedialDisplacement = leftApexMedialDisplacement if isLeft else -leftApexMedialDisplacement
            rotateLungAngleX = -rotateLeftLungX
            rotateLungAngleZ = rotateLeftLungZ if isLeft else -rotateLeftLungZ

            rotate_scale_lung(fieldmodule, coordinates, lungNodeset, disc_breadth, disc_height, oblique_slope_radians)

            if edgeSharpFactor != 0.0:
                sharpeningRidge(edgeSharpFactor, fieldmodule, coordinates, lungNodeset, halfBreadth)

            if apexSharpFactor != 0.0:
                sharpeningRidge(apexSharpFactor, fieldmodule, coordinates, lungNodeset, halfBreadth, isApex=True)

            if impression_depth_proportion > 0.0:
                form_mediastinal_surface(fieldmodule, coordinates, lungNodeset, disc_breadth, disc_height, disc_depth,
                                         impression_breadth_proportion, impression_height_proportion,
                                         impression_depth_proportion, lateral_shear_rate)

            if (diaphragm_angle_radians != 0.0) and (diaphragm_proportion > 0.0):
                form_diaphragm_surface(fieldmodule, coordinates, lungNodeset, diaphragm_angle_radians, bendZ,
                                       isLeft=isLeft)

            dorsalVentralXi = getDorsalVentralXiField(fieldmodule, coordinates, halfBreadth)
            if lungMedialCurvature != 0:
                bendingAroundZAxis(lungMedialCurvature, fieldmodule, coordinates, lungNodeset,
                                   stationaryPointXY=[0.0, 0.0],
                                   bias=lungMedialCurvatureBias,
                                   dorsalVentralXi=dorsalVentralXi)

            if rotateLungAngleX != 0.0:
                rotateLungs(rotateLungAngleX, fieldmodule, coordinates, lungNodeset, axis=1)

            if rotateLungAngleZ != 0.0:
                rotateLungs(rotateLungAngleZ, fieldmodule, coordinates, lungNodeset, axis=3)

            translateLungLocation(fieldmodule, coordinates, lungNodeset, spacing, bendZ)

            if apexMedialDisplacement != 0.0:
                medialShearRadian = math.atan(apexMedialDisplacement / height)
                tiltLungs(medialShearRadian, 0, 0, 0, fieldmodule, coordinates, lungNodeset)

            if forwardLeftRightApex != 0.0:
                ventralShearRadian = math.atan(forwardLeftRightApex / height)
                tiltLungs(0, ventralShearRadian, 0, 0, fieldmodule, coordinates, lungNodeset, isApex=True)

            smoothD3ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountLateral,
                                         elementsCountNormal, elementsCountOblique, elementsCountShell)
            smoothD1ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountNormal,
                                         elementsCountOblique, elementsCountShell)
            smoothD2ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountNormal,
                                         elementsCountOblique, elementsCountShell)

            if isOpenfissure:
                setValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                                  Node.VALUE_LABEL_D_DS3]
                nodeParameters = get_nodeset_field_parameters(lungNodeset, coordinates, setValueLabels)[1]
                nodeIdentifier = nodeParameters[-1][0] + 1
                # if numberOfLeftLung > 1:
                nodeIdentifier, copyIdentifiersLLU = disconnectFieldMeshGroupBoundaryNodes(
                    [coordinates], lowerLeftLungMeshGroup, upperLeftLungMeshGroup,
                    nodeIdentifier)
                nodeIdentifier, copyIdentifiersRLM = disconnectFieldMeshGroupBoundaryNodes(
                    [coordinates], lowerRightLungMeshGroup, middleRightLungMeshGroup,
                    nodeIdentifier)
                nodeIdentifier, copyIdentifiersRLU = disconnectFieldMeshGroupBoundaryNodes(
                    [coordinates], lowerRightLungMeshGroup, upperRightLungMeshGroup, nodeIdentifier)
                nodeIdentifier, copyIdentifiersRMU = disconnectFieldMeshGroupBoundaryNodes(
                    [coordinates], middleRightLungMeshGroup, upperRightLungMeshGroup,
                    nodeIdentifier)

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
        openFissures = options['Open fissures']

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

        # 1D edge markers
        mediastanum_group = {}
        mediastanum_exterior = {}
        mediastanum_2dgroup = {}
        mediastanum_border = {}

        mediastanumTerms = [
            "anterior mediastinum of left lung X",
            "anterior mediastinum of left lung Y",
            "anterior mediastinum of left lung Z",
            "anterior mediastinum of right lung X",
            "anterior mediastinum of right lung Y",
            "anterior mediastinum of right lung Z"
        ]

        for term in mediastanumTerms:
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [term, "None"])
            group2d = group.getGroup()
            group2d_exterior = fm.createFieldAnd(group2d, is_exterior)
            mediastanum_group[term] = group
            mediastanum_2dgroup[term] = group2d
            mediastanum_exterior[term] = group2d_exterior

        axisTerms = [
            "left lung X", "left lung Y", "left lung Z",
            "right lung X", "right lung Y", "right lung Z"
        ]

        for axisTerm in axisTerms:
            mediastanumTerm = f"anterior mediastinum of {axisTerm}"
            arbBorderTerm = f"anterior border of {axisTerm}"
            borderSide = "left" if "left lung" in axisTerm else "right"
            borderTerm = f"anterior border of {borderSide} lung"

            borderGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(borderTerm))

            if "lung X" in axisTerm:
                group2d_border = mediastanum_exterior[mediastanumTerm]

            elif "lung Y" in axisTerm:
                tempField = fm.createFieldAnd(mediastanum_exterior[mediastanumTerm], is_xi1_1)
                borderX_term = f"anterior border of {borderSide} lung X"
                borderX = mediastanum_border[borderX_term]
                group2d_border = fm.createFieldAnd(tempField, fm.createFieldNot(borderX))
                borderGroup.getMeshGroup(mesh1d).addElementsConditional(group2d_border)

            elif "lung Z" in axisTerm:
                tempField = fm.createFieldAnd(mediastanum_exterior[mediastanumTerm], is_xi2_0)
                borderX_term = f"anterior border of {borderSide} lung X"
                borderX = mediastanum_border[borderX_term]
                group2d_border = fm.createFieldAnd(tempField, fm.createFieldNot(borderX))
                borderGroup.getMeshGroup(mesh1d).addElementsConditional(group2d_border)

            mediastanum_border[arbBorderTerm] = group2d_border

        for group in mediastanum_group.values():
            annotationGroups.remove(group)

        # Arbitrary terms - are removed from the annotation groups later
        arbLobe_group = {}
        arbLobe_exterior = {}
        arbLobe_2dgroup = {}
        arbTerms = ["upper lobe of left lung dorsal", "upper lobe of right lung dorsal"]
        for arbTerm in arbTerms:
            group = findOrCreateAnnotationGroupForTerm(annotationGroups, region, [arbTerm, "None"])
            group2d = group.getGroup()
            group2d_exterior = fm.createFieldAnd(group2d, is_exterior)
            arbLobe_group.update({arbTerm: group})
            arbLobe_2dgroup.update({arbTerm: group2d})
            arbLobe_exterior.update({arbTerm: group2d_exterior})

        side_group = {}
        side_exterior = {}
        arbSideTerms = ["lateral left lung", "lateral right lung", "medial left lung", "medial right lung"]
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

            # medial and lateral in the subgroup
            # for sideTerm in ['lateral surface of ', 'medial surface of ']:
            for sideTerm in ['lateral surface of ']:
                surfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(sideTerm + term))
                if ('lateral' in sideTerm) and ('left' in term):
                    surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(fm.createFieldAnd(group2d_exterior, side_exterior["lateral left lung"]))
                # elif ('medial' in sideTerm) and ('left' in term):
                #     surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(fm.createFieldAnd(group2d_exterior, side_exterior["medial left lung"]))
                elif ('lateral' in sideTerm) and ('right' in term):
                    surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(fm.createFieldAnd(group2d_exterior, side_exterior["lateral right lung"]))
                # elif ('medial' in sideTerm) and ('right' in term):
                #     surfaceGroup.getMeshGroup(mesh2d).addElementsConditional(fm.createFieldAnd(group2d_exterior, side_exterior["medial right lung"]))

        # Fissures

        if not openFissures:
            fissureTerms = ["oblique fissure of right lung", "horizontal fissure of right lung"]
            if numberOfLeftLung > 1:
                fissureTerms.append("oblique fissure of left lung")
            lobeFissureTerms = ["oblique fissure of lower lobe of left lung",
                                "oblique fissure of upper lobe of left lung",
                                "oblique fissure of lower lobe of right lung",
                                "oblique fissure of middle lobe of right lung",
                                "oblique fissure of upper lobe of right lung",
                                "horizontal fissure of middle lobe of right lung",
                                "horizontal fissure of upper lobe of right lung"]
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
                            fissureGroup = fm.createFieldAnd(fissureGroup_temp, arbLobe_2dgroup['upper lobe of right lung dorsal'])
                        elif "oblique fissure of middle lobe of right lung" in lobeFissureTerm:
                            fissureGroup = fm.createFieldAnd(fissureGroup_temp, fm.createFieldNot(arbLobe_2dgroup['upper lobe of right lung dorsal']))
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
                annotationGroups,
                get_lung_term("horizontal fissure of right lung")).getGroup()
            obliqueFissureOfRightLungGroup = getAnnotationGroupForTerm(
                annotationGroups,
                get_lung_term("oblique fissure of right lung")).getGroup()
            obliqueFissureOfMiddleLobeOfRightLungGroup = getAnnotationGroupForTerm(
                annotationGroups,
                get_lung_term("oblique fissure of middle lobe of right lung")).getGroup()
            obliqueFissureOfUpperLobeOfRightLungGroup = getAnnotationGroupForTerm(
                annotationGroups,
                get_lung_term("oblique fissure of upper lobe of right lung")).getGroup()
            lobeSurfaceGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("lower lobe of right lung surface"))
            lobeSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueFissureOfRightLungGroup)
            lobeSurfaceGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("middle lobe of right lung surface"))
            lobeSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
                fm.createFieldOr(obliqueFissureOfMiddleLobeOfRightLungGroup,
                                 horizontalFissureOfRightLungGroup))
            lobeSurfaceGroup = getAnnotationGroupForTerm(
                annotationGroups, get_lung_term("upper lobe of right lung surface"))
            lobeSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(
                fm.createFieldOr(obliqueFissureOfUpperLobeOfRightLungGroup,
                                 horizontalFissureOfRightLungGroup))

        if openFissures:
            if numberOfLeftLung > 1:
                obliqueLowerLeft = fm.createFieldOr(fm.createFieldAnd(lobe_exterior['lower lobe of left lung surface'], is_xi2_1),
                                                    fm.createFieldAnd(lobe_exterior['lower lobe of left lung surface'], is_xi3_1))

                tempUpperLeft = fm.createFieldXor(fm.createFieldAnd(lobe_exterior['upper lobe of left lung surface'], is_exterior),
                                                     arbLobe_exterior['upper lobe of left lung dorsal'])

                obliqueUpperLeft = fm.createFieldOr(fm.createFieldAnd(lobe_exterior['upper lobe of left lung surface'], is_xi2_0),
                                                    tempUpperLeft)

                fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('oblique fissure of ' + 'lower lobe of left lung'))
                fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueLowerLeft)

                fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('oblique fissure of ' + 'upper lobe of left lung'))
                fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueUpperLeft)

                obliqueLeft = fm.createFieldOr(obliqueUpperLeft, obliqueLowerLeft)
                fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('oblique fissure of ' + 'left lung'))
                fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueLeft)

                leftLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('left lung surface'))
                leftLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(fm.createFieldAnd(lobe_exterior['left lung surface'], fm.createFieldNot(obliqueLeft)))

                for term in ("lower lobe of left lung surface", "upper lobe of left lung surface"):
                    leftLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(
                        annotationGroups, region, get_lung_term(term))
                    leftLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(lobe_exterior[term])
            else:
                leftLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('left lung surface'))
                leftLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(lobe_exterior['left lung surface'])

            obliqueLowerRight = fm.createFieldOr(
                fm.createFieldAnd(lobe_exterior['lower lobe of right lung surface'], is_xi2_1),
                fm.createFieldAnd(lobe_exterior['lower lobe of right lung surface'], is_xi3_1))

            tempObliqueMiddleRight1 = fm.createFieldAnd(
                fm.createFieldAnd(lobe_exterior['middle lobe of right lung surface'], is_xi2_0),
                fm.createFieldNot(boxTransition_exterior["box group"]))

            tempObliqueMiddleRight2 = fm.createFieldAnd(
                fm.createFieldAnd(lobe_exterior['middle lobe of right lung surface'], is_xi3_0),
                boxTransition_exterior["box group"])

            obliqueMiddleRight = fm.createFieldOr(tempObliqueMiddleRight1, tempObliqueMiddleRight2)

            tempObliqueUpperRight1 = fm.createFieldAnd(
                fm.createFieldAnd(lobe_exterior['upper lobe of right lung surface'], is_xi3_0),
                boxTransition_exterior["box group"])

            obliqueUpperRight = fm.createFieldOr(
                tempObliqueUpperRight1,
                fm.createFieldAnd(arbLobe_exterior['upper lobe of right lung dorsal'], is_xi2_0))

            obliqueRight = fm.createFieldOr(fm.createFieldOr(obliqueLowerRight, obliqueUpperRight), obliqueMiddleRight)

            tempHorizontalMiddleRight1 = fm.createFieldAnd(
                lobe_exterior['middle lobe of right lung surface'],
                fm.createFieldNot(obliqueMiddleRight))

            tempHorizontalMiddleRight2 = fm.createFieldAnd(lobe_exterior['middle lobe of right lung surface'], is_xi3_0)

            horizontalMiddleRight = fm.createFieldAnd(tempHorizontalMiddleRight1, fm.createFieldNot(tempHorizontalMiddleRight2))

            tempHorizontalUpperRight1 = fm.createFieldOr(
                fm.createFieldAnd(lobe_exterior['upper lobe of right lung surface'], is_xi1_0),
                fm.createFieldAnd(lobe_exterior['upper lobe of right lung surface'], is_xi1_1))

            horizontalUpperRight = fm.createFieldOr(
                tempHorizontalUpperRight1,
                fm.createFieldAnd(lobe_exterior['upper lobe of right lung surface'], is_xi2_1))

            horizontalRight = fm.createFieldOr(horizontalMiddleRight, horizontalUpperRight)

            fissureRight = fm.createFieldOr(obliqueRight, horizontalRight)

            rightLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('right lung surface'))
            rightLungSurface = fm.createFieldAnd(lobe_exterior['right lung surface'], fm.createFieldNot(fissureRight))
            rightLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(rightLungSurface)

            for term in ("lower lobe of right lung surface", "middle lobe of right lung surface",
                         "upper lobe of right lung surface"):
                rightLungSurfaceGroup = findOrCreateAnnotationGroupForTerm(
                    annotationGroups, region, get_lung_term(term))
                rightLungSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(lobe_exterior[term])

            fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('oblique fissure of ' + 'lower lobe of right lung'))
            fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueLowerRight)

            fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('oblique fissure of ' + 'middle lobe of right lung'))
            fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueMiddleRight)

            fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('oblique fissure of ' + 'upper lobe of right lung'))
            fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueUpperRight)

            fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('horizontal fissure of ' + 'middle lobe of right lung'))
            fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(horizontalMiddleRight)

            fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('horizontal fissure of ' + 'upper lobe of right lung'))
            fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(horizontalUpperRight)

            fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('oblique fissure of ' + 'right lung'))
            fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(obliqueRight)

            fissureSurfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term('horizontal fissure of ' + 'right lung'))
            fissureSurfaceGroup.getMeshGroup(mesh2d).addElementsConditional(horizontalRight)

        for key, value in arbLobe_group.items():
            annotationGroups.remove(value)

        for key, group in side_group.items():
            annotationGroups.remove(group)

        for key, group in boxTransition_group.items():
            annotationGroups.remove(group)

def rotate_scale_lung(fieldmodule, coordinates, lungNodesetGroup, disc_breadth, disc_height, oblique_slope_radians):
    """
    Rotate and scale lung disc to achieve oblique slope, height and breadth.
    Transformation tries to conserve original right angle between major and minor axes,
    and keep element sizes around to similar sizes.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param disc_breadth: Breadth of lung ellipse.
    :param disc_height: Height of lung ellipse.
    :param oblique_slope_radians: Slope of oblique fissure, where a positive angle tilts down ventrally.
    :return: None
    """
    x = fieldmodule.createFieldComponent(coordinates, 1)
    y = fieldmodule.createFieldComponent(coordinates, 2)
    z = fieldmodule.createFieldComponent(coordinates, 3)
    yz = fieldmodule.createFieldComponent(coordinates, [2, 3])
    r = fieldmodule.createFieldMagnitude(yz)
    theta = fieldmodule.createFieldAtan2(z, y)
    zero = fieldmodule.createFieldConstant(0.0)
    r_zero = fieldmodule.createFieldEqualTo(r, zero)
    safe_theta = fieldmodule.createFieldIf(r_zero, zero, theta)
    oblique_angle = fieldmodule.createFieldConstant(oblique_slope_radians)
    phi = safe_theta - oblique_angle
    cos_phi = fieldmodule.createFieldCos(phi)
    sin_phi = fieldmodule.createFieldSin(phi)

    scaling = 2.0 * (math.atan(disc_breadth / disc_height) - 0.25 * math.pi)
    cc = fieldmodule.createFieldConstant(scaling)
    delta_phi = cc * cos_phi * sin_phi
    new_theta = phi + delta_phi
    cos_new_theta = fieldmodule.createFieldCos(new_theta)
    sin_new_theta = fieldmodule.createFieldSin(new_theta)
    new_y = r * cos_new_theta * fieldmodule.createFieldConstant(disc_breadth)
    new_z = r * sin_new_theta * fieldmodule.createFieldConstant(disc_height)
    new_coordinates = fieldmodule.createFieldConcatenate([x, new_y, new_z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def form_mediastinal_surface(fieldmodule, coordinates, lungNodesetGroup, disc_breadth, disc_height, disc_depth,
                             impression_breadth_proportion, impression_height_proportion,
                             impression_depth_proportion, lateral_shear_rate):
    """
    Form mediastinal surface by draping over ellipsoid impression.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param disc_breadth: Breadth of lung disc ellipse.
    :param disc_height: Height of lung disc ellipse.
    :param disc_depth: Depth/thickness of lung disc half ellipsoid.
    :param impression_breadth_proportion: Impression breadth as proportion of disc breadth, circle chord.
    :param impression_height_proportion: Impression height as proportion of disc height, circle chord.
    :param impression_depth_proportion: Impression depth as proportion of disc thickness, giving radii.
    :param lateral_shear_rate: Rate of shear as proportion of depth into disc.
    :return: None.
    """
    # get scaling from ellipse to curve of ellipsoid impression
    b = 0.5 * impression_breadth_proportion * disc_breadth
    h = 0.5 * impression_height_proportion * disc_height
    d = impression_depth_proportion * disc_depth
    breadth_radius = 0.5 * (b * b + d * d) / d
    breadth_angle = 2.0 * math.asin(b / breadth_radius)
    breadth_cx = breadth_radius - d
    breadth_cy = 0.5 * disc_breadth - b
    height_radius = 0.5 * (h * h + d * d) / d
    height_angle = 2.0 * math.asin(h / height_radius)
    height_cx = height_radius - d
    height_cy = breadth_cy

    x = fieldmodule.createFieldComponent(coordinates, 1)
    y = fieldmodule.createFieldComponent(coordinates, 2)
    z = fieldmodule.createFieldComponent(coordinates, 3)
    b_angle_factor = fieldmodule.createFieldConstant(breadth_angle / disc_breadth)
    b_angle = fieldmodule.createFieldMultiply(y, b_angle_factor)
    h_angle_factor = fieldmodule.createFieldConstant(height_angle / disc_height)
    h_angle = fieldmodule.createFieldMultiply(z, h_angle_factor)
    bcx = fieldmodule.createFieldConstant(breadth_cx)
    dd = fieldmodule.createFieldConstant(d)
    bcy = fieldmodule.createFieldConstant(breadth_cy)
    rb = bcx - x + dd
    b_shear_rate = fieldmodule.createFieldConstant(lateral_shear_rate * disc_depth / breadth_radius)
    mod_b_angle = b_angle + b_shear_rate * x
    cosb = fieldmodule.createFieldCos(mod_b_angle)
    sinb = fieldmodule.createFieldSin(mod_b_angle)
    new_x = bcx - rb * cosb
    new_y = bcy + rb * sinb

    new_coordinates = fieldmodule.createFieldConcatenate([new_x, new_y, z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def form_diaphragm_surface(fieldmodule, coordinates, lungNodesetGroup, diaphragm_angle_radians, diaphragm_bend_z,
                           isLeft=True):
    """
    Form diaphragm surface by bending model up, phasing out to nothing at the apex.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param diaphragm_angle_radians: Angle of diaphragm from down mediastinal surface.
    :param diaphragm_bend_z: Z coordinate of centre of diaphragm bend.
    :param isLeft: True if left lung, False if right lung.
    :return: None.
    """
    x = fieldmodule.createFieldComponent(coordinates, 1)
    y = fieldmodule.createFieldComponent(coordinates, 2)
    z = fieldmodule.createFieldComponent(coordinates, 3)
    bend_z = fieldmodule.createFieldConstant(diaphragm_bend_z)
    minus_1 = fieldmodule.createFieldConstant(-1.0)
    nx = x * minus_1
    mz = z - bend_z
    nx_mz = fieldmodule.createFieldConcatenate([nx, mz])
    r = fieldmodule.createFieldMagnitude(nx_mz)

    theta = fieldmodule.createFieldAtan2(mz, nx)
    zero = fieldmodule.createFieldConstant(0.0)
    r_zero = fieldmodule.createFieldEqualTo(r, zero)
    safe_theta = fieldmodule.createFieldIf(r_zero, zero, theta)
    pi__2 = fieldmodule.createFieldConstant(0.5 * math.pi)
    delta_theta = (pi__2 - safe_theta) * fieldmodule.createFieldConstant(diaphragm_angle_radians / math.pi)
    new_theta = safe_theta + delta_theta
    cos_new_theta = fieldmodule.createFieldCos(new_theta) * minus_1 if isLeft else fieldmodule.createFieldCos(new_theta)
    sin_new_theta = fieldmodule.createFieldSin(new_theta)
    new_x = r * cos_new_theta
    new_z = bend_z + r * sin_new_theta

    new_coordinates = fieldmodule.createFieldConcatenate([new_x, y, new_z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def translateLungLocation(fm, coordinates, lungNodesetGroup, spaceFromCentre, bendZ):
    """
    Translates the lung mesh along the x- and z-axes to adjust its position within the anatomical space.
    :param fm: Field module used for creating and managing fields.
    :param coordinates: The coordinate field representing the initial position of the lung, with a circular profile in the y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing the nodes of the lung to be transformed.
    :param spaceFromCentre: The distance to translate the lung along the x-axis, allowing lateral positioning relative
    to the centerline (i.e. x = 0) of the body.
    :param bendZ: The amount to translate the lung along the z-axis, typically at a position where the diaphragm and the
    medial surface meets.
    :return: None
    """
    x = fm.createFieldComponent(coordinates, 1)
    y = fm.createFieldComponent(coordinates, 2)
    z = fm.createFieldComponent(coordinates, 3)
    spaceFromCentre = fm.createFieldConstant(spaceFromCentre)

    new_x = x + spaceFromCentre
    new_z = z - fm.createFieldConstant(bendZ)
    new_coordinates = fm.createFieldConcatenate([new_x, y, new_z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def rotateLungs(rotateAngle, fm, coordinates, lungNodesetGroup, axis):
    """
    Rotate the lung mesh at the center of the elements about z-axis
    :param rotateAngle: Angle of rotation in degrees.
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param axis: Axis of rotation.
    :return: None
    """
    rotateAngle = -rotateAngle / 180 * math.pi  # negative value due to right handed rule

    if axis == 1:
        rotateMatrix = fm.createFieldConstant([1.0, 0.0, 0.0,
                                               0.0, math.cos(rotateAngle), math.sin(rotateAngle),
                                               0.0, -math.sin(rotateAngle), math.cos(rotateAngle)])

    elif axis == 3:
        rotateMatrix = fm.createFieldConstant([math.cos(rotateAngle), math.sin(rotateAngle), 0.0,
                                               -math.sin(rotateAngle), math.cos(rotateAngle), 0.0,
                                               0.0, 0.0, 1.0])

    translate_coordinates = fm.createFieldMatrixMultiply(3, rotateMatrix, coordinates)

    fieldassignment = coordinates.createFieldassignment(translate_coordinates)
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


def bendingAroundZAxis(curvature, fm, coordinates, lungNodesetGroup, stationaryPointXY, bias=0.0, dorsalVentralXi=None):
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


def sharpeningRidge(sharpeningFactor, fm, coordinates, lungNodesetGroup, halfBreadth, isApex=False):
    """
    Linearly transforms the coordinates of the edge on the anterior side of the lung giving the effect of sharpening
    if isApex is False. If isApex is True, then the function transforms the coordinates of the edge on the posterior side
    giving the effect of blunting.
    :param sharpeningFactor: A value between 0 and 1, where 1 represents the maximum sharpness.
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param halfBreadth: Half breadth of lung.
    :param isApex: False if transforming the anterior edge, True if transforming the posterior edge of the lung.
    """
    # Transformation matrix = [ -k1y + 1, | [x,
    #                                 1, |  y,
    #                                 1] |  z]
    yLength = 0.1 * halfBreadth if isApex else 0.75 * halfBreadth
    offset = fm.createFieldConstant([0.0, yLength, 0.0])
    origin = fm.createFieldAdd(coordinates, offset)
    k1 = -sharpeningFactor / (halfBreadth * 1.75) if isApex else -sharpeningFactor / (halfBreadth * 1.75)
    scale = fm.createFieldConstant([0.0, -k1, 0.0]) if isApex else fm.createFieldConstant([0.0, k1, 0.0])
    # scale = fm.createFieldConstant([0.0, -k1, 0.0]) if isApex else fm.createFieldConstant([0.0, k1, 0.0])

    scaleFunction = fm.createFieldMultiply(origin, scale)
    constant = fm.createFieldConstant([0.0, 1.0, 1.0])
    constantFunction = fm.createFieldAdd(scaleFunction, constant)
    componentIndices = [3, 2, 3] if isApex else [2, 3, 3]
    # componentIndices = [3, 3, 2] if isApex else [2, 3, 3]

    transformation_matrix = fm.createFieldComponent(constantFunction, componentIndices)
    taper_coordinates = fm.createFieldMultiply(origin, transformation_matrix)
    translate_coordinates = fm.createFieldSubtract(taper_coordinates, offset)

    # if isApex:
    #     y = fm.createFieldComponent(coordinates, 2)
    #     isRightSide = fm.createFieldGreaterThan(y, fm.createFieldConstant(-0.05))
    #     translate_coordinates = fm.createFieldIf(isRightSide, coordinates, translate_coordinates)

    fieldassignment = coordinates.createFieldassignment(translate_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def tiltLungs(tiltApex_xAxis, tiltApex_yAxis, tiltDiap_yAxis, tiltDiap_xAxis, fm, coordinates, lungNodesetGroup,
              isApex=False):
    """
    Applies a shear transformation to the lung mesh to achieve tilting effects at both the apex and diaphragm.
    :param tiltApex_xAxis: Shear factor along the x-axis for the lung apex, tilting the top of the lung forward or backward.
    :param tiltApex_yAxis: Shear factor along the y-axis for the lung apex, tilting the top of the lung side-to-side.
    :param tiltDiap_yAxis: Shear factor along the y-axis for the diaphragm, controlling side-to-side tilt at the base.
    :param tiltDiap_xAxis: Shear factor along the x-axis for the diaphragm, controlling forward or backward tilt at the base.
    :param fm: Field module used to create and manage fields and transformations.
    :param coordinates: The coordinate field of the lung, initially with a circular profile in the y-z plane.
    :param lungNodesetGroup: The Zinc NodesetGroup containing the nodes to be transformed.
    :param isApex: If True, only transform coordinates above the lung base (apex region). If False, transform
    coordinates at or below the lung base (diaphragm region).
    """
    # FieldConstant - Matrix = [   x1,    x4, sh_zx,
    #                              x2,    x5, sh_zy,
    #                           sh_xz, sh_yz,    x9]
    sh_xz = tiltApex_xAxis
    sh_yz = tiltApex_yAxis
    sh_zy = tiltDiap_yAxis
    sh_zx = tiltDiap_xAxis

    shearMatrix = fm.createFieldConstant([1.0, 0.0, sh_xz, 0.0, 1.0, sh_yz, sh_zx, sh_zy, 1.0])
    newCoordinates = fm.createFieldMatrixMultiply(3, shearMatrix, coordinates)

    if tiltApex_yAxis != 0.0:
        z = fm.createFieldComponent(coordinates, 3)
        bendZ = fm.createFieldConstant(0.0)

        z_tilt = fm.createFieldGreaterThan(z, bendZ) if isApex else fm.createFieldLessThan(z, bendZ)
        newCoordinates = fm.createFieldIf(z_tilt, newCoordinates, coordinates)

    fieldassignment = coordinates.createFieldassignment(newCoordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def scaleNodes(fieldmodule, coordinates, lungNodesetGroup, disc_depth, scale_factor):
    """
    Redistribute lung nodes so that the outermost nodes remain in place while inner nodes are scaled outward.
    The scaling factor is 1 at the outer boundary and increases towards the center, causing nodes to shift outward.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in the y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param disc_depth: Depth of lung ellipse.
    :param scale_factor: Maximum scaling effect at the center (should be > 1 for expansion).
    """
    x = fieldmodule.createFieldComponent(coordinates, 1)
    y = fieldmodule.createFieldComponent(coordinates, 2)
    y_abs = fieldmodule.createFieldAbs(y)
    z = fieldmodule.createFieldComponent(coordinates, 3)

    zero = fieldmodule.createFieldConstant(0.0)
    one = fieldmodule.createFieldConstant(1.0)
    four = fieldmodule.createFieldConstant(4.0)

    tol = 1.0E-3

    aa = fieldmodule.createFieldConstant((disc_depth + disc_depth * tol) ** 2)
    b = fieldmodule.createFieldConstant(0.5)
    cc = fieldmodule.createFieldConstant((0.5 + 0.5 * tol) ** 2)

    exy = fieldmodule.createFieldSqrt((one - (x * x) / aa - (z * z) / cc))
    b_new = b * exy
    xi = y_abs / b_new
    xi2 = xi * xi
    xi3 = xi2 * xi

    beta = fieldmodule.createFieldConstant(scale_factor * 0.1)

    y_lt_zero = fieldmodule.createFieldLessThan(y, zero)

    # anterior side
    dy_a = four * (xi - xi2) * beta * exy
    new_y_a = y + fieldmodule.createFieldIf(exy, dy_a, zero)

    # posterior side
    dy_p = four * (xi3 - xi2) * beta * exy
    new_y_p = y + dy_p

    new_y = fieldmodule.createFieldIf(y_lt_zero, new_y_p, new_y_a)
    new_coordinates = fieldmodule.createFieldConcatenate([x, new_y, z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def resampleShellNodes(fieldmodule, coordinates, lungNodeset, elementsCountLateral, elementsCountNormal,
                       elementsCountOblique, elementsCountShell):
    """
    Resample shell nodes around the ellipsoid in the yz plane so that the shell nodes are aligned linearly.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in the y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    :param elementsCountShell: Number of elements across the shell layer.
    """
    elementsCountLateral -= 2 * elementsCountShell
    elementsCountNormal -= 2 * elementsCountShell
    elementsCountOblique -= 2 * elementsCountShell

    setValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                      Node.VALUE_LABEL_D_DS3]
    nodes = lungNodeset
    nodeParameters = get_nodeset_field_parameters(nodes, coordinates, setValueLabels)[1]
    sNodeIdentifier = nodeParameters[0][0]
    nodeParametersDict = dict(nodeParameters)
    elementsCountAround = getElementsCountAround(elementsCountNormal, elementsCountOblique)
    boundaryNids = getBoundaryNodeIdentifiers(elementsCountLateral, elementsCountNormal, elementsCountOblique,
                                              elementsCountAround, elementsCountShell, sNodeIdentifier=sNodeIdentifier)

    fieldcache = fieldmodule.createFieldcache()

    # resample across oblique direction
    normalMidIndex = (elementsCountNormal - 1) // 2
    for e3 in [0, -1]:
        for e2 in range(1, elementsCountNormal - 2):
            if (e3 == 0 and e2 < normalMidIndex) or (e3 == -1 and e2 >= normalMidIndex):
                sNodeIdentifier = [boundaryNids[0][e2][e3]]
                mNodeIdentifier = [boundaryNids[c][e2][e3] for c in range(1, elementsCountLateral - 2)]
                eNodeIdentifier = [boundaryNids[-1][e2][e3]]
            else:
                sNodeIdentifier = [boundaryNids[-1][e2][e3]]
                mNodeIdentifier = [boundaryNids[c][e2][e3] for c in range(-2, -(elementsCountLateral - 1), -1)]
                eNodeIdentifier = [boundaryNids[0][e2][e3]]
            snids = sNodeIdentifier + mNodeIdentifier + eNodeIdentifier

            bx, bd1 = zip(*[nodeParametersDict.get(id)[:2] for id in snids])
            bx, bd1 = [b[0] for b in bx], [d[0] for d in bd1]

            tx, td1 = sampleCubicHermiteCurves(bx, bd1, elementsCountLateral - 2, arcLengthDerivatives=True)[:2]
            if elementsCountLateral > 4:
                td1[1:-1] = smoothCubicHermiteDerivativesLine(tx[1:-1], td1[1:-1])

            for n in range(len(snids)):
                node = nodes.findNodeByIdentifier(snids[n])
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, tx[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, td1[n])

    lateralMidIndex = (elementsCountLateral - 1) // 2
    for e1 in range(1, elementsCountLateral - 2):
        for e2 in [0, -1]:
            for e3 in [0, -1]:
                sNodeIdentifier = [boundaryNids[e1][normalMidIndex][e3]]
                eNodeIdentifier = [boundaryNids[e1][e2][e3]]
                midRange = range(normalMidIndex - 1, 0, -1) if e2 == 0 else range(normalMidIndex + 1,
                                                                                  elementsCountNormal - 2)
                mNodeIdentifier = [boundaryNids[e1][c][e3] for c in midRange]
                snids = sNodeIdentifier + mNodeIdentifier + eNodeIdentifier

                bx, bd, n2List = [], [], []
                for n, id in enumerate(snids):
                    n2 = 1 if (e1 > lateralMidIndex and e2 == e3 and n == len(snids) - 1) or \
                              (e1 < lateralMidIndex and e2 != e3 and n == len(snids) - 1) else 2
                    n2List.append(n2)
                    bx.append(nodeParametersDict.get(id)[0][0])
                    bd.append(nodeParametersDict.get(id)[n2][0])
                bd[0] = mult(bd[0], -1) if e2 == 0 else bd[0]

                tx, td = sampleCubicHermiteCurves(bx, bd, len(bx) - 1, elementLengthStartEndRatio=0.7,
                                                  arcLengthDerivatives=True)[:2]

                for n, id in enumerate(snids):
                    n2 = n2List[n]
                    node = nodes.findNodeByIdentifier(id)
                    fieldcache.setNode(node)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, tx[n])
                    label = Node.VALUE_LABEL_D_DS1 if n2 == 1 else Node.VALUE_LABEL_D_DS2
                    coordinates.setNodeParameters(fieldcache, -1, label, 1, td[n])

    # resample across normal direction
    obliqueMidIndex = (elementsCountOblique - 1) // 2
    for e3 in range(1, elementsCountOblique - 2):
        if e3 == obliqueMidIndex:
            continue
        else:
            for e2 in [0, -1]:
                sNodeIdentifier = [boundaryNids[0][e2][e3]]
                eNodeIdentifier = [boundaryNids[-1][e2][e3]]
                mNodeIdentifier = [boundaryNids[c][e2][e3] for c in range(1, elementsCountLateral - 2)]
                snids = sNodeIdentifier + mNodeIdentifier + eNodeIdentifier

                bx, bd, n2List, dirList = [], [], [], []
                for n, id in enumerate(snids):
                    halfIndex = len(snids) // 2
                    isLeftHalf = n <= halfIndex
                    isRightHalf = n < halfIndex

                    if (e2 == 0 and e3 < obliqueMidIndex) or (e2 != 0 and e3 > obliqueMidIndex):
                        n2 = 1 if isLeftHalf else 2
                        dir = 1 if isLeftHalf else -1
                    else:
                        n2 = 2 if isRightHalf else 1
                        dir = 1 if isRightHalf else -1

                    n2List.append(n2)
                    dirList.append(dir)
                    bx.append(nodeParametersDict[id][0][0])
                    bd.append(mult(nodeParametersDict[id][n2][0], dir))

                tx, td = sampleCubicHermiteCurves(bx, bd, elementsCountLateral - 2, arcLengthDerivatives=True)[:2]

                for n, id in enumerate(snids):
                    n2, dir = n2List[n], dirList[n]
                    node = nodes.findNodeByIdentifier(id)
                    fieldcache.setNode(node)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, tx[n])
                    label = Node.VALUE_LABEL_D_DS1 if n2 == 1 else Node.VALUE_LABEL_D_DS2
                    coordinates.setNodeParameters(fieldcache, -1, label, 1, mult(td[n], dir))

    for e1 in range(1, elementsCountLateral - 2):
        for e2, e3 in [(0, 0), (0, -1), (-1, 0), (-1, -1)]:
            sNodeIdentifier = [boundaryNids[e1][e2][e3]]
            eNodeIdentifier = [boundaryNids[e1][e2][obliqueMidIndex]]
            midRange = range(e3 + 1, obliqueMidIndex) if e3 == 0 else range(elementsCountOblique - 3, obliqueMidIndex,
                                                                            -1)
            mNodeIdentifier = [boundaryNids[e1][e2][c] for c in midRange]
            snids = sNodeIdentifier + mNodeIdentifier + eNodeIdentifier

            bx, bd, n2List, dirList = [], [], [], []
            for n, id in enumerate(snids):
                isLast = (n == len(snids) - 1)
                if e1 == lateralMidIndex:
                    n2, dir = 2, (-1 if isLast and (e2 == e3) else 1)
                else:
                    isRightSide = e1 > lateralMidIndex
                    sameEdge = e2 == e3
                    if isLast:
                        n2 = 1
                        dir = -1 if (isRightSide and not sameEdge) or (not isRightSide and sameEdge) else 1
                    else:
                        n2 = 1 if (isRightSide and sameEdge) or (not isRightSide and not sameEdge) else 2
                        dir = 1

                n2List.append(n2)
                dirList.append(dir)
                bx.append(nodeParametersDict[id][0][0])
                bd.append(mult(nodeParametersDict[id][n2][0], dir))

            tx, td = sampleCubicHermiteCurves(bx, bd, len(bx) - 1, elementLengthStartEndRatio=1.3,
                                              arcLengthDerivatives=True)[:2]

            for n, id in enumerate(snids):
                n2, dir = n2List[n], dirList[n]
                node = nodes.findNodeByIdentifier(id)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, tx[n])
                label = Node.VALUE_LABEL_D_DS1 if n2 == 1 else Node.VALUE_LABEL_D_DS2
                coordinates.setNodeParameters(fieldcache, -1, label, 1, mult(td[n], dir))


def smoothD3ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountLateral, elementsCountNormal,
                                 elementsCountOblique, elementsCountShell):
    """
    Smooth derivatives of nodes connecting the box mesh to the shell mesh.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field.
    :param lungNodeset: Zinc NodesetGroup containing nodes to transform.
    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    :param elementsCountShell: Number of elements across the shell layer.
    """
    elementsCountLateral -= 2 * elementsCountShell
    elementsCountNormal -= 2 * elementsCountShell
    elementsCountOblique -= 2 * elementsCountShell

    setValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                      Node.VALUE_LABEL_D_DS3]
    nodes = lungNodeset
    nodeParameters = get_nodeset_field_parameters(nodes, coordinates, setValueLabels)[1]
    sNodeIdentifier = nodeParameters[0][0]
    nodeParametersDict = dict(nodeParameters)
    elementsCountAround = getElementsCountAround(elementsCountNormal, elementsCountOblique)

    nNodes = len(nodeParameters)
    nodeIds = getBoxRowNodeIdentifiers(elementsCountNormal, elementsCountOblique, nNodes, sNodeIdentifier)
    boundaryNodeIds = getBoundaryNodeIdentifiers(elementsCountLateral, elementsCountNormal, elementsCountOblique,
                                                 elementsCountAround, elementsCountShell, sNodeIdentifier)

    for e1 in [0, -1]:
        e1p = e1 + 1 if e1 == 0 else e1 - 1
        for e2 in range(elementsCountNormal - 1):
            for e3 in range(elementsCountOblique - 1):
                sid = nodeIds[e1][e2][e3]
                bid = nodeIds[e1p][e2][e3]
                bx = nodeParametersDict[bid][0][0]
                sx = nodeParametersDict[sid][0][0]

                shellFactor = 1.0  # 1.0 / elementsCountRim
                sd3 = mult(sub(sx, bx), shellFactor)

                fieldcache = fieldmodule.createFieldcache()

                node = nodes.findNodeByIdentifier(sid)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, sd3)

    for e1 in range(1, elementsCountLateral - 2):
        for e2 in range(elementsCountNormal - 1):
            for e3 in range(len(boundaryNodeIds[e1][e2])):
                sid = boundaryNodeIds[e1][e2][e3]

                if e2 in [0, elementsCountNormal - 2]:
                    bid = sid + (elementsCountOblique - 1) if e2 == 0 else sid - (elementsCountOblique - 1)
                else:
                    bid = sid + 1 if e3 == 0 else sid - 1

                bx = nodeParametersDict[bid][0][0]
                sx = nodeParametersDict[sid][0][0]

                shellFactor = 1.0  # 1.0 / elementsCountRim
                sd3 = mult(sub(sx, bx), shellFactor)

                fieldcache = fieldmodule.createFieldcache()
                node = nodes.findNodeByIdentifier(sid)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, sd3)


def smoothD1ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountNormal, elementsCountOblique,
                                 elementsCountShell):
    """
    Smooth derivatives of nodes connecting the box mesh to the shell mesh.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field.
    :param lungNodeset: Zinc NodesetGroup containing nodes to transform.
    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    :param elementsCountShell: Number of elements across the shell layer.
    """
    elementsCountNormal -= 2 * elementsCountShell
    elementsCountOblique -= 2 * elementsCountShell

    setValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                      Node.VALUE_LABEL_D_DS3]
    nodes = lungNodeset
    nodeParameters = get_nodeset_field_parameters(nodes, coordinates, setValueLabels)[1]
    sNodeIdentifier = nodeParameters[0][0]
    nodeParametersDict = dict(nodeParameters)

    nNodes = len(nodeParameters)
    nodeIds = getBoxRowNodeIdentifiers(elementsCountNormal, elementsCountOblique, nNodes, sNodeIdentifier)

    midNormalIndex = (elementsCountNormal - 1) // 2
    for e1 in [0, -1]:
        for e2 in range(1, elementsCountNormal - 2):
            snids = nodeIds[e1][e2]
            sx, sd1 = [], []
            if e1 == -1:
                dir = 1 if e2 < midNormalIndex else -1
            else:
                dir = -1 if e2 < midNormalIndex else 1
            for id in snids:
                sx.append(nodeParametersDict[id][0][0])
                sd1.append(mult(nodeParametersDict[id][1][0], dir))

            sd1 = smoothCubicHermiteDerivativesLine(sx, sd1,
                                                    fixStartDirection=True, fixEndDirection=True,
                                                    magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)

            fieldcache = fieldmodule.createFieldcache()
            for n, id in enumerate(snids):
                node = nodes.findNodeByIdentifier(id)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, mult(sd1[n], dir))


def smoothD2ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountNormal, elementsCountOblique,
                                 elementsCountShell):
    """
    Smooth derivatives of nodes connecting the box mesh to the shell mesh.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field.
    :param lungNodeset: Zinc NodesetGroup containing nodes to transform.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    :param elementsCountShell: Number of elements across the shell layer.
    """
    elementsCountNormal -= 2 * elementsCountShell
    elementsCountOblique -= 2 * elementsCountShell

    setValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                      Node.VALUE_LABEL_D_DS3]
    nodes = lungNodeset
    nodeParameters = get_nodeset_field_parameters(nodes, coordinates, setValueLabels)[1]
    sNodeIdentifier = nodeParameters[0][0]
    nodeParametersDict = dict(nodeParameters)
    elementsCountAround = getElementsCountAround(elementsCountNormal, elementsCountOblique)

    nNodes = len(nodeParameters)
    nodeIds = getBoxRowNodeIdentifiers(elementsCountNormal, elementsCountOblique, nNodes, sNodeIdentifier)

    midNormalIndex = (elementsCountNormal - 1) // 2
    midObliqueIndex = (elementsCountOblique - 1) // 2
    for e1 in [0, -1]:
        for e3 in range(1, elementsCountOblique - 2):
            snids = []
            sx, sd2 = [], []
            n2List, dirList = [], []
            for e2 in range(elementsCountNormal - 1):
                id = nodeIds[e1][e2][e3]
                if e1 == 0:
                    n2 = 1 if (e2 == 0 and e3 < midObliqueIndex) or (
                                e2 == elementsCountNormal - 2 and e3 > midObliqueIndex) else 2
                else:
                    n2 = 1 if (e2 == 0 and e3 > midObliqueIndex) or (
                                e2 == elementsCountNormal - 2 and e3 < midObliqueIndex) else 2
                dir = -1 if e2 < midNormalIndex else 1
                snids.append(id)
                n2List.append(n2)
                dirList.append(dir)
                sx.append(nodeParametersDict[id][0][0])
                sd2.append(mult(nodeParametersDict[id][n2][0], dir))

            sd2 = smoothCubicHermiteDerivativesLine(sx, sd2, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)

            fieldcache = fieldmodule.createFieldcache()
            for n, id in enumerate(snids):
                node = nodes.findNodeByIdentifier(id)
                n2 = n2List[n]
                dir = dirList[n]
                label = Node.VALUE_LABEL_D_DS1 if n2 == 1 else Node.VALUE_LABEL_D_DS2
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
                coordinates.setNodeParameters(fieldcache, -1, label, 1, mult(sd2[n], dir))


def getElementsCountAround(elementsCountNormal, elementsCountOblique):
    return (elementsCountNormal - 1) * 2 + (elementsCountOblique - 3) * 2


def getBoxRowNodeIdentifiers(elementsCountNormal, elementsCountOblique, nNodes, sNodeIdentifier=1):
    """
    Generates a list of node identifiers forming horizontal rows (y-direction) within a box core
    structure inside an ellipsoid.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    :param nNodes: Total number of nodes.
    :param sNodeIdentifier: Starting node identifier.
    :return: List of identifiers for box nodes grouped in rows.
    """
    boxRowNids = []
    eNodeIdentifier = nNodes + sNodeIdentifier - 1

    snid = sNodeIdentifier
    for _ in range(2):
        row = []
        nid = snid
        for _ in range(elementsCountNormal - 1):
            row.append([nid + i for i in range(elementsCountOblique - 1)])
            nid += elementsCountOblique - 1
        boxRowNids.append(row)
        snid = nid

    for e1 in [-2, -1]:
        row = []
        sNodeIdentifier = eNodeIdentifier - ((elementsCountOblique - 1) * (elementsCountNormal - 1)) * abs(e1) + 1
        for _ in range(elementsCountNormal - 1):
            row.append([sNodeIdentifier + i for i in range(elementsCountOblique - 1)])
            sNodeIdentifier += elementsCountOblique - 1
        boxRowNids.append(row)

    return boxRowNids


def getBoundaryNodeIdentifiers(elementsCountLateral, elementsCountNormal, elementsCountOblique, elementsCountAround,
                               elementsCountShell, sNodeIdentifier=1):
    """
    Generates a list of boundary node identifiers.
    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    :param elementsCountAround: Number of elements around the boundary.
    :param elementsCountShell:
    :param sNodeIdentifier: Starting node identifier (default is 1).
    :return: List of boundary node identifiers grouped by layers.
    """
    boundaryNids = []
    increment = (elementsCountNormal - 1) * (elementsCountOblique - 1)

    for e1 in range(elementsCountLateral - 1):
        boundaryNids.append([])
        nid = sNodeIdentifier
        for e3 in range(elementsCountNormal - 1):
            boundaryNids[e1].append([])
            if e3 == 0 or e3 == elementsCountNormal - 2:
                for e2 in range(elementsCountOblique - 1):
                    boundaryNids[e1][e3].append(nid)
                    nid += 1
            else:
                for e2 in range(2):
                    boundaryNids[e1][e3].append(nid)
                    increment_value = elementsCountOblique + (elementsCountShell * 2) \
                        if (e1 > 0 and e1 < elementsCountLateral - 2) else elementsCountOblique - 2
                    nid += increment_value

            if e1 > 0 and e1 < elementsCountLateral - 2:
                if e3 + 1 == 1 or e3 + 1 == elementsCountNormal - 2:
                    nid = boundaryNids[e1][e3][-1] + ((elementsCountOblique - 1) * (elementsCountShell + 1)) + 1
                else:
                    nid = boundaryNids[e1][e3][-1] + 1
            else:
                nid = boundaryNids[e1][e3][-1] + 1

        if e1 == 0:
            sNodeIdentifier += increment * (2 + elementsCountShell)
        elif 0 < e1 < elementsCountLateral - 3:
            sNodeIdentifier += increment + elementsCountAround
        else:
            sNodeIdentifier += increment * (elementsCountShell + 2) + elementsCountAround * (elementsCountShell + 1)

    return boundaryNids
