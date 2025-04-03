"""
Generates a lung scaffold by deforming a hemisphere.
"""
from cmlibs.maths.vectorops import magnitude, set_magnitude, mult, sub
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc import field
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.lung_terms import get_lung_term
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.derivativemoothing import DerivativeSmoothing
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine, sampleCubicHermiteCurvesSmooth, \
    sampleCubicHermiteCurves, DerivativeScalingMode
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
import math

from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters


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
            "Human 2",
            "Human 3",
            "Ellipsoid",
            "Teardrop"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        useParameterSetName = "Human 1" if (parameterSetName == "Default") else parameterSetName
        isHuman = "Human" in useParameterSetName
        if parameterSetName in ["Human 1", "Default"]:
            options = {
                "Half sphere": True,
                "Left lung only": True,
                "Scale": False,
                "Scale factor": 1.0,
                "Use serendipity": False,
                "Number of elements lateral": 4,
                "Number of elements normal": 8,
                "Number of elements oblique": 8,
                "Number of elements shell": 0,
                "Left-right lung spacing": 0.0,
                "Apex edge sharpness factor": 0.0,
                "Ventral edge sharpness factor": 0.0,
                "Left-right apex medial shear displacement": 0.0,
                "Left-right apex ventral shear displacement": 0.0,
                "Left-right base shear displacement": 0.0,
                "Diaphragm angle degrees": 60.0,
                "Diaphragm proportion": 0.3,
                "Disc breadth": 0.7,
                "Disc height": 1.0,
                "Disc depth": 0.35,
                "Impression breadth proportion": 0.8,
                "Impression height proportion": 0.8,
                "Impression depth proportion": 0.25,
                "Lateral shear rate": 0.75,
                "Oblique slope degrees": 45.0,
                "Medial curvature": 0.0,
                "Medial curvature bias": 0.0,
                "Ventral-medial rotation degrees": 0.0,
                "Refine": False,
                "Refine number of elements": 4,
            }
        elif parameterSetName == "Human 2":
            options = {
                "Half sphere": False,
                "Left lung only": False,
                "Scale": False,
                "Scale factor": 1.0,
                "Use serendipity": False,
                "Number of elements lateral": 4,
                "Number of elements normal": 8,
                "Number of elements oblique": 8,
                "Number of elements shell": 0,
                "Left-right lung spacing": 0.5,
                "Apex edge sharpness factor": 0.2,
                "Ventral edge sharpness factor": 0.8,
                "Left-right apex medial shear displacement": 0.2,
                "Left-right apex ventral shear displacement": -0.1,
                "Left-right base shear displacement": 0.0,
                "Diaphragm angle degrees": 30.0,
                "Diaphragm proportion": 0.3,
                "Disc breadth": 0.8,
                "Disc height": 1.1,
                "Disc depth": 0.25,
                "Impression breadth proportion": 0.8,
                "Impression height proportion": 0.8,
                "Impression depth proportion": 0.8,
                "Lateral shear rate": -2.0,
                "Oblique slope degrees": 45.0,
                "Medial curvature": 1.0,
                "Medial curvature bias": 0.0,
                "Ventral-medial rotation degrees": 0.0,
                "Refine": False,
                "Refine number of elements": 4,
            }
        elif parameterSetName == "Human 3":
            options = {
                "Half sphere": False,
                "Left lung only": False,
                "Scale": False,
                "Scale factor": 1.0,
                "Use serendipity": False,
                "Number of elements lateral": 6,
                "Number of elements normal": 8,
                "Number of elements oblique": 8,
                "Number of elements shell": 0,
                "Left-right lung spacing": 0.4,
                "Apex edge sharpness factor": 0.3,
                "Ventral edge sharpness factor": 0.9,
                "Left-right apex medial shear displacement": 0.2,
                "Left-right apex ventral shear displacement": -0.1,
                "Left-right base shear displacement": 0.0,
                "Diaphragm angle degrees": 30.0,
                "Diaphragm proportion": 0.3,
                "Disc breadth": 0.8,
                "Disc height": 1.2,
                "Disc depth": 0.25,
                "Impression breadth proportion": 0.8,
                "Impression height proportion": 0.8,
                "Impression depth proportion": 1.0,
                "Lateral shear rate": -1.0,
                "Oblique slope degrees": 45.0,
                "Medial curvature": 1.0,
                "Medial curvature bias": 0.0,
                "Ventral-medial rotation degrees": 5.0,
                "Refine": False,
                "Refine number of elements": 4,
            }
        elif parameterSetName == "Ellipsoid":
            options = {
                "Half sphere": False,
                "Left lung only": False,
                "Scale": False,
                "Scale factor": 1.0,
                "Use serendipity": False,
                "Number of elements lateral": 4,
                "Number of elements normal": 8,
                "Number of elements oblique": 8,
                "Number of elements shell": 0,
                "Left-right lung spacing": 1.0,
                "Apex edge sharpness factor": 0.0,
                "Ventral edge sharpness factor": 0.0,
                "Left-right apex medial shear displacement": 0.0,
                "Left-right apex ventral shear displacement": 0.0,
                "Left-right base shear displacement": 0.0,
                "Diaphragm angle degrees": 0.0,
                "Diaphragm proportion": 0.0,
                "Disc breadth": 0.8,
                "Disc height": 1.1,
                "Disc depth": 0.25,
                "Impression breadth proportion": 0.0,
                "Impression height proportion": 0.0,
                "Impression depth proportion": 0.0,
                "Lateral shear rate": 0.0,
                "Oblique slope degrees": 45.0,
                "Medial curvature": 0.0,
                "Medial curvature bias": 0.0,
                "Ventral-medial rotation degrees": 0.0,
                "Refine": False,
                "Refine number of elements": 4,
            }
        elif parameterSetName == "Teardrop":
            options = {
                "Half sphere": False,
                "Left lung only": True,
                "Scale": False,
                "Scale factor": 1.0,
                "Use serendipity": False,
                "Number of elements lateral": 4,
                "Number of elements normal": 4,
                "Number of elements oblique": 4,
                "Number of elements shell": 0,
                "Left-right lung spacing": 1.0,
                "Apex edge sharpness factor": 0.0,
                "Ventral edge sharpness factor": 0.95,
                "Left-right apex medial shear displacement": 0.0,
                "Left-right apex ventral shear displacement": 0.0,
                "Left-right base shear displacement": 0.0,
                "Diaphragm angle degrees": 0.0,
                "Diaphragm proportion": 0.0,
                "Disc breadth": 0.8,
                "Disc height": 1.1,
                "Disc depth": 0.25,
                "Impression breadth proportion": 0.0,
                "Impression height proportion": 0.0,
                "Impression depth proportion": 0.0,
                "Lateral shear rate": 0.0,
                "Oblique slope degrees": 45.0,
                # "Medial proportion": 0.0,
                "Medial curvature": 0.0,
                "Medial curvature bias": 0.0,
                "Ventral-medial rotation degrees": 0.0,
                "Refine": False,
                "Refine number of elements": 4,
            }
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Left lung only",
            "Scale",
            "Scale factor",
            "Use serendipity",
            "Number of elements lateral",
            "Number of elements normal",
            "Number of elements oblique",
            "Number of elements shell",
            "Left-right lung spacing",
            "Apex edge sharpness factor",
            "Ventral edge sharpness factor",
            "Left-right apex medial shear displacement",
            "Left-right apex ventral shear displacement",
            "Left-right base shear displacement",
            "Diaphragm angle degrees",
            "Diaphragm proportion",
            "Disc breadth",
            "Disc height",
            "Disc depth",
            "Impression breadth proportion",
            "Impression height proportion",
            "Impression depth proportion",
            "Lateral shear rate",
            "Oblique slope degrees",
            "Medial curvature",
            "Medial curvature bias",
            "Ventral-medial rotation degrees",
            "Refine",
            "Refine number of elements"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependentChanges = False
        if options["Number of elements lateral"] < 2:
            options["Number of elements lateral"] = 2
        for key in [
            "Number of elements normal",
            "Number of elements oblique"]:
            if options[key] < 4:
                options[key] = 4
            elif options[key] % 2:
                options[key] += 1

        maxShellElements = min(options["Number of elements lateral"],
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
            "Medial curvature bias"
        ]:
            if options[dimension] < 0.0:
                options[dimension] = 0.0
            elif options[dimension] > 1.0:
                options[dimension] = 1.0

        if options['Ventral-medial rotation degrees'] < -90.0:
            options['Ventral-medial rotation degrees'] = -90.0
        elif options['Ventral-medial rotation degrees'] > 90.0:
            options['Ventral-medial rotation degrees'] = 90.0

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
        showLeftLungOnly = options["Left lung only"]
        useHalfSphere = options["Half sphere"]
        useSerendipity = options["Use serendipity"]

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
        forwardLeftRightBase = options["Left-right base shear displacement"]
        diaphragm_angle_radians = math.radians(options["Diaphragm angle degrees"])
        diaphragm_proportion = options["Diaphragm proportion"]
        disc_breadth = options["Disc breadth"]
        disc_height = options["Disc height"]
        disc_depth = options["Disc depth"]
        impression_breadth_proportion = options["Impression breadth proportion"]
        impression_height_proportion = options["Impression height proportion"]
        impression_depth_proportion = options["Impression depth proportion"]
        lateral_shear_rate = options["Lateral shear rate"]
        oblique_slope_radians = math.radians(options["Oblique slope degrees"])
        leftLungMedialCurvature = options["Medial curvature"]
        lungMedialCurvatureBias = options["Medial curvature bias"]
        rotateLeftLung = options["Ventral-medial rotation degrees"]

        if useHalfSphere:
            rangeOfRequiredElements = [
                [elementsCountLateral, 2 * elementsCountLateral],
                [0, elementsCountOblique],
                [0, elementsCountNormal]
            ]
        else:
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

        boxGroup = AnnotationGroup(region, ("box group", ""))
        boxMeshGroup = boxGroup.getMeshGroup(mesh)
        boxNodesetGroup = boxGroup.getNodesetGroup(nodes)

        transitionGroup = AnnotationGroup(region, ("transition group", ""))
        transitionMeshGroup = transitionGroup.getMeshGroup(mesh)
        transitionNodesetGroup = transitionGroup.getNodesetGroup(nodes)

        lungGroup = AnnotationGroup(region, get_lung_term("lung"))
        leftLungGroup = AnnotationGroup(region, get_lung_term("left lung"))
        annotationGroups = [leftLungGroup, lungGroup]
        lungMeshGroup = lungGroup.getMeshGroup(mesh)
        leftLungMeshGroup = leftLungGroup.getMeshGroup(mesh)

        rightLungGroup = AnnotationGroup(region, get_lung_term("right lung"))
        rightLungMeshGroup = rightLungGroup.getMeshGroup(mesh)
        annotationGroups.append(rightLungGroup)

        leftLungNodesetGroup = leftLungGroup.getNodesetGroup(nodes)
        rightLungNodesetGroup = rightLungGroup.getNodesetGroup(nodes)

        centre = [0.0, 0.0, 0.0]
        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]
        # axes = [[disc_depth, 0.0, 0.0], [0.0, disc_breadth * 0.5, 0.0], [0.0, 0.0, disc_height * 0.5]]
        axes = [[disc_depth, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]

        if useHalfSphere:
            elementsCountAcross = [2 * elementsCountLateral, elementsCountOblique, elementsCountNormal]
        else:
            elementsCountAcross = [elementsCountLateral, elementsCountOblique, elementsCountNormal]
        shellProportion = 1.0

        leftLung, rightLung = 0, 1
        lungs = [leftLung] if showLeftLungOnly else [leftLung, rightLung]
        for i in lungs:
            if i == leftLung:
                meshGroups = [boxMeshGroup, transitionMeshGroup, leftLungMeshGroup]
            else:
                meshGroups = [boxMeshGroup, transitionMeshGroup, rightLungMeshGroup]

            sphere = SphereMesh(fieldmodule, coordinates, centre, axes, elementsCountAcross,
                                elementsCountShell, elementsCountTransition, shellProportion,
                                sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                                boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False,
                                useSerendipity=useSerendipity, meshGroups=meshGroups)

        annotationGroups = [boxGroup, transitionGroup, leftLungGroup, rightLungGroup]

        if options["Scale"]:
            scale_factor = options["Scale factor"]
            for i in lungs:
                isLeft = True if i == leftLung else False
                lungNodeset = leftLungNodesetGroup if isLeft else rightLungNodesetGroup
                scaleNodes2(region, fieldmodule, coordinates, lungNodeset, disc_breadth, disc_height, disc_depth,
                           scale_factor)
                if elementsCountLateral > 4:
                    resampleShellNodes(fieldmodule, coordinates, lungNodeset, elementsCountLateral,
                                       elementsCountNormal, elementsCountOblique)

        for i in lungs:
            isLeft = True if i == leftLung else False
            lungNodeset = leftLungNodesetGroup if isLeft else rightLungNodesetGroup
            halfBreadth = disc_breadth * 0.5
            height = disc_height
            spacing = -lungSpacing if i == leftLung else lungSpacing
            bendZ = (diaphragm_proportion - 0.5) * disc_height
            lungMedialCurvature = -leftLungMedialCurvature if isLeft else leftLungMedialCurvature
            apexMedialDisplacement = leftApexMedialDisplacement if isLeft else -leftApexMedialDisplacement
            rotateLungAngle = rotateLeftLung if isLeft else -rotateLeftLung

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
                form_diaphragm_surface(fieldmodule, coordinates, lungNodeset, diaphragm_angle_radians, bendZ, isLeft=isLeft)

            dorsalVentralXi = getDorsalVentralXiField(fieldmodule, coordinates, halfBreadth)
            if lungMedialCurvature != 0:
                bendingAroundZAxis(lungMedialCurvature, fieldmodule, coordinates, lungNodeset,
                                   stationaryPointXY=[0.0, 0.0],
                                   bias=lungMedialCurvatureBias,
                                   dorsalVentralXi=dorsalVentralXi)

            if rotateLungAngle != 0.0:
                rotateLungs(rotateLungAngle, fieldmodule, coordinates, lungNodeset)

            translateLungLocation(fieldmodule, coordinates, lungNodeset, spacing, bendZ)

            if apexMedialDisplacement != 0.0:
                medialShearRadian = math.atan(apexMedialDisplacement / height)
                tiltLungs(medialShearRadian, 0, 0, 0, fieldmodule, coordinates, lungNodeset)

            if forwardLeftRightApex != 0.0:
                ventralShearRadian = math.atan(forwardLeftRightApex / height)
                tiltLungs(0, ventralShearRadian, 0, 0, fieldmodule, coordinates, lungNodeset, isApex=True)

            if forwardLeftRightBase != 0.0:
                baseShearRadian = math.atan(forwardLeftRightBase / height)
                tiltLungs(0, baseShearRadian, 0, 0.05, fieldmodule, coordinates, lungNodeset, isApex=False)

            if options["Scale"]:
                smoothD3ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountLateral,
                                             elementsCountNormal, elementsCountOblique)
                smoothD1ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountNormal, elementsCountOblique)
                smoothD2ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountNormal, elementsCountOblique)

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


    # @classmethod
    # def defineFaceAnnotations(cls, region, options, annotationGroups):
    #     """
    #     Add face annotation groups from the highest dimension mesh.
    #     Must have defined faces and added subelements for highest dimension groups.
    #     :param region: Zinc region containing model.
    #     :param options: Dict containing options. See getDefaultOptions().
    #     :param annotationGroups: List of annotation groups for top-level elements.
    #     New face annotation groups are appended to this list.
    #     """
        # # numberOfLeftLung = options['Number of left lung lobes']
        # numberOfLeftLung = 2
        # # hasAccessoryLobe = options['Accessory lobe']
        # # openFissures = options['Open fissures']
        #
        # # create fissure groups
        # fm = region.getFieldmodule()
        # mesh1d = fm.findMeshByDimension(1)
        # mesh2d = fm.findMeshByDimension(2)
        #
        #
        #
        # # 1D Annotation
        # is_exterior = fm.createFieldIsExterior()
        # is_xi1_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0)
        # is_xi1_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1)
        # is_xi2_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0)
        # is_xi2_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1)
        # is_xi3_0 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0)
        # is_xi3_1 = fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1)
        #
        # # 1D edge markers
        # mediastanumTerms = ["left lung", "right lung"]
        # for mediastanumTerm in mediastanumTerms:
        #     mediastanumGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(mediastanumTerm))
        #     is_anteriorBorderGroup = mediastanumGroup.getGroup()
        #     is_mediastanumGroup_exterior = fm.createFieldAnd(is_anteriorBorderGroup, is_exterior)
        #     is_mediastanumGroup_exterior_xi1_0 = fm.createFieldAnd(is_mediastanumGroup_exterior, is_xi1_0)
        #     is_mediastanumGroup_exterior_xi1_01 = fm.createFieldAnd(is_mediastanumGroup_exterior_xi1_0, is_xi1_1)
        #     check1 = fm.createFieldAnd(is_mediastanumGroup_exterior, is_xi1_1)
        #     check2 = fm.createFieldAnd(is_mediastanumGroup_exterior, is_xi2_0)
        #     check3 = fm.createFieldAnd(is_mediastanumGroup_exterior, is_xi2_1)
        #     check4 = fm.createFieldAnd(is_mediastanumGroup_exterior, is_xi3_0)
        #     check5 = fm.createFieldAnd(is_mediastanumGroup_exterior, is_xi3_1)
        #     is_ridge = fm.createFieldAnd(is_mediastanumGroup_exterior_xi1_01, fm.createFieldNot(is_xi3_0))
        #     # if openFissures:
        #     #     is_ridge = fm.createFieldAnd(is_ridge, fm.createFieldNot(is_xi3_1))
        #     borderTerm = "anterior border of left lung" if "left lung" in mediastanumTerm else "anterior border of right lung"
        #     anteriorBorderLeftGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_lung_term(borderTerm))
        #     anteriorBorderLeftGroup.getMeshGroup(mesh1d).addElementsConditional(is_ridge)
        #     # anteriorBorderLeftGroup.getMeshGroup(mesh1d).addElementsConditional(check5)
        #
        #     # annotationGroups.remove(mediastanumGroup)


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


def form_diaphragm_surface(fieldmodule, coordinates, lungNodesetGroup, diaphragm_angle_radians, diaphragm_bend_z, isLeft=True):
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


def rotateLungs(rotateZ, fm, coordinates, lungNodesetGroup):
    """
    Rotate the lung mesh at the center of the elements about z-axis
    :param rotateZ:
    :param fm: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :return: None
    """
    # FieldConstant - Matrix = [   x1,    x4,    x7,
    #                              x2,    x5,    x8,
    #                              x3,    x6,    x9]

    rotateZ = -rotateZ / 180 * math.pi  # negative value due to right handed rule
    rotateZMatrix = fm.createFieldConstant([math.cos(rotateZ), math.sin(rotateZ), 0.0,
                                            -math.sin(rotateZ), math.cos(rotateZ), 0.0,
                                            0.0, 0.0, 1.0])
    translate_coordinates = fm.createFieldMatrixMultiply(3, rotateZMatrix, coordinates)

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
    k1 = sharpeningFactor / (halfBreadth * 1.75) if isApex else -sharpeningFactor / (halfBreadth * 1.75)
    # scale = fm.createFieldConstant([0.0, 0.0, k1]) if isApex else fm.createFieldConstant([0.0, k1, 0.0])
    scale = fm.createFieldConstant([0.0, -k1, 0.0]) if isApex else fm.createFieldConstant([0.0, k1, 0.0])


    scaleFunction = fm.createFieldMultiply(origin, scale)
    constant = fm.createFieldConstant([0.0, 1.0, 1.0])
    constantFunction = fm.createFieldAdd(scaleFunction, constant)
    # componentIndices = [3, 2, 3] if isApex else [2, 3, 3]
    componentIndices = [3, 3, 2] if isApex else [2, 3, 3]


    transformation_matrix = fm.createFieldComponent(constantFunction, componentIndices)
    taper_coordinates = fm.createFieldMultiply(origin, transformation_matrix)
    translate_coordinates = fm.createFieldSubtract(taper_coordinates, offset)

    if isApex:
        y = fm.createFieldComponent(coordinates, 2)
        isRightSide = fm.createFieldGreaterThan(y, fm.createFieldConstant(-0.05))
        translate_coordinates = fm.createFieldIf(isRightSide,coordinates, translate_coordinates)

    fieldassignment = coordinates.createFieldassignment(translate_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def tiltLungs(tiltApex_xAxis, tiltApex_yAxis, tiltDiap_yAxis, tiltDiap_xAxis, fm, coordinates, lungNodesetGroup, isApex=False):
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


def redistributeBoxNodes(fieldmodule, coordinates, lungNodesetGroup, disc_breadth, disc_height, scale_factor):
    """
    Redistribute lung disc nodes so that the outermost nodes remain in place while inner nodes are scaled outward.
    The scaling factor is 1 at the outer boundary and increases towards the center, causing nodes to shift outward.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in the y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param disc_breadth: Breadth of lung ellipse.
    :param disc_height: Height of lung ellipse.
    :param scale_factor: Maximum scaling effect at the center (should be > 1 for expansion).
    """
    x = fieldmodule.createFieldComponent(coordinates, 1)
    y = fieldmodule.createFieldComponent(coordinates, 2)
    z = fieldmodule.createFieldComponent(coordinates, 3)
    yz = fieldmodule.createFieldComponent(coordinates, [2, 3])
    r = fieldmodule.createFieldMagnitude(yz)
    theta = fieldmodule.createFieldAtan2(z, y)
    one = fieldmodule.createFieldConstant(1.0)

    cos_theta = fieldmodule.createFieldCos(theta)
    sin_theta = fieldmodule.createFieldSin(theta)

    cos2 = fieldmodule.createFieldMultiply(cos_theta, cos_theta)
    sin2 = fieldmodule.createFieldMultiply(sin_theta, sin_theta)
    ab = fieldmodule.createFieldConstant(disc_breadth * disc_height)
    aa = fieldmodule.createFieldConstant(disc_breadth ** 2)
    bb = fieldmodule.createFieldConstant(disc_height ** 2)

    denom = fieldmodule.createFieldSqrt((aa * sin2 + bb * cos2))
    max_r = ab / denom * fieldmodule.createFieldConstant(0.9)

    scaling = fieldmodule.createFieldConstant(scale_factor) + (one - r / max_r) * fieldmodule.createFieldConstant(0.9)

    new_y = y * scaling
    new_z = z * scaling

    new_coordinates = fieldmodule.createFieldConcatenate([x, new_y, new_z])

    # val1 = fieldmodule.createFieldConstant(one - r / max_r)
    # val2 = fieldmodule.createFieldConstant(0.9)
    # isGreater2 = fieldmodule.createFieldGreaterThan(fieldmodule.createFieldAbs(x), fieldmodule.createFieldConstant(0.25) * fieldmodule.createFieldConstant(0.9))
    # isGreater1 = fieldmodule.createFieldGreaterThan((r / max_r), val2)

    # isGreater = fieldmodule.createFieldAnd(isGreater1, isGreater2)
    # isBox = fieldmodule.createFieldIf(isGreater, )
    # new_coordinates = fieldmodule.createFieldIf(isGreater, coordinates, new_coordinates)

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()


def scaleNodes(region, fieldmodule, coordinates, lungNodesetGroup, disc_breadth, disc_height, disc_depth, scale_factor):
    """
    Redistribute lung nodes so that the outermost nodes remain in place while inner nodes are scaled outward.
    The scaling factor is 1 at the outer boundary and increases towards the center, causing nodes to shift outward.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in the y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param disc_depth: Depth of lung ellipse.
    :param disc_breadth: Breadth of lung ellipse.
    :param disc_height: Height of lung ellipse.
    :param scale_factor: Maximum scaling effect at the center (should be > 1 for expansion).
    """
    x = fieldmodule.createFieldComponent(coordinates, 1)
    y = fieldmodule.createFieldComponent(coordinates, 2)
    z = fieldmodule.createFieldComponent(coordinates, 3)
    one = fieldmodule.createFieldConstant(1.0)
    zero = fieldmodule.createFieldConstant(0.0)

    # calculate r and theta
    r = fieldmodule.createFieldSqrt((y * y) + (z * z))
    theta = fieldmodule.createFieldAtan2(z, y)
    r_zero = fieldmodule.createFieldEqualTo(r, zero)
    safe_theta = fieldmodule.createFieldIf(r_zero, zero, theta)

    # calculate new semi-axes
    a = fieldmodule.createFieldConstant(disc_depth)
    term = fieldmodule.createFieldSqrt(
        one - fieldmodule.createFieldMultiply(x, x) / fieldmodule.createFieldMultiply(a, a))

    b = fieldmodule.createFieldConstant(disc_breadth * 0.5)
    c = fieldmodule.createFieldConstant(disc_height * 0.5)
    b_new = b * term
    c_new = c * term

    # calculate maximum distance from the centre to the ellipsoid surface
    factor = fieldmodule.createFieldSqrt(one - ((x * x) / (a * a)))
    y_surface = b_new * factor * fieldmodule.createFieldCos(safe_theta)
    z_surface = c_new * factor * fieldmodule.createFieldSin(safe_theta)

    max_r = fieldmodule.createFieldSqrt(x * x + y_surface * y_surface + z_surface * z_surface)

    # threshold
    y_abs = fieldmodule.createFieldAbs(y)
    z_abs = fieldmodule.createFieldAbs(z)

    r_threshold = max_r * fieldmodule.createFieldConstant(0.9)
    # r_threshold = b_new

    isGreaterThanThreshold = fieldmodule.createFieldGreaterThan(r, r_threshold)
    # isGreaterThanThreshold = fieldmodule.createFieldGreaterThan(y_abs, r_threshold)

    isGreaterThanApex = fieldmodule.createFieldGreaterThan(fieldmodule.createFieldAbs(x), a * fieldmodule.createFieldConstant(0.8))

    isConditional = fieldmodule.createFieldOr(isGreaterThanThreshold, isGreaterThanApex)

    # calculate scaling factor
    scaling = fieldmodule.createFieldConstant(scale_factor) + (one - r / max_r) * fieldmodule.createFieldConstant(0.9)
    # scaling = (fieldmodule.createFieldConstant(scale_factor) + (one - r / max_r) +
    #            (fieldmodule.createFieldConstant(0.1)) * (fieldmodule.createFieldAbs(x) / a))

    # scaling = ((y - y_surface * fieldmodule.createFieldConstant(0.5)) * (y - y_surface * fieldmodule.createFieldConstant(0.5))
    #            + fieldmodule.createFieldConstant(scale_factor))

    new_y = fieldmodule.createFieldIf(isConditional, y * one, y * scaling)
    new_z = fieldmodule.createFieldIf(isConditional, z * one, z * scaling)

    # new_y = y * scaling
    # new_z = z * scaling

    new_coordinates = fieldmodule.createFieldConcatenate([x, new_y, new_z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.setNodeset(lungNodesetGroup)
    fieldassignment.assign()

    # smoothing = DerivativeSmoothing(region, coordinates, scalingMode=DerivativeScalingMode.HARMONIC_MEAN)
    # smoothing.smooth(updateDirections=True)


def scaleNodes2(region, fieldmodule, coordinates, lungNodesetGroup, disc_breadth, disc_height, disc_depth, scale_factor):
    """
    Redistribute lung nodes so that the outermost nodes remain in place while inner nodes are scaled outward.
    The scaling factor is 1 at the outer boundary and increases towards the center, causing nodes to shift outward.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in the y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param disc_depth: Depth of lung ellipse.
    :param disc_breadth: Breadth of lung ellipse.
    :param disc_height: Height of lung ellipse.
    :param scale_factor: Maximum scaling effect at the center (should be > 1 for expansion).
    """
    x = fieldmodule.createFieldComponent(coordinates, 1)
    y = fieldmodule.createFieldComponent(coordinates, 2)
    y_abs = fieldmodule.createFieldAbs(y)
    z = fieldmodule.createFieldComponent(coordinates, 3)

    zero = fieldmodule.createFieldConstant(0.0)
    one = fieldmodule.createFieldConstant(1.0)
    four = fieldmodule.createFieldConstant(4.0)

    tol = 0.01

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

    # smoothing = DerivativeSmoothing(region, coordinates, scalingMode=DerivativeScalingMode.HARMONIC_MEAN)
    # smoothing.smooth(updateDirections=True)


def resampleShellNodes(fieldmodule, coordinates, lungNodeset, elementsCountLateral, elementsCountNormal, elementsCountOblique):
    """
    Resample shell nodes around the ellipsoid in the yz plane so that the shell nodes are aligned linearly.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in the y-z plane.
    :param lungNodesetGroup: Zinc NodesetGroup containing nodes to transform.
    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    """
    setValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                      Node.VALUE_LABEL_D_DS3]
    nodes = lungNodeset
    nodeParameters = get_nodeset_field_parameters(nodes, coordinates, setValueLabels)[1]
    sNodeIdentifier = nodeParameters[0][0]
    nodeParametersDict = dict(nodeParameters)
    elementsCountAround = getElementsCountAround(elementsCountNormal, elementsCountOblique)
    boundaryNids = getBoundaryNodeIdentifiers(elementsCountLateral, elementsCountNormal, elementsCountOblique,
                                            elementsCountAround, sNodeIdentifier=sNodeIdentifier)

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

            # tx, td1 = sampleCubicHermiteCurves([bx[0], bx[-1]], [bd1[0], bd1[-1]], elementsCountLateral - 2, arcLengthDerivatives=True)[:2]
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
                midRange = range(normalMidIndex - 1, 0, -1) if e2 == 0 else range(normalMidIndex + 1, elementsCountNormal - 2)
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

                tx, td = sampleCubicHermiteCurves(bx, bd, len(bx) - 1, elementLengthStartEndRatio=0.7, arcLengthDerivatives=True)[:2]

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

                # tx, td = sampleCubicHermiteCurves([bx[0], bx[-1]], [bd[0], bd[-1]], elementsCountLateral - 2, arcLengthDerivatives=True)[:2]
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
            midRange = range(e3 + 1, obliqueMidIndex) if e3 == 0 else range(elementsCountOblique - 3, obliqueMidIndex, -1)
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

            tx, td = sampleCubicHermiteCurves(bx, bd, len(bx) - 1, elementLengthStartEndRatio=1.3, arcLengthDerivatives=True)[:2]

            for n, id in enumerate(snids):
                n2, dir = n2List[n], dirList[n]
                node = nodes.findNodeByIdentifier(id)
                fieldcache.setNode(node)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, tx[n])
                label = Node.VALUE_LABEL_D_DS1 if n2 == 1 else Node.VALUE_LABEL_D_DS2
                coordinates.setNodeParameters(fieldcache, -1, label, 1, mult(td[n], dir))


def smoothD3ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountLateral, elementsCountNormal, elementsCountOblique):
    """
    Smooth derivatives of nodes connecting the box mesh to the shell mesh.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field.
    :param lungNodeset: Zinc NodesetGroup containing nodes to transform.
    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    """
    setValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                      Node.VALUE_LABEL_D_DS3]
    nodes = lungNodeset
    nodeParameters = get_nodeset_field_parameters(nodes, coordinates, setValueLabels)[1]
    sNodeIdentifier = nodeParameters[0][0]
    nodeParametersDict = dict(nodeParameters)
    elementsCountAround = getElementsCountAround(elementsCountNormal, elementsCountOblique)

    nNodes = len(nodeParameters)
    nodeIds = getBoxRowNodeIdentifiers(elementsCountNormal, elementsCountOblique, nNodes, sNodeIdentifier)
    boundaryNodeIds = getBoundaryNodeIdentifiers(elementsCountLateral, elementsCountNormal, elementsCountOblique,elementsCountAround, sNodeIdentifier)

    for e1 in [0, -1]:
        e1p = e1 + 1 if e1 == 0 else e1 - 1
        for e2 in range(elementsCountNormal - 1):
            for e3 in range(elementsCountOblique - 1):
                sid = nodeIds[e1][e2][e3]
                bid = nodeIds[e1p][e2][e3]
                bx = nodeParametersDict[bid][0][0]
                sx = nodeParametersDict[sid][0][0]

                shellFactor = 1.0 # 1.0 / elementsCountRim
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


def smoothD1ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountNormal, elementsCountOblique):
    """
    Smooth derivatives of nodes connecting the box mesh to the shell mesh.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field.
    :param lungNodeset: Zinc NodesetGroup containing nodes to transform.
    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    """
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


def smoothD2ShellNodeDerivatives(fieldmodule, coordinates, lungNodeset, elementsCountNormal, elementsCountOblique):
    """
    Smooth derivatives of nodes connecting the box mesh to the shell mesh.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field.
    :param lungNodeset: Zinc NodesetGroup containing nodes to transform.
    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    """
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
                    n2 = 1 if (e2 == 0 and e3 < midObliqueIndex) or (e2 == elementsCountNormal - 2 and e3 > midObliqueIndex) else 2
                else:
                    n2 = 1 if (e2 == 0 and e3 > midObliqueIndex) or (e2 == elementsCountNormal - 2 and e3 < midObliqueIndex) else 2
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

            # snids = nodeIds[e1][e2]
            # sx, sd1 = [], []
            # if e1 == -1:
            #     dir = 1 if e2 < midNormalIndex else -1
            # else:
            #     dir = -1 if e2 < midNormalIndex else 1
            # for id in snids:
            #     sx.append(nodeParametersDict[id][0][0])
            #     sd1.append(mult(nodeParametersDict[id][1][0], dir))
            #
            # sd1 = smoothCubicHermiteDerivativesLine(sx, sd1, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
            #
            # fieldcache = fieldmodule.createFieldcache()
            # for n, id in enumerate(snids):
            #     node = nodes.findNodeByIdentifier(id)
            #     fieldcache.setNode(node)
            #     coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
            #     coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, mult(sd1[n], dir))


def getElementsCountAround(elementsCountNormal, elementsCountOblique):
    return (elementsCountNormal - 1) * 2 + (elementsCountOblique - 3) * 2


# def getBoxColumnNodeIdentifiers(elementsCountLateral, elementsCountNormal, elementsCountOblique):
#     """
#     Generates a list of node identifiers forming vertical columns (z-direction) within a box core
#     structure inside an ellipsoid.
#     :param elementsCountLateral: Number of elements along the lateral direction.
#     :param elementsCountNormal: Number of elements along the normal direction.
#     :param elementsCountOblique: Number of elements along the oblique direction.
#     :return: List of identifiers for box nodes grouped in columns.
#     """
#     boxColumnNids = []
#     sNodeIdentifier = 1
#     midLateral = elementsCountLateral // 2
#     midNormal = (elementsCountNormal - 2) // 2
#
#     for e1 in range(elementsCountLateral + 1):
#         increment = elementsCountOblique if e1 == midLateral else elementsCountOblique - 1
#         boxColumnNids.append([])
#         for e2 in range(elementsCountOblique - 1):
#             nids = []
#             nid = sNodeIdentifier
#             for e3 in range(elementsCountNormal - 1):
#                 nids.append(nid)
#                 nid += (increment + 1) if e1 == midLateral and e3 in {midNormal - 1, midNormal} else increment
#             boxColumnNids[e1].append(nids)
#             sNodeIdentifier += 1
#         sNodeIdentifier = nids[-1] + (elementsCountOblique if e1 in {midLateral, midLateral - 1} else 1)
#
#     return boxColumnNids


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


def getBoundaryNodeIdentifiers(elementsCountLateral, elementsCountNormal, elementsCountOblique, elementsCountAround, sNodeIdentifier=1):
    """
    Generates a list of boundary node identifiers.

    :param elementsCountLateral: Number of elements along the lateral direction.
    :param elementsCountNormal: Number of elements along the normal direction.
    :param elementsCountOblique: Number of elements along the oblique direction.
    :param elementsCountAround: Number of elements around the boundary.
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
                    increment_value = elementsCountOblique if (e1 > 0 and e1 < elementsCountLateral - 2) else elementsCountOblique - 2
                    nid += increment_value

            if e1 > 0 and e1 < elementsCountLateral - 2:
                if e3 + 1 == 1 or e3 + 1 == elementsCountNormal - 2:
                    nid = boundaryNids[e1][e3][-1] + elementsCountOblique
                else:
                    nid = boundaryNids[e1][e3][-1] + 1
            else:
                nid = boundaryNids[e1][e3][-1] + 1

        if e1 == 0:
            sNodeIdentifier += increment * 2
        elif 0 < e1 < elementsCountLateral - 3:
            sNodeIdentifier += increment + elementsCountAround
        else:
            sNodeIdentifier += increment * 2 + elementsCountAround

    return boundaryNids
