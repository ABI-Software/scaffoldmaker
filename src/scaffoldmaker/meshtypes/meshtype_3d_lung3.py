"""
Generates a lung scaffold by deforming a hemisphere.
"""
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
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
            "Human 1"]

    @classmethod
    def getDefaultOptions(cls, parameterSetName="Default"):
        useParameterSetName = "Human 1" if (parameterSetName == "Default") else parameterSetName
        isHuman = "Human" in useParameterSetName
        options = {
            "Number of elements lateral": 3,
            "Number of elements normal": 8,
            "Number of elements oblique": 8,
            "Number of elements shell": 0,
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
            "Refine": False,
            "Refine number of elements": 4,
        }
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Number of elements lateral",
            "Number of elements normal",
            "Number of elements oblique",
            "Number of elements shell",
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
            if options[dimension] <= 0.0:
                options[dimension] = 0.0
            elif options[dimension] > 1.0:
                options[dimension] = 1.0
        return dependentChanges

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        elementsCountLateral = options["Number of elements lateral"]
        elementsCountNormal = options["Number of elements normal"]
        elementsCountOblique = options["Number of elements oblique"]
        elementsCountShell = options["Number of elements shell"]
        elementsCountTransition = 1
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

        rangeOfRequiredElements = [
            [elementsCountLateral, 2 * elementsCountLateral],
            [0, elementsCountOblique],
            [0, elementsCountNormal]
        ]
        boxDerivatives = [1, 2, 3]
        sphereBoxDerivatives = [-boxDerivatives[0], boxDerivatives[1], boxDerivatives[2]]
        sphere_shape = SphereShape.SPHERE_SHAPE_FULL

        fieldmodule = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fieldmodule)

        mesh = fieldmodule.findMeshByDimension(3)
        boxGroup = AnnotationGroup(region, ("box group", ""))
        boxMeshGroup = boxGroup.getMeshGroup(mesh)
        transitionGroup = AnnotationGroup(region, ("transition group", ""))
        transitionMeshGroup = transitionGroup.getMeshGroup(mesh)
        meshGroups = [boxMeshGroup, transitionMeshGroup]
        annotationGroups = [boxGroup, transitionGroup]

        centre = [0.0, 0.0, 0.0]
        axis1 = [1.0, 0.0, 0.0]
        axis2 = [0.0, 1.0, 0.0]
        axis3 = [0.0, 0.0, 1.0]
        axes = [[disc_depth, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]
        elementsCountAcross = [2 * elementsCountLateral, elementsCountOblique, elementsCountNormal]
        shellProportion = 1.0
        sphere1 = SphereMesh(fieldmodule, coordinates, centre, axes, elementsCountAcross,
                             elementsCountShell, elementsCountTransition, shellProportion,
                             sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                             boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False,  meshGroups=meshGroups)

        rotate_scale_lung(fieldmodule, coordinates, disc_breadth, disc_height, oblique_slope_radians)

        if impression_depth_proportion > 0.0:
            form_mediastinal_surface(fieldmodule, coordinates, disc_breadth, disc_height, disc_depth,
                                     impression_breadth_proportion, impression_height_proportion,
                                     impression_depth_proportion, lateral_shear_rate)

        if (diaphragm_angle_radians != 0.0) and (diaphragm_proportion > 0.0):
            form_diaphragm_surface(fieldmodule, coordinates, diaphragm_angle_radians,
                                   (diaphragm_proportion - 0.5) * disc_height)

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


def rotate_scale_lung(fieldmodule, coordinates, disc_breadth, disc_height, oblique_slope_radians):
    """
    Rotate and scale lung disc to achieve oblique slope, height and breadth.
    Transformation tries to conserve original right angle between major and minor axes,
    and keep element sizes around to similar sizes.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param disc_breadth: Depth of lung ellipse.
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
    fieldassignment.assign()


def form_mediastinal_surface(fieldmodule, coordinates, disc_breadth, disc_height, disc_depth,
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
    fieldassignment.assign()


def form_diaphragm_surface(fieldmodule, coordinates, diaphragm_angle_radians, diaphragm_bend_z):
    """
    Form diaphragm surface by bending model up, phasing out to nothing at the apex.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param diaphragm_angle_radians: Angle of diaphragm from down mediastinal surface.
    :param diaphragm_bend_z: Z coordinate if centre of diaphragm bend.
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
    cos_new_theta = fieldmodule.createFieldCos(new_theta) * minus_1
    sin_new_theta = fieldmodule.createFieldSin(new_theta)
    new_x = r * cos_new_theta
    new_z = bend_z + r * sin_new_theta

    new_coordinates = fieldmodule.createFieldConcatenate([new_x, y, new_z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.assign()
