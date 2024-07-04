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
            "Depth": 0.8,
            "Height": 1.0,
            "Width": 0.4,
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
            "Depth",
            "Height",
            "Width",
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
        for dimension in ["Depth", "Height", "Width"]:
            if options[dimension] <= 0.0:
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
        depth = options["Depth"]
        height = options["Height"]
        width = options["Width"]
        obliqueSlopeRadians = math.radians(options["Oblique slope degrees"])

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
        axes = [[width, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]
        elementsCountAcross = [2 * elementsCountLateral, elementsCountOblique, elementsCountNormal]
        shellProportion = 1.0
        sphere1 = SphereMesh(fieldmodule, coordinates, centre, axes, elementsCountAcross,
                             elementsCountShell, elementsCountTransition, shellProportion,
                             sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                             boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False,  meshGroups=meshGroups)

        rotate_scale_lung(fieldmodule, coordinates, depth, height, obliqueSlopeRadians)

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


def rotate_scale_lung(fieldmodule, coordinates, depth, height, obliqueSlopeRadians):
    """
    Rotate and scale lung to achieve oblique slope and basic height and depth.
    Transformation tries to conserve original right angle between major and minor axes,
    and keep element sizes around to similar sizes.
    :param fieldmodule: Field module being worked with.
    :param coordinates: The coordinate field, initially circular in y-z plane.
    :param depth: Depth of lung ellipse.
    :param height: Height of lung ellipse.
    :param obliqueSlopeRadians: Slope of oblique fissure, where a positive angle tilts down ventrally.
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
    oblique_angle = fieldmodule.createFieldConstant(obliqueSlopeRadians)
    phi = safe_theta - oblique_angle
    cos_phi = fieldmodule.createFieldCos(phi)
    sin_phi = fieldmodule.createFieldSin(phi)

    scaling = 2.0 * (math.atan(depth / height) - 0.25 * math.pi)
    cc = fieldmodule.createFieldConstant(scaling)
    delta_phi = cc * cos_phi * sin_phi
    new_theta = phi + delta_phi
    cos_new_theta = fieldmodule.createFieldCos(new_theta)
    sin_new_theta = fieldmodule.createFieldSin(new_theta)
    depth_field = fieldmodule.createFieldConstant(depth)
    height_field = fieldmodule.createFieldConstant(height)
    new_y = r * cos_new_theta * depth_field
    new_z = r * sin_new_theta * height_field
    new_coordinates = fieldmodule.createFieldConcatenate([x, new_y, new_z])

    fieldassignment = coordinates.createFieldassignment(new_coordinates)
    fieldassignment.assign()
