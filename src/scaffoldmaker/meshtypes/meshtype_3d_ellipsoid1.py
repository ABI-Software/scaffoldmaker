"""
Generates a solid ellipsoid of hexahedral elements.
"""
import math
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.ellipsoidmesh import EllipsoidMesh
from scaffoldmaker.utils.meshrefinement import MeshRefinement


class MeshType_3d_ellipsoid1(Scaffold_base):
    """
    Generates a solid ellipsoid of hexahedral elements.
    """

    @classmethod
    def getName(cls):
        return "3D Ellipsoid 1"

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {
            "Number of elements across axis 1": 6,
            "Number of elements across axis 2": 4,
            "Number of elements across axis 3": 8,
            "Number of transition elements": 1,
            "Size x": 1.5,
            "Size y": 1.0,
            "Size z": 2.0,
            "Axis 1 rotation y degrees": 0.0,
            "Axis 3 relative rotation y degrees": -90.0,
            "Refine": False,
            "Refine number of elements": 4,
        }
        return options

    @classmethod
    def getOrderedOptionNames(cls):
        return [
            "Number of elements across axis 1",
            "Number of elements across axis 2",
            "Number of elements across axis 3",
            "Number of transition elements",
            "Size x",
            "Size y",
            "Size z",
            "Axis 1 rotation y degrees",
            "Axis 3 relative rotation y degrees",
            "Refine",
            "Refine number of elements"
        ]

    @classmethod
    def checkOptions(cls, options):
        dependent_changes = False
        max_transition_count = None
        for key in [
            "Number of elements across axis 1",
            "Number of elements across axis 2",
            "Number of elements across axis 3"
        ]:
            if options[key] < 4:
                options[key] = 4
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
            "Size x",
            "Size y",
            "Size z"
        ]:
            if options[key] <= 0.0:
                options[key] = 1.0

        rotation_ranges = [
            ("Axis 1 rotation y degrees", -90, 90.0),
            ("Axis 3 relative rotation y degrees", -150.0, -30.0)
        ]
        for key, min_value, max_value in rotation_ranges:
            if options[key] <= min_value:
                options[key] = min_value
            elif options[key] > max_value:
                options[key] = max_value

        for key in [
            "Refine number of elements"
        ]:
            if options[key] < 1:
                options[key] = 1

        return dependent_changes

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: empty list of AnnotationGroup, None
        """
        element_counts = [options[key] for key in [
            "Number of elements across axis 1", "Number of elements across axis 2", "Number of elements across axis 3"]]
        transition_element_count = options["Number of transition elements"]
        sizes = [options[key] for key in ["Size x", "Size y", "Size z"]]
        axis1_rotation_radians = math.radians(options["Axis 1 rotation y degrees"])
        axis3_rotation_radians = axis1_rotation_radians + math.radians(options["Axis 3 relative rotation y degrees"])

        fieldmodule = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fieldmodule)

        ellipsoid = EllipsoidMesh(element_counts, transition_element_count, sizes,
                                  axis1_rotation_radians, axis3_rotation_radians)
        ellipsoid.build()
        ellipsoid.generateMesh(fieldmodule, coordinates)

        return [], None

    @classmethod
    def refineMesh(cls, meshRefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshRefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        assert isinstance(meshRefinement, MeshRefinement)
        refineElementsCount = options["Refine number of elements"]
        meshRefinement.refineAllElementsCubeStandard3d(refineElementsCount, refineElementsCount, refineElementsCount)
