"""
Utility functions for generating solid bifurcation and trifurcation.
"""
import math

from enum import Enum

from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import vector, geometry
from scaffoldmaker.utils.cylindermesh import Ellipse2D, EllipseShape, CylinderCentralPath, CylinderShape, CylinderEnds,\
    CylinderMesh
from scaffoldmaker.utils.derivativemoothing import DerivativeSmoothing
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite, \
    smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.mirror import Mirror
from scaffoldmaker.utils.shieldmesh import ShieldMesh3D
from scaffoldmaker.utils.spheremesh import SphereMesh, SphereShape
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues


class PathNodes:
    """
    Coordinates for the start and end nodes of the centralpath.
    """
    def __init__(self, part1, radius, length, elements_count, attach_bottom=True):
        """
        Extends part1 scaffold with a cylinder with given radius and length and elements count.
        :param part1: Scaffold with shield structure at its end.
        :param attach_bottom: If true, the path starts from the top of the par1. Otherwise top of the
         path attaches to the bottom of the part1.
        """
        if attach_bottom:
            n3 = -1
            sc = 1
        else:
            n3 = 0
            sc = -1
        # get nodes parameters from part1 end.
        csh = part1.px[n3][elements_count[1] // 2][elements_count[0] // 2]
        d1sh = part1.pd2[n3][elements_count[1] // 2][elements_count[0] // 2]
        d2sh = part1.pd1[n3][elements_count[1] // 2][elements_count[0] // 2]
        d3sh = part1.pd3[n3][elements_count[1] // 2][elements_count[0] // 2]
        d2sh = vector.setMagnitude(d2sh, -radius[0][0])
        d3sh = vector.setMagnitude(d3sh, radius[0][1])
        cw = vector.addVectors([csh, vector.setMagnitude(d1sh, sc * length)], [1, 1])
        d1w = vector.setMagnitude(d1sh, length / elements_count[2])
        d2w = vector.setMagnitude(d2sh, radius[1][0])
        d3w = vector.setMagnitude(d3sh, radius[1][1])

        if attach_bottom:
            path_list = [[csh, d1sh, d2sh, [0.0, 0.0, 0.0], d3sh, [0.0, 0.0, 0.0]],
                         [cw, d1w, d2w, [0.0, 0.0, 0.0], d3w, [0.0, 0.0, 0.0]]]
        else:
            path_list = [[cw, d1w, d2w, [0.0, 0.0, 0.0], d3w, [0.0, 0.0, 0.0]],
                         [csh, d1sh, d2sh, [0.0, 0.0, 0.0], d3sh, [0.0, 0.0, 0.0]]]
        self.path_list = path_list

    def get_path_list(self):
        return self.path_list


class BranchCylinder:
    """
    Generates a cylinder on top of the given part.
    """
    def __init__(self, region, mesh, nodes, fieldmodule, coordinates, path_list, elements_count, part1,
                 attach_bottom=True):
        """
        Generate a cylinder to extend part1. see PathNodes class.
        :param region: Zinc region
        :param mesh:  Zinc mesh.
        :param nodes: Zinc nodeset.
        :param fieldmodule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :param path_list: node parameters on the cylinder centralpath.
        :param elements_count:
        :param part1: Scaffold with shield structure at its end.
        :param attach_bottom: If true, the cylinder starts from the top of the par1. Otherwise top of the
         cylinder attaches to the bottom of the part1
        """
        if attach_bottom:
            n3, n3p = 0, -1
            node_ranges = [1, elements_count[2]]
        else:
            n3, n3p = -1, 0
            node_ranges = [0, elements_count[2] - 1]

        # generate the cylinder
        centralPath = ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': len(path_list) - 1
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                 Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3],
                path_list)
        })
        cylinderCentralPath = CylinderCentralPath(region, centralPath, elements_count[2])
        cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL
        centre = path_list[0][0]
        base = CylinderEnds(elements_count[1], elements_count[0], 0, 1, 1.0,
                            centre, cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
                            cylinderCentralPath.minorRadii[0])
        cylinder = CylinderMesh(fieldmodule, coordinates, elements_count[2], base,
                                cylinderShape=cylinderShape,
                                cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False,
                                rangeOfRequiredElementsAlong=[-1, -1])

        # skip generating common nodes with the part1.
        cylinder.generateNodes(nodes, fieldmodule, coordinates, node_ranges)

        # add the common node parameters.
        cylinder_shield = cylinder.getShield()
        for n2 in range(elements_count[1] + 1):
            for n1 in range(elements_count[0] + 1):
                cylinder_shield.nodeId[n3][n2][n1] = part1.nodeId[n3p][n2][n1]
                cylinder_shield.px[n3][n2][n1] = part1.px[n3p][n2][n1]
                cylinder_shield.pd1[n3][n2][n1] = part1.pd1[n3p][n2][n1]
                cylinder_shield.pd2[n3][n2][n1] = part1.pd2[n3p][n2][n1]
                cylinder_shield.pd3[n3][n2][n1] = part1.pd3[n3p][n2][n1]

        cylinder.generateElements(mesh, fieldmodule, coordinates)
        self.cylinder = cylinder

    def get_cylinder(self):
        return self.cylinder


class BranchCap:
    """
    Create a cap to attach to the given scaffold.
    """
    def __init__(self, fieldmodule, coordinates, mesh, nodes, part1, radius):
        """
        A hemisphere scaffold to attach to the part1 scaffold.
        :param part1: Scaffold with shield structure at its end.
        """
        sphere_shape = SphereShape.SPHERE_SHAPE_FULL
        sphere_base = part1._ellipses[-1]
        sphere_centre = sphere_base.centre
        sphere_radius_3 = radius
        axes = [sphere_base.majorAxis, sphere_base.minorAxis,
                vector.setMagnitude(vector.crossproduct3(sphere_base.majorAxis, sphere_base.minorAxis),
                                    sphere_radius_3)]
        elementsCountAcross = [part1._elementsCountAcrossMajor, part1._elementsCountAcrossMinor, 4]
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [0, -1]]
        sphereBoxDerivatives = [1, 3, 2]

        sphere1 = SphereMesh(fieldmodule, coordinates, sphere_centre, axes, elementsCountAcross,
                             0, 1, 1.0,
                             sphereShape=sphere_shape, rangeOfRequiredElements=rangeOfRequiredElements,
                             boxDerivatives=sphereBoxDerivatives, useCrossDerivatives=False)

        hemisphere = ShieldMesh3D(elementsCountAcross, 0)

        # get hemisphere nodes from both cylinder end and top of the sphere and mix them
        btx = hemisphere.px
        btd1 = hemisphere.pd1
        btd2 = hemisphere.pd2
        btd3 = hemisphere.pd3

        hemisphere._boxDerivatives = sphere1._shield3D._boxDerivatives
        hemisphere._boxMapping = sphere1._shield3D._boxMapping
        hemisphere._box_deriv_mapping = sphere1._shield3D._box_deriv_mapping
        hemisphere._element_needs_scale_factor = sphere1._shield3D._element_needs_scale_factor
        hemisphere._xi_mapping = sphere1._shield3D._xi_mapping
        hemisphere._xi_signs = sphere1._shield3D._xi_signs

        for n3 in range(elementsCountAcross[2] + 1):
            for n2 in range(elementsCountAcross[0] + 1):
                for n1 in range(elementsCountAcross[1] + 1):
                    if n3 > elementsCountAcross[2] // 2:
                        if sphere1._shield3D.px[n3][n2][n1]:
                            # hemisphere.nodeId[n3][n2][n1] = sphere1._shield3D.nodeId[n3][n2][n1]
                            btx[n3][n2][n1] = sphere1._shield3D.px[n3][n2][n1]
                            btd1[n3][n2][n1] = sphere1._shield3D.pd1[n3][n2][n1]
                            btd2[n3][n2][n1] = sphere1._shield3D.pd2[n3][n2][n1]
                            btd3[n3][n2][n1] = sphere1._shield3D.pd3[n3][n2][n1]

                    # cylinder end
                    elif n3 == elementsCountAcross[2] // 2:
                        # find nodes on the triple line. Note that cylinder and sphere have a little bit different
                        # numbering for nodes on the triple line
                        n2c, n1c = n2, n1
                        if n2 < 1 and n1 == n2:
                            n1c = 1
                        elif n2 < 1 and n1 == elementsCountAcross[1] - n2:
                            n1c = elementsCountAcross[1] - 1
                        elif n2 > elementsCountAcross[1] - 1:
                            if n1 == elementsCountAcross[1] - n2:
                                n1c = 1
                            elif n1 == n2:
                                n1c = elementsCountAcross[1] - 1
                        hemisphere.nodeId[n3][n2][n1] = part1._shield.nodeId[-1][n2c][n1c]

        # generate hemisphere extra nodes.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [3, 4]]
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        nodeIdentifier = hemisphere.generateNodes(fieldmodule, coordinates, nodeIdentifier,
                                                  rangeOfRequiredElements)

        # generate hemisphere elements.
        rangeOfRequiredElements = [[0, elementsCountAcross[0]], [0, elementsCountAcross[1]], [2, 4]]
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        elementIdentifier = hemisphere.generateElements(fieldmodule, coordinates, elementIdentifier,
                                                        rangeOfRequiredElements)

        self.sphere = hemisphere

    def get_sphere(self):
        return self.sphere


class BranchType(Enum):
    RIGHT_ARM = 1
    NECK = 2
    LEFT_ARM = 3


class TrifurcationMesh:
    """
    Trifurcation mesh generator.
    """
    def __init__(self, fieldmodule, coordinates, region, torso_radius, left_arm_radius, right_arm_radius, neck_radius,
                 shoulder_height, neck_height, right_arm_angle, left_arm_angle, right_shoulder_length, armpit,
                 elements_count):
        """
        :param fieldmodule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :param region: Zinc region
        :param torso_radius: upper torso radius
        :param left_arm_radius:
        :param right_arm_radius:
        :param neck_radius:
        :param shoulder_height:
        :param neck_height:
        :param right_arm_angle:
        :param left_arm_angle:
        :param right_shoulder_length:
        :param armpit:
        :param elements_count:
        """

        # generate the mesh
        elementsCount = [2, 2, 5]
        self._elementsCount = elementsCount
        self._elements_count = elements_count
        self._region = region

        self.torso_radius = torso_radius
        self.left_arm_radius = left_arm_radius
        self.right_arm_radius = right_arm_radius
        self.neck_radius = neck_radius
        self.shoulder_height = shoulder_height
        self.neck_height = neck_height
        self.right_arm_angle = right_arm_angle
        self.left_arm_angle = left_arm_angle
        self.right_shoulder_length = right_shoulder_length
        self.armpit = armpit

        self._coordinates = coordinates
        self._fieldmodule = fieldmodule
        self._nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._mesh = fieldmodule.findMeshByDimension(3)

        self._torso_upper_part = None
        self._shoulder_part = None
        self._shoulder_part_left = None
        self._neck_part = None

        self.createTrifurcationMesh3d(fieldmodule, coordinates)

    def createTrifurcationMesh3d(self, fieldmodule, coordinates):
        """
        Create a trifurcation.
        :param fieldmodule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :return: Final values of nextNodeIdentifier, nextElementIdentifier.
        """
        nodes = self._nodes
        mesh = self._mesh
        base_c, shoulder_lc, shoulder_rc, neck_c = self._get_node_params()

        elementsCount = self._elements_count
        # elementsCount = [6,6,2]
        torso_upper_part = BaseLeg(elementsCount, base_c)
        self.generateNodes(nodes, fieldmodule, coordinates, torso_upper_part)
        self.generateElements(mesh, fieldmodule, coordinates, torso_upper_part)
        self._torso_upper_part = torso_upper_part

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
        shoulder_part = BaseLeg(self._elementsCount, nodeparams2)
=======
        shoulder_part = BaseLeg(self._elementsCount, nodeparams2, shoulder_right=True)
>>>>>>> e91ccef (Add range of elements along cylinder parameter. Add arms,neck and head.)
=======
        shoulder_part = BaseLeg(self._elementsCount, shoulder_lc, shoulder_right=True)
>>>>>>> 9619665 (Modify axes and clean the code and make some classes)
=======
        shoulder_part = BaseLeg(elementsCount, shoulder_lc, shoulder_right=True)
>>>>>>> c6c23f0 (Change middle box and fix derivatives)
        shoulder_part._shoulder = True
        self.remove_duplicate_nodes_from_shoulder(shoulder_part)
        self.generateNodes(nodes, fieldmodule, coordinates, shoulder_part)
        self.join_shoulder_to_bottom_part(shoulder_part, torso_upper_part)
        self.generateElements(mesh, fieldmodule, coordinates, shoulder_part)
        self._shoulder_part = shoulder_part

        shoulder_part_left = BaseLeg(elementsCount, shoulder_rc, shoulder_left=True)
        shoulder_part_left._shoulder_left = True
        self.remove_duplicate_nodes_from_shoulder(shoulder_part_left, 1)
        self.generateNodes(nodes, fieldmodule, coordinates, shoulder_part_left)
        self.join_shoulder_to_bottom_part(shoulder_part_left, torso_upper_part, 1)
        self.generateElements(mesh, fieldmodule, coordinates, shoulder_part_left)
        self._shoulder_part_left = shoulder_part_left

        neck_part = BaseLeg(elementsCount, neck_c)
        neck_part._neck = True
        self.remove_duplicate_nodes_from_neck(neck_part)
        self.generateNodes(nodes, fieldmodule, coordinates, neck_part)
        self.join_neck_to_torso(neck_part, shoulder_part, shoulder_part_left)
        self.generateElements(mesh, fieldmodule, coordinates, neck_part)
        self._neck_part = neck_part

        box_part = BoxPart([elementsCount[0] - 2, elementsCount[1], elementsCount[0] - 2], torso_upper_part,
                           shoulder_part, shoulder_part_left, neck_part)
        self.generateNodes(nodes, fieldmodule, coordinates, box_part)
        self.generateElements(mesh, fieldmodule, coordinates, box_part)

    def create_branch_cylinder(self, radius, length, number_of_elements, path_list=None,
                               part1=None, attach_bottom=True, branch_type=BranchType.LEFT_ARM):
        """
        Creates a cylinder attached to the part1
        :param radius:
        :param length:
        :param number_of_elements:
        :param path_list:
        :param part1:
        :param attach_bottom:
        :param branch_type:
        :return:
        """
        if branch_type == BranchType.LEFT_ARM:
            part1 = self._shoulder_part
        elif branch_type == BranchType.RIGHT_ARM:
            part1 = self._shoulder_part_left
        elif branch_type == BranchType.NECK:
            part1 = self._neck_part

        if not path_list:
            pn = PathNodes(part1, radius, length, number_of_elements, attach_bottom=attach_bottom)
            path_list = pn.get_path_list()
        bc = BranchCylinder(self._region, self._mesh, self._nodes, self._fieldmodule, self._coordinates,
                            path_list, number_of_elements, part1, attach_bottom=attach_bottom)
        cylinder = bc.get_cylinder()

        return cylinder

    def create_branch_cap(self, part1, radius):
        """
        Creates a cap attached to the part1
        :param part1:
        :param radius:
        :return:
        """
        cap = BranchCap(self._fieldmodule, self._coordinates, self._mesh, self._nodes, part1, radius)

        return cap

    def smooth_all_derivatives(self):
        smoothing = DerivativeSmoothing(self._region, self._coordinates)
        smoothing.smooth(True)
        del smoothing

    def _get_node_params(self):
        """
        Get node parameter for the landmarks
        :return:
        """
        class CylinderCurves:
            def __init__(self, bottom_curves, top_curves):
                self.curves = [bottom_curves, top_curves]

        class EllipseCurves:
            def __init__(self, x_curve1, d1_curve1, x_curve2, d1_curve2):
                x_curve2[0] = [c for c in x_curve1[0]]
                self.xc1 = x_curve1
                self.d1c1 = d1_curve1
                self.xc2 = x_curve2
                self.d1c2 = d1_curve2

        armpit = self.armpit

        elementsCountQuarter = [3,3,2]
        bc = EllipseCurves([[0.0, 0.0, 0.0], [self.torso_radius, 0.0, 0.0]],
                           [[self.torso_radius / elementsCountQuarter[0], 0.0, 0.0]]*2,
                           [[0.0, 0.0, 0.0], [0.0, self.torso_radius, 0.0]],
                           [[0.0, 1/elementsCountQuarter[1], 0.0]]*2)

        tc = EllipseCurves([[0.0, 0.0, 1.4], armpit],
                           [[1 / elementsCountQuarter[0], 0.0, 0.0], [0.5 * 0.7071, 0.0, -0.5 * 0.7071]],
                           [[0.0, 0.0, 1.4], [0.0, 1.0, 1.4]],
                           [[0.0, 1 / elementsCountQuarter[1], 0.0]]*2)

        base_c = CylinderCurves(bc, tc)

        x_shoulder_base_centre = [0.75, 0.0, self.shoulder_height]
        kv = [0.0, 1.0, 0.0]
        cv = [self.right_shoulder_length, 0.0, 0.0]
<<<<<<< HEAD
        cev = vector.rotateVectorAroundVector(cv, kv, self.right_arm_angle)
<<<<<<< HEAD
        ce = vector.addVectors([cev, [0.0, 0.0, self.shoulder_height]], [1, 1])
=======
        cevunit = vector.normalise(cev)
        ce = vector.addVectors([cev, x_shoulder_base_centre], [1, 1])
<<<<<<< HEAD
>>>>>>> e91ccef (Add range of elements along cylinder parameter. Add arms,neck and head.)
        x_shoulder_end_centre = ce
        # x_shoulder_end_curve1 = [1.7, 0.0, self.shoulder_height - self.left_arm_radius]
=======
>>>>>>> 9619665 (Modify axes and clean the code and make some classes)
        x_shoulder_end_curve1 = vector.addVectors([ce, vector.setMagnitude(vector.crossproduct3(kv, cev), self.right_arm_radius)], [1, 1])
        x_shoulder_end_curve2 = vector.addVectors([ce, vector.setMagnitude(kv, self.right_arm_radius)], [1, 1])
<<<<<<< HEAD
        # d1_shoulder_end_curve1 = [[0.0, 0.0, -self.left_arm_radius / self._elementsCount[1]], [0.0, 0.0, -self.left_arm_radius / self._elementsCount[1]]]
        d1_shoulder_end_curve1 = [vector.setMagnitude(vector.crossproduct3(kv, cev), self.right_arm_radius/self._elementsCount[1]),
                                  vector.setMagnitude(vector.crossproduct3(kv, cev), self.right_arm_radius/self._elementsCount[1])]
<<<<<<< HEAD
=======

>>>>>>> e91ccef (Add range of elements along cylinder parameter. Add arms,neck and head.)
        d1_shoulder_end_curve2 = [[0.0, self.left_arm_radius / self._elementsCount[0], 0.0], [0.0, self.left_arm_radius / self._elementsCount[0], 0.0]]

        nodeparams2 = [[x_shoulder_base_centre, x_shoulder_base_curve1, x_shoulder_base_curve2, d1_shoulder_base_curve1,
                       d1_shoulder_base_curve2],
                      [x_shoulder_end_centre, x_shoulder_end_curve1, x_shoulder_end_curve2, d1_shoulder_end_curve1,
                       d1_shoulder_end_curve2]]
=======
        d1_shoulder_end_curve1 = [vector.setMagnitude(vector.crossproduct3(kv, cev), self.right_arm_radius/self._elementsCount[1])]*2
=======
        cev = vector.rotateVectorAroundVector(cv, kv, self.left_arm_angle)
        cevunit = vector.normalise(cev)
        ce = vector.addVectors([cev, x_shoulder_base_centre], [1, 1])
        x_shoulder_end_curve1 = vector.addVectors(
            [ce, vector.setMagnitude(vector.crossproduct3(kv, cev), self.left_arm_radius)], [1, 1])
        x_shoulder_end_curve2 = vector.addVectors([ce, vector.setMagnitude(kv, self.left_arm_radius)], [1, 1])
<<<<<<< HEAD
        d1_shoulder_end_curve1 = [vector.setMagnitude(vector.crossproduct3(kv, cev), self.left_arm_radius/self._elementsCount[1])]*2
>>>>>>> 5931d2a (Separate arms from upper torso)
=======
        d1_shoulder_end_curve1 = [vector.setMagnitude(vector.crossproduct3(kv, cev),
                                                      self.left_arm_radius/self._elementsCount[1])]*2
>>>>>>> 3e2b739 (Minor changes)
        d1_shoulder_end_curve2 = [[0.0, self.left_arm_radius / self._elementsCount[0], 0.0]]*2

        bc = EllipseCurves([[0.75, 0.0, self.shoulder_height], armpit],
                           [[0.0, 0.0, -1 / self._elementsCount[1]], [0.5 * 0.7071, 0.0, -0.5 * 0.7071]],
                           [[0.75, 0.0, self.shoulder_height], [0.75, 1.0, self.shoulder_height]],
                           [[0.0, 1 / self._elementsCount[0], 0.0]]*2)
        tc = EllipseCurves([ce, x_shoulder_end_curve1], d1_shoulder_end_curve1,
                           [ce, x_shoulder_end_curve2], d1_shoulder_end_curve2)
        shoulder_lc = CylinderCurves(bc, tc)
>>>>>>> 9619665 (Modify axes and clean the code and make some classes)

        x_shoulder_base_centre = [-0.75, 0.0, self.shoulder_height]
        kv = [0.0, 1.0, 0.0]
        cv = [-self.right_shoulder_length, 0.0, 0.0]
        cev = vector.rotateVectorAroundVector(cv, kv, -self.right_arm_angle)
        ce = vector.addVectors([cev, x_shoulder_base_centre], [1, 1])
        x_shoulder_end_curve2 = vector.addVectors(
            [ce, vector.setMagnitude(vector.crossproduct3(cev, kv), self.right_arm_radius)], [1, 1])
        x_shoulder_end_curve1 = vector.addVectors([ce, vector.setMagnitude(kv, self.right_arm_radius)], [1, 1])
        d1_shoulder_end_curve2 = [vector.setMagnitude(vector.crossproduct3(cev, kv),
                                                      self.right_arm_radius/self._elementsCount[1])]*2
        d1_shoulder_end_curve1 = [[-0.0, self.right_arm_radius / self._elementsCount[0], 0.0]]*2

        bc = EllipseCurves([[-0.75, 0.0, self.shoulder_height], [-0.75, 1.0, self.shoulder_height]],
                           [[-0.0, 1 / self._elementsCount[0], 0.0]]*2,
                           [[-0.75, 0.0, self.shoulder_height], [-armpit[0], armpit[1], armpit[2]]],
                           [[-0.0, 0.0, -1 / self._elementsCount[1]], [-0.5 * 0.7071, 0.0, -0.5 * 0.7071]])
        tc = EllipseCurves([ce, x_shoulder_end_curve1], d1_shoulder_end_curve1,
                           [ce, x_shoulder_end_curve2], d1_shoulder_end_curve2)
        shoulder_rc = CylinderCurves(bc, tc)

        x_neck_base_centre = [0.0, 0.0, 2 * self.shoulder_height - 1.0 - self.right_arm_radius/2]
        x_neck_base_curve2 = [0.0, 1.0, 2 * self.shoulder_height - 1.0 - self.right_arm_radius/2]
        x_neck_base_curve1 = [1.2, 0.0, 2 * self.shoulder_height - 1.0]
        d1_neck_base_curve2 = [[0.0, 1 / self._elementsCount[0], 0.0]]*2
        d1_neck_base_curve1 = [[1 / self._elementsCount[0], 0.0, 0.0], [0.5, 0.0, 0.5]]
        x_neck_end_centre = [0.0, 0.0, self.neck_height]
        x_neck_end_curve2 = [0.0, self.neck_radius, self.neck_height]
        x_neck_end_curve1 = [self.neck_radius, 0.0, self.neck_height]
        d1_neck_end_curve2 = [[0.0, self.neck_radius / self._elementsCount[1], 0.0]]*2
        d1_neck_end_curve1 = [[self.neck_radius / self._elementsCount[0], 0.0, 0.0]]*2

        bc = EllipseCurves([x_neck_base_centre, x_neck_base_curve1], d1_neck_base_curve1,
                           [x_neck_base_centre, x_neck_base_curve2], d1_neck_base_curve2)
        tc = EllipseCurves([x_neck_end_centre, x_neck_end_curve1], d1_neck_end_curve1,
                           [x_neck_end_centre, x_neck_end_curve2], d1_neck_end_curve2)
        neck_c = CylinderCurves(bc, tc)

        return base_c, shoulder_lc, shoulder_rc, neck_c

    def join_to_torso(self, joining_torso, torso, shoulder_joint, bottom_joint):
        """
        Attach torso to shoulder
        :param joining_torso:
        :param shoulder_joint:
        :param bottom_joint:
        :return:
        """
        for n2 in range(joining_torso._elementsCount[1] + 1):
            for n1 in range(joining_torso._elementsCount[1] + 1):
                joining_torso.nodeId[1][n2][n1] = torso.nodeId[0][n2][4 - n1]
                joining_torso.px[1][n2][n1] = torso.px[0][n2][4 - n1]
                joining_torso.pd1[1][n2][n1] = torso.pd1[0][n2][4 - n1]
                joining_torso.pd2[1][n2][n1] = torso.pd2[0][n2][4 - n1]
                joining_torso.pd3[1][n2][n1] = torso.pd3[0][n2][4 - n1]

                if n1 <= joining_torso._elementsCount[0]//2:
                    joining_torso.nodeId[0][n2][n1] = shoulder_joint.nodeId[0][n2][joining_torso._elementsCount[0] - n1]
                    joining_torso.px[0][n2][n1] = shoulder_joint.px[0][n2][joining_torso._elementsCount[0] - n1]
                    joining_torso.pd1[0][n2][n1] = shoulder_joint.pd1[0][n2][joining_torso._elementsCount[0] - n1]
                    joining_torso.pd2[0][n2][n1] = shoulder_joint.pd2[0][n2][joining_torso._elementsCount[0] - n1]
                    joining_torso.pd3[0][n2][n1] = shoulder_joint.pd3[0][n2][joining_torso._elementsCount[0] - n1]
                else:
                    joining_torso.nodeId[0][n2][n1] = bottom_joint.nodeId[1][n2][n1]
                    joining_torso.px[0][n2][n1] = bottom_joint.px[1][n2][n1]
                    joining_torso.pd1[0][n2][n1] = bottom_joint.pd1[1][n2][n1]
                    joining_torso.pd2[0][n2][n1] = bottom_joint.pd2[1][n2][n1]
                    joining_torso.pd3[0][n2][n1] = bottom_joint.pd3[1][n2][n1]

    def joint_shoulder_joint_to_cylinder_and_box(self, shoulder_connecting_to_box, joining_box, cylinder_part,
                                                 cidxs, bidx):
        """

        :param shoulder_connecting_to_box:
        :param joining_box:
        :param cylinder_part:
        :param cidxs:
        :param bidx:
        :return:
        """
        for n2 in range(shoulder_connecting_to_box._elementsCount[1] + 1):
            for n1 in range(shoulder_connecting_to_box._elementsCount[0]//2,
                            shoulder_connecting_to_box._elementsCount[0] + 1):
                shoulder_connecting_to_box.nodeId[cidxs[0]][n2][n1] = cylinder_part.nodeId[cidxs[1]][n2][n1]
                shoulder_connecting_to_box.px[cidxs[0]][n2][n1] = cylinder_part.px[cidxs[1]][n2][n1]
                shoulder_connecting_to_box.pd1[cidxs[0]][n2][n1] = cylinder_part.pd1[cidxs[1]][n2][n1]
                shoulder_connecting_to_box.pd2[cidxs[0]][n2][n1] = cylinder_part.pd2[cidxs[1]][n2][n1]
                shoulder_connecting_to_box.pd3[cidxs[0]][n2][n1] = cylinder_part.pd3[cidxs[1]][n2][n1]

            shoulder_connecting_to_box.nodeId[bidx][n2][2] = joining_box.nodeId[1][n2][1]
            shoulder_connecting_to_box.px[bidx][n2][2] = joining_box.px[1][n2][1]
            shoulder_connecting_to_box.pd1[bidx][n2][2] = joining_box.pd1[1][n2][1]
            shoulder_connecting_to_box.pd2[bidx][n2][2] = joining_box.pd2[1][n2][1]
            shoulder_connecting_to_box.pd3[bidx][n2][2] = joining_box.pd3[1][n2][1]

    def join_box_to_bottom_and_shoulder(self, joining_box, bottom_part, shoulder_part):
        """

        :param bottom_part:
        :param shoulder_part:
        :return:
        """

        for n2 in range(bottom_part._elementsCount[1] + 1):
            for n1 in range(2):
                joining_box.nodeId[0][n2][n1] = bottom_part.nodeId[bottom_part._elementsCount[2]][n2][n1 + 1]
                joining_box.px[0][n2][n1] = bottom_part.px[bottom_part._elementsCount[2]][n2][n1 + 1]
                joining_box.pd1[0][n2][n1] = bottom_part.pd1[bottom_part._elementsCount[2]][n2][n1 + 1]
                joining_box.pd2[0][n2][n1] = bottom_part.pd2[bottom_part._elementsCount[2]][n2][n1 + 1]
                joining_box.pd3[0][n2][n1] = bottom_part.pd3[bottom_part._elementsCount[2]][n2][n1 + 1]
            joining_box.nodeId[1][n2][0] = shoulder_part.nodeId[0][n2][2]
            joining_box.px[1][n2][0] = shoulder_part.px[0][n2][2]
            joining_box.pd1[1][n2][0] = shoulder_part.pd1[0][n2][2]
            joining_box.pd2[1][n2][0] = shoulder_part.pd2[0][n2][2]
            joining_box.pd3[1][n2][0] = shoulder_part.pd3[0][n2][2]

    def remove_duplicate_nodes_from_shoulder(self, shoulder_part, c=0):
        """

        :param shoulder_part:
        :param c:
        :return:
        """
        def condition(n2, n1):
            if c:
                return n2 == 0 or n2 == 1
            else:
                return n1 == 0 or n1 == 1

        for n3 in range(1):
            for n2 in range(shoulder_part._elementsCount[0] + 1):
                for n1 in range(shoulder_part._elementsCount[1] + 1):
                    if condition(n2, n1):
                        shoulder_part.px[n3][n2][n1] = None
                        shoulder_part.pd1[n3][n2][n1] = None
                        shoulder_part.pd2[n3][n2][n1] = None
                        shoulder_part.pd3[n3][n2][n1] = None

    def join_shoulder_to_bottom_part(self, shoulder_part, bottom_part, c=0):
        """

        :param shoulder_part:
        :param bottom_part:
        :return:
        """
        def condition(n2, n1):
            if c:
                return n2 == 0 or n2 == 1
            else:
                return n1 == 0 or n1 == 1

        def index(n2, n1):
            if c:
                if n2 == 0 and n1 == 1:
                    return n1 - 1, shoulder_part._elementsCount[0] - n2 - 1
                if n2 == 0 and n1 == shoulder_part._elementsCount[0] - 1:
                    return n1 + 1, shoulder_part._elementsCount[0] - n2 - 1
                else:
                    return n1, shoulder_part._elementsCount[0] - n2
            else:
                return n2, n1

        for n3 in range(1):
            for n2 in range(shoulder_part._elementsCount[0] + 1):
                for n1 in range(shoulder_part._elementsCount[1] + 1):
                    if condition(n2, n1):
                        n2b, n1b = index(n2, n1)
                        n3b = bottom_part._elementsCount[2]
                        shoulder_part.nodeId[n3][n2][n1] = bottom_part.nodeId[n3b][n2b][n1b]
                        shoulder_part.px[n3][n2][n1] = bottom_part.px[n3b][n2b][n1b]
                        shoulder_part.pd1[n3][n2][n1] = bottom_part.pd1[n3b][n2b][n1b]
                        shoulder_part.pd2[n3][n2][n1] = bottom_part.pd2[n3b][n2b][n1b]
                        shoulder_part.pd3[n3][n2][n1] = bottom_part.pd3[n3b][n2b][n1b]

    def copyBaseLeg2Trifurcation(self, baseleg, idx):
        """

        :param baseleg:
        :param idx:
        :return:
        """
        for n3 in range(self._elementsCount[2]//2 + 1):
            for n2 in range(self._elementsCount[0] + 1):
                for n1 in range(self._elementsCount[1] + 1):
                    if idx == 1:
                        n3s = n3
                    elif idx == 2:
                        n3s = self._elementsCount[2]//2 + 2 + n3
                    self.px[n3s][n2][n1] = baseleg.px[n3][n2][n1]
                    self.pd1[n3s][n2][n1] = baseleg.pd1[n3][n2][n1]
                    self.pd2[n3s][n2][n1] = baseleg.pd2[n3][n2][n1]
                    self.pd3[n3s][n2][n1] = baseleg.pd3[n3][n2][n1]
                    if idx == 2 and n3 == 0:
                        if (n2 == 0 and n1 == 1) or (n2 == 1 and n1 == 1) or (n2 == 2 and n1 == 0) or (n2 == 2 and
                                                                                                       n1 == 1):
                            self.px[n3s][n2][n1] = None
                            self.pd1[n3s][n2][n1] = None
                            self.pd2[n3s][n2][n1] = None
                            self.pd3[n3s][n2][n1] = None

    def remove_duplicate_nodes_from_neck(self, neck_part):
        """

        :param neck_part:
        :return:
        """
        for n2 in range(neck_part._elementsCount[1] + 1):
            for n1 in range(neck_part._elementsCount[0] + 1):
                if n1 <= 1 or n1 >= neck_part._elementsCount[0] - 1:
                    neck_part.px[0][n2][n1] = None
                    neck_part.pd1[0][n2][n1] = None
                    neck_part.pd2[0][n2][n1] = None
                    neck_part.pd3[0][n2][n1] = None

    def join_neck_to_torso(self, neck_part, shoulder_part, shoulder_part_left):
        """

        :param neck_part:
        :param shoulder_part:
        :param shoulder_part_left:
        :return:
        """
        for n2 in range(neck_part._elementsCount[1] + 1):
            for n1 in range(neck_part._elementsCount[0] + 1):
                if n1 < 2:
                    n1s = neck_part._elementsCount[0] - n1
                    neck_part.nodeId[0][n2][n1] = shoulder_part.nodeId[0][n2][n1s]
                    neck_part.px[0][n2][n1] = shoulder_part.px[0][n2][n1s]
                    neck_part.pd1[0][n2][n1] = shoulder_part.pd1[0][n2][n1s]
                    neck_part.pd2[0][n2][n1] = shoulder_part.pd2[0][n2][n1s]
                    neck_part.pd3[0][n2][n1] = shoulder_part.pd3[0][n2][n1s]
                elif n1 >= neck_part._elementsCount[0] - 1:
                    if n2 == 0 or n2 == neck_part._elementsCount[1]:
                        if n1 == neck_part._elementsCount[0] - 1:
                            n2s = n1+1
                            if n2 == 0:
                                n1s = n2 + 1
                            if n2 == neck_part._elementsCount[1]:
                                n1s = n2 - 1
                    else:
                        n2s = n1
                        n1s = n2
                    if n1 == 4:
                        if n2 !=2:
                            continue
                    neck_part.nodeId[0][n2][n1] = shoulder_part_left.nodeId[0][n2s][n1s]
                    neck_part.px[0][n2][n1] = shoulder_part_left.px[0][n2s][n1s]
                    neck_part.pd1[0][n2][n1] = shoulder_part_left.pd1[0][n2s][n1s]
                    neck_part.pd2[0][n2][n1] = shoulder_part_left.pd2[0][n2s][n1s]
                    neck_part.pd3[0][n2][n1] = shoulder_part_left.pd3[0][n2s][n1s]

    def generateNodes(self, nodes, fieldModule, coordinates, part_structure):
        """
        Create cylinder nodes from coordinates.
        :param nodes: nodes from coordinates.
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self.topologygenerateNodes(fieldModule, coordinates, nodeIdentifier, part_structure)
        self._endNodeIdentifier = nodeIdentifier

    def generateElements(self, mesh, fieldModule, coordinates, part_structure):
        """
        Create cylinder elements from nodes.
        :param mesh:
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        self._startElementIdentifier = elementIdentifier
        elementIdentifier = self.topologygenerateElements(fieldModule, coordinates, elementIdentifier, part_structure,
                                                          [])
        self._endElementIdentifier = elementIdentifier

    def topologygenerateNodes(self, fieldmodule, coordinates, startNodeIdentifier, part_structure):
        """
        Create shield nodes from coordinates.
         """
        nodeIdentifier = startNodeIdentifier
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        cache = fieldmodule.createFieldcache()

        for n3 in range(part_structure._elementsCount[2] + 1):
            for n2 in range(part_structure._elementsCount[1] + 1):
                for n1 in range(part_structure._elementsCount[0] + 1):
                    if part_structure.px[n3][n2][n1]:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        part_structure.nodeId[n3][n2][n1] = nodeIdentifier
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1,
                                                      part_structure.px[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1,
                                                      part_structure.pd1[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1,
                                                      part_structure.pd2[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                                      part_structure.pd3[n3][n2][n1])
                        nodeIdentifier += 1

        return nodeIdentifier

    def topologygenerateElements(self, fieldmodule, coordinates, startElementIdentifier, part_structure, meshGroups=[]):
        """
        Create shield elements from nodes.
        :param fieldmodule: Zinc fieldmodule to create elements in.
        :param coordinates: Coordinate field to define.
        :param startElementIdentifier: First element identifier to use.
        :param meshGroups: Zinc mesh groups to add elements to.
        :return: next elementIdentifier.
         """
        elementIdentifier = startElementIdentifier
        useCrossDerivatives = False
        mesh = fieldmodule.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftNoCrossDerivatives()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # isEven = (self.elementsCountAcross % 2) == 0
        e1a = 0
        e1b = e1a + 1
        e1z = part_structure._elementsCount[0] - 1
        e1y = e1z - 1
        e2a = 0
        e2b = e2a + 1
        e2c = e2a + 2
        e2z = part_structure._elementsCount[1]-1
        e2y = e2z - 1
        # e2x = e2z - 2
        for e3 in range(part_structure._elementsCount[2]):
            for e2 in range(part_structure._elementsCount[1]):
                for e1 in range(part_structure._elementsCount[0]):
                    eft1 = eft
                    scalefactors = None
                    nids = [part_structure.nodeId[e3][e2][e1], part_structure.nodeId[e3][e2 + 1][e1],
                            part_structure.nodeId[e3+1][e2][e1], part_structure.nodeId[e3+1][e2 + 1][e1],
                            part_structure.nodeId[e3][e2][e1 + 1], part_structure.nodeId[e3][e2 + 1][e1 + 1],
                            part_structure.nodeId[e3+1][e2][e1 + 1], part_structure.nodeId[e3+1][e2 + 1][e1 + 1]]

                    if isinstance(part_structure, BoxPart):
                        if e3 == 0:
                            if e1 == 0:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e2 == 0:
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                elif e2 == part_structure._elementsCount[1]:
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])

                                if e2 == e2a or e2 == e2z:
                                    if e2 == e2a:
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                    elif e2 == e2z:
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                else:
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                            elif e1 == e1z:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if e2 == e2a:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                elif e2 == e2z:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])

                                else:
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                            else:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if e2 == e2a:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                if e2 == e2z:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                        elif e3 == part_structure._elementsCount[2] - 1:
                            if e1 == 0:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e2 == e2a:
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])

                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                elif e2 == e2z:
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                else:
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])

                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                            elif e1 == e1z:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if e2 == e2a:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                elif e2 == e2z:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, []), (Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                else:
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                            else:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if e2 == e2a:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                elif e2 == e2z:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])

                        else:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            if e1 == 0:
                                if e2 == e2a:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                elif e2 == e2z:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                else:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                            elif e1 == e1z:
                                if e2 == e2a:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                elif e2 == e2z:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                else:
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D2_DS1DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                            else:
                                if e2 == e2a:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                elif e2 == e2z:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])

                    elif part_structure._neck and e3 == 0:
                        if (e2 < e2b) or (e2 > e2y):
                            if (e1 < e1b) or (e1 > e1y):
                                continue  # no element due to triple point closure
                        if (e2 == e2a) or (e2 == e2z):
                            # bottom and top row elements
                            if e2 == e2a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e1 == e1a:
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                elif e1 == e1y:

                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, []),
                                                            (Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1]),
                                                            (Node.VALUE_LABEL_D_DS1, [])])

                                elif e1b < e1 < e1y:
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,

                                                           [(Node.VALUE_LABEL_D_DS1, [])])

                                else:
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    # if (e1 == e1b) or (e1 == e1y):
                                    if e1 == e1b:
                                        # map bottom triple point element
                                        if e1 == e1b:
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [1])])
                                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [])])
                                            if e3 == 0:
                                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                                       [(Node.VALUE_LABEL_D_DS3, [])])
                                            else:
                                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                                       [(Node.VALUE_LABEL_D_DS2, [1]),
                                                                        (Node.VALUE_LABEL_D_DS3, [])])
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [1])])
                            elif e2 == e2z:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map top triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                    elif e1 == e1y:
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                elif e1b < e1 < e1y:
                                    remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])

                        elif (e2 == e2b) or (e2 == e2y):
                            if (e1 <= e1a) or (e1 >= e1z):
                                if e1 == e1a:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    if e2 == e2b:
                                        if e3 == part_structure._elementsCount[2] // 2 + 1:
                                            e3r = e3-1  # to join upper leg with the lower leg.
                                            nids[0] = part_structure.nodeId[e3r][e2a][e1b]
                                            nids[2] = part_structure.nodeId[e3+1][e2a][e1b]
                                            nids[1] = part_structure.nodeId[e3r][e2 + 1][e1]
                                            nids[4] = part_structure.nodeId[e3r][e2][e1 + 1]
                                            nids[5] = part_structure.nodeId[e3r][e2 + 1][e1 + 1]
                                        else:
                                            nids[0] = part_structure.nodeId[e3][e2a][e1b]
                                            nids[2] = part_structure.nodeId[e3+1][e2a][e1b]
                                        if e3 == 0:
                                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, [1]),
                                                                    (Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, [1])])

                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                    elif e2 == e2y:
                                        nids[1] = part_structure.nodeId[e3][e2z+1][e1b]
                                        nids[3] = part_structure.nodeId[e3+1][e2z+1][e1b]
                                        if e3 == 0:
                                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, [1]),
                                                                    (Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                elif e1 == e1z:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    if e2 == e2b:
                                        nids[4] = part_structure.nodeId[e3][e2a][e1z]
                                        nids[6] = part_structure.nodeId[e3+1][e2a][e1z]
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [-1.0]
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                        if e3 == 0:
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D2_DS1DS2,
                                                                   [(Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS2, [])])
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1]),
                                                                (Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                    elif e2 == e2y:
                                        nids[5] = part_structure.nodeId[e3][e2z+1][e1z]
                                        nids[7] = part_structure.nodeId[e3+1][e2z+1][e1z]
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [-1.0]
                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS2, [1])])
                                        if e3 == 0:
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS2, [1])])
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D2_DS1DS2,
                                                                   [(Node.VALUE_LABEL_D_DS3, [])])

                                        else:
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS2, [1])])
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS2, [])])
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS3, [])])

                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                            elif e1 == e1b:
                                if e2 == e2b:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    if e3 == 0:
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS2, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                if e2 == e2y:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [1]),
                                                            (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                            elif e1 == e1y:
                                if e2 == e2b:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]

                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                if e2 == e2y:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                            (Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])

                        elif e2b < e2 < e2y:
                            if e1 == e1a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, [1]),
                                                        (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [])])
                                remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                            elif e1 == e1b:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                            elif e1 == e1y:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [])])
                            elif e1 == e1z:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D2_DS1DS2,
                                                       [(Node.VALUE_LABEL_D_DS1, [])])
                                remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])

                        else:
                            if e1 < e1a:
                                nids = [part_structure.nodeId[e3][e2 + 1][e1 + 1],
                                        part_structure.nodeId[e3][e2][e1 + 1],
                                        part_structure.nodeId[e3+1][e2 + 1][e1 + 1],
                                        part_structure.nodeId[e3+1][e2][e1 + 1],
                                        part_structure.nodeId[e3][e2 + 1][e1],
                                        part_structure.nodeId[e3][e2][e1],
                                        part_structure.nodeId[e3+1][e2 + 1][e1],
                                        part_structure.nodeId[e3+1][e2][e1]]
                            elif e1 == e1a:
                                # map left column elements
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])

                    else:
                        if (e2 < e2b) or (e2 > e2y):
                            if (e1 < e1b) or (e1 > e1y):
                                continue  # no element due to triple point closure
                            if (e2 < e2a) or (e2 > e2z):
                                if e2 < e2a:
                                    nids = [part_structure.nodeId[e3][e2+1][e1],
                                            part_structure.nodeId[e3][e2+1][e1+1],
                                            part_structure.nodeId[e3+1][e2+1][e1],
                                            part_structure.nodeId[e3+1][e2+1][e1+1],
                                            part_structure.nodeId[e3][e2][e1],
                                            part_structure.nodeId[e3][e2][e1+1],
                                            part_structure.nodeId[e3+1][e2][e1],
                                            part_structure.nodeId[e3+1][e2][e1+1]]
                                elif e2 > e2z:
                                    nids = [part_structure.nodeId[e3][e2][e1+1],
                                            part_structure.nodeId[e3][e2][e1],
                                            part_structure.nodeId[e3+1][e2][e1+1],
                                            part_structure.nodeId[e3+1][e2][e1],
                                            part_structure.nodeId[e3][e2+1][e1+1],
                                            part_structure.nodeId[e3][e2+1][e1],
                                            part_structure.nodeId[e3+1][e2+1][e1+1],
                                            part_structure.nodeId[e3+1][e2+1][e1]]
                            elif (e2 == e2a) or (e2 == e2z):
                                # bottom and top row elements
                                if e2 == e2a:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    if not part_structure._shoulder_left:
                                        if isinstance(part_structure, BaseLeg):
                                            if part_structure._shoulder:
                                                if e3 == 0 and e1 == e1b:
                                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                                           [(Node.VALUE_LABEL_D_DS2, []),
                                                                            (Node.VALUE_LABEL_D_DS1, [])])
                                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                                           [(Node.VALUE_LABEL_D_DS2, []),
                                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                        if e3 == 0:
                                            if e2 == 0 and (e1 <= e1y-1 and e1 >= e1b+1):
                                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS3,
                                                                       [(Node.VALUE_LABEL_D_DS1, [])])
                                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
                                                                       [(Node.VALUE_LABEL_D_DS2, [])])
                                    if (e1 == e1b) or (e1 == e1y):
                                        # map bottom triple point element
                                        if e1 == e1b:
                                            # if e3 != 2:
                                            if part_structure._shoulder_left:
                                                if e3 == 0:
                                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                                else:
                                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                                            (Node.VALUE_LABEL_D_DS3, [])])
                                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                                       [(Node.VALUE_LABEL_D_DS1, []),
                                                                        (Node.VALUE_LABEL_D_DS3, [])])
                                            else:
                                                if part_structure._shoulder:
                                                    if e3 == 0:
                                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                                    else:
                                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1,
                                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                                else:
                                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1,
                                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                                            (Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            if part_structure._shoulder_left:
                                                if e3 == 0:
                                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,
                                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                                else:
                                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1,
                                                                           [(Node.VALUE_LABEL_D_DS1, []),
                                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                                remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1,
                                                                       [(Node.VALUE_LABEL_D_DS1, []),
                                                                        (Node.VALUE_LABEL_D_DS3, [1])])
                                            else:
                                                remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1,
                                                                       [(Node.VALUE_LABEL_D_DS1, []),
                                                                        (Node.VALUE_LABEL_D_DS3, [1])])
                                elif e2 == e2z:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    if isinstance(part_structure, BaseLeg):
                                        if part_structure._shoulder and e3 == 0 and e1 == e1b:
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                    (Node.VALUE_LABEL_D_DS2, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                    if (e1 == e1b) or (e1 == e1y):
                                        # map top triple point element
                                        if e1 == e1b:
                                            remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [1])])
                                            if e3 == 0:
                                                if part_structure._shoulder:
                                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                                           [(Node.VALUE_LABEL_D_DS2, []),
                                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                            if part_structure._shoulder:
                                                if e3 == 0:
                                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [])])

                        elif (e2 == e2b) or (e2 == e2y):
                            if (e1 <= e1a) or (e1 >= e1z):
                                if e1 < e1a:
                                    e2r = e1
                                    if e2 == e2b:
                                        nids = [part_structure.nodeId[e3][e2c][e1+1],
                                                part_structure.nodeId[e3][e2r+1][e1b],
                                                part_structure.nodeId[e3+1][e2c][e1+1],
                                                part_structure.nodeId[e3+1][e2r+1][e1b],
                                                part_structure.nodeId[e3][e2c][e1],
                                                part_structure.nodeId[e3][e2r][e1b],
                                                part_structure.nodeId[e3+1][e2c][e1],
                                                part_structure.nodeId[e3+1][e2r][e1b]]
                                    if e2 == e2y:
                                        e2r = 2*part_structure._elementsCount[1] - e1-1
                                        nids = [part_structure.nodeId[e3][e2r][e1b],
                                                part_structure.nodeId[e3][e2y][e1+1],
                                                part_structure.nodeId[e3+1][e2r][e1b],
                                                part_structure.nodeId[e3+1][e2y][e1+1],
                                                part_structure.nodeId[e3][e2r+1][e1b],
                                                part_structure.nodeId[e3][e2y][e1],
                                                part_structure.nodeId[e3+1][e2r+1][e1b],
                                                part_structure.nodeId[e3+1][e2y][e1]]
                                elif e1 == e1a:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    if e2 == e2b:
                                        if e3 == part_structure._elementsCount[2] // 2 + 1:
                                            e3r = e3-1  # to join upper leg with the lower leg.
                                            nids[0] = part_structure.nodeId[e3r][e2a][e1b]
                                            nids[2] = part_structure.nodeId[e3+1][e2a][e1b]
                                            nids[1] = part_structure.nodeId[e3r][e2 + 1][e1]
                                            nids[4] = part_structure.nodeId[e3r][e2][e1 + 1]
                                            nids[5] = part_structure.nodeId[e3r][e2 + 1][e1 + 1]
                                        else:
                                            nids[0] = part_structure.nodeId[e3][e2a][e1b]
                                            nids[2] = part_structure.nodeId[e3+1][e2a][e1b]

                                        if part_structure._shoulder_left:
                                            if e3 == 0:
                                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                                       [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                        (Node.VALUE_LABEL_D_DS2, [])])
                                                remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                                       [(Node.VALUE_LABEL_D_DS2, [])])
                                                remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                                       [(Node.VALUE_LABEL_D_DS1, []),
                                                                        (Node.VALUE_LABEL_D_DS3, [1])])
                                            else:
                                                remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,
                                                                       [(Node.VALUE_LABEL_D_DS1, []),
                                                                        (Node.VALUE_LABEL_D_DS3, [])])
                                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            tripleN = [5, 7]
                                            remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [])])
                                    elif e2 == e2y:
                                        nids[1] = part_structure.nodeId[e3][e2z+1][e1b]
                                        nids[3] = part_structure.nodeId[e3+1][e2z+1][e1b]
                                        tripleN = [6, 8]
                                        remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    if e3 == 0:
                                        if part_structure._shoulder:
                                            remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [1])])
                                elif e1 == e1z:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    if e2 == e2b:
                                        nids[4] = part_structure.nodeId[e3][e2a][e1z]
                                        nids[6] = part_structure.nodeId[e3+1][e2a][e1z]
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [-1.0]
                                        if part_structure._shoulder_left:
                                            if e3 == 0:
                                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                                       [(Node.VALUE_LABEL_D_DS1, []),
                                                                        (Node.VALUE_LABEL_D_DS3, [])])
                                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1,
                                                                       [(Node.VALUE_LABEL_D_DS2, [])])
                                                remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1,
                                                                       [(Node.VALUE_LABEL_D_DS1, []),
                                                                        (Node.VALUE_LABEL_D_DS2, [])])
                                            else:
                                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                                       [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                        (Node.VALUE_LABEL_D_DS3, [])])
                                            remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                    (Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                    (Node.VALUE_LABEL_D_DS3, [])])
                                    elif e2 == e2y:
                                        nids[5] = part_structure.nodeId[e3][e2z+1][e1z]
                                        nids[7] = part_structure.nodeId[e3+1][e2z+1][e1z]
                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                elif e1 > e1z:
                                    e2r = part_structure._elementsCount[0] - e1
                                    if e2 == e2b:
                                        nids = [part_structure.nodeId[e3][e2r][e1z],
                                                part_structure.nodeId[e3][e2c][e1],
                                                part_structure.nodeId[e3+1][e2r][e1z],
                                                part_structure.nodeId[e3+1][e2c][e1],
                                                part_structure.nodeId[e3][e2r-1][e1z],
                                                part_structure.nodeId[e3][e2c][e1+1],
                                                part_structure.nodeId[e3+1][e2r-1][e1z],
                                                part_structure.nodeId[e3+1][e2c][e1+1]]
                                    elif e2 == e2y:
                                        e2r = e2z+e1-e1z
                                        nids[1] = part_structure.nodeId[e3][e2r][e1z]
                                        nids[3] = part_structure.nodeId[e3+1][e2r][e1z]
                                        nids[5] = part_structure.nodeId[e3][e2r+1][e1z]
                                        nids[7] = part_structure.nodeId[e3+1][e2r+1][e1z]
                            elif e1 == e1b:
                                if e2 == e2b:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    if part_structure._shoulder_left:
                                        if e3 == 0:
                                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                                    if part_structure._shoulder:
                                        if e3 == 0:
                                            setEftScaleFactorIds(eft1, [1], [])
                                            scalefactors = [-1.0]
                                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [1])])
                                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS2, [])])

                                if e2 == e2y:
                                    if part_structure._shoulder:
                                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                                        if e3 == 0:
                                            setEftScaleFactorIds(eft1, [1], [])
                                            scalefactors = [-1.0]
                                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                                   [(Node.VALUE_LABEL_D_DS2, []),
                                                                    (Node.VALUE_LABEL_D_DS3, [1])])
                                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS2, [])])

                            elif e1 == e1y:
                                if e2 == e2b:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    if part_structure._shoulder_left:
                                        if e3 == 0:
                                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS3,
                                                                   [(Node.VALUE_LABEL_D_DS1, [])])

                            else:
                                if part_structure._shoulder_left:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    if e3 == 0 and e2 <= e2b:
                                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])

                        else:
                            if e1 < e1a:
                                nids = [part_structure.nodeId[e3][e2 + 1][e1 + 1],
                                        part_structure.nodeId[e3][e2][e1 + 1],
                                        part_structure.nodeId[e3+1][e2 + 1][e1 + 1],
                                        part_structure.nodeId[e3+1][e2][e1 + 1],
                                        part_structure.nodeId[e3][e2 + 1][e1],
                                        part_structure.nodeId[e3][e2][e1],
                                        part_structure.nodeId[e3+1][e2 + 1][e1],
                                        part_structure.nodeId[e3+1][e2][e1]]
                            elif e1 == e1a:
                                # map left column elements
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                if e3 == 0 and part_structure._shoulder:
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, []),
                                                            (Node.VALUE_LABEL_D_DS3, [1])])

                            elif e1 == e1b:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if e3 == 0 and part_structure._shoulder:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, []),
                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, []),
                                                            (Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])

                    if not all(nids):
                        continue

                    if eft1 is not eft:
                        elementtemplate1.defineField(coordinates, -1, eft1)
                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                    else:
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7
                    part_structure.elementId[e3][e2][e1] = elementIdentifier
                    elementIdentifier += 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

        return elementIdentifier


class BaseLeg:
    """
    Base case for creating a child
    """
<<<<<<< HEAD
<<<<<<< HEAD
    def __init__(self, elementsCount, nodeparams):
=======
    def __init__(self, elementsCount, nodeparams, shoulder_right=False, shoulder_left=False):
>>>>>>> e91ccef (Add range of elements along cylinder parameter. Add arms,neck and head.)
=======
    def __init__(self, elementsCount, cyl_curves, shoulder_right=False, shoulder_left=False):
>>>>>>> 9619665 (Modify axes and clean the code and make some classes)
        """

        :param elementsCount:
        :param cyl_curves:
        :param shoulder_right:
        :param shoulder_left:
        """
        self._elementsCount = [elementsCount[0]//2, elementsCount[1]//2, elementsCount[2]]
        self._shoulder = False
        self._shoulder_left = False
        self._neck = False

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                   for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                       for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])]
                          for c in range(elementsCount[2])]
        self.generateBaseLeg(cyl_curves.curves)

        bc = cyl_curves.curves[0]
        tc = cyl_curves.curves[1]
        n, d = geometry.get_plane_normal_vector_and_distance(bc.xc1[0], bc.xc1[1], tc.xc1[0])
        plane = [n[0], n[1], n[2], d]
        mirrore = Mirror(plane)
        if shoulder_left:
            mirror0 = Mirror([0.0, 0.0, 1.0, bc.xc1[0][2]])
        for n2 in range(elementsCount[1]//2+1, elementsCount[1]+1):
            for n3 in range(elementsCount[2]+1):
                if shoulder_left and n3 == 0:
                    mirror = mirror0
                else:
                    mirror = mirrore
                for n1 in range(elementsCount[0]//2+1):
                    n3q = n3
                    n2q = elementsCount[1] - n2
                    n1q = n1
                    if self.px[n3q][n2q][n1q]:
                        self.px[n3][n2][n1] = mirror.mirrorImageOfPoint(self.px[n3q][n2q][n1q])
                        self.pd1[n3][n2][n1] = mirror.reverseMirrorVector(self.pd1[n3q][n2q][n1q])
                        self.pd2[n3][n2][n1] = mirror.mirrorVector(self.pd2[n3q][n2q][n1q])
                        self.pd3[n3][n2][n1] = mirror.mirrorVector(self.pd3[n3q][n2q][n1q])

        n, d = geometry.get_plane_normal_vector_and_distance(bc.xc1[0], bc.xc2[1], tc.xc1[0])
        plane = [n[0], n[1], n[2], d]
<<<<<<< HEAD
        mirror = Mirror(plane)
        for n2 in range(elementsCount[1]+1):
            for n3 in range(elementsCount[2]+1):
=======
        mirrore = Mirror(plane)
        if shoulder_right:
            mirror0 = Mirror([0.0, 0.0, 1.0, bc.xc1[0][2]])
        for n2 in range(elementsCount[1]+1):
            for n3 in range(elementsCount[2]+1):
                if shoulder_right and n3 == 0:
                    mirror = mirror0
                else:
                    mirror = mirrore
>>>>>>> e91ccef (Add range of elements along cylinder parameter. Add arms,neck and head.)
                for n1 in range(elementsCount[0]//2+1, elementsCount[0]+1):
                    n3q = n3
                    n2q = n2
                    n1q = elementsCount[0] - n1
                    if self.px[n3q][n2q][n1q]:
                        self.px[n3][n2][n1] = mirror.mirrorImageOfPoint(self.px[n3q][n2q][n1q])
                        if n2 == 0 or n2 == elementsCount[1] or n1 == elementsCount[0]:
                            self.pd1[n3][n2][n1] = mirror.reverseMirrorVector(self.pd1[n3q][n2q][n1q])
                            self.pd3[n3][n2][n1] = mirror.mirrorVector(self.pd3[n3q][n2q][n1q])
                        else:
                            self.pd1[n3][n2][n1] = mirror.mirrorVector(self.pd1[n3q][n2q][n1q])
                            self.pd3[n3][n2][n1] = mirror.reverseMirrorVector(self.pd3[n3q][n2q][n1q])
                        self.pd2[n3][n2][n1] = mirror.mirrorVector(self.pd2[n3q][n2q][n1q])
<<<<<<< HEAD
=======

        # sample in between ellipses
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3
        elementsCountOut = 2
        for n2 in range(elementsCount[1] + 1):
            for n1 in range(elementsCount[0] + 1):
                if btx[0][n2][n1]:
                    tx, td2, pe, pxi, psf = sampleCubicHermiteCurves(
                        [btx[0][n2][n1], btx[self._elementsCount[2]][n2][n1]],
                        [btd2[0][n2][n1], btd2[self._elementsCount[2]][n2][n1]],
                        elementsCountOut)
                    td1 = interpolateSampleCubicHermite(
                        [btd1[0][n2][n1], btd1[self._elementsCount[2]][n2][n1]],
                        [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
                    td3 = interpolateSampleCubicHermite(
                        [btd3[0][n2][n1], btd3[self._elementsCount[2]][n2][n1]],
                        [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]

                    for n3 in range(1, self._elementsCount[2]):
                        n3i = n3
                        self.px[n3][n2][n1] = tx[n3i]
                        self.pd1[n3][n2][n1] = td1[n3i]
                        self.pd2[n3][n2][n1] = td2[n3i]
                        self.pd3[n3][n2][n1] = td3[n3i]

>>>>>>> 9619665 (Modify axes and clean the code and make some classes)
        self._elementsCount = elementsCount

    def generateBaseLeg(self, surface_curves):
        """
        Generate base leg that is a cylinder generated from cylinder ends.
        :return:
        """
        bottom_curves = surface_curves[0]
        top_curves = surface_curves[1]
        self.generate_surface(bottom_curves, 0)
        self.generate_surface(top_curves, self._elementsCount[2])
        self.generateMiddleLevels()
        self.smoothd2()

    def generate_surface(self, curves, n3):
        """

        :return:
        """
        centre = curves.xc1[0]
        txc1, td1c1 = self.generate1DPath(curves.xc1, curves.d1c1, self._elementsCount[0])
        txc2, td1c2 = self.generate1DPath(curves.xc2, curves.d1c2, self._elementsCount[1])
        ellipse = self.generateSurfaceUsingTwoCurves(centre, txc1, td1c1, txc2, td1c2)
        self.copyEllipseNodesToTrifurcation(ellipse, n3)

    def generateMiddleLevels(self):
        """

        :return:
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3
<<<<<<< HEAD
        # generate the armpit curves.
        elementsCountOut = 2
        txcc1, td1cc1, pec, pxic, psfc = sampleCubicHermiteCurves([btx[0][self._elementsCount[1]][self._elementsCount[0]],
                                                                btx[self._elementsCount[2]][self._elementsCount[1]][self._elementsCount[0]]],
                                                               [btd2[0][self._elementsCount[1]][self._elementsCount[0]],
                                                                btd2[self._elementsCount[2]][self._elementsCount[1]][self._elementsCount[0]]],
                                                               elementsCountOut)
        txec1, td1ec1, pec1, pxic1, psfc1 = sampleCubicHermiteCurves([btx[0][self._elementsCount[1]][0], btx[self._elementsCount[2]][self._elementsCount[1]][0]],
                                                               [btd2[0][self._elementsCount[1]][0], btd2[self._elementsCount[2]][self._elementsCount[1]][0]],
                                                               elementsCountOut)
        txec2, td1ec2, pec2, pxic2, psfc2 = sampleCubicHermiteCurves([btx[0][0][self._elementsCount[0]], btx[self._elementsCount[2]][0][self._elementsCount[0]]],
                                                               [btd2[0][0][self._elementsCount[0]], btd2[self._elementsCount[2]][0][self._elementsCount[0]]],
                                                               elementsCountOut)

        tdlcc1 = interpolateSampleCubicHermite([btd3[0][self._elementsCount[1]][self._elementsCount[0]],
                                                btd3[self._elementsCount[2]][self._elementsCount[1]][self._elementsCount[0]]],
                                               [[0.0, 0.0, 0.0]] * 2, pec, pxic, psfc)[0]
        tdlcc2 = interpolateSampleCubicHermite([btd1[0][self._elementsCount[1]][self._elementsCount[0]],
                                                btd1[self._elementsCount[2]][self._elementsCount[1]][self._elementsCount[0]]],
                                               [[0.0, 0.0, 0.0]] * 2, pec, pxic, psfc)[0]
        tdlec1 = interpolateSampleCubicHermite([btd3[0][self._elementsCount[1]][0],
                                                btd3[self._elementsCount[2]][self._elementsCount[1]][0]],
                                               [[0.0, 0.0, 0.0]] * 2, pec1, pxic1, psfc1)[0]
        tdlec2 = interpolateSampleCubicHermite([btd3[0][0][self._elementsCount[0]],
                                                btd3[self._elementsCount[2]][0][self._elementsCount[0]]],
                                               [[0.0, 0.0, 0.0]] * 2, pec2, pxic2, psfc2)[0]

        for n3 in range(1, self._elementsCount[2]):
            centre = txcc1[n3]
            txc1, td1c1 = self.generate1DPath([centre, txec1[n3]], [[-c for c in tdlcc1[n3]], tdlec1[n3]], self._elementsCount[0])
            txc2, td1c2 = self.generate1DPath([centre, txec2[n3]], [[-c for c in tdlcc2[n3]], tdlec2[n3]], self._elementsCount[1])
            ellipse = self.generateSurfaceUsingTwoCurves(centre, txc1, td1c1, txc2, td1c2)
            self.copyEllipseNodesToBifurcation(ellipse, n3)
=======

        # generate the armpit curves.
        elementsCountOut = 2

        # sample in between ellipses
        for n2 in range(self._elementsCount[1] + 1):
            for n1 in range(self._elementsCount[0] + 1):
                if btx[0][n2][n1]:
                    tx, td2, pe, pxi, psf = sampleCubicHermiteCurves(
                        [btx[0][n2][n1], btx[self._elementsCount[2]][n2][n1]],
                        [btd2[0][n2][n1], btd2[self._elementsCount[2]][n2][n1]],
                        elementsCountOut)
                    td1 = interpolateSampleCubicHermite(
                        [btd1[0][n2][n1], btd1[self._elementsCount[2]][n2][n1]],
                        [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
                    td3 = interpolateSampleCubicHermite(
                        [btd3[0][n2][n1], btd3[self._elementsCount[2]][n2][n1]],
                        [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]

                    for n3 in range(1, self._elementsCount[2]):
                        n3i = n3
                        self.px[n3][n2][n1] = tx[n3i]
                        self.pd1[n3][n2][n1] = td1[n3i]
                        self.pd2[n3][n2][n1] = td2[n3i]
                        self.pd3[n3][n2][n1] = td3[n3i]
>>>>>>> 9619665 (Modify axes and clean the code and make some classes)

    def smoothd2(self):
        """

        :return:
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3
        # smooth d2
        for n2 in range(self._elementsCount[1] + 1):
            for n1 in range(self._elementsCount[0] + 1):
                if btx[0][n2][n1]:
                    nx = []
                    nd1 = []
                    for n3 in range(self._elementsCount[2] + 1):
                        nx.append(btx[n3][n2][n1])
                        nd1.append(btd2[n3][n2][n1])
                    td2 = smoothCubicHermiteDerivativesLine(nx, nd1)
                    for n3 in range(self._elementsCount[2] + 1):
                        btd2[n3][n2][n1] = td2[n3]

    def generate1DPath(self, nx, nd1, elementsCountOut):
        """
        Given end nodes generate 1d path
        :return:
        """
        tx, td1, pe, pxi, psf = sampleCubicHermiteCurves(nx, nd1, elementsCountOut)
        return tx, td1

    def generateSurfaceUsingTwoCurves(self, centre, txc1, td1c1, txc2, td1c2):
        """
        Get major and minor curves to generate the rounded surface.
        :return:
        """

        # self._coreMinorRadii.append(
        #     (1 - self._shellProportion * self._elementsCountAcrossShell / elementsMinor) * self._minorRadii[n3])
        # self._coreMajorRadii.append(
        #     (1 - self._shellProportion * self._elementsCountAcrossShell / elementsMajor) * self._majorRadii[n3])

        # ratio = rx/self.elementsCountAcrossShell if self.elementsCountAcrossShell > 0 else 0
        # majorAxis = [d*(1 - ratio*(1-self.coreMajorRadius/self.majorRadius)) for d in self.majorAxis]
        # minorAxis = [d*(1 - ratio*(1-self.coreMinorRadius/self.minorRadius)) for d in self.minorAxis]
        centre_d1c1 = [c for c in td1c1[0]]
        centre_d1c2 = [c for c in td1c2[0]]
        majorAxis = [c*self._elementsCount[0] for c in centre_d1c2]
        minorAxis = [-c*self._elementsCount[1] for c in centre_d1c1]

        elementsCountAcrossShell = 0
        elementsCountAcrossTransition = 1
        shellProportion = 1.0
        coreMajorRadius = 1.0
        coreMinorRadius = 1.0
        elementsCountAround = self._elementsCount[0] + self._elementsCount[1] - 2
        ellipse = Ellipse2D(centre, majorAxis, minorAxis, 2*self._elementsCount[0], 2*self._elementsCount[1],
                            elementsCountAcrossShell, elementsCountAcrossTransition,
                            shellProportion, coreMajorRadius, coreMinorRadius,
                            ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

        shield = ellipse.getShield()
        ellipse.generateBase1DMesh(0)

        # Modify x and d1 for nodes on the end of curve 1 and 2.
        ellipse.px[self._elementsCount[1]][0] = txc1[-1]
        ellipse.pd3[self._elementsCount[1]][0] = td1c1[-1]
        ellipse.px[0][self._elementsCount[0]] = txc2[-1]
        ellipse.pd3[0][self._elementsCount[0]] = td1c2[-1]

        # Bend the surface according to its new curves
        # Modify x and d1 for the surface perimeter
        nx = []
        nd1 = []
        for n in range(elementsCountAround + 1):
            n1, n2 = shield.convertRimIndex(n)
            xi = n/elementsCountAround
            # ellipse.px[n2][n1][2] = (1 - xi) * ellipse.px[self._elementsCount[1]][0][2] + xi * \
            #                         ellipse.px[0][self._elementsCount[0]][2]
            x_p = vector.addVectors([vector.scaleVector(ellipse.px[self._elementsCount[1]][0], 1 - xi),
                                     vector.scaleVector(ellipse.px[0][self._elementsCount[0]], xi)], [1, 1])
            delta = [-ellipse.px[0][self._elementsCount[0]][c] + ellipse.px[self._elementsCount[1]][0][c]
                     for c in range(3)]
            normal = vector.normalise(vector.crossproduct3(td1c1[0], td1c2[0]))
            delta_p = vector.addVectors([x_p, ellipse.px[self._elementsCount[1]][self._elementsCount[0]]], [1, -1])
            delta_p_normalmag = vector.dotproduct(delta_p, normal)
            delta_p_normal = vector.scaleVector(normal, delta_p_normalmag)
            delnormag = vector.dotproduct(delta, normal)
            delnor = vector.scaleVector(normal, delnormag)
            if 0 < n < elementsCountAround:
                # ellipse.px[n2][n1] = [ellipse.px[n2][n1][c] +(1-xi) * delnor[c] for c in range(3)]
                ellipse.px[n2][n1] = vector.addVectors([ellipse.px[n2][n1], delta_p_normal], [1, 1])
            nx.append(ellipse.px[n2][n1])
            nd1.append(ellipse.pd1[n2][n1])
        td1 = smoothCubicHermiteDerivativesLine(nx, nd1)
        for n in range(1, elementsCountAround):
            n1, n2 = shield.convertRimIndex(n)
            ellipse.pd1[n2][n1] = td1[n]

        # for n in range(elementsCountAround + 1, 2*elementsCountAround+1):
        #     n1, n2 = shield.convertRimIndex(n)
        #     n1q1, n2q1 = shield.convertRimIndex(2*elementsCountAround - n)
        #     ellipse.px[n2][n1] = [-ellipse.px[n2q1][n1q1][0], ellipse.px[n2q1][n1q1][1], ellipse.px[n2q1][n1q1][2]]
        #     ellipse.pd1[n2][n1] = [ellipse.pd1[n2q1][n1q1][0], -ellipse.pd1[n2q1][n1q1][1], -ellipse.pd1[n2q1][n1q1][2]]
        #     ellipse.pd2[n2][n1] = [ellipse.pd2[n2q1][n1q1][0], ellipse.pd2[n2q1][n1q1][1], ellipse.pd2[n2q1][n1q1][2]]
        #     ellipse.pd3[n2][n1] = [-ellipse.pd3[n2q1][n1q1][0], ellipse.pd3[n2q1][n1q1][1], ellipse.pd3[n2q1][n1q1][2]]

        # recalculate the inner nodes.
        ellipse.setRimNodes()
        rscx, rscd1, rscd2, rscd3 = ellipse.createMirrorCurve()
        rscx2 = [ellipse.px[c][self._elementsCount[0]] for c in range(self._elementsCount[1]+1)]
        ellipse.createRegularRowCurves(rscx2, rscd1, rscd3)
        # ellipse.createRegularRowCurves(ellipse.px[:self._elementsCount[1]+1][self._elementsCount[0]], rscd1, rscd3)
        ellipse.createRegularColumnCurves()
        shield.getTriplePoints(0)
        ellipse.smoothTriplePointsCurves()
        ellipse.smoothTransitionRims()
        if ellipse.ellipseShape == EllipseShape.Ellipse_SHAPE_FULL:
            ellipse.generateNodesForUpperHalf()

        # smooth triple1
        nx = [ellipse.px[1][1], ellipse.px[0][1]]
        nd1 = [[-ellipse.pd1[1][1][c]-ellipse.pd3[1][1][c] for c in range(3)], ellipse.pd3[0][1]]
        td3 = smoothCubicHermiteDerivativesLine(nx, nd1)

        ellipse.pd3[0][1] = td3[1]

        for n2 in range(self._elementsCount[1] + 1):
            for n1 in range(self._elementsCount[0] + 1):
                if ellipse.px[n2][n1]:
                    ellipse.pd2[n2][n1] = vector.crossproduct3(ellipse.pd3[n2][n1], ellipse.pd1[n2][n1])

        return ellipse

    def copyEllipseNodesToTrifurcation(self, ellipse, n3):
        """
        Copy ellipse nodes to trifurcation
        :param ellipse:
        :return:
        """

        for n2 in range(self._elementsCount[1] + 1):
            for n1 in range(self._elementsCount[0] + 1):
                n2e, n1e = n2, n1
                self.px[n3][n2][n1] = ellipse.px[n2e][n1e]
                self.pd1[n3][n2][n1] = ellipse.pd1[n2e][n1e]
                self.pd2[n3][n2][n1] = ellipse.pd2[n2e][n1e]
                self.pd3[n3][n2][n1] = ellipse.pd3[n2e][n1e]


class BifurcationMesh:
    """
    Bifurction mesh generator.
    """

    def __init__(self, fieldModule, coordinates, region, centre, radii, right_leg_angle, left_leg_angle, part1=None):
        """
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        """
        # generate the mesh
        elementsCount = [2, 2, 5]
        self._elementsCount = elementsCount
        self._region = region
        self._centre = centre
        self._radii = radii
        self.right_leg_angle = right_leg_angle
        self.left_leg_angle = left_leg_angle
        self._part1 = part1

        self._coordinates = coordinates
        self._fieldmodule = fieldModule
        self._nodes = fieldModule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._mesh = fieldModule.findMeshByDimension(3)

        self.createBifurcationMesh3d(fieldModule, coordinates)

    def createBifurcationMesh3d(self, fieldmodule, coordinates):
        """
        Create a trifurcation.
        :param fieldmodule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :return: Final values of nextNodeIdentifier, nextElementIdentifier.
        """
        # assert (self._elementsCountAlong > 0), 'createCylinderMesh3d:  Invalid number of along elements'
        # assert (self._elementsCountAcrossMinor > 3), 'createCylinderMesh3d: Invalid number of across elements'
        # assert (self._elementsCountAcrossMinor % 2 == 0), 'createCylinderMesh3d: number of across elements' \
        #                                                   ' is not an even number'
        # assert (self._elementsCountAcrossMajor > 1), 'createCylinderMesh3d: Invalid number of up elements'
        # assert (self._cylinderShape in [self._cylinderShape.CYLINDER_SHAPE_FULL,
        #                                 self._cylinderShape.CYLINDER_SHAPE_LOWER_HALF]), \
        #     'createCylinderMesh3d: Invalid cylinder mode.'
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldmodule.findMeshByDimension(3)

        elementsCountAcrossShell = 0
        elementsCountAcrossTransition = 1
        shellProportion = 1.0
        coreMajorRadius = 1.0
        coreMinorRadius = 1.0
        # elementsCountAround = self._elementsCount[0] + self._elementsCount[1] - 2

        xa = [1, 0, 0]
        ya = [0, 1, 0]
        za = [0, 0, 1]
        radius1 = 1.6
        radius2 = 1.3
        radius3 = 0.8
        majorAxis = vector.scaleVector(ya, radius2)
        minorAxis = vector.scaleVector(xa, -radius1)
        centre1 = vector.addVectors([self._centre, vector.setMagnitude(za, -0.0)], [1, 1])

        ellipse1 = Ellipse2D(centre1, majorAxis, minorAxis, 2*self._elementsCount[0], 2*self._elementsCount[1],
                             elementsCountAcrossShell, elementsCountAcrossTransition,
                             shellProportion, coreMajorRadius, coreMinorRadius,
                             ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

        majorAxis = vector.scaleVector(ya, radius2)
        minorAxis = vector.scaleVector(za, -radius3)

        ellipse2 = Ellipse2D(centre1, majorAxis, minorAxis, 2*self._elementsCount[0], 2*self._elementsCount[1],
                             elementsCountAcrossShell, elementsCountAcrossTransition,
                             shellProportion, coreMajorRadius, coreMinorRadius,
                             ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

        elementsCount = [ellipse2.elementsCountAcrossMinor, ellipse2.elementsCountAcrossMajor, 1]

        shield = ellipse2.getShield()
        for n2 in range(elementsCount[1] + 1):
            for n1 in range(elementsCount[0] + 1):
                if n1 <= elementsCount[1]//2:
                    shield.px[0][n2][n1] = None

        radius1 = 0.8
        radius2 = 0.8
        hip_length = 1.6
        # angle = math.pi/2*0.7
        angle = self.left_leg_angle
        hip_distance = 1.4  # distance from belly button.

        kv = [0.0, 1.0, 0.0]
        x_hip_joint = vector.addVectors([centre1,
                                         vector.setMagnitude(vector.rotateVectorAroundVector(xa, kv, math.pi/4),
                                                             hip_distance)], [1, 1])
        Lv = vector.setMagnitude(vector.rotateVectorAroundVector(xa, kv, angle), hip_length)
        # Lv = vector.setMagnitude(vector.rotateVectorAroundVector(xa, kv, angle), hip_length)
        centre = vector.addVectors([x_hip_joint, [0, 0, 0]], [1, 1])
        # centre = vector.addVectors([centre1, Lv], [1, 1])
        minorAxis = vector.setMagnitude(vector.vectorRejection(vector.scaleVector(xa, -1), Lv), radius1)
        majorAxis = vector.setMagnitude(vector.crossproduct3(Lv, minorAxis), radius2)

        ellipse3 = Ellipse2D(centre, majorAxis, minorAxis, 2*self._elementsCount[0], 2*self._elementsCount[1],
                             elementsCountAcrossShell, elementsCountAcrossTransition,
                             shellProportion, coreMajorRadius, coreMinorRadius,
                             ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

        hip_left = Hip(elementsCount, 'left')

        pr = [hip_left.nodeId, hip_left.px, hip_left.pd1, hip_left.pd2, hip_left.pd3]
        p1 = [ellipse1.nodeId[0], ellipse1.px, ellipse1.pd1, ellipse1.pd2, ellipse1.pd3]
        p1v2 = [self._part1._shield.nodeId[0], self._part1._shield.px[0],
                self._part1._shield.pd1[0], self._part1._shield.pd2[0], self._part1._shield.pd3[0]]
        p2 = [ellipse2.nodeId[0], ellipse2.px, ellipse2.pd1, ellipse2.pd2, ellipse2.pd3]
        p3 = [ellipse3.nodeId[0], ellipse3.px, ellipse3.pd1, ellipse3.pd2, ellipse3.pd3]
        for n3 in range(elementsCount[2] + 1):
            for n2 in range(elementsCount[1] + 1):
                for n1 in range(elementsCount[0] + 1):
                    for i in range(5):
                        if n3 == 0:
                            pr[i][n3][n2][n1] = p3[i][n2][n1]
                        elif n3 == elementsCount[2]:
                            if n1 <= elementsCount[0]//2:
                                if i != 1:
                                    pr[i][n3][n2][n1] = p1v2[i][n2][n1]
                            else:
                                pr[i][n3][n2][n1] = p2[i][n2][n1]

        self.generateNodes(nodes, fieldmodule, coordinates, hip_left)
        # Add the common nodeIds
        for n3 in range(elementsCount[2] + 1):
            for n2 in range(elementsCount[1] + 1):
                for n1 in range(elementsCount[0] + 1):
                    for i in range(5):
                        if n3 == elementsCount[2]:
                            if n1 <= elementsCount[0]//2:
                                pr[i][n3][n2][n1] = p1v2[i][n2][n1]
        self.generateElements(mesh, fieldmodule, coordinates, hip_left)

        x_hip_joint = vector.addVectors([centre1,
                                         vector.setMagnitude(vector.rotateVectorAroundVector(
                                             vector.scaleVector(xa, -1), kv, -math.pi/4),
                                                             hip_distance)], [1, 1])
        Lv = vector.setMagnitude(vector.rotateVectorAroundVector(vector.scaleVector(xa, -1), kv, -self.right_leg_angle),
                                 hip_length)
        centre = vector.addVectors([x_hip_joint, [0, 0, 0]], [1, 1])
        # centre = vector.addVectors([centre1, Lv], [1, 1])
        minorAxis = vector.setMagnitude(vector.vectorRejection(vector.scaleVector(xa, -1), Lv), radius1)
        majorAxis = vector.setMagnitude(vector.crossproduct3(Lv, minorAxis), radius2)

        ellipse4 = Ellipse2D(centre, majorAxis, minorAxis, 2*self._elementsCount[0], 2*self._elementsCount[1],
                             elementsCountAcrossShell, elementsCountAcrossTransition,
                             shellProportion, coreMajorRadius, coreMinorRadius,
                             ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

        hip_right = Hip(elementsCount, 'right')

        pl = [hip_right.nodeId, hip_right.px, hip_right.pd1, hip_right.pd2, hip_right.pd3]
        p1 = [ellipse1.nodeId[0], ellipse1.px, ellipse1.pd1, ellipse1.pd2, ellipse1.pd3]
        p2 = [ellipse2.nodeId[0], ellipse2.px, ellipse2.pd1, ellipse2.pd2, ellipse2.pd3]
        p3 = [ellipse4.nodeId[0], ellipse4.px, ellipse4.pd1, ellipse4.pd2, ellipse4.pd3]
        for n3 in range(elementsCount[2] + 1):
            for n2 in range(elementsCount[1] + 1):
                for n1 in range(elementsCount[0] + 1):
                    for i in range(5):
                        if n3 == 0:
                            pl[i][n3][n2][n1] = p3[i][n2][n1]
                        elif n3 == elementsCount[2]:
                            if n1 > elementsCount[0]//2:
                                if i!= 1:
                                    pl[i][n3][n2][n1] = p1v2[i][n2][n1]
                            else:
                                if i != 1:
                                    n2r, n1r = n2, elementsCount[0] - n1
                                    pl[i][n3][n2][n1] = pr[i][n3][n2r][n1r]

        self.generateNodes(nodes, fieldmodule, coordinates, hip_right)
        for n3 in range(elementsCount[2] + 1):
            for n2 in range(elementsCount[1] + 1):
                for n1 in range(elementsCount[0] + 1):
                    for i in range(5):
                        if n3 == elementsCount[2]:
                            if n1 > elementsCount[0]//2:
                                pl[i][n3][n2][n1] = p1v2[i][n2][n1]
        self.generateElements(mesh, fieldmodule, coordinates, hip_right)

        pn = PathNodes(hip_left, [[radius1]*2, [radius1*0.8]*2], 6.0, [4, 4, 7], attach_bottom=False)
        path_list = pn.get_path_list()
        bc = BranchCylinder(self._region, self._mesh, self._nodes, self._fieldmodule, self._coordinates,
                            path_list, [4, 4, 7], hip_left, attach_bottom=False)
        cylinder = bc.get_cylinder()

        pn = PathNodes(hip_right, [[radius1]*2, [radius1*0.8]*2], 6.0, [4, 4, 7], attach_bottom=False)
        path_list = pn.get_path_list()
        bc = BranchCylinder(self._region, self._mesh, self._nodes, self._fieldmodule, self._coordinates,
                            path_list, [4, 4, 7], hip_right, attach_bottom=False)
        cylinder = bc.get_cylinder()

    def generateNodes(self, nodes, fieldModule, coordinates, part_structure):
        """
        Create cylinder nodes from coordinates.
        :param nodes: nodes from coordinates.
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self.topologygenerateNodes(fieldModule, coordinates, nodeIdentifier, part_structure)
        self._endNodeIdentifier = nodeIdentifier

    def topologygenerateNodes(self, fieldmodule, coordinates, startNodeIdentifier, part_structure):
        """
        Create shield nodes from coordinates.
         """
        nodeIdentifier = startNodeIdentifier
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        cache = fieldmodule.createFieldcache()

        for n3 in range(1 + 1):  # TODO should change to number of elements along, major, minor.
            for n2 in range(4 + 1):
                for n1 in range(4 + 1):
                    if part_structure.px[n3][n2][n1]:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        part_structure.nodeId[n3][n2][n1] = nodeIdentifier
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1,
                                                      part_structure.px[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1,
                                                      part_structure.pd1[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1,
                                                      part_structure.pd2[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1,
                                                      part_structure.pd3[n3][n2][n1])
                        nodeIdentifier += 1

        return nodeIdentifier

    def generateElements(self, mesh, fieldModule, coordinates, part_structure):
        """
        Create cylinder elements from nodes.
        :param mesh:
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        self._startElementIdentifier = elementIdentifier
        elementIdentifier = self.topologygenerateElements(fieldModule, coordinates, elementIdentifier, part_structure,
                                                          [])
        self._endElementIdentifier = elementIdentifier

    def topologygenerateElements(self, fieldmodule, coordinates, startElementIdentifier, part_structure, meshGroups=[]):
        """
        Create shield elements from nodes.
        :param fieldmodule: Zinc fieldmodule to create elements in.
        :param coordinates: Coordinate field to define.
        :param startElementIdentifier: First element identifier to use.
        :param meshGroups: Zinc mesh groups to add elements to.
        :return: next elementIdentifier.
         """
        elementIdentifier = startElementIdentifier
        useCrossDerivatives = False
        mesh = fieldmodule.findMeshByDimension(3)

        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftNoCrossDerivatives()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementtemplate1 = mesh.createElementtemplate()
        elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)

        # isEven = (self.elementsCountAcross % 2) == 0
        e1a = 0
        e1b = e1a + 1
        e1z = part_structure._elementsCount[0] - 1
        e1y = e1z - 1
        e2a = 0
        e2b = e2a + 1
        e2c = e2a + 2
        e2z = part_structure._elementsCount[1]-1
        e2y = e2z - 1
        # e2x = e2z - 2
        for e3 in range(part_structure._elementsCount[2]):
            for e2 in range(part_structure._elementsCount[1]):
                for e1 in range(part_structure._elementsCount[0]):
                    eft1 = eft
                    scalefactors = None
                    nids = [part_structure.nodeId[e3][e2][e1], part_structure.nodeId[e3][e2 + 1][e1],
                            part_structure.nodeId[e3+1][e2][e1], part_structure.nodeId[e3+1][e2 + 1][e1],
                            part_structure.nodeId[e3][e2][e1 + 1], part_structure.nodeId[e3][e2 + 1][e1 + 1],
                            part_structure.nodeId[e3+1][e2][e1 + 1], part_structure.nodeId[e3+1][e2 + 1][e1 + 1]]

                    if (e2 < e2b) or (e2 > e2y):
                        if (e1 < e1b) or (e1 > e1y):
                            continue  # no element due to triple point closure
                        if (e2 < e2a) or (e2 > e2z):
                            if e2 < e2a:
                                nids = [part_structure.nodeId[e3][e2+1][e1], part_structure.nodeId[e3][e2+1][e1+1],
                                        part_structure.nodeId[e3+1][e2+1][e1], part_structure.nodeId[e3+1][e2+1][e1+1],
                                        part_structure.nodeId[e3][e2][e1], part_structure.nodeId[e3][e2][e1+1],
                                        part_structure.nodeId[e3+1][e2][e1],  part_structure.nodeId[e3+1][e2][e1+1]]
                            elif e2 > e2z:
                                nids = [part_structure.nodeId[e3][e2][e1+1], part_structure.nodeId[e3][e2][e1],
                                        part_structure.nodeId[e3+1][e2][e1+1], part_structure.nodeId[e3+1][e2][e1],
                                        part_structure.nodeId[e3][e2+1][e1+1], part_structure.nodeId[e3][e2+1][e1],
                                        part_structure.nodeId[e3+1][e2+1][e1+1], part_structure.nodeId[e3+1][e2+1][e1]]
                        elif (e2 == e2a) or (e2 == e2z):
                            # bottom and top row elements
                            if e2 == e2a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if part_structure._side == 'left':
                                    # if e3,e2,e1 = 0, 0, 2
                                    if e1 == e1y:
                                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [])])

                                if part_structure._side == 'left':
                                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    if e1 == e1y:
                                        remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                elif part_structure._side == 'right':
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])

                                if (e1 == e1b) or (e1 == e1y):
                                    # map bottom triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                    else:

                                        remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                            elif e2 == e2z:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if isinstance(part_structure, BaseLeg):
                                    if part_structure._shoulder and e3 == 0 and e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS2, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                # TODO if e3 == -1
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                if part_structure._side == 'left':
                                    if e1 == e1y:
                                        remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [])])
                                elif part_structure._side == 'right':
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D2_DS1DS2,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [2, 6, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS3, [])])

                                if (e1 == e1b) or (e1 == e1y):
                                    # map top triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [])])

                    elif (e2 == e2b) or (e2 == e2y):
                        if (e1 <= e1a) or (e1 >= e1z):
                            if e1 < e1a:
                                e2r = e1
                                if e2 == e2b:
                                    nids = [part_structure.nodeId[e3][e2c][e1+1],
                                            part_structure.nodeId[e3][e2r+1][e1b],
                                            part_structure.nodeId[e3+1][e2c][e1+1],
                                            part_structure.nodeId[e3+1][e2r+1][e1b],
                                            part_structure.nodeId[e3][e2c][e1],
                                            part_structure.nodeId[e3][e2r][e1b],
                                            part_structure.nodeId[e3+1][e2c][e1],
                                            part_structure.nodeId[e3+1][e2r][e1b]]
                                if e2 == e2y:
                                    e2r = 2*part_structure._elementsCount[1] - e1-1
                                    nids = [part_structure.nodeId[e3][e2r][e1b],
                                            part_structure.nodeId[e3][e2y][e1+1],
                                            part_structure.nodeId[e3+1][e2r][e1b],
                                            part_structure.nodeId[e3+1][e2y][e1+1],
                                            part_structure.nodeId[e3][e2r+1][e1b],
                                            part_structure.nodeId[e3][e2y][e1],
                                            part_structure.nodeId[e3+1][e2r+1][e1b],
                                            part_structure.nodeId[e3+1][e2y][e1]]
                            elif e1 == e1a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e2 == e2b:
                                    if e3 == part_structure._elementsCount[2] // 2 + 1:
                                        e3r = e3-1  # to join upper leg with the lower leg.
                                        nids[0] = part_structure.nodeId[e3r][e2a][e1b]
                                        nids[2] = part_structure.nodeId[e3+1][e2a][e1b]
                                        nids[1] = part_structure.nodeId[e3r][e2 + 1][e1]
                                        nids[4] = part_structure.nodeId[e3r][e2][e1 + 1]
                                        nids[5] = part_structure.nodeId[e3r][e2 + 1][e1 + 1]
                                    else:
                                        nids[0] = part_structure.nodeId[e3][e2a][e1b]
                                        nids[2] = part_structure.nodeId[e3+1][e2a][e1b]
                                    tripleN = [5, 7] if part_structure._side == 'left' else [5]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2y:
                                    nids[1] = part_structure.nodeId[e3][e2z+1][e1b]
                                    nids[3] = part_structure.nodeId[e3+1][e2z+1][e1b]
                                    tripleN = [6, 8] if part_structure._side == 'left' else [6]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1])])
                                if part_structure._side == 'right':
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    if e2 == e2b:
                                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, []),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                    elif e2 == e2y:
                                        remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                (Node.VALUE_LABEL_D_DS3, [1])])
                                elif part_structure._side == 'left':
                                    remapEftNodeValueLabel(eft1, [1, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS3, [1])])
                            elif e1 == e1z:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if e2 == e2b:
                                    nids[4] = part_structure.nodeId[e3][e2a][e1z]
                                    nids[6] = part_structure.nodeId[e3+1][e2a][e1z]
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [1]),
                                                            (Node.VALUE_LABEL_D_DS3, [])])
                                    if part_structure._side == 'left':
                                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                elif e2 == e2y:
                                    nids[5] = part_structure.nodeId[e3][e2z+1][e1z]
                                    nids[7] = part_structure.nodeId[e3+1][e2z+1][e1z]
                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    # if e3,e2,e1 = 0 2 3
                                    if part_structure._side == 'left':
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [-1.0]
                                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                            elif e1 > e1z:
                                e2r = part_structure._elementsCount[0] - e1
                                if e2 == e2b:
                                    nids = [part_structure.nodeId[e3][e2r][e1z],
                                            part_structure.nodeId[e3][e2c][e1],
                                            part_structure.nodeId[e3+1][e2r][e1z],
                                            part_structure.nodeId[e3+1][e2c][e1],
                                            part_structure.nodeId[e3][e2r-1][e1z],
                                            part_structure.nodeId[e3][e2c][e1+1],
                                            part_structure.nodeId[e3+1][e2r-1][e1z],
                                            part_structure.nodeId[e3+1][e2c][e1+1]]
                                elif e2 == e2y:
                                    e2r = e2z+e1-e1z
                                    nids[1] = part_structure.nodeId[e3][e2r][e1z]
                                    nids[3] = part_structure.nodeId[e3+1][e2r][e1z]
                                    nids[5] = part_structure.nodeId[e3][e2r+1][e1z]
                                    nids[7] = part_structure.nodeId[e3+1][e2r+1][e1z]
                        elif e1 == e1b:
                            if e2 == e2b:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if part_structure._side == 'right':
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                            else:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if part_structure._side == 'right':
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])

                        elif e1 == e1y:
                            if e2 == e2b or e2 == e2y:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if part_structure._side == 'left':
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])

                    else:
                        if e1 < e1a:
                            nids = [part_structure.nodeId[e3][e2 + 1][e1 + 1],
                                    part_structure.nodeId[e3][e2][e1 + 1],
                                    part_structure.nodeId[e3+1][e2 + 1][e1 + 1],
                                    part_structure.nodeId[e3+1][e2][e1 + 1],
                                    part_structure.nodeId[e3][e2 + 1][e1],
                                    part_structure.nodeId[e3][e2][e1],
                                    part_structure.nodeId[e3+1][e2 + 1][e1],
                                    part_structure.nodeId[e3+1][e2][e1]]
                        elif e1 == e1a:
                            # map left column elements
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [-1.0]
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS3, [1])])
                            if e3 == 0 and part_structure._shoulder:
                                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1])])

                        elif e1 == e1b:
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            if e3 == 0 and part_structure._shoulder:
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS2, [])])

                    if not all(nids):
                        continue

                    if eft1 is not eft:
                        elementtemplate1.defineField(coordinates, -1, eft1)
                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                    else:
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                    result2 = element.setNodesByIdentifier(eft1, nids)
                    if scalefactors:
                        result3 = element.setScaleFactors(eft1, scalefactors)
                    else:
                        result3 = 7
                    part_structure.elementId[e3][e2][e1] = elementIdentifier
                    elementIdentifier += 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

        return elementIdentifier


class Hip:
    def __init__(self, elementsCount, side):
        self._elementsCount = elementsCount
        self._side = side

        elementsCount = elementsCount

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                   for c in range(elementsCount[2] + 1)]
        self._px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                       for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])]
                          for c in range(elementsCount[2])]


class BoxPart:
    def __init__(self, elementsCount, torso, shoulder_left, shoulder_right, neck):
        # self._elementsCount = [2, 4, 2]
        # self._elementsCount = [4, 6, 4]
        self._elementsCount = elementsCount
        self._joining_box = True

        # elementsCount = [2, 4, 2]
        # elementsCount = [4, 6, 4]
        elementsCount = elementsCount

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                   for c in range(elementsCount[2] + 1)]
        self._px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                    for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)]
                       for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])]
                          for c in range(elementsCount[2])]

        torso_elements_count = [len(torso.px[0][0])-1, len(torso.px[0])-1, len(torso.px)-1]
        p = [self.nodeId, self._px, self.pd1, self.pd2, self.pd3]
        p1 = [torso.nodeId, torso.px, torso.pd1, torso.pd2, torso.pd3]
        p2 = [neck.nodeId, neck.px, neck.pd1, neck.pd2, neck.pd3]
        p3 = [shoulder_left.nodeId, shoulder_left.px, shoulder_left.pd1, shoulder_left.pd2, shoulder_left.pd3]
        p4 = [shoulder_right.nodeId, shoulder_right.px, shoulder_right.pd1, shoulder_right.pd2, shoulder_right.pd3]

        for n3 in range(elementsCount[2] + 1):
            for n2 in range(elementsCount[1] + 1):
                for n1 in range(elementsCount[0] + 1):
                    for i in range(5):
                        if (n2 == 0 or n2 == elementsCount[1]) and n1 == 0:
                            n2p, n1p = n2, n1 + 1
                        elif (n2 == 0 or n2 == elementsCount[1]) and n1 == torso_elements_count[0] - 2:
                            n2p, n1p = n2, n1 + 1
                        else:
                            n2p, n1p = n2, n1 + 1

                        if n3 == 0:
                            p[i][n3][n2][n1] = p1[i][torso_elements_count[2]][n2p][n1p]
                        elif n3 == elementsCount[2]:
                            p[i][elementsCount[2]][n2][n1] = p2[i][0][n2p][n1p]
                        else:
                            if n1 == 0:
                                p[i][n3][n2][0] = p3[i][0][n2][n3+1]
                            elif n1 == elementsCount[0]:
                                p[i][n3][n2][elementsCount[0]] = p4[i][0][n3+1][n2]

        elementsCountOut = elementsCount[2]
        for n2 in range(elementsCount[1] + 1):
            for n1 in range(1, elementsCount[0]):
                nx = [self._px[0][n2][n1], self._px[elementsCount[2]][n2][n1]]
                nd1 = [self.pd2[0][n2][n1], self.pd2[elementsCount[2]][n2][n1]]
                tx, td2, pe, pxi, psf = sampleCubicHermiteCurves(nx, nd1, elementsCountOut)
                td1 = interpolateSampleCubicHermite(
                    [self.pd1[0][n2][n1], self.pd1[elementsCount[2]][n2][n1]],
                    [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
                td3 = interpolateSampleCubicHermite(
                    [self.pd3[0][n2][n1], self.pd3[elementsCount[2]][n2][n1]],
                    [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
                for n3 in range(1, elementsCount[2]):
                    self.px[n3][n2][n1] = tx[n3]
                    self._px[n3][n2][n1] = tx[n3]
                    self.pd1[n3][n2][n1] = td1[n3]
                    self.pd2[n3][n2][n1] = td2[n3]
                    self.pd3[n3][n2][n1] = td3[n3]
