"""
Utility functions for generating a 3-D solid bifurcation.
"""
import copy
from enum import Enum
from scaffoldmaker.utils import vector, geometry
import math
from opencmiss.zinc.field import Field
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from scaffoldmaker.utils.shieldmesh import ShieldMesh2D, ShieldShape2D, ShieldRimDerivativeMode
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite, \
    smoothCubicHermiteDerivativesLine, interpolateSampleLinear
from opencmiss.zinc.node import Node
from scaffoldmaker.utils.mirror import Mirror
from scaffoldmaker.meshtypes.meshtype_1d_path1 import extractPathParametersFromRegion
from scaffoldmaker.utils.cylindermesh import Ellipse2D, EllipseShape, CylinderCentralPath, CylinderShape, CylinderEnds, CylinderMesh
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from opencmiss.zinc.element import Element
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues


class BifurcationMesh:
    """
    Bifurction mesh generator.
    """

    def __init__(self, fieldModule, coordinates, region, torso_radius, left_arm_radius, right_arm_radius, neck_radius,
                 shoulder_height, neck_height, right_arm_angle, right_arm_length):
        """
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        """
        # generate the mesh
        elementsCount = [2, 2, 5]
        self._elementsCount = elementsCount
        self._region = region

        self.torso_radius = torso_radius
        self.left_arm_radius = left_arm_radius
        self.right_arm_radius = right_arm_radius
        self.neck_radius = neck_radius
        self.shoulder_height = shoulder_height
        self.neck_height = neck_height
        self.right_arm_angle = right_arm_angle
        self.right_arm_length = right_arm_length


        self.createBifurcationMesh3d(fieldModule, coordinates)

    def createBifurcationMesh3d(self, fieldmodule, coordinates):
        """
        Create a bifurcation.
        :param fieldModule: Zinc fieldModule to create elements in.
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

        nodeparams1, nodeparams2, nodeparams3, nodeparams4 = self._get_node_params()

        bottom_part = BaseLeg(self._elementsCount, nodeparams1)
        self.generateNodes(nodes, fieldmodule, coordinates, bottom_part)
        self.generateElements(mesh, fieldmodule, coordinates, bottom_part)

        shoulder_part = BaseLeg(self._elementsCount, nodeparams2)
        shoulder_part._shoulder = True
        self.remove_duplicate_nodes_from_shoulder(shoulder_part)
        self.generateNodes(nodes, fieldmodule, coordinates, shoulder_part)
        self.join_shoulder_to_bottom_part(shoulder_part, bottom_part)
        self.generateElements(mesh, fieldmodule, coordinates, shoulder_part)

        joining_box = JoiningBox([1, 4, 1], self.shoulder_height)
        self.generateNodes(nodes, fieldmodule, coordinates, joining_box)
        self.join_box_to_bottom_and_shoulder(joining_box, bottom_part, shoulder_part)
        self.generateElements(mesh, fieldmodule, coordinates, joining_box)

        shoulder_part_left = BaseLeg(self._elementsCount, nodeparams3)
        shoulder_part_left._shoulder_left = True
        self.remove_duplicate_nodes_from_shoulder(shoulder_part_left, 1)
        self.generateNodes(nodes, fieldmodule, coordinates, shoulder_part_left)
        self.join_shoulder_to_bottom_part(shoulder_part_left, bottom_part, 1)
        self.generateElements(mesh, fieldmodule, coordinates, shoulder_part_left)

        neck_part = BaseLeg(self._elementsCount, nodeparams4)
        neck_part._neck = True
        self.remove_duplicate_nodes_from_neck(neck_part)
        self.generateNodes(nodes, fieldmodule, coordinates, neck_part)
        self.join_neck_to_torso(neck_part, shoulder_part, shoulder_part_left)
        self.generateElements(mesh, fieldmodule, coordinates, neck_part)


        joining_box_right = JoiningBoxRight([1, 4, 1], bottom_part, shoulder_part_left, joining_box)
        self.generateElements(mesh, fieldmodule, coordinates, joining_box_right)

        joining_box_3 = JoiningBox3([1, 4, 1], joining_box, neck_part)
        self.generateElements(mesh, fieldmodule, coordinates, joining_box_3)

        joining_box_4 = JoiningBox4([1, 4, 1], joining_box_3, neck_part, joining_box_right)
        self.generateElements(mesh, fieldmodule, coordinates, joining_box_4)


        # shoulder_connecting_to_box = CylinderConnectingToBox(shoulder_part, [0, 0], -1)
        # self.generateNodes(nodes, fieldmodule, coordinates, shoulder_connecting_to_box)
        # self.joint_shoulder_joint_to_cylinder_and_box(shoulder_connecting_to_box, joining_box, shoulder_part, [1, 0], 0)
        # self.generateElements(mesh, fieldmodule, coordinates, shoulder_connecting_to_box)
        #
        # bottom_connecting_to_box = CylinderConnectingToBox(bottom_part, [1, 2], 1)
        # self.generateNodes(nodes, fieldmodule, coordinates, bottom_connecting_to_box)
        # self.joint_shoulder_joint_to_cylinder_and_box(bottom_connecting_to_box, joining_box, bottom_part, [0, 2], 1)
        # self.generateElements(mesh, fieldmodule, coordinates, bottom_connecting_to_box)
        #
        # centre = [-6.580413981734434e-01, 5.756093176338770e-02, 2.797218065767146e+00]
        # outer_point = [4.417e-01, 4.174e-02, 3.897e+00]
        #
        # p1 = [-1.12e+00, 9.29e-02, 3.22e+00]
        # p2 =  [-1.62e+00, 1.061e-01, 3.72e+00]
        # p3 = [-2.12e+00, 1.14e-01, 4.26e+00]
        # p4 = [-2.605e+00, 1.22e-01, 4.76e+00]
        # d11 = vector.addVectors([p1, centre], [1, -1])
        # d21 = vector.addVectors([outer_point, centre], [1, -1])
        # d31 = vector.setMagnitude(vector.crossproduct3(d11, d21), -1.0)


        # centralPath = ScaffoldPackage(MeshType_1d_path1, {
        #     'scaffoldSettings': {
        #         'Coordinate dimensions': 3,
        #         'D2 derivatives': True,
        #         'D3 derivatives': True,
        #         'Length': 3.0,
        #         'Number of elements': 4
        #     },
        #     'meshEdits': exnodeStringFromNodeValues(
        #         [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
        #          Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
        #             [[-6.58e-01, 5.75e-02, 2.797e+00], d11, d31, [0.0, 0.0, 0.0], d21, [0.0, 0.0, 0.0]],
        #             [[-1.12e+00, 9.29e-02, 3.22e+00], d11, d31, [0.0, 0.0, 0.0], d21, [0.0, 0.0, 0.0]],
        #             [[-1.62e+00, 1.061e-01, 3.72e+00], d11, d31, [0.0, 0.0, 0.0], d21, [0.0, 0.0, 0.0]],
        #             [[-2.12e+00, 1.14e-01, 4.26e+00], d11, d31, [0.0, 0.0, 0.0], d21, [0.0, 0.0, 0.0]],
        #             [[-2.605e+00, 1.22e-01, 4.76e+00], d11, d31, [0.0, 0.0, 0.0], d21, [0.0, 0.0, 0.0]]
        #         ])
        # })
        #
        # cylinderCentralPath = CylinderCentralPath(self._region, centralPath, 5)
        #
        # cylinderShape = CylinderShape.CYLINDER_SHAPE_FULL
        #
        # base = CylinderEnds(4, 4, 0, 1, 1.0,
        #                     centre, cylinderCentralPath.alongAxis[0], cylinderCentralPath.majorAxis[0],
        #                     cylinderCentralPath.minorRadii[0])
        # torso_cylinder = CylinderMesh(fieldmodule, coordinates, 5, base,
        #                          cylinderShape=cylinderShape,
        #                          cylinderCentralPath=cylinderCentralPath, useCrossDerivatives=False)

        # joining_torso = JoiningTorso([4, 4, 1])
        # self.join_to_torso(joining_torso, torso_cylinder._shield, shoulder_connecting_to_box, bottom_connecting_to_box)
        # print(joining_torso.nodeId)
        # self.generateElements(mesh, fieldmodule, coordinates, joining_torso)
        #
        #
        # torso_part = torso_cylinder._shield

    def _get_node_params(self):
        """

        :return:
        """
        armpit = [1.2, 0.0, 1.0]
        x_bottom_base_centre = [0.0, 0.0, 0.0]
        x_bottom_base_curve1 = [self.torso_radius, 0.0, 0.0]
        x_bottom_base_curve2 = [0.0, self.torso_radius, 0.0]
        d1_bottom_base_curve1 = [[self.torso_radius/self._elementsCount[0], 0.0, 0.0],
                                 [self.torso_radius/self._elementsCount[0], 0.0, 0.0]]
        d1_bottom_base_curve2 = [[0.0, 1/self._elementsCount[1], 0.0],
                                 [0.0, 1/self._elementsCount[1], 0.0]]
        x_bottom_end_centre = [0.0, 0.0, 1.4]
        x_bottom_end_curve1 = armpit
        x_bottom_end_curve2 = [0.0, 1.0, 1.4]
        d1_bottom_end_curve1 = [[1/self._elementsCount[0], 0.0, 0.0],
                                [0.5*0.7071, 0.0, -0.5*0.7071]]
        d1_bottom_end_curve2 = [[0.0, 1/self._elementsCount[1], 0.0],
                                [0.0, 1/self._elementsCount[1], 0.0]]

        nodeparams1 = [[x_bottom_base_centre, x_bottom_base_curve1, x_bottom_base_curve2, d1_bottom_base_curve1,
                       d1_bottom_base_curve2],
                      [x_bottom_end_centre, x_bottom_end_curve1, x_bottom_end_curve2, d1_bottom_end_curve1,
                       d1_bottom_end_curve2]]

        x_shoulder_base_centre = [0.75, 0.0, self.shoulder_height]
        x_shoulder_base_curve1 = armpit
        x_shoulder_base_curve2 = [0.75, 1.0, self.shoulder_height]
        d1_shoulder_base_curve1 = [[0.0, 0.0, -1 / self._elementsCount[1]], [0.5 * 0.7071, 0.0, -0.5 * 0.7071]]
        d1_shoulder_base_curve2 = [[0.0, 1 / self._elementsCount[0], 0.0], [0.0, 1 / self._elementsCount[0], 0.0]]
        kv = [0.0, 1.0, 0.0]
        cv = [self.right_arm_length, 0.0, 0.0]
        cev = vector.rotateVectorAroundVector(cv, kv, self.right_arm_angle)
        ce = vector.addVectors([cev, [0.0, 0.0, self.shoulder_height]], [1, 1])
        x_shoulder_end_centre = ce
        # x_shoulder_end_curve1 = [1.7, 0.0, self.shoulder_height - self.left_arm_radius]
        x_shoulder_end_curve1 = vector.addVectors([ce, vector.setMagnitude(vector.crossproduct3(kv, cev), self.right_arm_radius)], [1, 1])
        # x_shoulder_end_curve2 = [1.7, self.left_arm_radius, self.shoulder_height]
        x_shoulder_end_curve2 = vector.addVectors([ce, vector.setMagnitude(kv, self.right_arm_radius)], [1, 1])
        # d1_shoulder_end_curve1 = [[0.0, 0.0, -self.left_arm_radius / self._elementsCount[1]], [0.0, 0.0, -self.left_arm_radius / self._elementsCount[1]]]
        d1_shoulder_end_curve1 = [vector.setMagnitude(vector.crossproduct3(kv, cev), self.right_arm_radius/self._elementsCount[1]),
                                  vector.setMagnitude(vector.crossproduct3(kv, cev), self.right_arm_radius/self._elementsCount[1])]
        d1_shoulder_end_curve2 = [[0.0, self.left_arm_radius / self._elementsCount[0], 0.0], [0.0, self.left_arm_radius / self._elementsCount[0], 0.0]]

        nodeparams2 = [[x_shoulder_base_centre, x_shoulder_base_curve1, x_shoulder_base_curve2, d1_shoulder_base_curve1,
                       d1_shoulder_base_curve2],
                      [x_shoulder_end_centre, x_shoulder_end_curve1, x_shoulder_end_curve2, d1_shoulder_end_curve1,
                       d1_shoulder_end_curve2]]

        x_shoulder_base_centre = [-0.75, 0.0, self.shoulder_height]
        x_shoulder_base_curve2 = [-armpit[0], armpit[1], armpit[2]]
        x_shoulder_base_curve1 = [-0.75, 1.0, self.shoulder_height]
        d1_shoulder_base_curve2 = [[-0.0, 0.0, -1 / self._elementsCount[1]], [-0.5 * 0.7071, 0.0, -0.5 * 0.7071]]
        d1_shoulder_base_curve1 = [[-0.0, 1 / self._elementsCount[0], 0.0], [-0.0, 1 / self._elementsCount[0], 0.0]]
        x_shoulder_end_centre = [-1.7, 0.0, self.shoulder_height]
        x_shoulder_end_curve2 = [-1.7, 0.0, self.shoulder_height - self.right_arm_radius]
        x_shoulder_end_curve1 = [-1.7, self.right_arm_radius, self.shoulder_height]
        d1_shoulder_end_curve2 = [[-0.0, 0.0, -self.right_arm_radius / self._elementsCount[1]], [0.0, 0.0, -self.right_arm_radius / self._elementsCount[1]]]
        d1_shoulder_end_curve1 = [[-0.0, self.right_arm_radius / self._elementsCount[0], 0.0], [0.0, self.right_arm_radius / self._elementsCount[0], 0.0]]

        nodeparams3 = [[x_shoulder_base_centre, x_shoulder_base_curve1, x_shoulder_base_curve2, d1_shoulder_base_curve1,
                       d1_shoulder_base_curve2],
                      [x_shoulder_end_centre, x_shoulder_end_curve1, x_shoulder_end_curve2, d1_shoulder_end_curve1,
                       d1_shoulder_end_curve2]]

        x_neck_base_centre = [0.0, 0.0, 2 * self.shoulder_height - 1.0 - self.right_arm_radius/2]
        x_neck_base_curve2 = [0.0, 1.0, 2 * self.shoulder_height - 1.0 - self.right_arm_radius/2]
        x_neck_base_curve1 = [1.2, 0.0, 2 * self.shoulder_height - 1.0]
        d1_neck_base_curve2 = [[0.0, 1 / self._elementsCount[0], 0.0], [0.0, 1 / self._elementsCount[0], 0.0]]
        d1_neck_base_curve1 = [[1 / self._elementsCount[0], 0.0, 0.0], [0.5, 0.0, 0.5]]
        x_neck_end_centre = [0.0, 0.0, self.neck_height]
        x_neck_end_curve2 = [0.0, self.neck_radius, self.neck_height]
        x_neck_end_curve1 = [self.neck_radius, 0.0, self.neck_height]
        d1_neck_end_curve2 = [[0.0, self.neck_radius / self._elementsCount[1], 0.0], [0.0, self.neck_radius / self._elementsCount[1], 0.0]]
        d1_neck_end_curve1 = [[self.neck_radius / self._elementsCount[0], 0.0, 0.0], [self.neck_radius / self._elementsCount[0], 0.0, 0.0]]

        nodeparams4 = [[x_neck_base_centre, x_neck_base_curve1, x_neck_base_curve2, d1_neck_base_curve1,
                       d1_neck_base_curve2],
                      [x_neck_end_centre, x_neck_end_curve1, x_neck_end_curve2, d1_neck_end_curve1,
                       d1_neck_end_curve2]]

        return nodeparams1, nodeparams2, nodeparams3, nodeparams4

    def join_to_torso(self, joining_torso, torso, shoulder_joint, bottom_joint):
        """

        :param joining_torso:
        :param shoulder_part:
        :param bottom_part:
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

    def joint_shoulder_joint_to_cylinder_and_box(self, shoulder_connecting_to_box, joining_box, cylinder_part, cidxs, bidx):
        """

        :param shoulder_connecting_to_box:
        :param joining_box:
        :param shoulder_part:
        :return:
        """
        for n2 in range(shoulder_connecting_to_box._elementsCount[1] + 1):
            for n1 in range(shoulder_connecting_to_box._elementsCount[0]//2, shoulder_connecting_to_box._elementsCount[0] + 1):
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
        :param bottom_part:
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
            # if c:
            #     n2, n1 = n1, n2
            # return (n2 == 0 and n1 == 1) or (n2 == 1 and n1 == 1) or (n2 == 2 and n1 == 0) or\
            #                 (n2 == 2 and n1 == 1) or (n2 == 3 and n1 == 1) or (n2 == 4 and n1 == 1)
            if c:
                return n2 == 0 or n2 == 1
            else:
                return n1 == 0 or n1 == 1

        def index(n2, n1):
            if c:
                if n2 == 0 and n1 == 1:
                    return n1 - 1, 4 - n2 - 1
                if n2 == 0 and n1 == 3:
                    return n1 + 1, 4 - n2 - 1
                else:
                    return n1, 4 - n2
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

    def copyBaseLeg2Bifurcation(self, baseleg, idx):
        """

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
                        if (n2 == 0 and n1 == 1) or (n2 == 1 and n1 == 1) or (n2 == 2 and n1 == 0) or (n2 == 2 and n1 == 1):
                            self.px[n3s][n2][n1] = None
                            self.pd1[n3s][n2][n1] = None
                            self.pd2[n3s][n2][n1] = None
                            self.pd3[n3s][n2][n1] = None


    def remove_duplicate_nodes_from_neck(self, neck_part):
        """

        :param shoulder_part:
        :param bottom_part:
        :return:
        """
        for n2 in range(neck_part._elementsCount[1] + 1):
            for n1 in range(neck_part._elementsCount[0] + 1):
                if n1 != 2:
                    neck_part.px[0][n2][n1] = None
                    neck_part.pd1[0][n2][n1] = None
                    neck_part.pd2[0][n2][n1] = None
                    neck_part.pd3[0][n2][n1] = None

    def join_neck_to_torso(self, neck_part, shoulder_part, shoulder_part_left):
        """

        :param shoulder_part:
        :param bottom_part:
        :return:
        """
        for n2 in range(neck_part._elementsCount[1] + 1):
            for n1 in range(neck_part._elementsCount[0] + 1):
                if n1 < 2:
                    n1s = 4 - n1
                    neck_part.nodeId[0][n2][n1] = shoulder_part.nodeId[0][n2][n1s]
                    neck_part.px[0][n2][n1] = shoulder_part.px[0][n2][n1s]
                    neck_part.pd1[0][n2][n1] = shoulder_part.pd1[0][n2][n1s]
                    neck_part.pd2[0][n2][n1] = shoulder_part.pd2[0][n2][n1s]
                    neck_part.pd3[0][n2][n1] = shoulder_part.pd3[0][n2][n1s]
                elif n1 > 2:
                    if n2 == 0 or n2 == 4:
                        if n1 == 3:
                            n2s = n1+1
                            if n2 == 0:
                                n1s = n2 + 1
                            if n2 == 4:
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
        elementIdentifier = self.topologygenerateElements(fieldModule, coordinates, elementIdentifier, part_structure, [])
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
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, part_structure.px [n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, part_structure.pd1[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, part_structure.pd2[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, part_structure.pd3[n3][n2][n1])
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
                    nids = [ part_structure.nodeId[e3][e2][e1], part_structure.nodeId[e3][e2 + 1][e1],
                             part_structure.nodeId[e3+1][e2][e1], part_structure.nodeId[e3+1][e2 + 1][e1],
                             part_structure.nodeId[e3][e2][e1 + 1], part_structure.nodeId[e3][e2 + 1][e1 + 1],
                             part_structure.nodeId[e3+1][e2][e1 + 1], part_structure.nodeId[e3+1][e2 + 1][e1 + 1] ]

                    if isinstance(part_structure, JoiningBox):
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [-1.0]
                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                        remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                        if e2 == e2a or e2 == e2z:
                            if e2 == e2a:
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                            elif e2 == e2z:
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                                remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                        else:
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                    elif isinstance(part_structure, JoiningBoxRight):
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [-1.0]
                        if e2 == e2a:
                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D2_DS1DS2, [(Node.VALUE_LABEL_D_DS2, [])])
                        elif e2 == e2z:
                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                            remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D2_DS1DS2, [(Node.VALUE_LABEL_D_DS2, [])])

                        else:
                            remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                            remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D2_DS1DS2, [(Node.VALUE_LABEL_D_DS2, [])])

                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])

                    elif isinstance(part_structure, JoiningBox3):
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [-1.0]
                        if e2 == e2a:
                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1])])

                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                        elif e2 == e2z:
                            remapEftNodeValueLabel(eft1, [1, 3, 5, 6], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 3, 5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                        else:
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS3,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])

                    elif isinstance(part_structure, JoiningBox4):
                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        scalefactors = [-1.0]
                        if e2 == e2a:
                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])

                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D2_DS1DS2, [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                        elif e2 == e2z:
                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                            remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D2_DS1DS2, [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, []), (Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                        else:
                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D2_DS1DS2, [])])
                            remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                            remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D2_DS1DS2, [(Node.VALUE_LABEL_D_DS2, [])])

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
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                elif e1 == e1y:

                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, []), (Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1]), (Node.VALUE_LABEL_D_DS1, [])])

                                else:
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1])])
                                    if (e1 == e1b) or (e1 == e1y):
                                        # map bottom triple point element
                                        if e1 == e1b:
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                            remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
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
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                    elif e1 == e1y:
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
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

                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2,[(Node.VALUE_LABEL_D_DS3, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2,[(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    elif e2 == e2y:
                                        nids[1] = part_structure.nodeId[e3][e2z+1][e1b]
                                        nids[3] = part_structure.nodeId[e3+1][e2z+1][e1b]
                                        remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                    remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                elif e1 == e1z:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    if e2 == e2b:
                                        nids[4] = part_structure.nodeId[e3][e2a][e1z]
                                        nids[6] = part_structure.nodeId[e3+1][e2a][e1z]
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [ -1.0 ]
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                                    elif e2 == e2y:
                                        nids[5] = part_structure.nodeId[e3][e2z+1][e1z]
                                        nids[7] = part_structure.nodeId[e3+1][e2z+1][e1z]
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [-1.0]
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                            elif e1 == e1b:
                                if e2 == e2b:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                if e2 == e2y:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                            elif e1 == e1y:
                                if e2 == e2b:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]

                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                                if e2 == e2y:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])

                        else:
                            if e1 < e1a:
                                nids = [ part_structure.nodeId[e3][e2 + 1][e1 + 1], part_structure.nodeId[e3][e2][e1 + 1], part_structure.nodeId[e3+1][e2 + 1][e1 + 1], part_structure.nodeId[e3+1][e2][e1 + 1],
                                         part_structure.nodeId[e3][e2 + 1][e1], part_structure.nodeId[e3][e2][e1], part_structure.nodeId[e3+1][e2 + 1][e1], part_structure.nodeId[e3+1][e2][e1]]
                            elif e1 == e1a:
                                # map left column elements
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [ -1.0 ]
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS3, [1])])

                    else:
                        if (e2 < e2b) or (e2 > e2y):
                            if (e1 < e1b) or (e1 > e1y):
                                continue  # no element due to triple point closure
                            if (e2 < e2a) or (e2 > e2z):
                                if e2 < e2a:
                                    nids = [part_structure.nodeId[e3][e2+1][e1], part_structure.nodeId[e3][e2+1][e1+1], part_structure.nodeId[e3+1][e2+1][e1], part_structure.nodeId[e3+1][e2+1][e1+1],
                                            part_structure.nodeId[e3][e2][e1], part_structure.nodeId[e3][e2][e1+1], part_structure.nodeId[e3+1][e2][e1],  part_structure.nodeId[e3+1][e2][e1+1]]
                                elif e2 > e2z:
                                    nids = [part_structure.nodeId[e3][e2][e1+1], part_structure.nodeId[e3][e2][e1], part_structure.nodeId[e3+1][e2][e1+1], part_structure.nodeId[e3+1][e2][e1],
                                            part_structure.nodeId[e3][e2+1][e1+1], part_structure.nodeId[e3][e2+1][e1], part_structure.nodeId[e3+1][e2+1][e1+1], part_structure.nodeId[e3+1][e2+1][e1]]
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
                                                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                               [(Node.VALUE_LABEL_D_DS3, [1])])
                                        remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                    if (e1 == e1b) or (e1 == e1y):
                                        # map bottom triple point element
                                        if e1 == e1b:
                                            # if e3 != 2:
                                            #     remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                            if part_structure._shoulder_left:
                                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                                                remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                                                remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                                            else:
                                                if part_structure._shoulder:
                                                    if e3 == 0:
                                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                                                        remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                                        remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                                    else:
                                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                                else:
                                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS1,
                                                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            if part_structure._shoulder_left:
                                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                                                remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [1])])
                                                remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                                                remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
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
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [1])])
                                    remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS3, [])])
                                    if (e1 == e1b) or (e1 == e1y):
                                        # map top triple point element
                                        if e1 == e1b:
                                            remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                            if part_structure._shoulder:
                                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS1, []),(Node.VALUE_LABEL_D_DS3, [])])

                        elif (e2 == e2b) or (e2 == e2y):
                            if (e1 <= e1a) or (e1 >= e1z):
                                if e1 < e1a:
                                    e2r = e1
                                    if e2 == e2b:
                                        nids = [part_structure.nodeId[e3][e2c][e1+1], part_structure.nodeId[e3][e2r+1][e1b], part_structure.nodeId[e3+1][e2c][e1+1], part_structure.nodeId[e3+1][e2r+1][e1b],
                                                part_structure.nodeId[e3][e2c][e1], part_structure.nodeId[e3][e2r][e1b], part_structure.nodeId[e3+1][e2c][e1], part_structure.nodeId[e3+1][e2r][e1b]]
                                    if e2 == e2y:
                                        e2r = 2*part_structure._elementsCount[1] - e1-1
                                        nids = [part_structure.nodeId[e3][e2r][e1b], part_structure.nodeId[e3][e2y][e1+1], part_structure.nodeId[e3+1][e2r][e1b], part_structure.nodeId[e3+1][e2y][e1+1],
                                                part_structure.nodeId[e3][e2r+1][e1b], part_structure.nodeId[e3][e2y][e1], part_structure.nodeId[e3+1][e2r+1][e1b], part_structure.nodeId[e3+1][e2y][e1]]
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
                                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                            remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                            remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            tripleN = [5, 7]
                                            remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    elif e2 == e2y:
                                        nids[1] = part_structure.nodeId[e3][e2z+1][e1b]
                                        nids[3] = part_structure.nodeId[e3+1][e2z+1][e1b]
                                        tripleN = [6, 8]
                                        remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                    remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                                elif e1 == e1z:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    if e2 == e2b:
                                        setEftScaleFactorIds(eft1, [1], [])
                                        nids[4] = part_structure.nodeId[e3][e2a][e1z]
                                        nids[6] = part_structure.nodeId[e3+1][e2a][e1z]
                                        setEftScaleFactorIds(eft1, [1], [])
                                        scalefactors = [ -1.0 ]
                                        if part_structure._shoulder_left:
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                            remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                            if e3 == 0:
                                                remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                            remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                        else:
                                            remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                    elif e2 == e2y:
                                        nids[5] = part_structure.nodeId[e3][e2z+1][e1z]
                                        nids[7] = part_structure.nodeId[e3+1][e2z+1][e1z]
                                        remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e1 > e1z:
                                    e2r = part_structure._elementsCount[0] - e1
                                    if e2 == e2b:
                                        nids = [part_structure.nodeId[e3][e2r][e1z], part_structure.nodeId[e3][e2c][e1], part_structure.nodeId[e3+1][e2r][e1z], part_structure.nodeId[e3+1][e2c][e1],
                                                part_structure.nodeId[e3][e2r-1][e1z], part_structure.nodeId[e3][e2c][e1+1], part_structure.nodeId[e3+1][e2r-1][e1z], part_structure.nodeId[e3+1][e2c][e1+1]]
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
                                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                                    if part_structure._shoulder:
                                        remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])
                                if e2 == e2y:
                                    if part_structure._shoulder:
                                        eft1 = tricubichermite.createEftNoCrossDerivatives()
                                        remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [])])

                            elif e1 == e1y:
                                if e2 == e2b:
                                    eft1 = tricubichermite.createEftNoCrossDerivatives()
                                    if part_structure._shoulder_left:
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS1, [])])

                        else:
                            if e1 < e1a:
                                nids = [ part_structure.nodeId[e3][e2 + 1][e1 + 1], part_structure.nodeId[e3][e2][e1 + 1], part_structure.nodeId[e3+1][e2 + 1][e1 + 1], part_structure.nodeId[e3+1][e2][e1 + 1],
                                         part_structure.nodeId[e3][e2 + 1][e1], part_structure.nodeId[e3][e2][e1], part_structure.nodeId[e3+1][e2 + 1][e1], part_structure.nodeId[e3+1][e2][e1]]
                            elif e1 == e1a:
                                # map left column elements
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [ -1.0 ]
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS3, [1])])

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
                    #print('create element shield', elementIdentifier, result2, result3, nids)
                    part_structure.elementId[e3][e2][e1] = elementIdentifier
                    elementIdentifier += 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

        return elementIdentifier


class BaseLeg:
    """
    Base case for creating a child
    """
    def __init__(self, elementsCount, nodeparams):
        """

        :param fieldmodule:
        :param coordinates:
        :param mesh:
        :param nodes:
        """
        elementsCount = [4, 4, 2]
        self._elementsCount = [elementsCount[0]//2, elementsCount[1]//2, elementsCount[2]]
        self.nodeparams = nodeparams
        self._shoulder = False
        self._shoulder_left = False
        self._neck = False

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])] for c in range(elementsCount[2])]
        self.generateBaseLeg(nodeparams)

        n, d = geometry.get_plane_normal_vector_and_distance(nodeparams[0][0], nodeparams[0][1], nodeparams[1][0])
        plane = [n[0], n[1], n[2], d]
        mirror = Mirror(plane)
        for n2 in range(elementsCount[1]//2+1, elementsCount[1]+1):
            for n3 in range(elementsCount[2]+1):
                for n1 in range(elementsCount[0]//2+1):
                    n3q = n3
                    n2q = elementsCount[1] - n2
                    n1q = n1
                    if self.px[n3q][n2q][n1q]:
                        self.px[n3][n2][n1] = mirror.mirrorImageOfPoint(self.px[n3q][n2q][n1q])
                        self.pd1[n3][n2][n1] = mirror.reverseMirrorVector(self.pd1[n3q][n2q][n1q])
                        self.pd2[n3][n2][n1] = mirror.mirrorVector(self.pd2[n3q][n2q][n1q])
                        self.pd3[n3][n2][n1] = mirror.mirrorVector(self.pd3[n3q][n2q][n1q])

        n, d = geometry.get_plane_normal_vector_and_distance(nodeparams[0][0], nodeparams[0][2], nodeparams[1][0])
        plane = [n[0], n[1], n[2], d]
        mirror = Mirror(plane)
        for n2 in range(elementsCount[1]+1):
            for n3 in range(elementsCount[2]+1):
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
        self._elementsCount = elementsCount

    def generateBaseLeg(self, nodeparams):
        """
        Generate base leg that is a cylinder generated from cylinder ends.
        :return:
        """
        bottomnodeparams = nodeparams[0]
        topnodeparams = nodeparams[1]
        self.genetateBottomSurface(bottomnodeparams)
        self.generateTopSurface(topnodeparams)
        self.generateMiddleLevels()
        self.smoothd2()

    def genetateBottomSurface(self, bottomnodeparams):
        """
        Use major and minor curves to generate the ellipse
        :return:
        """
        centre, xec1, xec2, d1c1, d1c2 = bottomnodeparams[0], bottomnodeparams[1], bottomnodeparams[2], bottomnodeparams[3], bottomnodeparams[4]
        txc1, td1c1 = self.generate1DPath([centre, xec1], d1c1, self._elementsCount[0])
        txc2, td1c2 = self.generate1DPath([centre, xec2], d1c2, self._elementsCount[1])
        ellipse = self.generateSurfaceUsingTwoCurves(centre, txc1, td1c1, txc2, td1c2)
        self.copyEllipseNodesToBifurcation(ellipse, 0)

    def generateTopSurface(self, topnodeparams):
        """

        :return:
        """
        centre, xec1, xec2, d1c1, d1c2 = topnodeparams[0], topnodeparams[1], topnodeparams[2], topnodeparams[3], topnodeparams[4]
        txc1, td1c1 = self.generate1DPath([centre, xec1], d1c1, self._elementsCount[0])
        txc2, td1c2 = self.generate1DPath([centre, xec2], d1c2, self._elementsCount[1])
        ellipse = self.generateSurfaceUsingTwoCurves(centre, txc1, td1c1, txc2, td1c2)
        self.copyEllipseNodesToBifurcation(ellipse, self._elementsCount[2])

    def generateMiddleLevels(self):
        """

        :return:
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3
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
        majorAxis = [c*self._elementsCount[0] for c in td1c2[0]]
        minorAxis = [-c*self._elementsCount[1] for c in td1c1[0]]
        elementsCountAcrossShell = 0
        elementsCountAcrossTransition = 1
        shellProportion = 1.0
        coreMajorRadius = 1.0
        coreMinorRadius = 1.0
        elementsCountAround = self._elementsCount[0] + self._elementsCount[1] - 2
        ellipse = Ellipse2D(centre, majorAxis, minorAxis, 2*self._elementsCount[0], 2*self._elementsCount[1],
                            elementsCountAcrossShell, elementsCountAcrossTransition,
                            shellProportion, coreMajorRadius, coreMinorRadius, ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

        shield = ellipse.getShield()
        ellipse.generateBase1DMesh(0)

        ellipse.px[self._elementsCount[1]][0] = txc1[-1]
        ellipse.pd3[self._elementsCount[1]][0] = td1c1[-1]
        ellipse.px[0][self._elementsCount[0]] = txc2[-1]
        ellipse.pd3[0][self._elementsCount[0]] = td1c2[-1]

        nx = []
        nd1 = []
        for n in range(elementsCountAround + 1):
            n1, n2 = shield.convertRimIndex(n)
            xi = n/elementsCountAround
            # ellipse.px[n2][n1][2] = (1 - xi) * ellipse.px[self._elementsCount[1]][0][2] + xi * \
            #                         ellipse.px[0][self._elementsCount[0]][2]
            x_p = vector.addVectors([vector.scaleVector(ellipse.px[self._elementsCount[1]][0], 1 - xi),
                                     vector.scaleVector(ellipse.px[0][self._elementsCount[0]], xi)], [1, 1])
            delta = [-ellipse.px[0][self._elementsCount[0]][c] + ellipse.px[self._elementsCount[1]][0][c] for c in range(3)]
            normal = vector.normalise(vector.crossproduct3(td1c1[0], td1c2[0]))
            delta_p = vector.addVectors([x_p, ellipse.px[self._elementsCount[1]][self._elementsCount[0]]], [1, -1])
            delta_p_normalmag = vector.dotproduct(delta_p, normal)
            delta_p_normal = vector.scaleVector(normal, delta_p_normalmag)
            delnormag = vector.dotproduct(delta, normal)
            delnor = vector.scaleVector(normal, delnormag)
            if 0<n<elementsCountAround:
                # ellipse.px[n2][n1] = [ellipse.px[n2][n1][c] +(1-xi) * delnor[c] for c in range(3)]
                ellipse.px[n2][n1] = vector.addVectors([ellipse.px[n2][n1], delta_p_normal], [1, 1])
            nx.append(ellipse.px[n2][n1])
            nd1.append(ellipse.pd1[n2][n1])
        td1 = smoothCubicHermiteDerivativesLine(nx, nd1)
        for n in range(1, elementsCountAround):
            n1, n2 = shield.convertRimIndex(n)
            ellipse.pd1[n2][n1] = td1[n]


        ellipse.setRimNodes()
        rscx, rscd1, rscd2, rscd3 = ellipse.createMirrorCurve()
        ellipse.createRegularRowCurves(ellipse.px[:self._elementsCount[1]+1][self._elementsCount[0]], rscd1, rscd3)
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

    def copyEllipseNodesToBifurcation(self, ellipse, n3):
        """
        Copy ellipse nodes to bifurcation
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


class JoiningBox:
    def __init__(self, elementsCount, height):
        self._elementsCount = elementsCount
        self._joining_box = True

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])] for c in range(elementsCount[2])]

        self.px[1][0][1] = [0.0, 1.0, height]
        self.px[1][1][1] = [0.0, 0.5, height]
        self.px[1][2][1] = [0.0, 0.0, height]
        self.px[1][3][1] = [0.0, -0.5, height]
        self.px[1][4][1] = [0.0, -1.0, height]
        for n2 in range(elementsCount[1]+1):
            self.pd1[1][n2][1] = [0.0, -0.5, 0.0]
            self.pd2[1][n2][1] = [0.5, 0.0, 0.0]
            self.pd3[1][n2][1] = [0.0, 0.0, 0.7]


class CylinderConnectingToBox:
    """

    """
    def __init__(self, cylinder_part, idxs, sign):
        elementsCount = [cylinder_part._elementsCount[0], cylinder_part._elementsCount[1], 1]
        self._elementsCount = elementsCount

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])] for c in range(elementsCount[2])]

        for n2 in range(elementsCount[1] + 1):
            for n1 in range(elementsCount[0]//2 + 1, elementsCount[0] + 1):
                if cylinder_part.px[idxs[0]][n2][n1]:
                    self.px[idxs[0]][n2][n1] = vector.addVectors([cylinder_part.px[idxs[1]][n2][n1],
                                                            cylinder_part.pd2[idxs[1]][n2][n1]], [1, sign])
                    self.pd1[idxs[0]][n2][n1] = cylinder_part.pd1[idxs[1]][n2][n1]
                    self.pd2[idxs[0]][n2][n1] = cylinder_part.pd2[idxs[1]][n2][n1]
                    self.pd3[idxs[0]][n2][n1] = cylinder_part.pd3[idxs[1]][n2][n1]


class JoiningTorso:
    def __init__(self, elementsCount):
        self._elementsCount = elementsCount
        self._joining_torso = True

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])] for c in range(elementsCount[2])]


class BifurcationMeshCrotch:
    """
    Bifurction mesh generator.
    """

    def __init__(self, fieldModule, coordinates, region):
        """
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        """
        # generate the mesh
        elementsCount = [2, 2, 5]
        self._elementsCount = elementsCount
        self._region = region

        self.createBifurcationMesh3d(fieldModule, coordinates)

    def createBifurcationMesh3d(self, fieldmodule, coordinates):
        """
        Create a bifurcation.
        :param fieldModule: Zinc fieldModule to create elements in.
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


class CrotchEllipses:
    def __init__(self):
        """

        """


class JoiningBoxRight:
    def __init__(self, elementsCount, bottom_part, shoulder_part_left, joining_box):
        self._elementsCount = elementsCount
        self._joining_box = True

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])] for c in range(elementsCount[2])]

        for n2 in range(elementsCount[1]+1):
            self.nodeId[1][n2][0] = joining_box.nodeId[1][n2][1]
            self.px[1][n2][0] = joining_box.px[1][n2][1]
            self.pd1[1][n2][0] = joining_box.pd1[1][n2][1]
            self.pd2[1][n2][0] = joining_box.pd2[1][n2][1]
            self.pd3[1][n2][0] = joining_box.pd3[1][n2][1]

            self.nodeId[0][n2][0] = joining_box.nodeId[0][n2][1]
            self.px[0][n2][0] = joining_box.px[0][n2][1]
            self.pd1[0][n2][0] = joining_box.pd1[0][n2][1]
            self.pd2[0][n2][0] = joining_box.pd2[0][n2][1]
            self.pd3[0][n2][0] = joining_box.pd3[0][n2][1]

            self.nodeId[0][n2][1] = bottom_part.nodeId[2][n2][3]
            self.px[0][n2][1] = bottom_part.px[2][n2][3]
            self.pd1[0][n2][1] = bottom_part.pd1[2][n2][3]
            self.pd2[0][n2][1] = bottom_part.pd2[2][n2][3]
            self.pd3[0][n2][1] = bottom_part.pd3[2][n2][3]

            n2s = 2
            n1s = n2
            self.nodeId[1][n2][1] = shoulder_part_left.nodeId[0][n2s][n1s]
            self.px[1][n2][1] = shoulder_part_left.px[0][n2s][n1s]
            self.pd1[1][n2][1] = shoulder_part_left.pd1[0][n2s][n1s]
            self.pd2[1][n2][1] = shoulder_part_left.pd2[0][n2s][n1s]
            self.pd3[1][n2][1] = shoulder_part_left.pd3[0][n2s][n1s]


class JoiningBox3:
    def __init__(self, elementsCount, joining_box, neck_part):
        self._elementsCount = elementsCount
        self._joining_box = True

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])] for c in range(elementsCount[2])]

        for n2 in range(elementsCount[1]+1):
            self.nodeId[0][n2][0] = joining_box.nodeId[1][n2][0]
            self.px[0][n2][0] = joining_box.px[1][n2][0]
            self.pd1[0][n2][0] = joining_box.pd1[1][n2][0]
            self.pd2[0][n2][0] = joining_box.pd2[1][n2][0]
            self.pd3[0][n2][0] = joining_box.pd3[1][n2][0]
            self.nodeId[0][n2][1] = joining_box.nodeId[1][n2][1]
            self.px[0][n2][1] = joining_box.px[1][n2][1]
            self.pd1[0][n2][1] = joining_box.pd1[1][n2][1]
            self.pd2[0][n2][1] = joining_box.pd2[1][n2][1]
            self.pd3[0][n2][1] = joining_box.pd3[1][n2][1]

            self.nodeId[1][n2][0] = neck_part.nodeId[0][n2][1]
            self.px[1][n2][0] = neck_part.px[0][n2][1]
            self.pd1[1][n2][0] = neck_part.pd1[0][n2][1]
            self.pd2[1][n2][0] = neck_part.pd2[0][n2][1]
            self.pd3[1][n2][0] = neck_part.pd3[0][n2][1]

            self.nodeId[1][n2][1] = neck_part.nodeId[0][n2][2]
            self.px[1][n2][1] = neck_part.px[0][n2][2]
            self.pd1[1][n2][1] = neck_part.pd1[0][n2][2]
            self.pd2[1][n2][1] = neck_part.pd2[0][n2][2]
            self.pd3[1][n2][1] = neck_part.pd3[0][n2][2]


class JoiningBox4:
    def __init__(self, elementsCount, joining_box_3, neck_part, joining_box_right):
        self._elementsCount = elementsCount
        self._joining_box = True
        self._neck_part = False

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])] for c in range(elementsCount[2])]

        for n2 in range(elementsCount[1]+1):
            self.nodeId[0][n2][0] = joining_box_3.nodeId[0][n2][1]
            self.px[0][n2][0] = joining_box_3.px[0][n2][1]
            self.pd1[0][n2][0] = joining_box_3.pd1[0][n2][1]
            self.pd2[0][n2][0] = joining_box_3.pd2[0][n2][1]
            self.pd3[0][n2][0] = joining_box_3.pd3[0][n2][1]

            self.nodeId[1][n2][0] = joining_box_3.nodeId[1][n2][1]
            self.px[1][n2][0] = joining_box_3.px[1][n2][1]
            self.pd1[1][n2][0] = joining_box_3.pd1[1][n2][1]
            self.pd2[1][n2][0] = joining_box_3.pd2[1][n2][1]
            self.pd3[1][n2][0] = joining_box_3.pd3[1][n2][1]

            self.nodeId[0][n2][1] = joining_box_right.nodeId[1][n2][1]
            self.px[0][n2][1] = joining_box_right.px[1][n2][1]
            self.pd1[0][n2][1] = joining_box_right.pd1[1][n2][1]
            self.pd2[0][n2][1] = joining_box_right.pd2[1][n2][1]
            self.pd3[0][n2][1] = joining_box_right.pd3[1][n2][1]

            self.nodeId[1][n2][1] = neck_part.nodeId[0][n2][3]
            self.px[1][n2][1] = neck_part.px[0][n2][3]
            self.pd1[1][n2][1] = neck_part.pd1[0][n2][3]
            self.pd2[1][n2][1] = neck_part.pd2[0][n2][3]
            self.pd3[1][n2][1] = neck_part.pd3[0][n2][3]





