"""
Utility functions for generating a solid spheroid.
"""

from enum import Enum
from scaffoldmaker.utils import vector
import math
from opencmiss.zinc.field import Field
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from scaffoldmaker.utils.shieldmesh import ShieldMesh3D, ShieldShape3D
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite, \
    smoothCubicHermiteDerivativesLine, interpolateSampleLinear, interpolateCubicHermite
from scaffoldmaker.utils.mirror import Mirror
from scaffoldmaker.utils.cylindermesh import Ellipse2D, EllipseShape


class SphereShape(Enum):
    SPHERE_SHAPE_FULL = 1
    SPHERE_SHAPE_HALF_NNP = 2    # NNP is a3>=3
    SPHERESHIELD_SHAPE_OCTANT_PPP = 3


class SphereMesh:
    """
    Sphere mesh generator.
    """

    def __init__(self, fieldModule, coordinates, centre, axes, elementsCountAcross,
                 elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                 sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False):
        """
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :param centre, axes: centre and axes of the sphere.
        :param elementsCountAcross: [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3] Total number of elements
         across the sphere axes.
        :param elementsCountAcrossShell, elementsCountAcrossTransition: Total number of elements across each axis
         consists of regular elements in the middle cube, transition elements from cube to a sphere (core boundary)
          and shell elements around it. Shell nodes and derivatives are similar to the core boundary and don't need
           remapping. The topology of the shield structure is extended to 3D with a quadruple points.
        :param sphereShape: A value from enum sphereMode specifying. Octant_PPP for example, is the octant in axis1>=0
        axis2>=0 and axis3>=0
        """

        self._axes = axes
        self._radius = [vector.magnitude(axis) for axis in axes]
        self._coreRadius = []
        self._shield = None
        self._elementsCount = elementsCountAcross
        self._elementsCountAcrossShell = elementsCountAcrossShell
        self._elementsCountAcrossTransition = elementsCountAcrossTransition
        self._elementsCountAcrossRim = self._elementsCountAcrossShell + self._elementsCountAcrossTransition - 1
        self._shellProportion = shellProportion
        self._elementsCountAround12 = 2 * (self._elementsCount[0] + self._elementsCount[1] -
                                         4*(self._elementsCountAcrossRim + 1))
        self._startNodeIdentifier = 1
        self._startElementIdentifier = 1
        self._endNodeIdentifier = 1
        self._endElementIdentifier = 1
        self._sphereShape = sphereShape

        self._useCrossDerivatives = useCrossDerivatives

        self._centre = centre

        for i in range(3):
            elementsAxis = elementsCountAcross[i] - elementsCountAcrossShell * (1 - shellProportion)
            self._coreRadius.append(
                (1 - shellProportion * elementsCountAcrossShell / elementsAxis) * self._radius[i])

        # generate the mesh
        self.createSphereMesh3d(fieldModule, coordinates)

    def createSphereMesh3d(self, fieldModule, coordinates):
        """
        Create a sphere mesh based on the shield topology.
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :return: Final values of nextNodeIdentifier, nextElementIdentifier.
        """
        for i in range(3):
            assert (self._elementsCount[i] > 1), 'createSphereMesh3d:  Invalid number of elements'
            # assert (self._elementsCount[i] % 2 == 0), 'createSphereMesh3d: number of across elements' \
            #                                               ' is not an even number'

        assert (self._sphereShape in [self._sphereShape.SPHERE_SHAPE_FULL,
                                        self._sphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP]), \
            'createSphereMesh3d: Invalid sphere mode.'

        nodes = fieldModule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldModule.findMeshByDimension(3)

        elementsCountRim = self._elementsCountAcrossRim

        shieldMode = ShieldShape3D.SHIELD_SHAPE_FULL if self._sphereShape is self._sphereShape.SPHERE_SHAPE_FULL \
            else ShieldShape3D.SHIELD_SHAPE_OCTANT_PPP

        self._shield3D = ShieldMesh3D(self._elementsCount, elementsCountRim, shieldMode=shieldMode)

        self.calculateBoundaryEllipses()
        self.generateNodes(nodes, fieldModule, coordinates)
        self.generateElements(mesh, fieldModule, coordinates)

    def calculateBoundaryEllipses(self):
        """

        :return:
        """

        centre = self._centre
        elementsCountAcrossShell = self._elementsCountAcrossShell
        elementsCountAcrossTransition = self._elementsCountAcrossTransition
        shellProportion = self._shellProportion

        ellipseAxes = [[self._axes[0], self._axes[1], self._axes[2]],
                       [[-c for c in self._axes[2]], self._axes[1], self._axes[0]],
                       [self._axes[0], [-c for c in self._axes[2]], self._axes[1]]]

        coreRadius = [(self._coreRadius[0], self._coreRadius[1]),
                      (self._coreRadius[2], self._coreRadius[1]),
                      (self._coreRadius[0], self._coreRadius[2])]

        elementsCount = [(2 * self._elementsCount[0], 2 * self._elementsCount[1]),
                         (2 * self._elementsCount[2], 2 * self._elementsCount[1]),
                         (2 * self._elementsCount[0], 2 * self._elementsCount[2])]

        for i in range(3):
            majorAxis = ellipseAxes[i][0]
            minorAxis = ellipseAxes[i][1]
            alongAxis = ellipseAxes[i][2]
            elementsCountAcrossMajor = elementsCount[i][0]
            elementsCountAcrossMinor = elementsCount[i][1]
            coreMajorRadius = coreRadius[i][0]
            coreMinorRadius = coreRadius[i][1]
            ellipse = Ellipse2D(centre, majorAxis, minorAxis,
                         elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                         elementsCountAcrossTransition, shellProportion, coreMajorRadius, coreMinorRadius,
                         ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

            self.copyEllipsesNodesToShieldNodes(ellipse, i)
        self.createAdditionalPointsForIncreasingElementsCount()

    def copyEllipsesNodesToShieldNodes(self, ellipse, ellipsenumber):
        """
        Copy coordinates and derivatives of ellipse to shield.
        :param n3: the index number of ellipse along the central path.
        """

        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        shield = ellipse.getShield()

        # Modify the shield12 to get only the quarter that you want. make others None so you can generate the nodes.
        if ellipsenumber == 0:
            for n2 in range(self._elementsCount[0] + 1):  # TODO modify this to number of elements.
                for n1 in range(self._elementsCount[1] + 1):
                    # only first quadrant is needed  TODO we need to modify this to number of elements half of the elments across each direction.
                    # btd2[0][n2][n1] = [c for c in alongAxis]
                    if n2 == 0 and n1 == self._elementsCount[1] - 1:
                        n1s = n1 + 1
                    else:
                        n1s = n1
                    n1e = n1 + self._elementsCount[1]
                    if shield.px[0][n2][n1e]:
                        btx[0][n2][n1s] = shield.px[0][n2][n1e]
                        btd1[0][n2][n1s] = shield.pd1[0][n2][n1e]
                        btd2[0][n2][n1s] = shield.pd2[0][n2][n1e]
                        btd3[0][n2][n1s] = shield.pd3[0][n2][n1e]
        elif ellipsenumber == 1:
            shield.remap_derivatives([2, 3], circleMapping=None)
            n2s = self._elementsCount[0]
            for n3 in range(self._elementsCount[2] + 1):
                for n1 in range(self._elementsCount[1] + 1):
                    # btd2[n3][n2s][n1] = [c for c in alongAxis]
                    if n3 == self._elementsCount[2] and n1 == self._elementsCount[1] - 1:
                        n1s = n1 + 1
                    else:
                        n1s = n1
                    n2e = n3 + self._elementsCount[2]
                    n1e = n1 + self._elementsCount[1]
                    if n3 == 0:
                        if n1 == self._elementsCount[1]:
                            btd2[n3][n2s][n1s] = [c for c in shield.pd1[0][n2e][n1e]]
                        else:
                            btd2[n3][n2s][n1s] = [c for c in shield.pd2[0][n2e][n1e]]
                    else:
                        if shield.px[0][n2e][n1 + self._elementsCount[1]]:
                            btx[n3][n2s][n1s] = [c for c in shield.px[0][n2e][n1e]]
                            btd1[n3][n2s][n1s] = [c for c in shield.pd1[0][n2e][n1e]]
                            btd2[n3][n2s][n1s] = [c for c in shield.pd2[0][n2e][n1e]]
                            btd3[n3][n2s][n1s] = [c for c in shield.pd3[0][n2e][n1e]]
        elif ellipsenumber == 2:
            shield.remap_derivatives([1, -2], circleMapping=None)
            n1s = 0
            for n3 in range(self._elementsCount[2] + 1):
                for n2 in range(self._elementsCount[0] + 1):

                    if n2 == 0 and n3 == self._elementsCount[2] - 1:
                        n3s = self._elementsCount[2]
                    else:
                        n3s = n3
                    n1e = self._elementsCount[2] - n3
                    if n3 == 0:
                        if n2 == 0:
                            btd2[n3][n2][n1s] = [-c for c in shield.pd1[0][n2][n1e]]
                        elif 0 < n2 < self._elementsCount[0]:
                            btd2[n3][n2][n1s] = [c for c in shield.pd2[0][n2][n1e]]
                    else:
                        if n2 < self._elementsCount[0]:
                            # btd2[n3][n2][n1s] = [c for c in alongAxis]
                            if shield.px[0][n2][n1e]:
                                btx[n3s][n2][n1s] = [c for c in shield.px[0][n2][n1e]]
                                btd1[n3s][n2][n1s] = [c for c in shield.pd1[0][n2][n1e]]
                                btd2[n3s][n2][n1s] = [c for c in shield.pd2[0][n2][n1e]]
                                btd3[n3s][n2][n1s] = [c for c in shield.pd3[0][n2][n1e]]
                        else:
                            if n3 == self._elementsCount[2]:
                                btd2[n3][n2][n1s] = [c for c in shield.pd1[0][n2][n1e]]
                            else:
                                btd1[n3][n2][n1s] = [c for c in shield.pd1[0][n2][n1e]]

    def createAdditionalPointsForIncreasingElementsCount(self):
        """

        :return:
        """

        self.calculate_surface_quadruple_point()
        self.sample_triple_curves_on_sphere()
        self.sample_regular_curves_on_sphere()
        self.create_interior_nodes()
        self.smooth_derivatives_to_surface()

    def calculate_surface_quadruple_point(self):
        """
        Calculate coordinates and derivatives of the quadruple point on the surface, where 3 hex elements merge.
        :return:
        """

        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n1z = self._elementsCount[1]
        n1y = n1z - 1
        n2z = self._elementsCount[0]
        n3z = self._elementsCount[2]
        n3y = n3z - 1

        radius = self._radius[0]  # TODO need to be changed for spheroid

        elementsAroundEllipse12 = self._elementsCount[0] + self._elementsCount[1] - 2
        radiansAroundEllipse12 = math.pi / 2
        radiansPerElementAroundEllipse12 = radiansAroundEllipse12 / elementsAroundEllipse12
        elementsAroundEllipse13 = self._elementsCount[0] + self._elementsCount[2] - 2
        radiansAroundEllipse13 = math.pi / 2
        radiansPerElementAroundEllipse13 = radiansAroundEllipse13 / elementsAroundEllipse13

        theta_2 = n3y * radiansPerElementAroundEllipse13
        theta_3 = n1y * radiansPerElementAroundEllipse12
        phi_3 = calculate_azimuth(theta_3, theta_2)
        # We assume it is a sphere not a spheroid for now. TODO Use the relations for spheroid instead
        # ratio = -0.1 * (min(self._elementsCount) - 2) + 1 if self._elementsCount[0] <= 2 else 0.2
        ratio = 1
        # local_x = intersection_of_two_great_circles_on_sphere(btx[0][0][n1y-1], btx[n3z][n2z][n1z], btx[0][2][n1z], btx[n3z][0][0])
        local_x = spherical_to_cartesian(radius, theta_3, ratio * phi_3 + (1-ratio)*math.pi/2)

        x = local_to_global_coordinates(local_x, self._axes, self._centre)

        a1, a2, a3 = local_orthogonal_unit_vectors(x, self._axes[2])
        n3r, n2r, n1r = self.get_triple_curves_end_node_parameters(1, index_output=True)
        btx[n3r][n2r][n1r] = x
        btd1[n3r][n2r][n1r] = a1  # initialise
        btd2[n3r][n2r][n1r] = a2  # initialise
        btd3[n3r][n2r][n1r] = a3

    def sample_triple_curves_on_sphere(self):
        """
        Sample points on the triple curves of the 'quadruple point' on the sphere surface.
        :return:
        """
        n1z = self._elementsCount[1]
        n2z = self._elementsCount[0]
        n3z = self._elementsCount[2]

        # sample on curve 1 of the triple curves and smooth the end derivatives.
        n3r1, n2r1, n1r1 = self.get_triple_curves_end_node_parameters(1, index_output=True)
        n3r2, n2r2, n1r2 = self.get_triple_curves_end_node_parameters(1, cx=1, index_output=True)
        self.sample_curves_between_two_nodes_on_sphere([n3r1, n2r1, n1r1], [n3r2, n2r2, n1r2], self._elementsCount[0] - 1,
                                                       [1, None], [1], [-2])
        # curve 2
        n3r1, n2r1, n1r1 = self.get_triple_curves_end_node_parameters(1, cx=2, index_output=True)
        n3r2, n2r2, n1r2 = self.get_triple_curves_end_node_parameters(1, index_output=True)
        self.sample_curves_between_two_nodes_on_sphere([n3r1, n2r1, n1r1], [n3r2, n2r2, n1r2], self._elementsCount[1] - 1,
                                                       [2], [-2], [None, -2])
        # curve 3.
        n3r1, n2r1, n1r1 = self.get_triple_curves_end_node_parameters(1, cx=3, index_output=True)
        n3r2, n2r2, n1r2 = self.get_triple_curves_end_node_parameters(1, index_output=True)
        self.sample_curves_between_two_nodes_on_sphere([n3r1, n2r1, n1r1], [n3r2, n2r2, n1r2], self._elementsCount[2] - 1,
                                                       [2], [2], [None, None])

    def sample_regular_curves_on_sphere(self):
        """
        Create all other nodes on the sphere except the nodes on the triple curves and ellipses.
        :return:
        """
        n2a = 1
        n1z = self._elementsCount[1]
        n2z = self._elementsCount[0]
        n3z = self._elementsCount[2]

        # regular curves crossing curve 1
        for n2 in range(2, self._elementsCount[0]):
            # bottom right
            self.sample_curves_between_two_nodes_on_sphere([0, n2, n1z], [n3z, n2, n1z], self._elementsCount[2] - 1,
                                                           [2], [2], [None, None])
            # top
            self.sample_curves_between_two_nodes_on_sphere([n3z, n2, 0], [n3z, n2, n1z], self._elementsCount[1] - 1,
                                                           [2], [-2], [None, -2])

        # regular curves crossing curve 2
        for n1 in range(1, self._elementsCount[1] - 1):
            # bottom left. Top is done before.
            self.sample_curves_between_two_nodes_on_sphere([0, 0, n1], [n3z, 0, n1], self._elementsCount[2] - 1,
                                                           [2], [2], [None, None])

        # smooth regular curves crossing curve 1
        for n2 in range(n2a + 1, n2z):
            self.smooth_derivatives_regular_surface_curve(2, n2, [[2], [None, None]], [[-2], [-2]], [[None, -2], [-2]])

        # smooth regular curves crossing curve 2
        for n1 in range(1, n1z - 1):
            self.smooth_derivatives_regular_surface_curve(1, n1, [[2], [None, None]], [[2], [1]], [[None, 1], [-2]])

        # smooth regular curves crossing curve 3
        for n3 in range(1, self._elementsCount[2] - 1):
            self.smooth_derivatives_regular_surface_curve(3, n3, [[2], [None, None]], [[1], [1]], [[None, 1], [-2]])

    def create_interior_nodes(self):
        """

        :return:
        """

        self.calculate_interior_quadruple_point()
        self.sample_interior_curves()
        self.smooth_regular_interior_curves()

    def sample_curves_between_two_nodes_on_sphere(self, id1, id2, elementsOut, dStart, dbetween, dEnd):
        """
        samples curves on the sphere surface between two points given by their indexes.
        :param id1, id2: [n3,n2,n1] for the first and second points.
        :param dStart, dBetween, dEnd: Specifies the derivatives that are used for this curve at the beginning, end and
         in between. e.g. dStart=[2, -1, None] means d2 for the first node, -1 for the second node and skip the third one.
        :return:
        """
        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        # Find what index is constant
        if id1[0] != id2[0]:
            varying_index = 3
            elementsCount = self._elementsCount[2]
        elif id1[1] != id2[1]:
            varying_index = 2
            elementsCount = self._elementsCount[0]
        elif id1[2] != id2[2]:
            varying_index = 1
            elementsCount = self._elementsCount[1]

        else:
            raise ValueError("None of n1, n2, or n3 is constant. Only on the constant curves.")

        btd = {1: btd1, 2: btd2, 3: btd3}
        idi = {0: id1[0], 1: id1[1], 2: id1[2]}

        nx, nd1 = sample_curves_on_sphere(btx[id1[0]][id1[1]][id1[2]], btx[id2[0]][id2[1]][id2[2]], elementsOut)

        nit = 0
        for ni in range(elementsCount + 1):
            idi[3 - varying_index] = ni

            if ni < len(dStart):
                if dStart[ni]:
                    btd[dStart[ni]][idi[0]][idi[1]][idi[2]] = nd1[nit] if dStart[ni] > 0 else vector.scaleVector(
                        nd1[nit], -1)
                    nit += 1
            elif ni > elementsCount - len(dEnd):
                nie = ni - elementsCount + len(dEnd) - 1
                if dEnd[nie]:
                    btd[abs(dEnd[nie])][idi[0]][idi[1]][idi[2]] = nd1[nit] if dEnd[nie] > 0 else vector.scaleVector(
                        nd1[nit], -1)
                    nit += 1
            else:
                btx[idi[0]][idi[1]][idi[2]] = nx[nit]

                a1, a2, a3 = local_orthogonal_unit_vectors(nx[nit], self._axes[2])
                btd1[idi[0]][idi[1]][idi[2]] = a1  # initialise
                btd2[idi[0]][idi[1]][idi[2]] = a2  # initialise
                btd3[idi[0]][idi[1]][idi[2]] = a3  # initialise

                btd[abs(dbetween[0])][idi[0]][idi[1]][idi[2]] = nd1[nit] if dbetween[0] > 0 else vector.scaleVector(
                    nd1[nit], -1)
                nit += 1

    def smooth_derivatives_regular_surface_curve(self, constant_index, nc, dStart, dBetween, dEnd):
        """
        Smooth derivatives for each constant index curve. e.g. n2 = 3
        :param constant_index: Specifies n1, n2 or n3 is constant.
        :param nc: Index that is constant across the curve.
        :param dStart, dBetween, dEnd: See sample_curves_between_two_nodes_on_sphere. The difference here is
         the values are given for two curves that connect one end to the other end of the sphere surface.
        :return:
        """

        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n1z = self._elementsCount[1]
        n2z = self._elementsCount[0]
        n3z = self._elementsCount[2]

        if constant_index == 1:
            elementsCount = [self._elementsCount[2], self._elementsCount[0]]
        elif constant_index == 2:
            elementsCount = [self._elementsCount[1], self._elementsCount[2]]
        elif constant_index == 3:
            elementsCount = [self._elementsCount[1], self._elementsCount[0]]

        btd = {1: btd1, 2: btd2, 3: btd3}

        tx = []
        td = []
        for se in range(2):
            for ni in range(elementsCount[se] + 1):
                if constant_index == 1:
                    ids = [ni, 0, nc] if se == 0 else [n3z, ni, nc]
                elif constant_index == 2:
                    ids = [n3z, nc, ni] if se == 0 else [elementsCount[se] - ni, nc, n1z]
                elif constant_index == 3:
                    ids = [nc, 0, ni] if se == 0 else [nc, ni, n1z]

                if ni < len(dStart[se]):
                    if dStart[se][ni]:
                        tx.append(btx[ids[0]][ids[1]][ids[2]])
                        if dStart[0][ni] > 0:
                            td.append(btd[abs(dStart[se][ni])][ids[0]][ids[1]][ids[2]])
                        else:
                            td.append(vector.scaleVector(btd[abs(dStart[se][ni])][ids[0]][ids[1]][ids[2]], -1))
                elif ni > elementsCount[se] - len(dEnd[se]):
                    nie = ni - elementsCount[se] + len(dEnd[se]) - 1
                    if dEnd[se][nie]:
                        tx.append(btx[ids[0]][ids[1]][ids[2]])
                        if dEnd[se][nie] > 0:
                            td.append(btd[abs(dEnd[se][nie])][ids[0]][ids[1]][ids[2]])
                        else:
                            td.append(vector.scaleVector(btd[abs(dEnd[se][nie])][ids[0]][ids[1]][ids[2]], -1))
                else:
                    tx.append(btx[ids[0]][ids[1]][ids[2]])
                    if dBetween[se][0] > 0:
                        td.append(btd[abs(dBetween[se][0])][ids[0]][ids[1]][ids[2]])
                    else:
                        td.append(vector.scaleVector(btd[abs(dBetween[se][0])][ids[0]][ids[1]][ids[2]], -1))

        td = smoothCubicHermiteDerivativesLine(tx, td, fixStartDirection=True, fixEndDirection=True)

        nit = 0
        for se in range(2):
            for ni in range(elementsCount[se] + 1):
                if constant_index == 1:
                    ids = [ni, 0, nc] if se == 0 else [n3z, ni, nc]
                elif constant_index == 2:
                    ids = [n3z, nc, ni] if se == 0 else [elementsCount[se] - ni, nc, n1z]
                elif constant_index == 3:
                    ids = [nc, 0, ni] if se == 0 else [nc, ni, n1z]

                if ni < len(dStart[se]):
                    if dStart[se][ni]:
                        if dStart[0][ni] > 0:
                            btd[abs(dStart[se][ni])][ids[0]][ids[1]][ids[2]] = td[nit]
                        else:
                            btd[abs(dStart[se][ni])][ids[0]][ids[1]][ids[2]] = vector.scaleVector(td[nit], -1)
                        nit += 1
                elif ni > elementsCount[se] - len(dEnd[se]):
                    nie = ni - elementsCount[se] + len(dEnd[se]) - 1
                    if dEnd[se][nie]:
                        if dEnd[se][nie] > 0:
                            btd[abs(dEnd[se][nie])][ids[0]][ids[1]][ids[2]] = td[nit]
                        else:
                            btd[abs(dEnd[se][nie])][ids[0]][ids[1]][ids[2]] = vector.scaleVector(td[nit], -1)
                        nit += 1
                else:
                    if dBetween[se][0] > 0:
                        btd[abs(dBetween[se][0])][ids[0]][ids[1]][ids[2]] = td[nit]
                    else:
                        btd[abs(dBetween[se][0])][ids[0]][ids[1]][ids[2]] = vector.scaleVector(td[nit], -1)
                    nit += 1

    def calculate_interior_quadruple_point(self):
        """

        :return:
        """
        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n1z = self._elementsCount[1]
        n1y = n1z - 1
        n3z = self._elementsCount[2]
        n3y = n3z - 1
        n2z = self._elementsCount[0]

        if self._elementsCount[2] == min(self._elementsCount):
            cx = 3
        elif self._elementsCount[1] == min(self._elementsCount):
            cx = 2
        else:
            cx = 1
        n3r0, n2r0, n1r0 = self.get_triple_curves_end_node_parameters(0, cx=cx, index_output=True)
        n3r, n2r, n1r = self.get_triple_curves_end_node_parameters(1, cx=cx, index_output=True)

        ts = vector.magnitude(vector.addVectors([btx[n3r0][n2r0][n1r0], btx[n3r][n2r][n1r]], [1, -1]))
        ra = vector.magnitude(btx[n3z][0][n1z])
        x = vector.scaleVector(btx[n3z][0][n1z], (1 - ts/ra))
        n3r0, n2r0, n1r0 = self.get_triple_curves_end_node_parameters(0, index_output=True)
        n3r1, n2r1, n1r1 = self.get_triple_curves_end_node_parameters(1, index_output=True)
        btx[n3r0][n2r0][n1r0] = x
        btd1[n3r0][n2r0][n1r0] = [-(btx[n3r1][n2r1][n1r1][0] - btx[n3r0][n2r0][n1r0][0]), 0.0, 0.0]
        btd2[n3r0][n2r0][n1r0] = [0.0, 0.0, (btx[n3r1][n2r1][n1r1][2] - btx[n3r0][n2r0][n1r0][2])]
        btd3[n3r0][n2r0][n1r0] = [0.0, (btx[n3r1][n2r1][n1r1][1] - btx[n3r0][n2r0][n1r0][1]), 0.0]

    def sample_interior_curves(self):
        """

        :return:
        """
        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n1z = self._elementsCount[1]
        n1y = n1z - 1
        n3z = self._elementsCount[2]
        n3y = n3z - 1
        n2z = self._elementsCount[0]

        # btd = {1: btd1, 2: btd2, 3: btd3}


        #
        # n3s = [[n3y, n3y, 0, 0], [n3y, n3y, 0, 0], [0, n3y, 0, n3y]]
        # n2s = [[1, 1, 1, 1], [1, n2z, 1, n2z], [1, 1, 1, 1]]
        # n1s = [[0, n1y, 0, n1y], [n1y, n1y, n1y, n1y], [n1y, n1y, 0, 0]]
        #
        # n3x2 = n3y
        # n3d1 = 0
        # n2x1 = 1
        # n2d1 = 1
        # n1x2 = n1y
        # for nic in range(3):
        #     n3x1 = 0 if nic == 2 else n3y
        #     n3d2 = n3y if nic == 2 else 0
        #     n2x2 = n2z if nic == 1 else 1
        #     n2d2 = n2z if nic == 1 else 1
        #     n1x1 = 0 if nic == 1 else n1y
        #     n1d1 = n1y if nic == 1 else 0
        #     n1d2 = 0 if nic == 2 else n1y
        #
        #     tx, td = sampleCubicHermiteCurves([btx[n3x1][n2x1][n1x1], btx[n3x2][n2x2][n1x2]],
        #                                        [btd[nic][n3d1][n2d1][n1d1], btd[nic][n3d2][n2d2][n1d2]], self._elementsCount[nic] - 1)[:2]
        #
        #     for ni in range(ni0, self._elementsCount[nic] - 1 + ni0)

        tx, td1 = sampleCubicHermiteCurves([btx[n3y][1][n1y], btx[n3y][n2z][n1y]],
                                           [btd1[0][1][n1y], btd1[0][n2z][n1y]], self._elementsCount[0] - 1)[:2]

        for n2 in range(2, self._elementsCount[0]):
            btx[n3y][n2][n1y] = tx[n2-1]
            btd1[n3y][n2][n1y] = td1[n2-1]
            btd2[n3y][n2][n1y] = [0.0, 0.0, (btx[n3z][n2][n1z][2] - btx[n3y][n2][n1y][2])]
            btd3[n3y][n2][n1y] = [0.0, (btx[n3z][n2][n1z][1] - btx[n3y][n2][n1y][1]), 0.0]

        # curve 2 and parallel curves TODO change all [1] to n3y.
        for n2 in range(1, self._elementsCount[0]):
            tx, td3 = sampleCubicHermiteCurves([btx[n3y][n2][0], btx[n3y][n2][n1y]],
                                               [btd3[0][n2][0], btd3[0][n2][n1y]], self._elementsCount[1] - 1)[:2]

            for n1 in range(1, self._elementsCount[1] - 1):
                btx[n3y][n2][n1] = tx[n1]
                btd3[n3y][n2][n1] = td3[n1]

        for n2 in range(1, self._elementsCount[0]):
            for n1 in range(1, self._elementsCount[1] - 1):
                if n2 == 1:
                    btd1[n3y][n2][n1] = [btx[n3y][1][n1][0] - btx[n3z][0][n1][0], 0.0, 0.0]
                    btd2[n3y][n2][n1] = [0.0, 0.0, -btx[n3y][1][n1][2] + btx[n3z][0][n1][2]]
                else:
                    btd1[n3y][n2][n1] = vector.addVectors([btx[n3y][n2][n1], btx[n3y][n2+1][n1]], [-1, 1])
                    btd2[n3y][n2][n1] = vector.addVectors([btx[n3y][n2][n1], btx[0][n2][n1]], [1, -1])




        # sample along curve0_3
        for n2 in range(1, self._elementsCount[0]):
            for n1 in range(1, self._elementsCount[1]):
                tx, td2 = sampleCubicHermiteCurves([btx[0][n2][n1], btx[n3y][n2][n1]],
                                                   [btd2[0][n2][0], btd2[n3y][n2][0]], self._elementsCount[2]-1)[:2]

                for n3 in range(1, self._elementsCount[2] - 1):
                    btx[n3][n2][n1] = tx[n3]
                    btd2[n3][n2][n1] = td2[n3]

        for n3 in range(1, self._elementsCount[2] - 1):
            for n2 in range(1, self._elementsCount[0]):
                for n1 in range(1, self._elementsCount[1]):
                    if n2 == 1 and n1 == n1y:
                        btd1[n3][n2][n1] = [btx[n3][n2][n1][0] - btx[n3][n2-1][n1+1][0], 0.0, 0.0]
                        btd3[n3][n2][n1] = [0.0, btx[n3][n2-1][n1+1][1] - btx[n3][n2][n1][1], 0.0]
                    else:
                        btd1[n3][n2][n1] = vector.addVectors([btx[n3][n2+1][n1], btx[n3][n2][n1]], [1, -1])
                        btd3[n3][n2][n1] = vector.addVectors([btx[n3][n2][n1+1], btx[n3][n2][n1]], [1, -1])

    def smooth_regular_interior_curves(self):
        """

        :return:
        """
        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n1z = self._elementsCount[1]
        n1y = n1z - 1
        n3z = self._elementsCount[2]
        n3y = n3z - 1
        n2z = self._elementsCount[0]

        # smooth d1 in regular 1
        if self._elementsCount[0] >= 3:
            for n3 in range(1, self._elementsCount[2]):
                for n1 in range(1, self._elementsCount[1]):
                    tx = []
                    td1 = []
                    for n2 in range(1, self._elementsCount[0]+1):
                        tx.append(btx[n3][n2][n1])
                        td1.append(btd1[n3][n2][n1])
                    td1 = smoothCubicHermiteDerivativesLine(tx, td1, fixEndDirection=True)
                    for n2 in range(1, self._elementsCount[0]+1):
                        btd1[n3][n2][n1] = td1[n2-1]
        else:
            for n3 in range(1, self._elementsCount[2]):
                for n1 in range(1, self._elementsCount[1]):
                    btd1[n3][1][n1] = vector.addVectors([btx[n3][2][n1], btx[n3][1][n1]], [1, -1])
                    btd1[n3][2][n1] = vector.setMagnitude(btd1[n3][2][n1], vector.magnitude(btd1[n3][1][n1]))

        # smooth d3 in regular
        if self._elementsCount[1] >= 3:
            for n3 in range(1, self._elementsCount[2]):
                for n2 in range(1, self._elementsCount[0]):
                    tx = []
                    td3 = []
                    for n1 in range(self._elementsCount[1]):
                        tx.append(btx[n3][n2][n1])
                        td3.append(btd3[n3][n2][n1])

                    td3 = smoothCubicHermiteDerivativesLine(tx, td3, fixStartDirection=True)

                    for n1 in range(self._elementsCount[1]):
                        btd3[n3][n2][n1] = td3[n1]
        else:
            for n3 in range(1, self._elementsCount[2]):
                for n2 in range(1, self._elementsCount[0]):
                    btd3[n3][n2][1] = vector.addVectors([btx[n3][n2][1], btx[n3][n2][0]], [1, -1])
                    btd3[n3][n2][0] = vector.setMagnitude(btd3[n3][n2][0], vector.magnitude(btd3[n3][n2][1]))

        # regular curves d2
        for n2 in range(1, self._elementsCount[0]):
            for n1 in range(1, self._elementsCount[1]):
                if self._elementsCount[2] >= 3:
                    tx = []
                    td2 = []
                    for n3 in range(self._elementsCount[2]):
                        tx.append(btx[n3][n2][n1])
                        td2.append(btd2[n3][n2][n1])
                    td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixStartDirection=True)
                    for n3 in range(self._elementsCount[2]):
                        btd2[n3][n2][n1] = td2[n3]
                else:
                    btd2[1][n2][n1] = vector.addVectors([btx[1][n2][n1], btx[0][n2][n1]], [1, -1])
                    btd2[0][n2][n1] = vector.setMagnitude(btd2[0][n2][n1], vector.magnitude(btd2[1][n2][n1]))

    def get_triple_curves_end_node_parameters(self, rx, cx=None, index_output=False):
        """
        :param cx:
        :param rx:
        :return:
        """
        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n3z = self._elementsCount[2]
        n2z = self._elementsCount[0]
        n1z = self._elementsCount[1]
        n3y = n3z - 1
        n1y = n1z - 1

        if not cx:
            n3r = n3y + rx
            n2r = 1 - rx
            n1r = n1y + rx
        else:
            if cx == 1:
                n3r = n3y + rx
                n2r = n2z
                n1r = n1y + rx
            elif cx == 2:
                n3r = n3y + rx
                n2r = 1 - rx
                n1r = 0
            elif cx == 3:
                n3r = 0
                n2r = 1 - rx
                n1r = n1y + rx
            else:
                raise ValueError("curve index must be 1,2 or 3.")

        if index_output:
            return n3r, n2r, n1r
        else:
            return btx[n3r][n2r][n1r], btd1[n3r][n2r][n1r], btd2[n3r][n2r][n1r], btd3[n3r][n2r][n1r]

    def smooth_derivatives_to_surface(self):
        '''
        Smooth derivatives leading to quadruple point where 3 hex elements merge.
        :param n3: Index of through-wall coordinates to use.
        '''
        n3z = self._elementsCount[2]
        n1z = self._elementsCount[1]

        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        for n3 in range(1, self._elementsCount[2] + 1):
            for n2 in range(self._elementsCount[0]):
                for n1 in range(1, self._elementsCount[1]+1):
                    if self.on_sphere(n3, n2, n1):
                        # find indices of the neighbour node inside and contribution of its derivatives.
                        n3r = n3 - 1 if n3 == n3z else n3
                        n2r = n2 + 1 if n2 == 0 else n2
                        n1r = n1 - 1 if n1 == n1z else n1
                        co = [-1, 1, 1]
                        co[0] = -1 if n2 == 0 else 0
                        co[1] = 1 if n3 == n3z else 0
                        co[2] = 1 if n1 == n1z else 0

                        tx = []
                        td3 = []
                        tx.append(btx[n3r][n2r][n1r])
                        td3.append(
                            [(co[0]*btd1[n3r][n2r][n1r][c] + co[1]*btd2[n3r][n2r][n1r][c] + co[2]*btd3[n3r][n2r][n1r][c]) for c in range(3)])

                        tx.append(btx[n3][n2][n1])
                        td3.append(btd3[n3][n2][n1])

                        td3 = smoothCubicHermiteDerivativesLine(tx, td3, fixStartDirection=True, fixEndDirection=True)
                        btd3[n3][n2][n1] = td3[1]

    def on_sphere(self, n3, n2, n1):
        """
        Check if the given point is on the sphere.
        :param n3, n2, n1: node indexes in data structure [n3][n2][n1]
        :return: True, if it is on the sphere.
        """
        n3z = self._elementsCount[2]
        n1z = self._elementsCount[1]

        btx = self._shield3D.px

        return (n3 == n3z or n2 == 0 or n1 == n1z) and btx[n3][n2][n1]

    def generateNodes(self, nodes, fieldModule, coordinates):
        """
        Create cylinder nodes from coordinates.
        :param nodes: nodes from coordinates.
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self._shield3D.generateNodes(fieldModule, coordinates, nodeIdentifier)
        self._endNodeIdentifier = nodeIdentifier

    def generateElements(self, mesh, fieldModule, coordinates):
        """
        Create cylinder elements from nodes.
        :param mesh:
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        self._startElementIdentifier = elementIdentifier
        elementIdentifier = self._shield3D.generateElements(fieldModule, coordinates, elementIdentifier, [])
        self._endElementIdentifier = elementIdentifier


def calculate_azimuth(theta, theta_p):
    """
    Given polar angles of a point on the sphere surfaces, calculate the azimuth angle.
    :param theta: polar angle. In orthonormal coordinate system (axis1, axis2, axis3) with right-hand rule,
    theta is angle between common axis and point projection on plane of theta. In case theta=theta_3 and theta_p = theta_1, theta is between axis2 and projection
    :param theta_p: polar angle wrt other direction.
    :return: Azimuth angle.
    """
    return math.atan(1/(math.tan(theta_p)*math.cos(theta)))


def local_orthogonal_unit_vectors(x, axis3):
    """
    Find local orthogonal unit vectors for a point on a sphere
    :param x: coordinates of the point.
    :param axis3: The third axis in Cartesian coordinate system (axis1, axis2, axis3)
    :return: e1, e2, e3. Unit vectors. e3 is normal to the boundary, e2 is in (e3, axis3) plane and e1 normal to them.
    """
    e3 = vector.normalise(x)
    e2 = vector.vectorRejection(axis3, e3)
    e2 = vector.normalise(e2)
    e1 = vector.crossproduct3(e2, e3)

    return e1, e2, e3


def calculate_arc_length(x1, x2):
    """
    Calculate the arc length between points x1 and x2.
    :param x1, x2: points coordinates.
    :return: arc length
    """
    radius = vector.magnitude(x1)
    angle = vector.angleBetweenVectors(x1, x2)
    return radius * angle


def sample_curves_on_sphere(x1, x2, elementsOut):
    """

    :param x1, x2: points coordinates.
    :param elementsOut:
    :return:
    """
    deltax = vector.addVectors([x1, x2], [-1, 1])
    normal = vector.crossproduct3(x1, deltax)
    angle = vector.angleBetweenVectors(x1, x2)
    anglePerElement = angle/elementsOut
    arcLengthPerElement = calculate_arc_length(x1, x2)/elementsOut

    nx = []
    nd1 = []
    for n1 in range(elementsOut + 1):
        radiansAcross = n1 * anglePerElement
        x = vector.rotateVectorAroundVector(x1, normal, radiansAcross)
        d1 = vector.setMagnitude(vector.crossproduct3(normal, x), arcLengthPerElement)
        nx.append(x)
        nd1.append(d1)

    return nx, nd1


def spherical_to_cartesian(r, theta, phi):
    """
    :param r: Radius.
    :param theta: in radians.
    :param phi: azimuth angle in radians
    :return: x=[x1, x2, x3] coordinates.
    """
    return [r*math.sin(phi)*math.cos(theta), r*math.sin(phi)*math.sin(theta), r*math.cos(phi)]


def cartesian_to_spherical(x):
    """
    :return: [r, theta, phi].
    """
    r = vector.magnitude(x)
    theta = math.atan2(x[1], x[0])
    phi = math.acos(x[2]/r)
    return r, theta, phi


def local_to_global_coordinates(local_x, local_axes, local_origin=None):
    """
    Get global coordinates of a point with local coordinates x = [x1, x2, x3] and axes of local coordinate system.
    :param local_x: Coordinates in local coordinates system as a list of 3 components.
    :param local_origin: Origin of local coordinates system specified as a list of 3 components wrt global coordinates system.
    :param local_axes: Axes of local coordinates system, specified as a list of list 3X3 with respect to global coordinates system.
    :return: Global coordinate system.
    """
    if local_origin is None:
        local_origin = [0.0, 0.0, 0.0]
    return vector.addVectors([vector.addVectors(local_axes, local_x), local_origin])


def intersection_of_two_great_circles_on_sphere(p1, q1, p2, q2):
    """
    Find the intersection between arcs P1Q1 and P2Q2 on sphere.
    :param p1, q1, p2, q2: arcs extremities coordinates.
    :return: Point Sx, intersection between the arcs.
    """
    normal_to_plane_OP1Q1 = vector.crossproduct3(p1, q1)
    normal_to_plane_OP2Q2 = vector.crossproduct3(p2, q2)

    planes_intersection_vector = vector.crossproduct3(normal_to_plane_OP1Q1, normal_to_plane_OP2Q2)
    if vector.magnitude(planes_intersection_vector) == 0:
        sx = None
    else:
        sx = vector.setMagnitude(planes_intersection_vector, vector.magnitude(p1))
        p1q1_angle = vector.angleBetweenVectors(p1, q1)
        p1s_angle = vector.angleBetweenVectors(p1, sx)
        p2s_angle = vector.angleBetweenVectors(p2, sx)
        if p1s_angle > p1q1_angle or p2s_angle > p1q1_angle:
            sx = vector.scaleVector(sx, -1)

    return sx


def point_projection_on_sphere(p1, radius):
    """
    Find closest point to p1 on the sphere.
    :param p1: point.
    :return:
    """
    return vector.setMagnitude(p1, radius)
