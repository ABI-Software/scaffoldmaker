"""
Utility functions for generating a solid spheroid.
"""

from enum import Enum
from scaffoldmaker.utils import vector, geometry
import math
from opencmiss.zinc.field import Field
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from scaffoldmaker.utils.shieldmesh import ShieldMesh3D, ShieldShape3D, ShieldRimDerivativeMode
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite, \
    smoothCubicHermiteDerivativesLine, interpolateSampleLinear, interpolateCubicHermite
from opencmiss.zinc.node import Node
from scaffoldmaker.utils.mirror import Mirror
from scaffoldmaker.meshtypes.meshtype_1d_path1 import extractPathParametersFromRegion
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

        self.calculateBoundaryElipses()
        self.generateNodes(nodes, fieldModule, coordinates)
        self.generateElements(mesh, fieldModule, coordinates)

    def calculateBoundaryElipses(self):
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
        self.sample_triple_curves()
        self._shield3D.getQuadruplePoint2()
        self.fixD2DerivativesOnTheEllipses()
        n3z = self._elementsCount[2]
        self._shield3D.smoothDerivativesToSurfaceQuadruple(n3z)
        self.smoothDerivativesToSurface()

    def calculate_surface_quadruple_point(self):
        """
        Calculate coordinates and derivatives of points where 3 hex elements merge.
        :return:
        """

        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n1z = self._elementsCount[1]
        n3z = self._elementsCount[2]

        radius = self._radius[0]  # TODO need to be changed for spheroid

        elementsAroundEllipse12 = self._elementsCount[0] + self._elementsCount[1] - 2
        radiansAroundEllipse12 = math.pi / 2
        radiansPerElementAroundEllipse12 = radiansAroundEllipse12 / elementsAroundEllipse12
        elementsAroundEllipse13 = self._elementsCount[0] + self._elementsCount[2] - 2
        radiansAroundEllipse13 = math.pi / 2
        radiansPerElementAroundEllipse13 = radiansAroundEllipse13 / elementsAroundEllipse13

        for n in range(min(self._elementsCount[0], self._elementsCount[1]) - 1):
            theta_2 = (n+1) * radiansPerElementAroundEllipse13
            theta_3 = (self._elementsCount[1] - 1) * radiansPerElementAroundEllipse12
            phi_3 = calculate_azimuth(theta_3, theta_2)
            # We assume it is a sphere not a spheroid for now. TODO Use the relations for spheroid instead
            x = spherical_to_cartesian(radius, theta_3, phi_3)

            a1, a2, a3 = local_orthogonal_unit_vectors(x, self._axes[2])
            n1 = self._elementsCount[1] - n if n == 0 else self._elementsCount[1] - n - 1
            n2 = n if n == 0 else n+1
            btx[n3z][n2][n1] = x
            btd1[n3z][n2][n1] = a1  # initialise
            btd2[n3z][n2][n1] = a2  # initialise
            btd3[n3z][n2][n1] = a3

    def sample_triple_curves(self):
        """
        Sample points on the triple curves of quadruple point on the sphere surface.
        :return:
        """
        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n3a = 0
        n1z = self._elementsCount[1]
        n2z = self._elementsCount[0]
        n3z = self._elementsCount[2]

        # sample on curve 1 of the triple curves and smooth the end derivatives.
        nx, nd1 = sample_curves_sphere(btx[n3z][0][n1z], btx[n3z][n2z][n1z], self._elementsCount[0] - 1)
        for nc in range(self._elementsCount[0]):
            if nc == 0:
                btd1[n3z][0][n1z] = nd1[0]
            elif nc == n2z - 1:
                btd2[n3z][n2z][n1z] = vector.scaleVector(nd1[-1], -1)
            else:
                btx[n3z][nc + 1][n1z] = nx[nc]
                btd1[n3z][nc + 1][n1z] = nd1[nc]

        # smooth d2 curve
        for n2 in range(n2a + 1, n2z):
            a1, a2, a3 = local_orthogonal_unit_vectors(btx[n3z][n2][n1z], self._axes[2])
            btd3[n3z][n2][n1z] = a3
            tx = []
            td2 = []
            for n1 in range(self._elementsCount[1] - 1):
                tx.append(btx[n3z][n2][n1])
                td2.append(btd2[n3z][n2][n1])
            tx.append(btx[n3z][n2][n1z])
            td2.append(vector.crossproduct3(btd3[n3z][n2][n1z],btd1[n3z][n2][n1z]))
            for n3 in range(self._elementsCount[2] - 2, -1, -1):
                tx.append(btx[n3][n2][n1z])
                td2.append(vector.scaleVector(btd2[n3][n2][n1z], -1))

            td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixStartDirection=True, fixEndDirection=True)

            btd2[n3z][n2][0] = td2[0]
            for n1 in range(1, self._elementsCount[1] - 1):
                btd2[n3z][n2][n1] = vector.scaleVector(td2[n1], -1)
            btd2[n3z][n2][n1z] = vector.scaleVector(td2[self._elementsCount[1] - 1], -1)
            for n3 in range(self._elementsCount[2] - 2, -1, -1):
                btd2[n3][n2][n1z] = vector.scaleVector(td2[-1 - n3], -1)

        # curve 2
        nx, nd1 = sample_curves_sphere(btx[n3z][0][0], btx[n3z][0][n1z], self._elementsCount[1] - 1)
        for n1 in range(self._elementsCount[1]+1):
            if n1 == n1z-1:
                continue
            if n1 == 0:
                btd2[n3z][0][0] = nd1[0]
            elif n1 == n1z:
                btd2[n3z][0][n1z] = vector.scaleVector(nd1[-1], -1)
            else:
                btx[n3z][0][n1] = nx[n1]
                btd2[n3z][0][n1] = vector.scaleVector(nd1[n1], -1)

        # smooth d1 curve
        for n1 in range(1, n1z - 1):
            a1, a2, a3 = local_orthogonal_unit_vectors(btx[n3z][0][n1], self._axes[2])
            btd3[n3z][0][n1] = a3
            tx = []
            td1 = []
            tx.append(btx[0][0][n1])
            td1.append(btd2[0][0][n1])
            for n3 in range(1, self._elementsCount[2] - 1):
                tx.append(btx[n3][0][n1])
                td1.append(btd1[n3][0][n1])
            tx.append(btx[n3z][0][n1])
            td1.append(vector.crossproduct3(btd2[n3z][0][n1], btd3[n3z][0][n1]))
            for n2 in range(2, self._elementsCount[0]):
                tx.append(btx[n3z][n2][n1])
                td1.append(btd1[n3z][n2][n1])
            tx.append(btx[n3z][n2z][n1])
            td1.append(vector.scaleVector(btd2[n3z][n2z][n1], -1))

            # for n2 in range(self._elementsCount[1] - 1):
            #     tx.append(btx[n3z][n2][n1])
            #     td1.append(btd2[n3z][n2][n1])
            # tx.append(btx[n3z][n2][n1z])
            # td1.append(vector.crossproduct3(btd3[n3z][n2][n1z],btd1[n3z][n2][n1z]))
            # for n3 in range(self._elementsCount[2] - 2, -1, -1):
            #     tx.append(btx[n3][n2][n1z])
            #     td1.append(vector.scaleVector(btd2[n3][n2][n1z], -1))

            td1 = smoothCubicHermiteDerivativesLine(tx, td1, fixStartDirection=True, fixEndDirection=True)
            btd2[0][0][n1] = td1[0]
            for n3 in range(1, self._elementsCount[2] - 1):
                btd1[n3][0][n1] = td1[n3]
            btd1[n3z][0][n1] = td1[self._elementsCount[2] - 1]
            for n2 in range(2, self._elementsCount[0]):
                btd1[n3z][n2][n1] = td1[self._elementsCount[2] + n2-2]
            btd2[n3z][n2z][n1] = vector.scaleVector(td1[-1], -1)



            # btd2[n3z][n2][0] = td1[0]
            # for n1 in range(1, self._elementsCount[1] - 1):
            #     btd2[n3z][n2][n1] = td1[n1]
            # btd2[n3z][n2][n1z] = vector.scaleVector(td1[1], -1)
            # for n3 in range(self._elementsCount[2] - 2, -1, -1):
            #     btd2[n3][n2][n1z] = vector.scaleVector(td1[-1 - n3], -1)





        # sample on curve 3 of the triple curves and smooth the end derivatives.
        arcLength = calculate_arc_length(btx[0][0][n1z], btx[n3z][0][n1z])
        btd2[0][0][n1z] = vector.setMagnitude(btd2[0][0][n1z], arcLength)

    def fixD2DerivativesOnTheEllipses(self):
        """

        :return:
        """
        n3a = 0
        n3b = self._elementsCount[2] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2b = self._elementsCount[0]
        n1a = 0
        n1b = self._elementsCount[1] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)

        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        n1 = self._elementsCount[1] - 1
        for n2 in range(self._elementsCount[0] + 1):
            if n2a <= n2 < self._elementsCount[0]:
                btd3[n3b][n2][n1a] = [btx[n3b][n2][n1b][c] - btx[n3b][n2][n1a][c] for c in range(3)]
                btd2[0][n2][n1] = [btx[1][n2][n1][c] - btx[0][n2][n1][c] for c in range(3)]

    def smoothDerivativesToSurface(self):
        '''
        Smooth derivatives leading to quadruple point where 3 hex elements merge.
        :param n3: Index of through-wall coordinates to use.
        '''
        n3a = 0
        n3z = self._elementsCount[2]
        n3b = self._elementsCount[2] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2b = self._elementsCount[0]
        n1a = 0
        n1b = self._elementsCount[1] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n1z = self._elementsCount[1]

        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        for n2 in range(self._elementsCount[0]):
            for n1 in range(1, self._elementsCount[1]+1):
                if btx[n3z][n2][n1]:
                    tx = []
                    td3 = []
                    if n2 == 0:
                        if n1 == n1z:
                            n2r = n2 + 1
                            n1r = n1 - 1
                            co = [-1, 1, 1]
                        else:
                            n2r = n2 + 1
                            n1r = n1
                            co = [-1, 1, 0]
                    else:
                        if n1 == n1z:
                            n2r = n2
                            n1r = n1 - 1
                            co = [0, 1, 1]
                        else:
                            n2r = n2
                            n1r = n1
                            co = [0, 1, 0]

                    tx.append(btx[n3z-1][n2r][n1r])
                    td3.append(
                        [(co[0]*btd1[n3z-1][n2r][n1r][c] + co[1]*btd2[n3z-1][n2r][n1r][c] + co[2]*btd3[n3z-1][n2r][n1r][c]) for c in range(3)])
                    for n3 in range(n3z, n3z + 1):
                        tx.append(btx[n3][n2][n1])
                        td3.append(btd3[n3][n2][n1])
                    td3 = smoothCubicHermiteDerivativesLine(tx, td3, fixStartDirection=True, fixEndDirection=True)
                    for n3 in range(n3z, n3z + 1):
                        btd3[n3][n2][n1] = td3[n3 - n3z + 1]

        # for n2 in range(self._elementsCount[0] + 1):
        #     if 1 < n2 < self._elementsCount[0]:
        #         tx = []
        #         td3 = []
        #         for n3 in range(n3b, n3b+2):
        #             n1 = self._elementsCount[1] - 2 + n3
        #             tx.append(btx[n3][n2][n1])
        #             if n3 == n3b:
        #                 td3.append(vector.addVectors(btd2[n3][n2][n1], btd3[n3][n2][n1]))
        #             else:
        #                 td3.append(btd3[n3][n2][n1])
        #
        #         td3 = smoothCubicHermiteDerivativesLine(tx, td3, fixStartDirection=True, fixEndDirection=True)
        #
        #         for nc3 in range(1, self._elementsCountAcrossShell + self._elementsCountAcrossTransition + 1):
        #             n3 = nc3 + n3b
        #             n1 = n3
        #             btd3[n3][n2][n1] = td3[nc3]

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


def sample_curves_sphere(x1, x2, elementsOut):
    """

    :param x1, x2: points coordinates.
    :param elementsOut:
    :return:
    """
    deltax = vector.addVectors(x1, x2, -1, 1)
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


