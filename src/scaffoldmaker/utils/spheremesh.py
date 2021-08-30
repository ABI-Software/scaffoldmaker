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
    smoothCubicHermiteDerivativesLine, interpolateSampleLinear
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

        self._shield3D = ShieldMesh3D(self._elementsCount, elementsCountRim,
                 shieldMode=ShieldShape3D.SHIELD_SHAPE_OCTANT_PPP)

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

            self.copyEllipsesNodesToShieldNodes(ellipse, alongAxis, i)


        self.createAdditionalPointsForIncreasingElementsCount()



    def copyEllipsesNodesToShieldNodes(self, ellipse, alongAxis, ellipsenumber):
        """
        Copy coordinates and derivatives of ellipse to shield.
        :param n3: the index number of ellipse along the central path.
        """

        shield = ellipse.getShield()

        # Modify the shield12 to get only the quarter that you want. make others None so you can generate the nodes.
        if ellipsenumber == 0:
            for n2 in range(self._elementsCount[0] + 1):  # TODO modify this to number of elements.
                for n1 in range(self._elementsCount[1] + 1):
                    # only first quadrant is needed  TODO we need to modify this to number of elements half of the elments across each direction.
                    self._shield3D.pd2[0][n2][n1] = [c for c in alongAxis]
                    if n2 == 0 and n1 == self._elementsCount[1] - 1:
                        n1s = n1 + 1
                    else:
                        n1s = n1
                    n1e = n1 + self._elementsCount[1]
                    if shield.px[0][n2][n1e]:
                        self._shield3D.px[0][n2][n1s] = [c for c in shield.px[0][n2][n1e]]
                        self._shield3D.pd1[0][n2][n1s] = [c for c in shield.pd1[0][n2][n1e]]
                        self._shield3D.pd3[0][n2][n1s] = [c for c in shield.pd3[0][n2][n1e]]
        elif ellipsenumber == 1:
            n2s = self._elementsCount[0]
            for n3 in range(self._elementsCount[2] + 1):
                for n1 in range(self._elementsCount[1] + 1):
                    self._shield3D.pd2[n3][n2s][n1] = [c for c in alongAxis]
                    if n3 == self._elementsCount[2] and n1 == self._elementsCount[1] - 1:
                        n1s = n1 + 1
                    else:
                        n1s = n1
                    n2e = n3 + self._elementsCount[2]
                    n1e = n1 + self._elementsCount[1]
                    if n3 == 0:
                        self._shield3D.pd2[n3][n2s][n1s] = [c for c in shield.pd1[0][n2e][n1e]]
                    else:
                        if shield.px[0][n2e][n1 + self._elementsCount[1]]:
                            self._shield3D.px[n3][n2s][n1s] = [c for c in shield.px[0][n2e][n1e]]
                            self._shield3D.pd1[n3][n2s][n1s] = [c for c in shield.pd1[0][n2e][n1e]]
                            self._shield3D.pd3[n3][n2s][n1s] = [c for c in shield.pd3[0][n2e][n1e]]
        elif ellipsenumber == 2:
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
                            self._shield3D.pd2[n3][n2][n1s] = [c for c in shield.pd1[0][n2][n1e]]
                        elif 0 < n2 < self._elementsCount[0]:
                            self._shield3D.pd2[n3][n2][n1s] = [-c for c in shield.pd3[0][n2][n1e]]
                    else:
                        if n2 < self._elementsCount[0]:
                            self._shield3D.pd2[n3][n2][n1s] = [c for c in alongAxis]
                            if shield.px[0][n2][n1e]:
                                self._shield3D.px[n3s][n2][n1s] = [c for c in shield.px[0][n2][n1e]]
                                self._shield3D.pd1[n3s][n2][n1s] = [c for c in shield.pd1[0][n2][n1e]]
                                self._shield3D.pd3[n3s][n2][n1s] = [c for c in shield.pd3[0][n2][n1e]]
                        else:
                            self._shield3D.pd2[n3][n2][n1s] = [c for c in shield.pd1[0][n2][n1e]]

    def createAdditionalPointsForIncreasingElementsCount(self):
        """

        :return:
        """

        self.calculate_regular_nodes()
        self._shield3D.getQuadruplePoint()
        self.fixD2DerivativesOnTheEllipses()
        self.remapDerivativesOnTheEllipses()

        # reverse d2 for 0,0,0
        self._shield3D.pd2[0][0][0] = [-c for c in self._shield3D.pd2[0][0][0]]

        self.calculate_surface_nodes()
        self.smooth_triple_curves()

    def calculate_regular_nodes(self):
        """

        :return:
        """
        n3a = 0
        n3b = self._elementsCount[2] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2b = self._elementsCount[0]
        n1a = 0
        n1b = self._elementsCount[1] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)

        if self._elementsCount[0] == 3:
            self.CreateRegularNodeOnBoundary()
            n3 = 0
            n1 = 0
            for n2 in range(n2a + 1, n2b):
                # First create the node for regular node. We are using the neighbour nodes x,y,z to do it but  TODO eventully we need to use the sampling between boundary nodes.
                x = [self._shield3D.px[n3][n2][n1 + 1][0], self._shield3D.px[n3 + 1][n2 + 1][n1 + 1][1],
                     self._shield3D.px[n3 + 1][n2][n1][2]]
                self._shield3D.px[n3 + 1][n2][n1 + 1] = x
                self._shield3D.pd1[n3 + 1][n2][n1 + 1] = [
                    (self._shield3D.px[n3 + 1][n2 + 1][n1 + 1][c] - self._shield3D.px[n3 + 1][n2][n1 + 1][c]) for c in
                    range(3)]
                self._shield3D.pd2[n3 + 1][n2][n1 + 1] = [
                    -(self._shield3D.px[n3][n2][n1 + 1][c] - self._shield3D.px[n3 + 1][n2][n1 + 1][c]) for c in
                    range(3)]
                self._shield3D.pd3[n3 + 1][n2][n1 + 1] = [
                    -(self._shield3D.px[n3 + 1][n2][n1][c] - self._shield3D.px[n3 + 1][n2][n1 + 1][c]) for c in
                    range(3)]
                # using new point obtain the d2 derivatives of the neighbour nodes
                self._shield3D.pd2[n3 + 1][n2][n1] = [
                    self._shield3D.px[n3 + 1][n2][n1 + 1][c] - self._shield3D.px[n3 + 1][n2][n1][c] for c in range(3)]
                self._shield3D.pd2[n3][n2][n1 + 1] = [
                    self._shield3D.px[n3 + 1][n2][n1 + 1][c] - self._shield3D.px[n3][n2][n1 + 1][c] for c in range(3)]
                self._shield3D.pd2[n3 + 1][n2 + 1][n1 + 1] = [
                    self._shield3D.px[n3 + 1][n2 + 1][n1 + 1][c] - self._shield3D.px[n3][n2 + 1][n1 + 1][c] for c in
                    range(3)]

            self.remapRegularCurvesOnEllipses()

    def CreateRegularNodeOnBoundary(self):
        """

        :return:
        """
        n3a = 0
        n3b = self._elementsCount[2] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2b = self._elementsCount[0]
        n1a = 0
        n1b = self._elementsCount[1] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)

        # create the other two nodes to create the elements.
        radius = self._radius[0]  # TODO need to be changed for spheroid

        elementsAroundEllipse12 = self._elementsCount[0] + self._elementsCount[1] - 2
        radiansAroundEllipse12 = math.pi / 2
        radiansPerElementAroundEllipse12 = radiansAroundEllipse12 / elementsAroundEllipse12
        elementsAroundEllipse13 = self._elementsCount[0] + self._elementsCount[2] - 2
        radiansAroundEllipse13 = math.pi / 2
        radiansPerElementAroundEllipse13 = radiansAroundEllipse13 / elementsAroundEllipse13

        n3 = 1
        n1 = 0
        for n2 in range(n2a + 1, n2b):
            theta_1 = math.pi / 4  # TODO actually it should change with the number of elements.
            theta_3 = 2*radiansPerElementAroundEllipse12
            phi_3 = calculate_azimuth(math.pi/2 - theta_3, theta_1)
            # We assume it is a sphere not a spheroid for now. TODO Use the relations for spheroid instead
            x = [radius * math.sin(phi_3) * math.cos(theta_3), radius * math.sin(phi_3) * math.sin(theta_3),
                 radius * math.cos(phi_3)]

            a1, a2, a3 = local_orthogonal_unit_vectors(x, self._axes[2])
            self._shield3D.px[n3 + 1][n2][n1 + 2] = x
            self._shield3D.pd1[n3 + 1][n2][n1 + 2] = a1
            self._shield3D.pd2[n3 + 1][n2][n1 + 2] = a2
            self._shield3D.pd3[n3 + 1][n2][n1 + 2] = a3

            arcLength = radius * phi_3
            self._shield3D.pd2[n3 + 1][n2][n1 + 2] = vector.setMagnitude(self._shield3D.pd2[n3 + 1][n2][n1 + 2], arcLength)

    def remapRegularCurvesOnEllipses(self):
        """

        :return:
        """
        n3a = 0
        n3b = self._elementsCount[2] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2b = self._elementsCount[0]
        n1a = 0
        n1b = self._elementsCount[1] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)

        # We need to remap the derivatives as we want d1 to be in -axis1 direction, d2, axis3 and d3, axis2. We need to think about the directions as well. e.g.
        # d1, d2, d3 could be simply axis1,axis2 and axis3. However this one was chosen so it could be consistent with the cylinder mesh that makes joining them
        # a little easier. TODO d1,d2,d3 directions in regular elements.
        # for nodes on the ellipse2, d2 is -d3 and d3 is d2. TODO we can add another method to ellipse so it gives us what we expect for derivatives.
        n3 = 0
        n1 = 0
        for n2 in range(n2a + 1, n2b):
            temp = self._shield3D.pd2[n3 + 1][n2][n1]
            self._shield3D.pd2[n3 + 1][n2][n1] = [-c for c in self._shield3D.pd3[n3 + 1][n2][n1]]
            self._shield3D.pd3[n3 + 1][n2][n1] = temp
            # for ellipse 1, swap d2 and d1
            temp = self._shield3D.pd2[n3 + 1][n2 + 1][n1 + 1]
            self._shield3D.pd2[n3 + 1][n2 + 1][n1 + 1] = [c for c in self._shield3D.pd1[n3 + 1][n2 + 1][n1 + 1]]
            self._shield3D.pd1[n3 + 1][n2 + 1][n1 + 1] = temp

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

        self._shield3D.pd2[n3b][n2a][n1a] = [self._shield3D.px[n3b][n2a][n1b][c] - self._shield3D.px[n3b][n2a][n1a][c]
                                             for c in range(3)]
        self._shield3D.pd2[0][1][1] = [self._shield3D.px[1][1][1][c] - self._shield3D.px[0][1][1][c] for c in range(3)]
        self._shield3D.pd2[n3b][n2b][n1b] = [
            -(self._shield3D.px[n3b][n2a][n1b][c] - self._shield3D.px[n3b][n2b][n1b][c]) for c in range(3)]

    def remapDerivativesOnTheEllipses(self):
        """

        :return:
        """
        n3a = 0
        n3b = self._elementsCount[2] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2b = self._elementsCount[0]
        n1a = 0
        n1b = self._elementsCount[1] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)

        temp = self._shield3D.pd2[1][1][0]
        self._shield3D.pd2[1][1][0] = [-c for c in self._shield3D.pd3[1][1][0]]
        self._shield3D.pd3[1][1][0] = [c for c in temp]
        temp = self._shield3D.pd2[n3b][n2b][n1b]
        self._shield3D.pd2[n3b][n2b][n1b] = [c for c in self._shield3D.pd1[n3b][n2b][n1b]]
        self._shield3D.pd1[n3b][n2b][n1b] = temp
        temp = self._shield3D.pd2[n3b][n2b][n1a]
        self._shield3D.pd2[n3b][n2b][n1a] = [c for c in self._shield3D.pd1[n3b][n2b][n1a]]
        self._shield3D.pd1[n3b][n2b][n1a] = temp

    def calculate_surface_nodes(self):
        """

        :return:
        """
        n3a = 0
        n3b = self._elementsCount[2] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2b = self._elementsCount[0]
        n1a = 0
        n1b = self._elementsCount[1] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)

        radius = self._radius[0]  # TODO need to be changed for spheroid

        elementsAroundEllipse12 = self._elementsCount[0] + self._elementsCount[1] - 2
        radiansAroundEllipse12 = math.pi / 2
        radiansPerElementAroundEllipse12 = radiansAroundEllipse12 / elementsAroundEllipse12
        elementsAroundEllipse13 = self._elementsCount[0] + self._elementsCount[2] - 2
        radiansAroundEllipse13 = math.pi / 2
        radiansPerElementAroundEllipse13 = radiansAroundEllipse13 / elementsAroundEllipse13

        theta_2 = radiansPerElementAroundEllipse13  # TODO actually it should change with the number of elements.
        theta_3 = radiansPerElementAroundEllipse12
        phi_3 = calculate_azimuth(theta_3, theta_2)
        # We assume it is a sphere not a spheroid for now. TODO Use the relations for spheroid instead
        x = [radius * math.sin(phi_3) * math.cos(theta_3), radius * math.sin(phi_3) * math.sin(theta_3),
             radius * math.cos(phi_3)]

        a1, a2, a3 = local_orthogonal_unit_vectors(x, self._axes[2])
        self._shield3D.px[n3b + 1][0][n1b + 1] = x
        self._shield3D.pd1[n3b + 1][0][n1b + 1] = a1
        self._shield3D.pd2[n3b + 1][0][n1b + 1] = a2
        self._shield3D.pd3[n3b + 1][0][n1b + 1] = a3

        arcLength = calculate_arc_length(x, self._shield3D.px[n3a][0][n1b + 1], radius)
        self._shield3D.pd2[n3b + 1][0][n1b + 1] = vector.setMagnitude(self._shield3D.pd2[n3b + 1][0][n1b + 1], arcLength)
        self._shield3D.smoothDerivativesToSurfaceQuadruple(n3b+1)

        self._shield3D.pd1[n3b + 1][0][n1b + 1] = [
            self._shield3D.px[n3b + 1][n2b][n1b + 1][c] - self._shield3D.px[n3b + 1][0][n1b + 1][c] for c in range(3)]
        # self._shield3D.pd2[n3b + 1][0][n1b + 1] = [
        #     -(self._shield3D.px[0][0][n1b + 1][c] - self._shield3D.px[n3b + 1][0][n3b + 1][c]) for c in range(3)]
        # self._shield3D.pd3[n3b + 1][0][n1b + 1] = [
        #     -(self._shield3D.px[1][1][1][c] - self._shield3D.px[n3b + 1][0][n1b + 1][c]) for c in range(3)]



        n1a = 1
        n2a = 1
        n3a = 1
        self._shield3D.pd1[n3a][n2a][n1a] = [-(self._shield3D.px[n3a+1][n2a-1][n1a+1][0] - self._shield3D.px[n3a][n2a][n1a][0]), 0.0, 0.0]
        self._shield3D.pd2[n3a][n2a][n1a] = [0.0, 0.0, (self._shield3D.px[n3a+1][n2a-1][n1a+1][2] - self._shield3D.px[n3a][n2a][n1a][2])]
        self._shield3D.pd3[n3a][n2a][n1a] = [0.0, (self._shield3D.px[n3a+1][n2a-1][n1a+1][1] - self._shield3D.px[n3a][n2a][n1a][1]), 0.0]

    def smooth_triple_curves(self):
        """

        :return:
        """
        n3a = 0
        n3b = self._elementsCount[2] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2a = (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)
        n2b = self._elementsCount[0]
        n1a = 0
        n1b = self._elementsCount[1] - (self._elementsCountAcrossShell + self._elementsCountAcrossTransition)

        radius = 1.0
        theta = math.atan(math.sqrt(1.0 / 2.0))
        arclength = radius * theta
        self._shield3D.pd2[0][0][n1b + 1] = vector.setMagnitude(self._shield3D.pd2[0][0][n1b + 1], arclength)
        self._shield3D.pd2[n3b + 1][0][0] = vector.setMagnitude(self._shield3D.pd2[n3b + 1][0][0], arclength)
        self._shield3D.pd2[n3b + 1][n2b][n1b + 1] = vector.setMagnitude(self._shield3D.pd2[n3b + 1][n2b][n1b + 1],
                                                                        arclength)
        self._shield3D.pd2[n3b + 1][0][n1b + 1] = vector.vectorRejection(self._shield3D.pd2[n3b + 1][0][n1b + 1],
                                                                         self._shield3D.pd3[n3b + 1][0][n1b + 1])
        self._shield3D.pd2[n3b + 1][0][n1b + 1] = vector.setMagnitude(self._shield3D.pd2[n3b + 1][0][n1b + 1],
                                                                      arclength)
        self._shield3D.pd1[n3b + 1][0][n1b + 1] = vector.crossproduct3(self._shield3D.pd2[n3b + 1][0][n1b + 1],
                                                                       self._shield3D.pd3[n3b + 1][0][n1b + 1])
        d2mag = vector.magnitude(self._shield3D.pd2[n3b + 1][0][n1b + 1])
        thetad1ncurve = math.pi / 6
        self._shield3D.pd1[n3b + 1][0][n1b + 1] = vector.setMagnitude(self._shield3D.pd1[n3b + 1][0][n1b + 1],
                                                                      d2mag / math.tan(thetad1ncurve))

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


def calculate_arc_length(x1, x2, radius):
    """
    Calculate the arc length between points x1 and x2.
    :param x1, x2: points coordinates.
    :param radius: sphere radius.
    :return: arc length
    """
    angle = vector.angleBetweenVectors(x1, x2)
    return radius * angle
