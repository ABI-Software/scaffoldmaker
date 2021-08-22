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
            assert (self._elementsCount[i] % 2 == 0), 'createSphereMesh3d: number of across elements' \
                                                          ' is not an even number'

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

        self.calculateBoundaryElipses(fieldModule, coordinates)
        # self.generateNodes(nodes, fieldModule, coordinates)
        # self.generateElements(mesh, fieldModule, coordinates)

    def calculateBoundaryElipses(self, fieldModule, coordinates):
        """

        :return:
        """

        centre = self._centre
        elementsCountAcrossMajor = 2 * self._elementsCount[1]
        elementsCountAcrossMinor = 2 * self._elementsCount[0]
        elementsCountAcrossShell = self._elementsCountAcrossShell
        elementsCountAcrossTransition = self._elementsCountAcrossTransition
        shellProportion = self._shellProportion

        ellipseAxes = [[self._axes[0], self._axes[1], self._axes[2]],
                       [[-c for c in self._axes[2]], self._axes[1], self._axes[0]],
                       [self._axes[0], [-c for c in self._axes[2]], self._axes[1]]]

        coreRadius = [(self._coreRadius[0], self._coreRadius[1]),
                      (self._coreRadius[2], self._coreRadius[1]),
                      (self._coreRadius[0], self._coreRadius[2])]

        for i in range(3):
            majorAxis = ellipseAxes[i][0]
            minorAxis = ellipseAxes[i][1]
            alongAxis = ellipseAxes[i][2]
            coreMajorRadius = coreRadius[i][0]
            coreMinorRadius = coreRadius[i][1]
            ellipse = Ellipse2D(centre, majorAxis, minorAxis,
                         elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                         elementsCountAcrossTransition, shellProportion, coreMajorRadius, coreMinorRadius,
                         ellipseShape=EllipseShape.Ellipse_SHAPE_FULL)

            self.copyEllipsesNodesToShieldNodes(ellipse, alongAxis, i)

        self._shield3D.getQudaruplePoint()
        self._shield3D.pd2[1][1][0] = [self._shield3D.px[1][1][1][c] - self._shield3D.px[1][1][0][c] for c in range(3)]
        self._shield3D.pd2[0][1][1] = [self._shield3D.px[1][1][1][c] - self._shield3D.px[0][1][1][c] for c in range(3)]
        self._shield3D.pd2[1][2][1] = [-(self._shield3D.px[1][1][1][c] - self._shield3D.px[1][2][1][c]) for c in range(3)]

        temp = self._shield3D.pd2[1][1][0]
        self._shield3D.pd2[1][1][0] = [-c for c in self._shield3D.pd3[1][1][0]]
        self._shield3D.pd3[1][1][0] = [c for c in temp]
        temp = self._shield3D.pd2[1][2][1]
        self._shield3D.pd2[1][2][1] = [c for c in self._shield3D.pd1[1][2][1]]
        self._shield3D.pd1[1][2][1] = temp
        temp = self._shield3D.pd2[1][2][0]
        self._shield3D.pd2[1][2][0] = [c for c in self._shield3D.pd1[1][2][0]]
        self._shield3D.pd1[1][2][0] = temp

        x = [0.57735026, 0.57735026, 0.57735026]
        self._shield3D.px[2][0][2] = x
        self._shield3D.pd1[2][0][2] = [self._shield3D.px[2][2][2][c] - self._shield3D.px[2][0][2][c] for c in range(3)]
        self._shield3D.pd2[2][0][2] = [-(self._shield3D.px[0][0][2][c] - self._shield3D.px[2][0][2][c]) for c in range(3)]
        self._shield3D.pd3[2][0][2] = [-(self._shield3D.px[1][1][1][c] - self._shield3D.px[2][0][2][c]) for c in range(3)]

        radius = 1.0
        theta = math.atan(math.sqrt(1.0/2.0))
        arclength = radius*theta
        self._shield3D.pd2[0][0][2] = vector.setMagnitude(self._shield3D.pd2[0][0][2], arclength)
        self._shield3D.pd2[1][0][0] = vector.setMagnitude(self._shield3D.pd2[1][0][0], arclength)
        self._shield3D.pd2[2][2][2] = vector.setMagnitude(self._shield3D.pd2[2][2][2], arclength)
        self._shield3D.pd2[2][0][2] = vector.vectorRejection(self._shield3D.pd2[2][0][2], self._shield3D.pd3[2][0][2])
        self._shield3D.pd2[2][0][2] = vector.setMagnitude(self._shield3D.pd2[2][0][2], arclength)
        self._shield3D.pd1[2][0][2] = vector.crossproduct3(self._shield3D.pd2[2][0][2], self._shield3D.pd3[2][0][2])
        self._shield3D.pd1[2][0][2] = vector.setMagnitude(self._shield3D.pd1[2][0][2], arclength)


        self._shield3D.generateNodes(fieldModule, coordinates, 1)
        self._shield3D.generateElements(fieldModule, coordinates, 1)

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
                    n1e = n1 + 2
                    if shield.px[0][n2][n1 + self._elementsCount[1]]:
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
                    n2e = n3 + 2
                    n1e = n1 + 2
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

                    if n2 == self._elementsCount[2] and n3 == self._elementsCount[2] - 1:
                        n3s = n3 + 1
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



    def generateNodes(self, nodes, fieldModule, coordinates):
        """
        Create cylinder nodes from coordinates.
        :param nodes: nodes from coordinates.
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self._shield.generateNodes(fieldModule, coordinates, nodeIdentifier)
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
        elementIdentifier = self._shield.generateElements(fieldModule, coordinates, elementIdentifier, [])
        self._endElementIdentifier = elementIdentifier
