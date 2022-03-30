"""
Utility functions for generating a solid sphere/spheroid (ellipsoid in general).
"""

from enum import Enum

import math

from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from opencmiss.zinc.field import Field
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.cylindermesh import Ellipse2D, EllipseShape
from scaffoldmaker.utils.shieldmesh import ShieldMesh3D


class SphereShape(Enum):
    SPHERE_SHAPE_FULL = 1
    SPHERE_SHAPE_HALF_AAP = 2    # AAP is a hemisphere where x_a3>=0
    SPHERESHIELD_SHAPE_OCTANT_PPP = 3  # positive axis1, positive axis2, positive axis3
    SPHERESHIELD_SHAPE_OCTANT_NPP = 4
    SPHERESHIELD_SHAPE_OCTANT_PNP = 5
    SPHERESHIELD_SHAPE_OCTANT_NNP = 6
    SPHERESHIELD_SHAPE_OCTANT_PPN = 7
    SPHERESHIELD_SHAPE_OCTANT_NPN = 8
    SPHERESHIELD_SHAPE_OCTANT_PNN = 9
    SPHERESHIELD_SHAPE_OCTANT_NNN = 10


class SphereMesh:
    """
    Sphere mesh generator.
    """

    def __init__(self, fieldModule, coordinates, centre, axes, elementsCountAcross,
                 elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                 sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, rangeOfRequiredElements=None,
                 boxDerivatives=None, useCrossDerivatives=False,  meshGroups=[]):
        """
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :param centre, axes: centre and axes of the sphere.
        :param elementsCountAcross: [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3]
        Total number of elements across the sphere axes.
        :param elementsCountAcrossShell, elementsCountAcrossTransition: Total number of elements across each axis
         consists of regular elements in the middle cube, transition elements from cube to a sphere (core boundary)
          and shell elements around it. Shell nodes and derivatives are similar to the core boundary and don't need
           remapping. The topology of the shield structure is extended to 3D with quadruple points.
        :param sphereShape: A value from enum SphereShape specifying the shape of sphere. Octant_PPP, for example is
         the octant in positive axis1, positive axis2 and positive axis3. SPHERE_SHAPE_HALF_AAP is a hemisphere on
          the positive side of axis3.
        :param rangeOfRequiredElements: Specifies the range of elements required to be created. It can be used to
         create the part of sphere required. If None or same as elementsCountAcross the whole part will be created.
        :param boxDerivatives: It is a list of [deriv1,deriv2,deriv3]. It is used to set the derivative directions.
        default is [1, 3, 2] which means it makes -d1, d3 and d2 in direction of axis1, axis2 and axis3. To make
        d1, d2 and d3 directions in direction of axis1, axis2 and axis3 use [-1, 2, 3].
        """

        self._axes = axes
        self._radius = [vector.magnitude(axis) for axis in axes]
        self._coreRadius = []
        self._shield = None
        self._elementsCount = [ec - 2 * elementsCountAcrossShell for ec in elementsCountAcross]
        self._elementsCount_total = elementsCountAcross
        self._elementsCountAcrossShell = elementsCountAcrossShell
        self._elementsCountAcrossTransition = elementsCountAcrossTransition
        self._elementsCountAcrossRim = self._elementsCountAcrossShell + self._elementsCountAcrossTransition - 1
        self._shellProportion = shellProportion

        self._startNodeIdentifier = 1
        self._startElementIdentifier = 1
        self._endNodeIdentifier = 1
        self._endElementIdentifier = 1
        self._sphereShape = sphereShape
        self._rangeOfRequiredElements = rangeOfRequiredElements

        self._useCrossDerivatives = useCrossDerivatives

        self._centre = centre

        self._boxDerivatives = boxDerivatives if boxDerivatives else [1, 3, 2]
        self._meshGroups = meshGroups
        self._shield3D = None

        for i in range(3):
            # Make a sphere of radius 1 first, then add the shell layers on top of it and resize it back to R=1.0
            elementsAxis = elementsCountAcross[i] - elementsCountAcrossShell * (1 - 1.0)
            self._coreRadius.append(
                (1 - 1.0 * elementsCountAcrossShell / elementsAxis) * self._radius[i])

        # generate the mesh
        self.createSphereMesh3d(fieldModule, coordinates)

    def createSphereMesh3d(self, fieldModule, coordinates):
        """
        Create a sphere mesh based on the shield topology.
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :return: Final values of nextNodeIdentifier, nextElementIdentifier.
        """
        if 'OCTANT' in str(self._sphereShape):
            minimum = [2, 2, 2]
            era = 0
        elif 'SPHERE_SHAPE_HALF_AAP' in str(self._sphereShape):
            minimum = [4, 4, 2]
            era = 2
        else:
            minimum = [4, 4, 4]
            era = 3

        for i in range(3):
            assert (self._elementsCount[i] >= minimum[i]), 'createSphereMesh3d:  Invalid number of elements'

        for i in range(era):
            assert (self._elementsCount[i] % 2 == 0), 'createSphereMesh3d: number of across elements' \
                                                          ' is not an even number'

        assert (self._sphereShape in [SphereShape.SPHERE_SHAPE_FULL,
                                      SphereShape.SPHERE_SHAPE_HALF_AAP,
                                      SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP]), \
            'createSphereMesh3d: Invalid sphere mode.'

        elementsCountRim = self._elementsCountAcrossRim

        self._shield3D = ShieldMesh3D(self._elementsCount, elementsCountRim, box_derivatives=self._boxDerivatives)

        nodes = fieldModule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldModule.findMeshByDimension(3)

        elementsCountAcrossTransition = self._elementsCountAcrossTransition

        # order of octants is important because common nodes will overwrite
        OctantVariationsAll = [SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPN, SphereShape.SPHERESHIELD_SHAPE_OCTANT_PNN,
                               SphereShape.SPHERESHIELD_SHAPE_OCTANT_NNN, SphereShape.SPHERESHIELD_SHAPE_OCTANT_NPN,
                               SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, SphereShape.SPHERESHIELD_SHAPE_OCTANT_PNP,
                               SphereShape.SPHERESHIELD_SHAPE_OCTANT_NNP, SphereShape.SPHERESHIELD_SHAPE_OCTANT_NPP]

        if self._sphereShape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP:
            OctantVariationsType = [SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP]
        elif self._sphereShape == SphereShape.SPHERE_SHAPE_HALF_AAP:
            OctantVariationsType = OctantVariationsAll[4:9]
        else:
            OctantVariationsType = OctantVariationsAll

        # Make a sphere of radius 1 first, then add the shell layers on top of it and resize it back to R=1.0
        for octantType in OctantVariationsType:
            axes, elementsCountAcross, boxDerivatives = self.get_octant_axes_and_elements_count(octantType)
            octant = OctantMesh(self._centre, axes, elementsCountAcross,
                                0, elementsCountAcrossTransition, 1.0,
                                sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP,
                                useCrossDerivatives=False, boxDerivatives=boxDerivatives)
            self.copy_octant_nodes_to_sphere_shield(octant, octantType)

        self.modify_octant_common_nodes()
        if self._elementsCountAcrossShell > 0:
            self.add_shell_layers()
        self.sphere_to_spheroid()

        self.generateNodes(nodes, fieldModule, coordinates)
        self.generateElements(mesh, fieldModule, coordinates)

    def get_octant_axes_and_elements_count(self, octant_shape):
        """
        Get the octants axes, elements count across and box mapping for 8 octants in 8 different regions.
        :param octant_shape: A value from enum SphereShape specifying the shape of sphere
        :return: axes, elementsCountAcross, boxDerivatives for the octant.
        """
        boxDerivativesV1 = self._boxDerivatives
        boxDerivativesV2 = [self._boxDerivatives[1], self._boxDerivatives[0], self._boxDerivatives[2]]
        axesV1 = [self._axes[0], self._axes[1], self._axes[2]]
        axesV2 = [self._axes[1], self._axes[0], self._axes[2]]
        elementsCountAcrossV1 = [c for c in self._elementsCount]
        elementsCountAcrossV2 = [self._elementsCount[1], self._elementsCount[0], self._elementsCount[2]]

        if octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP:
            octant_number = 1
        elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PNP:
            octant_number = 2
        elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NNP:
            octant_number = 3
        elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NPP:
            octant_number = 4
        elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPN:
            octant_number = 5
        elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PNN:
            octant_number = 6
        elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NNN:
            octant_number = 7
        elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NPN:
            octant_number = 8
        else:
            raise ValueError("Not implemented.")

        axesSignsV1 = [1 if str(octant_shape).split('_')[-1][c] == 'P' else -1 for c in range(3)]
        axesSignsV2 = [axesSignsV1[1], axesSignsV1[0], axesSignsV1[2]]

        signs = self._shield3D.get_octant_signs(octant_number)
        if octant_number in [2, 4, 5, 7]:
            boxDerivatives = [boxDerivativesV2[c] * signs[c] for c in range(3)]
            elementsCountAcross = elementsCountAcrossV2
            axes = [vector.scaleVector(axesV2[c], axesSignsV2[c]) for c in range(3)]
        else:
            boxDerivatives = [boxDerivativesV1[c] * signs[c] for c in range(3)]
            elementsCountAcross = elementsCountAcrossV1
            axes = [vector.scaleVector(axesV1[c], axesSignsV1[c]) for c in range(3)]

        if self._sphereShape == SphereShape.SPHERE_SHAPE_HALF_AAP:
            elementsCountAcross[0] = elementsCountAcross[0]//2
            elementsCountAcross[1] = elementsCountAcross[1]//2
        elif self._sphereShape == SphereShape.SPHERE_SHAPE_FULL:
            elementsCountAcross[0] = elementsCountAcross[0] // 2
            elementsCountAcross[1] = elementsCountAcross[1] // 2
            elementsCountAcross[2] = elementsCountAcross[2] // 2

        axes = [vector.normalise(v) for v in axes]
        return axes, elementsCountAcross, boxDerivatives

    def copy_octant_nodes_to_sphere_shield(self, octant, octant_shape):
        """
        Copy octant nodes to sphere shield.
        :param octant: octant instance.
        :param octant_shape: A value from enum SphereShape specifying the shape of sphere
        :return:
        """
        n3a = 0
        n3z = self._elementsCount[2]
        n2a = 0
        n2z = self._elementsCount[0]
        n1a = 0
        n1z = self._elementsCount[1]
        if self._sphereShape == SphereShape.SPHERE_SHAPE_FULL:
            if octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP:
                octant_number = 1
                n3a = self._elementsCount[2]//2
                n2z = self._elementsCount[0]//2
                n1a = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PNP:
                octant_number = 2
                n3a = self._elementsCount[2]//2
                n2z = self._elementsCount[0]//2
                n1z = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NNP:
                octant_number = 3
                n3a = self._elementsCount[2]//2
                n2a = self._elementsCount[0]//2
                n1z = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NPP:
                octant_number = 4
                n3a = self._elementsCount[2]//2
                n2a = self._elementsCount[0]//2
                n1a = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPN:
                octant_number = 5
                n3z = self._elementsCount[2]//2
                n2z = self._elementsCount[0]//2
                n1a = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PNN:
                octant_number = 6
                n3z = self._elementsCount[2]//2
                n2z = self._elementsCount[0]//2
                n1z = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NNN:
                octant_number = 7
                n3z = self._elementsCount[2]//2
                n2a = self._elementsCount[0]//2
                n1z = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NPN:
                octant_number = 8
                n3z = self._elementsCount[2]//2
                n2a = self._elementsCount[0]//2
                n1a = self._elementsCount[1]//2
            else:
                raise ValueError("Not implemented.")
        elif self._sphereShape == SphereShape.SPHERE_SHAPE_HALF_AAP:
            if octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP:
                octant_number = 1
                n2z = self._elementsCount[0]//2
                n1a = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PNP:
                octant_number = 2
                n2z = self._elementsCount[0]//2
                n1z = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NNP:
                octant_number = 3
                n2a = self._elementsCount[0]//2
                n1z = self._elementsCount[1]//2
            elif octant_shape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_NPP:
                octant_number = 4
                n2a = self._elementsCount[0]//2
                n1a = self._elementsCount[1]//2
        elif self._sphereShape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP:
            octant_number = 1
        else:
            raise ValueError("Not implemented.")

        shield3D = octant.get_shield()
        for n3 in range(n3a, n3z + 1):
            for n2 in range(n2a, n2z + 1):
                for n1 in range(n1a, n1z + 1):
                    n3o, n2o, n1o = self._shield3D.get_octant_node_index(octant_number, n3, n2, n1)
                    self._shield3D.px[n3][n2][n1] = shield3D.px[n3o][n2o][n1o]
                    self._shield3D.pd1[n3][n2][n1] = shield3D.pd1[n3o][n2o][n1o]
                    self._shield3D.pd2[n3][n2][n1] = shield3D.pd2[n3o][n2o][n1o]
                    self._shield3D.pd3[n3][n2][n1] = shield3D.pd3[n3o][n2o][n1o]

    def modify_octant_common_nodes(self):
        """

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

        # modify pole on highest z.
        if self._sphereShape == SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP:
            btd1[n3z][n2z][0] = [-c for c in btd1[n3z][n2z][0]]
            btd2[n3z][n2z][0] = [-c for c in btd2[n3z][n2z][0]]
        elif self._sphereShape == SphereShape.SPHERE_SHAPE_HALF_AAP:
            temp = btd1[n3z][n2z//2][n1z//2]
            btd1[n3z][n2z//2][n1z//2] = btd2[n3z][n2z//2][n1z//2]
            btd2[n3z][n2z//2][n1z//2] = [-c for c in temp]
        elif self._sphereShape == SphereShape.SPHERE_SHAPE_FULL:
            temp = btd1[n3z][n2z//2][n1z//2]
            btd1[n3z][n2z//2][n1z//2] = btd2[n3z][n2z//2][n1z//2]
            btd2[n3z][n2z//2][n1z//2] = [-c for c in temp]

    def add_shell_layers(self):
        """
        Add shell layers on top of the sphere. Calculates shell element thickness using shell proportion and adds
        shell elements. Note that Radius without the shell elements is 1.0. R= 1.0 + thickness X number of shell layers.
        :return:
        """

        def on_sphere(n3, n2, n1):
            """
            Check if the given point is on the sphere.
            :param n3, n2, n1: node indexes in data structure [n3][n2][n1]
            :return: True, if it is on the sphere.
            """
            n3z = self._elementsCount[2]
            n2z = self._elementsCount[0]
            n1z = self._elementsCount[1]

            x = self._shield3D.px

            return (n3 == n3z or n3 == 0 or n2 == 0 or n2 == n2z or n1 == 0 or n1 == n1z) and x[n3][n2][n1]

        shield3D_with_shell = ShieldMesh3D(self._elementsCount_total, self._elementsCountAcrossShell,
                                           box_derivatives=self._boxDerivatives)

        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        shell_proportion = self._shellProportion
        thickness = [2*1.0/ne * shell_proportion for ne in self._elementsCount]
        elementsCountShell = self._elementsCountAcrossShell

        for n3 in range(self._elementsCount[2] + 1):
            for n2 in range(self._elementsCount[0] + 1):
                for n1 in range(self._elementsCount[1] + 1):
                    # new numbering due to introducing shell layers
                    n3n, n2n, n1n = n3 + elementsCountShell, n2 + elementsCountShell, n1 + elementsCountShell
                    # generate rim/shell nodes. Note, sphere radius without shell is 1.0. To keep it 1.0 after adding
                    # shell elements, we need to shrink it back. See ratio in sphere_to_spheroid.
                    if on_sphere(n3, n2, n1):
                        for n in range(elementsCountShell):
                            n3s, n2s, n1s = n3n, n2n, n1n
                            nl = n + 1
                            if n3 == 0:
                                n3s = elementsCountShell - nl
                            if n2 == 0:
                                n2s = elementsCountShell - nl
                            if n1 == 0:
                                n1s = elementsCountShell - nl
                            if n3 == self._elementsCount[2]:
                                n3s = n3n + nl
                            if n2 == self._elementsCount[0]:
                                n2s = n2n + nl
                            if n1 == self._elementsCount[1]:
                                n1s = n1n + nl

                            shield3D_with_shell.px[n3s][n2s][n1s] =\
                                [(1 + nl * thickness[c]) * btx[n3][n2][n1][c] for c in range(3)]
                            shield3D_with_shell.pd1[n3s][n2s][n1s] = \
                                [(1 + nl * thickness[c]) * btd1[n3][n2][n1][c] for c in range(3)]
                            shield3D_with_shell.pd2[n3s][n2s][n1s] =\
                                [(1 + nl * thickness[c]) * btd2[n3][n2][n1][c] for c in range(3)]
                            shield3D_with_shell.pd3[n3s][n2s][n1s] =\
                                [shell_proportion * btd3[n3][n2][n1][c] for c in range(3)]
                    # Add the rest of the nodes.
                    if btx[n3][n2][n1]:
                        shield3D_with_shell.px[n3n][n2n][n1n] = btx[n3][n2][n1]
                        shield3D_with_shell.pd1[n3n][n2n][n1n] = btd1[n3][n2][n1]
                        shield3D_with_shell.pd2[n3n][n2n][n1n] = btd2[n3][n2][n1]
                        shield3D_with_shell.pd3[n3n][n2n][n1n] = btd3[n3][n2][n1]

        self._shield3D = shield3D_with_shell

    def sphere_to_spheroid(self):
        """
        Using the radius in each direction,transform the sphere to ellipsoid.
        :return:
        """
        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        ratio = [1/(1 + 2/ne * self._shellProportion * self._elementsCountAcrossShell) for ne in self._elementsCount]

        for n3 in range(self._elementsCount_total[2] + 1):
            for n2 in range(self._elementsCount_total[0] + 1):
                for n1 in range(self._elementsCount_total[1] + 1):
                    if btx[n3][n2][n1]:
                        btx[n3][n2][n1] = [ratio[c] * self._radius[c] * btx[n3][n2][n1][c] for c in range(3)]
                        btd1[n3][n2][n1] = [ratio[c] * self._radius[c] * btd1[n3][n2][n1][c] for c in range(3)]
                        btd2[n3][n2][n1] = [ratio[c] * self._radius[c] * btd2[n3][n2][n1][c] for c in range(3)]
                        btd3[n3][n2][n1] = [ratio[c] * self._radius[c] * btd3[n3][n2][n1][c] for c in range(3)]

    def generateNodes(self, nodes, fieldModule, coordinates):
        """
        Create cylinder nodes from coordinates.
        :param nodes: nodes from coordinates.
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self._shield3D.generateNodes(fieldModule, coordinates, nodeIdentifier,
                                                      self._rangeOfRequiredElements)
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
        elementIdentifier = self._shield3D.generateElements(fieldModule, coordinates, elementIdentifier,
                                                            self._rangeOfRequiredElements, self._meshGroups)
        self._endElementIdentifier = elementIdentifier


class OctantMesh:
    """
    Octant mesh generator.
    """

    def __init__(self, centre, axes, elementsCountAcross,
                 elementsCountAcrossShell, elementsCountAcrossTransition, shellProportion,
                 sphereShape=SphereShape.SPHERESHIELD_SHAPE_OCTANT_PPP, useCrossDerivatives=False, boxDerivatives=None):
        """
        :param centre, axes: centre and axes of the octant.
        :param elementsCountAcross: [elementsCountAcrossAxis1, elementsCountAcrossAxis2, elementsCountAcrossAxis3] Total
         number of elements across the octant axes.
        :param elementsCountAcrossShell, elementsCountAcrossTransition: Total number of elements across each axis
         consists of regular elements in the middle cube, transition elements from cube to a sphere (core boundary)
          and shell elements around it. Shell nodes and derivatives are similar to the core boundary and don't need
           remapping. The topology of the shield structure is extended to 3D with a quadruple points.
        :param sphereShape: A value from enum SphereShape specifying one of the 8 octant regions. Octant_PPP for
         example, is the octant in positive axis1, positive axis2 and positive axis3.
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

        self._sphereShape = sphereShape

        self._useCrossDerivatives = useCrossDerivatives

        self._centre = centre

        self._boxDerivatives = boxDerivatives if boxDerivatives else [1, 3, 2]
        self._shield3D = None

        for i in range(3):
            elementsAxis = elementsCountAcross[i] - elementsCountAcrossShell * (1 - 1.0)
            self._coreRadius.append(
                (1 - 1.0 * elementsCountAcrossShell / elementsAxis) * self._radius[i])

        # generate the mesh
        self.createOctantMesh3d()

    def createOctantMesh3d(self):
        """
        Create an octant mesh based on the shield topology.
        """
        elementsCountRim = self._elementsCountAcrossRim

        self._shield3D = ShieldMesh3D(self._elementsCount, elementsCountRim)

        self.create_boundary_ellipses_nodes()
        self.create_surface_and_interior_nodes()
        self._shield3D.set_derivatives(self._boxDerivatives)

    def create_boundary_ellipses_nodes(self):
        """
        Use ellipse class to create 3 ellipses needed for the octant.
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

        # set derivatives in up and right direction for each ellipse in square region and the circle around.
        # We need to do this so the octant derivatives comes out as we want.
        squareDerivatives = [[1, 3], [2, 3], [1, -2]]
        circleDerivatives = [[1, 3], [2, 3], [-2, 3]]

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

            self.copy_ellipses_nodes_to_shield_nodes(ellipse, i+1, squareDerivatives[i], circleDerivatives[i])

    def copy_ellipses_nodes_to_shield_nodes(self, ellipse, ellipse_number, squareDerivatives, circleDerivatives):
        """
        Copy coordinates and derivatives of ellipse to shield.
        """
        btx = self._shield3D.px
        btd1 = self._shield3D.pd1
        btd2 = self._shield3D.pd2
        btd3 = self._shield3D.pd3

        shield = ellipse.getShield()

        # Modify the shield to get only the quarter that you want. make others None so you can generate the nodes.
        shield.remap_derivatives(squareDerivatives, circleMapping=circleDerivatives)
        for n3 in range(self._elementsCount[2] + 1):
            for n2 in range(self._elementsCount[0] + 1):
                for n1 in range(self._elementsCount[1] + 1):
                    n3s, n2s, n1s = n3, n2, n1
                    n3e, n2e, n1e = 0, n2, n1
                    if ellipse_number == 1:
                        if n3 > 0:
                            continue
                        if n2 == 0 and n1 == self._elementsCount[1] - 1:
                            n1s = n1 + 1
                        n1e = n1 + self._elementsCount[1]
                    elif ellipse_number == 2:
                        if n2 < self._elementsCount[0]:
                            continue
                        if n3 == self._elementsCount[2] and n1 == self._elementsCount[1] - 1:
                            n1s = n1 + 1
                        n2e = n3 + self._elementsCount[2]
                        n1e = n1 + self._elementsCount[1]
                        if n3 == 0:
                            if n1 == self._elementsCount[1]:
                                btd2[n3][n2s][n1s] = shield.pd2[0][n2e][n1e]
                            else:
                                btd2[n3][n2s][n1s] = shield.pd2[0][n2e][n1e]
                            continue
                    elif ellipse_number == 3:
                        if n1 > 0:
                            continue
                        if n2 == 0 and n3 == self._elementsCount[2] - 1:
                            n3s = self._elementsCount[2]
                        n1e = self._elementsCount[2] - n3
                        if n3 == 0:
                            if n2 == 0:
                                btd2[n3][n2][n1s] = shield.pd2[0][n2][n1e]
                            elif 0 < n2 < self._elementsCount[0]:
                                btd2[n3][n2][n1s] = shield.pd2[0][n2][n1e]
                            continue
                        else:
                            if n2 == self._elementsCount[0]:
                                if n3 == self._elementsCount[2]:
                                    btd1[n3][n2][n1s] = shield.pd2[0][n2][n1e]
                                else:
                                    btd1[n3][n2][n1s] = shield.pd1[0][n2][n1e]
                                continue

                    if shield.px[n3e][n2e][n1e]:
                        btx[n3s][n2s][n1s] = shield.px[n3e][n2e][n1e]
                        btd1[n3s][n2s][n1s] = shield.pd1[n3e][n2e][n1e]
                        btd2[n3s][n2s][n1s] = shield.pd2[n3e][n2e][n1e]
                        btd3[n3s][n2s][n1s] = shield.pd3[n3e][n2e][n1e]

    def create_surface_and_interior_nodes(self):
        """
        Create octant exterior surface nodes and interior nodes
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

        elementsAroundEllipse12 = self._elementsCount[0] + self._elementsCount[1] - 2
        radiansAroundEllipse12 = math.pi / 2
        radiansPerElementAroundEllipse12 = radiansAroundEllipse12 / elementsAroundEllipse12
        elementsAroundEllipse13 = self._elementsCount[0] + self._elementsCount[2] - 2
        radiansAroundEllipse13 = math.pi / 2
        radiansPerElementAroundEllipse13 = radiansAroundEllipse13 / elementsAroundEllipse13

        theta_2 = n3y * radiansPerElementAroundEllipse13
        theta_3 = n1y * radiansPerElementAroundEllipse12
        phi_3 = calculate_azimuth(theta_3, theta_2)
        # ratio = -0.1 * (min(self._elementsCount) - 2) + 1 if self._elementsCount[0] <= 2 else 0.2
        ratio = 1
        # local_x = intersection_of_two_great_circles_on_sphere(btx[0][0][n1y-1],
        # btx[n3z][n2z][n1z], btx[0][2][n1z], btx[n3z][0][0])
        local_x = spherical_to_cartesian(1.0, theta_3, ratio * phi_3 + (1-ratio)*math.pi/2)

        x = local_to_global_coordinates(local_x, self._axes, self._centre)

        a1, a2, a3 = local_orthogonal_unit_vectors(x, self._axes[2], self._centre)
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
        # sample on curve 1 of the triple curves and smooth the end derivatives.
        n3r1, n2r1, n1r1 = self.get_triple_curves_end_node_parameters(1, index_output=True)
        n3r2, n2r2, n1r2 = self.get_triple_curves_end_node_parameters(1, cx=1, index_output=True)
        self.sample_curves_between_two_nodes_on_sphere([n3r1, n2r1, n1r1], [n3r2, n2r2, n1r2],
                                                       self._elementsCount[0] - 1, [1, None], [1], [1])
        # curve 2
        n3r1, n2r1, n1r1 = self.get_triple_curves_end_node_parameters(1, cx=2, index_output=True)
        n3r2, n2r2, n1r2 = self.get_triple_curves_end_node_parameters(1, index_output=True)
        self.sample_curves_between_two_nodes_on_sphere([n3r1, n2r1, n1r1], [n3r2, n2r2, n1r2],
                                                       self._elementsCount[1] - 1, [1], [-2], [None, -2])
        # curve 3.
        n3r1, n2r1, n1r1 = self.get_triple_curves_end_node_parameters(1, cx=3, index_output=True)
        n3r2, n2r2, n1r2 = self.get_triple_curves_end_node_parameters(1, index_output=True)
        self.sample_curves_between_two_nodes_on_sphere([n3r1, n2r1, n1r1], [n3r2, n2r2, n1r2],
                                                       self._elementsCount[2] - 1, [2], [2], [None, None])

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
                                                           [1], [-2], [None, -2])

        # regular curves crossing curve 2
        for n1 in range(1, self._elementsCount[1] - 1):
            # bottom left. Top is done before.
            self.sample_curves_between_two_nodes_on_sphere([0, 0, n1], [n3z, 0, n1], self._elementsCount[2] - 1,
                                                           [2], [2], [None, None])

        # smooth regular curves crossing curve 1
        for n2 in range(n2a + 1, n2z):
            self.smooth_derivatives_regular_surface_curve(2, n2, [[1], [None, None]], [[-2], [-2]], [[None, -2], [-2]])

        # smooth regular curves crossing curve 2
        for n1 in range(1, n1z - 1):
            self.smooth_derivatives_regular_surface_curve(1, n1, [[2], [None, None]], [[2], [1]], [[None, 1], [1]])

        # smooth regular curves crossing curve 3
        for n3 in range(1, self._elementsCount[2] - 1):
            self.smooth_derivatives_regular_surface_curve(3, n3, [[1], [None, None]], [[1], [1]], [[None, 1], [1]])

    def create_interior_nodes(self):
        """
        Create box nodes.
        :return:
        """
        self.calculate_interior_quadruple_point()
        self.sample_interior_curves()
        self.smooth_regular_interior_curves()

    def smooth_derivatives_to_surface(self):
        """
         Smooth derivatives leading to quadruple point where 3 hex elements merge.
        :return:
        """
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
                            [(co[0]*btd1[n3r][n2r][n1r][c] + co[1]*btd2[n3r][n2r][n1r][c] +
                              co[2]*btd3[n3r][n2r][n1r][c]) for c in range(3)])

                        tx.append(btx[n3][n2][n1])
                        td3.append(btd3[n3][n2][n1])

                        td3 = smoothCubicHermiteDerivativesLine(tx, td3, fixStartDirection=True, fixEndDirection=True)
                        btd3[n3][n2][n1] = td3[1]

    def sample_curves_between_two_nodes_on_sphere(self, id1, id2, elementsOut, dStart, dbetween, dEnd):
        """
        Samples curves on the sphere surface between two points given by their indexes.
        :param id1, id2: [n3,n2,n1] for the first and second points.
        :param dStart, dBetween, dEnd: Specifies the derivatives that are used for this curve at the beginning, end and
         in between. e.g. dStart=[2, -1, None] means d2 for the first node, -d1 for the second node and skip the third
          one.
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

        nx, nd1 = sample_curves_on_sphere(btx[id1[0]][id1[1]][id1[2]], btx[id2[0]][id2[1]][id2[2]], self._centre,
                                          elementsOut)

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

                a1, a2, a3 = local_orthogonal_unit_vectors(nx[nit], self._axes[2], self._centre)
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
        Calculate coordinates and derivatives of the quadruple point inside, where 3 hex elements merge.
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
        ra = vector.addVectors([btx[n3z][0][n1z], self._centre], [1, -1])
        radius = vector.magnitude(ra)
        local_x = vector.scaleVector(ra, (1 - ts/radius))
        x = vector.addVectors([local_x, self._centre], [1, 1])
        n3r0, n2r0, n1r0 = self.get_triple_curves_end_node_parameters(0, index_output=True)
        n3r1, n2r1, n1r1 = self.get_triple_curves_end_node_parameters(1, index_output=True)
        btx[n3r0][n2r0][n1r0] = x
        btd1[n3r0][n2r0][n1r0] = [-(btx[n3r1][n2r1][n1r1][0] - btx[n3r0][n2r0][n1r0][0]), 0.0, 0.0]
        btd2[n3r0][n2r0][n1r0] = [0.0, 0.0, (btx[n3r1][n2r1][n1r1][2] - btx[n3r0][n2r0][n1r0][2])]
        btd3[n3r0][n2r0][n1r0] = [0.0, (btx[n3r1][n2r1][n1r1][1] - btx[n3r0][n2r0][n1r0][1]), 0.0]

    def sample_interior_curves(self):
        """
        Sample box curves.
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

        tx, td1 = sampleCubicHermiteCurves([btx[n3y][1][n1y], btx[n3y][n2z][n1y]],
                                           [btd1[0][1][n1y], btd1[0][n2z][n1y]], self._elementsCount[0] - 1)[:2]

        for n2 in range(2, self._elementsCount[0]):
            btx[n3y][n2][n1y] = tx[n2-1]
            btd1[n3y][n2][n1y] = td1[n2-1]
            btd2[n3y][n2][n1y] = [0.0, 0.0, (btx[n3z][n2][n1z][2] - btx[n3y][n2][n1y][2])]
            btd3[n3y][n2][n1y] = [0.0, (btx[n3z][n2][n1z][1] - btx[n3y][n2][n1y][1]), 0.0]

        # curve 2 and parallel curves
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
        Smooth box curves.
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
        Find the indexes or node parameters for the 6 end nodes of unique curves of triple curves on the surface and
        inside.
        if cx is not given, it returns the quadruple points identified by rx.
        :param cx: curve index. Curve 1 connects quadruples to ellipse 23.
         Similarly, curve 2, connects quadruples to ellipse 13
        :param rx: 1 means on the exterior surface, 0 in one below.
        :return: indexes of the point or the node parameters.
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

    def get_shield(self):
        return self._shield3D


def calculate_azimuth(theta, theta_p):
    """
    Given polar angles of a point on the sphere surfaces, calculate the azimuth angle.
    :param theta: polar angle. In orthonormal coordinate system (axis1, axis2, axis3) with right-hand rule,
    theta is angle between common axis and point projection on plane of theta. In case theta=theta_3 and
    theta_p = theta_1, theta is between axis2 and projection
    :param theta_p: polar angle wrt other direction.
    :return: Azimuth angle.
    """
    return math.atan(1/(math.tan(theta_p)*math.cos(theta)))


def local_orthogonal_unit_vectors(x, axis3, origin):
    """
    Find local orthogonal unit vectors for a point on a sphere
    :param x: coordinates of the point.
    :param axis3: The third axis in Cartesian coordinate system (axis1, axis2, axis3)
    :return: e1, e2, e3. Unit vectors. e3 is normal to the boundary, e2 is in (e3, axis3) plane and e1 normal to them.
    """
    r = vector.addVectors([x, origin], [1, -1])
    e3 = vector.normalise(r)
    e2 = vector.vectorRejection(axis3, e3)
    e2 = vector.normalise(e2)
    e1 = vector.crossproduct3(e2, e3)

    return e1, e2, e3


def calculate_arc_length(x1, x2, origin):
    """
    Calculate the arc length between points x1 and x2.
    :param x1, x2: points coordinates.
    :return: arc length
    """
    r1 = vector.addVectors([x1, origin], [1, -1])
    r2 = vector.addVectors([x2, origin], [1, -1])
    radius = vector.magnitude(r1)
    angle = vector.angleBetweenVectors(r1, r2)
    return radius * angle


def sample_curves_on_sphere(x1, x2, origin, elementsOut):
    """

    :param x1, x2: points coordinates.
    :param elementsOut:
    :return:
    """
    r1 = vector.addVectors([x1, origin], [1, -1])
    r2 = vector.addVectors([x2, origin], [1, -1])
    deltax = vector.addVectors([r1, r2], [-1, 1])
    normal = vector.crossproduct3(r1, deltax)
    angle = vector.angleBetweenVectors(r1, r2)
    anglePerElement = angle/elementsOut
    arcLengthPerElement = calculate_arc_length(x1, x2, origin)/elementsOut

    nx = []
    nd1 = []
    for n1 in range(elementsOut + 1):
        radiansAcross = n1 * anglePerElement
        r = vector.rotateVectorAroundVector(r1, normal, radiansAcross)
        x = vector.addVectors([r, origin], [1, 1])
        d1 = vector.setMagnitude(vector.crossproduct3(normal, r), arcLengthPerElement)
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
    :param local_origin: Origin of local coordinates system specified as a list of 3 components wrt global coordinates
     system.
    :param local_axes: Axes of local coordinates system, specified as a list of list 3X3 with respect to global
     coordinates system.
    :return: Global coordinate system.
    """
    if local_origin is None:
        local_origin = [0.0, 0.0, 0.0]
    normalised_axes = [vector.normalise(v) for v in local_axes]
    return vector.addVectors([vector.addVectors(normalised_axes, local_x), local_origin])


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
