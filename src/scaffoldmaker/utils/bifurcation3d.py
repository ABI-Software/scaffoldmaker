"""
Utility functions for generating a 3-D solid bifurcation.
"""

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
from scaffoldmaker.utils.cylindermesh import Ellipse2D, EllipseShape
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds
from opencmiss.zinc.element import Element


class BifurcationMesh:
    """
    Bifurction mesh generator.
    """

    def __init__(self, fieldModule, coordinates):
        """
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        """
        # generate the mesh
        elementsCount = [2, 2, 5]
        self._elementsCount = elementsCount

        # self.px = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.pd1 = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.pd2 = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.pd3 = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.nodeId = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 2)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 2)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 2)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 2)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 2)]

        self.elementId = [[[None] * (elementsCount[0]+1) for c in range(elementsCount[1])] for c in range(elementsCount[2])]
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

        nodeparams1 = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [[1/self._elementsCount[0], 0.0, 0.0], [1/self._elementsCount[0], 0.0, 0.0]],
                       [[0.0, 1/self._elementsCount[1], 0.0], [0.0, 1/self._elementsCount[1], 0.0]]],
                      [[0.0, 0.0, 1.4], [1.2, 0.0, 1.0], [0.0, 1.0, 1.4], [[1/self._elementsCount[0], 0.0, 0.0], [0.5*0.7071, 0.0, -0.5*0.7071]],
                       [[0.0, 1/self._elementsCount[1], 0.0], [0.0, 1/self._elementsCount[1], 0.0]]]]

        nodeparams2 = [[[0.5, 0.0, 2.2], [1.2, 0.0, 1.0], [0.5, 1.0, 2.2], [[0.0, 0.0, -1 / self._elementsCount[1]], [0.5 * 0.7071, 0.0, -0.5 * 0.7071]],
                        [[0.0, 1 / self._elementsCount[0], 0.0], [0.0, 1 / self._elementsCount[0], 0.0]]],
                       [[1.7, 0.0, 2.2], [1.7, 0.0, 1.2], [1.7, 1.0, 2.2], [[0.0, 0.0, -1 / self._elementsCount[1]], [0.0, 0.0, -1 / self._elementsCount[1]]],
                        [[0.0, 1 / self._elementsCount[0], 0.0], [0.0, 1 / self._elementsCount[0], 0.0]]]]

        baseleg1 = BaseLeg(self._elementsCount, nodeparams1)
        self.copyBaseLeg2Bifurcation(baseleg1, 1)

        self.px[self._elementsCount[2]//2+1][0][self._elementsCount[0]] = [0.0, 1.0, 2.2]
        self.px[self._elementsCount[2]//2+1][1][self._elementsCount[0]] = [0.0, 0.5, 2.2]
        self.px[self._elementsCount[2]//2+1][2][self._elementsCount[0]] = [0.0, 0.0, 2.2]
        self.pd1[self._elementsCount[2]//2+1][0][self._elementsCount[0]] = [0.0, -0.5, 0.0]
        self.pd1[self._elementsCount[2]//2+1][1][self._elementsCount[0]] = [0.0, -0.5, 0.0]
        self.pd1[self._elementsCount[2]//2+1][2][self._elementsCount[0]] = [0.0, -0.5, 0.0]
        self.pd2[self._elementsCount[2]//2+1][0][self._elementsCount[0]] = [0.5, 0.0, 0.0]
        self.pd2[self._elementsCount[2]//2+1][1][self._elementsCount[0]] = [0.5, 0.0, 0.0]
        self.pd2[self._elementsCount[2]//2+1][2][self._elementsCount[0]] = [0.5, 0.0, 0.0]
        self.pd3[self._elementsCount[2]//2+1][0][self._elementsCount[0]] = [0.0, 0.0, 0.7]
        self.pd3[self._elementsCount[2]//2+1][1][self._elementsCount[0]] = [0.0, 0.0, 0.7]
        self.pd3[self._elementsCount[2]//2+1][2][self._elementsCount[0]] = [0.0, 0.0, 0.7]

        baseleg2 = BaseLeg(self._elementsCount, nodeparams2)
        self.copyBaseLeg2Bifurcation(baseleg2, 2)

        # self.generateBaseLeg(fieldModule, coordinates, mesh, nodes)
        self.generateNodes(nodes, fieldmodule, coordinates)
        self.generateElements(mesh, fieldmodule, coordinates)



        # elementsCount = [2, 2, 5]
        # self.px = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.pd1 = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.pd2 = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.pd3 = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.nodeId = [[[None] * (elementsCount[0] + 2) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        # self.elementId = [[[None] * (elementsCount[0]+1) for c in range(elementsCount[1])] for c in range(elementsCount[2])]
        #
        # nodeparams1 = [[[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [[0.0, -1/self._elementsCount[1], 0.0], [0.0, -1/self._elementsCount[1], 0.0]],
        #                [[1/self._elementsCount[0], 0.0, 0.0], [1/self._elementsCount[0], 0.0, 0.0]]],
        #               [[0.0, 0.0, 1.4], [0.0, -1.0, 1.4], [1.2, 0.0, 1.0],[[0.0, -1/self._elementsCount[1], 0.0], [0.0, -1/self._elementsCount[1], 0.0]],
        #                [[1/self._elementsCount[0], 0.0, 0.0], [0.5*0.7071, 0.0, -0.5*0.7071]]]]
        #
        # nodeparams2 = [[[0.5, 0.0, 2.2], [0.5, -1.0, 2.2], [1.2, 0.0, 1.0],[[0.0, -1 / self._elementsCount[1], 0.0], [0.0, -1 / self._elementsCount[1], 0.0]],
        #                [[0.0, 0.0, -1 / self._elementsCount[1]], [0.5 * 0.7071, 0.0, -0.5 * 0.7071]]],
        #                [[1.7, 0.0, 2.2], [1.7, -1.0, 2.2], [1.7, 0.0, 1.2], [[0.0, -1 / self._elementsCount[0], 0.0], [0.0, -1 / self._elementsCount[0], 0.0]],
        #                 [[0.0, 0.0, -1 / self._elementsCount[1]], [0.0, 0.0, -1 / self._elementsCount[1]]],]]
        #
        #
        # baseleg1 = BaseLeg(self._elementsCount, nodeparams1)
        # self.copyBaseLeg2Bifurcation(baseleg1, 1)
        #
        # # self.px[self._elementsCount[2]//2+1][0][self._elementsCount[0] + 1] = [0.0, 1.0, 2.2]
        # # self.px[self._elementsCount[2]//2+1][1][self._elementsCount[0] + 1] = [0.0, 0.5, 2.2]
        # # self.px[self._elementsCount[2]//2+1][2][self._elementsCount[0] + 1] = [0.0, 0.0, 2.2]
        # # self.pd1[self._elementsCount[2]//2+1][0][self._elementsCount[0] + 1] = [0.0, -0.5, 0.0]
        # # self.pd1[self._elementsCount[2]//2+1][1][self._elementsCount[0] + 1] = [0.0, -0.5, 0.0]
        # # self.pd1[self._elementsCount[2]//2+1][2][self._elementsCount[0] + 1] = [0.0, -0.5, 0.0]
        # # self.pd2[self._elementsCount[2]//2+1][0][self._elementsCount[0] + 1] = [0.5, 0.0, 0.0]
        # # self.pd2[self._elementsCount[2]//2+1][1][self._elementsCount[0] + 1] = [0.5, 0.0, 0.0]
        # # self.pd2[self._elementsCount[2]//2+1][2][self._elementsCount[0] + 1] = [0.5, 0.0, 0.0]
        # # self.pd3[self._elementsCount[2]//2+1][0][self._elementsCount[0] + 1] = [0.0, 0.0, 0.7]
        # # self.pd3[self._elementsCount[2]//2+1][1][self._elementsCount[0] + 1] = [0.0, 0.0, 0.7]
        # # self.pd3[self._elementsCount[2]//2+1][2][self._elementsCount[0] + 1] = [0.0, 0.0, 0.7]
        #
        # baseleg2 = BaseLeg(self._elementsCount, nodeparams2)
        # self.copyBaseLeg2Bifurcation(baseleg2, 2)

        # self.generateBaseLeg(fieldModule, coordinates, mesh, nodes)
        # self.generateNodes(nodes, fieldmodule, coordinates)
        # self.generateElements(mesh, fieldmodule, coordinates)

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

    def generateNodes(self, nodes, fieldModule, coordinates):
        """
        Create cylinder nodes from coordinates.
        :param nodes: nodes from coordinates.
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self.topologygenerateNodes(fieldModule, coordinates, nodeIdentifier)
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
        elementIdentifier = self.topologygenerateElements(fieldModule, coordinates, elementIdentifier, [])
        self._endElementIdentifier = elementIdentifier

    def topologygenerateNodes(self, fieldmodule, coordinates, startNodeIdentifier):
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

        for n3 in range(self._elementsCount[2] + 2):
            for n2 in range(self._elementsCount[1] + 1):
                for n1 in range(self._elementsCount[0] + 1):
                    if self.px[n3][n2][n1]:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        self.nodeId[n3][n2][n1] = nodeIdentifier
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, self.px [n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, self.pd1[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, self.pd2[n3][n2][n1])
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, self.pd3[n3][n2][n1])
                        nodeIdentifier += 1

        return nodeIdentifier

    def topologygenerateElements(self, fieldmodule, coordinates, startElementIdentifier, meshGroups=[]):
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
        e1z = 2*self._elementsCount[0] - 1
        e1y = e1z - 1
        e2a = 0
        e2b = e2a + 1
        e2c = e2a + 2
        e2z = 2*self._elementsCount[1]-1
        e2y = e2z - 1
        # e2x = e2z - 2
        for e3 in range(self._elementsCount[2]):
            for e2 in range(self._elementsCount[1]):
                for e1 in range(self._elementsCount[0]):
                    eft1 = eft
                    scalefactors = None
                    if e3 >= 3:
                        e3t = e3 + 1
                    else:
                        e3t = e3
                    nids = [ self.nodeId[e3t][e2][e1], self.nodeId[e3t][e2 + 1][e1], self.nodeId[e3t+1][e2][e1], self.nodeId[e3t+1][e2 + 1][e1],
                             self.nodeId[e3t][e2][e1 + 1], self.nodeId[e3t][e2 + 1][e1 + 1], self.nodeId[e3t+1][e2][e1 + 1], self.nodeId[e3t+1][e2 + 1][e1 + 1] ]

                    if (e2 < e2b) or (e2 > e2y):
                        if (e1 < e1b) or (e1 > e1y):
                            continue  # no element due to triple point closure
                        if (e2 < e2a) or (e2 > e2z):
                            if e2 < e2a:
                                nids = [self.nodeId[e3][e2+1][e1], self.nodeId[e3][e2+1][e1+1], self.nodeId[e3+1][e2+1][e1], self.nodeId[e3+1][e2+1][e1+1],
                                        self.nodeId[e3][e2][e1], self.nodeId[e3][e2][e1+1], self.nodeId[e3+1][e2][e1],  self.nodeId[e3+1][e2][e1+1]]
                            elif e2 > e2z:
                                nids = [self.nodeId[e3][e2][e1+1], self.nodeId[e3][e2][e1], self.nodeId[e3+1][e2][e1+1], self.nodeId[e3+1][e2][e1],
                                        self.nodeId[e3][e2+1][e1+1], self.nodeId[e3][e2+1][e1], self.nodeId[e3+1][e2+1][e1+1], self.nodeId[e3+1][e2+1][e1]]
                        elif (e2 == e2a) or (e2 == e2z):
                            # bottom and top row elements
                            if e2 == e2a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e3 == self._elementsCount[2] // 2 + 1 and e1 == e1b:
                                    e3r = e3t-1
                                    nids[0] = self.nodeId[e3r][e2][e1]
                                    nids[1] = self.nodeId[e3r][e2 + 1][e1]

                                if e3 == 2:
                                    remapEftNodeValueLabel(eft1, [3, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 5], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS1, [])])
                                elif e3 == 3:
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS1, [])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS3,
                                                               [(Node.VALUE_LABEL_D_DS1, [])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 5, 7], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [])])

                                else:
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                                           [(Node.VALUE_LABEL_D_DS3, [1])])
                                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [])])

                                if (e1 == e1b) or (e1 == e1y):
                                    # map bottom triple point element
                                    if e1 == e1b:
                                        if e3 != 2:
                                            remapEftNodeValueLabel(eft1, [ 2, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                    else:
                                        remapEftNodeValueLabel(eft1, [ 6, 8 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [] ), ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            elif e2 == e2z:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [1])])
                                remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS3, [])])
                                if (e1 == e1b) or (e1 == e1y):
                                    # map top triple point element
                                    if e1 == e1b:
                                        remapEftNodeValueLabel(eft1, [1, 3], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [1])])
                                    else:
                                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1,[(Node.VALUE_LABEL_D_DS1, []),(Node.VALUE_LABEL_D_DS3, [])])

                    elif (e2 == e2b) or (e2 == e2y):
                        if (e1 <= e1a) or (e1 >= e1z):
                            if e1 < e1a:
                                e2r = e1
                                if e2 == e2b:
                                    nids = [self.nodeId[e3][e2c][e1+1], self.nodeId[e3][e2r+1][e1b], self.nodeId[e3+1][e2c][e1+1], self.nodeId[e3+1][e2r+1][e1b],
                                            self.nodeId[e3][e2c][e1], self.nodeId[e3][e2r][e1b], self.nodeId[e3+1][e2c][e1], self.nodeId[e3+1][e2r][e1b]]
                                if e2 == e2y:
                                    e2r = 2*self._elementsCount[1] - e1-1
                                    nids = [self.nodeId[e3][e2r][e1b], self.nodeId[e3][e2y][e1+1], self.nodeId[e3+1][e2r][e1b], self.nodeId[e3+1][e2y][e1+1],
                                            self.nodeId[e3][e2r+1][e1b], self.nodeId[e3][e2y][e1], self.nodeId[e3+1][e2r+1][e1b], self.nodeId[e3+1][e2y][e1]]
                            elif e1 == e1a:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                scalefactors = [-1.0]
                                if e2 == e2b:
                                    if e3 == self._elementsCount[2] // 2 + 1:
                                        e3r = e3-1  # to join upper leg with the lower leg.
                                        nids[0] = self.nodeId[e3r][e2a][e1b]
                                        nids[2] = self.nodeId[e3+1][e2a][e1b]
                                        nids[1] = self.nodeId[e3r][e2 + 1][e1]
                                        nids[4] = self.nodeId[e3r][e2][e1 + 1]
                                        nids[5] = self.nodeId[e3r][e2 + 1][e1 + 1]
                                    else:
                                        nids[0] = self.nodeId[e3][e2a][e1b]
                                        nids[2] = self.nodeId[e3+1][e2a][e1b]
                                    tripleN = [5, 7]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                elif e2 == e2y:
                                    nids[1] = self.nodeId[e3][e2z+1][e1b]
                                    nids[3] = self.nodeId[e3+1][e2z+1][e1b]
                                    tripleN = [6, 8]
                                    remapEftNodeValueLabel(eft1, tripleN, Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS1, [ ( Node.VALUE_LABEL_D_DS1, [1] ) ])
                                remapEftNodeValueLabel(eft1, [ 1, 2, 3, 4 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS3, [1] ) ])
                            elif e1 == e1z:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if e2 == e2b:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    nids[4] = self.nodeId[e3][e2a][e1z]
                                    nids[6] = self.nodeId[e3+1][e2a][e1z]
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [ -1.0 ]
                                    remapEftNodeValueLabel(eft1, [ 1, 3 ], Node.VALUE_LABEL_D_DS3, [ ( Node.VALUE_LABEL_D_DS1, [1] ), ( Node.VALUE_LABEL_D_DS3, [] ) ])
                                elif e2 == e2y:
                                    nids[5] = self.nodeId[e3][e2z+1][e1z]
                                    nids[7] = self.nodeId[e3+1][e2z+1][e1z]
                                    remapEftNodeValueLabel(eft1, [2, 4], Node.VALUE_LABEL_D_DS3,[(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                            elif e1 > e1z:
                                e2r = self._elementsCount[0] - e1
                                if e2 == e2b:
                                    nids = [self.nodeId[e3][e2r][e1z], self.nodeId[e3][e2c][e1], self.nodeId[e3+1][e2r][e1z], self.nodeId[e3+1][e2c][e1],
                                            self.nodeId[e3][e2r-1][e1z], self.nodeId[e3][e2c][e1+1], self.nodeId[e3+1][e2r-1][e1z], self.nodeId[e3+1][e2c][e1+1]]
                                elif e2 == e2y:
                                    e2r = e2z+e1-e1z
                                    nids[1] = self.nodeId[e3][e2r][e1z]
                                    nids[3] = self.nodeId[e3+1][e2r][e1z]
                                    nids[5] = self.nodeId[e3][e2r+1][e1z]
                                    nids[7] = self.nodeId[e3+1][e2r+1][e1z]
                        elif e1 == e1b:
                            if e2 == e2b:
                                eft1 = tricubichermite.createEftNoCrossDerivatives()
                                if e3 == self._elementsCount[2] // 2 + 1:
                                    e3r = e3 - 1
                                    nids[0] = self.nodeId[e3r][e2][e1]
                                    nids[1] = self.nodeId[e3r][e2 + 1][e1]
                                    remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS3,
                                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                                if e3 == 2:
                                    setEftScaleFactorIds(eft1, [1], [])
                                    scalefactors = [-1.0]
                                    remapEftNodeValueLabel(eft1, [1, 2], Node.VALUE_LABEL_D_DS2,
                                                           [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                                    remapEftNodeValueLabel(eft1, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                    remapEftNodeValueLabel(eft1, [3, 4, 7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                                    # remapEftNodeValueLabel(eft1, [3, 4], Node.VALUE_LABEL_D2_DS1DS2, [(Node.VALUE_LABEL_D_DS2, [1])])





                    else:
                        if e1 < e1a:
                            nids = [ self.nodeId[e3][e2 + 1][e1 + 1], self.nodeId[e3][e2][e1 + 1], self.nodeId[e3+1][e2 + 1][e1 + 1], self.nodeId[e3+1][e2][e1 + 1],
                                     self.nodeId[e3][e2 + 1][e1], self.nodeId[e3][e2][e1], self.nodeId[e3+1][e2 + 1][e1], self.nodeId[e3+1][e2][e1]]
                        elif e1 == e1a:
                            # map left column elements
                            eft1 = tricubichermite.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            scalefactors = [ -1.0 ]
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 2, 3, 4], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS3, [1])])

                    if e3 == 2:
                        if e1 == e1b:
                            # joining elements
                            if e2 == 0:
                                nids = [15, 17, 22, 24, 16, 18, 23, 25]
                                remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [3], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])

                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [])])
                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS3,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])

                                remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS3, [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                                remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D2_DS1DS2, [(Node.VALUE_LABEL_D_DS3, [])])





                            elif e2 == 1:
                                nids = [17, 20, 24, 26, 18, 21, 25, 27]

                    print(e3, e2, e1)
                    print(nids)

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
                    self.elementId[e3][e2][e1] = elementIdentifier
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
        elementsCount = [2, 2, 2]
        self._elementsCount = elementsCount
        self.nodeparams = nodeparams

        self.px = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd1 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd2 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.pd3 = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.nodeId = [[[None] * (elementsCount[0] + 1) for c in range(elementsCount[1] + 1)] for c in range(elementsCount[2] + 1)]
        self.elementId = [[[None] * elementsCount[0] for c in range(elementsCount[1])] for c in range(elementsCount[2])]
        self.generateBaseLeg(nodeparams)

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

        # self.generateNodes(nodes, fieldmodule, coordinates)
        # self.generateElements(mesh, fieldmodule, coordinates)

        # self.genetateBottomSurface([1.7, 0.0, 2.2], [1.7, 0.0, 1.2], [1.7, 1.0, 2.2],
        #                            [[0.0, 0.0, -1 / self._elementsCount[1]], [0.0, 0.0, -1 / self._elementsCount[1]]],
        #                            [[0.0, 1 / self._elementsCount[0], 0.0], [0.0, 1 / self._elementsCount[0], 0.0]])
        # self.generateTopSurface([0.5, 0.0, 2.2], [1.2, 0.0, 1.0], [0.5, 1.0, 2.2],
        #                         [[0.0, 0.0, -1 / self._elementsCount[1]], [0.5 * 0.7071, 0.0, -0.5 * 0.7071]],
        #                         [[0.0, 1 / self._elementsCount[0], 0.0], [0.0, 1 / self._elementsCount[0], 0.0]])
        # self.generateMiddleLevels()
        # self.smoothd2()

        # self.generateNodes(nodes, fieldmodule, coordinates)
        # self.generateElements(mesh, fieldmodule, coordinates)

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
            delta = [ellipse.px[0][self._elementsCount[0]][c] - ellipse.px[self._elementsCount[1]][0][c] for c in range(3)]
            normal = vector.normalise(vector.crossproduct3(td1c1[0], td1c2[0]))
            delnormag = vector.dotproduct(delta, normal)
            delnor = [delnormag*c for c in normal]
            if 0<n<elementsCountAround:
                ellipse.px[n2][n1] = [ellipse.px[n2][n1][c] -(1-xi) * delnor[c] for c in range(3)]
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
