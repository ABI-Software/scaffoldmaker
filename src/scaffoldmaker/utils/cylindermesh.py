"""
Utility functions for generating a generalised 3-D solid cylinder (extruded ellipse/circle). The radii and ellipses
 orientation can be controlled using a central path subscaffold. It can be used to generate
a solid truncated cone. It also can be used for transition from a 2D base to another base (e.g., ellipse to a circle).
"""

import math
from enum import Enum

from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.meshtype_1d_path1 import extractPathParametersFromRegion
from scaffoldmaker.utils import vector, geometry
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite, \
    smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.mirror import Mirror
from scaffoldmaker.utils.shieldmesh import ShieldMesh2D, ShieldShape2D, ShieldRimDerivativeMode


class CylinderShape(Enum):
    CYLINDER_SHAPE_FULL = 1  # full cylinder is created
    CYLINDER_SHAPE_LOWER_HALF = 2  # lower half cylinder


class EllipseShape(Enum):
    Ellipse_SHAPE_FULL = 1  # full ellipse is created
    Ellipse_SHAPE_LOWER_HALF = 2  # lower half ellipse


class CylinderType(Enum):
    CYLINDER_STRAIGHT = 1  # regular cylinder
    CYLINDER_TAPERED = 2  # cylinder radius changes along the cylinder axis


class ConeBaseProgression(Enum):
    GEOMETRIC_PROGRESSION = 1  # geometric sequence decrease for major radius of bases
    ARITHMETIC_PROGRESSION = 2  # arithmetic sequence decrease for major radius of bases


class CylinderEnds:
    """
    Stores base ellipse parameters.
    """

    def __init__(self, elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell=0,
                 elementsCountAcrossTransition=1, shellProportion=1.0,
                 centre=None, alongAxis=None, majorAxis=None, minorRadius=None):
        """
        :param elementsCountAcrossMajor: Number of elements across major axis. Must be at least 2 + elementsCountRim for
         half and 4 + elementsCountRim for full cylinder.
        :param elementsCountAcrossMinor: Number of elements across minor axis.
        :param elementsCountAcrossShell: Number of elements across shell.
        :param elementsCountAcrossTransition: Number of elements between core boundary and inner square.
        :param shellProportion: Ratio of thickness of each layer in shell wrt thickness of each layer in core.
        :param centre: Centre of the ellipse.
        :param alongAxis: The cylinder axis that the base is extruded along.
        :param majorAxis: The major axis of the base. Should be perpendicular to alongAxis
        :param minorRadius: The minor radius of the ellipse.
        """
        self._centre = centre
        self._alongAxis = alongAxis
        self._majorAxis = majorAxis
        self._minorRadius = minorRadius
        if alongAxis:
            self._minorAxis = vector.setMagnitude(vector.crossproduct3(alongAxis, majorAxis), minorRadius)
        self._elementsCountAcrossMinor = elementsCountAcrossMinor
        self._elementsCountAcrossMajor = elementsCountAcrossMajor
        self._elementsCountAcrossShell = elementsCountAcrossShell
        self._elementsCountAcrossTransition = elementsCountAcrossTransition
        self._shellProportion = shellProportion
        self._majorRadius = vector.magnitude(majorAxis)
        self.px = None
        self.pd1 = None
        self.pd2 = None
        self.pd3 = None


class Tapered:
    """
    Stores parameters for making a tapered cylinder.
    """

    def __init__(self, majorRatio=1.0, majorProgressionMode=ConeBaseProgression.GEOMETRIC_PROGRESSION,
                 minorRatio=1.0, minorProgressionMode=ConeBaseProgression.GEOMETRIC_PROGRESSION):
        """
        :param ratio: radius common ratio increment along cylinder axis.
        :param progressionMode: controls the change in radius along cylinder axis. r_n+1 = r_n+ratio for arithmetic,
        r_n+1 = r_n*ratio for geometric progression.
        """
        self.majorRatio = majorRatio
        self.majorProgressionMode = majorProgressionMode
        self.minorRatio = minorRatio
        self.minorProgressionMode = minorProgressionMode


class CylinderCentralPath:
    """
    Stores ellipses parameters a long the central path.
    """

    def __init__(self, region, centralPath, elementsCount):
        """
        :param region: Zinc region to define model in.
        :param centralPath: Central path subscaffold comes from meshtype_1d_path1 and used to calculate ellipse radii.
        :param elementsCount: Number of elements needs to be sampled along the central path.
        """
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd3, cd12, cd13 = extractPathParametersFromRegion(tmpRegion,
                                                                        [Node.VALUE_LABEL_VALUE,
                                                                         Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                                                                         Node.VALUE_LABEL_D_DS3,
                                                                         Node.VALUE_LABEL_D2_DS1DS2,
                                                                         Node.VALUE_LABEL_D2_DS1DS3])
        del tmpRegion
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i], ',', cd12[i], ',', cd3[i], ',', cd13[i], '],')

        sx, sd1, se, sxi, ssf = sampleCubicHermiteCurves(cx, cd1, elementsCount)
        sd2, sd12 = interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        sd3, sd13 = interpolateSampleCubicHermite(cd3, cd13, se, sxi, ssf)

        self.centres = sx
        self.majorRadii = [vector.magnitude(a) for a in sd2]
        self.majorAxis = sd2
        self.minorRadii = [vector.magnitude(a) for a in sd3]
        self.minorAxis = sd3
        self.alongAxis = sd1


class CylinderMesh:
    """
    Cylinder mesh generator. Extrudes an ellipse/circle.
    """

    def __init__(self, fieldModule, coordinates, elementsCountAlong, base=None, end=None,
                 cylinderShape=CylinderShape.CYLINDER_SHAPE_FULL, rangeOfRequiredElements=None,
                 tapered=None, cylinderCentralPath=None, useCrossDerivatives=False , meshGroupsElementsAlong=[], meshGroups=[]):
        """
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :param base: Cylinder base ellipse. It is an instance of class CylinderEnds.
        :param end: Cylinder end ellipse. It is an instance of class CylinderEnds.
        :param elementsCountAlong: Number of elements along the cylinder axis.
        :param cylinderShape: A value from enum CylinderMode specifying.
        :param rangeOfRequiredElements: Specifies the range of elements required to be created. It can be used to
         create the part of cylinder required. If None or same as elementsCountAcross the whole part will be created.
        """

        self._centres = None
        self._majorAxis = None
        self._minorAxis = None
        self._majorRadii = None
        self._minorRadii = None
        self._coreMajorRadii = []
        self._coreMinorRadii = []
        self._base = base
        self._end = end
        self._shield = None
        self._elementsCountAcrossMinor = base._elementsCountAcrossMinor
        self._elementsCountAcrossMajor = base._elementsCountAcrossMajor
        self._elementsCountUp = base._elementsCountAcrossMajor // 2 \
            if cylinderShape == CylinderShape.CYLINDER_SHAPE_FULL else base._elementsCountAcrossMajor
        self._elementsCountAcrossShell = base._elementsCountAcrossShell
        self._elementsCountAcrossTransition = base._elementsCountAcrossTransition
        self._elementsCountAcrossRim = self._elementsCountAcrossShell + self._elementsCountAcrossTransition - 1
        self._shellProportion = base._shellProportion
        self._elementsCountAlong = elementsCountAlong
        self._elementsCountAround = 2 * (self._elementsCountAcrossMajor+self._elementsCountAcrossMinor -
                                         4*(self._elementsCountAcrossRim + 1))
        self._startNodeIdentifier = 1
        self._startElementIdentifier = 1
        self._endNodeIdentifier = 1
        self._endElementIdentifier = 1
        self._cylinderShape = cylinderShape
        self._cylinderType = CylinderType.CYLINDER_STRAIGHT
        if (tapered is not None) or cylinderCentralPath:
            self._cylinderType = CylinderType.CYLINDER_TAPERED
            self._tapered = tapered
        self._useCrossDerivatives = useCrossDerivatives
        if rangeOfRequiredElements:
            self._rangeOfRequiredElements = rangeOfRequiredElements
        else:
            self._rangeOfRequiredElements = [
                [0, self._elementsCountAcrossMajor],
                [0, self._elementsCountAcrossMinor],
                [0, self._elementsCountAlong],
            ]

        self._meshGroups = meshGroups
        self._meshGroupsElementsAlong = meshGroupsElementsAlong
        self._cylinderCentralPath = cylinderCentralPath
        if cylinderCentralPath:
            self.calculateEllipseParams(cylinderCentralPath=self._cylinderCentralPath)
            self._base = CylinderEnds(base._elementsCountAcrossMajor, base._elementsCountAcrossMinor,
                                      base._elementsCountAcrossShell, base._elementsCountAcrossTransition,
                                      base._shellProportion, self._centres[0],
                                      None, self._majorAxis[0], self._minorRadii[0])
        else:
            self._length = vector.magnitude(base._alongAxis)
            arcLengthAlong = self._length / elementsCountAlong
            self.calculateEllipseParams(arcLengthAlong, cylinderCentralPath=self._cylinderCentralPath)

        # generate the mesh
        self.createCylinderMesh3d(fieldModule, coordinates)

    def createCylinderMesh3d(self, fieldModule, coordinates):
        """
        Create an extruded shape (ellipse/circle) mesh. Currently limited to ellipse or circle base with the alongAxis
        perpendicular to the base.
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :return: Final values of nextNodeIdentifier, nextElementIdentifier.
        """
        assert (self._elementsCountAlong > 0), 'createCylinderMesh3d:  Invalid number of along elements'
        assert (self._elementsCountAcrossMinor > 3), 'createCylinderMesh3d: Invalid number of across elements'
        assert (self._elementsCountAcrossMinor % 2 == 0), 'createCylinderMesh3d: number of across elements' \
                                                          ' is not an even number'
        assert (self._elementsCountAcrossMajor > 1), 'createCylinderMesh3d: Invalid number of up elements'
        assert (self._cylinderShape in [self._cylinderShape.CYLINDER_SHAPE_FULL,
                                        self._cylinderShape.CYLINDER_SHAPE_LOWER_HALF]), \
            'createCylinderMesh3d: Invalid cylinder mode.'

        nodes = fieldModule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldModule.findMeshByDimension(3)

        elementsCountRim = self._elementsCountAcrossRim

        shieldMode = ShieldShape2D.SHIELD_SHAPE_FULL if self._cylinderShape is self._cylinderShape.CYLINDER_SHAPE_FULL \
            else ShieldShape2D.SHIELD_SHAPE_LOWER_HALF
        ellipseShape = EllipseShape.Ellipse_SHAPE_FULL \
            if self._cylinderShape is self._cylinderShape.CYLINDER_SHAPE_FULL else EllipseShape.Ellipse_SHAPE_LOWER_HALF
        self._shield = ShieldMesh2D(self._elementsCountAcrossMinor, self._elementsCountAcrossMajor, elementsCountRim,
                                  None, self._elementsCountAlong, shieldMode,
                                  shieldType=ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND)

        # generate ellipses mesh along cylinder axis
        n3Count = 0 if self._cylinderType == CylinderType.CYLINDER_STRAIGHT else self._elementsCountAlong
        self._ellipses = []
        for n3 in range(n3Count + 1):
            ellipse = Ellipse2D(self._centres[n3], self._majorAxis[n3], self._minorAxis[n3],
                                self._elementsCountAcrossMajor, self._elementsCountAcrossMinor, self._elementsCountAcrossShell,
                                self._elementsCountAcrossTransition,
                                self._shellProportion, self._coreMajorRadii[n3], self._coreMinorRadii[n3],
                                ellipseShape=ellipseShape)
            self._ellipses.append(ellipse)
            self.copyEllipsesNodesToShieldNodes(n3)

        for n3 in range(n3Count + 1):
            self.calculateD2Derivatives(n3, n3Count)

        if self._cylinderType == CylinderType.CYLINDER_TAPERED:
            self.smoothd2Derivatives()

        # The other ellipses for a straight cylinder.
        if self._cylinderType == CylinderType.CYLINDER_STRAIGHT:
            arcLengthAlong = vector.magnitude(self._base._alongAxis) / self._elementsCountAlong
            for n2 in range(self._elementsCountUp + 1):
                for n3 in range(self._elementsCountAlong + 1):
                    for n1 in range(self._elementsCountAcrossMinor + 1):
                        if self._shield.px[0][n2][n1]:
                            temx = [self._shield.px[0][n2][n1][c] + n3 * arcLengthAlong *
                                    vector.normalise(self._base._alongAxis)[c] for c in range(3)]
                            self._shield.px[n3][n2][n1] = temx
                            self._shield.pd1[n3][n2][n1] = self._shield.pd1[0][n2][n1]
                            self._shield.pd2[n3][n2][n1] = self._shield.pd2[0][n2][n1]
                            self._shield.pd3[n3][n2][n1] = self._shield.pd3[0][n2][n1]

        self.generateNodes(nodes, fieldModule, coordinates)
        self.generateElements(mesh, fieldModule, coordinates)

        if self._end is None:
            self._end = CylinderEnds(self._elementsCountAcrossMajor, self._elementsCountAcrossMinor,
                                     self._elementsCountAcrossShell, self._elementsCountAcrossTransition,
                                     self._shellProportion,
                                     self._centres[-1], self._shield.pd2[-1][0][1],
                                     vector.setMagnitude(self._base._majorAxis, self._majorRadii[-1]),
                                     self._minorRadii[-1])
        self.setEndsNodes()

    def calculateD2Derivatives(self, n3, n3Count):
        """
        calculate d2 derivatives.
        :param n3: Index of along cylinder axis coordinates to use
        :param n3Count: number of bases to create coordinates for.
        """
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3
        # get ellipse d2 and next ellipse x, d2
        for n2 in range(self._elementsCountAcrossMajor + 1):
            for n1 in range(self._elementsCountAcrossMinor + 1):
                if btd1[n3][n2][n1]:
                    n3n = n3 if (n3 < n3Count) else n3 - 1
                    btd2[n3][n2][n1] = [(btx[n3n + 1][n2][n1][c] - btx[n3n][n2][n1][c]) for c in range(3)]

    def smoothd2Derivatives(self):
        """
        smooth d2 derivatives using initial values calculated by calculateD2Derivatives
        """
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3
        for n2 in range(self._elementsCountAcrossMajor + 1):
            for n1 in range(self._elementsCountAcrossMinor + 1):
                td2 = []
                tx = []
                if btx[0][n2][n1]:
                    for n3 in range(self._elementsCountAlong + 1):
                        tx.append(btx[n3][n2][n1])
                        td2.append(btd2[n3][n2][n1])
                    td2 = smoothCubicHermiteDerivativesLine(tx, td2, fixStartDirection=True)
                    for n3 in range(self._elementsCountAlong + 1):
                        btd2[n3][n2][n1] = td2[n3]

    def setEndsNodes(self):
        """
        sets ellipse coordinates, derivatives and node ids.
        """
        self._base.px = self._shield.px[0]
        self._base.pd1 = self._shield.pd1[0]
        self._base.pd2 = self._shield.pd2[0]
        self._base.pd3 = self._shield.pd3[0]
        self._end.px = self._shield.px[-1]
        self._end.pd1 = self._shield.pd1[-1]
        self._end.pd2 = self._shield.pd2[-1]
        self._end.pd3 = self._shield.pd3[-1]

    def calculateEllipseParams(self, arcLengthAlong=None, cylinderCentralPath=None):
        """
        Calculate the ellipses major and minor radii, major and minor axis and ellipses centres.
        :param arcLengthAlong: arc length along the cylinder axis. Only if cylinderCentralPath is false.
        :param cylinderCentralPath: Stores radii and centres of the ellipses along the cylinder length
        """
        if not cylinderCentralPath:
            self._centres = [self._base._centre for _ in range(self._elementsCountAlong+1)]
            self._majorAxis = [self._base._majorAxis for _ in range(self._elementsCountAlong+1)]
            self._minorAxis = [self._base._minorAxis for _ in range(self._elementsCountAlong+1)]
            self._majorRadii = [self._base._majorRadius for _ in range(self._elementsCountAlong+1)]
            self._minorRadii = [self._base._minorRadius for _ in range(self._elementsCountAlong+1)]
        if cylinderCentralPath:
            self._centres = cylinderCentralPath.centres
            self._majorAxis = cylinderCentralPath.majorAxis
            self._minorAxis = cylinderCentralPath.minorAxis
            self._majorRadii = cylinderCentralPath.majorRadii
            self._minorRadii = cylinderCentralPath.minorRadii

            elementsMinor = self._elementsCountAcrossMinor//2 - self._elementsCountAcrossShell*(1-self._shellProportion)
            elementsMajor = self._elementsCountUp - self._elementsCountAcrossShell*(1-self._shellProportion)
            for n3 in range(self._elementsCountAlong + 1):
                self._coreMinorRadii.append((1-self._shellProportion*self._elementsCountAcrossShell/elementsMinor)*self._minorRadii[n3])
                self._coreMajorRadii.append((1-self._shellProportion*self._elementsCountAcrossShell/elementsMajor)*self._majorRadii[n3])

        if self._cylinderType == CylinderType.CYLINDER_TAPERED and not cylinderCentralPath:
            centre = self._base._centre
            majorRadius = self._base._majorRadius
            minorRadius = self._base._minorRadius
            for n3 in range(1, self._elementsCountAlong+1):
                majorRadius, majorAxis = computeNextRadius(majorRadius, self._base._majorAxis,
                                                           self._tapered.majorRatio,
                                                           self._tapered.majorProgressionMode)
                minorRadius, minorAxis = computeNextRadius(minorRadius, self._base._minorAxis,
                                                           self._tapered.minorRatio,
                                                           self._tapered.minorProgressionMode)
                centre = computeNextCentre(centre, arcLengthAlong, self._base._alongAxis)
                self._centres[n3] = centre
                self._majorAxis[n3] = majorAxis
                self._minorAxis[n3] = minorAxis
                self._majorRadii[n3] = majorRadius
                self._minorRadii[n3] = minorRadius

    def copyEllipsesNodesToShieldNodes(self, n3):
        """
        Copy coordinates and derivatives of ellipses to shield.
        :param n3: the index number of ellipse along the central path.
        """
        self._shield.px[n3] = self._ellipses[n3].px
        self._shield.pd1[n3] = self._ellipses[n3].pd1
        self._shield.pd2[n3] = self._ellipses[n3].pd2
        self._shield.pd3[n3] = self._ellipses[n3].pd3

    def generateNodes(self, nodes, fieldModule, coordinates):
        """
        Create cylinder nodes from coordinates.
        :param nodes: nodes from coordinates.
        :param fieldModule: Zinc fieldmodule to create nodes in. Uses DOMAIN_TYPE_NODES.
        :param coordinates: Coordinate field to define.
        """
        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self._shield.generateNodes(fieldModule, coordinates, nodeIdentifier,
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
        elementIdentifier = self._shield.generateElements(fieldModule, coordinates, elementIdentifier,
                                                          self._rangeOfRequiredElements,
                                                          self._meshGroupsElementsAlong, self._meshGroups)
        self._endElementIdentifier = elementIdentifier

    def getElementsCountAround(self):
        return self._elementsCountAround

    def getElementIdentifiers(self):
        return self._shield.elementId


class Ellipse2D:
    """
    Generate a 2D ellipse.
    """

    def __init__(self, centre, majorAxis, minorAxis,
                 elementsCountAcrossMajor, elementsCountAcrossMinor, elementsCountAcrossShell,
                 elementsCountAcrossTransition, shellProportion, coreMajorRadius, coreMinorRadius,
                 ellipseShape=EllipseShape.Ellipse_SHAPE_FULL):
        """
        :param centre: Ellipse centre.
        :param majorAxis: A vector for ellipse major axis.
        :param minorAxis: Ellipse minor axis.
        :param elementsCountAcrossMajor:
        :param elementsCountAcrossMinor:
        :param ellipseShape: The shape of the ellipse which can be full or lower half.
        """
        self.centre = centre
        self.majorAxis = majorAxis
        self.minorAxis = minorAxis
        self.majorRadius = vector.magnitude(majorAxis)
        self.minorRadius = vector.magnitude(minorAxis)
        self.coreMajorRadius = coreMajorRadius
        self.coreMinorRadius = coreMinorRadius
        elementsCountRim = elementsCountAcrossShell + elementsCountAcrossTransition - 1
        shieldMode = ShieldShape2D.SHIELD_SHAPE_FULL if ellipseShape is EllipseShape.Ellipse_SHAPE_FULL\
            else ShieldShape2D.SHIELD_SHAPE_LOWER_HALF
        shield = ShieldMesh2D(elementsCountAcrossMinor, elementsCountAcrossMajor, elementsCountRim,
                            None, 1, shieldMode,
                            shieldType=ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND)
        self.elementsCountAcrossMajor = elementsCountAcrossMajor
        self.elementsCountAround = shield.elementsCountAroundFull
        self.elementsCountUp = shield.elementsCountUp
        self.elementsCountAcrossMinor = elementsCountAcrossMinor
        self.elementsCountAcrossShell = elementsCountAcrossShell
        self.elementsCountAcrossTransition = elementsCountAcrossTransition
        self.elementsCountAcrossRim = elementsCountRim
        self.shellProportion = shellProportion
        self.nodeId = shield.nodeId
        self.px = shield.px[0]
        self.pd1 = shield.pd1[0]
        self.pd2 = shield.pd2[0]
        self.pd3 = shield.pd3[0]
        self.__shield = shield
        self.ellipseShape = ellipseShape
        # generate the ellipse
        self.generate2DEllipseMesh()

    def generate2DEllipseMesh(self):
        """
        Generates a 2d ellipse using shield structure in shieldmesh.
        """
        self.generateBase1DMesh(0)
        if self.elementsCountAcrossShell > 0:
            self.generateBase1DMesh(self.elementsCountAcrossShell)
        self.setRimNodes()
        rscx, rscd1, rscd2, rscd3 = self.createMirrorCurve()
        self.createRegularRowCurves(rscx, rscd1, rscd3)
        self.createRegularColumnCurves()
        self.__shield.getTriplePoints(0)
        self.smoothTriplePointsCurves()
        self.smoothTransitionRims()
        if self.ellipseShape == EllipseShape.Ellipse_SHAPE_FULL:
            self.generateNodesForUpperHalf()
        self.calculateD2()

    def generateBase1DMesh(self, rx):
        """
        Generate nodes around the perimeter of the ellipse.
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3

        ratio = rx/self.elementsCountAcrossShell if self.elementsCountAcrossShell > 0 else 0
        majorAxis = [d*(1 - ratio*(1-self.coreMajorRadius/self.majorRadius)) for d in self.majorAxis]
        minorAxis = [d*(1 - ratio*(1-self.coreMinorRadius/self.minorRadius)) for d in self.minorAxis]
        majorRadius = vector.magnitude(majorAxis)

        nx, nd1 = createEllipsePerimeter(
            self.centre, majorAxis, minorAxis, self.elementsCountAround, majorRadius)
        nte = normalToEllipse(self.majorAxis, self.minorAxis)

        tbx, tbd1, tbd2, tbd3 = [], [], [], []
        for n in range(self.elementsCountAround + 1):
            tbx.append(nx[n])
            tbd1.append(nd1[n])
            tbd2.append(nte)
            tbd3.append(vector.normalise(vector.crossproduct3(tbd1[n], nte)))

        for n in range(self.elementsCountAround + 1):
            n1, n2 = self.__shield.convertRimIndex(n, rx)
            btx[n2][n1] = tbx[n]
            btd1[n2][n1] = tbd1[n]
            btd2[n2][n1] = tbd2[n]
            btd3[n2][n1] = tbd3[n]

    def setRimNodes(self):
        """
        Set nodes around the shell outer layer and core boundary in order needed for creating a shield mesh.
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3
        #
        elementsCountShell = self.elementsCountAcrossShell
        for n in range(self.elementsCountAround + 1):
            if self.elementsCountAcrossShell > 1:
                n1, n2 = self.__shield.convertRimIndex(n)
                n1c, n2c = self.__shield.convertRimIndex(n, self.elementsCountAcrossShell)
                tx, td3, pe, pxi, psf = sampleCubicHermiteCurves([btx[n2c][n1c], btx[n2][n1]],
                                                                 [btd3[n2c][n1c], btd3[n2][n1]], elementsCountShell,
                                                                 arcLengthDerivatives=True)

                for rx in range(1, self.elementsCountAcrossShell):
                    n1rx, n2rx = self.__shield.convertRimIndex(n, rx)
                    btx[n2rx][n1rx] = tx[-rx + elementsCountShell]
                    btd1[n2rx][n1rx] = btd1[n2][n1]
                    btd2[n2rx][n1rx] = btd2[n2][n1]
                    btd3[n2rx][n1rx] = btd3[n2][n1]

        for rx in range(1, self.elementsCountAcrossShell):
            self.__shield.smoothDerivativesAroundRim(0, n3d=None, rx=rx)

    def createMirrorCurve(self):
        """
        generate coordinates and derivatives for the mirror curve
        :return: Coordinates and derivatives for the mirror curve
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3

        n2a = self.elementsCountAcrossShell
        rcx = []
        tmdx = btx[n2a][self.elementsCountAcrossMinor // 2]
        tmdd3 = btd3[n2a][self.elementsCountAcrossMinor // 2]
        tmux = [
            0.5 * (btx[self.elementsCountUp][0][c] + btx[self.elementsCountUp][self.elementsCountAcrossMinor]
            [c]) for c in range(3)]
        rcx.append(tmdx)
        rcx.append(tmux)
        rcd3 = [vector.setMagnitude(tmdd3, -1), vector.setMagnitude(tmdd3, -1)]
        rscx, rscd1 = sampleCubicHermiteCurves(rcx, rcd3, self.elementsCountUp - n2a, arcLengthDerivatives=True)[0:2]

        # get d2, d3
        rscd2 = []
        rscd3 = []
        for n in range(len(rscx)):
            d3 = vector.normalise(
                [btx[self.elementsCountUp][self.elementsCountAcrossMinor][c] - btx[self.elementsCountUp][0]
                [c] for c in range(3)])
            d2 = vector.normalise(vector.crossproduct3(d3, rscd1[n]))
            rscd2.append(d2)
            rscd3.append(d3)

        return rscx, rscd1, rscd2, rscd3

    def createRegularRowCurves(self, rscx, rscd1, rscd3):
        """
        generate curves along regular rows using the mirror curve obtained from createMirrorCurve.
        :param rscx: Coordinates of the nodes on the middle curve.
        :param rscd1: d1 derivatives of the nodes on the middle curve.
        :param rscd3: d3 derivatives of the nodes on the middle curve.
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3

        elementsCountRim = self.elementsCountAcrossRim
        n2a = self.elementsCountAcrossShell
        n2b = elementsCountRim
        n2d = n2b + 2
        n2m = self.elementsCountUp
        n1a = self.elementsCountAcrossShell
        n1b = elementsCountRim
        n1z = self.elementsCountAcrossMinor - self.elementsCountAcrossShell
        for n2 in range(n2d, n2m + 1):
            txm, td3m, pe, pxi, psf = sampleCubicHermiteCurves(
                [btx[n2][n1a], rscx[n2 - n2a], btx[n2][n1z]],
                [vector.setMagnitude(btd3[n2][n1a], -1.0), rscd3[n2 - n2a], btd3[n2][n1z]],
                self.elementsCountAcrossMinor-2*self.elementsCountAcrossShell, arcLengthDerivatives=True)
            td1m = interpolateSampleCubicHermite([[-btd1[n2][n1a][c] for c in range(3)], rscd1[n2 - n2a],
                                                  btd1[n2][n1z]], [[0.0, 0.0, 0.0]] * 3, pe, pxi, psf)[0]

            tx = []
            td3 = []
            for n1 in range(self.elementsCountAcrossMinor + 1):
                if n1 <= n1a:
                    tx.append(btx[n2][n1])
                    td3.append([-btd3[n2][n1][c] for c in range(3)])
                elif (n1 > n1a) and (n1 < n1z):
                    tx.append(txm[n1 - n1a])
                    td3.append(td3m[n1 - n1a])
                else:
                    tx.append(btx[n2][n1])
                    td3.append((btd3[n2][n1]))

            td3 = smoothCubicHermiteDerivativesLine(tx, td3, fixStartDirection=True, fixEndDirection=True)

            for n1 in range(self.elementsCountAcrossMinor + 1):
                if n1 <= n1a:
                    btd3[n2][n1] = [-td3[n1][c] for c in range(3)]
                elif (n1 > n1a) and (n1 < n1z):
                    btx[n2][n1] = tx[n1]
                    if n2 == n2m:
                        if n1 <= n1b:
                            btd1[n2][n1] = vector.setMagnitude(btd1[n2m][-1], -1.0)
                        else:
                            btd1[n2][n1] = vector.setMagnitude(btd1[n2m][-1], 1.0)
                    else:
                        if n1 <= n1b:
                            btd1[n2][n1] = [-d for d in td1m[n1 - n1a]]
                        else:
                            btd1[n2][n1] = td1m[n1 - n1a]
                    if n1 <= n1b:
                        btd3[n2][n1] = [-d for d in td3[n1]]
                    else:
                        btd3[n2][n1] = td3[n1]
                else:
                    btd3[n2][n1] = td3[n1]

    def createRegularColumnCurves(self):
        """
        up regular columns of shield: get d1, initial d3 below regular rows
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3

        elementsCountRim = self.elementsCountAcrossRim
        n1d = elementsCountRim + 2
        n1x = self.elementsCountAcrossMinor - elementsCountRim - 1
        n2a = self.elementsCountAcrossShell
        n2b = elementsCountRim
        n2d = 2 + n2b
        n2m = self.elementsCountUp
        for n1 in range(n1d, n1x):
            txm, td1m, pe, pxi, psf = sampleCubicHermiteCurves(
                [btx[n2a][n1], btx[n2d][n1]], [[-btd3[n2a][n1][c] for c in range(3)], btd1[n2d][n1]],
                2 + self.elementsCountAcrossTransition - 1, arcLengthDerivatives=True)
            td3m = interpolateSampleCubicHermite([btd1[n2a][n1], btd3[n2d][n1]], [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]

            tx = []
            td1 = []
            for n2 in range(n2m + 1):
                if n2 <= n2a:
                    tx.append(btx[n2][n1])
                    td1.append([-btd3[n2][n1][c] for c in range(3)])
                elif (n2 > n2a) and (n2 < n2d):
                    tx.append(txm[n2 - n2a])
                    td1.append(td1m[n2 - n2a])
                else:
                    tx.append(btx[n2][n1])
                    td1.append(btd1[n2][n1])

            td1 = smoothCubicHermiteDerivativesLine(tx, td1, fixStartDirection=True, fixEndDirection=True)

            for n2 in range(n2m + 1):
                if n2 <= n2a:
                    btd3[n2][n1] = [-td1[n2][c] for c in range(3)]
                elif (n2 > n2a) and (n2 < n2d):
                    btx[n2][n1] = tx[n2]
                    if n2 <= n2b:
                        btd1[n2][n1] = td3m[n2 - n2a]
                        btd3[n2][n1] = [-d for d in td1[n2]]
                    else:
                        btd1[n2][n1] = td1[n2]
                        btd3[n2][n1] = td3m[n2 - n2a]
                else:
                    btd1[n2][n1] = td1[n2]

    def smoothTriplePointsCurves(self):
        """
        Smooth row and column curves passing triple points (i.e., row 1 and columns 1 and -2).
        """
        btx = self.px
        btd1 = self.pd1
        btd2 = self.pd2
        btd3 = self.pd3

        n1c = 1 + self.elementsCountAcrossRim
        n1y = self.elementsCountAcrossMinor - self.elementsCountAcrossRim
        n1x = n1y - 1
        n2a = self.elementsCountAcrossShell
        n2c = 1 + self.elementsCountAcrossRim
        n2m = self.elementsCountUp

        # smooth shield row n2c
        btd3[n2c][n1c:n1y] = smoothCubicHermiteDerivativesLine(btx[n2c][n1c:n1y], btd3[n2c][n1c:n1y])

        # smooth Shield columns n1c, n1x
        for n1 in [n1c, n1x]:
            tx = []
            td1 = []
            for n2 in range(n2c, n2m + 1):
                tx.append(btx[n2][n1])
                td1.append(btd1[n2][n1])
            td1 = smoothCubicHermiteDerivativesLine(tx, td1, fixEndDirection=True, fixStartDerivative=True)
            for n in range(n2m-n2c+1):
                btd1[n + n2c][n1] = td1[n]

        # sample nodes to triple points
        for n1 in [n1c, n1x]:
            c1 = 1 if n1 == n1c else -1
            txm, td3m, pe, pxi, psf = sampleCubicHermiteCurves(
                [btx[n2a][n1], btx[n2c][n1]], [[-btd3[n2a][n1][c] for c in range(3)],
                                               [btd1[n2c][n1][c] + c1 * btd3[n2c][n1][c] for c in range(3)]],
                1 + self.elementsCountAcrossTransition - 1, arcLengthDerivatives=True)
            td1m = interpolateSampleCubicHermite([btd1[n2a][n1], [-c1*d for d in btd1[n2c][n1]]], [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]

            for n2 in range(n2a + 1, n2c):
                btx[n2][n1] = txm[n2 - n2a]
                btd1[n2][n1] = td1m[n2 - n2a]
                btd3[n2][n1] = [-d for d in td3m[n2 - n2a]]

        # smooth
        self.__shield.smoothDerivativesToTriplePoints(0, fixAllDirections=True)

    def smoothTransitionRims(self):
        """
        smooth rims in the transition zone
        """
        n2a = self.elementsCountAcrossShell
        n2b = self.elementsCountAcrossRim
        for rx in range(n2a + 1, n2b + 1):
            self.__shield.smoothDerivativesAroundRim(0, n3d=None, rx=rx)

    def generateNodesForUpperHalf(self):
        """
        Generates coordinates and derivatives for the upper half by mirroring the lower half nodes and derivatives.
         It keeps the d1 direction.
        It uses mirrorPlane: plane ax+by+cz=d in form of [a,b,c,d]
        """
        mirrorPlane = [-d for d in self.majorAxis] + [-vector.dotproduct(self.majorAxis, self.centre)]
        mirror = Mirror(mirrorPlane)
        for n2 in range(self.elementsCountUp):
            for n1 in range(self.elementsCountAcrossMinor + 1):
                if self.px[n2][n1]:
                    self.px[2 * self.elementsCountUp - n2][n1] = mirror.mirrorImageOfPoint(
                        self.px[n2][n1])
                    self.pd1[2 * self.elementsCountUp - n2][n1] = mirror.reverseMirrorVector(
                        self.pd1[n2][n1])
                    self.pd3[2 * self.elementsCountUp - n2][n1] = mirror.mirrorVector(
                        self.pd3[n2][n1])

    def calculateD2(self):
        """
        :return:
        """
        btx = self.px
        btd2 = self.pd2

        nte = normalToEllipse(self.majorAxis, self.minorAxis)
        for n2 in range(self.elementsCountAcrossMajor + 1):
            for n1 in range(self.elementsCountAcrossMinor + 1):
                if btx[n2][n1]:
                    btd2[n2][n1] = nte

    def getShield(self):
        return self.__shield


def createEllipsePerimeter(centre, majorAxis, minorAxis, elementsCountAround, height):
    """
    Generate a set of points and derivatives for an ellipse
    starting at pole majorAxis from centre.
    :param elementsCountAround: Number of elements around.
    :param centre: Centre of full ellipse.
    :param majorAxis: Vector in direction of starting major radius, magnitude is ellipse major radius.
    :param minorAxis: Vector normal to major axis, magnitude is ellipse minor axis length.
    :param height: Height of arc of ellipsoid from starting point along majorAxis.
    :return: Lists nx, nd1. Ordered fastest around, starting at major radius.
    """
    nx = []
    nd1 = []
    magMajorAxis = vector.magnitude(majorAxis)
    magMinorAxis = vector.magnitude(minorAxis)
    unitMajorAxis = vector.normalise(majorAxis)
    unitMinorAxis = vector.normalise(minorAxis)
    useHeight = min(max(0.0, height), 2.0 * magMajorAxis)

    totalRadians = math.acos((magMajorAxis - useHeight)/magMajorAxis)
    radians = 0.0

    arcLengthUp = geometry.getEllipseArcLength(magMajorAxis, magMinorAxis, radians, totalRadians, 'integrate')
    elementsCountUp = elementsCountAround // 2
    elementArcLength = arcLengthUp / elementsCountUp

    radians = geometry.updateEllipseAngleByArcLength(magMajorAxis, magMinorAxis, radians, -arcLengthUp, method='Newton')
    for n1 in range(2 * elementsCountUp + 1):
        cosRadians = math.cos(radians)
        sinRadians = math.sin(radians)
        nx.append(
            [(centre[c] + cosRadians * majorAxis[c] + sinRadians * minorAxis[c]) for c in range(3)])

        ndab = vector.setMagnitude([-sinRadians * magMajorAxis, cosRadians * magMinorAxis], elementArcLength)
        nd1.append(
            [(ndab[0] * unitMajorAxis[c] + ndab[1] * unitMinorAxis[c]) for c in range(3)])
        radians = geometry.updateEllipseAngleByArcLength(magMajorAxis, magMinorAxis, radians, elementArcLength,
                                                         method='Newton')
    return nx, nd1


def normalToEllipse(v1, v2):
    """
    Find unit normal vector of an ellipse using two vectors in the ellipse. The direction is v1xv2
    :param v1: vector 1.
    :param v2: vector 2.
    :return:
    """
    nte = vector.normalise(vector.crossproduct3(v1, v2))
    return nte


def computeNextRadius(radius, axis, ratio, progression):
    """
    calculate next radius based on the progression method. r_n+1=r_n*ratio for geometric. r_n+1=r_ratio for arithmetic.
    :param radius: radius (major or minor) along the central path.
    :param axis: major or minor axis along the central path.
    :param ratio: common ratio (common difference) for changing the next radius.
    :param progression: arithmetic or geometric.
    :return: next radius and axis.
    """
    if progression == ConeBaseProgression.GEOMETRIC_PROGRESSION:
        radius = radius * ratio
    elif progression == ConeBaseProgression.ARITHMETIC_PROGRESSION:
        radius += ratio
    axis = vector.setMagnitude(axis, radius)
    return radius, axis


def computeNextCentre(centre, arcLength, axis):
    """
    compute next centre coordinate
    :param axis:
    :param arcLength: the length to go forward.
    :param centre: the start centre.
    :return: next centre coordinates.
    """
    centre = [centre[c]+arcLength * vector.normalise(axis)[c] for c in range(3)]
    return centre
