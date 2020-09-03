"""
Utility functions for generating a generalised 3-D solid cylinder (extruded ellipse/circle). It can be used to generate
a solid truncated cone. It also can be used for transition from a 2D base to another base (e.g., ellipse to a circle).
"""

from enum import Enum
from scaffoldmaker.utils import vector, geometry
import math
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.utils.zinc.finiteelement import getMaximumNodeIdentifier, getMaximumElementIdentifier
from scaffoldmaker.utils.shieldmesh import ShieldMesh, ShieldShape, ShieldRimDerivativeMode
from scaffoldmaker.utils.interpolation import sampleCubicHermiteCurves, interpolateSampleCubicHermite, \
    smoothCubicHermiteDerivativesLine


class CylinderShape(Enum):
    CYLINDER_SHAPE_FULL = 1  # full cylinder is created
    CYLINDER_SHAPE_LOWER_HALF = 2  # lower half cylinder


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

    def __init__(self, elementsCountAcrossMajor, elementsCountAcrossMinor,
                 centre, alongAxis, majorAxis, minorRadius):
        """
        :param elementsCountAcrossMajor: Number of elements across major axis. Must be at least 2 + elementsCountRim for
         half and 4 + elementsCountRim for full cylinder.
        :param elementsCountAcrossMinor: Number of elements across minor axis.
        :param centre: Centre of the ellipse.
        :param alongAxis: The cylinder axis that the base is extruded along.
        :param majorAxis: The major axis of the base. Should be perpendicular to alongAxis
        :param minorRadius: The minor radius of the ellipse.
        """
        self._centre = centre
        self._alongAxis = alongAxis
        self._majorAxis = majorAxis
        self._minorRadius = minorRadius
        self._minorAxis = vector.setMagnitude(vector.crossproduct3(alongAxis, majorAxis), minorRadius)
        self._elementsCountAcrossMinor = elementsCountAcrossMinor
        self._elementsCountAcrossMajor = elementsCountAcrossMajor
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


class CylinderMesh:
    """
    Cylinder mesh generator. Extrudes an ellipse/circle.
    """

    def __init__(self, fieldModule, coordinates, base, elementsCountAlong, end=None,
                 cylinderShape=CylinderShape.CYLINDER_SHAPE_FULL,
                 tapered=None, useCrossDerivatives=False):
        """
        :param fieldModule: Zinc fieldModule to create elements in.
        :param coordinates: Coordinate field to define.
        :param base: Cylinder base ellipse. It is an instance of class CylinderEnds.
        :param end: Cylinder end ellipse. It is an instance of class CylinderEnds.
        :param elementsCountAlong: Number of elements along the cylinder axis.
        :param cylinderShape: A value from enum CylinderMode specifying.
        """
        self._base = base
        self._end = end
        self._shield = None
        self._elementsCountAcrossMinor = base._elementsCountAcrossMinor
        self._elementsCountAcrossMajor = base._elementsCountAcrossMajor
        self._elementsCountUp = base._elementsCountAcrossMajor // 2 \
            if cylinderShape == CylinderShape.CYLINDER_SHAPE_FULL else base._elementsCountAcrossMajor
        self._elementsCountAlong = elementsCountAlong
        self._startNodeIdentifier = 1
        self._startElementIdentifier = 1
        self._endNodeIdentifier = 1
        self._endElementIdentifier = 1
        self._cylinderShape = cylinderShape
        self._cylinderType = CylinderType.CYLINDER_STRAIGHT
        if tapered is not None:
            self._cylinderType = CylinderType.CYLINDER_TAPERED
            self._tapered = tapered
        self._useCrossDerivatives = useCrossDerivatives
        self._length = vector.magnitude(base._alongAxis)
        self._majorRadius = vector.magnitude(base._majorAxis)
        self._basesCentres = [None for _ in range(elementsCountAlong + 1)]
        self._basesCentres[0] = self._base._centre
        # generate the mesh
        self.createCylinderMesh3d(fieldModule, coordinates)

    @staticmethod
    def createCylinderBaseMesh2D(centre, majorAxis, minorAxis, elementsCountAround, height):
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
        totalRadians = geometry.getEllipseRadiansToX(magMajorAxis, 0.0, magMajorAxis - useHeight,
                                                     initialTheta=0.5 * math.pi * useHeight / magMajorAxis)
        radians = 0.0
        arcLengthUp = geometry.getEllipseArcLength(magMajorAxis, magMinorAxis, radians, totalRadians)
        elementsCountUp = elementsCountAround // 2
        elementArcLength = arcLengthUp / elementsCountUp
        radians = geometry.updateEllipseAngleByArcLength(magMajorAxis, magMinorAxis, radians, -arcLengthUp)
        for n1 in range(2 * elementsCountUp + 1):
            cosRadians = math.cos(radians)
            sinRadians = math.sin(radians)
            nx.append(
                [(centre[c] + cosRadians * majorAxis[c] + sinRadians * minorAxis[c]) for c in range(3)])

            ndab = vector.setMagnitude([-sinRadians * magMajorAxis, cosRadians * magMinorAxis], elementArcLength)
            nd1.append(
                [(ndab[0] * unitMajorAxis[c] + ndab[1] * unitMinorAxis[c]) for c in range(3)])
            radians = geometry.updateEllipseAngleByArcLength(magMajorAxis, magMinorAxis, radians, elementArcLength)
        return nx, nd1

    def generateBasesMesh(self, majorRadius, elementsCountAround, arcLengthAlong, minorAxis):
        """
        generate bases of the truncated cone along the cone axis.
        :param minorAxis:
        :param arcLengthAlong:
        :param majorRadius: major radius of the cone ellipse base.
        :param elementsCountAround: major radius of the cone ellipse base.
        :return:
        """
        self._majorRadii = []
        self._minorRadii = []
        nx, nd1 = self.createCylinderBaseMesh2D(
            self._base._centre, self._base._majorAxis, minorAxis, elementsCountAround, majorRadius)
        majorRadius1 = majorRadius
        self._majorRadii.append(majorRadius1)
        minorRadius1 = vector.magnitude(minorAxis)
        self._minorRadii.append(minorRadius1)
        tnx, tnd1, tnd2, tnd3 = [], [], [], []
        self._basesCentres = [self._base._centre for _ in range(self._elementsCountAlong + 1)]
        for n3 in range(self._elementsCountAlong + 1):
            tbx, tbd1, tbd2, tbd3 = [], [], [], []
            for n in range(elementsCountAround + 1):
                tbx.append(nx[n])
                tbd1.append(nd1[n])
                tbd2.append([arcLengthAlong * vector.normalise(self._base._alongAxis)[c] for c in range(3)])
                tbd3.append(vector.normalise(vector.crossproduct3(tbd1[n], tbd2[n])))
            tnx.append(tbx)
            tnd1.append(tbd1)
            tnd2.append(tbd2)
            tnd3.append(tbd3)
            if (self._cylinderType == CylinderType.CYLINDER_TAPERED) and (n3 < self._elementsCountAlong):
                nx = nd1 = []
                if self._tapered.majorProgressionMode == ConeBaseProgression.GEOMETRIC_PROGRESSION:
                    majorRadius1 = majorRadius1 * self._tapered.majorRatio
                elif self._tapered.majorProgressionMode == ConeBaseProgression.ARITHMETIC_PROGRESSION:
                    majorRadius1 += self._tapered.majorRatio
                majorAxis1 = vector.setMagnitude(self._base._majorAxis, majorRadius1)
                if self._tapered.minorProgressionMode == ConeBaseProgression.GEOMETRIC_PROGRESSION:
                    minorRadius1 = minorRadius1 * self._tapered.minorRatio
                elif self._tapered.minorProgressionMode == ConeBaseProgression.ARITHMETIC_PROGRESSION:
                    minorRadius1 += self._tapered.minorRatio
                minorAxis1 = vector.setMagnitude(self._base._minorAxis, minorRadius1)
                baseC = [self._base._centre[c] + (n3 + 1) * arcLengthAlong * vector.normalise(self._base._alongAxis)[c]
                         for c in range(3)]
                self._basesCentres[n3 + 1] = baseC
                self._majorRadii.append(majorRadius1)
                self._minorRadii.append(minorRadius1)
                nx, nd1 = self.createCylinderBaseMesh2D(
                    baseC, majorAxis1, minorAxis1, elementsCountAround, majorRadius1)

        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        for n3 in range(self._elementsCountAlong + 1):
            for n in range(elementsCountAround + 1):
                n1, n2 = self._shield.convertRimIndex(n)
                tx, td1, td2, td3 = tnx[n3], tnd1[n3], tnd2[n3], tnd3[n3]
                btx[n3][n2][n1] = tx[n]
                if n2 > self._shield.elementsCountRim:  # regular rows
                    btd1[n3][n2][n1] = td1[n]
                    btd3[n3][n2][n1] = td3[n]
                if n2 >= 2:
                    btd3[n3][n2][n1] = vector.setMagnitude(minorAxis, vector.dotproduct(td3[n], minorAxis))
                else:  # around rim
                    btd1[n3][n2][n1] = td1[n]
                    btd3[n3][n2][n1] = td3[n]
                btd2[n3][n2][n1] = td2[n]

    def createMirrorCurve(self, n3):
        """
        generate coordinates and derivatives for the mirror curve
        :param n3: Index of along cylinder axis coordinates to use
        :return: Coordinates and derivatives for the mirror curve
        """
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        rcx = []
        tmdx = btx[n3][0][self._elementsCountAcrossMinor // 2]
        tmdd3 = btd3[n3][0][self._elementsCountAcrossMinor // 2]
        tmux = [
            0.5 * (btx[n3][self._elementsCountUp][0][c] + btx[n3][self._elementsCountUp][self._elementsCountAcrossMinor]
            [c]) for c in range(3)]
        rcx.append(tmdx)
        rcx.append(tmux)
        rcd3 = [vector.setMagnitude(tmdd3, -1), vector.setMagnitude(tmdd3, -1)]
        rscx, rscd1 = sampleCubicHermiteCurves(rcx, rcd3, self._shield.elementsCountUp, lengthFractionStart=1,
                                               arcLengthDerivatives=True)[0:2]

        # get d2, d3
        rscd2 = []
        rscd3 = []
        for n in range(len(rscx)):
            d3 = vector.normalise(
                [btx[n3][self._elementsCountUp][self._elementsCountAcrossMinor][c] - btx[n3][self._elementsCountUp][0]
                [c] for c in range(3)])
            d2 = vector.normalise(vector.crossproduct3(d3, rscd1[n]))
            rscd2.append(d2)
            rscd3.append(d3)

        return rscx, rscd1, rscd2, rscd3

    def createRegularRowCurves(self, n3, rscx, rscd1, rscd3):
        """
        generate curves along regular rows using the mirror curve obtained from createMirrorCurve.
        :param rscx: Coordinates of the nodes on the middle curve.
        :param rscd1: d1 derivatives of the nodes on the middle curve.
        :param rscd3: d3 derivatives of the nodes on the middle curve.
        :param n3: Index of along cylinder axis coordinates to use
        """
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        elementsCountRim = 0
        for n2 in range(elementsCountRim + 2, self._elementsCountUp + 1):
            btx[n3][n2], btd3[n3][n2], pe, pxi, psf = sampleCubicHermiteCurves(
                [btx[n3][n2][0], rscx[n2], btx[n3][n2][-1]],
                [vector.setMagnitude(btd3[n3][n2][0], -1.0), rscd3[n2], btd3[n3][n2][-1]],
                self._elementsCountAcrossMinor, lengthFractionStart=1, lengthFractionEnd=1, arcLengthDerivatives=True)
            btd1[n3][n2] = \
                interpolateSampleCubicHermite([[-btd1[n3][n2][0][c] for c in range(3)], rscd1[n2], btd1[n3][n2][-1]],
                                              [[0.0, 0.0, 0.0]] * 3, pe, pxi, psf)[0]
            btd3[n3][n2][0] = [-btd3[n3][n2][0][c] for c in range(3)]
            btd1[n3][n2][0] = [-btd1[n3][n2][0][c] for c in range(3)]

    def createRegularColumnCurves(self, n3):
        """
        up regular columns of shield: get d1, initial d3 below regular rows
        """
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        for n1 in range(2, self._elementsCountAcrossMinor - 1):
            tx, td1, pe, pxi, psf = sampleCubicHermiteCurves(
                [btx[n3][0][n1], btx[n3][2][n1]], [[-btd3[n3][0][n1][c] for c in range(3)], btd1[n3][2][n1]], 2,
                lengthFractionStart=1, arcLengthDerivatives=True)
            for n2 in range(3, self._elementsCountUp + 1):
                tx.append(btx[n3][n2][n1])
                td1.append(btd1[n3][n2][n1])
            td1 = smoothCubicHermiteDerivativesLine(tx, td1, fixStartDirection=True)
            td3 = \
                interpolateSampleCubicHermite([btd1[n3][0][n1], btd3[n3][2][n1]], [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[
                    0]
            for n2 in range(self._elementsCountUp + 1):
                if n2 < 2:
                    btx[n3][n2][n1] = tx[n2]
                    if n2 == 0:
                        btd3[n3][n2][n1] = [-td1[0][c] for c in range(3)]
                    else:
                        btd3[n3][n2][n1] = td3[n2]
                if n2 == 0:
                    btd1[n3][n2][n1] = td3[n2]
                else:
                    btd1[n3][n2][n1] = td1[n2]

    def smoothTriplePointsCurves(self, n3, n2b, n1b, m1a):
        """
        Smooth row and column curves passing triple points (i.e., row 1 and columns 1 and -2).
        """
        btx = self._shield.px
        btd1 = self._shield.pd1
        btd2 = self._shield.pd2
        btd3 = self._shield.pd3

        # smooth shield row 1
        btd3[n3][n2b][n1b:m1a] = smoothCubicHermiteDerivativesLine(btx[n3][n2b][n1b:m1a], btd3[n3][n2b][n1b:m1a])

        # smooth Shield columns 1, -2
        for n1 in [1, -2]:
            tx = []
            td1 = []
            for n2 in range(1, self._elementsCountUp + 1):
                tx.append(btx[n3][n2][n1])
                td1.append(btd1[n3][n2][n1])
            td1 = smoothCubicHermiteDerivativesLine(tx, td1)
            for n in range(self._elementsCountUp):
                btd1[n3][n + 1][n1] = td1[n]

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
        # get base d2 and next base x, d2
        for n2 in range(self._elementsCountUp + 1):
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
        for n2 in range(self._elementsCountUp + 1):
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
        assert (self._elementsCountAcrossMajor > 2), 'createCylinderMesh3d: Invalid number of up elements'
        assert (self._cylinderShape in [self._cylinderShape.CYLINDER_SHAPE_FULL,
                                        self._cylinderShape.CYLINDER_SHAPE_LOWER_HALF]), \
            'createCylinderMesh3d: Invalid cylinder mode.'
        plane = [-d for d in self._base._majorAxis] + [-vector.dotproduct(self._base._majorAxis, self._base._centre)]

        nodes = fieldModule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        mesh = fieldModule.findMeshByDimension(3)

        # create the base ellipse
        minorAxis = vector.setMagnitude(vector.crossproduct3(self._base._alongAxis, self._base._majorAxis),
                                        self._base._minorRadius)
        majorRadius = vector.magnitude(self._base._majorAxis)
        elementsCountAround = 2 * (self._elementsCountUp - 2) + self._elementsCountAcrossMinor

        # the bottom curve node coordinates and derivatives
        arcLengthAlong = vector.magnitude(self._base._alongAxis) / self._elementsCountAlong
        elementsCountRim = 0

        shieldMode = ShieldShape.SHIELD_SHAPE_FULL if self._cylinderShape is self._cylinderShape.CYLINDER_SHAPE_FULL \
            else ShieldShape.SHIELD_SHAPE_LOWER_HALF
        self._shield = ShieldMesh(self._elementsCountAcrossMinor, self._elementsCountAcrossMajor, elementsCountRim,
                                  None,
                                  self._elementsCountAlong, shieldMode,
                                  shieldType=ShieldRimDerivativeMode.SHIELD_RIM_DERIVATIVE_MODE_AROUND)

        # generate bases mesh along cylinder axis
        self.generateBasesMesh(majorRadius, elementsCountAround, arcLengthAlong, minorAxis)

        n3Count = 0 if self._cylinderType == CylinderType.CYLINDER_STRAIGHT else self._elementsCountAlong
        for n3 in range(n3Count + 1):
            rscx, rscd1, rscd2, rscd3 = self.createMirrorCurve(n3)
            self.createRegularRowCurves(n3, rscx, rscd1, rscd3)
            self.createRegularColumnCurves(n3)
            self._shield.getTriplePoints(n3)
            n1b = 1
            m1a = self._shield.elementsCountAcross
            m1b = m1a - 1
            m1c = m1a - 2
            n2b = 1
            self.smoothTriplePointsCurves(n3, n2b, n1b, m1a)
            self._shield.smoothDerivativesToTriplePoints(n3, fixAllDirections=True)

        for n3 in range(n3Count + 1):
            self.calculateD2Derivatives(n3, n3Count)

        if self._cylinderType == CylinderType.CYLINDER_TAPERED:
            self.smoothd2Derivatives()

        # The other bases.
        if self._cylinderType == CylinderType.CYLINDER_STRAIGHT:
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

        #################
        # Create nodes
        #################

        # nodeTemplate = nodes.createNodetemplate()
        # nodeTemplate.defineField(coordinates)
        # nodeTemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        # nodeTemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        # nodeTemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        # nodeTemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        self._startNodeIdentifier = nodeIdentifier
        nodeIdentifier = self._shield.generateNodes(fieldModule, coordinates, nodeIdentifier, plane)
        self._endNodeIdentifier = nodeIdentifier

        #################
        # Create elements
        #################

        elementIdentifier = max(1, getMaximumElementIdentifier(mesh) + 1)
        self._startElementIdentifier = elementIdentifier
        elementIdentifier = self._shield.generateElements(fieldModule, coordinates, elementIdentifier, [])
        self._endElementIdentifier = elementIdentifier

        if self._end is None:
            self._end = CylinderEnds(self._elementsCountAcrossMajor, self._elementsCountAcrossMinor,
                                     self._basesCentres[-1], self._shield.pd2[-1][0][1],
                                     vector.setMagnitude(self._base._majorAxis, self._majorRadii[-1]),
                                     self._minorRadii[-1])
        self.setEndsNodes()
