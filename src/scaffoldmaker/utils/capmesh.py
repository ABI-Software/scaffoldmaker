"""
Specialisation of Tube Network Mesh for building 3-D cap mesh.
"""

import math

from cmlibs.maths.vectorops import magnitude, sub, add, set_magnitude, normalize, rotate_vector_around_vector, cross, \
    angle, mult, div
from cmlibs.zinc.element import Element
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.eft_utils import determineCubicHermiteSerendipityEft
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLoop, smoothCubicHermiteDerivativesLine, \
    sampleCubicHermiteCurves, interpolateSampleCubicHermite, DerivativeScalingMode
from scaffoldmaker.utils.spheremesh import calculate_arc_length, local_to_global_coordinates, spherical_to_cartesian, \
    calculate_azimuth


class CapMesh:
    """
    Cap mesh generator.
    """

    def __init__(self, elementsCountAround, elementsCountCoreBoxMajor, elementsCountCoreBoxMinor,
                 elementsCountThroughShell, elementsCountTransition, networkPathParameters, isCap, isCore):
        """
        :param elementsCountAround: Number of elements around this segment.
        :param elementsCountCoreBoxMajor: Number of elements across core box major axis.
        :param elementsCountCoreBoxMinor: Number of elements across core box minor axis.
        :param elementsCountThroughShell: Number of elements between inner and outer tube.
        :param elementsCountTransition: Number of elements across transition zone between core box elements and
        rim elements.
        :param networkPathParameters: List containing path parameters of a tube network.
        :param isCap: List [startCap, endCap] with boolean values. True if the tube segment requires a cap at the
        start of a segment, or at the end of a segment, respectively. [True, True] if the segment requires cap at both
        ends.
        :param isCore: True for tube network with a solid core, False for regular tube network.
        """
        self._isCap = isCap
        self._isCore = isCore

        self._elementsCountAround = elementsCountAround
        self._elementsCountCoreBoxMajor = elementsCountCoreBoxMajor
        self._elementsCountCoreBoxMinor = elementsCountCoreBoxMinor
        self._elementsCountThroughShell = elementsCountThroughShell
        self._elementsCountTransition = elementsCountTransition

        self._networkPathParameters = networkPathParameters
        self._tubeBoxCoordinates = None  # tube box coordinates
        self._tubeTransitionCoordinates = None  # tube transition coordinates
        self._tubeShellCoordinates = None  # tube rim coordinates

        self._isStartCap = None
        self._generateData = None

        self._boxExtCoordinates = None
        # coordinates and derivatives for box nodes extended from the tube segment
        # list[startCap, endCap][x, d1, d2, d3][nAcrossMajor][nAcrossMinor]
        self._transitionExtCoordinates = None
        self._shellExtCoordinates = None
        # coordinates and derivatives for shell nodes extended from the tube segment
        # list[startCap, endCap][x, d1, d2, d3][nThroughWall][nAround]
        self._boxExtNodeIds = None
        self._rimExtNodeIds = None

        self._boxCoordinates = None
        # list[startCap, endCap][[x, d1, d2, d3][nAcrossMajor][nAcrossMinor]
        self._shellCoordinates = None
        # list[startCap, endCap][x, d1, d2, d3][nThroughWall][apex, rim][nAround] if the tube is without the solid core.
        # list[startCap, endCap][[x, d1, d2, d3][nThroughWall][nAcrossMajor][nAcrossMinor] if the tube is with the core.
        self._startCapNodeIds = None
        # capNodeIds that form the cap at the start of a tube segment.
        # list[nThroughWall][apex, rim] if the tube is without the core.
        # list[nThroughWall][nAcrossMajor][nAcrossMinor] if the tube is with the core.
        self._endCapNodeIds = None
        # capNodeIds that form the cap at the end of a tube segment.
        # list structure is identical to startCapNodeIds.
        self._startCapElementIds = None
        # elementIds that form the cap at the start of a tube segment.
        # list[nThroughWall][apex, rim] if the tube is without the core.
        # list[box, rim]: [box] and [rim] sublists have different structures.
        # [box][core, transition, shield][nAcrossMajor][nAcrossMinor]
        # [rim][base, shield][nAround]
        self._endCapElementIds = None
        # elementIds that form the cap at the end of a tube segment.
        self._startExtElementIds = None
        self._endExtElementIds = None

        # annotation groups created if core:
        self._coreGroup = None
        self._shellGroup = None

    def _extendTubeEnds(self):
        """
        Add additional tube sections with smaller element size along the tube at either ends of the tube with smaller
        D2 derivatives. This function is to minimise the effect of large difference in D2 derivatives between the cap
        mesh and the tube mesh.
        """
        isStartCap = self._isStartCap
        idx = 0 if isStartCap else -1

        layoutD1 = self._networkPathParameters[0][1][idx]
        ext = self._getExtensionLength()
        unitVector = normalize(layoutD1)
        signValue = -1 if isStartCap else 1

        coreBoxMajorNodesCount = self._elementsCountCoreBoxMajor + 1
        coreBoxMinorNodesCount = self._elementsCountCoreBoxMinor + 1

        boxCoordinates, transitionCoordinates, shellCoordinates = [], [], []
        for nx in range(4):
            boxCoordinates.append([])
            transitionCoordinates.append([])
            shellCoordinates.append([])
            if self._isCore:
                for m in range(coreBoxMajorNodesCount):
                    boxCoordinates[nx].append([])
                    for n in range(coreBoxMinorNodesCount):
                        boxCoordinates[nx][m].append([])
                if self._elementsCountTransition > 1:
                    for n3 in range(self._elementsCountTransition - 1):
                        transitionCoordinates[nx].append([])
            for n3 in range(self._elementsCountThroughShell + 1):
                shellCoordinates[nx].append([])

        if self._isCore:
            for m in range(coreBoxMajorNodesCount):
                xList, d2List = [], []
                x = self._tubeBoxCoordinates[0][idx][m]
                for n in range(coreBoxMinorNodesCount):
                    tx = add(x[n], set_magnitude(unitVector, ext * signValue))
                    td2 = mult(sub(tx, x[n]), signValue)
                    xList.append(tx)
                    d2List.append(td2)
                boxCoordinates[0][m] = xList
                boxCoordinates[1][m] = self._tubeBoxCoordinates[1][idx][m]
                boxCoordinates[2][m] = d2List
                boxCoordinates[3][m] = self._tubeBoxCoordinates[3][idx][m]

            if self._elementsCountTransition > 1:
                for n3 in range(self._elementsCountTransition - 1):
                    xList, d2List = [], []
                    x = self._tubeTransitionCoordinates[0][idx][n3]
                    for nx in range(self._elementsCountAround):
                        tx = add(x[nx], set_magnitude(unitVector, ext * signValue))
                        td2 = mult(sub(tx, x[nx]), signValue)
                        xList.append(tx)
                        d2List.append(td2)
                    transitionCoordinates[0][n3] = xList
                    transitionCoordinates[1][n3] = self._tubeTransitionCoordinates[1][idx][n3]
                    transitionCoordinates[2][n3] = d2List
                    transitionCoordinates[3][n3] = self._tubeTransitionCoordinates[3][idx][n3]

        for n3 in range(self._elementsCountThroughShell + 1):
            xList, d2List = [], []
            x = self._tubeShellCoordinates[0][idx][n3]
            for nx in range(self._elementsCountAround):
                tx = add(x[nx], set_magnitude(unitVector, ext * signValue))
                td2 = mult(sub(tx, x[nx]), signValue)
                xList.append(tx)
                d2List.append(td2)
            shellCoordinates[0][n3] = xList
            shellCoordinates[1][n3] = self._tubeShellCoordinates[1][idx][n3]
            shellCoordinates[2][n3] = d2List
            shellCoordinates[3][n3] = self._tubeShellCoordinates[3][idx][n3]

        if self._isCore:
            self._boxExtCoordinates = [None, None] if self._boxExtCoordinates is None else self._boxExtCoordinates
            self._boxExtCoordinates[idx] = boxCoordinates
            if self._elementsCountTransition > 1:
                self._transitionExtCoordinates = [None, None] if self._transitionExtCoordinates is None \
                    else self._transitionExtCoordinates
                self._transitionExtCoordinates[idx] = transitionCoordinates

        self._shellExtCoordinates = [None, None] if self._shellExtCoordinates is None else self._shellExtCoordinates
        self._shellExtCoordinates[idx] = shellCoordinates

    def _remapCapCoordinates(self):
        """
        Remap box and rim coordinates of the cap nodes based on the scale of tube extension.
        """
        isStartCap = self._isStartCap
        idx = 0 if isStartCap else -1

        layoutD1 = self._networkPathParameters[0][1][idx]
        ext = self._getExtensionLength()
        unitVector = normalize(layoutD1)
        signValue = -1 if isStartCap else 1

        coreBoxMajorNodesCount = self._getNodesCountCoreBoxMajor()
        coreBoxMinorNodesCount = self._getNodesCountCoreBoxMinor()
        nodesCountRim = self._getNodesCountRim()

        if self._isCore:
            for m in range(coreBoxMajorNodesCount):
                xList = []
                x = self._boxCoordinates[idx][0][m]
                for n in range(coreBoxMinorNodesCount):
                    tx = add(x[n], set_magnitude(unitVector, ext * signValue))
                    xList.append(tx)
                self._boxCoordinates[idx][0][m] = xList
            for n3 in range(nodesCountRim):
                for m in range(coreBoxMajorNodesCount):
                    xList = []
                    x = self._shellCoordinates[idx][0][n3][m]
                    for n in range(coreBoxMinorNodesCount):
                        tx = add(x[n], set_magnitude(unitVector, ext * signValue))
                        xList.append(tx)
                    self._shellCoordinates[idx][0][n3][m] = xList

    def _sampleCapCoordinatesWithoutCore(self):
        """
        Calculates coordinates and derivatives for the cap elements. It first calculates the coordinates for the apex
        nodes, and then calculates the coordinates for rim nodes on the shell surface.
        Used when the solid core is inactive.
        """
        self._shellCoordinates = [None, None] if self._shellCoordinates is None else self._shellCoordinates

        isStartCap = self._isStartCap
        idx = 0 if isStartCap else -1
        signValue = 1 if isStartCap else -1
        pathParameters = self._networkPathParameters[idx]

        outerRadius = self._getOuterShellRadius()
        shellThickness = self._getShellThickness()
        ext = self._getExtensionLength()
        centre = add(pathParameters[0][idx], set_magnitude(pathParameters[1][idx], ext * -signValue))
        outerWidth = outerLength = outerRadius
        innerWidth = innerLength = outerRadius - shellThickness

        elementLengthRatioEquatorApex = 1.0
        lengthRatio = 1.0

        bOuter = 2.0 / (1.0 + elementLengthRatioEquatorApex / lengthRatio)
        aOuter = 1.0 - bOuter
        bInner = 2.0 / (1.0 + elementLengthRatioEquatorApex / lengthRatio)
        aInner = 1.0 - bInner

        elementsCountUp = 2
        radiansPerElementAround = 2.0 * math.pi / self._elementsCountAround
        positionOuterArray = [(0, 0)] * elementsCountUp
        positionInnerArray = [(0, 0)] * elementsCountUp
        radiansUpOuterArray = [0] * elementsCountUp
        radiansUpInnerArray = [0] * elementsCountUp
        vector2OuterArray = [(0, 0)] * elementsCountUp
        vector2InnerArray = [(0, 0)] * elementsCountUp

        for n2 in range(2):
            xi = n2 * 2 / (2 * elementsCountUp)
            nxiOuter = aOuter * xi * xi + bOuter * xi
            dnxiOuter = 2.0 * aOuter * xi + bOuter
            radiansUpOuterArray[n2] = radiansUpOuter = nxiOuter * math.pi * 0.5
            dRadiansUpOuter = dnxiOuter * math.pi / (2 * elementsCountUp)
            cosRadiansUpOuter = math.cos(radiansUpOuter)
            sinRadiansUpOuter = math.sin(radiansUpOuter)
            positionOuterArray[n2] = [outerWidth * sinRadiansUpOuter, -outerLength * cosRadiansUpOuter]
            vector2OuterArray[n2] = (outerWidth * cosRadiansUpOuter * dRadiansUpOuter,
                                     outerLength * sinRadiansUpOuter * dRadiansUpOuter)

            nxiInner = aInner * xi * xi + bInner * xi
            dnxiInner = 2.0 * aInner * xi + bInner
            radiansUpInnerArray[n2] = radiansUpInner = nxiInner * math.pi * 0.5
            dRadiansUpInner = dnxiInner * math.pi / (2 * elementsCountUp)
            cosRadiansUpInner = math.cos(radiansUpInner)
            sinRadiansUpInner = math.sin(radiansUpInner)
            positionInnerArray[n2] = [innerWidth * sinRadiansUpInner, -innerLength * cosRadiansUpInner]
            vector2InnerArray[n2] = (innerWidth * cosRadiansUpInner * dRadiansUpInner,
                                     innerLength * sinRadiansUpInner * dRadiansUpInner)

        xList, d1List, d2List, d3List, rList = [], [], [], [], []
        elementsCountThroughShell = self._elementsCountThroughShell
        for n3 in range(elementsCountThroughShell + 1):
            for lst in [xList, d1List, d2List, d3List]:
                lst.append([])
            for n2 in range(2):
                n3_fraction = n3 / elementsCountThroughShell
                positionOuter = positionOuterArray[n2]
                positionInner = positionInnerArray[n2]
                position = [positionOuter[0] * n3_fraction + positionInner[0] * (1.0 - n3_fraction),
                            positionOuter[1] * n3_fraction + positionInner[1] * (1.0 - n3_fraction)]
                vector2Outer = vector2OuterArray[n2]
                vector2Inner = vector2InnerArray[n2]
                vector2 = [vector2Outer[0] * n3_fraction + vector2Inner[0] * (1.0 - n3_fraction),
                           vector2Outer[1] * n3_fraction + vector2Inner[1] * (1.0 - n3_fraction)]
                vector3 = [(positionOuter[0] - positionInner[0]) / elementsCountThroughShell,
                           (positionOuter[1] - positionInner[1]) / elementsCountThroughShell]
                # calculate coordinates
                if n2 == 0:  # apex
                    x = apex = add(pathParameters[0][idx],
                                   set_magnitude(pathParameters[1][idx], (position[1] - ext) * signValue))
                    d1 = set_magnitude(pathParameters[4][idx], vector2[0] * signValue)
                    d2 = set_magnitude(pathParameters[2][idx], vector2[0])
                    d3 = set_magnitude(pathParameters[1][idx], vector3[1] * signValue)
                    for lst, value in zip([xList, d1List, d2List, d3List], [x, d1, d2, d3]):
                        lst[-1].append(value)
                else:
                    refAxis = normalize(pathParameters[4][idx])
                    rotateAngle = n2 * (math.pi / 2) / elementsCountUp * signValue
                    radius = innerLength + (outerLength - innerLength) * n3 / elementsCountThroughShell
                    rList.append(radius)
                    tx = rotate_vector_around_vector(normalize(sub(apex, centre)), refAxis, -rotateAngle)
                    tx = set_magnitude(tx, radius)
                    for n1 in range(self._elementsCountAround):
                        radiansAround = n1 * radiansPerElementAround
                        cosRadiansAround = math.cos(radiansAround)
                        sinRadiansAround = math.sin(radiansAround)

                        rx = rotate_vector_around_vector(tx, pathParameters[1][idx], radiansAround)
                        rx = add(rx, centre)
                        d1 = [0.0, position[0] * -sinRadiansAround * radiansPerElementAround * signValue,
                              position[0] * cosRadiansAround * radiansPerElementAround * signValue]
                        d2 = [vector2[1], vector2[0] * cosRadiansAround,
                              vector2[0] * sinRadiansAround]
                        d3 = [vector3[1], vector3[0] * cosRadiansAround * signValue,
                              vector3[0] * sinRadiansAround * signValue]
                        for lst, value in zip([xList, d1List, d2List, d3List], [rx, d1, d2, d3]):
                            lst[-1].append(value)

        xCoordinates = [[] for _ in range(4)]
        for n, value in zip(range(4), [xList, d1List, d2List, d3List]):
            for n3 in range(self._elementsCountThroughShell + 1):
                xCoordinates[n].append([])
                xCoordinates[n][n3] = [value[n3][0], value[n3][1:]]
        self._shellCoordinates[idx] = xCoordinates

        # transform sphere to spheroid
        for n3 in range(elementsCountThroughShell + 1):
            radii = self._getTubeRadii(centre, n3, idx)
            oRadii = [1.0, rList[n3], rList[n3]]
            ratio = self._getRatioBetweenTwoRadii(radii, oRadii)
            self._sphereToSpheroid(n3, ratio, centre)
        # smooth derivatives
        for n3 in range(elementsCountThroughShell + 1):
            xList = self._shellCoordinates[idx][0]
            d1List = self._shellCoordinates[idx][1]
            sd1 = smoothCubicHermiteDerivativesLoop(xList[n3][1], d1List[n3][1])
            d1List[n3][1] = sd1
            for n1 in range(self._elementsCountAround):
                radiansAround = n1 * radiansPerElementAround
                x = xList[n3][1][n1]
                xStart, xEnd = xList[n3][0], self._shellExtCoordinates[idx][0][n3][n1]
                nx = [xStart, x, xEnd] if isStartCap else [xEnd, x, xStart]
                d2List = self._shellCoordinates[idx][2]
                d2Start = rotate_vector_around_vector(d2List[n3][0], pathParameters[1][idx], radiansAround)
                d2Start = set_magnitude(d2Start, magnitude(d2List[n3][0]) * signValue)
                d2End = set_magnitude(self._shellExtCoordinates[idx][2][n3][n1],
                                      magnitude(self._shellExtCoordinates[idx][1][n3][n1]))
                d2 = d2List[n3][1][n1]
                nd = [d2Start, d2, d2End] if isStartCap else [d2End, d2, d2Start]
                sd2 = smoothCubicHermiteDerivativesLine(nx, nd, fixStartDerivative=True, fixEndDerivative=True)
                d2List[n3][1][n1] = sd2[1]
        for n1 in range(self._elementsCountAround):
            xList = self._shellCoordinates[idx][0]
            d3List = self._shellCoordinates[idx][3]
            nx = [xList[n3][1][n1] for n3 in range(elementsCountThroughShell + 1)]
            nd = [d3List[n3][1][n1] for n3 in range(elementsCountThroughShell + 1)]
            sd3 = smoothCubicHermiteDerivativesLine(nx, nd)
            for n3 in range(elementsCountThroughShell + 1):
                d3List[n3][1][n1] = sd3[n3]

    def _sampleCapCoordinatesWithCore(self, s):
        """
        Blackbox function for calculating coordinates and derivatives for the cap elements.
        It first calculates the coordinates for shell nodes, then calculates for box nodes.
        nodes, and then calculates the coordinates for rim nodes on the shell surface.
        Used when the solid core is active.
        :param s: Index for isCap list. 0 indicates start cap and 1 indicates end cap.
        """
        self._isStartCap = isStartCap = True if self._isCap[0] and s == 0 else False
        idx = 0 if isStartCap else -1
        centre = self._networkPathParameters[0][0][idx]

        self._extendTubeEnds()  # extend tube end
        # shell nodes
        nodesCountRim = self._getNodesCountRim()
        for n3 in range(nodesCountRim):
            ox = self._getRimExtCoordinatesAround(n3)[0]
            radius = self._getRadius(ox)
            radii = self._getTubeRadii(centre, n3, idx)  # radii for spheroid
            oRadii = [1.0, radius, radius]  # original radii used to create the sphere
            ratio = self._getRatioBetweenTwoRadii(radii, oRadii)
            # ratio between original radii for the sphere and the new radii for spheroid
            self._calculateMajorAndMinorNodesCoordinates(n3, centre, ratio)
            self._calculateShellQuadruplePoints(n3, centre, radius)
            self._calculateShellRegularNodeCoordinates(n3, centre)
            self._sphereToSpheroid(n3, ratio, centre)
        self._determineShellDerivatives()
        # box nodes
        self._calculateBoxQuadruplePoints(centre)
        self._calculateBoxMajorAndMinorNodes()
        self._determineBoxDerivatives()

        self._remapCapCoordinates()
        self._smoothDerivatives()

    def _createShellCoordinatesList(self):
        """
        Creates an empty list for storing rim coordinates. Only applies when the solid core is active.
        """
        self._shellCoordinates = [] if self._shellCoordinates is None else self._shellCoordinates
        elementsCountRim = self._getElementsCountRim()
        for s in range(2):
            self._shellCoordinates.append([] if self._isCap[s] else None)
            if self._shellCoordinates[s] is not None:
                for nx in range(4):
                    self._shellCoordinates[s].append([])
                    for n3 in range(elementsCountRim):
                        self._shellCoordinates[s][nx].append([])
                        self._shellCoordinates[s][nx][n3] = [[] for _ in range(self._elementsCountCoreBoxMajor + 1)]
                        for m in range(self._elementsCountCoreBoxMajor + 1):
                            self._shellCoordinates[s][nx][n3][m] = \
                                [None for _ in range(self._elementsCountCoreBoxMinor + 1)]

    def _getOuterShellRadius(self):
        """
        Calculates the radius of an outer shell. It takes the average of a half-distance between two opposing nodes on
        the outer shell of a tube segment.
        :return: Radius of the cap shell.
        """
        ox = self._tubeShellCoordinates[0][0][-1] if self._isStartCap else self._tubeShellCoordinates[0][-1][-1]
        radii = [magnitude(sub(ox[i], ox[i + self._elementsCountAround // 2])) / 2 for i in
                 range(self._elementsCountAround // 2)]
        return sum(radii) / len(radii)

    def _getRadius(self, ox):
        """
        Calculates the radius of a shell. It takes the average of a half-distance between two opposing nodes around a
        tube segment.
        :param ox: Coordinates of shell nodes around a tube segment.
        :return: Radius of the cap shell.
        """
        radii = [magnitude(sub(ox[i], ox[i + self._elementsCountAround // 2])) / 2 for i in
                 range(self._elementsCountAround // 2)]
        return sum(radii) / len(radii)

    def _getShellThickness(self):
        """
        Calculates the thickness of a shell, based on the thickness of a tube segment at either ends.
        It takes the average of a distance between the outer and the inner node pair around the rim of a tube segment.
        :return: Thickness of the cap shell.
        """
        ix = self._tubeShellCoordinates[0][0][0] if self._isStartCap else self._tubeShellCoordinates[0][-1][0]
        ox = self._tubeShellCoordinates[0][0][-1] if self._isStartCap else self._tubeShellCoordinates[0][-1][-1]
        shellThicknesses = [magnitude(sub(ox[i], ix[i])) for i in range(self._elementsCountAround)]
        return sum(shellThicknesses) / len(shellThicknesses)

    def _getExtensionLength(self):
        """
        Calculates the length of extended tube segment. Currently set to half of the outer tube radius.
        :return: Length of extended tube segment.
        """
        outerRadius = self._getOuterShellRadius()
        return outerRadius / 2

    def _calculateMajorAndMinorNodesCoordinates(self, n3, centre, ratio):
        """
        Calculates coordinates and derivatives for major and minor axis nodes on the surface of a cap shell by rotating
        the major and minor axis nodes on the rim of a tube segment.
        :param n3: Node index from inner to outer rim.
        :param centre: Centre coordinates of a tube segment at either ends.
        :param ratio: List of ratios between original circular radii and new radii if the tube is non-circular.
        [x-axis, major axis, minor axis]. The values should equal 1.0 if the tube cross-section is circular.
        """
        idx = 0 if self._isStartCap else -1
        layoutD2 = self._networkPathParameters[0][2][idx]
        layoutD3 = self._networkPathParameters[0][4][idx]

        elementsCountAcrossMajor = self._elementsCountCoreBoxMajor + 2
        elementsCountAcrossMinor = self._elementsCountCoreBoxMinor + 2

        refAxis = normalize(layoutD2)
        rotateAngle = (math.pi / elementsCountAcrossMinor) if self._isStartCap else \
            -(math.pi / elementsCountAcrossMinor)
        minorAxisNodesCoordinates = [[], []]  # [startCap, endCap]
        n1 = self._elementsCountAround * 3 // 4
        ix = self._getTubeRimCoordinates(n1, idx, n3)
        for n in range(1, elementsCountAcrossMinor):
            for nx in [0, 1]:
                vi = sub(ix[nx], centre) if nx == 0 else mult(ix[nx], -1)
                vi = div(vi, ratio[2])
                vr = rotate_vector_around_vector(vi, refAxis, n * rotateAngle)
                vr = add(vr, centre) if nx == 0 else vr
                minorAxisNodesCoordinates[nx].append(vr)

        refAxis = normalize(layoutD3)
        rotateAngle = (math.pi / elementsCountAcrossMajor) if self._isStartCap else \
            -(math.pi / elementsCountAcrossMajor)
        majorAxisNodesCoordinates = [[], []]  # [startCap, endCap]
        ix = self._getTubeRimCoordinates(0, idx, n3)
        for m in range(1, elementsCountAcrossMajor):
            for nx in [0, 1]:  # [x, d1]
                vi = sub(ix[nx], centre) if nx == 0 else ix[nx]
                vi = div(vi, ratio[1])
                vr = rotate_vector_around_vector(vi, refAxis, m * rotateAngle)
                vr = add(vr, centre) if nx == 0 else vr
                majorAxisNodesCoordinates[nx].append(vr)

        midMajorIndex = elementsCountAcrossMajor // 2 - 1
        midMinorIndex = elementsCountAcrossMinor // 2 - 1
        for n in range(elementsCountAcrossMinor - 1):
            for i in [0, 1]:
                nx = [0, 1][i]
                if self._shellCoordinates[idx][nx][n3][midMajorIndex][n] is None:
                    self._shellCoordinates[idx][nx][n3][midMajorIndex][n] = minorAxisNodesCoordinates[i][n]
        for m in range(elementsCountAcrossMajor - 1):
            for i in [0, 1]:
                nx = [0, 2][i]
                if self._shellCoordinates[idx][nx][n3][m][midMinorIndex] is None:
                    self._shellCoordinates[idx][nx][n3][m][midMinorIndex] = majorAxisNodesCoordinates[i][m]

        # derivatives
        for n in range(elementsCountAcrossMinor - 1):
            if self._shellCoordinates[idx][2][n3][midMajorIndex][n] is None:
                self._shellCoordinates[idx][2][n3][midMajorIndex][n] = majorAxisNodesCoordinates[1][midMajorIndex]
        tx = self._shellCoordinates[idx][0][n3][midMajorIndex]
        td2 = self._shellCoordinates[idx][2][n3][midMajorIndex]
        self._shellCoordinates[idx][2][n3][midMajorIndex] = smoothCubicHermiteDerivativesLine(tx, td2)

        for m in range(elementsCountAcrossMajor - 1):
            if self._shellCoordinates[idx][1][n3][m][midMinorIndex] is None:
                self._shellCoordinates[idx][1][n3][m][midMinorIndex] = minorAxisNodesCoordinates[1][midMinorIndex]
        tx = [self._shellCoordinates[idx][0][n3][m][midMinorIndex] for m in range(elementsCountAcrossMajor - 1)]
        td1 = [self._shellCoordinates[idx][1][n3][m][midMinorIndex] for m in range(elementsCountAcrossMajor - 1)]
        sd1 = smoothCubicHermiteDerivativesLine(tx, td1)
        for m in range(elementsCountAcrossMajor - 1):
            self._shellCoordinates[idx][1][n3][m][midMinorIndex] = sd1[m]

    def _calculateShellQuadruplePoints(self, n3, centre, radius):
        """
        Calculate coordinates and derivatives of the quadruple point on the surface, where 3 hex elements merge.
        :param n3: Node index from inner to outer rim.
        :param centre: Centre coordinates of a tube segment at either ends.
        :param radius: Shell radius.
        """
        idx = 0 if self._isStartCap else -1

        layoutD1 = self._networkPathParameters[0][1][idx]
        layoutD2 = self._networkPathParameters[0][2][idx]
        layoutD3 = self._networkPathParameters[0][4][idx]

        axesList = [[mult(layoutD1, -1), layoutD2, mult(layoutD3, -1)],
                    [layoutD2, mult(layoutD1, -1), layoutD3],
                    [mult(layoutD2, -1), mult(layoutD1, -1), mult(layoutD3, -1)],
                    [mult(layoutD1, -1), mult(layoutD2, -1), layoutD3]]

        elementsCountUp = 2
        elementsCountAcrossMajor = self._elementsCountCoreBoxMajor + 2
        elementsCountAcrossMinor = self._elementsCountCoreBoxMinor + 2
        signValue = 1 if self._isStartCap else -1
        for counter, (m, n) in enumerate([(0, 0), (0, -1), (-1, 0), (-1, -1)]):
            elementsCount = (
                [elementsCountUp, elementsCountAcrossMajor // 2, elementsCountAcrossMinor // 2] if m == n else
                [elementsCountAcrossMajor // 2, elementsCountUp, elementsCountAcrossMinor // 2]
            )

            radiansPerElementAroundEllipse12 = math.pi / (2 * (elementsCount[0] + elementsCount[1] - 2))
            radiansPerElementAroundEllipse13 = math.pi / (2 * (elementsCount[0] + elementsCount[2] - 2))

            theta_2 = (elementsCount[2] - 1) * radiansPerElementAroundEllipse13
            theta_3 = (elementsCount[1] - 1) * radiansPerElementAroundEllipse12
            phi_3 = calculate_azimuth(theta_3, theta_2)

            local_x = spherical_to_cartesian(radius, theta_3, phi_3)
            c = counter if self._isStartCap else -(counter + 1)
            axes = [mult(axis, signValue) for axis in axesList[c]]
            x = local_to_global_coordinates(local_x, axes, centre)

            self._shellCoordinates[idx][0][n3][m][n] = x

    def _calculateShellRegularNodeCoordinates(self, n3, centre):
        """
        Calculate coordinates and derivatives of all other shell nodes on the cap surface.
        :param n3: Node index from inner to outer rim.
        :param centre: Centre coordinates of a tube segment at either ends.
        """
        idx = 0 if self._isStartCap else -1
        elementsCountAcrossMajor = self._elementsCountCoreBoxMajor + 2
        elementsCountAcrossMinor = self._elementsCountCoreBoxMinor + 2
        midMajorIndex = elementsCountAcrossMajor // 2 - 1
        midMinorIndex = elementsCountAcrossMinor // 2 - 1

        elementsOut = self._elementsCountCoreBoxMinor // 2
        for m in range(self._elementsCountCoreBoxMajor + 1):
            for n in [0, -1]:
                x1 = self._shellCoordinates[idx][0][n3][m][n]
                x2 = self._shellCoordinates[idx][0][n3][m][midMinorIndex]
                if x1 is None:
                    continue
                nx, nd1 = (self._sampleCurvesOnSphere(x1, x2, centre, elementsOut) if n == 0 else
                           self._sampleCurvesOnSphere(x2, x1, centre, elementsOut))
                start, end = (
                    (n + 1, midMinorIndex) if n == 0 else (midMinorIndex + 1, self._elementsCountCoreBoxMinor))
                for c in range(start, end):
                    idx_c = c % (len(nx) - 1)
                    self._shellCoordinates[idx][0][n3][m][c] = nx[idx_c]
                    self._shellCoordinates[idx][1][n3][m][c] = nd1[idx_c]
                    self._shellCoordinates[idx][2][n3][m][c] = [0, 0, 0]

        elementsOut = self._elementsCountCoreBoxMajor // 2
        for n in range(self._elementsCountCoreBoxMinor + 1):
            for m in [0, -1]:
                x1 = self._shellCoordinates[idx][0][n3][m][n]
                x2 = self._shellCoordinates[idx][0][n3][midMajorIndex][n]
                if x1 is None:
                    continue
                nx, nd2 = (self._sampleCurvesOnSphere(x1, x2, centre, elementsOut) if m == 0 else
                           self._sampleCurvesOnSphere(x2, x1, centre, elementsOut))
                start, end = (
                    (m + 1, midMajorIndex) if m == 0 else (midMajorIndex + 1, self._elementsCountCoreBoxMajor))
                for c in range(start, end):
                    idx_c = c % (len(nx) - 1)
                    self._shellCoordinates[idx][0][n3][c][n] = nx[idx_c]
                    self._shellCoordinates[idx][1][n3][c][n] = [0, 0, 0]
                    self._shellCoordinates[idx][2][n3][c][n] = nd2[idx_c]

    def _sphereToSpheroid(self, n3, ratio, centre):
        """
        Transform the sphere to ellipsoid using the radius in each direction.
        :param n3: Node index from inner to outer rim.
        :param ratio: List of ratios between original circular radii and new radii if the tube is non-circular.
        [x-axis, major axis, minor axis]. The values should equal 1.0 if the tube cross-section is circular.
        :param centre: Centre coordinates of a tube segment at either ends.
        """
        idx = 0 if self._isStartCap else -1

        # rotation angles to use in case when the cap is tilted from xyz axes.
        layoutD2 = normalize(self._networkPathParameters[0][2][idx])
        layoutD3 = normalize(self._networkPathParameters[0][4][idx])
        thetaD2 = angle(layoutD2, [0.0, 1.0, 0.0])
        thetaD3 = angle(layoutD3, [0.0, 0.0, 1.0])

        mCount = self._elementsCountCoreBoxMajor + 1 if self._isCore else 1
        nCount = self._elementsCountCoreBoxMinor + 1 if self._isCore else self._elementsCountAround
        # process shell coordinates
        for m in range(mCount):
            mp = m if self._isCore else 1
            for n in range(nCount):
                btx = self._shellCoordinates[idx][0][n3][mp][n]
                btx = sub(btx, centre)
                # apply forward transformations
                for vec, theta in [(layoutD3, thetaD2), (layoutD2, thetaD3)]:
                    btx = rotate_vector_around_vector(btx, vec, theta)
                # scale by ratios
                btx = [ratio[c] * btx[c] for c in range(3)]
                # apply inverse transformations
                for vec, theta in [(layoutD2, -thetaD3), (layoutD3, -thetaD2)]:
                    btx = rotate_vector_around_vector(btx, vec, theta)
                # update shell coordinates
                self._shellCoordinates[idx][0][n3][mp][n] = add(btx, centre)

    def _getRatioBetweenTwoRadii(self, radii, oRadii):
        """
        Calculates the ratio between the original radius of a sphere and the new radius of an ellipsoid.
        :param radii: List of new radius in each direction.
        :param oRadii: List of original radius in each direction.
        :return: List of ratio between two radii in x, y, and z-direction.
        """
        return [radii[c] / oRadii[c] for c in range(3)]

    def _getTubeRadii(self, centre, n3, idx):
        """
        Calculates the radius of a tube segment in major axis and minor axis.
        :param centre:
        :param n3: Node index from inner to outer rim.
        :param idx: 0 if calculating for the start segment, -1 if calculating for the end segment.
        :return: List of radii in major and minor axes. The radius in x-direction is set to 1.0 by default because the
        radius in this direction is constant.
        """
        n1m, n1n = 0, self._elementsCountAround // 4
        ixm, ixn = (self._getTubeRimCoordinates(n, idx, n3)[0] for n in [n1m, n1n])
        majorRadius, minorRadius = (magnitude(sub(coord, centre)) for coord in [ixm, ixn])
        # if majorRadius > minorRadius:
        #     xRadius = majorRadius / minorRadius
        # elif majorRadius < minorRadius:
        #     xRadius = minorRadius / majorRadius
        # else:
        #     xRadius = 1.0
        return [1.0, majorRadius, minorRadius]

    def _determineShellDerivatives(self):
        """
        Compute d1, d2, and d3 derivatives for the shell nodes.
        """
        idx = 0 if self._isStartCap else -1
        nodesCountRim = self._getNodesCountRim()
        for n3 in range(nodesCountRim):
            for m in [0, -1]:
                for n in [0, -1]:
                    # initialise derivatives for quadruple points
                    mp = m + 1 if m == 0 else m - 1
                    np = n + 1 if n == 0 else n - 1
                    d1 = sub(self._shellCoordinates[idx][0][n3][mp][n], self._shellCoordinates[idx][0][n3][m][n])
                    self._shellCoordinates[idx][1][n3][m][n] = mult(d1, -1) if m == -1 else d1
                    d2 = sub(self._shellCoordinates[idx][0][n3][m][np], self._shellCoordinates[idx][0][n3][m][n])
                    self._shellCoordinates[idx][2][n3][m][n] = mult(d2, -1) if n == -1 else d2

        for n3 in range(nodesCountRim):
            for m in range(self._elementsCountCoreBoxMajor + 1):
                signValue = 1 if self._isStartCap else -1
                tx = self._shellCoordinates[idx][0][n3][m]
                td2 = self._shellCoordinates[idx][2][n3][m]
                sd2 = smoothCubicHermiteDerivativesLine(tx, td2)
                for n in range(self._elementsCountCoreBoxMinor + 1):
                    self._shellCoordinates[idx][2][n3][m][n] = mult(sd2[n], signValue)
            for n in range(self._elementsCountCoreBoxMinor + 1):
                tx = [self._shellCoordinates[idx][0][n3][m][n] for m in range(self._elementsCountCoreBoxMajor + 1)]
                td1 = [self._shellCoordinates[idx][1][n3][m][n] for m in range(self._elementsCountCoreBoxMajor + 1)]
                sd1 = smoothCubicHermiteDerivativesLine(tx, td1)
                for m in range(self._elementsCountCoreBoxMajor + 1):
                    self._shellCoordinates[idx][1][n3][m][n] = sd1[m]

        elementsCountRim = self._elementsCountThroughShell + self._elementsCountTransition - 1
        for m in range(self._elementsCountCoreBoxMajor + 1):
            for n in range(self._elementsCountCoreBoxMinor + 1):
                otx = self._shellCoordinates[idx][0][-1][m][n]
                itx = self._shellCoordinates[idx][0][0][m][n]
                shellFactor = 1.0 / elementsCountRim
                sd3 = mult(sub(otx, itx), shellFactor)
                for n3 in range(nodesCountRim):
                    self._shellCoordinates[idx][3][n3][m][n] = sd3

    def _calculateBoxQuadruplePoints(self, centre):
        """
        Calculate coordinates and derivatives of the quadruple point for the box elements, where 3 hex elements merge.
        :param centre: Centre coordinates of a tube segment at either ends.
        """
        idx = 0 if self._isStartCap else -1
        capBoxCoordinates = []
        for nx in range(4):
            capBoxCoordinates.append([])
            capBoxCoordinates[nx] = [[] for _ in range(self._elementsCountCoreBoxMajor + 1)]
            for m in range(self._elementsCountCoreBoxMajor + 1):
                capBoxCoordinates[nx][m] = [None for _ in range(self._elementsCountCoreBoxMinor + 1)]

        boxCoordinates = self._tubeBoxCoordinates[0][idx]
        rimCoordinates = self._tubeTransitionCoordinates[0][idx][0] if self._elementsCountTransition > 1 \
            else self._tubeShellCoordinates[0][idx][0]
        capCoordinates = self._shellCoordinates[idx][0][0]

        elementsCountAround = self._elementsCountAround
        nodesCountAcrossMinorHalf = len(boxCoordinates[0]) // 2
        triplePointIndexesList = []
        for n in range(0, elementsCountAround, elementsCountAround // 2):
            triplePointIndexesList.append((n - nodesCountAcrossMinorHalf) % elementsCountAround)
            triplePointIndexesList.append(n + nodesCountAcrossMinorHalf)
        triplePointIndexesList[-2], triplePointIndexesList[-1] = triplePointIndexesList[-1], triplePointIndexesList[-2]

        for counter, (m, n) in enumerate([(0, 0), (0, -1), (-1, 0), (-1, -1)]):
            tpIndex = triplePointIndexesList[counter]
            x1, x2, x3 = boxCoordinates[m][n], rimCoordinates[tpIndex], capCoordinates[m][n]
            ts = magnitude(sub(x1, x2))
            ra = sub(x3, centre)
            radius = magnitude(ra)
            local_x = mult(ra, (1 - ts / radius))
            capBoxCoordinates[0][m][n] = add(local_x, centre)
            capBoxCoordinates[1][m][n] = self._tubeBoxCoordinates[1][idx][m][n]
            capBoxCoordinates[3][m][n] = self._tubeBoxCoordinates[3][idx][m][n]

        if self._boxCoordinates is None:
            self._boxCoordinates = [None] * 2
        self._boxCoordinates[idx] = capBoxCoordinates

    def _calculateBoxMajorAndMinorNodes(self):
        """
        Calculate coordinates and derivatives for box nodes along the central major and minor axes.
        """
        idx = 0 if self._isStartCap else -1
        midMajorIndex = self._elementsCountCoreBoxMajor // 2
        midMinorIndex = self._elementsCountCoreBoxMinor // 2

        # box side nodes
        for m in [0, -1]:
            nx = [self._boxCoordinates[idx][0][m][n] for n in [0, -1]]
            nd1 = [self._boxCoordinates[idx][1][m][n] for n in [0, -1]]
            nd3 = [self._boxCoordinates[idx][3][m][n] for n in [0, -1]]
            tx, td3, pe, pxi, psf = sampleCubicHermiteCurves(nx, nd3, self._elementsCountCoreBoxMinor, arcLengthDerivatives=True)
            td1 = interpolateSampleCubicHermite(nd1, [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
            for n in range(1, self._elementsCountCoreBoxMinor):
                self._boxCoordinates[idx][0][m][n] = tx[n]
                self._boxCoordinates[idx][1][m][n] = td1[n]
                self._boxCoordinates[idx][3][m][n] = td3[n]
        for n in [0, -1]:
            nx = [self._boxCoordinates[idx][0][m][n] for m in [0, -1]]
            nd1 = [self._boxCoordinates[idx][1][m][n] for m in [0, -1]]
            nd3 = [self._boxCoordinates[idx][3][m][n] for m in [0, -1]]
            tx, td1, pe, pxi, psf = sampleCubicHermiteCurves(nx, nd1, self._elementsCountCoreBoxMajor, arcLengthDerivatives=True)
            td3 = interpolateSampleCubicHermite(nd3, [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
            for m in range(1, self._elementsCountCoreBoxMajor):
                self._boxCoordinates[idx][0][m][n] = tx[m]
                self._boxCoordinates[idx][1][m][n] = td1[m]
                self._boxCoordinates[idx][3][m][n] = td3[m]

        # box major and minor nodes
        nx = [self._boxCoordinates[idx][0][midMajorIndex][n] for n in [0, -1]]
        nd1 = [self._boxCoordinates[idx][1][midMajorIndex][n] for n in [0, -1]]
        nd3 = [self._boxCoordinates[idx][3][midMajorIndex][n] for n in [0, -1]]
        tx, td3, pe, pxi, psf = sampleCubicHermiteCurves(nx, nd3, self._elementsCountCoreBoxMinor, arcLengthDerivatives=True)
        td1 = interpolateSampleCubicHermite(nd1, [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
        for n in range(1, self._elementsCountCoreBoxMinor):
            self._boxCoordinates[idx][0][midMajorIndex][n] = tx[n]
            self._boxCoordinates[idx][1][midMajorIndex][n] = td1[n]
            self._boxCoordinates[idx][3][midMajorIndex][n] = td3[n]

        nx = [self._boxCoordinates[idx][0][m][midMinorIndex] for m in [0, -1]]
        nd1 = [self._boxCoordinates[idx][1][m][midMinorIndex] for m in [0, -1]]
        tx, td1, pe, pxi, psf = sampleCubicHermiteCurves(nx, nd1, self._elementsCountCoreBoxMajor, arcLengthDerivatives=True)
        td3 = interpolateSampleCubicHermite(nd3, [[0.0, 0.0, 0.0]] * 2, pe, pxi, psf)[0]
        for m in range(1, self._elementsCountCoreBoxMajor):
            self._boxCoordinates[idx][0][m][midMinorIndex] = tx[m]
            self._boxCoordinates[idx][1][m][midMinorIndex] = td1[m]
            self._boxCoordinates[idx][3][m][midMinorIndex] = td3[m]

        # remaining nodes
        for m in range(self._elementsCountCoreBoxMajor):
            for n in range(self._elementsCountCoreBoxMinor):
                if self._boxCoordinates[idx][0][m][n] is None:
                    nx = [self._boxCoordinates[idx][0][i][n] for i in [0, midMajorIndex, -1]]
                    nd1 = [self._boxCoordinates[idx][1][i][n] for i in [0, midMajorIndex, -1]]
                    nd3 = [self._boxCoordinates[idx][3][i][n] for i in [0, midMajorIndex, -1]]
                    tx, td1, pe, pxi, psf = sampleCubicHermiteCurves(nx, nd1, self._elementsCountCoreBoxMajor, arcLengthDerivatives=True)
                    td3 = interpolateSampleCubicHermite(nd3, [[0.0, 0.0, 0.0]] * 3, pe, pxi, psf)[0]
                    for mi in range(1, self._elementsCountCoreBoxMajor):
                        self._boxCoordinates[idx][0][mi][n] = tx[mi]
                        self._boxCoordinates[idx][1][mi][n] = td1[mi]
                        self._boxCoordinates[idx][3][mi][n] = td3[mi]

        # smooth derivatives
        for m in range(self._elementsCountCoreBoxMajor + 1):
            nx = self._boxCoordinates[idx][0][m]
            nd3 = self._boxCoordinates[idx][3][m]
            sd3 = smoothCubicHermiteDerivativesLine(nx, nd3)
            for n in range(self._elementsCountCoreBoxMinor + 1):
                self._boxCoordinates[idx][3][m][n] = sd3[n]
        for n in range(self._elementsCountCoreBoxMinor + 1):
            nx = [self._boxCoordinates[idx][0][m][n] for m in range(self._elementsCountCoreBoxMajor + 1)]
            nd1 = [self._boxCoordinates[idx][1][m][n] for m in range(self._elementsCountCoreBoxMajor + 1)]
            sd1 = smoothCubicHermiteDerivativesLine(nx, nd1)
            for m in range(self._elementsCountCoreBoxMajor + 1):
                self._boxCoordinates[idx][1][m][n] = sd1[m]

    def _determineBoxDerivatives(self):
        """
        Calculate d2 derivatives of box nodes.
        """
        idx = 0 if self._isStartCap else -1
        signValue = 1 if self._isStartCap else -1
        for m in range(self._elementsCountCoreBoxMajor + 1):
            for n in range(self._elementsCountCoreBoxMinor + 1):
                otx = self._tubeBoxCoordinates[0][idx][m][n]
                itx = self._boxCoordinates[idx][0][m][n]
                d2 = mult(sub(otx, itx), signValue)
                self._boxCoordinates[idx][2][m][n] = d2

    def _smoothDerivatives(self):
        """
        Smooths derivatives to eliminate zero Jacobian contours at the cap-tube joint.
        """
        nodesCountCoreBoxMajor = self._getNodesCountCoreBoxMajor()
        nodesCountCoreBoxMinor = self._getNodesCountCoreBoxMinor()
        nodesCountRim = self._getNodesCountRim()
        nloop = self._getNodesCountRim()

        idx = 0 if self._isStartCap else -1
        if self._isCap[idx]:
            boxCoordinates = self._boxCoordinates[idx]
            boxExtCoordinates = self._boxExtCoordinates[idx]
            shellCoordinates = self._shellCoordinates[idx]
            tubeBoxCoordinates = self._tubeBoxCoordinates
            boxBoundaryNodeToBoxId = self._createBoxBoundaryNodeIdsList()

            # smooth derivatives along d2 for box
            for m in range(nodesCountCoreBoxMajor):
                for n in range(nodesCountCoreBoxMinor):
                    nx = [boxCoordinates[0][m][n], boxExtCoordinates[0][m][n], tubeBoxCoordinates[0][idx][m][n]]
                    nd2 = [boxCoordinates[2][m][n], boxExtCoordinates[2][m][n], tubeBoxCoordinates[2][idx][m][n]]
                    sd2 = smoothCubicHermiteDerivativesLine(nx, nd2, fixAllDirections=True,
                                                            magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
                    boxCoordinates[2][m][n] = sd2[0]
                    boxExtCoordinates[2][m][n] = sd2[1]

            # smooth derivatives along d3 for shell
            for m in range(nodesCountCoreBoxMajor):
                for n in range(nodesCountCoreBoxMinor):
                    bx = boxCoordinates[0][m][n]
                    bd3 = sub(shellCoordinates[0][0][m][n], bx)
                    nx = [bx] + [shellCoordinates[0][n3][m][n] for n3 in range(nloop)]
                    nd3 = [bd3] + [shellCoordinates[3][n3][m][n] for n3 in range(nloop)]
                    sd3 = smoothCubicHermiteDerivativesLine(nx, nd3, fixAllDirections=True,
                                                            magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
                    for n3 in range(nloop):
                        shellCoordinates[3][n3][m][n] = sd3[n3 + 1]

            # smooth transition/shell ext coordinates along d2
            for n2 in range(self._elementsCountAround):
                m, n = boxBoundaryNodeToBoxId[n2]
                for n3 in range(nodesCountRim - 1):
                    tx = self._getRimExtCoordinatesAround(n3)[0][n2]
                    bx = shellCoordinates[0][n3][m][n]
                    bd2 = sub(bx, tx)
                    nx = [bx, tx]
                    nd2 = [bd2, self._getRimExtCoordinatesAround(n3)[2][n2]]
                    sd2 = smoothCubicHermiteDerivativesLine(nx, nd2, fixAllDirections=True,
                                                            magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
                    if self._elementsCountTransition - 1 > n3:
                        self._transitionExtCoordinates[idx][2][n3][n2] = sd2[1]
                    else:
                        n3p = n3 - 1 - self._elementsCountTransition if self._elementsCountTransition > 1 else n3
                        self._shellExtCoordinates[idx][2][n3p][n2] = sd2[1]

            # smooth shell ext coordinates along d3
            for n2 in range(self._elementsCountAround):
                m, n = boxBoundaryNodeToBoxId[n2]
                bx = boxCoordinates[0][m][n]
                tx = [self._getRimExtCoordinatesAround(n3)[0][n2] for n3 in range(nodesCountRim)]
                bd3 = sub(bx, tx[0])
                td3 = [self._getRimExtCoordinatesAround(n3)[3][n2] for n3 in range(nodesCountRim)]
                nx = [bx] + tx
                nd3 = [bd3] + td3
                sd3 = smoothCubicHermiteDerivativesLine(nx, nd3, fixAllDirections=True,
                                                        magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
                self._shellExtCoordinates[idx][3][0][n2] = sd3[1]

    def _createBoxBoundaryNodeIdsList(self):
        """
        Creates a list (in a circular format similar to other rim node id lists) of core box node ids that are
        located at the boundary of the core. This list is used to easily stitch inner rim nodes with box nodes.
        Used specifically for solid core at the junction.
        :return: A list of box node ids stored in a circular format, and a lookup list that translates indexes used in
        boxBoundaryNodeIds list to indexes that can be used in boxCoordinates list.
        """
        boxBoundaryNodeToBoxId = []
        elementsCountCoreBoxMajor = self._elementsCountCoreBoxMajor
        elementsCountCoreBoxMinor = self._elementsCountCoreBoxMinor
        coreBoxMajorNodesCount = elementsCountCoreBoxMajor + 1
        coreBoxMinorNodesCount = elementsCountCoreBoxMinor + 1

        for n3 in range(coreBoxMajorNodesCount):
            if n3 == 0 or n3 == coreBoxMajorNodesCount - 1:
                n1List = list(range(coreBoxMinorNodesCount)) if n3 == 0 else (
                    list(range(coreBoxMinorNodesCount - 1, -1, -1)))
                for n1 in n1List:
                    boxBoundaryNodeToBoxId.append([n3, n1])
            else:
                for n1 in [-1, 0]:
                    boxBoundaryNodeToBoxId.append([n3, n1])

        start = elementsCountCoreBoxMajor - 2
        idx = elementsCountCoreBoxMinor + 2
        for n in range(int(start), -1, -1):
            boxBoundaryNodeToBoxId.append(boxBoundaryNodeToBoxId.pop(idx + 2 * n))

        nloop = elementsCountCoreBoxMinor // 2
        for _ in range(nloop):
            boxBoundaryNodeToBoxId.insert(len(boxBoundaryNodeToBoxId), boxBoundaryNodeToBoxId.pop(0))

        return boxBoundaryNodeToBoxId

    def _createBoundaryNodeIdsList(self, nodeIds):
        """
        Creates a list (in a circular format similar to other rim node id lists) of box node ids that are
        located at the boundary of the box component.
        This list is used to easily stitch inner rim nodes with box nodes.
        :param nodeIds: List of box node ids to be rearranged.
        :return: A list of box node ids stored in a circular format, and a lookup list that translates indexes used in
        boxBoundaryNodeIds list to indexes that can be used in boxCoordinates list.
        """
        capBoundaryNodeIds, capBoundaryNodeToBoxId = [], []
        boxElementsCountRow = self._elementsCountCoreBoxMajor + 1
        boxElementsCountColumn = self._elementsCountCoreBoxMinor + 1
        for n3 in range(boxElementsCountRow):
            if n3 == 0 or n3 == boxElementsCountRow - 1:
                ids = nodeIds[n3] if n3 == 0 else nodeIds[n3][::-1]
                n1List = list(range(boxElementsCountColumn)) if n3 == 0 else (
                    list(range(boxElementsCountColumn - 1, -1, -1)))
                capBoundaryNodeIds += [ids[c] for c in range(boxElementsCountColumn)]
                for n1 in n1List:
                    capBoundaryNodeToBoxId.append([n3, n1])
            else:
                for n1 in [-1, 0]:
                    capBoundaryNodeIds.append(nodeIds[n3][n1])
                    capBoundaryNodeToBoxId.append([n3, n1])

        start = self._elementsCountCoreBoxMajor - 2
        idx = self._elementsCountCoreBoxMinor + 2
        for n in range(int(start), -1, -1):
            capBoundaryNodeIds.append(capBoundaryNodeIds.pop(idx + 2 * n))
            capBoundaryNodeToBoxId.append(capBoundaryNodeToBoxId.pop(idx + 2 * n))

        nloop = self._elementsCountCoreBoxMinor // 2
        for _ in range(nloop):
            capBoundaryNodeIds.insert(len(capBoundaryNodeIds), capBoundaryNodeIds.pop(0))
            capBoundaryNodeToBoxId.insert(len(capBoundaryNodeToBoxId), capBoundaryNodeToBoxId.pop(0))

        return capBoundaryNodeIds, capBoundaryNodeToBoxId

    def _getBoxCoordinates(self, m, n):
        """
        :param m: Index along the major axis.
        :param n: Index along the minor axis.
        :return: x[], d1[], d2[], and d3[] of box nodes in the cap mesh.
        """
        idx = 0 if self._isStartCap else -1
        return [self._boxCoordinates[idx][0][m][n],
                self._boxCoordinates[idx][1][m][n],
                self._boxCoordinates[idx][2][m][n],
                self._boxCoordinates[idx][3][m][n]]

    def _getBoxExtCoordinates(self, m, n):
        """
        :param m: Index along the major axis.
        :param n: Index along the minor axis.
        :return: x[], d1[], d2[], and d3[] of box nodes extended from the tube segment.
        """
        idx = 0 if self._isStartCap else -1
        return [self._boxExtCoordinates[idx][0][m][n],
                self._boxExtCoordinates[idx][1][m][n],
                self._boxExtCoordinates[idx][2][m][n],
                self._boxExtCoordinates[idx][3][m][n]]

    def _getRimExtCoordinates(self, n1, n3):
        """
        :param n1: Index around rim.
        :param n3: Index from inner to outer rim.
        :return: x[], d1[], d2[], and d3[] of rim nodes extended from the tube segment.
        """
        idx = 0 if self._isStartCap else -1
        transitionNodeCount = (len(self._transitionExtCoordinates[idx][0])
                               if (self._transitionExtCoordinates and self._transitionExtCoordinates[idx]) else 0)

        if n3 < transitionNodeCount:
            return [self._transitionExtCoordinates[idx][0][n3][n1],
                    self._transitionExtCoordinates[idx][1][n3][n1],
                    self._transitionExtCoordinates[idx][2][n3][n1],
                    self._transitionExtCoordinates[idx][3][n3][n1]]
        sn3 = n3 - transitionNodeCount
        return [self._shellExtCoordinates[idx][0][sn3][n1],
                self._shellExtCoordinates[idx][1][sn3][n1],
                self._shellExtCoordinates[idx][2][sn3][n1],
                self._shellExtCoordinates[idx][3][sn3][n1]]

    def _getRimExtCoordinatesAround(self, n3):
        """
        :param n3: Index from inner to outer rim.
        :return: x[], d1[], d2[], and d3[] of rim nodes extended from the tube segment.
        """
        idx = 0 if self._isStartCap else -1
        transitionNodeCount = (len(self._transitionExtCoordinates[idx][0])
                               if (self._transitionExtCoordinates and self._transitionExtCoordinates[idx]) else 0)

        if n3 < transitionNodeCount:
            return [self._transitionExtCoordinates[idx][0][n3],
                    self._transitionExtCoordinates[idx][1][n3],
                    self._transitionExtCoordinates[idx][2][n3],
                    self._transitionExtCoordinates[idx][3][n3]]
        sn3 = n3 - transitionNodeCount
        return [self._shellExtCoordinates[idx][0][sn3],
                self._shellExtCoordinates[idx][1][sn3],
                self._shellExtCoordinates[idx][2][sn3],
                self._shellExtCoordinates[idx][3][sn3]]

    def _getRimCoordinatesWithCore(self, m, n, n3):
        """
        Get coordinates and derivatives for cap rim. Only applies when core option is active.
        :param m: Index across major axis.
        :param n: Index across minor axis.
        :param n3: Index along the tube.
        :return: cap rim coordinates and derivatives for points at n3, m and n.
        """
        idx = 0 if self._isStartCap else -1
        return [self._shellCoordinates[idx][0][n3][m][n],
                self._shellCoordinates[idx][1][n3][m][n],
                self._shellCoordinates[idx][2][n3][m][n],
                self._shellCoordinates[idx][3][n3][m][n]]

    def _getTubeBoxCoordinates(self, m, n):
        """
        Get coordinates and derivatives for tube box.
        :param m: Index across major axis.
        :param n: Index across minor axis.
        :return: Tube box coordinates and derivatives for points at n2, m and n.
        """
        idx = 0 if self._isStartCap else -1
        return [self._tubeBoxCoordinates[0][idx][m][n],
                self._tubeBoxCoordinates[1][idx][m][n],
                self._tubeBoxCoordinates[2][idx][m][n],
                self._tubeBoxCoordinates[3][idx][m][n]]

    def _getTubeRimCoordinates(self, n1, n2, n3):
        """
        Get rim parameters at a point.
        :param n1: Node index around.
        :param n2: Node index along segment.
        :param n3: Node index from inner to outer rim.
        :return: Tube rim coordinates and derivatives for points at n1, n2 and n3.
        """
        transitionNodeCount = (len(self._tubeTransitionCoordinates[0][0])
                               if (self._tubeTransitionCoordinates and self._tubeTransitionCoordinates[0]) else 0)
        if n3 < transitionNodeCount and self._isCore:
            return [self._tubeTransitionCoordinates[0][n2][n3][n1],
                    self._tubeTransitionCoordinates[1][n2][n3][n1],
                    self._tubeTransitionCoordinates[2][n2][n3][n1],
                    self._tubeTransitionCoordinates[3][n2][n3][n1]]
        sn3 = n3 - transitionNodeCount
        return [self._tubeShellCoordinates[0][n2][sn3][n1],
                self._tubeShellCoordinates[1][n2][sn3][n1],
                self._tubeShellCoordinates[2][n2][sn3][n1],
                self._tubeShellCoordinates[3][n2][sn3][n1]]

    def _getTriplePointIndexes(self):
        """
        Get a node ID at triple points (special four corners) of the solid core.
        :return: A list of circular (n1) indexes used to identify triple points.
        """
        elementsCountAround = self._elementsCountAround
        nodesCountAcrossMinorHalf = self._getNodesCountCoreBoxMinor() // 2
        triplePointIndexesList = []

        for n in range(0, elementsCountAround, elementsCountAround // 2):
            triplePointIndexesList.append(n + nodesCountAcrossMinorHalf)
            triplePointIndexesList.append((n - nodesCountAcrossMinorHalf) % elementsCountAround)

        return triplePointIndexesList

    def _getTriplePointLocation(self, e1, isShell=False):
        """
        Determines the location of a specific triple point relative to the solid core box.
        There are four locations: Top left (location = 1); top right (location = -1); bottom left (location = 2);
        and bottom right (location = -2). Location is None if not located at any of the four specified locations.
        :param isShell: True if the triple point is located on the shell layer, False if located on the core box.
        :return: Location identifier.
        """
        em = self._elementsCountCoreBoxMinor // 2
        eM = self._elementsCountCoreBoxMajor // 2
        ec = self._elementsCountAround // 4

        lftColumnElements = list(range(0, ec - eM)) + list(range(3 * ec + eM, self._elementsCountAround))
        topRowElements = list(range(ec - eM, ec + eM))
        rhtColumnElements = list((range(2 * ec - em, 2 * ec + em)))
        btmRowElements = list(range(3 * ec - eM, 3 * ec + eM))

        ni = len(lftColumnElements) // 2
        if e1 == topRowElements[0] or e1 == lftColumnElements[ni - 1]:
            location = 1  # "TopLeft"
        elif e1 == topRowElements[-1] or e1 == rhtColumnElements[0]:
            location = -1  # "TopRight"
        elif e1 == btmRowElements[-1] or e1 == lftColumnElements[ni]:
            location = 2  # "BottomLeft"
        elif e1 == btmRowElements[0] or e1 == rhtColumnElements[-1]:
            location = -2  # "BottomRight"
        else:
            location = 0

        if isShell and not self._isStartCap and location != 0:
            location = location + 1 if location > 0 else location - 1
            if abs(location) > 2:
                location = location - 2 if location > 0 else location + 2

        return location

    def _getBoxBoundaryLocation(self, m, n):
        """
        Determines: 1. Where the node is located on the box boundary. There are four side locations: Top, bottom, left,
        and right; and 2. the location of a triple point relative to the solid core box. There are four locations:
        Top left, top right, bottom left, and bottom right.
        Location is 0 if not located at any of the four specified locations.
        :return: Side location identifier and triple point location identifier.
        """
        mEnd = self._getNodesCountCoreBoxMajor() - 1
        nEnd = self._getNodesCountCoreBoxMinor() - 1

        m = m + self._getNodesCountCoreBoxMajor() if m < 0 else m
        n = n + self._getNodesCountCoreBoxMinor() if n < 0 else n

        if 0 < m < mEnd:
            location = 1 if n == nEnd else -1 if n == 0 else 0
        elif 0 < n < nEnd:
            location = 2 if m == 0 else -2 if m == mEnd else 0
        else:
            location = 0

        tpLocation = 0
        if location == 0:
            if n == nEnd:
                tpLocation = 1 if m == 0 else -1  # Top Left or Top Right
            elif n == 0:
                tpLocation = 2 if m == 0 else -2  # Bottom Left or Bottom Right

        return location, tpLocation

    def _getNodesCountCoreBoxMajor(self):
        idx = -1 if self._boxCoordinates[0] is None else 0
        return len(self._boxCoordinates[idx][0])

    def _getNodesCountCoreBoxMinor(self):
        idx = -1 if self._boxCoordinates[0] is None else 0
        return len(self._boxCoordinates[idx][0][0])

    def _getNodesCountRim(self):
        nodesCountRim = self._elementsCountThroughShell + self._elementsCountTransition
        return nodesCountRim

    def _getElementsCountRim(self):
        elementsCountRim = max(1, self._elementsCountThroughShell)
        if self._isCore:
            elementsCountRim += self._elementsCountTransition
        return elementsCountRim

    def _sampleCurvesOnSphere(self, x1, x2, origin, elementsOut):
        """
        Sample coordinates and d1 derivatives of
        :param x1, x2: Coordinates of points 1 and 2 on the spherical surface of cap mesh.
        :param origin: Centre point coordinates.
        :param elementsOut: The number of elements required between points 1 and 2.
        :return: Lists of sampled x and d1 between points 1 and 2.
        """
        r1, r2 = sub(x1, origin), sub(x2, origin)
        deltax = sub(r2, r1)
        normal = cross(r1, deltax)
        theta = angle(r1, r2)
        anglePerElement = theta / elementsOut
        arcLengthPerElement = calculate_arc_length(x1, x2, origin) / elementsOut

        nx, nd1 = [], []
        for n1 in range(elementsOut + 1):
            radiansAcross = n1 * anglePerElement
            r = rotate_vector_around_vector(r1, normal, radiansAcross)
            nx.append(add(r, origin))
            nd1.append(set_magnitude(cross(normal, r), arcLengthPerElement))

        return nx, nd1

    def _generateNodesWithoutCore(self):
        """
        Blackbox function for generating cap nodes. Used only when the tube segment does not have a core.
        """
        generateData = self._generateData
        coordinates = generateData.getCoordinates()
        fieldcache = generateData.getFieldcache()
        nodes = generateData.getNodes()
        nodetemplate = generateData.getNodetemplate()

        nodesCountShell = len(self._tubeShellCoordinates[0][0])
        capNodeIds = []
        idx = 0 if self._isStartCap else -1
        for n2 in range(2):
            capNodeIds.append([])
            for n3 in range(nodesCountShell):
                capNodeIds[n2].append([])
                if n2 == 0:  # apex
                    rx, rd1, rd2, rd3 = (self._shellCoordinates[idx][i][n3][n2] for i in range(4))
                    nodeIdentifier = generateData.nextNodeIdentifier()
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    fieldcache.setNode(node)
                    for nodeValue, rValue in zip([Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                  Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3],
                                                 [rx, rd1, rd2, rd3]):
                        coordinates.setNodeParameters(fieldcache, -1, nodeValue, 1, rValue)
                    capNodeIds[n2][n3].append(nodeIdentifier)
                else:
                    for n1 in range(self._elementsCountAround):
                        rx, rd1, rd2, rd3 = (self._shellCoordinates[idx][i][n3][n2][n1] for i in range(4))
                        nodeIdentifier = generateData.nextNodeIdentifier()
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        for nodeValue, rValue in zip([Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                      Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3],
                                                     [rx, rd1, rd2, rd3]):
                            coordinates.setNodeParameters(fieldcache, -1, nodeValue, 1, rValue)
                        capNodeIds[n2][n3].append(nodeIdentifier)

        if self._isStartCap:
            self._startCapNodeIds = capNodeIds
        else:
            self._endCapNodeIds = capNodeIds

    def _generateNodesWithCore(self):
        """
        Blackbox function for generating cap nodes. Used only when the tube segment has a core.
        """
        generateData = self._generateData
        coordinates = generateData.getCoordinates()
        fieldcache = generateData.getFieldcache()
        nodes = generateData.getNodes()
        nodetemplate = generateData.getNodetemplate()

        nodesCountCoreBoxMajor = self._getNodesCountCoreBoxMajor()
        nodesCountCoreBoxMinor = self._getNodesCountCoreBoxMinor()
        nloop = self._getNodesCountRim() + 1
        capNodeIds = []
        idx = 0 if self._isStartCap else -1
        for n3 in range(nloop):
            if n3 == 0:
                capNodeIds.append([])
                for m in range(nodesCountCoreBoxMajor):
                    capNodeIds[n3].append([])
                    for n in range(nodesCountCoreBoxMinor):
                        rx, rd1, rd2, rd3 = (self._boxCoordinates[idx][i][m][n] for i in range(4))
                        nodeIdentifier = generateData.nextNodeIdentifier()
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        for nodeValue, rValue in zip([Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                      Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3],
                                                     [rx, rd1, rd2, rd3]):
                            coordinates.setNodeParameters(fieldcache, -1, nodeValue, 1, rValue)
                        capNodeIds[n3][m].append(nodeIdentifier)
            else:
                capNodeIds.append([])
                for m in range(nodesCountCoreBoxMajor):
                    n3p = n3 - 1
                    capNodeIds[n3].append([])
                    for n in range(nodesCountCoreBoxMinor):
                        rx, rd1, rd2, rd3 = (self._shellCoordinates[idx][i][n3p][m][n] for i in range(4))
                        nodeIdentifier = generateData.nextNodeIdentifier()
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        for nodeValue, rValue in zip([Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                      Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3],
                                                     [rx, rd1, rd2, rd3]):
                            coordinates.setNodeParameters(fieldcache, -1, nodeValue, 1, rValue)
                        capNodeIds[n3][m].append(nodeIdentifier)

        if self._isStartCap:
            self._startCapNodeIds = capNodeIds
        else:
            self._endCapNodeIds = capNodeIds

    def _generateExtendedTubeNodes(self):
        """
        Blackbox function for generating tube nodes extended from the original tube segment.
        """
        idx = 0 if self._isStartCap else -1
        generateData = self._generateData
        coordinates = generateData.getCoordinates()
        fieldcache = generateData.getFieldcache()
        nodes = generateData.getNodes()
        nodetemplate = generateData.getNodetemplate()

        # create core box nodes
        self._boxExtNodeIds = [None, None] if self._boxExtNodeIds is None else self._boxExtNodeIds
        if self._isCore:
            self._boxExtNodeIds[idx] = []
            nodesCountCoreBoxMajor = self._getNodesCountCoreBoxMajor()
            nodesCountAcrossMinor = self._getNodesCountCoreBoxMinor()
            for n3 in range(nodesCountCoreBoxMajor):
                self._boxExtNodeIds[idx].append([])
                rx, rd1, rd2, rd3 = [self._boxExtCoordinates[idx][i][n3] for i in range(4)]
                for n1 in range(nodesCountAcrossMinor):
                    nodeIdentifier = generateData.nextNodeIdentifier()
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    fieldcache.setNode(node)
                    for nodeValue, rValue in zip([Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                  Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3],
                                                 [rx[n1], rd1[n1], rd2[n1], rd3[n1]]):
                        coordinates.setNodeParameters(fieldcache, -1, nodeValue, 1, rValue)
                    self._boxExtNodeIds[idx][n3].append(nodeIdentifier)

        # create rim nodes and transition nodes (if there are more than 1 layer of transition)
        nodesCountRim = self._getNodesCountRim()
        elementsCountTransition = self._elementsCountTransition
        self._rimExtNodeIds = [None, None] if self._rimExtNodeIds is None else self._rimExtNodeIds
        self._rimExtNodeIds[idx] = []
        for n3 in range(nodesCountRim):
            n3p = n3 - (elementsCountTransition - 1)
            tx = self._transitionExtCoordinates[idx] if self._isCore and elementsCountTransition > 1 and n3 < (
                    elementsCountTransition - 1) else self._shellExtCoordinates[idx]
            rx, rd1, rd2, rd3 = [tx[i][n3 if self._isCore and elementsCountTransition > 1 and n3 < (
                    elementsCountTransition - 1) else n3p] for i in range(4)]
            ringNodeIds = []
            for n1 in range(self._elementsCountAround):
                nodeIdentifier = generateData.nextNodeIdentifier()
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                fieldcache.setNode(node)
                for nodeValue, rValue in zip([Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                              Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3],
                                             [rx[n1], rd1[n1], rd2[n1], rd3[n1]]):
                    coordinates.setNodeParameters(fieldcache, -1, nodeValue, 1, rValue)
                ringNodeIds.append(nodeIdentifier)
            self._rimExtNodeIds[idx].append(ringNodeIds)

    def _generateElementsWithoutCore(self, elementsCountRim, annotationMeshGroups):
        """
        Blackbox function for generating cap elements. Used only when the tube segment does not have a core.
        :param elementsCountRim: Number of elements through the rim.
        :param annotationMeshGroups: List of all annotated mesh groups.
        """
        generateData = self._generateData
        coordinates = generateData.getCoordinates()
        mesh = generateData.getMesh()
        elementtemplateStd, eftStd = generateData.getStandardElementtemplate()
        eftfactory = eftfactory_tricubichermite(mesh, False)
        isStartCap = self._isStartCap

        if isStartCap:
            capNodeIds = self._startCapNodeIds
            self._startCapElementIds = [] if self._startCapElementIds is None else self._startCapElementIds
        else:
            capNodeIds = self._endCapNodeIds
            self._endCapElementIds = [] if self._endCapElementIds is None else self._endCapElementIds
        capElementIds = []
        for e3 in range(elementsCountRim):
            capElementIds.append([])
            for e2 in range(2):
                capElementIds[e3].append([])
                if e2 == 0:
                    for e1 in range(self._elementsCountAround):
                        e1p = (e1 + 1) % self._elementsCountAround
                        nids = []
                        for n3 in [e3, e3 + 1]:
                            if isStartCap:
                                nids += [capNodeIds[0][n3][0], capNodeIds[1][n3][e1], capNodeIds[1][n3][e1p]]
                            else:
                                nids += [capNodeIds[1][n3][e1], capNodeIds[1][n3][e1p], capNodeIds[0][n3][0]]
                        elementIdentifier = generateData.nextElementIdentifier()
                        va, vb = e1, (e1 + 1) % self._elementsCountAround
                        if isStartCap:
                            eftCap = eftfactory.createEftShellPoleBottom(va * 100, vb * 100)
                        else:
                            eftCap = eftfactory.createEftShellPoleTop(va * 100, vb * 100)
                        elementtemplateCap = mesh.createElementtemplate()
                        elementtemplateCap.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                        elementtemplateCap.defineField(coordinates, -1, eftCap)
                        element = mesh.createElement(elementIdentifier, elementtemplateCap)
                        element.setNodesByIdentifier(eftCap, nids)

                        # set general linear map coefficients
                        radiansPerElementAround = math.pi * 2.0 / self._elementsCountAround
                        radiansAround = e1 * radiansPerElementAround if isStartCap else math.pi + e1 * radiansPerElementAround
                        radiansAroundNext = ((e1 + 1) % radiansPerElementAround) * radiansPerElementAround if isStartCap \
                            else math.pi + ((e1 + 1) % self._elementsCountAround) * radiansPerElementAround
                        scalefactors = [
                            1.0,
                            math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                            math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                            math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                            math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
                        ]
                        if not isStartCap:
                            for s in [0, 1, 4, 7, 10]:
                                scalefactors[s] *= -1
                        element.setScaleFactors(eftCap, scalefactors)
                        for annotationMeshGroup in annotationMeshGroups:
                            annotationMeshGroup.addElement(element)
                        capElementIds[e3][e2].append(elementIdentifier)
                else:
                    idx = 0 if isStartCap else -1
                    for e1 in range(self._elementsCountAround):
                        e1p = (e1 + 1) % self._elementsCountAround
                        nids = []
                        for n3 in [e3, e3 + 1]:
                            nids += [capNodeIds[e2][n3][e1], capNodeIds[e2][n3][e1p],
                                     self._rimExtNodeIds[idx][n3][e1], self._rimExtNodeIds[idx][n3][e1p]]
                            if not isStartCap:
                                for a in [nids]:
                                    a[-4], a[-2] = a[-2], a[-4]
                                    a[-3], a[-1] = a[-1], a[-3]
                        elementIdentifier = generateData.nextElementIdentifier()
                        element = mesh.createElement(elementIdentifier, elementtemplateStd)
                        element.setNodesByIdentifier(eftStd, nids)
                        capElementIds[e3][e2].append(elementIdentifier)

        if isStartCap:
            self._startCapElementIds = capElementIds
        else:
            self._endCapElementIds = capElementIds

    def _generateElementsWithCore(self, annotationMeshGroups):
        """
        Blackbox function for generating cap elements. Used only when the tube segment has a core.
        """
        isStartCap = self._isStartCap
        idx = 0 if isStartCap else -1
        elementsCountAround = self._elementsCountAround
        elementsCountCoreBoxMinor = self._elementsCountCoreBoxMinor
        elementsCountCoreBoxMajor = self._elementsCountCoreBoxMajor
        elementsCountRim = self._getElementsCountRim()

        generateData = self._generateData
        coordinates = generateData.getCoordinates()
        mesh = generateData.getMesh()
        elementtemplateStd, eftStd = generateData.getStandardElementtemplate()

        nodeLayoutTransition = generateData.getNodeLayoutTransition()
        nodeLayoutCapTransition = generateData.getNodeLayoutCapTransition()

        if isStartCap:
            capNodeIds = self._startCapNodeIds
            self._startCapElementIds = [] if self._startCapElementIds is None else self._startCapElementIds
        else:
            capNodeIds = self._endCapNodeIds
            self._endCapElementIds = [] if self._endCapElementIds is None else self._endCapElementIds
        boxExtNodeIds = self._boxExtNodeIds[idx]
        rimExtNodeIds = self._rimExtNodeIds[idx]
        triplePointIndexesList = self._getTriplePointIndexes()

        capElementIds = []
        boxBoundaryNodeIds, boxBoundaryNodeToCapIndex = [], []
        rimBoundaryNodeIds, rimBoundaryNodeToCapIndex = [], []
        for n2 in range(elementsCountRim + 1):
            if n2 == 0:
                boxBoundaryNodeIds.append(self._createBoundaryNodeIdsList(capNodeIds[n2])[0])
                boxBoundaryNodeToCapIndex.append(self._createBoundaryNodeIdsList(capNodeIds[n2])[1])
            else:
                rimBoundaryNodeIds.append(self._createBoundaryNodeIdsList(capNodeIds[n2])[0])
                rimBoundaryNodeToCapIndex.append(self._createBoundaryNodeIdsList(capNodeIds[n2])[1])

        # box & shield
        # box
        boxElementIds = []
        for e3 in range(elementsCountCoreBoxMajor):
            boxElementIds.append([])
            e3p = e3 + 1
            for e1 in range(elementsCountCoreBoxMinor):
                nids, nodeParameters, nodeLayouts = [], [], []
                for n1 in [e1, e1 + 1]:
                    nids += [capNodeIds[0][e3][n1], capNodeIds[0][e3p][n1],
                             boxExtNodeIds[e3][n1], boxExtNodeIds[e3p][n1]]
                    if not isStartCap:
                        for a in [nids]:
                            a[-4], a[-2] = a[-2], a[-4]
                            a[-3], a[-1] = a[-1], a[-3]
                elementIdentifier = generateData.nextElementIdentifier()
                element = mesh.createElement(elementIdentifier, elementtemplateStd)
                element.setNodesByIdentifier(eftStd, nids)
                for annotationMeshGroup in annotationMeshGroups:
                    annotationMeshGroup.addElement(element)
                boxElementIds[e3].append(elementIdentifier)
        capElementIds.append(boxElementIds)

        # box shield elements (elements joining the box and the shell elements)
        boxshieldElementIds = []
        for e3 in range(elementsCountCoreBoxMajor):
            boxshieldElementIds.append([])
            e3p = e3 + 1
            for e1 in range(elementsCountCoreBoxMinor):
                nids, nodeParameters, nodeLayouts = [], [], []
                elementIdentifier = generateData.nextElementIdentifier()
                for n1 in [e1, e1 + 1]:
                    for n3 in [e3, e3p]:
                        nids += [capNodeIds[1][n3][n1]]
                        nodeParameter = self._getRimCoordinatesWithCore(n3, n1, 0)
                        nodeParameters.append(nodeParameter)
                        nodeLayouts.append(nodeLayoutCapTransition)
                    for n3 in [e3, e3p]:
                        boxLocation, tpLocation = self._getBoxBoundaryLocation(n3, n1)
                        nid = capNodeIds[0][n3][n1]
                        nids += [nid]
                        nodeParameter = self._getBoxCoordinates(n3, n1)
                        nodeParameters.append(nodeParameter)
                        nodeLayoutCapBoxShield = generateData.getNodeLayoutCapBoxShield(boxLocation, isStartCap)
                        nodeLayoutCapBoxShieldTriplePoint = generateData.getNodeLayoutCapBoxShieldTriplePoint(tpLocation, isStartCap)
                        if nid in boxBoundaryNodeIds[0]:
                            nodeLayouts.append(
                                nodeLayoutCapBoxShield if tpLocation == 0 else nodeLayoutCapBoxShieldTriplePoint)
                        else:
                            nodeLayouts.append(None)
                    if not isStartCap:
                        for a in [nids, nodeParameters, nodeLayouts]:
                            a[-4], a[-2] = a[-2], a[-4]
                            a[-3], a[-1] = a[-1], a[-3]
                eft, scalefactors = determineCubicHermiteSerendipityEft(mesh, nodeParameters, nodeLayouts)
                elementtemplate = mesh.createElementtemplate()
                elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                elementtemplate.defineField(coordinates, -1, eft)
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, nids)
                if scalefactors:
                    element.setScaleFactors(eft, scalefactors)
                for annotationMeshGroup in annotationMeshGroups:
                    annotationMeshGroup.addElement(element)
                boxshieldElementIds[e3].append(elementIdentifier)
        capElementIds.append(boxshieldElementIds)

        # shield
        for e3 in range(elementsCountRim - 1):
            e3p = e3 + 1
            shieldElementIds = []
            for e2 in range(self._elementsCountCoreBoxMajor):
                e2p = e2 + 1
                shieldElementIds.append([])
                for e1 in range(elementsCountCoreBoxMinor):
                    e1p = e1 + 1
                    nids = []
                    for n3 in [e3p, e3p + 1]:
                        nids += [capNodeIds[n3][e2][e1], capNodeIds[n3][e2p][e1],
                                 capNodeIds[n3][e2][e1p], capNodeIds[n3][e2p][e1p]]
                        if not isStartCap:
                            for a in [nids]:
                                a[-4], a[-2] = a[-2], a[-4]
                                a[-3], a[-1] = a[-1], a[-3]
                    elementIdentifier = generateData.nextElementIdentifier()
                    element = mesh.createElement(elementIdentifier, elementtemplateStd)
                    element.setNodesByIdentifier(eftStd, nids)
                    for annotationMeshGroup in annotationMeshGroups:
                        annotationMeshGroup.addElement(element)
                    shieldElementIds[e2].append(elementIdentifier)
            capElementIds.append(shieldElementIds)

        if isStartCap:
            self._startCapElementIds.append(capElementIds)
        else:
            self._endCapElementIds.append(capElementIds)

        # rim
        capElementIds = []
        # box transition
        ringElementIds = []
        boxExtBoundaryNodeIds, boxExtBoundaryNodestoBoxIds = self._createBoundaryNodeIdsList(boxExtNodeIds)
        for e1 in range(elementsCountAround):
            nids, nodeParameters, nodeLayouts = [], [], []
            n1p = (e1 + 1) % self._elementsCountAround
            boxLocation = self._getTriplePointLocation(e1)
            shellLocation = self._getTriplePointLocation(e1, isShell=True)
            nodeLayoutTransitionTriplePoint = generateData.getNodeLayoutTransitionTriplePoint(boxLocation)
            nodeLayoutCapShellTransitionTriplePoint = generateData.getNodeLayoutCapShellTriplePoint(shellLocation)
            for n3 in [0, 1]:
                for n1 in [e1, n1p]:
                    nid = boxBoundaryNodeIds[n3][n1] if n3 == 0 else rimBoundaryNodeIds[n3 - 1][n1]
                    nids += [nid]
                    mi, ni = boxBoundaryNodeToCapIndex[n3][n1] if n3 == 0 else rimBoundaryNodeToCapIndex[n3 - 1][n1]
                    location, tpLocation = self._getBoxBoundaryLocation(mi, ni)
                    nodeLayoutCapBoxShield = generateData.getNodeLayoutCapBoxShield(location, isStartCap)
                    nodeLayoutCapBoxShieldTriplePoint = generateData.getNodeLayoutCapBoxShieldTriplePoint(tpLocation, isStartCap)
                    if n3 == 0:
                        nodeParameter = self._getBoxCoordinates(mi, ni)
                        nodeLayout = nodeLayoutCapBoxShieldTriplePoint if n1 in triplePointIndexesList else (
                            nodeLayoutCapBoxShield)
                    else:
                        nodeParameter = self._getRimCoordinatesWithCore(mi, ni, 0)
                        nodeLayout = nodeLayoutCapShellTransitionTriplePoint if n1 in triplePointIndexesList else (
                            nodeLayoutCapTransition)
                    nodeParameters.append(nodeParameter)
                    nodeLayouts.append(nodeLayout)
                for n1 in [e1, n1p]:
                    if n3 == 0:
                        nid = boxExtBoundaryNodeIds[n1]
                        mi, ni = boxExtBoundaryNodestoBoxIds[n1]
                        nodeParameter = self._getBoxExtCoordinates(mi, ni)
                    else:
                        nid = rimExtNodeIds[0][n1]
                        nodeParameter = self._getRimExtCoordinates(n1, 0)
                    nids += [nid]
                    nodeParameters.append(nodeParameter)
                    nodeLayouts.append(nodeLayoutTransitionTriplePoint if n1 in triplePointIndexesList and n3 == 0
                                       else nodeLayoutTransition)
                if not isStartCap:
                    for a in [nids, nodeParameters, nodeLayouts]:
                        a[-4], a[-2] = a[-2], a[-4]
                        a[-3], a[-1] = a[-1], a[-3]
            eft, scalefactors = determineCubicHermiteSerendipityEft(mesh, nodeParameters, nodeLayouts)
            elementtemplate = mesh.createElementtemplate()
            elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            elementtemplate.defineField(coordinates, -1, eft)
            elementIdentifier = generateData.nextElementIdentifier()
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, nids)
            if scalefactors:
                element.setScaleFactors(eft, scalefactors)
            for annotationMeshGroup in annotationMeshGroups:
                annotationMeshGroup.addElement(element)
            ringElementIds.append(elementIdentifier)
        capElementIds.append(ringElementIds)

        # shell
        triplePointIndexesList = self._getTriplePointIndexes()
        for e3 in range(elementsCountRim - 1):
            rimElementIds = []
            for e1 in range(self._elementsCountAround):
                nids, nodeParameters, nodeLayouts = [], [], []
                e1p = (e1 + 1) % self._elementsCountAround
                location = self._getTriplePointLocation(e1, isShell=True)
                nodeLayoutCapTransition = generateData.getNodeLayoutCapTransition()
                nodeLayoutCapShellTriplePoint = generateData.getNodeLayoutCapShellTriplePoint(location)
                for n3 in [e3, e3 + 1]:
                    for n1 in [e1, e1p]:
                        nids += [rimBoundaryNodeIds[n3][n1]]
                        mi, ni = rimBoundaryNodeToCapIndex[n3][n1]
                        nodeParameter = self._getRimCoordinatesWithCore(mi, ni, n3)
                        nodeParameters.append(nodeParameter)
                        nodeLayouts.append(nodeLayoutCapShellTriplePoint if n1 in triplePointIndexesList else
                                           nodeLayoutCapTransition)
                    for n1 in [e1, e1p]:
                        nids += [rimExtNodeIds[n3][n1]]
                        nodeParameter = self._getRimExtCoordinates(n1, n3)
                        nodeParameters.append(nodeParameter)
                        nodeLayouts.append(None)
                    if not isStartCap:
                        for a in [nids, nodeParameters, nodeLayouts]:
                            a[-4], a[-2] = a[-2], a[-4]
                            a[-3], a[-1] = a[-1], a[-3]
                elementIdentifier = generateData.nextElementIdentifier()
                eft, scalefactors = determineCubicHermiteSerendipityEft(mesh, nodeParameters, nodeLayouts)
                elementtemplate = mesh.createElementtemplate()
                elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                elementtemplate.defineField(coordinates, -1, eft)
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, nids)
                if scalefactors:
                    element.setScaleFactors(eft, scalefactors)
                for annotationMeshGroup in annotationMeshGroups:
                    annotationMeshGroup.addElement(element)
                rimElementIds.append(elementIdentifier)
            capElementIds.append(rimElementIds)

        if isStartCap:
            self._startCapElementIds.append(capElementIds)
        else:
            self._endCapElementIds.append(capElementIds)

    def _generateExtendedTubeElements(self, tubeBoxNodeIds, tubeRimNodeIds, annotationMeshGroups):
        """
        Blackbox function for generating extended tube elements.
        :param tubeBoxNodeIds: List of tube box nodes.
        :param tubeRimNodeIds: List of tube rim nodes.
        """
        isStartCap = self._isStartCap
        idx = 0 if isStartCap else -1
        scalingMode = 3 if isStartCap else 4
        generateData = self._generateData
        coordinates = generateData.getCoordinates()
        mesh = generateData.getMesh()
        elementtemplateStd, eftStd = generateData.getStandardElementtemplate()
        boxExtElementIds, rimExtElementIds = [], []

        boxExtNodeIds = self._boxExtNodeIds[idx]
        rimExtNodeIds = self._rimExtNodeIds[idx]

        if isStartCap:
            self._startExtElementIds = [] if self._startExtElementIds is None else self._startExtElementIds
        else:
            self._endExtElementIds = [] if self._endExtElementIds is None else self._endExtElementIds

        if self._isCore:
            boxExtBoundaryNodeIds, boxExtBoundaryNodesToBoxIds = self._createBoundaryNodeIdsList(boxExtNodeIds)
            tubeBoxBoundaryNodeIds, tubeBoxBoundaryNodesToBoxIds = self._createBoundaryNodeIdsList(tubeBoxNodeIds[idx])
            # create box elements
            boxElementIds = []
            for e3 in range(self._elementsCountCoreBoxMajor):
                boxElementIds.append([])
                e3p = e3 + 1
                for e1 in range(self._elementsCountCoreBoxMinor):
                    nids = []
                    for n1 in [e1, e1 + 1]:
                        nids += [boxExtNodeIds[e3][n1], boxExtNodeIds[e3p][n1],
                                 tubeBoxNodeIds[idx][e3][n1], tubeBoxNodeIds[idx][e3p][n1]]
                        if not isStartCap:
                            for a in [nids]:
                                a[-4], a[-2] = a[-2], a[-4]
                                a[-3], a[-1] = a[-1], a[-3]
                    elementIdentifier = generateData.nextElementIdentifier()
                    element = mesh.createElement(elementIdentifier, elementtemplateStd)
                    element.setNodesByIdentifier(eftStd, nids)
                    for annotationMeshGroup in annotationMeshGroups:
                        annotationMeshGroup.addElement(element)
                    boxElementIds[e3].append(elementIdentifier)
            boxExtElementIds.append(boxElementIds)

            # create core transition elements first layer after box
            triplePointIndexesList = self._getTriplePointIndexes()
            boxElementIds = []
            for e1 in range(self._elementsCountAround):
                nids, nodeParameters, nodeLayouts = [], [], []
                n1p = (e1 + 1) % self._elementsCountAround
                location = self._getTriplePointLocation(e1)
                nodeLayoutTransition = generateData.getNodeLayoutTransition()
                nodeLayoutTransitionTriplePoint = generateData.getNodeLayoutTransitionTriplePoint(location)
                for n2 in [0, 1]:
                    for n1 in [e1, n1p]:
                        if n2 == 0:
                            nid = boxExtBoundaryNodeIds[n1]
                            mi, ni = boxExtBoundaryNodesToBoxIds[n1]
                            nodeParameter = self._getBoxExtCoordinates(mi, ni)
                        else:
                            nid = tubeBoxBoundaryNodeIds[n1]
                            mi, ni = tubeBoxBoundaryNodesToBoxIds[n1]
                            nodeParameter = self._getTubeBoxCoordinates(mi, ni)
                        nids += [nid]
                        nodeParameters.append(nodeParameter)
                        nodeLayouts.append(nodeLayoutTransitionTriplePoint if n1 in triplePointIndexesList else
                                           nodeLayoutTransition)
                if not isStartCap:
                    for a in [nids, nodeParameters, nodeLayouts]:
                        a[-4], a[-2] = a[-2], a[-4]
                        a[-3], a[-1] = a[-1], a[-3]
                for n2 in [0, 1]:
                    for n1 in [e1, n1p]:
                        if n2 == 0:
                            nid = rimExtNodeIds[0][n1]
                            nodeParameter = self._getRimExtCoordinates(n1, 0)
                        else:
                            nid = tubeRimNodeIds[idx][0][n1]
                            nodeParameter = self._getTubeRimCoordinates(n1, idx, 0)
                        nids += [nid]
                        nodeParameters.append(nodeParameter)
                        nodeLayouts.append(None)
                if not isStartCap:
                    for a in [nids, nodeParameters, nodeLayouts]:
                        a[-4], a[-2] = a[-2], a[-4]
                        a[-3], a[-1] = a[-1], a[-3]
                eft, scalefactors = determineCubicHermiteSerendipityEft(mesh, nodeParameters, nodeLayouts)
                if self._elementsCountTransition == 1:
                    eft, scalefactors = generateData.resolveEftCoreBoundaryScaling(
                        eft, scalefactors, nodeParameters, nids, scalingMode)
                elementtemplate = mesh.createElementtemplate()
                elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                elementtemplate.defineField(coordinates, -1, eft)
                elementIdentifier = generateData.nextElementIdentifier()
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, nids)
                if scalefactors:
                    element.setScaleFactors(eft, scalefactors)
                for annotationMeshGroup in annotationMeshGroups:
                    annotationMeshGroup.addElement(element)
                boxElementIds.append(elementIdentifier)
            boxExtElementIds.append(boxElementIds)

        # create regular rim elements - all elements outside first transition layer
        elementsCountRim = self._getElementsCountRim()
        elementsCountRimRegular = elementsCountRim - 1 if self._isCore else elementsCountRim
        nTransition = elementsCountRimRegular - self._elementsCountThroughShell
        for e3 in range(elementsCountRimRegular):
            if e3 < nTransition:
                boxElementIds = []
            else:
                rimExtElementIds.append([])
            lastTransition = self._isCore and (e3 == (self._elementsCountTransition - 2))
            for e1 in range(self._elementsCountAround):
                elementtemplate, eft = elementtemplateStd, eftStd
                n1p = (e1 + 1) % self._elementsCountAround
                nids = []
                for n3 in [e3, e3 + 1]:
                    nids += [rimExtNodeIds[n3][e1], rimExtNodeIds[n3][n1p],
                             tubeRimNodeIds[idx][n3][e1], tubeRimNodeIds[idx][n3][n1p]]
                    if not isStartCap:
                        for a in [nids]:
                            a[-4], a[-2] = a[-2], a[-4]
                            a[-3], a[-1] = a[-1], a[-3]
                elementIdentifier = generateData.nextElementIdentifier()
                scalefactors = []
                if lastTransition:
                    # get node parameters for computing scale factors
                    nodeParameters = []
                    for n3 in (e3, e3 + 1):
                        for n2 in (0, 1):
                            for n1 in (e1, n1p):
                                if n2 == 0:
                                    nodeParameter = self._getRimExtCoordinates(n1, n3)
                                else:
                                    nodeParameter = self._getTubeRimCoordinates(n1, idx, n3)
                                nodeParameters.append(nodeParameter)
                        if not isStartCap:
                            for a in [nodeParameters]:
                                a[-4], a[-2] = a[-2], a[-4]
                                a[-3], a[-1] = a[-1], a[-3]
                    eft = generateData.createElementfieldtemplate()
                    eft, scalefactors = generateData.resolveEftCoreBoundaryScaling(
                        eft, scalefactors, nodeParameters, nids, scalingMode)
                    elementtemplateTransition = mesh.createElementtemplate()
                    elementtemplateTransition.setElementShapeType(Element.SHAPE_TYPE_CUBE)
                    elementtemplateTransition.defineField(coordinates, -1, eft)
                    elementtemplate = elementtemplateTransition
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft, nids)
                if scalefactors:
                    element.setScaleFactors(eft, scalefactors)
                for annotationMeshGroup in annotationMeshGroups:
                    annotationMeshGroup.addElement(element)
                if e3 < nTransition:
                    boxElementIds.append(elementIdentifier)
                else:
                    rimExtElementIds[-1].append(elementIdentifier)
            if e3 < nTransition:
                boxExtElementIds.append(boxElementIds)

        if isStartCap:
            self._startExtElementIds.append(boxExtElementIds)
            self._startExtElementIds.append(rimExtElementIds)
        else:
            self._endExtElementIds.append(boxExtElementIds)
            self._endExtElementIds.append(rimExtElementIds)

    def sampleCoordinates(self, tubeBoxCoordinates, tubeTransitionCoordinates, tubeShellCoordinates):
        """
        Sample cap coordinates.
        """
        self._tubeBoxCoordinates = tubeBoxCoordinates  # tube box coordinates
        self._tubeTransitionCoordinates = tubeTransitionCoordinates  # tube transition coordinates
        self._tubeShellCoordinates = tubeShellCoordinates  # tube rim coordinates

        if self._isCore:
            self._createShellCoordinatesList()
        for s in range(2):
            self._isStartCap = True if self._isCap[0] and s == 0 else False
            if self._isCap[s]:
                if self._isCore:
                    self._sampleCapCoordinatesWithCore(s)
                else:
                    self._extendTubeEnds()
                    self._sampleCapCoordinatesWithoutCore()
            else:
                continue

        return self._boxExtCoordinates, self._transitionExtCoordinates, self._shellExtCoordinates

    def generateNodes(self, generateData, isStartCap=True, isCore=False):
        """
        Blackbox function for generating cap and extended tube nodes.
        :param generateData: Class object from TubeNetworkMeshGenerateData.
        :param isStartCap: True if generating a cap mesh at the start of a tube segment, False if generating at the end
        of a tube segment.
        :param isCore: True for generating a solid core inside the tube, False for regular tube network.
        """
        self._isStartCap = isStartCap
        self._generateData = generateData
        if isStartCap:
            self._generateNodesWithCore() if isCore else self._generateNodesWithoutCore()
            self._generateExtendedTubeNodes()
        else:
            self._generateExtendedTubeNodes()
            self._generateNodesWithCore() if isCore else self._generateNodesWithoutCore()

    def generateElements(self, elementsCountRim, tubeBoxNodeIds, tubeRimNodeIds, annotationMeshGroups,
                         isStartCap=True, isCore=False):
        """
        Blackbox function for generating cap and extended tube elements.
        :param elementsCountRim: Number of elements through the rim.
        :param tubeBoxNodeIds: List of tube box nodes.
        :param tubeRimNodeIds: List of tube rim nodes.
        :param annotationMeshGroups: List of all annotated mesh groups.
        :param isStartCap: True if generating a cap mesh at the start of a tube segment, False if generating at the end
        of a tube segment.
        :param isCore: True for generating a solid core inside the tube, False for regular tube network.
        """
        self._isStartCap = isStartCap
        if isCore:
            self._generateElementsWithCore(annotationMeshGroups)
            self._generateExtendedTubeElements(tubeBoxNodeIds, tubeRimNodeIds, annotationMeshGroups)
        else:
            self._generateElementsWithoutCore(elementsCountRim, annotationMeshGroups)
            self._generateExtendedTubeElements(tubeBoxNodeIds, tubeRimNodeIds, annotationMeshGroups)

    def _addElementsFromIdentifiers(self, mesh, meshGroup, identifiers, tRange=None, mRange=None, nRange=None,
                                    e2Range=None):
        """
        Adds elements to the mesh group based on a structured list of element identifiers.
        :param mesh: The master mesh object used to look up elements by identifier.
        :param meshGroup: Zinc MeshGroup to add elements to.
        :param identifiers: A nested list of element identifiers.
        :param tRange: Range of transition indices.
        :param mRange: Range along the major axis.
        :param nRange: Range along the minor axis.
        :param e2Range: Range around the tube, for ring/circular indexing.
        """
        for t in tRange:
            if mRange is not None and nRange is not None:
                for m in mRange:
                    for n in nRange:
                        elementIdentifier = identifiers[t][m][n]
                        element = mesh.findElementByIdentifier(elementIdentifier)
                        meshGroup.addElement(element)
            elif e2Range is not None:
                for e2 in e2Range:
                    elementIdentifier = identifiers[t][e2]
                    element = mesh.findElementByIdentifier(elementIdentifier)
                    meshGroup.addElement(element)

    def addBoxElementsToMeshGroup(self, meshGroup, e1Range=None, e2Range=None, e3Range=None):
        """
        Add ranges of box elements to mesh group.
        :param meshGroup: Zinc MeshGroup to add elements to.
        :param e1Range: Range between start and limit element indexes in major / d2 direction.
        :param e2Range: Range between start and limit element indexes around the tube.
        :param e3Range: Range between start and limit element indexes in minor / d3 direction.
        """
        elementsCountAround = self._elementsCountAround
        elementsCountCoreBoxMajor = self._elementsCountCoreBoxMajor
        elementsCountCoreBoxMinor = self._elementsCountCoreBoxMinor
        elementsCountTransition = self._elementsCountTransition

        mesh = meshGroup.getMasterMesh()
        e1Range = range(elementsCountCoreBoxMajor) if e1Range is None else e1Range
        e2Range = range(elementsCountAround) if e2Range is None else e2Range
        e3Range = range(elementsCountCoreBoxMinor) if e3Range is None else e3Range

        for i, isCap in enumerate(self._isCap):
            if not isCap:
                continue

            isStart = (i == 0)
            capElementIds = self._startCapElementIds if isStart else self._endCapElementIds
            extElementIds = self._startExtElementIds if isStart else self._endExtElementIds

            # Cap base block (structured m, t, n)
            self._addElementsFromIdentifiers(mesh, meshGroup, capElementIds[0],
                                             tRange=range(elementsCountTransition + 1), mRange=e1Range, nRange=e3Range)

            # Cap tube elements (t, e2)
            self._addElementsFromIdentifiers(mesh, meshGroup, capElementIds[1],
                                             tRange=range(elementsCountTransition), e2Range=e2Range)

            # Extension base block (just one layer at t=0)
            self._addElementsFromIdentifiers(mesh, meshGroup, extElementIds[0][0:1],
                                             tRange=[0], mRange=e1Range, nRange=e3Range)

            # Extension tube (t > 0, e2)
            self._addElementsFromIdentifiers(mesh, meshGroup, extElementIds[0],
                                             tRange=range(1, elementsCountTransition + 1), e2Range=e2Range)

    def addShellElementsToMeshGroup(self, meshGroup, e1Range=None, e2Range=None, e3Range=None):
        """
        Add ranges of shell elements to mesh group.
        :param meshGroup: Zinc MeshGroup to add elements to.
        :param e1Range: Range between start and limit element indexes in major / d2 direction.
        :param e2Range: Range between start and limit element indexes around the tube.
        :param e3Range: Range between start and limit element indexes in minor / d3 direction.
        """
        elementsCountAround = self._elementsCountAround
        elementsCountCoreBoxMajor = self._elementsCountCoreBoxMajor
        elementsCountCoreBoxMinor = self._elementsCountCoreBoxMinor
        elementsCountTransition = self._elementsCountTransition
        elementsCountThroughShell = self._elementsCountThroughShell

        mesh = meshGroup.getMasterMesh()
        e1Range = range(elementsCountCoreBoxMajor) if e1Range is None else e1Range
        e2Range = range(elementsCountAround) if e2Range is None else e2Range
        e3Range = range(elementsCountCoreBoxMinor) if e3Range is None else e3Range

        for i, isCap in enumerate(self._isCap):
            if not isCap:
                continue

            isStart = (i == 0)
            capElementIds = self._startCapElementIds if isStart else self._endCapElementIds
            extElementIds = self._startExtElementIds if isStart else self._endExtElementIds

            tCapShellRange = range(elementsCountTransition + 1, elementsCountTransition + 1 + elementsCountThroughShell)
            self._addElementsFromIdentifiers(mesh, meshGroup, capElementIds[0],
                                             tRange=tCapShellRange, mRange=e1Range, nRange=e3Range)

            tTubeShellRange = range(elementsCountTransition, elementsCountTransition + elementsCountThroughShell)
            self._addElementsFromIdentifiers(mesh, meshGroup, capElementIds[1], tRange=tTubeShellRange, e2Range=e2Range)

            tExtShellRange = range(elementsCountThroughShell)
            self._addElementsFromIdentifiers(mesh, meshGroup, extElementIds[1], tRange=tExtShellRange, e2Range=e2Range)

    def addAllElementsToMeshGroup(self, meshGroup):
        """
        Add all elements in the segment to mesh group.
        :param meshGroup: Zinc MeshGroup to add elements to.
        """
        self.addBoxElementsToMeshGroup(meshGroup)
        self.addShellElementsToMeshGroup(meshGroup)

    def addSideD2ElementsToMeshGroup(self, side: bool, meshGroup):
        """
        Add elements to the mesh group on side of +d2 or -d2, often matching left and right.
        Only works with even numbers around and phase starting at +d2.
        :param side: False for +d2 direction, True for -d2 direction.
        :param meshGroup: Zinc MeshGroup to add elements to.
        """
        e2Start = (self._elementsCountAround // 4) if side else -((self._elementsCountAround + 2) // 4)
        e2Limit = e2Start + (self._elementsCountAround // 2)
        e2Limit = e2Limit + 1 if (self._elementsCountAround % 4) == 2 else e2Limit
        e2Range = range(e2Start, e2Limit)
        if self._isCore:
            e1Start = (self._elementsCountCoreBoxMajor // 2) if side else 0
            e1Limit = self._elementsCountCoreBoxMajor if side else ((self._elementsCountCoreBoxMajor + 1) // 2)
            e1Range = range(e1Start, e1Limit)
            self.addBoxElementsToMeshGroup(meshGroup, e1Range=e1Range, e2Range=e2Range)
        else:
            e1Range = None
        self.addShellElementsToMeshGroup(meshGroup, e1Range=e1Range, e2Range=e2Range)

    def addSideD3ElementsToMeshGroup(self, side: bool, meshGroup):
        """
        Add elements to the mesh group on side of +d3 or -d3, often matching anterior/ventral and posterior/dorsal.
        Only works with even numbers around and phase starting at +d2.
        :param side: False for +d3 direction, True for -d3 direction.
        :param meshGroup: Zinc MeshGroup to add elements to.
        """
        e2Start = (self._elementsCountAround // 2) if side else 0
        e2Limit = e2Start + (self._elementsCountAround // 2)
        e2Range = range(e2Start, e2Limit)
        if self._isCore:
            e3Start = 0 if side else (self._elementsCountCoreBoxMinor // 2)
            e3Limit = ((self._elementsCountCoreBoxMinor + 1) // 2) if side else self._elementsCountCoreBoxMinor
            e3Range = range(e3Start, e3Limit)
            self.addBoxElementsToMeshGroup(meshGroup, e2Range=e2Range, e3Range=e3Range)
        else:
            e3Range = None
        self.addShellElementsToMeshGroup(meshGroup, e2Range=e2Range, e3Range=e3Range)
