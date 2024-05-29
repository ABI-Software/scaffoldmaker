"""
Specialisation of Network Mesh for building 2-D and 3-D tube mesh networks.
"""
from cmlibs.maths.vectorops import add, cross, dot, magnitude, mult, normalize, sub
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.eft_utils import remapEftLocalNodes, remapEftNodeValueLabelVersion, setEftScaleFactorIds
from scaffoldmaker.utils.interpolation import DerivativeScalingMode, evaluateCoordinatesOnCurve, \
    smoothCubicHermiteDerivativesLoop, smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.networkmesh import NetworkMesh, NetworkMeshBuilder, NetworkMeshGenerateData, \
    NetworkMeshJunction, NetworkMeshSegment, getPathRawTubeCoordinates, pathValueLabels, resampleTubeCoordinates
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters
import math


class TubeNetworkMeshGenerateData(NetworkMeshGenerateData):
    """
    Data for passing to TubeNetworkMesh generateMesh functions.
    """

    def __init__(self, region, meshDimension, islinearThroughWall, isShowTrimSurfaces,
            coordinateFieldName="coordinates",
                 startNodeIdentifier=1, startElementIdentifier=1):
        """
        :param islinearThroughWall: Callers should only set if 3-D with no core.
        :param isShowTrimSurfaces: Tells junction generateMesh to make 2-D trim surfaces.
        """
        super(TubeNetworkMeshGenerateData, self).__init__(
            region, meshDimension, coordinateFieldName, startNodeIdentifier, startElementIdentifier)
        self._isLinearThroughWall = islinearThroughWall
        self._isShowTrimSurfaces = isShowTrimSurfaces

        # get node template for standard and cross nodes
        self._nodetemplate = self._nodes.createNodetemplate()
        self._nodetemplate.defineField(self._coordinates)
        self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if (meshDimension == 3) and not islinearThroughWall:
            self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        # get element template and eft for standard case
        self._standardElementtemplate = self._mesh.createElementtemplate()
        self._standardElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE if (meshDimension == 3)
                                                     else Element.SHAPE_TYPE_SQUARE)
        elementbasis = self._fieldmodule.createElementbasis(
            meshDimension, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY)
        if (meshDimension == 3) and islinearThroughWall:
            elementbasis.setFunctionType(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        self._standardEft = self._mesh.createElementfieldtemplate(elementbasis)
        self._standardElementtemplate.defineField(self._coordinates, -1, self._standardEft)

    def getStandardEft(self):
        return self._standardEft

    def getStandardElementtemplate(self):
        return self._standardElementtemplate

    def getNodetemplate(self):
        return self._nodetemplate

    def islinearThroughWall(self):
        return self._isLinearThroughWall

    def isShowTrimSurfaces(self):
        return self._isShowTrimSurfaces

class TubeNetworkMeshSegment(NetworkMeshSegment):

    def __init__(self, networkSegment, pathParametersList, elementsCountAround, elementsCountThroughWall):
        """
        :param networkSegment: NetworkSegment this is built from.
        :param pathParametersList: [pathParameters] if 2-D or [outerPathParameters, innerPathParameters] if 3-D
        :param elementsCountAround: Number of elements around this segment.
        :param elementsCountThroughWall: Number of elements between inner and outer tube if 3-D, 1 if 2-D.
        """
        super(TubeNetworkMeshSegment, self).__init__(networkSegment, pathParametersList)
        self._elementsCountAround = elementsCountAround
        assert elementsCountThroughWall > 0
        self._elementsCountThroughWall = elementsCountThroughWall
        self._rawTubeCoordinatesList = []
        self._rawTrackSurfaceList = []
        for pathParameters in pathParametersList:
            px, pd1, pd2, pd12 = getPathRawTubeCoordinates(pathParameters, elementsCountAround)
            self._rawTubeCoordinatesList.append((px, pd1, pd2, pd12))
            nx, nd1, nd2, nd12 = [], [], [], []
            for i in range(len(px)):
                nx += px[i]
                nd1 += pd1[i]
                nd2 += pd2[i]
                nd12 += pd12[i]
            self._rawTrackSurfaceList.append(TrackSurface(len(px[0]), len(px) - 1, nx, nd1, nd2, nd12, loop1=True))
        # list[pathsCount][4] of sx, sd1, sd2, sd12; all [nAlong][nAround]:
        self._sampledTubeCoordinates = [None for p in range(self._pathsCount)]
        self._rimCoordinates = None
        self._rimNodeIds = None

    def getRawTrackSurface(self, pathIndex=0):
        return self._rawTrackSurfaceList[pathIndex]

    def sample(self, targetElementLength):
        trimSurfaces = [self._junctions[j].getTrimSurfaces(self) for j in range(2)]
        elementsCountAlong = max(1, math.ceil(self._length / targetElementLength))
        if ((elementsCountAlong == 1) and
                (self._junctions[0].getSegmentsCount() > 2) and
                (self._junctions[1].getSegmentsCount() > 2)):
            elementsCountAlong = 2  # at least 2 segments if bifurcating at both ends
        elif self._isLoop and (elementsCountAlong < 2):
            elementsCountAlong = 2
        for p in range(self._pathsCount):
            self._sampledTubeCoordinates[p] = resampleTubeCoordinates(
                self._rawTubeCoordinatesList[p], elementsCountAlong,
                startSurface=trimSurfaces[0][p], endSurface=trimSurfaces[1][p])

        if self._dimension == 2:
            self._rimCoordinates = self._sampledTubeCoordinates
        else:
            wallFactor = 1.0 / self._elementsCountThroughWall
            ox, od1, od2 = self._sampledTubeCoordinates[0][0:3]
            ix, id1, id2 = self._sampledTubeCoordinates[1][0:3]
            rx, rd1, rd2, rd3 = [], [], [], []
            for n2 in range(elementsCountAlong + 1):
                for r in (rx, rd1, rd2, rd3):
                    r.append([])
                otx, otd1, otd2 = ox[n2], od1[n2], od2[n2]
                itx, itd1, itd2 = ix[n2], id1[n2], id2[n2]
                wd3 = [mult(sub(otx[n1], itx[n1]), wallFactor) for n1 in range(self._elementsCountAround)]
                for n3 in range(self._elementsCountThroughWall + 1):
                    oFactor = n3 / self._elementsCountThroughWall
                    iFactor = 1.0 - oFactor
                    for r in (rx, rd1, rd2, rd3):
                        r[n2].append([])
                    for n1 in range(self._elementsCountAround):
                        if n3 == 0:
                            x, d1, d2 = itx[n1], itd1[n1], itd2[n1]
                        elif n3 == self._elementsCountThroughWall:
                            x, d1, d2 = otx[n1], otd1[n1], otd2[n1]
                        else:
                            x = add(mult(itx[n1], iFactor), mult(otx[n1], oFactor))
                            d1 = add(mult(itd1[n1], iFactor), mult(otd1[n1], oFactor))
                            d2 = add(mult(itd2[n1], iFactor), mult(otd2[n1], oFactor))
                        d3 = wd3[n1]
                        for r, value in zip((rx, rd1, rd2, rd3), (x, d1, d2, d3)):
                            r[n2][n3].append(value)
            self._rimCoordinates = rx, rd1, rd2, rd3
        self._rimNodeIds = [None] * (elementsCountAlong + 1)

    @classmethod
    def blendSampledCoordinates(cls, segment1 , nodeIndexAlong1, segment2, nodeIndexAlong2):
        nodesCountAround = segment1._elementsCountAround
        nodesCountRim = len(segment1._rimCoordinates[0][0])
        if ((nodesCountAround != segment2._elementsCountAround) or
                (nodesCountRim != len(segment2._rimCoordinates[0][0]))):
            return  # can't blend unless these match

        # blend rim coordinates
        for n3 in range(nodesCountRim):
            s1d2 = segment1._rimCoordinates[2][nodeIndexAlong1][n3]
            s2d2 = segment2._rimCoordinates[2][nodeIndexAlong2][n3]
            for n1 in range(nodesCountAround):
                # harmonic mean magnitude
                s1d2Mag = magnitude(s1d2[n1])
                s2d2Mag = magnitude(s2d2[n1])
                d2Mag = 2.0 / ((1.0 / s1d2Mag) + (1.0 / s2d2Mag))
                d2 = mult(s1d2[n1], d2Mag / s1d2Mag)
                s1d2[n1] = d2
                s2d2[n1] = d2

    def generateMesh(self, generateData: TubeNetworkMeshGenerateData):
        elementsCountAlong = len(self._rimCoordinates[0]) - 1
        elementsCountRim = max(1, len(self._rimCoordinates[0][0]) - 1)
        coordinates = generateData.getCoordinates()
        fieldcache = generateData.getFieldcache()
        startSkipCount = 1 if (self._junctions[0].getSegmentsCount() > 2) else 0
        endSkipCount = 1 if (self._junctions[1].getSegmentsCount() > 2) else 0

        # create nodes
        nodes = generateData.getNodes()
        islinearThroughWall = generateData.islinearThroughWall()
        nodetemplate = generateData.getNodetemplate()
        for n2 in range(elementsCountAlong + 1):
            if (n2 < startSkipCount) or (n2 > elementsCountAlong - endSkipCount):
                self._rimNodeIds[n2] = None
                continue
            if self._rimNodeIds[n2]:
                continue
            if self._isLoop and (n2 == elementsCountAlong):
                self._rimNodeIds[n2] = self._rimNodeIds[0]
                continue

            # create rim nodes
            self._rimNodeIds[n2] = []
            for n3 in range(elementsCountRim + 1):
                rx = self._rimCoordinates[0][n2][n3]
                rd1 = self._rimCoordinates[1][n2][n3]
                rd2 = self._rimCoordinates[2][n2][n3]
                rd3 = None if islinearThroughWall else self._rimCoordinates[3][n2][n3]
                ringNodeIds = []
                for n1 in range(self._elementsCountAround):
                    nodeIdentifier = generateData.nextNodeIdentifier()
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    fieldcache.setNode(node)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx[n1])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1[n1])
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2[n1])
                    if rd3:
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3[n1])
                    ringNodeIds.append(nodeIdentifier)
                self._rimNodeIds[n2].append(ringNodeIds)

        # create elements
        annotationMeshGroups = generateData.getAnnotationMeshGroups(self._annotationTerms)
        mesh = generateData.getMesh()
        elementtemplateStd = generateData.getStandardElementtemplate()
        eftStd = generateData.getStandardEft()
        for e2 in range(startSkipCount, elementsCountAlong - endSkipCount):
            e2p = e2 + 1

            # create rim elements
            for e3 in range(elementsCountRim):
                for e1 in range(self._elementsCountAround):
                    e1p = (e1 + 1) % self._elementsCountAround
                    nids = []
                    for n3 in [e3, e3 + 1] if (self._dimension == 3) else [0]:
                        nids += [self._rimNodeIds[e2][n3][e1], self._rimNodeIds[e2][n3][e1p],
                                 self._rimNodeIds[e2p][n3][e1], self._rimNodeIds[e2p][n3][e1p]]
                    elementIdentifier = generateData.nextElementIdentifier()
                    element = mesh.createElement(elementIdentifier, elementtemplateStd)
                    element.setNodesByIdentifier(eftStd, nids)
                    for annotationMeshGroup in annotationMeshGroups:
                        annotationMeshGroup.addElement(element)


class TubeNetworkMeshJunction(NetworkMeshJunction):
    """
    Describes junction between multiple tube segments, some in, some out.
    """

    def __init__(self, inNetworkSegments: list, outNetworkSegments: list, networkMeshSegments):
        """
        :param inNetworkSegments: List of inward segments.
        :param outNetworkSegments: List of outward segments.
        :param networkMeshSegments: dict NetworkSegment -> TubeNetworkMeshSegment.
        """
        super(TubeNetworkMeshJunction, self).__init__(
            inNetworkSegments, outNetworkSegments, networkMeshSegments)
        # following are calculated in determineCrossIndexes()
        # GRC check following:
        self._connectingCoordinateRings = [[]] * self._segmentsCount  # second row of coordinates from end, made into nodes
        self._endCoordinateRings = [[]] * self._segmentsCount  # row of coordinates at end, NOT made into nodes
        self._aroundCounts = [0] * self._segmentsCount
        self._connectionCounts = [0] * self._segmentsCount
        self._aCrossIndexes = None
        self._bCrossIndexes = None
        pathsCount = self._segments[0].getPathsCount()
        self._trimSurfaces = [[None for p in range(pathsCount)] for s in range(self._segmentsCount)]
        self._calculateTrimSurfaces()

    def _calculateTrimSurfaces(self):
        """
        Calculate surfaces for trimming adjacent segments so they can smoothly transition at junction.
        Algorithm gets 6 intersection points of longitudinal lines around each segment with other tubes or the end.
        These are joined to make a partial cone radiating back to the centre of the junction.
        Longitudinal lines start on the edge adjacent to the other segment most normal to it.
        """
        if self._segmentsCount < 3:
            return
        pathsCount = self._segments[0].getPathsCount()
        # get directions at end of segments' paths:
        dirEnd = [[] for s in range(self._segmentsCount)]
        for s in range(self._segmentsCount):
            endIndex = -1 if self._segmentsIn[s] else 0
            for p in range(pathsCount):
                pathParameters = self._segments[s].getPathParameters(p)
                dirEnd[s].append(normalize(pathParameters[1][endIndex]))

        trimPointsCountAround = 6
        for s in range(self._segmentsCount):
            endIndex = -1 if self._segmentsIn[s] else 0
            for p in range(pathsCount):
                # get index of other segment most normal to this segment
                osMax = None
                osDotMax = 0.0
                for os in range(self._segmentsCount):
                    if os == s:
                        continue
                    osDot = math.fabs(dot(dirEnd[s][p], dirEnd[os][p]))
                    if osDot > osDotMax:
                        osMax = os
                        osDotMax = osDot
                pathParameters = self._segments[s].getPathParameters(p)
                # get components of osMax end direction aligned with d2, d3 and compute initial phase
                xEnd = pathParameters[0][endIndex]  # used as centre trim surface radiates out from
                d2End = pathParameters[2][endIndex]
                d3End = pathParameters[4][endIndex]
                endEllipseNormal = normalize(cross(d2End, d3End))
                osDirEnd = dirEnd[osMax][p]
                if self._segmentsIn[osMax]:
                    osDirEnd = [-d for d in osDirEnd]
                dx = dot(osDirEnd, d2End)
                dy = dot(osDirEnd, d3End)
                phaseAngle = math.atan2(dy, dx)
                lx, ld1, ld2, ld12 = getPathRawTubeCoordinates(
                    pathParameters, trimPointsCountAround, radius=1.0, phaseAngle=phaseAngle)
                pointsCountAlong = len(pathParameters[0])

                # get coordinates and directions of intersection points of longitudinal lines and other track surfaces
                rx = []
                rd1 = []
                trim = False
                for n1 in range(trimPointsCountAround):
                    # proportion1 = n1 / trimPointsCountAround
                    cx = [lx[n2][n1] for n2 in range(pointsCountAlong)]
                    cd1 = [ld1[n2][n1] for n2 in range(pointsCountAlong)]
                    cd2 = [ld2[n2][n1] for n2 in range(pointsCountAlong)]
                    cd12 = [ld12[n2][n1] for n2 in range(pointsCountAlong)]
                    x = lx[endIndex][n1]
                    d1 = ld1[endIndex][n1]
                    maxProportionFromEnd = 0.0
                    # find intersection point with other segments which is furthest from end
                    for os in range(self._segmentsCount):
                        if os == s:
                            continue
                        otherSegment = self._segments[os]
                        otherTrackSurface = otherSegment.getRawTrackSurface(p)
                        otherSurfacePosition, curveLocation, isIntersection = \
                            otherTrackSurface.findNearestPositionOnCurve(
                                cx, cd2, loop=False, sampleEnds=False, sampleHalf=2 if self._segmentsIn[s] else 1)
                        if isIntersection:
                            proportion2 = (curveLocation[0] + curveLocation[1]) / (pointsCountAlong - 1)
                            proportionFromEnd = math.fabs(proportion2 - (1.0 if self._segmentsIn[s] else 0.0))
                            if proportionFromEnd > maxProportionFromEnd:
                                trim = True
                                x, d2 = evaluateCoordinatesOnCurve(cx, cd2, curveLocation, loop=False, derivative=True)
                                d1 = evaluateCoordinatesOnCurve(cd1, cd12, curveLocation, loop=False)
                                n = cross(d1, d2)  # normal to this surface
                                ox, od1, od2 = otherTrackSurface.evaluateCoordinates(
                                    otherSurfacePosition, derivatives=True)
                                on = cross(od1, od2)  # normal to other surface
                                d1 = cross(n, on)
                                maxProportionFromEnd = proportionFromEnd
                    # ensure d1 directions go around in same direction as loop
                    if dot(endEllipseNormal, cross(sub(x, xEnd), d1)) < 0.0:
                        d1 = [-d for d in d1]
                    rx.append(x)
                    rd1.append(d1)

                if trim:
                    rd1 = smoothCubicHermiteDerivativesLoop(rx, rd1, fixAllDirections=True,
                                                            magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
                    rd2 = [sub(rx[n1], xEnd) for n1 in range(trimPointsCountAround)]
                    rd12 = smoothCurveSideCrossDerivatives(rx, rd1, [rd2], loop=True)[0]
                    nx = []
                    nd1 = []
                    nd2 = []
                    nd12 = []
                    for factor in (0.5, 1.5):
                        for n1 in range(trimPointsCountAround):
                            d2 = sub(rx[n1], xEnd)
                            x = add(xEnd, mult(d2, factor))
                            d1 = mult(rd1[n1], factor)
                            d12 = mult(rd12[n1], factor)
                            nx.append(x)
                            nd1.append(d1)
                            nd2.append(d2)
                            nd12.append(d12)
                    trimSurface = TrackSurface(trimPointsCountAround, 1, nx, nd1, nd2, nd12, loop1=True)
                    self._trimSurfaces[s][p] = trimSurface

    def getTrimSurfaces(self, segment):
        """
        :param segment: TubeNetworkMeshSegment which must join at junction.
        :return: List of trim surfaces for paths of segment at junction.
        """
        return self._trimSurfaces[self._segments.index(segment)]

    def sample(self, targetElementLength):
        """
        Blend sampled d2 derivatives across 2-segment junctions with the same version.
        Sample junction coordinates between second-from-end segment coordinates.
        :param targetElementLength: Ignored here as always 2 elements across junction.
        """
        if self._segmentsCount == 2:
            TubeNetworkMeshSegment.blendSampledCoordinates(
                self._segments[0], -1 if self._segmentsIn[0] else 0,
                self._segments[1], -1 if self._segmentsIn[1] else 0)

    def generateMesh(self, generateData: TubeNetworkMeshGenerateData):
        if generateData.isShowTrimSurfaces():
            dimension = generateData.getMeshDimension()
            nodeIdentifier, elementIdentifier = generateData.getNodeElementIdentifiers()
            faceIdentifier = elementIdentifier if (dimension == 2) else None
            for s in range(self._segmentsCount):
                for trimSurface in self._trimSurfaces[s]:
                    if trimSurface:
                        nodeIdentifier, faceIdentifier = \
                            trimSurface.generateMesh(generateData.getRegion(), nodeIdentifier, faceIdentifier)
            if dimension == 2:
                elementIdentifier = faceIdentifier
            generateData.setNodeElementIdentifiers(nodeIdentifier, elementIdentifier)

class TubeNetworkMeshBuilder(NetworkMeshBuilder):

    def __init__(self, networkMesh: NetworkMesh, targetElementDensityAlongLongestSegment: float,
                 defaultElementsCountAround: int, elementsCountThroughWall: int,
                 layoutAnnotationGroups: list = [], annotationElementsCountsAround: list = [],
                 annotationElementsCountsAcrossMajor: list = []):
        super(TubeNetworkMeshBuilder, self).__init__(
            networkMesh, targetElementDensityAlongLongestSegment, layoutAnnotationGroups)
        self._defaultElementsCountAround = defaultElementsCountAround
        self._elementsCountThroughWall = elementsCountThroughWall
        self._annotationElementsCountsAround = annotationElementsCountsAround
        self._annotationElementsCountsAcrossMajor = annotationElementsCountsAcrossMajor
        layoutFieldmodule = self._layoutRegion.getFieldmodule()
        self._layoutInnerCoordinates = layoutFieldmodule.findFieldByName("inner coordinates").castFiniteElement()
        if not self._layoutInnerCoordinates.isValid():
            self._layoutInnerCoordinates = None

    def createSegment(self, networkSegment):
        pathParametersList = [get_nodeset_path_ordered_field_parameters(
            self._layoutNodes, self._layoutCoordinates, pathValueLabels,
            networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())]
        if self._layoutInnerCoordinates:
            pathParametersList.append(get_nodeset_path_ordered_field_parameters(
                self._layoutNodes, self._layoutInnerCoordinates, pathValueLabels,
                networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions()))
        elementsCountAround = self._defaultElementsCountAround
        i = 0
        for layoutAnnotationGroup in self._layoutAnnotationGroups:
            if i >= len(self._annotationElementsCountsAround):
                break
            if self._annotationElementsCountsAround[i] > 0:
                if networkSegment.hasLayoutElementsInMeshGroup(layoutAnnotationGroup.getMeshGroup(self._layoutMesh)):
                    elementsCountAround = self._annotationElementsCountsAround[i]
                    break
            i += 1
        return TubeNetworkMeshSegment(networkSegment, pathParametersList, elementsCountAround,
                                      self._elementsCountThroughWall)

    def createJunction(self, inNetworkSegments, outNetworkSegments, networkMeshSegments):
        """
        :param inNetworkSegments: List of input segments.
        :param outNetworkSegments: List of output segments.
        :param networkMeshSegments: dict NetworkSegment -> NetworkMeshSegment-derived object.
        :return: A TubeNetworkMeshJunction.
        """
        return TubeNetworkMeshJunction(inNetworkSegments, outNetworkSegments, networkMeshSegments)
