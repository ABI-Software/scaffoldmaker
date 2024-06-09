"""
Specialisation of Network Mesh for building 2-D and 3-D tube mesh networks.
"""
from cmlibs.maths.vectorops import add, cross, dot, magnitude, mult, normalize, sub
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.eft_utils import determineCubicHermiteSerendipityEft, HermiteNodeLayoutManager
from scaffoldmaker.utils.interpolation import DerivativeScalingMode, evaluateCoordinatesOnCurve, \
    interpolateCubicHermite, interpolateCubicHermiteDerivative, interpolateLagrangeHermiteDerivative, \
    smoothCubicHermiteDerivativesLoop, smoothCurveSideCrossDerivatives
from scaffoldmaker.utils.networkmesh import NetworkMesh, NetworkMeshBuilder, NetworkMeshGenerateData, \
    NetworkMeshJunction, NetworkMeshSegment, getPathRawTubeCoordinates, pathValueLabels, resampleTubeCoordinates
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters
import copy
import math


class TubeNetworkMeshGenerateData(NetworkMeshGenerateData):
    """
    Data for passing to TubeNetworkMesh generateMesh functions.
    """

    def __init__(self, region, meshDimension, isLinearThroughWall, isShowTrimSurfaces,
            coordinateFieldName="coordinates", startNodeIdentifier=1, startElementIdentifier=1):
        """
        :param isLinearThroughWall: Callers should only set if 3-D with no core.
        :param isShowTrimSurfaces: Tells junction generateMesh to make 2-D trim surfaces.
        """
        super(TubeNetworkMeshGenerateData, self).__init__(
            region, meshDimension, coordinateFieldName, startNodeIdentifier, startElementIdentifier)
        self._isLinearThroughWall = isLinearThroughWall
        self._isShowTrimSurfaces = isShowTrimSurfaces

        # get node template for standard and cross nodes
        self._nodetemplate = self._nodes.createNodetemplate()
        self._nodetemplate.defineField(self._coordinates)
        self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if (meshDimension == 3) and not isLinearThroughWall:
            self._nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        # get element template and eft for standard case
        self._standardElementtemplate = self._mesh.createElementtemplate()
        self._standardElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE if (meshDimension == 3)
                                                     else Element.SHAPE_TYPE_SQUARE)
        elementbasis = self._fieldmodule.createElementbasis(
            meshDimension, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY)
        if (meshDimension == 3) and isLinearThroughWall:
            elementbasis.setFunctionType(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        self._standardEft = self._mesh.createElementfieldtemplate(elementbasis)
        self._standardElementtemplate.defineField(self._coordinates, -1, self._standardEft)

        self._nodeLayoutManager = HermiteNodeLayoutManager()

    def getStandardEft(self):
        return self._standardEft

    def getStandardElementtemplate(self):
        return self._standardElementtemplate

    def getNodeLayoutManager(self):
        return self._nodeLayoutManager
    def getNodetemplate(self):
        return self._nodetemplate

    def isLinearThroughWall(self):
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

    def getElementsCountAround(self):
        return self._elementsCountAround

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
            # copy first sampled tube coordinates, but insert single-entry 'n3' index after n2
            self._rimCoordinates = (
                [[ring] for ring in self._sampledTubeCoordinates[0][0]],
                [[ring] for ring in self._sampledTubeCoordinates[0][1]],
                [[ring] for ring in self._sampledTubeCoordinates[0][2]],
                None)
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

    def getSampledElementsCountAlong(self):
        return len(self._sampledTubeCoordinates[0][0])

    def getSampledTubeCoordinatesRing(self, pathIndex, nodeIndexAlong):
        """
        Get a ring of sampled coordinates at the supplied node index.
        :param pathIndex: 0 for outer/primary, 1 or -1 for inner/secondary.
        :param nodeIndexAlong: Node index from 0 to self._elementsCountAlong, or negative to count from end.
        :return: sx[nAround]
        """
        return self._sampledTubeCoordinates[pathIndex][0][nodeIndexAlong]

    def getElementsCountRim(self):
        return max(1, len(self._rimCoordinates[0][0]) - 1)

    def getNodesCountRim(self):
        return len(self._rimCoordinates[0][0])

    def getRimCoordinatesListAlong(self, n1, n2List, n3):
        """
        Get list of parameters for n2 indexes along segment at given n1, n3.
        :param n1: Node index around segment.
        :param n2List: List of node indexes along segment.
        :param n3: Node index from inner to outer rim.
        :return: [x[], d1[], d2[], d3[]]. d3[] may be None
        """
        paramsList = []
        for i in range(4):
            params = []
            for n2 in n2List:
                params.append(self._rimCoordinates[i][n2][n3][n1] if self._rimCoordinates[i] else None)
            paramsList.append(params)
        return paramsList

    def getRimCoordinates(self, n1, n2, n3):
        """
        Get rim parameters at a point.
        :param n1: Node index around.
        :param n2: Node index along segment.
        :param n3: Node index from inner to outer rim.
        :return: x, d1, d2, d3
        """
        return (self._rimCoordinates[0][n2][n3][n1],
                self._rimCoordinates[1][n2][n3][n1],
                self._rimCoordinates[2][n2][n3][n1],
                self._rimCoordinates[3][n2][n3][n1] if self._rimCoordinates[3] else None)

    def getRimNodeId(self, n1, n2, n3):
        """
        Get a rim node ID for a point.
        :param n1: Node index around.
        :param n2: Node index along segment.
        :param n3: Node index from inner to outer rim.
        :return: Node identifier.
        """
        return self._rimNodeIds[n2][n3][n1]

    def getRimNodeIdsSlice(self, n2):
        """
        Get slice of rim node IDs.
        :param n2: Node index along segment, including negative indexes from end.
        :return: Node IDs arrays through wall and around, or None if not set.
        """
        return self._rimNodeIds[n2]

    def generateMesh(self, generateData: TubeNetworkMeshGenerateData, n2Only=None):
        """
        :param n2Only: If set, create nodes only for that single n2 index along. Must be >= 0!
        """
        elementsCountAlong = len(self._rimCoordinates[0]) - 1
        elementsCountRim = self.getElementsCountRim()
        coordinates = generateData.getCoordinates()
        fieldcache = generateData.getFieldcache()
        startSkipCount = 1 if (self._junctions[0].getSegmentsCount() > 2) else 0
        endSkipCount = 1 if (self._junctions[1].getSegmentsCount() > 2) else 0

        # create nodes
        nodes = generateData.getNodes()
        isLinearThroughWall = generateData.isLinearThroughWall()
        nodetemplate = generateData.getNodetemplate()
        for n2 in range(elementsCountAlong + 1) if (n2Only is None) else [n2Only]:
            if (n2 < startSkipCount) or (n2 > elementsCountAlong - endSkipCount):
                self._rimNodeIds[n2] = None
                continue
            if self._rimNodeIds[n2]:
                continue
            # get shared nodes from single adjacent segment, including loop on itself
            # only handles one in, one out
            if n2 == 0:
                if self._junctions[0].getSegmentsCount() == 2:
                    segments = self._junctions[0].getSegments()
                    rimNodeIds = segments[0].getRimNodeIdsSlice(-1)
                    if rimNodeIds:
                        self._rimNodeIds[n2] = rimNodeIds
                        continue
            if n2 == elementsCountAlong:
                if self._junctions[1].getSegmentsCount() == 2:
                    segments = self._junctions[1].getSegments()
                    rimNodeIds = segments[1].getRimNodeIdsSlice(0)
                    if rimNodeIds:
                        self._rimNodeIds[n2] = rimNodeIds
                        continue

            # create rim nodes
            nodesCountRim = len(self._rimCoordinates[0][0])
            self._rimNodeIds[n2] = []
            for n3 in range(nodesCountRim):
                rx = self._rimCoordinates[0][n2][n3]
                rd1 = self._rimCoordinates[1][n2][n3]
                rd2 = self._rimCoordinates[2][n2][n3]
                rd3 = None if isLinearThroughWall else self._rimCoordinates[3][n2][n3]
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

        if n2Only is not None:
            return

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

    def __init__(self, inSegments: list, outSegments: list):
        """
        :param inSegments: List of inward TubeNetworkMeshSegment.
        :param outSegments: List of outward TubeNetworkMeshSegment.
        """
        super(TubeNetworkMeshJunction, self).__init__(inSegments, outSegments)
        pathsCount = self._segments[0].getPathsCount()
        self._trimSurfaces = [[None for p in range(pathsCount)] for s in range(self._segmentsCount)]
        self._calculateTrimSurfaces()
        # rim indexes are issued for interior points connected to 2 or more segment node indexes
        # based on the outer surface, and reused through the wall
        self._rimIndexToSegmentNodeList = []  # list[rim index] giving list[(segment number, node index around)]
        self._segmentNodeToRimIndex = []  # list[segment number][node index around] to rimIndex
        # rim coordinates sampled in the junction are indexed by n3 (through the wall) and 'rim index'
        self._rimCoordinates = None  # if set, (rx[], rd1[], rd2[], rd3[]) each over [n3][rim index]
        self._rimNodeIds = None  # if set, nodeIdentifier[n3][rim index]

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

    def _sampleMidPoint(self, segmentsParameterLists):
        """
        Get mid-point coordinates and derivatives within junction from 2 or more segments' parameters.
        :param segmentsParameterLists: List over segment indexes s of [x, d1, d2, d3], each with 2 last parameters.
        d3 will be None for 2-D of bicubic-linear.
        :return: Mid-point x, d1, d2, d3. Derivative magnitudes will need smoothing.`
        """
        segmentsIn = [dot(sub(params[0][1], params[0][0]), params[2][1]) > 0.0 for params in segmentsParameterLists]
        segmentsCount = len(segmentsIn)
        assert segmentsCount > 1
        d3Defined = None not in segmentsParameterLists[0][3]
        # for each segment get inward parameters halfway between last 2 parameters
        xi = 0.5
        hx = []
        hd1 = []
        hd2 = []
        hd3 = [] if d3Defined else None
        for s in range(segmentsCount):
            params = segmentsParameterLists[s]
            hd2m = [params[2][i] if segmentsIn[s] else [-d for d in params[2][i]] for i in range(2)]
            hx.append(interpolateCubicHermite(params[0][0], hd2m[0], params[0][1], hd2m[1], xi))
            hd1.append(mult(add(params[1][0], params[1][1]), 0.5 if segmentsIn[s] else -0.5))
            hd2.append(interpolateCubicHermiteDerivative(params[0][0], hd2m[0], params[0][1], hd2m[1], xi))
            if d3Defined:
                hd3.append(mult(add(params[3][0], params[3][1]), 0.5))
        # get lists of mid-point parameters for all segment permutations
        mx = []
        md1 = []
        md2 = []
        md3 = [] if d3Defined else None
        xi = 0.5
        for s1 in range(segmentsCount - 1):
            for s2 in range(s1 + 1, segmentsCount):
                hd2s2 = [-d for d in hd2[s2]]
                mx.append(interpolateCubicHermite(hx[s1], hd2[s1], hx[s2], hd2s2, xi))
                md1.append(mult(add(hd1[s1], [-d for d in hd1[s2]]), 0.5 if segmentsIn[s1] else -0.5))
                md2.append(interpolateCubicHermiteDerivative(hx[s1], hd2[s1], hx[s2], hd2s2, xi))
                if d3Defined:
                    md3.append(mult(add(hd3[s1], hd3[s2]), 0.5))
        if len(segmentsIn) == 2:
            if not segmentsIn[0]:
                md2[0] = [-d for d in md2[0]]
            return mx[0], md1[0], md2[0], md3[0] if d3Defined else None
        cx = [sum(x[c] for x in mx) / segmentsCount for c in range(3)]
        cd3 = [sum(d3[c] for d3 in md3) for c in range(3)] if d3Defined else None
        ns12 = [0.0, 0.0, 0.0]
        for m in range(len(md1)):
            cp12 = normalize(cross(md1[m], md2[m]))
            if d3Defined and (dot(cp12, cd3) < 0.0):
                cp12 = [-d for d in cp12]
            for c in range(3):
                ns12[c] += cp12[c]
        ns12 = normalize(ns12)
        # maintain symmetry of bifurcations
        if (segmentsIn == [True, True, False]) or (segmentsIn == [False, False, True]):
            si = (0, 1)
        elif (segmentsIn == [True, False, True]) or (segmentsIn == [False, True, False]):
            si = (2, 0)
        else:
            si = (1, 2)
        params = [segmentsParameterLists[s] for s in si]
        td = [interpolateLagrangeHermiteDerivative(
            cx, params[s][0][0], [-d for d in params[s][2][0]] if segmentsIn[si[s]] else params[s][2][0], 0.0)
            for s in range(2)]
        if segmentsIn.count(True) == 2:
            # reverse so matches inward directions
            td = [[-c for c in d] for d in td]
        cd1, cd2 = (sub(d, mult(ns12, dot(d, ns12))) for d in td)
        if dot(cross(cd1, cd2), ns12) < 0.0:
            cd1, cd2 = cd2, cd1
        return cx, cd1, cd2, cd3

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
            return

        aroundCounts = [segment.getElementsCountAround() for segment in self._segments]
        rimIndexesCount = 0

        if self._segmentsCount == 3:
            # numbers of elements directly connecting pairs of segments
            connectionCounts = [(aroundCounts[s] + aroundCounts[s - 2] - aroundCounts[s - 1]) // 2 for s in range(3)]
            for s in range(3):
                if aroundCounts[s] != (connectionCounts[s - 1] + connectionCounts[s]):
                    print("Can't make tube bifurcation between elements counts around", aroundCounts)
                    return

            self._rimIndexToSegmentNodeList = []  # list[rim index] giving list[(segment index, node index around)]
            self._segmentNodeToRimIndex = [[None] * aroundCounts[s] for s in range(self._segmentsCount)]
            for s1 in range(3):
                s2 = (s1 + 1) % 3
                startNodeIndex1 = (aroundCounts[s1] - connectionCounts[s1]) // 2
                startNodeIndex2 = connectionCounts[s1] // -2
                for n in range(connectionCounts[s1] + 1):
                    nodeIndex1 = startNodeIndex1 + (n if self._segmentsIn[s1] else (connectionCounts[s1] - n))
                    if self._segmentNodeToRimIndex[s1][nodeIndex1] is None:
                        nodeIndex2 = startNodeIndex2 + ((connectionCounts[s1] - n) if self._segmentsIn[s2] else n)
                        if self._segmentNodeToRimIndex[s2][nodeIndex2] is None:
                            rimIndex = rimIndexesCount
                            # keep list in order from lowest s
                            self._rimIndexToSegmentNodeList.append(
                                [[s1, nodeIndex1], [s2, nodeIndex2]] if (s1 < s2) else
                                [[s2, nodeIndex2], [s1, nodeIndex1]])
                            self._segmentNodeToRimIndex[s2][nodeIndex2] = rimIndex
                            rimIndexesCount += 1
                        else:
                            rimIndex = self._segmentNodeToRimIndex[s2][nodeIndex2]
                            # keep list in order from lowest s
                            segmentNodeList = self._rimIndexToSegmentNodeList[rimIndex]
                            index = 0
                            for i in range(len(segmentNodeList)):
                                if s1 < segmentNodeList[i][0]:
                                    break
                                index += 1
                            segmentNodeList.insert(index, [s1, nodeIndex1])
                        self._segmentNodeToRimIndex[s1][nodeIndex1] = rimIndex

        if not rimIndexesCount:
            return

        # get node indexes giving lowest sum of distances between adjoining points on outer sampled tubes
        permutationCount = 1
        for count in aroundCounts:
            permutationCount *= count
        minIndexes = None
        minSum = None
        indexes = [0] * self._segmentsCount
        rings = [self._segments[s].getSampledTubeCoordinatesRing(0, -1 if self._segmentsIn[s] else 0)
                 for s in range(self._segmentsCount)]
        for p in range(permutationCount):
            sum = 0.0
            for rimIndex in range(rimIndexesCount):
                segmentNodeList = self._rimIndexToSegmentNodeList[rimIndex]
                sCount = len(segmentNodeList)
                for i in range(sCount - 1):
                    s1, n1 = segmentNodeList[i]
                    nodeIndex1 = (n1 + indexes[s1]) % aroundCounts[s1]
                    x1 = rings[s1][nodeIndex1]
                    for j in range(s1 + 1, sCount):
                        s2, n2 = segmentNodeList[j]
                        nodeIndex2 = (n2 + indexes[s2]) % aroundCounts[s2]
                        x2 = rings[s2][nodeIndex2]
                        sum += magnitude([x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]])
            if (minSum is None) or (sum < minSum):
                minIndexes = copy.copy(indexes)
                minSum = sum
            # permute through indexes:
            for s in range(self._segmentsCount):
                indexes[s] += 1
                if indexes[s] < aroundCounts[s]:
                    break
                indexes[s] = 0

        # offset node indexes by minIndexes
        for rimIndex in range(rimIndexesCount):
            segmentNodeList = self._rimIndexToSegmentNodeList[rimIndex]
            for segmentNode in segmentNodeList:
                s, n = segmentNode
                nodeIndex = (n + minIndexes[s]) % aroundCounts[s]
                self._segmentNodeToRimIndex[s][nodeIndex] = rimIndex
                segmentNode[1] = nodeIndex

        # sample rim coordinates
        nodesCountRim = self._segments[0].getNodesCountRim()
        rx, rd1, rd2, rd3 = [
            [[None] * rimIndexesCount for _ in range(nodesCountRim)] for i in range(4)]
        self._rimCoordinates = (rx, rd1, rd2, rd3)
        for n3 in range(nodesCountRim):
            for rimIndex in range(rimIndexesCount):
                segmentNodeList = self._rimIndexToSegmentNodeList[rimIndex]
                # segments have been ordered from lowest to highest s index
                segmentsParameterLists = []
                for s, n1 in segmentNodeList:
                    segmentsParameterLists.append(
                        self._segments[s].getRimCoordinatesListAlong(
                            n1, [-2, -1] if self._segmentsIn[s] else [1, 0], n3))
                rx[n3][rimIndex], rd1[n3][rimIndex], rd2[n3][rimIndex], rd3[n3][rimIndex] = \
                    self._sampleMidPoint(segmentsParameterLists)

        # smooth rim coordinates

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

        if self._segmentsCount < 3:
            return

        rimIndexesCount = len(self._rimIndexToSegmentNodeList)
        nodesCountRim = self._segments[0].getNodesCountRim()
        elementsCountRim = max(1, nodesCountRim - 1)
        if self._rimCoordinates:
            self._rimNodeIds = [[None] * rimIndexesCount for _ in range(nodesCountRim)]

        coordinates = generateData.getCoordinates()
        fieldcache = generateData.getFieldcache()
        nodes = generateData.getNodes()
        nodetemplate = generateData.getNodetemplate()
        isLinearThroughWall = generateData.isLinearThroughWall()
        mesh = generateData.getMesh()
        meshDimension = generateData.getMeshDimension()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(
            Element.SHAPE_TYPE_CUBE if (meshDimension == 3) else Element.SHAPE_TYPE_SQUARE)
        d3Defined = (meshDimension == 3) and not isLinearThroughWall

        nodeLayoutManager = generateData.getNodeLayoutManager()

        # nodes and elements are generated in order of segments
        for s in range(self._segmentsCount):
            segment = self._segments[s]
            n2 = (segment.getSampledElementsCountAlong() - 2) if self._segmentsIn[s] else 1
            segment.generateMesh(generateData, n2Only=n2)

            elementsCountAround = segment.getElementsCountAround()

            if self._rimCoordinates:
                # create rim nodes
                for n3 in range(nodesCountRim):
                    rx = self._rimCoordinates[0][n3]
                    rd1 = self._rimCoordinates[1][n3]
                    rd2 = self._rimCoordinates[2][n3]
                    rd3 = self._rimCoordinates[3][n3] if d3Defined else None
                    layerNodeIds = self._rimNodeIds[n3]
                    for n1 in range(elementsCountAround):
                        rimIndex = self._segmentNodeToRimIndex[s][n1]
                        nodeIdentifier = self._rimNodeIds[n3][rimIndex]
                        if nodeIdentifier is not None:
                            continue
                        nodeIdentifier = generateData.nextNodeIdentifier()
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        fieldcache.setNode(node)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, rx[rimIndex])
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, rd1[rimIndex])
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, rd2[rimIndex])
                        if rd3:
                            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, rd3[rimIndex])
                        layerNodeIds[rimIndex] = nodeIdentifier

            if self._rimCoordinates:
                # create rim elements
                annotationMeshGroups = generateData.getAnnotationMeshGroups(segment.getAnnotationTerms())
                eftList = [None] * elementsCountAround
                scalefactorsList = [None] * elementsCountAround
                for e3 in range(elementsCountRim):
                    for e1 in range(elementsCountAround):
                        n1p = (e1 + 1) % elementsCountAround
                        nids = []
                        nodeParameters = []
                        nodeLayouts = []
                        for n3 in [e3, e3 + 1] if (meshDimension == 3) else [e3]:
                            for n1 in [e1, n1p]:
                                nids.append(self._segments[s].getRimNodeId(n1, n2, n3))
                                if e3 == 0:
                                    rimCoordinates = self._segments[s].getRimCoordinates(n1, n2, n3)
                                    nodeParameters.append(rimCoordinates if d3Defined else
                                        (rimCoordinates[0], rimCoordinates[1], rimCoordinates[2], None))
                                    nodeLayouts.append(None)
                            for n1 in [e1, n1p]:
                                rimIndex = self._segmentNodeToRimIndex[s][n1]
                                nids.append(self._rimNodeIds[n3][rimIndex])
                                if e3 == 0:
                                    nodeParameters.append(
                                        (self._rimCoordinates[0][n3][rimIndex],
                                         self._rimCoordinates[1][n3][rimIndex],
                                         self._rimCoordinates[2][n3][rimIndex],
                                         self._rimCoordinates[3][n3][rimIndex] if d3Defined else None))
                                    segmentNodesCount = len(self._rimIndexToSegmentNodeList[rimIndex])
                                    nodeLayouts.append(
                                        nodeLayoutManager.getNodeLayout6Way12(d3Defined) if (segmentNodesCount == 3)
                                        else nodeLayoutManager.getNodeLayoutRegularPermuted(d3Defined))
                            if not self._segmentsIn[s]:
                                for a in [nids, nodeParameters, nodeLayouts] if (e3 == 0) else [nids]:
                                    a[-4], a[-2] = a[-2], a[-4]
                                    a[-3], a[-1] = a[-1], a[-3]
                        # exploit efts being same through the wall
                        eft = eftList[e1]
                        scalefactors = scalefactorsList[e1]
                        if not eft:
                            eft, scalefactors = determineCubicHermiteSerendipityEft(mesh, nodeParameters, nodeLayouts)
                            eftList[e1] = eft
                            scalefactorsList[e1] = scalefactors
                        elementtemplate.defineField(coordinates, -1, eft)
                        elementIdentifier = generateData.nextElementIdentifier()
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        element.setNodesByIdentifier(eft, nids)
                        if scalefactors:
                            element.setScaleFactors(eft, scalefactors)
                        for annotationMeshGroup in annotationMeshGroups:
                            annotationMeshGroup.addElement(element)


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

    def createJunction(self, inSegments, outSegments):
        """
        :param inSegments: List of inward TubeNetworkMeshSegment.
        :param outSegments: List of outward TubeNetworkMeshSegment.
        :return: A TubeNetworkMeshJunction.
        """
        return TubeNetworkMeshJunction(inSegments, outSegments)
