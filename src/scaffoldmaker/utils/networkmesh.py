"""
Utility class for defining network meshes from 1-D connectivity and lateral axes, with continuity control.
"""

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from cmlibs.maths.vectorops import cross, magnitude, normalize, sub
from scaffoldmaker.utils.interpolation import interpolateSampleCubicHermite, sampleCubicHermiteCurvesSmooth,\
    smoothCubicHermiteDerivativesLoop
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.vector import vectorRejection

import math
import sys


class NetworkNode:
    """
    Describes a single node in a network, storing number of versions etc.
    """

    def __init__(self, nodeIdentifier):
        self._nodeIdentifier = nodeIdentifier
        self._inSegments = []  # segments leading in i.e. ending on this node
        self._outSegments = []  # segments heading out i.e. starting on this node
        self._interiorSegment = None
        self._versionsCount = 0
        self._versionsUsed = set()  # set of version numbers actually used by elements
        self._posX = None  # integer x-coordinate in default layout
        self._x = [None] * 3  # coordinates in default layout

    def getNodeIdentifier(self):
        return self._nodeIdentifier

    def getInteriorSegment(self):
        """
        :return: NetworkSegment this node is interior to i.e. not start nor end, or None if at start or end.
        """
        return self._interiorSegment

    def setInteriorSegment(self, interiorSegment):
        """
        :param interiorSegment: NetworkSegment to set node as interior in, or None to clear.
        """
        self._interiorSegment = interiorSegment

    def addInSegment(self, segment):
        self._inSegments.append(segment)

    def getInSegments(self):
        return self._inSegments

    def addOutSegment(self, segment):
        self._outSegments.append(segment)

    def getOutSegments(self):
        return self._outSegments

    def getVersionsCount(self):
        return self._versionsCount

    def defineVersion(self, nodeVersion):
        self._versionsCount = max(self._versionsCount, nodeVersion)
        self._versionsUsed.add(nodeVersion)

    def checkVersions(self) -> bool:
        """
        :return: True if all node version numbers are in use.
        """
        versionsOK = True
        for nodeVersion in range(1, self._versionsCount + 1):
            if nodeVersion not in self._versionsUsed:
                print("Warning: Node", self._nodeIdentifier, "does not use version", nodeVersion, "of", self._versionsCount,
                      file=sys.stderr)
                versionsOK = False
        return versionsOK

    def getPosX(self):
        return self._posX

    def setPosX(self, posX):
        self._posX = posX

    def getX(self):
        return self._x

    def setX(self, x: list):
        self._x = x


class NetworkSegment:
    """
    Describes a segment of a network between junctions as a sequence of nodes with node derivative versions.
    """

    def __init__(self, networkNodes: list, nodeVersions: list):
        """
        :param networkNodes: List of NetworkNodes from start to end. Must be at least 2.
        :param nodeVersions: List of node versions to use for derivatives at network nodes.
        """
        assert isinstance(networkNodes, list) and (len(networkNodes) > 1) and (len(nodeVersions) == len(networkNodes))
        self._networkNodes = networkNodes
        self._nodeVersions = nodeVersions
        self._elementIdentifiers = [None] * (len(networkNodes) - 1)
        for networkNode in networkNodes[1:-1]:
            networkNode.setInteriorSegment(self)

    def getNetworkNodes(self):
        """
        :return: List of NetworkNode in order along segment.
        """
        return self._networkNodes

    def getNodeIdentifiers(self):
        """
        :return: List of node identifiers in order along segment.
        """
        return [networkNode.getNodeIdentifier() for networkNode in self._networkNodes]

    def getNodeVersions(self):
        """
        :return: List of node version number in node order along segment.
        """
        return self._nodeVersions

    def getElementIdentifiers(self):
        """
        :return: List of element identifiers along segment. Note: identifiers are all None until mesh is generated
        """
        return self._elementIdentifiers

    def setElementIdentifier(self, index, elementIdentifier):
        self._elementIdentifiers[index] = elementIdentifier

    def hasLayoutElementsInMeshGroup(self, meshGroup):
        """
        NetworkSegment must have elementIdentifiers first!
        :param meshGroup: 1D Mesh group to check membership of.
        :return: True if any element in segment is in 1D mesh group.
        """
        for elementIdentifier in self._elementIdentifiers:
            element = meshGroup.findElementByIdentifier(elementIdentifier)
            if element.isValid():
                return True
        return False

    def isCyclic(self):
        """
        Determine by advancing from end node through each out segment recursively whether network cycles back to
        this segment.
        """
        return False  # not implemented, assume not cyclic

    def split(self, splitNetworkNode):
        """
        Split segment to finish at splitNetworkNode, returning remainder as a new NetworkSegment.
        :param splitNetworkNode: NetworkNode to split segment at.
        :return:  New segment after split.
        """
        index = self._networkNodes.index(splitNetworkNode, 1, -1)  # throws exception if not an interior node
        splitNetworkNode.setInteriorSegment(None)
        nextSegment = NetworkSegment(self._networkNodes[index:], self._nodeVersions[index:])
        self._networkNodes = self._networkNodes[:index + 1]
        self._nodeVersions = self._nodeVersions[:index + 1]
        return nextSegment


class NetworkMesh:
    """
    Defines a 1-D network with lateral axes, and utility functions for fleshing into higher dimensional meshes.
    """

    def __init__(self, structureString):
        self._networkNodes = {}
        self._networkSegments = []
        self.build(structureString)
        self._region = None  # set when generated

    def build(self, structureString):
        """
        Set up network structure from structure string description.
        :param structureString: Comma-separated sequences of dash-separated node identifiers defining connectivity and
        continuity throughout network.
        Each listing of a given node identifier is either followed by a dot . and the version number, or version 1 is
        assumed. Each version gives a full set of all derivatives.
        Each sequence is a single contiguous segment until it is split by another sequence referencing an interior
        node.
        Example: "1-2-4-5-6,3-4.2-7,5-8" makes the following structure:
                  7
                 /
            1-2-4-5-6
               /   \
              3     8
        where C1 continuity is maintained along each sequence of connected nodes, e.g. 1-2-4-5-6 all use version 1 of
        d1 and side derivatives; 3-4.2-7 uses version 2 at node 4, and 5-8 uses version 1 at node 5 so both output
        segments at node 5 leave with the same longitudinal derivative and use the same side derivatives.
        """
        self._networkNodes = {}
        self._networkSegments = []
        sequenceStrings = structureString.split(",")
        for sequenceString in sequenceStrings:
            nodeIdentifiers = []
            nodeVersions = []
            nodeVersionStrings = sequenceString.split("-")
            if len(nodeVersionStrings) < 2:
                print("Network mesh: Ignoring empty or single node sequence", sequenceString, file=sys.stderr)
                continue
            try:
                for nodeVersionString in nodeVersionStrings:
                    values = [int(s) for s in nodeVersionString.split(".", maxsplit=1)]
                    nodeIdentifiers.append(values[0])
                    nodeVersions.append(values[1] if (len(values) == 2) else 1)
            except ValueError:
                print("Network mesh: Skipping invalid sequence", sequenceString, file=sys.stderr)
                continue
            sequenceNodes = []
            sequenceVersions = []
            for n in range(len(nodeIdentifiers)):
                nodeIdentifier = nodeIdentifiers[n]
                nodeVersion = nodeVersions[n]
                networkNode = existingNetworkNode = self._networkNodes.get(nodeIdentifier)
                if networkNode:
                    interiorSegment = networkNode.getInteriorSegment()
                    if interiorSegment:
                        nextSegment = interiorSegment.split(networkNode)
                        index = self._networkSegments.index(interiorSegment) + 1
                        self._networkSegments.insert(index, nextSegment)
                else:
                    networkNode = NetworkNode(nodeIdentifier)
                    self._networkNodes[nodeIdentifier] = networkNode
                networkNode.defineVersion(nodeVersion)
                sequenceNodes.append(networkNode)
                sequenceVersions.append(nodeVersion)
                if (len(sequenceNodes) > 1) and (existingNetworkNode or (nodeIdentifier == nodeIdentifiers[-1])):
                    networkSegment = NetworkSegment(sequenceNodes, sequenceVersions)
                    self._networkSegments.append(networkSegment)
                    sequenceNodes = sequenceNodes[-1:]
                    sequenceVersions = sequenceVersions[-1:]

        # warn about nodes without all versions in use
        for networkNode in self._networkNodes.values():
            networkNode.checkVersions()

        # set up in/out connectivity
        for networkSegment in self._networkSegments:
            segmentNodes = networkSegment.getNetworkNodes()
            segmentNodes[0].addOutSegment(networkSegment)
            segmentNodes[-1].addInSegment(networkSegment)

        # assign integer posX coordinates
        for _ in range(len(self._networkSegments)):  # limit total iterations so no endless loop
            changeCount = 0
            for networkSegment in self._networkSegments:
                segmentNodes = networkSegment.getNetworkNodes()
                posX = segmentNodes[0].getPosX()
                if posX is None:
                    posX = 0
                for node in segmentNodes:
                    existingPosX = node.getPosX()
                    if (existingPosX is None) or ((existingPosX < posX) and (
                            (node != segmentNodes[-1]) or not networkSegment.isCyclic())):
                        node.setPosX(posX)
                        changeCount += 1
                    posX += 1
            if changeCount == 0:
                break

        maxPosX = -1
        for node in self._networkNodes.values():
            maxPosX = max(maxPosX, node.getPosX())

        posNodes = []
        for posX in range(0, maxPosX + 1):
            posNodes.append([])
        for node in self._networkNodes.values():
            posNodes[node.getPosX()].append(node)
        maxCountY = 0
        for posX in range(0, maxPosX + 1):
            maxCountY = max(maxCountY, len(posNodes[posX]))
        for posX in range(0, maxPosX + 1):
            xx = float(posX)
            nodes = posNodes[posX]
            countY = len(nodes)
            rangeY = float(countY - 1)
            for iy in range(countY):
                x = [xx, rangeY * (-0.5 + iy / rangeY) if (countY > 1) else 0.0, 0.0]
                nodes[iy].setX(x)

    def getNetworkNodes(self):
        """
        :return: dict mapping node identifier to NetworkNode
        """
        return self._networkSegments

    def getNetworkSegments(self):
        """
        :return: List of segments.
        """
        return self._networkSegments

    def getRegion(self):
        """
        :return: Zinc Region containing generated network layout, or None if not generated
        """
        return self._region

    def create1DLayoutMesh(self, region):
        """
        Expects region Fieldmodule ChangeManager to be in effect.
        :param region: Zinc region to create network layout in.
        """
        self._region = region
        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplates = [None]
        fieldcache = fieldmodule.createFieldcache()
        derivativeValueLabels = [
            Node.VALUE_LABEL_D_DS1,
            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
            Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]
        nodetemplate = None
        for nodeIdentifier, networkNode in self._networkNodes.items():
            versionsCount = networkNode.getVersionsCount()
            # print("Node", nodeIdentifier, "versions", versionsCount)
            if versionsCount < len(nodetemplates):
                nodetemplate = nodetemplates[versionsCount]
            else:
                for v in range(len(nodetemplates), versionsCount + 1):
                    nodetemplate = nodes.createNodetemplate()
                    nodetemplate.defineField(coordinates)
                    for valueLabel in derivativeValueLabels:
                        nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, v)
                    nodetemplates.append(nodetemplate)
            node = nodes.createNode(networkNode.getNodeIdentifier(), nodetemplate)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, networkNode.getX())
            d1 = []
            for nodeVersion in range(1, versionsCount + 1):
                interiorSegment = networkNode.getInteriorSegment()
                if interiorSegment:
                    segmentNodes = interiorSegment.getNetworkNodes()
                    index = segmentNodes.index(networkNode)
                    prevNetworkNode = segmentNodes[index - 1]
                    nextNetworkNode = segmentNodes[index + 1]
                else:
                    prevNetworkNode = None
                    for inNetworkSegment in networkNode.getInSegments():
                        if inNetworkSegment.getNodeVersions()[-1] == nodeVersion:
                            prevNetworkNode = inNetworkSegment.getNetworkNodes()[-2]
                            break
                    nextNetworkNode = None
                    for outNetworkSegment in networkNode.getOutSegments():
                        if outNetworkSegment.getNodeVersions()[0] == nodeVersion:
                            nextNetworkNode = outNetworkSegment.getNetworkNodes()[1]
                            break
                if prevNetworkNode or nextNetworkNode:
                    if prevNetworkNode and nextNetworkNode:
                        d1 = sub(nextNetworkNode.getX(), prevNetworkNode.getX())
                    elif prevNetworkNode:
                        d1 = sub(networkNode.getX(), prevNetworkNode.getX())
                    else:
                        d1 = sub(nextNetworkNode.getX(), networkNode.getX())
                    d1 = normalize(d1)
                    d3 = [0.0, 0.0, 0.1]
                    d2 = cross(d3, normalize(d1))
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, nodeVersion, d1)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, nodeVersion, d2)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, nodeVersion, d3)
                else:
                    print("Warning: No data to define derivative version", nodeVersion, "at node", nodeIdentifier, ".",
                          file=sys.stderr)

        mesh = fieldmodule.findMeshByDimension(1)
        elementIdentifier = 1
        elementtemplates = {}  # dict (startVersion, endVersion) -> Elementtemplate
        for networkSegment in self._networkSegments:
            segmentNodes = networkSegment.getNetworkNodes()
            segmentNodeVersions = networkSegment.getNodeVersions()
            segmentElementCount = len(segmentNodes) - 1
            for e in range(segmentElementCount):
                startVersion = segmentNodeVersions[e]
                endVersion = segmentNodeVersions[e + 1]
                # print("Element", nextElementIdentifier, "versions", (startVersion, endVersion))
                elementtemplate_eft = elementtemplates.get((startVersion, endVersion))
                if elementtemplate_eft:
                    elementtemplate, eft = elementtemplate_eft
                else:
                    elementtemplate = mesh.createElementtemplate()
                    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
                    elementbasis = fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
                    eft = mesh.createElementfieldtemplate(elementbasis)
                    if startVersion != 1:
                        eft.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, startVersion)
                    if endVersion != 1:
                        eft.setTermNodeParameter(4, 1, 2, Node.VALUE_LABEL_D_DS1, endVersion)
                    elementtemplate.defineField(coordinates, -1, eft)
                    elementtemplates[(startVersion, endVersion)] = elementtemplate, eft
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(
                    eft, [segmentNodes[e].getNodeIdentifier(), segmentNodes[e + 1].getNodeIdentifier()])
                networkSegment.setElementIdentifier(e, elementIdentifier)
                elementIdentifier += 1


def getPathRawTubeCoordinates(pathParameters, elementsCountAround, radius=1.0, phaseAngle=0.0):
    """
    Generate coordinates around and along a tube in parametric space around the path parameters,
    at xi2^2 + xi3^2 = radius at the same density as path parameters.
    :param pathParameters: List over nodes of 6 parameters vectors [cx, cd1, cd2, cd12, cd3, cd13] giving
    coordinates cx along path centre, derivatives cd1 along path, cd2 and cd3 giving side vectors,
    and cd12, cd13 giving rate of change of side vectors. Parameters have 3 components.
    Same format as output of zinc_utils get_nodeset_path_ordered_field_parameters().
    :param elementsCountAround: Number of elements & nodes to create around tube. First location is at +d2.
    :param radius: Radius of tube in xi space.
    :param phaseAngle: Starting angle around ellipse, where 0.0 is at d2, pi/2 is at d3.
    :return: px[][], pd1[][], pd2[][], pd12[][] with first index in range(pointsCountAlong),
    second inner index in range(elementsCountAround)
    """
    assert len(pathParameters) == 6
    pointsCountAlong = len(pathParameters[0])
    assert pointsCountAlong > 1
    assert len(pathParameters[0][0]) == 3

    # sample around circle in xi space, later smooth and re-sample to get even spacing in geometric space
    ellipsePointCount = 16
    aroundScale = 2.0 * math.pi / ellipsePointCount
    sxi = []
    sdxi = []
    angleBetweenPoints = 2.0 * math.pi / ellipsePointCount
    for q in range(ellipsePointCount):
        theta = phaseAngle + q * angleBetweenPoints
        xi2 = radius * math.cos(theta)
        xi3 = radius * math.sin(theta)
        sxi.append([xi2, xi3])
        dxi2 = -xi3 * aroundScale
        dxi3 = xi2 * aroundScale
        sdxi.append([dxi2, dxi3])

    px = []
    pd1 = []
    pd2 = []
    pd12 = []
    for p in range(pointsCountAlong):
        cx, cd1, cd2, cd12, cd3, cd13 = [cp[p] for cp in pathParameters]
        tx = []
        td1 = []
        for q in range(ellipsePointCount):
            xi2 = sxi[q][0]
            xi3 = sxi[q][1]
            x = [(cx[c] + xi2 * cd2[c] + xi3 * cd3[c]) for c in range(3)]
            tx.append(x)
            dxi2 = sdxi[q][0]
            dxi3 = sdxi[q][1]
            d1 = [(dxi2 * cd2[c] + dxi3 * cd3[c]) for c in range(3)]
            td1.append(d1)
        # smooth to get reasonable derivative magnitudes
        td1 = smoothCubicHermiteDerivativesLoop(tx, td1, fixAllDirections=True)
        # resample to get evenly spaced points around loop, temporarily adding start point to end
        ex, ed1, pe, pxi, psf = sampleCubicHermiteCurvesSmooth(tx + tx[:1], td1 + td1[:1], elementsCountAround)
        exi, edxi = interpolateSampleCubicHermite(sxi + sxi[:1], sdxi + sdxi[:1], pe, pxi, psf)
        ex.pop()
        ed1.pop()
        exi.pop()
        edxi.pop()

        # check closeness of x(exi[i]) to ex[i]
        # find nearest xi2, xi3 if above is finite error
        # a small but non-negligible error, but results look fine so not worrying
        # dxi = []
        # for i in range(len(ex)):
        #     xi2 = exi[i][0]
        #     xi3 = exi[i][1]
        #     xi23 = xi2 * xi3
        #     x = [(cx[c] + xi2 * cd2[c] + xi3 * cd3[c]) for c in range(3)]
        #     dxi.append(sub(x, ex[i]))
        # print("error", p, "=", [magnitude(v) for v in dxi])

        # calculate d2, d12 at exi
        ed2 = []
        ed12 = []
        for i in range(len(ex)):
            xi2 = exi[i][0]
            xi3 = exi[i][1]
            d2 = [(cd1[c] + xi2 * cd12[c] + xi3 * cd13[c]) for c in range(3)]
            ed2.append(d2)
            dxi2 = edxi[i][0]
            dxi3 = edxi[i][1]
            d12 = [(dxi2 * cd12[c] + dxi3 * cd13[c]) for c in range(3)]
            ed12.append(d12)

        px.append(ex)
        pd1.append(ed1)
        pd2.append(ed2)
        pd12.append(ed12)

    return px, pd1, pd2, pd12


def resampleTubeCoordinates(rawTubeCoordinates, elementsCountAlong,
                            startSurface: TrackSurface=None, endSurface: TrackSurface=None):
    """
    Generate new tube coordinates evenly spaced along raw tube coordinates, optionally
    starting or ending at intersections at start/end track surfaces.
    :param rawTubeCoordinates: (px, pd1, pd2, pd12) returned by getPathRawTubeCoordinates().
    :param elementsCountAlong: Number of elements in resampled coordinates.
    :param startSurface: Optional TrackSurface specifying start of tube at intersection with it.
    :param endSurface: Optional TrackSurface specifying end of tube at intersection with it.
    :return: sx[][], sd1[][], sd2[][], sd12[][] with first index in range(elementsCountAlong + 1),
    second inner index in range(elementsCountAround)
    """
    px, pd1, pd2, pd12 = rawTubeCoordinates
    pointsCountAlong = len(px)
    elementsCountAround = len(px[0])
    # resample at even spacing along
    sx = [[None] * elementsCountAround for _ in range(elementsCountAlong + 1)]
    sd1 = [[None] * elementsCountAround for _ in range(elementsCountAlong + 1)]
    sd2 = [[None] * elementsCountAround for _ in range(elementsCountAlong + 1)]
    sd12 = [[None] * elementsCountAround for _ in range(elementsCountAlong + 1)]
    for q in range(elementsCountAround):
        cx = [px[p][q] for p in range(pointsCountAlong)]
        cd1 = [pd1[p][q] for p in range(pointsCountAlong)]
        cd2 = [pd2[p][q] for p in range(pointsCountAlong)]
        cd12 = [pd12[p][q] for p in range(pointsCountAlong)]
        startCurveLocation = None
        if startSurface:
            startSurfacePosition, startCurveLocation, startIntersects = startSurface.findNearestPositionOnCurve(cx, cd2)
            if not startIntersects:
                startCurveLocation = None
        endCurveLocation = None
        if endSurface:
            endSurfacePosition, endCurveLocation, endIntersects = endSurface.findNearestPositionOnCurve(cx, cd2)
            if not endIntersects:
                endCurveLocation = None
        qx, qd2, pe, pxi, psf = sampleCubicHermiteCurvesSmooth(
            cx, cd2, elementsCountAlong, startLocation=startCurveLocation, endLocation=endCurveLocation)
        qd1, qd12 = interpolateSampleCubicHermite(cd1, cd12, pe, pxi, psf)
        # swizzle
        for p in range(elementsCountAlong + 1):
            sx[p][q] = qx[p]
            sd1[p][q] = qd1[p]
            sd2[p][q] = qd2[p]
            sd12[p][q] = qd12[p]

    # recalculate d1 around intermediate rings, but still in plane
    # normally looks fine, but d1 derivatives are wavy when very distorted
    pStart = 0 if startSurface else 1
    pLimit = elementsCountAlong + 1 if endSurface else elementsCountAlong
    for p in range(pStart, pLimit):
        # first smooth to get d1 with new directions not tangential to surface
        td1 = smoothCubicHermiteDerivativesLoop(sx[p], sd1[p])
        # constraint to be tangential to surface
        td1 = [vectorRejection(td1[q], normalize(cross(sd1[p][q], sd2[p][q]))) for q in range(elementsCountAround)]
        # smooth magnitudes only
        sd1[p] = smoothCubicHermiteDerivativesLoop(sx[p], td1, fixAllDirections=True)

    return sx, sd1, sd2, sd12
