"""
Utility class for defining network meshes from 1-D connectivity and lateral axes, with continuity control.
"""

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.maths.vectorops import cross, normalize, sub

import sys


class NetworkNode:
    """
    Describes a single node in a network, storing number of versions etc.
    """

    def __init__(self, nodeIdentifier):
        self._nodeIdentifier = nodeIdentifier
        self._inSegments = []
        self._outSegments = []
        self._interiorSegment = None
        self._versionCount = 1
        self._posX = None  # integer x coordinate in default layout
        self._x = [None] * 3  # coordinates in default layout

    def __hash__(self):
        return hash(self._nodeIdentifier)

    def __eq__(self, other):
        if not isinstance(other, type(self)): return NotImplemented
        return self._nodeIdentifier == other._nodeIdentifier

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

    def getVersionCount(self):
        return self._versionCount

    def incrementVersionCount(self):
        self._versionCount += 1

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
    Describes a segment of a network between junctions as a sequence of nodes with start, end derivative versions.
    """

    def __init__(self, networkNodes, startNodeVersion=None, endNodeVersion=None):
        """
        :param networkNodes: List of NetworkNodes from start to end. Must be at least 2.
        :param startNodeVersion: Version number of derivatives at start node >= 1, or None to use current version
        count of start node.
        :param endNodeVersion: Version number of derivatives at end node >= 1, or None to use current version
        count of end node.
        """
        assert isinstance(networkNodes, list) and (len(networkNodes) > 1)
        assert (startNodeVersion is None) or (startNodeVersion >= 1)
        assert (endNodeVersion is None) or (endNodeVersion >= 1)
        assert isinstance(networkNodes, list) and (len(networkNodes) > 1)
        self._networkNodes = networkNodes
        for networkNode in networkNodes[1:-1]:
            networkNode.setInteriorSegment(self)
        self._startNodeVersion = startNodeVersion if startNodeVersion else networkNodes[0].getVersionCount()
        self._endNodeVersion = endNodeVersion if endNodeVersion else networkNodes[-1].getVersionCount()

    def getNetworkNodes(self):
        return self._networkNodes

    def getStartNodeVersion(self):
        return self._startNodeVersion

    def getEndNodeVersion(self):
        return self._endNodeVersion

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
        nextSegment = NetworkSegment(self._networkNodes[index:], 1, self._endNodeVersion)
        self._networkNodes = self._networkNodes[:index + 1]
        self._endNodeVersion = 1
        return nextSegment


class NetworkMesh:
    """
    Defines a 1-D network with lateral axes, and utility functions for fleshing into 1-D and higher dimensional meshes.
    Only supports a subset of possible networks, e.g. limited to directed acyclic graphs, with fleshing limited to
    bifurcating junctions only.
    """

    def __init__(self, structureString):
        """
        Set up network structure from string description.
        :param structureString: Comma-separated sequences of dash-separated node identifiers defining connectivity and
        continuity throughout network. Each listing of a given node identifier adds an addition version of d1 and side
        derivatives for use in that context. A single node identifier without connections divides a segment at that
        point so annotations can be separate.
        Example: "1-2-4-5-6,3-4-7,5-8" makes the following structure:
                  7
                 /
            1-2-4-5-6
               /   \
              3     8
        where C1 continuity is maintained along each sequence of connected nodes, e.g. 1-2-4-5-6 all use version 1 of
        d1 and side derivatives; 3-4-7 uses version 2 at node 4, and 5-8 uses version 2 at node 5.
        """
        self._networkNodes = {}
        self._networkSegments = []
        sequenceStrings = structureString.split(",")
        for sequenceString in sequenceStrings:
            nodeIdentifiers = [int(s) for s in sequenceString.split("-")]
            sequenceNodes = []
            for nodeIdentifier in nodeIdentifiers:
                # if nodeIdentifiers.count(nodeIdentifier) > 1:
                #     print("Ignoring sequence " + sequenceString + " due to repeated node", nodeIdentifier, \
                #           file=sys.stderr)
                #     break
                networkNode = existingNetworkNode = self._networkNodes.get(nodeIdentifier)
                if networkNode:
                    # future: check whether making a cycle
                    interiorSegment = networkNode.getInteriorSegment()
                    if interiorSegment:
                        nextSegment = interiorSegment.split(networkNode)
                        index = self._networkSegments.index(interiorSegment) + 1
                        self._networkSegments.insert(index, nextSegment)
                    if len(nodeIdentifiers) == 1:
                        if not interiorSegment:
                            print("Cannot split segment at non-interior node", nodeIdentifier, file=sys.stderr)
                        break
                    networkNode.incrementVersionCount()
                else:
                    if len(nodeIdentifiers) == 1:
                        print("Missing split node", nodeIdentifier, file=sys.stderr)
                        break
                    networkNode = NetworkNode(nodeIdentifier)
                    self._networkNodes[nodeIdentifier] = networkNode
                sequenceNodes.append(networkNode)
                if (len(sequenceNodes) > 1) and (existingNetworkNode or (nodeIdentifier == nodeIdentifiers[-1])):
                    networkSegment = NetworkSegment(sequenceNodes)
                    self._networkSegments.append(networkSegment)
                    sequenceNodes = sequenceNodes[-1:]

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
        rangeY = float(maxCountY - 1)
        for posX in range(0, maxPosX + 1):
            xx = float(posX)
            nodes = posNodes[posX]
            countY = len(nodes)
            for iy in range(countY):
                x = [xx, rangeY * (-0.5 + (iy + 0.5) / countY), 0.0]
                nodes[iy].setX(x)


    def create1DLayoutMesh(self, region):
        """
        Expects region Fieldmodule ChangeManager to be in effect.
        :param region: Zinc region to create network layout in.
        """
        fieldmodule = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fieldmodule)

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplates = [None]
        fieldcache = fieldmodule.createFieldcache()

        for nodeIdentifier, networkNode in self._networkNodes.items():
            versionCount = networkNode.getVersionCount()
            print("Node", nodeIdentifier, "versions", versionCount)
            if versionCount < len(nodetemplates):
                nodetemplate = nodetemplates[versionCount]
            else:
                for v in range(len(nodetemplates), versionCount + 1):
                    nodetemplate = nodes.createNodetemplate()
                    nodetemplate.defineField(coordinates)
                    for valueLabel in (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                                       Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3):
                        nodetemplate.setValueNumberOfVersions(coordinates, -1, valueLabel, versionCount)
                    nodetemplates.append(nodetemplate)
            node = nodes.createNode(networkNode.getNodeIdentifier(), nodetemplate)
            fieldcache.setNode(node)
            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, networkNode.getX())
            d1 = []
            for version in range(1, versionCount + 1):
                interiorSegment = networkNode.getInteriorSegment()
                if interiorSegment:
                    segmentNodes = interiorSegment.getNetworkNodes()
                    index = segmentNodes.index(networkNode)
                    prevNetworkNode = segmentNodes[index - 1]
                    nextNetworkNode = segmentNodes[index + 1]
                else:
                    prevNetworkNode = None
                    for inNetworkSegment in networkNode.getInSegments():
                        if inNetworkSegment.getEndNodeVersion() == version:
                            prevNetworkNode = inNetworkSegment.getNetworkNodes()[-2]
                            break
                    nextNetworkNode = None
                    for outNetworkSegment in networkNode.getOutSegments():
                        if outNetworkSegment.getStartNodeVersion() == version:
                            nextNetworkNode = outNetworkSegment.getNetworkNodes()[1]
                            break
                if prevNetworkNode and nextNetworkNode:
                    d1 = sub(nextNetworkNode.getX(), prevNetworkNode.getX())
                elif prevNetworkNode:
                    d1 = sub(networkNode.getX(), prevNetworkNode.getX())
                else:
                    d1 = sub(nextNetworkNode.getX(), networkNode.getX())
                d1 = normalize(d1)
                d2 = [-0.25 * d1[1], 0.25 * d1[0], 0.0]
                d3 = cross(d1, d2)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, version, d1)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, version, d2)
                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, version, d3)

        mesh = fieldmodule.findMeshByDimension(1)
        nextElementIdentifier = 1
        elementtemplates = {}  # dict (startVersion, endVersion) -> Elementtemplate
        for networkSegment in self._networkSegments:
            segmentNodes = networkSegment.getNetworkNodes()
            segmentElementCount = len(segmentNodes) - 1
            for e in range(segmentElementCount):
                startVersion = networkSegment.getStartNodeVersion() if (e == 0) else 1
                endVersion = networkSegment.getEndNodeVersion() if (e == (segmentElementCount - 1)) else 1
                print("Element", nextElementIdentifier, "versions", (startVersion, endVersion))
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
                element = mesh.createElement(nextElementIdentifier, elementtemplate)
                element.setNodesByIdentifier(
                    eft, [segmentNodes[e].getNodeIdentifier(), segmentNodes[e + 1].getNodeIdentifier()])
                nextElementIdentifier += 1
