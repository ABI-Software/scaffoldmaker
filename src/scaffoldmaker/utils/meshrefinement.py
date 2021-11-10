'''
Class for refining a mesh from one region to another.
'''
from __future__ import division

import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, findOrCreateFieldNodeGroup, \
    findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.result import RESULT_OK as ZINC_OK
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.utils.octree import Octree


class MeshRefinement:
    '''
    Class for refining a mesh from one region to another.
    '''

    def __init__(self, sourceRegion, targetRegion, sourceAnnotationGroups = []):
        '''
        Assumes targetRegion is empty.
        :param sourceAnnotationGroups: List of AnnotationGroup for source mesh in sourceRegion.
        A copy containing the refined elements is created by the MeshRefinement.
        '''
        self._sourceRegion = sourceRegion
        self._sourceFm = sourceRegion.getFieldmodule()
        self._sourceCache = self._sourceFm.createFieldcache()
        self._sourceCoordinates = findOrCreateFieldCoordinates(self._sourceFm)
        # get range of source coordinates for octree range
        self._sourceFm.beginChange()
        sourceNodes = self._sourceFm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        minimumsField = self._sourceFm.createFieldNodesetMinimum(self._sourceCoordinates, sourceNodes)
        result, minimums = minimumsField.evaluateReal(self._sourceCache, 3)
        assert result == ZINC_OK, 'MeshRefinement failed to get minimum coordinates'
        maximumsField = self._sourceFm.createFieldNodesetMaximum(self._sourceCoordinates, sourceNodes)
        result, maximums = maximumsField.evaluateReal(self._sourceCache, 3)
        assert result == ZINC_OK, 'MeshRefinement failed to get maximum coordinates'
        xrange = [ (maximums[i] - minimums[i]) for i in range(3) ]
        edgeTolerance = 0.5*(max(xrange))
        if edgeTolerance == 0.0:
            edgeTolerance = 1.0
        minimums = [ (minimums[i] - edgeTolerance) for i in range(3) ]
        maximums = [ (maximums[i] + edgeTolerance) for i in range(3) ]
        minimumsField = None
        maximumsField = None
        self._sourceMesh = self._sourceFm.findMeshByDimension(3)
        self._sourceNodes = self._sourceFm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._sourceElementiterator = self._sourceMesh.createElementiterator()
        self._octree = Octree(minimums, maximums)

        self._targetRegion = targetRegion
        self._targetFm = targetRegion.getFieldmodule()
        self._targetFm.beginChange()
        self._targetCache = self._targetFm.createFieldcache()
        self._targetCoordinates = findOrCreateFieldCoordinates(self._targetFm)

        self._targetNodes = self._targetFm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._nodetemplate = self._targetNodes.createNodetemplate()
        self._nodetemplate.defineField(self._targetCoordinates)

        self._targetMesh = self._targetFm.findMeshByDimension(3)
        self._targetBasis = self._targetFm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        self._targetEft = self._targetMesh.createElementfieldtemplate(self._targetBasis)
        self._targetElementtemplate = self._targetMesh.createElementtemplate()
        self._targetElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = self._targetElementtemplate.defineField(self._targetCoordinates, -1, self._targetEft)

        self._nodeIdentifier = 1
        self._elementIdentifier = 1
        # prepare annotation group map
        self._sourceAnnotationGroups = sourceAnnotationGroups
        self._annotationGroups = []
        self._sourceAndTargetMeshGroups = []
        self._sourceAndTargetNodesetGroups = []
        for sourceAnnotationGroup in sourceAnnotationGroups:
            targetAnnotationGroup = AnnotationGroup(self._targetRegion, sourceAnnotationGroup.getTerm())
            self._annotationGroups.append(targetAnnotationGroup)
            # assume have only highest dimension element or node/point annotation groups:
            if sourceAnnotationGroup.hasMeshGroup(self._sourceMesh):
                self._sourceAndTargetMeshGroups.append( ( sourceAnnotationGroup.getMeshGroup(self._sourceMesh), targetAnnotationGroup.getMeshGroup(self._targetMesh)) )
            else:
                self._sourceAndTargetNodesetGroups.append( ( sourceAnnotationGroup.getNodesetGroup(self._sourceNodes), targetAnnotationGroup.getNodesetGroup(self._targetNodes)) )

        # prepare element -> marker point list map
        self.elementMarkerMap = {}
        sourceMarkerGroup = findOrCreateFieldGroup(self._sourceFm, "marker")
        sourceMarkerName = findOrCreateFieldStoredString(self._sourceFm, name="marker_name")
        sourceMarkerLocation = findOrCreateFieldStoredMeshLocation(self._sourceFm, self._sourceMesh, name="marker_location")
        sourceMarkerNodes = findOrCreateFieldNodeGroup(sourceMarkerGroup, sourceNodes).getNodesetGroup()
        nodeIter = sourceMarkerNodes.createNodeiterator()
        node = nodeIter.next()
        while node.isValid():
            self._sourceCache.setNode(node)
            element, xi = sourceMarkerLocation.evaluateMeshLocation(self._sourceCache, self._sourceMesh.getDimension())
            if element.isValid():
                elementIdentifier = element.getIdentifier()
                markerName = sourceMarkerName.evaluateString(self._sourceCache)
                markerList = self.elementMarkerMap.get(elementIdentifier)
                if not markerList:
                    markerList = []
                    self.elementMarkerMap[elementIdentifier] = markerList
                markerList.append( ( markerName, xi, node.getIdentifier() ) )
            node = nodeIter.next()
        if self.elementMarkerMap:
            self._targetMarkerGroup = findOrCreateFieldGroup(self._targetFm, "marker")
            self._targetMarkerName = findOrCreateFieldStoredString(self._targetFm, name="marker_name")
            self._targetMarkerLocation = findOrCreateFieldStoredMeshLocation(self._targetFm, self._targetMesh, name="marker_location")
            self._targetMarkerNodes = findOrCreateFieldNodeGroup(self._targetMarkerGroup, self._targetNodes).getNodesetGroup()
            self._targetMarkerTemplate = self._targetMarkerNodes.createNodetemplate()
            self._targetMarkerTemplate.defineField(self._targetMarkerName)
            self._targetMarkerTemplate.defineField(self._targetMarkerLocation)


    def __del__(self):
        self._sourceFm.endChange()
        self._targetFm.endChange()


    def getAnnotationGroups(self):
        return self._annotationGroups


    def refineElementCubeStandard3d(self, sourceElement, numberInXi1, numberInXi2, numberInXi3,
            addNewNodesToOctree=True, shareNodeIds=None, shareNodeCoordinates=None):
        '''
        Refine cube sourceElement to numberInXi1*numberInXi2*numberInXi3 linear cube
        sub-elements, evenly spaced in xi.
        :param addNewNodesToOctree: If True (default) add newly created nodes to
        octree to be found when refining later elements. Set to False when nodes are at the
        same location and not intended to be shared.
        :param shareNodeIds, shareNodeCoordinates: Arrays of identifiers and coordinates of
        nodes which may be shared in refining this element. If supplied, these are preferentially
        used ahead of points in the octree. Used to control merging with known nodes, e.g.
        those returned by this function for elements which used addNewNodesToOctree=False.
        :return: Node identifiers, node coordinates used in refinement of sourceElement.
        '''
        assert (shareNodeIds and shareNodeCoordinates) or (not shareNodeIds and not shareNodeCoordinates), \
            'refineElementCubeStandard3d.  Must supply both of shareNodeIds and shareNodeCoordinates, or neither'
        shareNodesCount = len(shareNodeIds) if shareNodeIds else 0
        meshGroups = []
        for sourceAndTargetMeshGroup in self._sourceAndTargetMeshGroups:
            if sourceAndTargetMeshGroup[0].containsElement(sourceElement):
                meshGroups.append(sourceAndTargetMeshGroup[1])
        # create nodes
        nids = []
        nx = []
        xi = [ 0.0, 0.0, 0.0 ]
        tol = self._octree._tolerance
        for k in range(numberInXi3 + 1):
            kExterior = (k == 0) or (k == numberInXi3)
            xi[2] = k/numberInXi3
            for j in range(numberInXi2 + 1):
                jExterior = kExterior or (j == 0) or (j == numberInXi2)
                xi[1] = j/numberInXi2
                for i in range(numberInXi1 + 1):
                    iExterior = jExterior or (i == 0) or (i == numberInXi1)
                    xi[0] = i/numberInXi1
                    self._sourceCache.setMeshLocation(sourceElement, xi)
                    result, x = self._sourceCoordinates.evaluateReal(self._sourceCache, 3)
                    # only exterior points are ever common:
                    nodeId = None
                    if iExterior:
                        if shareNodeIds:
                            for n in range(shareNodesCount):
                                if (math.fabs(shareNodeCoordinates[n][0] - x[0]) <= tol) and \
                                   (math.fabs(shareNodeCoordinates[n][1] - x[1]) <= tol) and \
                                   (math.fabs(shareNodeCoordinates[n][2] - x[2]) <= tol):
                                    nodeId = shareNodeIds[n]
                                    break
                        if nodeId is None:
                            nodeId = self._octree.findObjectByCoordinates(x)
                    if nodeId is None:
                        node = self._targetNodes.createNode(self._nodeIdentifier, self._nodetemplate)
                        self._targetCache.setNode(node)
                        result = self._targetCoordinates.setNodeParameters(self._targetCache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        nodeId = self._nodeIdentifier
                        if iExterior and addNewNodesToOctree:
                            self._octree.addObjectAtCoordinates(x, nodeId)
                        self._nodeIdentifier += 1
                    nids.append(nodeId)
                    nx.append(x)
        # create elements
        startElementIdentifier = self._elementIdentifier
        for k in range(numberInXi3):
            ok = (numberInXi2 + 1)*(numberInXi1 + 1)
            for j in range(numberInXi2):
                oj = (numberInXi1 + 1)
                for i in range(numberInXi1):
                    bni = k*ok + j*oj + i
                    element = self._targetMesh.createElement(self._elementIdentifier, self._targetElementtemplate)
                    enids = [ nids[bni     ], nids[bni      + 1], nids[bni      + oj], nids[bni      + oj + 1],
                              nids[bni + ok], nids[bni + ok + 1], nids[bni + ok + oj], nids[bni + ok + oj + 1] ]
                    result = element.setNodesByIdentifier(self._targetEft, enids)
                    #if result != ZINC_OK:
                    #print('Element', self._elementIdentifier, result, enids)
                    self._elementIdentifier += 1

                    for meshGroup in meshGroups:
                        meshGroup.addElement(element)

        # re-map any markers embedded in the source element
        if self.elementMarkerMap:
            markerList = self.elementMarkerMap.get(sourceElement.getIdentifier())
            if markerList:
                numberInXi = [ numberInXi1, numberInXi2, numberInXi3 ]
                elementOffset = [ 1, numberInXi1, numberInXi1*numberInXi2 ]
                targetXi = [ 0.0 ]*3
                for marker in markerList:
                    markerName, sourceXi, sourceNodeIdentifier = marker
                    sourceNode = self._sourceNodes.findNodeByIdentifier(sourceNodeIdentifier)
                    node = self._targetMarkerNodes.createNode(self._nodeIdentifier, self._targetMarkerTemplate)
                    self._targetCache.setNode(node)
                    self._targetMarkerName.assignString(self._targetCache, markerName)
                    # determine which sub-element, targetXi that sourceXi maps to
                    targetElementIdentifier = startElementIdentifier
                    for i in range(3):
                        targetXi[i] = sourceXi[i]*numberInXi[i]
                        el = int(targetXi[i])
                        if el < numberInXi[i]:
                            targetXi[i] -= el
                        else:
                            el = numberInXi[i] - 1
                            targetXi[i] = 1.0
                        targetElementIdentifier += el*elementOffset[i]
                    targetElement = self._targetMesh.findElementByIdentifier(targetElementIdentifier)
                    result = self._targetMarkerLocation.assignMeshLocation(self._targetCache, targetElement, targetXi)
                    self._nodeIdentifier += 1
                    # add new node to matching annotation groups the previous one was in
                    for sourceAndTargetNodesetGroup in self._sourceAndTargetNodesetGroups:
                        if sourceAndTargetNodesetGroup[0].containsNode(sourceNode):
                            sourceAndTargetNodesetGroup[1].addNode(node)

        return nids, nx


    def refineAllElementsCubeStandard3d(self, numberInXi1, numberInXi2, numberInXi3):
        element = self._sourceElementiterator.next()
        while element.isValid():
            self.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            element = self._sourceElementiterator.next()
