'''
Class for refining a mesh from one region to another.
Created on April 4, 2018

@author: Richard Christie
'''

from scaffoldmaker.utils.octree import Octree
from scaffoldmaker.utils.zinc_utils import *
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.result import RESULT_OK as ZINC_OK

class MeshRefinement:
    '''
    Class for refining a mesh from one region to another.
    '''

    def __init__(self, sourceRegion, targetRegion):
        '''
        Assumes targetRegion is empty.
        '''
        self._sourceRegion = sourceRegion
        self._sourceFm = sourceRegion.getFieldmodule()
        self._sourceCache = self._sourceFm.createFieldcache()
        self._sourceCoordinates = getOrCreateCoordinateField(self._sourceFm)
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
        self._sourceFm.endChange()
        self._sourceMesh = self._sourceFm.findMeshByDimension(3)
        self._octree = Octree(minimums, maximums)

        self._targetRegion = targetRegion
        self._targetFm = targetRegion.getFieldmodule()
        self._targetFm.beginChange()
        self._targetCache = self._targetFm.createFieldcache()
        self._targetCoordinates = getOrCreateCoordinateField(self._targetFm)

        self._targetNodes = self._targetFm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        self._nodetemplate = self._targetNodes.createNodetemplate()
        self._nodetemplate.defineField(self._targetCoordinates)

        self._targetMesh = self._targetFm.findMeshByDimension(3)
        self._targetBasis = self._targetFm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
        self._targetEft = self._targetMesh.createElementfieldtemplate(self._targetBasis)
        self._targetElementtemplate = self._targetMesh.createElementtemplate()
        self._targetElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = self._targetElementtemplate.defineField(self._targetCoordinates, -1, self._targetEft)

        self._nodetemplate.defineField(self._targetCoordinates)

        self._nodeIdentifier = 1
        self._elementIdentifier = 1

    def __del__(self):
        self._targetFm.endChange()

    def refineElementCubeStandard3d(self, sourceElement, numberInXi1, numberInXi2, numberInXi3):
        # create nodes
        nids = []
        xi = [ 0.0, 0.0, 0.0 ]
        for k in range(numberInXi3 + 1):
            xi[2] = k/numberInXi3
            for j in range(numberInXi2 + 1):
                xi[1] = j/numberInXi2
                for i in range(numberInXi1 + 1):
                    xi[0] = i/numberInXi1
                    self._sourceCache.setMeshLocation(sourceElement, xi)
                    result, x = self._sourceCoordinates.evaluateReal(self._sourceCache, 3)
                    nodeId = self._octree.findObjectByCoordinates(x)
                    if nodeId is None:
                        node = self._targetNodes.createNode(self._nodeIdentifier, self._nodetemplate)
                        self._targetCache.setNode(node)
                        result = self._targetCoordinates.setNodeParameters(self._targetCache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        nodeId = self._nodeIdentifier
                        self._octree.addObjectAtCoordinates(x, nodeId)
                        self._nodeIdentifier += 1
                    nids.append(nodeId)
        # create elements
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

    def refineAllElementsCubeStandard3d(self, numberInXi1, numberInXi2, numberInXi3):
        elementIter = self._sourceMesh.createElementiterator()
        element = elementIter.next()
        while element.isValid():
            self.refineElementCubeStandard3d(element, numberInXi1, numberInXi2, numberInXi3)
            element = elementIter.next()
