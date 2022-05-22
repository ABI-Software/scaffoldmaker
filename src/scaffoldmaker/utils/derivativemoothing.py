'''
Class for globally smoothing field derivatives.
'''
from __future__ import division

import math

from opencmiss.maths.vectorops import magnitude
from opencmiss.utils.zinc.field import findOrCreateFieldGroup
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from scaffoldmaker.utils.interpolation import DerivativeScalingMode, getCubicHermiteArcLength, interpolateHermiteLagrangeDerivative, interpolateLagrangeHermiteDerivative
from scaffoldmaker.utils.vector import setMagnitude


class EdgeCurve:
    '''
    A description of a 1D Hermite curve mapped from an element edge
    '''

    def __init__(self, expressions):
        self._expressions = expressions
        self._arcLength = None
        self._lastArcLength = 0.0
        self._parameters = None

    def evaluateArcLength(self, nodes, field, fieldcache):
        componentsCount = field.getNumberOfComponents()
        self._parameters = []
        for expression in self._expressions:
            sumx = None
            for term in expression:
                nodeIdentifier, nodeValueLabel, nodeVersion, scaleFactor = term
                fieldcache.setNode(nodes.findNodeByIdentifier(nodeIdentifier))
                result, x = field.getNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, componentsCount)
                if scaleFactor:
                    x = [ v*scaleFactor for v in x ]
                if sumx:
                    sumx = [ (sumx[c] + x[c]) for c in range(componentsCount) ]
                else:
                    sumx = x
            self._parameters.append(sumx)
        self._arcLength = getCubicHermiteArcLength(*self._parameters)
        return self._arcLength

    def getArcLength(self):
        '''
        Caller must have called evaluateArcLength!
        '''
        return self._arcLength

    def getLastArcLength(self):
        return self._lastArcLength

    def updateLastArcLength(self):
        self._lastArcLength = self._arcLength

    def getDelta(self):
        '''
        Caller must have called evaluateArcLength!
        '''
        x1 = self._parameters[0]
        x2 = self._parameters[2]
        componentsCount = len(x1)
        return [ (x2[c] - x1[c]) for c in range(componentsCount) ]

    def getExpression(self, expressionIndex):
        '''
        Caller must have called evaluateArcLength!
        :param parameterIndex: 0 = start x, 1 = start d, 2 = end x, 3 = end d
        '''
        return self._expressions[expressionIndex]

    def getParameter(self, parameterIndex):
        '''
        Caller must have called evaluateArcLength!
        :param parameterIndex: 0 = start x, 1 = start d, 2 = end x, 3 = end d
        '''
        return self._parameters[parameterIndex]

class DerivativeSmoothing:
    '''
    Class for globally smoothing field derivatives.
    '''

    cubeEdgeLocalNodes = [
        [ [ 1, 2 ], [ 3, 4 ], [ 5, 6 ], [ 7, 8 ] ],
        [ [ 1, 3 ], [ 2, 4 ], [ 5, 7 ], [ 6, 8 ] ],
        [ [ 1, 5 ], [ 2, 6 ], [ 3, 7 ], [ 4, 8 ] ] ]
    squareEdgeLocalNodes = [
        [ [ 1, 2 ], [ 3, 4 ] ],
        [ [ 1, 3 ], [ 2, 4 ] ] ]

    def __init__(self, region, field, selectionGroupName=None, scalingMode=DerivativeScalingMode.ARITHMETIC_MEAN, editGroupName=None):
        '''
        :param selectionGroupName: Optional name of group to limit smoothing to nodes in group.
        Only element edges including those nodes are included in smoothing.
        '''
        self._region = region
        self._field = field.castFiniteElement()
        componentsCount = self._field.getNumberOfComponents()
        assert self._field.isValid()
        self._selectionGroupName = selectionGroupName
        self._selectionGroup = None
        self._selectionNodes = None
        self._scalingMode = scalingMode
        self._editGroupName = editGroupName
        self._fieldmodule = self._region.getFieldmodule()
        for dimension in range(3,0,-1):
            self._mesh = self._fieldmodule.findMeshByDimension(dimension)
            if self._mesh.getSize() > 0:
                break
        self._nodes = self._fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # edited nodes are added to the edit nodeset group, if group is supplied
        if editGroupName:
            editGroup = findOrCreateFieldGroup(self._fieldmodule, editGroupName, managed=True)
            editNodeGroup = editGroup.getFieldNodeGroup(self._nodes)
            if not editNodeGroup.isValid():
                editNodeGroup = editGroup.createFieldNodeGroup(self._nodes)
            self._editNodesetGroup = editNodeGroup.getNodesetGroup()
        else:
            self._editNodesetGroup = None
        # edges define curves with 4 expressions for x1, d1, x2, d2, followed by the arcLength
        # each expression is a list of terms.
        # each term contains a list of global node id, value label, version, element identifier, scale factor or None
        # edges are mapped from sorted start and end node
        self._edgesMap = {}
        # map global nodeid, derivative, version to list of EdgeCurve
        self._derivativeMap = {}
        if selectionGroupName:
            self._selectionGroup = self._fieldmodule.findFieldByName(selectionGroupName).castGroup()
            if not self._selectionGroup.isValid():
                print('DerivativeSmoothing: Selection group not found')
                return
            self._selectionNodes = self._selectionGroup.getFieldNodeGroup(self._nodes).getNodesetGroup()
            if (not self._selectionNodes.isValid()) or (self._selectionNodes.getSize() == 0):
                print('DerivativeSmoothing: No nodes selected for smoothing')
                return
        self.addElementEdges()

    def addElementEdges(self):
        '''
        Compile edge and derivative maps over elements.
        '''
        dimension = self._mesh.getDimension()
        elementIter = self._mesh.createElementiterator()
        element = elementIter.next()
        while element.isValid():
            eft = element.getElementfieldtemplate(self._field, -1)
            useElement = True
            if not eft.isValid():
                eft1 = element.getElementfieldtemplate(self._field, 1)
                if eft1.isValid():
                    print('DerivativeSmoothing: Skip element', element.getIdentifier(), 'with per-component element field, which is not implemented')
                #else:
                #    print('DerivativeSmoothing: Skip element', element.getIdentifier(), 'as field is not defined on it')
                useElement = False
            elif self._selectionNodes:
                # check element contains at least one nodes in selection:
                localNodesCount = eft.getNumberOfLocalNodes()
                for n in range(1, localNodesCount + 1):
                    node = element.getNode(eft, n)
                    if self._selectionNodes.containsNode(node):
                        break
                else:
                    #print('DerivativeSmoothing: Skip element', element.getIdentifier(), 'as it contains no selected nodes')
                    useElement = False
            if useElement:
                elementbasis = eft.getElementbasis()
                functionTypes = [ elementbasis.getFunctionType(d) for d in range(1, dimension + 1) ]
                fn = 1
                basisLocalNodeCount = elementbasis.getNumberOfNodes()
                basisLocalNodeFunctions = [ 1 ]
                for n in range(1, basisLocalNodeCount):
                    basisLocalNodeFunctions.append(basisLocalNodeFunctions[n - 1] + elementbasis.getNumberOfFunctionsPerNode(n));
                useDirections = []
                derivativeOffsets = []
                derivativeOffset = 1
                for functionType in functionTypes:
                    if functionType in (Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE,):  # Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY):
                        useDirections.append(True)
                        if (functionType == Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE) and (len(derivativeOffsets) == 2):
                            derivativeOffset += 1
                        derivativeOffsets.append(derivativeOffset)
                        derivativeOffset += 1
                    elif functionType == Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE:
                        useDirections.append(False)
                        derivativeOffsets.append(None)
                shapeType = element.getShapeType()
                if len(useDirections) != dimension:
                    print('DerivativeSmoothing: Skip element', element.getIdentifier(), 'with unsupported basis function types', functionTypes)
                elif shapeType == Element.SHAPE_TYPE_CUBE:
                    for d in range(dimension):
                        if useDirections[d]:
                            for localNodes in self.cubeEdgeLocalNodes[d]:
                                self.addEdgeCurve(element, eft, basisLocalNodeFunctions, localNodes, derivativeOffsets[d])
                elif shapeType == Element.SHAPE_TYPE_SQUARE:
                    for d in range(dimension):
                        if useDirections[d]:
                            for localNodes in self.squareEdgeLocalNodes[d]:
                                self.addEdgeCurve(element, eft, basisLocalNodeFunctions, localNodes, derivativeOffsets[d])
                elif shapeType == Element.SHAPE_TYPE_LINE:
                    if useDirections[0]:
                        self.addEdgeCurve(element, eft, basisLocalNodeFunctions, [ 1, 2 ], derivativeOffsets[0])
                else:
                    print('DerivativeSmoothing: Element', element.getIdentifier(), 'has unsupported shape type')
            element = elementIter.next()

    def addEdgeCurve(self, element, eft, basisLocalNodeFunctions, basisLocalNodes, derivativeOffset):
        '''
        Add the EdgeCurve object representing the Hermite curve between the local nodes.
        Reused existing curve for matching pairs of global nodes.
        Build map of node-derivative-version to list of edges.
        :param basisLocalNodes: Local node indexes of start and end for the raw basis.
        :param derivativeOffset: Offset from value parameter at node to derivative function.
        '''
        #print('element', element.getIdentifier(), 'fns', basisLocalNodeFunctions, 'basisLocalNodes', basisLocalNodes, 'derivativeOffset', derivativeOffset)
        locationBasisFunctions = [ basisLocalNodeFunctions[n - 1] for n in basisLocalNodes ]
        edgeNodes = [ element.getNode(eft, eft.getTermLocalNodeIndex(fn, 1)) for fn in locationBasisFunctions ]
        if self._selectionNodes:
            for node in edgeNodes:
                if self._selectionNodes.containsNode(node):
                    break
            else:
                return  # edge does not contain any nodes in selection
        nodeIdentifiers = [ node.getIdentifier() for node in edgeNodes ]
        if -1 in nodeIdentifiers:
            print('Non-nodal element', element.getIdentifier(), 'basis nodes', basisLocalNodes)
            return
        if nodeIdentifiers[0] == nodeIdentifiers[1]:
            #print('Collapsed edge on element', element.getIdentifier(), 'basis nodes', basisLocalNodes)
            return
        sortedNodeIdentifiers = tuple(sorted(nodeIdentifiers))
        existingEdgeCurve = self._edgesMap.get(sortedNodeIdentifiers)
        if existingEdgeCurve:
            #print('Found existing edge for global nodes', sortedNodeIdentifiers, 'on element', element.getIdentifier(), 'basis nodes', basisLocalNodes)
            return
        # element basis function numbers for x1, d1, x2, d2 of Hermite curve
        functionNumbers = [ locationBasisFunctions[0], locationBasisFunctions[0] + derivativeOffset,
                            locationBasisFunctions[1], locationBasisFunctions[1] + derivativeOffset ]
        expressions = []
        scaleFactors = []
        scaleFactorsCount = eft.getNumberOfLocalScaleFactors()
        if scaleFactorsCount > 0:
            # handle zinc returning single value as a scalar, change to list for consistency
            result, scaleFactors = element.getScaleFactors(eft, scaleFactorsCount)
            if not isinstance(scaleFactors, list):
                scaleFactors = [scaleFactors]
        for fn in functionNumbers:
            termCount = eft.getFunctionNumberOfTerms(fn)
            expression = []
            for t in range(1, termCount + 1):
                nodeIdentifier = element.getNode(eft, eft.getTermLocalNodeIndex(fn, t)).getIdentifier()
                valueLabel = eft.getTermNodeValueLabel(fn, t)
                nodeVersion = eft.getTermNodeVersion(fn, t)
                totalScaleFactor = 1.0
                if scaleFactorsCount > 0:
                    scaleFactorIndexesCount = eft.getTermScaling(fn, t, 0)[0]
                    if scaleFactorIndexesCount:
                        result, scaleFactorIndexes = eft.getTermScaling(fn, t, scaleFactorIndexesCount)
                        if not isinstance(scaleFactorIndexes, list):
                            scaleFactorIndexes = [scaleFactorIndexes]
                        for scaleFactorIndex in scaleFactorIndexes:
                            totalScaleFactor *= scaleFactors[scaleFactorIndex - 1]
                term = [ nodeIdentifier, valueLabel, nodeVersion, totalScaleFactor ]
                expression.append(term)
            expressions.append(expression)
        edge = EdgeCurve(expressions)
        self._edgesMap[sortedNodeIdentifiers] = edge

        # map from node derivative version to list of edges it is on
        for expressionIndex in (1, 3):
            expression = expressions[expressionIndex]
            if len(expression) == 1:  # handle single term only
                nodeIdentifier, nodeValueLabel, nodeVersion, totalScaleFactor = expression[0]
                if self._selectionNodes:
                    if not self._selectionNodes.containsNode(self._nodes.findNodeByIdentifier(nodeIdentifier)):
                        continue
                derivativeKey = (nodeIdentifier, nodeValueLabel, nodeVersion)
                derivativeEdge = (edge, expressionIndex, totalScaleFactor)
                derivativeEdges = self._derivativeMap.get(derivativeKey)
                if derivativeEdges:
                    derivativeEdges.append(derivativeEdge)
                else:
                    self._derivativeMap[derivativeKey] = [ derivativeEdge ]

    def smooth(self, updateDirections=False, maxIterations=10, arcLengthTolerance=1.0E-6):
        '''
        :param maxIterations: Maximum iterations before stopping if not converging.
        :param arcLengthTolerance: Ratio of difference in arc length from last iteration
        divided by current arc length under which convergence is achieved. Required to
        be met by every element edge.
        '''
        if not self._derivativeMap:
            return  # no nodes being smoothed
        componentsCount = self._field.getNumberOfComponents()
        with ChangeManager(self._fieldmodule):
            fieldcache = self._fieldmodule.createFieldcache()
            for iter in range(maxIterations + 1):
                converged = True
                for edge in self._edgesMap.values():
                    lastArcLength = edge.getLastArcLength()
                    arcLength = edge.evaluateArcLength(self._nodes, self._field, fieldcache)
                    edge.updateLastArcLength()
                    if (math.fabs(arcLength - lastArcLength)/arcLength) > arcLengthTolerance:
                        converged = False
                if converged:
                    print('Derivative smoothing: Converged after', iter, 'iterations.')
                    break
                elif (iter == maxIterations):
                    print('Derivative smoothing: Stopping after', maxIterations, 'iterations without converging.')
                    break
                for derivativeKey, derivativeEdges in self._derivativeMap.items():
                    edgeCount = len(derivativeEdges)
                    if edgeCount > 1:
                        nodeIdentifier, nodeValueLabel, nodeVersion = derivativeKey
                        fieldcache.setNode(self._nodes.findNodeByIdentifier(nodeIdentifier))
                        if updateDirections:
                            x = [ 0.0 for _ in range(componentsCount) ]
                        else:
                            result, x = self._field.getNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, componentsCount)
                        mag = 0.0
                        for derivativeEdge in derivativeEdges:
                            edge, expressionIndex, totalScaleFactor = derivativeEdge
                            arcLength = edge.getArcLength()
                            if updateDirections:
                                delta = edge.getDelta()
                                if totalScaleFactor < 0.0:
                                    delta = [ -d for d in delta ]
                                for c in range(componentsCount):
                                    x[c] += delta[c] / arcLength
                            if self._scalingMode == DerivativeScalingMode.ARITHMETIC_MEAN:
                                mag += arcLength/math.fabs(totalScaleFactor)
                            else: # self._scalingMode == DerivativeScalingMode.HARMONIC_MEAN
                                mag += math.fabs(totalScaleFactor)/arcLength
                        if self._scalingMode == DerivativeScalingMode.ARITHMETIC_MEAN:
                            mag /= edgeCount
                        else: # self._scalingMode == DerivativeScalingMode.HARMONIC_MEAN
                            mag = edgeCount/mag
                        if (mag <= 0.0):
                            print('Node', nodeIdentifier, 'value', nodeValueLabel, 'version', nodeVersion, \
                                  'has negative mag', mag)
                        x = setMagnitude(x, mag)
                        result = self._field.setNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, x)
                for derivativeKey, derivativeEdges in self._derivativeMap.items():
                    edgeCount = len(derivativeEdges)
                    if edgeCount == 1:
                        # boundary smoothing over single edge
                        nodeIdentifier, nodeValueLabel, nodeVersion = derivativeKey
                        edge, expressionIndex, totalScaleFactor = derivativeEdges[0]
                        # re-evaluate arc length so parameters are up-to-date for other end
                        arcLength = edge.evaluateArcLength(self._nodes, self._field, fieldcache)
                        fieldcache.setNode(self._nodes.findNodeByIdentifier(nodeIdentifier))  # since changed by evaluateArcLength
                        otherExpressionIndex = 3 if (expressionIndex == 1) else 1
                        otherd = edge.getParameter(otherExpressionIndex)
                        if updateDirections:
                            thisx = edge.getParameter(expressionIndex - 1)
                            otherx = edge.getParameter(otherExpressionIndex - 1)
                            bothEndsOnBoundary = False
                            otherExpression = edge.getExpression(otherExpressionIndex)
                            if len(otherExpression) == 1:
                                otherNodeIdentifier, otherValueLabel, otherNodeVersion, otherTotalScaleFactor = otherExpression[0]
                                otherDerivativeKey = (otherNodeIdentifier, otherValueLabel, otherNodeVersion)
                                otherDerivativeEdges = self._derivativeMap.get(otherDerivativeKey)
                                bothEndsOnBoundary = (otherDerivativeEdges is not None) and (len(otherDerivativeEdges) == 1)
                            if bothEndsOnBoundary:
                                if expressionIndex == 1:
                                    x = [ (otherx[c] - thisx[c]) for c in range(componentsCount) ]
                                else:
                                    x = [ (thisx[c] - otherx[c]) for c in range(componentsCount) ]
                            else:
                                if expressionIndex == 1:
                                    x = interpolateLagrangeHermiteDerivative(thisx, otherx, otherd, 0.0)
                                else:
                                    x = interpolateHermiteLagrangeDerivative(otherx, otherd, thisx, 1.0)
                            x = [ d/totalScaleFactor for d in x ]
                        else:
                            result, x = self._field.getNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, componentsCount)
                        othermag = magnitude(otherd)
                        mag = (2.0*arcLength - othermag)/math.fabs(totalScaleFactor)
                        if (mag <= 0.0):
                            print('Derivative smoothing: Node', nodeIdentifier, 'label', nodeValueLabel, 'version', nodeVersion, 'has negative magnitude', mag)
                        x = setMagnitude(x, mag)
                        result = self._field.setNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, x)
            # record modified nodes while ChangeManager is in effect
            if self._editNodesetGroup:
                for derivativeKey in self._derivativeMap:
                    self._editNodesetGroup.addNode((self._nodes.findNodeByIdentifier(derivativeKey[0])))
            del fieldcache
