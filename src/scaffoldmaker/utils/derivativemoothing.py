'''
Class for globally smoothing field derivatives.
'''
from __future__ import division
import math
from opencmiss.utils.maths.vectorops import add, magnitude
from opencmiss.utils.zinc.field import findOrCreateFieldGroup
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from opencmiss.zinc.result import RESULT_OK as ZINC_OK
from scaffoldmaker.utils.interpolation import DerivativeScalingMode, getCubicHermiteArcLength
from scaffoldmaker.utils.vector import setMagnitude


class EdgeCurve:
    '''
    A description of a 1D Hermite curve mapped from an element edge
    '''

    def __init__(self, expressions):
        self._expressions = expressions
        self._arcLength = 0.0

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

    def __init__(self, region, field, groupName = None, scalingMode = DerivativeScalingMode.ARITHMETIC_MEAN, editGroupName=None):
        '''
        :param groupName: Optional name of group to limit
        '''
        self._region = region
        self._field = field.castFiniteElement()
        componentsCount = self._field.getNumberOfComponents()
        assert self._field.isValid()
        self._groupName = groupName
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
        if groupName:
            self._group = self._fieldmodule.findFieldByName(groupName).castGroup()
            if not self._group.isValid():
                print('DerivativeSmoothing: Group', groupName, 'not found')
                return
            self._mesh = self._group.getFieldElementGroup(self._mesh).getMeshGroup()
            self._nodes = self._group.getFieldNodeGroup(self._nodes).getNodesetGroup()
        # edges define curves with 4 expressions for x1, d1, x2, d2, followed by the arcLength
        # each expression is a list of terms.
        # each term contains a list of global node id, value label, version, elementIdentifier, scale factor or None
        # edges are mapped from sorted start and end node
        self._edgesMap = {}
        # map global nodeid, derivative, version to list of EdgeCurve
        self._derivativeMap = {}
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
            if not eft.isValid():
                print('Ignore element',element.getIdentifier())
            if eft.isValid():
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
                    if functionType in (Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY):
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
                    print('DerivativeSmoothing: Element', element.getIdentifier(), 'has unsupported basis function types', functionTypes)
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
        globalNodes = [ element.getNode(eft, eft.getTermLocalNodeIndex(fn, 1)) for fn in locationBasisFunctions ]
        nodeIdentifiers = [ node.getIdentifier() for node in globalNodes ]
        if -1 in nodeIdentifiers:
            print('Non-nodal element', element.getIdentifier(), 'basis nodes', basisLocalNodes)
            return
        if nodeIdentifiers[0] == nodeIdentifiers[1]:
            #print('Collapsed edge on element', element.getIdentifier(), 'basis nodes', basisLocalNodes)
            return
        sortednodeIdentifiers = tuple(sorted(nodeIdentifiers))
        existingEdgeCurve = self._edgesMap.get(sortednodeIdentifiers)
        if existingEdgeCurve:
            #print('Found existing edge for global nodes', sortednodeIdentifiers, 'on element', element.getIdentifier(), 'basis nodes', basisLocalNodes)
            return
        # element basis function numbers for x1, d1, x2, d2 of Hermite curve
        functionNumbers = [ locationBasisFunctions[0], locationBasisFunctions[0] + derivativeOffset,
                            locationBasisFunctions[1], locationBasisFunctions[1] + derivativeOffset ]
        expressions = []
        scaleFactorsCount = eft.getNumberOfLocalScaleFactors()
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
        self._edgesMap[sortednodeIdentifiers] = edge

        # map from node derivative version to list of edges it is on
        for expressionIndex in (1, 3):
            expression = expressions[expressionIndex]
            if len(expression) == 1:  # handle single term only
                nodeIdentifier, nodeValueLabel, nodeVersion, scaleFactor = expression[0]
                derivativeKey = (nodeIdentifier, nodeValueLabel, nodeVersion)
                derivativeEdge = (edge, expressionIndex, math.fabs(scaleFactor))
                derivativeEdges = self._derivativeMap.get(derivativeKey)
                if derivativeEdges:
                    derivativeEdges.append(derivativeEdge)
                else:
                    self._derivativeMap[derivativeKey] = [ derivativeEdge ]
                if self._editNodesetGroup:
                    self._editNodesetGroup.addNode(self._nodes.findNodeByIdentifier(nodeIdentifier))  # so client know which nodes are modified

    def smooth(self):
        fieldcache = self._fieldmodule.createFieldcache()
        for edge in self._edgesMap.values():
            edge.evaluateArcLength(self._nodes, self._field, fieldcache)
        componentsCount = self._field.getNumberOfComponents()
        with ChangeManager(self._fieldmodule):
            for derivativeKey, derivativeEdges in self._derivativeMap.items():
                edgeCount = len(derivativeEdges)
                if edgeCount > 1:
                    nodeIdentifier, nodeValueLabel, nodeVersion = derivativeKey
                    fieldcache.setNode(self._nodes.findNodeByIdentifier(nodeIdentifier))
                    result, x = self._field.getNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, componentsCount)
                    mag = 0.0
                    for derivativeEdge in derivativeEdges:
                        edge, expressionIndex, totalScaleFactor = derivativeEdge
                        arcLength = edge.getArcLength()
                        if self._scalingMode == DerivativeScalingMode.ARITHMETIC_MEAN:
                            mag += arcLength/totalScaleFactor
                        else: # self._scalingMode == DerivativeScalingMode.HARMONIC_MEAN
                            mag += totalScaleFactor/arcLength
                    if self._scalingMode == DerivativeScalingMode.ARITHMETIC_MEAN:
                        mag /= edgeCount
                    else: # self._scalingMode == DerivativeScalingMode.HARMONIC_MEAN
                        mag = edgeCount/mag
                    if (mag <= 0.0):
                        print('Node', nodeIdentifier, 'label', label, 'version', version, 'has negative mag', mag)
                    x = setMagnitude(x, mag)
                    result = self._field.setNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, x)
            for derivativeKey, derivativeEdges in self._derivativeMap.items():
                edgeCount = len(derivativeEdges)
                if edgeCount == 1:
                    # boundary smoothing over single edge
                    nodeIdentifier, nodeValueLabel, nodeVersion = derivativeKey
                    fieldcache.setNode(self._nodes.findNodeByIdentifier(nodeIdentifier))
                    result, x = self._field.getNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, componentsCount)
                    edge, expressionIndex, totalScaleFactor = derivativeEdges[0]
                    # re-evaluate arc length so parameters are up-to-date for other end
                    arcLength = edge.evaluateArcLength(self._nodes, self._field, fieldcache)
                    otherd = edge.getParameter(3 if (expressionIndex == 1) else 1)
                    othermag = magnitude(otherd)
                    mag = (2.0*arcLength - othermag)/totalScaleFactor
                    if (mag <= 0.0):
                        print('Node', nodeIdentifier, 'label', nodeValueLabel, 'version', nodeVersion, 'has negative mag', mag)
                    x = setMagnitude(x, mag)
                    fieldcache.setNode(self._nodes.findNodeByIdentifier(nodeIdentifier))  # need to set again as changed node in edge.evaluateArcLength
                    result = self._field.setNodeParameters(fieldcache, -1, nodeValueLabel, nodeVersion, x)
