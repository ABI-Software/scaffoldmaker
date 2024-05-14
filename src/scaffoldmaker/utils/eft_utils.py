'''
Utility functions for element field templates shared by mesh generators.
'''
from cmlibs.maths.vectorops import add, dot, magnitude, mult, sub
from cmlibs.zinc.element import Elementbasis, Elementfieldtemplate
from cmlibs.zinc.node import Node
from cmlibs.zinc.result import RESULT_OK
import math


def getEftTermScaling(eft, functionIndex, termIndex):
    '''
    Convenience function to get the scale factor indexes scaling a term as a list.
    :eft:  Element field template to query.
    :functionIndex:  Function index from 1 to number of basis functions.
    :termIndex:  Term index from 1 to function number of terms.
    :return:  List containing any scale factor indexes scaling term. Empty if unscaled.
    '''
    scaleFactorCount, scaleFactorIndexes = eft.getTermScaling(functionIndex, termIndex, 1)
    if scaleFactorCount < 0:
        if eft.getNumberOfLocalScaleFactors() > 0:
            print('getEftTermScaling function', functionIndex, ' term', termIndex, ' scaleFactorCount ', scaleFactorCount)
        return []
    if scaleFactorCount == 0:
        return []
    if scaleFactorCount == 1:
        return [ scaleFactorIndexes ]
    scaleFactorCount, scaleFactorIndexes = eft.getTermScaling(functionIndex, termIndex, scaleFactorCount)
    return scaleFactorIndexes


def mapEftFunction1Node1Term(eft, function, localNode, valueLabel, version, scaleFactors):
    '''
    Set function of eft to map valueLabel, version from localNode with scaleFactors
    '''
    eft.setFunctionNumberOfTerms(function, 1)
    eft.setTermNodeParameter(function, 1, localNode, valueLabel, version)
    eft.setTermScaling(function, 1, scaleFactors)


def mapEftFunction1Node2Terms(eft, function, localNode, valueLabel1, version1, scaleFactors1, valueLabel2, version2, scaleFactors2):
    '''
    Set function of eft to sum 2 terms the respective valueLabels, versions from localNode with scaleFactors
    '''
    eft.setFunctionNumberOfTerms(function, 2)
    eft.setTermNodeParameter(function, 1, localNode, valueLabel1, version1)
    eft.setTermScaling(function, 1, scaleFactors1)
    eft.setTermNodeParameter(function, 2, localNode, valueLabel2, version2)
    eft.setTermScaling(function, 2, scaleFactors2)


def remapEftLocalNodes(eft, newNodeCount, localNodeIndexes):
    '''
    Remaps current local nodes to the new local ids, changing number of local nodes.
    Assumes node parameters are in use.
    :param localNodeIds: new local node identifiers starting at 1 for each current local node id - 1 referenced.
    '''
    functionCount = eft.getNumberOfFunctions()
    for f in range(1, functionCount + 1):
        termCount = eft.getFunctionNumberOfTerms(f)
        for t in range(1, termCount + 1):
            eft.setTermNodeParameter(f, t, localNodeIndexes[eft.getTermLocalNodeIndex(f, t) - 1], eft.getTermNodeValueLabel(f, t), eft.getTermNodeVersion(f, t))
    eft.setNumberOfLocalNodes(newNodeCount)


def remapEftNodeValueLabel(eft, localNodeIndexes, fromValueLabel, expressionTerms):
    '''
    Remap all uses of the given valueLabels to the expressionTerms.
    Note: Assumes valueLabel is currently single term and unscaled!
    :param localNodeIndexes:  List of local node indexes >= 1 to remap at.
    :param fromValueLabel:  Node value label to be remapped.
    :param expressionTerms: List of (valueLabel, scaleFactorIndexesList ) to remap to.
        e.g. [ (Node.VALUE_LABEL_D_DS2, []), (Node.VALUE_LABEL_D_DS3, [5, 6]) ]
    '''
    functionCount = eft.getNumberOfFunctions()
    for f in range(1, functionCount + 1):
        if eft.getFunctionNumberOfTerms(f) == 1:
            localNodeIndex = eft.getTermLocalNodeIndex(f, 1)
            if ((localNodeIndex in localNodeIndexes) and (eft.getTermNodeValueLabel(f, 1) == fromValueLabel) and
                    (not getEftTermScaling(eft, f, 1))):
                termCount = len(expressionTerms)
                eft.setFunctionNumberOfTerms(f, termCount)
                version = eft.getTermNodeVersion(f, 1)
                for t in range(1, termCount + 1):
                    expressionTerm = expressionTerms[t - 1]
                    eft.setTermNodeParameter(f, t, localNodeIndex, expressionTerm[0], version)
                    if expressionTerm[1]:
                        eft.setTermScaling(f, t, expressionTerm[1])


def remapEftNodeValueLabelVersion(eft, localNodeIndexes, fromValueLabel, expressionTerms):
    '''
    Remap all uses of the given valueLabels to the expressionTerms.
    Note: Assumes valueLabel is currently single term and unscaled!
    :param localNodeIndexes:  List of local node indexes >= 1 to remap at.
    :param fromValueLabel:  Node value label to be remapped.
    :param expressionTerms: List of (valueLabel, version, scaleFactorIndexesList) to remap to.
        e.g. [ (Node.VALUE_LABEL_D_DS2, 1, []), (Node.VALUE_LABEL_D_DS3, 2, [5, 6]) ]
    '''
    functionCount = eft.getNumberOfFunctions()
    for f in range(1, functionCount + 1):
        if eft.getFunctionNumberOfTerms(f) == 1:
            localNodeIndex = eft.getTermLocalNodeIndex(f, 1)
            if ((localNodeIndex in localNodeIndexes) and (eft.getTermNodeValueLabel(f, 1) == fromValueLabel) and
                    (not getEftTermScaling(eft, f, 1))):
                termCount = len(expressionTerms)
                eft.setFunctionNumberOfTerms(f, termCount)
                for t in range(1, termCount + 1):
                    expressionTerm = expressionTerms[t - 1]
                    eft.setTermNodeParameter(f, t, localNodeIndex, expressionTerm[0], expressionTerm[1])
                    if expressionTerm[2]:
                        eft.setTermScaling(f, t, expressionTerm[2])


def remapEftNodeValueLabelsVersion(eft, localNodeIndexes, valueLabels, version):
    '''
    Remap all uses of the given valueLabels to use the version.
    :param localNodeIndexes:  List of local node indexes >= 1 to remap at.
    :param valueLabels:  List of node value labels to be remapped.
    :param version: Version >= 1
    '''
    functionCount = eft.getNumberOfFunctions()
    for f in range(1, functionCount + 1):
        termCount = eft.getFunctionNumberOfTerms(f)
        for t in range(1, termCount + 1):
            localNodeIndex = eft.getTermLocalNodeIndex(f, t)
            valueLabel = eft.getTermNodeValueLabel(f, t)
            if (localNodeIndex in localNodeIndexes) and (valueLabel in valueLabels):
                result = eft.setTermNodeParameter(f, t, localNodeIndex, valueLabel, version)
                #print('remap result', result)


def remapEftNodeValueLabelWithNodes(eft, localNodeIndex, fromValueLabel, expressionTerms):
    '''
    Remap all uses of the given valueLabel to the expressionTerms.
    Note: Assumes valueLabel is currently single term and unscaled!
    :param localNodeIndex:  Local node index >= 1 to remap at.
    :param fromValueLabel:  Node value label to be remapped.
    :param expressionTerms: List of (localNodeIndex, valueLabel, scaleFactorIndexesList ) to remap to.
        e.g. [ (6, Node.VALUE_LABEL_VALUE, [1]), (7, Node.VALUE_LABEL_VALUE, []) ]
    '''
    functionCount = eft.getNumberOfFunctions()
    for f in range(1, functionCount + 1):
        if eft.getFunctionNumberOfTerms(f) == 1:
            if (eft.getTermLocalNodeIndex(f, 1) == localNodeIndex) and (eft.getTermNodeValueLabel(f, 1) == fromValueLabel) and (not getEftTermScaling(eft, f, 1)):
                termCount = len(expressionTerms)
                eft.setFunctionNumberOfTerms(f, termCount)
                version = eft.getTermNodeVersion(f, 1)
                for t in range(1, termCount + 1):
                    expressionTerm = expressionTerms[t - 1]
                    eft.setTermNodeParameter(f, t, expressionTerm[0], expressionTerm[1], version)
                    if expressionTerm[2]:
                        eft.setTermScaling(f, t, expressionTerm[2])


def scaleEftNodeValueLabels(eft, localNodeIndexes, valueLabels, addScaleFactorIndexes):
    '''
    Multiply all uses of the given ValueLabels by scale factor at scaleFactorIndex.
    Handles general maps.
    :param localNodeIndexes:  List of local node indexes >= 1 to scale value label.
    :param valueLabels:  List of node value label to be scaled.
    :param addScaleFactorIndexes: List of valid scale factor indexes >= 1 to append to existing scale factors.
    '''
    functionCount = eft.getNumberOfFunctions()
    for f in range(1, functionCount + 1):
        termCount = eft.getFunctionNumberOfTerms(f)
        for t in range(1, termCount + 1):
            if eft.getTermLocalNodeIndex(f, t) in localNodeIndexes:
                if eft.getTermNodeValueLabel(f, t) in valueLabels:
                    scaleFactorIndexes = getEftTermScaling(eft, f, t)
                    scaleFactorIndexes += addScaleFactorIndexes
                    eft.setTermScaling(f, t, scaleFactorIndexes)


def setEftScaleFactorIds(eft, globalScaleFactorIds, nodeScaleFactorIds):
    '''
    Set general followed by node scale factor identifiers.
    '''
    eft.setNumberOfLocalScaleFactors(len(globalScaleFactorIds) + len(nodeScaleFactorIds))
    s = 1
    for id in globalScaleFactorIds:
        eft.setScaleFactorType(s, Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL)
        eft.setScaleFactorIdentifier(s, id)
        s += 1
    for id in nodeScaleFactorIds:
        eft.setScaleFactorType(s, Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL)
        eft.setScaleFactorIdentifier(s, id)
        s += 1


def createEftElementSurfaceLayer(elementIn, eftIn, eftfactory, eftStd, removeNodeValueLabel=None):
    """
    Create eft for layer above xi3 = 1, from either tricubic hermite or bicubic hermite linear.
    Extracts mapping information from top xi3=1 surface of elementIn with eftIn.
    Returns eft for new element and scalefactors for it.
    :param elementIn: The element to build on.
    :param eftIn: The eft for the base element field. Assumed validated.
    :param eftfactory: The eft factory to create new eft from, if needed.
    Either eftfactory_tricubichermite or eftfactory_bicubichermitelinear.
    :param eftStd: The standard eft for the new element as created by eftfactory.createEftNoCrossDerivatives.
    Assumed unscaled, without cross derivatives, max 1 version and not collapsed. Never modified.
    :param removeNodeValueLabel: If not None, remove mappings from the value label in newly created layer.
    :return: eftStd or new eft if remapped, scalefactors
    """
    eft = eftStd
    scalefactors = None
    scaleFactorCountIn = eftIn.getNumberOfLocalScaleFactors()
    scaleFactorsUsed = [False]*scaleFactorCountIn  # flag scale factors used on top surface of element
    functionCountIn = eftIn.getNumberOfFunctions()
    halfFunctionCountIn = functionCountIn // 2
    elementBasis = eftfactory.getElementbasis()
    functionCountOut = elementBasis.getNumberOfFunctions()
    halfFunctionCountOut = functionCountOut // 2
    assert functionCountOut == eftStd.getNumberOfFunctions()
    assert (functionCountIn in (32, 64)) and (functionCountOut in (32, 64)) and (functionCountIn > functionCountOut)
    assert eftIn.getNumberOfLocalNodes() == 8
    assert eftStd.getNumberOfLocalNodes() == 8
    functionsPerNodeIn = functionCountIn // 8
    functionsPerNodeOut = functionCountOut // 8
    for i in range(halfFunctionCountOut):
        f = 1 + i
        fIn = halfFunctionCountIn + 1 + ((i // functionsPerNodeOut) * functionsPerNodeIn) + (i % functionsPerNodeOut)
        termCountIn = eftIn.getFunctionNumberOfTerms(fIn)
        noRemap = True
        newt = 0
        for j in range(termCountIn):
            t = 1 + j
            localNodeIndex = eftIn.getTermLocalNodeIndex(fIn, t) - 4
            nodeValueLabel = eftIn.getTermNodeValueLabel(fIn, t)
            nodeVersion = eftIn.getTermNodeVersion(fIn, t)
            scalingCount, scalingIndexes = eftIn.getTermScaling(fIn, t, 0)
            if scalingCount > 0:
                scalingIndexes = eftIn.getTermScaling(fIn, t, scalingCount)[1]
                #if (scalingCount > 1) or (scalingIndexes != 1):
                #    print("Scaling", scalingIndexes)
                if scalingCount == 1:
                    scalingIndexes = [scalingIndexes]
                for scalingIndex in scalingIndexes:
                    scaleFactorsUsed[scalingIndex - 1] = True
            if ((t > 1) or (scalingCount > 0)
                    or (eft.getTermLocalNodeIndex(f, 1) != localNodeIndex)
                    or (eft.getTermNodeValueLabel(f, 1) != nodeValueLabel)
                    or (eft.getTermNodeVersion(f, 1) != nodeVersion)):
                noRemap = False
                if eft is eftStd:
                    eft = eftfactory.createEftNoCrossDerivatives()
                    # initially use same number of scale factors as eftIn, with default type and identifiers
                    # later set type and identifiers, duplicating node scale factor indexes for top row
                    eft.setNumberOfLocalScaleFactors(scaleFactorCountIn)
                if t == 2:
                    eft.setFunctionNumberOfTerms(f, termCountIn)
                eft.setTermNodeParameter(f, t, localNodeIndex, nodeValueLabel, nodeVersion)
                if scalingCount > 0:
                    eft.setTermScaling(f, t, scalingIndexes)
                if nodeValueLabel != removeNodeValueLabel:
                    newt += 1
                    if newt > 1:
                        eft.setFunctionNumberOfTerms(f + halfFunctionCountOut, newt)
                    # top surface is the same, with node offset
                    assert localNodeIndex <= 4
                    assert RESULT_OK == eft.setTermNodeParameter(f + halfFunctionCountOut, newt, localNodeIndex + 4,
                                                                 nodeValueLabel, nodeVersion), localNodeIndex + 4
                    if scalingCount > 0:
                        eft.setTermScaling(f + halfFunctionCountOut, newt, scalingIndexes)
            else:
                newt = t
        assert noRemap or (newt > 0), "Element " + str(elementIn.getIdentifier()) + " f " + str(f) + \
            " terms " + str(termCountIn)
    if True in scaleFactorsUsed:
        # get used scale factors, reorder indexes
        result, scalefactorsIn = elementIn.getScaleFactors(eftIn, scaleFactorCountIn)
        if isinstance(scalefactorsIn, float):
            scalefactorsIn = [scalefactorsIn]
        scalefactors = []
        oldToNewScaleFactorIndex = {}
        globalScaleFactorIds = []
        for i in range(scaleFactorCountIn):
            if scaleFactorsUsed[i]:
                sf = i + 1
                scaleFactorType = eftIn.getScaleFactorType(sf)
                if scaleFactorType == Elementfieldtemplate.SCALE_FACTOR_TYPE_GLOBAL_GENERAL:
                    globalScaleFactorIds.append(eftIn.getScaleFactorIdentifier(sf))
                    scalefactors.append(scalefactorsIn[i])
                    oldToNewScaleFactorIndex[sf] = len(scalefactors)
                else:
                    # only other type current supported here:
                    assert scaleFactorType == Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL
        nodeScaleFactorIds = []
        for i in range(scaleFactorCountIn):
            if scaleFactorsUsed[i]:
                sf = i + 1
                scaleFactorType = eftIn.getScaleFactorType(sf)
                if scaleFactorType == Elementfieldtemplate.SCALE_FACTOR_TYPE_NODE_GENERAL:
                    nodeScaleFactorIds.append(eftIn.getScaleFactorIdentifier(sf))
                    scalefactors.append(scalefactorsIn[i])
                    oldToNewScaleFactorIndex[sf] = len(scalefactors)
        # must duplicate node scale factors as separate for top layer
        setEftScaleFactorIds(eft, globalScaleFactorIds, nodeScaleFactorIds + nodeScaleFactorIds)
        globalScaleFactorCount = len(globalScaleFactorIds)
        nodeScaleFactorCount = len(nodeScaleFactorIds)
        if nodeScaleFactorCount:
            scalefactors += scalefactors[-nodeScaleFactorCount:]
        # remap scale factor indexes
        for f in range(1, functionCountOut + 1):
            termCount = eft.getFunctionNumberOfTerms(f)
            for t in range(1, termCount + 1):
                scalingCount = eft.getTermScaling(f, t, 0)[0]
                if scalingCount > 0:
                    oldScalingIndexes = eft.getTermScaling(f, t, scalingCount)[1]
                    if scalingCount == 1:
                        oldScalingIndexes = [oldScalingIndexes]
                    scalingIndexes = [oldToNewScaleFactorIndex[index] for index in oldScalingIndexes]
                    if nodeScaleFactorCount and (f > halfFunctionCountOut):
                        # remap top row node scale factor indexes to extra set
                        for i in range(scalingCount):
                            if scalingIndexes[i] > globalScaleFactorCount:
                                scalingIndexes[i] += nodeScaleFactorCount
                    if scalingIndexes != oldScalingIndexes:
                        eft.setTermScaling(f, t, scalingIndexes)
    elif eft is not eftStd:
        eft.setNumberOfLocalScaleFactors(0)
    if eft is not eftStd:
        assert eft.validate(), "Element " + str(elementIn.getIdentifier()) + " eft not validated"

    return eft, scalefactors


def determineTricubicHermiteEft(mesh, nodeParameters, nodeDerivativeFixedWeights=None, serendipity=False,
                                mapCrossDerivatives=False):
    """
    Determine the tricubic Hermite normal/serendipity element field template for
    interpolating node parameters at corners of a cube, by matching deltas
    between corners with node derivatives.
    Node lists use zinc ordering which varies nodes fastest over lower element coordinate.
    :param mesh: A Zinc mesh of dimension 3.
    :param nodeParameters: List over 8 local nodes in Zinc ordering of 4 parameter
    vectors x, d1, d2, d3 each with 3 components.
    :param nodeDerivativeFixedWeights: Optional list over 8 local nodes in Zinc ordering of
    list of up to 3 derivatives each containing weights for global d1, d2, d3 to fix that
    derivative, or None to use default search.
    Example: [[], [None, [0.0, -1.0, 1.0]], [], [], [], [], [], []]
    Local nodes 1, 3-8 use the standard search for all derivatives.
    Local node 2 forces d/dxi2 = -d2 + d3, d/dxi1 and d/dxi3 use the standard search.
    :param serendipity: Set to True for serendipity basis.
    :param mapCrossDerivatives: For non-serendipity Hermite basis, map cross derivatives
    as appropriate for first derivatives. For example, if d/dxi1 = d2 and d/dxi2 = d3,
    d2/dxi1dxi2 should map to -d23. If any element derivative is a sum of node derivatives,
    eliminate cross derivatives which use that direction.
    :return: eft, scale factors list [-1.0] or None. Returned eft can be further modified.
    """
    if nodeDerivativeFixedWeights is None:
        nodeDerivativeFixedWeights = [[]] * 8
    assert mesh.getDimension() == 3
    assert len(nodeParameters) == 8
    assert len(nodeParameters[0]) == 4
    assert len(nodeDerivativeFixedWeights) == 8
    assert not (serendipity and mapCrossDerivatives)
    delta12 = sub(nodeParameters[1][0], nodeParameters[0][0])
    delta34 = sub(nodeParameters[3][0], nodeParameters[2][0])
    delta56 = sub(nodeParameters[5][0], nodeParameters[4][0])
    delta78 = sub(nodeParameters[7][0], nodeParameters[6][0])
    delta13 = sub(nodeParameters[2][0], nodeParameters[0][0])
    delta24 = sub(nodeParameters[3][0], nodeParameters[1][0])
    delta57 = sub(nodeParameters[6][0], nodeParameters[4][0])
    delta68 = sub(nodeParameters[7][0], nodeParameters[5][0])
    delta15 = sub(nodeParameters[4][0], nodeParameters[0][0])
    delta26 = sub(nodeParameters[5][0], nodeParameters[1][0])
    delta37 = sub(nodeParameters[6][0], nodeParameters[2][0])
    delta48 = sub(nodeParameters[7][0], nodeParameters[3][0])
    deltas = [
        [delta12, delta13, delta15],
        [delta12, delta24, delta26],
        [delta34, delta13, delta37],
        [delta34, delta24, delta48],
        [delta56, delta57, delta15],
        [delta56, delta68, delta26],
        [delta78, delta57, delta37],
        [delta78, delta68, delta48]
    ]
    fieldmodule = mesh.getFieldmodule()
    tricubicHermiteBasis = fieldmodule.createElementbasis(
        3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY if serendipity
        else Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    eft = mesh.createElementfieldtemplate(tricubicHermiteBasis)
    scalefactors = None
    derivativeLabels = [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3]
    crossDerivativeLabels = [Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3,
                             Node.VALUE_LABEL_D2_DS2DS3, Node.VALUE_LABEL_D3_DS1DS2DS3]
    for n in range(8):
        ln = n + 1
        derivativeExpressionTerms = []
        for ed in range(3):
            if mapCrossDerivatives:
                functionNumber = n * 8 + derivativeLabels[ed]
            else:
                functionNumber = n * 4 + ed + 2
            expressionTerms = []
            if nodeDerivativeFixedWeights[n] and (ed < len(nodeDerivativeFixedWeights[n])):
                weights = nodeDerivativeFixedWeights[n][ed]
                if weights is not None:
                    for nd in range(len(weights)):
                        weight = weights[nd]
                        if weight:
                            expressionTerms.append((derivativeLabels[nd], [1] if (weight < 0.0) else []))
            if not expressionTerms:
                delta = deltas[n][ed]  # delta side coordinates for nominal derivative direction
                magDelta = magnitude(delta)
                greatestSimilarity = 0.0  # can be negative in which can
                derivativeLabel = None
                for nd in range(3):
                    derivative = nodeParameters[n][nd + 1]
                    magDerivative = magnitude(derivative)
                    cosineSimilarity = dot(derivative, delta) / (magDerivative * magDelta)
                    magnitudeSimilarity = math.exp(-math.fabs((magDerivative - magDelta) / magDelta))
                    similarity = cosineSimilarity * magnitudeSimilarity
                    if math.fabs(similarity) > math.fabs(greatestSimilarity):
                        greatestSimilarity = similarity
                        derivativeLabel = derivativeLabels[nd]
                expressionTerms.append((derivativeLabel, [1] if (greatestSimilarity < 0.0) else []))
            termCount = len(expressionTerms)
            if termCount != 1:
                eft.setFunctionNumberOfTerms(functionNumber, termCount)
            for t in range(termCount):
                term = t + 1
                eft.setTermNodeParameter(functionNumber, term, ln, expressionTerms[t][0], 1)
                scaling = expressionTerms[t][1]
                if scaling:
                    if not scalefactors:
                        # assumes only
                        setEftScaleFactorIds(eft, [1], [])
                        scalefactors = [-1.0]
                    eft.setTermScaling(functionNumber, term, scaling)
            if mapCrossDerivatives:
                derivativeExpressionTerms.append(expressionTerms)
        if mapCrossDerivatives:
            for cd in range(4):
                functionNumber = n * 8 + crossDerivativeLabels[cd]
                derivativeIndexes = \
                    [0, 1] if (cd == 0) else \
                    [0, 2] if (cd == 1) else \
                    [1, 2] if (cd == 2) else \
                    [0, 1, 2]
                sign = 1.0
                crossDerivativeLabel = Node.VALUE_LABEL_VALUE
                for ed in derivativeIndexes:
                    if len(derivativeExpressionTerms[ed]) != 1:
                        crossDerivativeLabel = Node.VALUE_LABEL_VALUE
                        break
                    crossDerivativeLabel += derivativeExpressionTerms[ed][0][0] - Node.VALUE_LABEL_VALUE
                    if derivativeExpressionTerms[ed][0][1]:
                        sign *= -1.0
                if crossDerivativeLabel == Node.VALUE_LABEL_VALUE:
                    eft.setFunctionNumberOfTerms(functionNumber, 0)
                else:
                    eft.setTermNodeParameter(functionNumber, 1, ln, crossDerivativeLabel, 1)
                    if sign < 0.0:
                        eft.setTermScaling(functionNumber, 1, [1])
    return eft, scalefactors