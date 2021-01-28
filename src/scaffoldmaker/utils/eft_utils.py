'''
Utility functions for element field templates shared by mesh generators.
'''
from opencmiss.zinc.element import Elementfieldtemplate

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
            if (localNodeIndex in localNodeIndexes) and (eft.getTermNodeValueLabel(f, 1) == fromValueLabel) and (not getEftTermScaling(eft, f, 1)):
                termCount = len(expressionTerms)
                eft.setFunctionNumberOfTerms(f, termCount)
                version = eft.getTermNodeVersion(f, 1)
                for t in range(1, termCount + 1):
                    expressionTerm = expressionTerms[t - 1]
                    eft.setTermNodeParameter(f, t, localNodeIndex, expressionTerm[0], version)
                    if expressionTerm[1]:
                        eft.setTermScaling(f, t, expressionTerm[1])

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
