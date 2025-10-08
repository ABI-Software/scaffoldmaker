#%%
# Number of elements per segment. Used to calculate the number of nodes per segment. 
humanElementCounts = {
    'headElementsCount': 3, 
    'neckElementsCount': 2, 
    'brachiumElementsCount': 4, 
    'antebrachiumElementsCount': 2, 
    'handElementsCount': 1, 
    'thoraxElementsCount': 3, 
    'abdomenElementsCount': 4, 
    'upperLegElementsCount': 3,
    'lowerLegElementsCount': 2,
    'footElementsCount': 2
}


def createSegment(nodeCount:int, segmentStart:str, nodeIdentifier:int, endSegment=False, version=0):
    """
    Construct a segment of the human network node

    :param nodeCount: Number of nodes to add.
    :param networkLayout: String containing the current network layout.
    :param nodeIdentifier: Integer denoting the current node.
    :param endSegment: If true, adds a comma at the end of the segment.
    :param version: If > 0, adds version number on the last node of the segment.
    :return networklayout: String containing the updated layout.
    :return nodeIdentifier: The updated nodeIdentifier after adding the segment.
    """
    for i in range(nodeCount):
        segmentLayout = segmentStart + str(nodeIdentifier)
        if i == nodeCount - 1:
            if endSegment:
                if version == 0:
                    segmentConnector = ','
                else:
                    segmentConnector = '.' + str(version) + ','
                segmentLayout = segmentLayout + segmentConnector
            else:
                segmentLayout = segmentLayout + '-'
                nodeIdentifier += 1
        else:
            segmentLayout = segmentLayout + '-'
            nodeIdentifier += 1 
    return segmentLayout, nodeIdentifier

def constructNetworkLayoutStructure(humanElementCounts:dict):
    """
    Construct the network layout of the human wholebody scaffold. 
    The network layout consists of the following segments: 
    head, neck, thorax, abdomen, right/left arm, right/left leg. 
    Arms are subdivided into brachium, antebrachium and hand. 
    Legs are subdivided into upper leg, lower leg and foot. 

    :param humanElementCounts: Dictionary containing the number of elements 
        corresponding to each segment. 
    :return humanNetworkLayout: String containing the network layout
    """
    nodeIdentifier = 1
    # Head
    headNetworkLayout = str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    headNetworkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['headElementsCount'], 
        headNetworkLayout, nodeIdentifier, endSegment=True)
    # Neck
    necknNetworkLayout = str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    necknNetworkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['neckElementsCount'], 
        necknNetworkLayout, nodeIdentifier, endSegment=True, version=1)
    neckJointNode = nodeIdentifier
    # Thorax 
    thoraxNetworkLayout = str(nodeIdentifier) + '.1-'
    nodeIdentifier += 1 
    thoraxNetworkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['thoraxElementsCount'], 
        thoraxNetworkLayout, nodeIdentifier, endSegment=True)
    # Abdomen 
    abdomenNetworkLayout = str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    abdomenNetworkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['abdomenElementsCount'], 
        abdomenNetworkLayout, nodeIdentifier, endSegment=True, version=1)
    pelvisNodeJoint = nodeIdentifier
    # Arms
    armNetworkLayouts = []
    for i in range(2):
        version = 2 if (i == 0) else 3 #Left is 2, right is 3 
        armNetworkLayout = str(neckJointNode) + '.' + str(version) + '-'
        nodeIdentifier += 1  
        # Brachium 
        armNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['brachiumElementsCount'], 
            armNetworkLayout, nodeIdentifier)
        # Antebrachium 
        armNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['antebrachiumElementsCount'], 
            armNetworkLayout, nodeIdentifier)
        # Hand
        armNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['handElementsCount'], 
            armNetworkLayout, nodeIdentifier, endSegment=True)
        armNetworkLayouts.append(armNetworkLayout)
    #Legs 
    legNetworkLayouts = []
    for i in range(2):
        version = 2 if (i == 0) else 3 #Left is 2, right is 3 
        legNetworkLayout = str(pelvisNodeJoint) + '.' + str(version) + '-'
        nodeIdentifier += 1
        # Upper leg
        legNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['upperLegElementsCount'], 
            legNetworkLayout, nodeIdentifier)
        # Lower leg 
        legNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['lowerLegElementsCount'], 
            legNetworkLayout, nodeIdentifier)
        # Feet 
        legNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['footElementsCount'], 
            legNetworkLayout, nodeIdentifier, endSegment=True)
        legNetworkLayouts.append(legNetworkLayout)
    humanNetworkLayout = headNetworkLayout + necknNetworkLayout + \
        armNetworkLayouts[0] + armNetworkLayouts[1] + thoraxNetworkLayout + \
        abdomenNetworkLayout + legNetworkLayouts[0] + legNetworkLayouts[1]
    humanNetworkLayout = humanNetworkLayout[:-1] #Remove a comma at the end
    return humanNetworkLayout

# humanNetworkLayout.replace(',', ',\n').splitlines()

