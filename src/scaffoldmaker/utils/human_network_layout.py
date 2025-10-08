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


def addNodesToLayout(nodeCount:int, networkLayout:str, nodeIdentifier:int, endSegment=False, version=0):
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
        networkLayout = networkLayout + str(nodeIdentifier)
        if i == nodeCount - 1:
            if endSegment:
                if version == 0:
                    segmentConnector = ','
                else:
                    segmentConnector = '.' + str(version) + ','
                networkLayout = networkLayout + segmentConnector
            else:
                networkLayout = networkLayout + '-'
                nodeIdentifier += 1
        else:
            networkLayout = networkLayout + '-'
            nodeIdentifier += 1 
    return networkLayout, nodeIdentifier

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
    humanNetWorkLayout = ""
    nodeIdentifier = 1
    # Head
    humanNetWorkLayout = humanNetWorkLayout + str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['headElementsCount'], 
        humanNetWorkLayout, nodeIdentifier, endSegment=True)
    # Neck
    humanNetWorkLayout = humanNetWorkLayout + str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['neckElementsCount'], 
        humanNetWorkLayout, nodeIdentifier, endSegment=True, version=1)
    neckJointNode = nodeIdentifier
    # Thorax 
    humanNetWorkLayout = humanNetWorkLayout + str(nodeIdentifier) + '.1-'
    nodeIdentifier += 1 
    humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['thoraxElementsCount'], 
        humanNetWorkLayout, nodeIdentifier, endSegment=True)
    # Abdomen 
    humanNetWorkLayout = humanNetWorkLayout + str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['abdomenElementsCount'], 
        humanNetWorkLayout, nodeIdentifier, endSegment=True, version=1)
    pelvisNodeJoint = nodeIdentifier
    # Arms
    for i in range(2):
        version = 2 if (i == 0) else 3 #Left is 2, right is 3 
        humanNetWorkLayout = humanNetWorkLayout + str(neckJointNode) + '.' + str(version) + '-'
        nodeIdentifier += 1  
        # Brachium 
        humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['brachiumElementsCount'], 
            humanNetWorkLayout, nodeIdentifier)
        # Antebrachium 
        humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['antebrachiumElementsCount'], 
            humanNetWorkLayout, nodeIdentifier)
        # Hand
        humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['handElementsCount'], 
            humanNetWorkLayout, nodeIdentifier, endSegment=True)
    #Legs 
    for i in range(2):
        version = 2 if (i == 0) else 3 #Left is 2, right is 3 
        humanNetWorkLayout = humanNetWorkLayout + str(pelvisNodeJoint) + '.' + str(version) + '-'
        nodeIdentifier += 1
        # Upper leg
        humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['upperLegElementsCount'], 
            humanNetWorkLayout, nodeIdentifier)
        # Lower leg 
        humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['lowerLegElementsCount'], 
            humanNetWorkLayout, nodeIdentifier)
        # Feet 
        humanNetWorkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['footElementsCount'], 
            humanNetWorkLayout, nodeIdentifier, endSegment=True)
    return humanNetWorkLayout

# humanNetWorkLayout.replace(',', ',\n').splitlines()

