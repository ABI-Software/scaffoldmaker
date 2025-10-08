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
    humanNetworkLayout = ""
    nodeIdentifier = 1
    # Head
    humanNetworkLayout = humanNetworkLayout + str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    humanNetworkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['headElementsCount'], 
        humanNetworkLayout, nodeIdentifier, endSegment=True)
    # Neck
    humanNetworkLayout = humanNetworkLayout + str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    humanNetworkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['neckElementsCount'], 
        humanNetworkLayout, nodeIdentifier, endSegment=True, version=1)
    neckJointNode = nodeIdentifier
    # Thorax 
    humanNetworkLayout = humanNetworkLayout + str(nodeIdentifier) + '.1-'
    nodeIdentifier += 1 
    humanNetworkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['thoraxElementsCount'], 
        humanNetworkLayout, nodeIdentifier, endSegment=True)
    # Abdomen 
    humanNetworkLayout = humanNetworkLayout + str(nodeIdentifier) + '-'
    nodeIdentifier += 1 
    humanNetworkLayout, nodeIdentifier = addNodesToLayout(
        humanElementCounts['abdomenElementsCount'], 
        humanNetworkLayout, nodeIdentifier, endSegment=True, version=1)
    pelvisNodeJoint = nodeIdentifier
    # Arms
    for i in range(2):
        version = 2 if (i == 0) else 3 #Left is 2, right is 3 
        humanNetworkLayout = humanNetworkLayout + str(neckJointNode) + '.' + str(version) + '-'
        nodeIdentifier += 1  
        # Brachium 
        humanNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['brachiumElementsCount'], 
            humanNetworkLayout, nodeIdentifier)
        # Antebrachium 
        humanNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['antebrachiumElementsCount'], 
            humanNetworkLayout, nodeIdentifier)
        # Hand
        humanNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['handElementsCount'], 
            humanNetworkLayout, nodeIdentifier, endSegment=True)
    #Legs 
    for i in range(2):
        version = 2 if (i == 0) else 3 #Left is 2, right is 3 
        humanNetworkLayout = humanNetworkLayout + str(pelvisNodeJoint) + '.' + str(version) + '-'
        nodeIdentifier += 1
        # Upper leg
        humanNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['upperLegElementsCount'], 
            humanNetworkLayout, nodeIdentifier)
        # Lower leg 
        humanNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['lowerLegElementsCount'], 
            humanNetworkLayout, nodeIdentifier)
        # Feet 
        humanNetworkLayout, nodeIdentifier = addNodesToLayout(
            humanElementCounts['footElementsCount'], 
            humanNetworkLayout, nodeIdentifier, endSegment=True)
    humanNetworkLayout = humanNetworkLayout[:-1] #Remove a comma at the end
    return humanNetworkLayout

# humanNetWorkLayout.replace(',', ',\n').splitlines()

