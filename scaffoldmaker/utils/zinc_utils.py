'''
Utility functions for easing use of Zinc API.
Created on Jan 4, 2018

@author: Richard Christie
'''

def getElementNodeIdentifiers(element, eft):
    '''
    Get identifiers of all nodes used by eft in element.
    '''
    nodeIdentifiers = []
    nodeCount = eft.getNumberOfLocalNodes()
    for n in range(nodeCount):
        node = element.getNode(eft, n + 1)
        nodeIdentifiers.append(node.getIdentifier())
    return nodeIdentifiers

def getMaximumNodeIdentifier(nodeset):
    """
    :return: Maximum node identifier in nodeset or -1 if none.
    """
    maximumNodeId = -1
    nodeiterator = nodeset.createNodeiterator()
    node = nodeiterator.next()
    while node.isValid():
        id = node.getIdentifier()
        if id > maximumNodeId:
            maximumNodeId = id
        node = nodeiterator.next()
    return maximumNodeId

def getMaximumElementIdentifier(mesh):
    """
    :return: Maximum element identifier in mesh or -1 if none.
    """
    maximumElementId = -1
    elementiterator = mesh.createElementiterator()
    element = elementiterator.next()
    while element.isValid():
        id = element.getIdentifier()
        if id > maximumElementId:
            maximumElementId = id
        element = elementiterator.next()
    return maximumElementId
   