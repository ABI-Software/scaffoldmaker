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
