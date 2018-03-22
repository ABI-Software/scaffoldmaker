'''
Utility functions for easing use of Zinc API.
Created on Jan 4, 2018

@author: Richard Christie
'''

from opencmiss.zinc.field import Field

def getOrCreateCoordinateField(fieldmodule, name='coordinates', componentsCount=3):
    '''
    Finds or creates a rectangular cartesian coordinate field.
    New field has component names, 'x', 'y', 'z'.
    Raises exception if existing field of name is not finite element type or has incorrect attributes.
    :param fieldmodule:  Zinc fieldmodule to find or create field in.
    :param name:  Name of field to find or create.
    :param componentsCount: Number of components / dimension of field.
    '''
    assert (componentsCount > 0) and (componentsCount <= 3), 'getOrCreateCoordinateField.  Dimensions must be from 1 to 3'
    fieldmodule.beginChange()
    coordinates = fieldmodule.findFieldByName(name)
    if coordinates.isValid():
        coordinates = coordinates.castFiniteElement()
        assert coordinates.isValid(), 'getOrCreateCoordinateField.  Existing field \'' + name + '\' is not finite element type'
        assert coordinates.getNumberOfComponents() == componentsCount, 'getOrCreateCoordinateField.  Existing field \'' + name + '\' does not have ' + str(componentsCount) + ' components'
        assert coordinates.getCoordinateSystemType() == Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN, 'getOrCreateCoordinateField.  Existing field \'' + name + '\' is not rectangular Cartesian'
        return coordinates
    coordinates = fieldmodule.createFieldFiniteElement(componentsCount)
    coordinates.setName(name)
    coordinates.setManaged(True)
    coordinates.setTypeCoordinate(True)
    coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
    for c in range(componentsCount):
        coordinates.setComponentName(c + 1, ['x', 'y', 'z'][c])
    fieldmodule.endChange()
    return coordinates

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
   