'''
Utility functions for easing use of Zinc API.
Created on Jan 4, 2018

@author: Richard Christie
'''

from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.utils.interpolation import *
import scaffoldmaker.utils.vector as vector

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
    coordinates = fieldmodule.findFieldByName(name)
    if coordinates.isValid():
        coordinates = coordinates.castFiniteElement()
        assert coordinates.isValid(), 'getOrCreateCoordinateField.  Existing field \'' + name + '\' is not finite element type'
        assert coordinates.getNumberOfComponents() == componentsCount, 'getOrCreateCoordinateField.  Existing field \'' + name + '\' does not have ' + str(componentsCount) + ' components'
        assert coordinates.getCoordinateSystemType() == Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN, 'getOrCreateCoordinateField.  Existing field \'' + name + '\' is not rectangular Cartesian'
        return coordinates
    fieldmodule.beginChange()
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

def interpolateNodesCubicHermite(cache, coordinates, xi, normal_scale, \
        node1, derivative1, scale1, cross_derivative1, cross_scale1, \
        node2, derivative2, scale2, cross_derivative2, cross_scale2):
    """
    Interpolates position and first derivative with cubic Hermite basis.
    Interpolates cross derivative linearly.
    :param cache: Field cache to evaluate in.
    :param coordinates: Coordinates field.
    :param xi: Element coordinate to interpolate at.
    :param normal_scale: Magnitude of normal derivative to return.
    :param node1, node2: Start and end nodes.
    :param derivative1, derivative2: Node value label for derivatives.
    :param scale1, scale2: Real value scaling derivatives, to reverse if needed.
    :param cross_derivative1, cross_derivative2: Node value label for cross derivatives.
    :param cross_scale1, cross_scale2: Real value scaling cross_derivatives, to reverse if needed.
    :return: x, dx_ds, dx_ds_cross, dx_ds_normal
    """
    cache.setNode(node1)
    result, v1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
    result, d1 = coordinates.getNodeParameters(cache, -1, derivative1, 1, 3 )
    result, d1c = coordinates.getNodeParameters(cache, -1, cross_derivative1, 1, 3 )
    d1 = [ scale1*d for d in d1 ]
    d1c = [ cross_scale1*d for d in d1c ]
    cache.setNode(node2)
    result, v2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
    result, d2 = coordinates.getNodeParameters(cache, -1, derivative2, 1, 3 )
    result, d2c = coordinates.getNodeParameters(cache, -1, cross_derivative2, 1, 3 )
    d2 = [ scale2*d for d in d2 ]
    d2c = [ cross_scale2*d for d in d2c ]

    arcLength = computeCubicHermiteArcLength(v1, d1, v2, d2, True)
    mag = arcLength/vector.magnitude(d1)
    d1 = [ mag*d for d in d1 ]
    mag = arcLength/vector.magnitude(d2)
    d2 = [ mag*d for d in d2 ]

    xr = 1.0 - xi
    x = list(interpolateCubicHermite(v1, d1, v2, d2, xi))
    dx_ds = interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi)
    scale = min(xi, xr)
    dx_ds = [ scale*d for d in dx_ds ]
    dx_ds_cross = [ (xr*d1c[c] + xi*d2c[c]) for c in range(3) ]

    radialVector = vector.normalise(vector.crossproduct3(dx_ds_cross, dx_ds))
    dx_ds_normal = [ normal_scale*d for d in radialVector ]

    return x, dx_ds, dx_ds_cross, dx_ds_normal

def computeNodeDerivativeHermiteLagrange(cache, coordinates, node1, derivative1, scale1, node2, scale2):
    """
    Computes the derivative at node2 from quadratic Hermite-Lagrange interpolation of
    node1 value and derivative1 to node2 value.
    :param cache: Field cache to evaluate in.
    :param coordinates: Coordinates field.
    :param node1, node2: Start and end nodes.
    :param derivative1: Node value label for derivative at node1.
    :param scale1, scale2: Scaling to apply to derivatives at nodes, e.g. -1.0 to reverse.
    :return: dx_dxi at node2
    """
    cache.setNode(node1)
    result, v1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
    result, d1 = coordinates.getNodeParameters(cache, -1, derivative1, 1, 3 )
    d1 = [ d*scale1 for d in d1 ]
    cache.setNode(node2)
    result, v2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3 )
    xi = 1.0
    #phi1 = 1 - xi*xi
    #phi2 = xi - xi*xi
    #phi3 = xi*xi
    dphi1 = -2.0*xi
    dphi2 = 1 - 2.0*xi
    dphi3 = 2.0*xi
    d2 = [ (v1[c]*dphi1 + d1[c]*dphi2 + v2[c]*dphi3) for c in range(3) ]
    d2 = [ d*scale2 for d in d2 ]
    return d2
