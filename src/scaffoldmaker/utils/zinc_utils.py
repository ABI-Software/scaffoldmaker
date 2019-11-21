'''
Utility functions for easing use of Zinc API.
'''

from opencmiss.zinc.context import Context
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node, Nodeset
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector

class ZincCacheChanges:
    """
    Context manager for ensuring beginChange, endChange always called on
    supplied object, even with exceptions.
    Usage:
    with ZincCacheChanges(object):
        # make multiple changes to object or objects it owns
    """

    def __init__(self, object):
        """
        :param object: Zinc object with beginChange/endChange methods.
        """
        self._object = object

    def __enter__(self):
        self._object.beginChange()
        return self

    def __exit__(self, *args):
        self._object.endChange()


def getOrCreateCoordinateField(fieldmodule, name='coordinates', componentsCount=3):
    '''
    Finds or creates a rectangular cartesian coordinate field.
    New field has component names: 'x', 'y', 'z'.
    Raises exception if existing field of name is not finite element type or has incorrect attributes.
    :param fieldmodule:  Zinc fieldmodule to find or create field in.
    :param name:  Name of field to find or create.
    :param componentsCount: Number of components / dimension of field, from 1 to 3.
    '''
    assert (componentsCount > 0) and (componentsCount <= 3), 'getOrCreateCoordinateField.  Dimensions must be from 1 to 3'
    coordinates = fieldmodule.findFieldByName(name)
    if coordinates.isValid():
        coordinates = coordinates.castFiniteElement()
        assert coordinates.isValid(), 'getOrCreateCoordinateField.  Existing field \'' + name + '\' is not finite element type'
        assert coordinates.getNumberOfComponents() == componentsCount, 'getOrCreateCoordinateField.  Existing field \'' + name + '\' does not have ' + str(componentsCount) + ' components'
        assert coordinates.getCoordinateSystemType() == Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN, 'getOrCreateCoordinateField.  Existing field \'' + name + '\' is not rectangular Cartesian'
        return coordinates
    with ZincCacheChanges(fieldmodule):
        coordinates = fieldmodule.createFieldFiniteElement(componentsCount)
        coordinates.setName(name)
        coordinates.setManaged(True)
        coordinates.setTypeCoordinate(True)
        coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
        for c in range(componentsCount):
            coordinates.setComponentName(c + 1, ['x', 'y', 'z'][c])
    return coordinates

def getOrCreateFibreField(fieldmodule, name='fibres', componentsCount=3):
    '''
    Finds or creates a fibre field.
    New field has component names: 'fibre angle', 'imbrication angle', 'sheet angle'.
    Raises exception if existing field of name is not finite element type or has incorrect attributes.
    :param fieldmodule:  Zinc fieldmodule to find or create field in.
    :param name:  Name of field to find or create.
    :param componentsCount: Number of components of field, from 1 to 3.
    '''
    assert (componentsCount > 0) and (componentsCount <= 3), 'getOrCreateFibreField.  Dimensions must be from 1 to 3'
    fibres = fieldmodule.findFieldByName(name)
    if fibres.isValid():
        fibres = fibres.castFiniteElement()
        assert fibres.isValid(), 'getOrCreateFibreField.  Existing field \'' + name + '\' is not finite element type'
        assert fibres.getNumberOfComponents() == componentsCount, 'getOrCreateFibreField.  Existing field \'' + name + '\' does not have ' + str(componentsCount) + ' components'
        assert fibres.getCoordinateSystemType() == Field.COORDINATE_SYSTEM_TYPE_FIBRE, 'getOrCreateFibreField.  Existing field \'' + name + '\' is not fibre'
        return fibres
    with ZincCacheChanges(fieldmodule):
        fibres = fieldmodule.createFieldFiniteElement(componentsCount)
        fibres.setName(name)
        fibres.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_FIBRE)
        for c in range(componentsCount):
            fibres.setComponentName(c + 1, ['fibre angle', 'imbrication angle', 'sheet angle'][c])
    return fibres

def getOrCreateTextureCoordinateField(fieldmodule, name='texture coordinates', componentsCount=3):
    '''
    Finds or creates a rectangular cartesian texture coordinate field.
    New field has component names: 'u', 'v', 'w'.
    Raises exception if existing field of name is not finite element type or has incorrect attributes.
    :param fieldmodule:  Zinc fieldmodule to find or create field in.
    :param name:  Name of field to find or create.
    :param componentsCount: Number of components / dimension of field, from 1 to 3.
    '''
    assert (componentsCount > 0) and (componentsCount <= 3), 'getOrCreateTextureCoordinateField.  Dimensions must be from 1 to 3'
    coordinates = fieldmodule.findFieldByName(name)
    if coordinates.isValid():
        coordinates = coordinates.castFiniteElement()
        assert coordinates.isValid(), 'getOrCreateTextureCoordinateField.  Existing field \'' + name + '\' is not finite element type'
        assert coordinates.getNumberOfComponents() == componentsCount, 'getOrCreateTextureCoordinateField.  Existing field \'' + name + '\' does not have ' + str(componentsCount) + ' components'
        assert coordinates.getCoordinateSystemType() == Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN, 'getOrCreateTextureCoordinateField.  Existing field \'' + name + '\' is not rectangular Cartesian'
        return coordinates
    with ZincCacheChanges(fieldmodule):
        coordinates = fieldmodule.createFieldFiniteElement(componentsCount)
        coordinates.setName(name)
        coordinates.setManaged(True)
        coordinates.setTypeCoordinate(True)
        coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
        for c in range(componentsCount):
            coordinates.setComponentName(c + 1, ['u', 'v', 'w'][c])
    return coordinates

def getOrCreateFlatCoordinateField(fieldmodule, name='flat coordinates', componentsCount=3):
    '''
    Finds or creates a rectangular cartesian texture coordinate field.
    New field has component names: 'x', 'y', 'z'.
    Raises exception if existing field of name is not finite element type or has incorrect attributes.
    :param fieldmodule:  Zinc fieldmodule to find or create field in.
    :param name:  Name of field to find or create.
    :param componentsCount: Number of components / dimension of field, from 1 to 3.
    '''
    assert (componentsCount > 0) and (componentsCount <= 3), 'getOrCreateFlatCoordinateField.  Dimensions must be from 1 to 3'
    coordinates = fieldmodule.findFieldByName(name)
    if coordinates.isValid():
        coordinates = coordinates.castFiniteElement()
        assert coordinates.isValid(), 'getOrCreateFlatCoordinateField.  Existing field \'' + name + '\' is not finite element type'
        assert coordinates.getNumberOfComponents() == componentsCount, 'getOrCreateFlatCoordinateField.  Existing field \'' + name + '\' does not have ' + str(componentsCount) + ' components'
        assert coordinates.getCoordinateSystemType() == Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN, 'getOrCreateFlatCoordinateField.  Existing field \'' + name + '\' is not rectangular Cartesian'
        return coordinates
    with ZincCacheChanges(fieldmodule):
        coordinates = fieldmodule.createFieldFiniteElement(componentsCount)
        coordinates.setName(name)
        coordinates.setManaged(True)
        coordinates.setTypeCoordinate(True)
        coordinates.setCoordinateSystemType(Field.COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN)
        for c in range(componentsCount):
            coordinates.setComponentName(c + 1, ['x', 'y', 'z'][c])
    return coordinates

def getOrCreateElementXiField(fieldmodule, name='element_xi', mesh=None):
    '''
    Finds or creates a stored mesh location field for storing locations in the
    supplied mesh e.g. for defining on annotation points with mesh locations.
    Raises exception if existing field of name is not stored mesh location type.
    Note can't currently verify existing field stores locations in the supplied mesh.
    :param fieldmodule:  Zinc fieldmodule to find or create field in.
    :param name:  Name of field to find or create.
    :param mesh:  Mesh to store locations in.
    '''
    if mesh is None:
        mesh = fieldmodule.findMeshByDimension(3)
    assert mesh.isValid(), 'getOrCreateElementXiField.  Invalid mesh'
    elementXiField = fieldmodule.findFieldByName(name)
    if elementXiField.isValid():
        elementXiField = elementXiField.castStoredMeshLocation()
        assert elementXiField.isValid(), 'getOrCreateElementXiField.  Existing field \'' + name + '\' is not stored mesh location type'
        return elementXiField
    with ZincCacheChanges(fieldmodule):
        elementXiField = fieldmodule.createFieldStoredMeshLocation(mesh)
        elementXiField.setName(name)
        elementXiField.setManaged(True)
    return elementXiField

def getOrCreateLabelField(fieldmodule, name='label'):
    '''
    Finds or creates a stored string field for defining labels on nodes, e.g. annotation points.
    Raises exception if existing field of name is not string-valued.
    Note can't currently distinguish stored string from constant string fields.
    :param fieldmodule:  Zinc fieldmodule to find or create field in.
    :param name:  Name of field to find or create.
    '''
    labelField = fieldmodule.findFieldByName(name)
    if labelField.isValid():
        assert labelField.getValueType() == Field.VALUE_TYPE_STRING, 'getOrCreateLabelField.  Existing field \'' + name + '\' is not string valued'
        return labelField
    with ZincCacheChanges(fieldmodule):
        labelField = fieldmodule.createFieldStoredString()
        labelField.setName(name)
        labelField.setManaged(True)
    return labelField

def getOrCreateGroupField(fieldmodule, name):
    '''
    Finds or creates a Group field of the supplied name.
    Raises exception if existing field of name is not a group.
    :param fieldmodule:  Zinc fieldmodule to find or create field in.
    :param name:  Name of field to find or create.
    '''
    group = fieldmodule.findFieldByName(name)
    if group.isValid():
        group = group.castGroup()
        assert group.isValid(), 'getOrCreateGroupField.  Existing field \'' + name + '\' is not a group type'
        return group
    with ZincCacheChanges(fieldmodule):
        group = fieldmodule.createFieldGroup()
        group.setName(name)
        group.setManaged(True)
    return group

def getOrCreateNodesetGroup(group, nodeset):
    '''
    Gets or creates the NodesetGroup for the supplied nodeset in group.
    :param group:  Zinc FieldGroup.
    :param nodeset:  A nodeset from group region to get or create subgroup of.
    '''
    nodeGroup = group.getFieldNodeGroup(nodeset)
    if not nodeGroup.isValid():
        nodeGroup = group.createFieldNodeGroup(nodeset)
    return nodeGroup.getNodesetGroup()

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

def getElementNodeIdentifiers4Node(element, eft):
    '''
    Get 4 node identifiers for an element with 4 basis nodes, handling
    collapses e.g. where eft has fewer nodes. Asserts basis has 4 nodes.
    :param element: Element to query.
    :param eft: Element field template nodes are stored for in element.
    :return: List of 4 local node identifiers.
    '''
    elementbasis = eft.getElementbasis()
    basisNodesCount = elementbasis.getNumberOfNodes()
    assert basisNodesCount == 4, 'getElementNodeIdentifiers4Node:  Element ' + str(element.getIdentifier()) + ' is not using a 4 node basis'
    nodeIdentifiers = []
    fn = 1
    for n in range(basisNodesCount):
        ln = eft.getTermLocalNodeIndex(fn, 1)
        nodeIdentifiers.append(element.getNode(eft, ln).getIdentifier())
        fn += elementbasis.getNumberOfFunctionsPerNode(n + 1)
    return nodeIdentifiers

def getElementNodeIdentifiers8Node(element, eft):
    '''
    Get 8 node identifiers for an element with 8 basis nodes, handling
    collapses e.g. where eft has fewer nodes. Asserts basis has 8 nodes.
    :param element: Element to query.
    :param eft: Element field template nodes are stored for in element.
    :return: List of 8 local node identifiers.
    '''
    elementbasis = eft.getElementbasis()
    basisNodesCount = elementbasis.getNumberOfNodes()
    assert basisNodesCount == 8, 'getElementNodeIdentifiers8Node:  Element ' + str(element.getIdentifier()) + ' is not using an 8 node basis'
    nodeIdentifiers = []
    fn = 1
    for n in range(basisNodesCount):
        ln = eft.getTermLocalNodeIndex(fn, 1)
        nodeIdentifiers.append(element.getNode(eft, ln).getIdentifier())
        fn += elementbasis.getNumberOfFunctionsPerNode(n + 1)
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

def evaluateFieldRange(field : Field, nodeset : Nodeset):
    """
    :return: minimums, maximums (as lists for components).
    """
    fieldmodule = field.getFieldmodule()
    componentsCount = field.getNumberOfComponents()
    with ZincCacheChanges(fieldmodule):
        fieldNodesetMinimum = fieldmodule.createFieldNodesetMinimum(field, nodeset)
        fieldNodesetMaximum = fieldmodule.createFieldNodesetMaximum(field, nodeset)
        fieldcache = fieldmodule.createFieldcache()
        result, minimums = fieldNodesetMinimum.evaluateReal(fieldcache, componentsCount)
        assert result == RESULT_OK
        result, maximums = fieldNodesetMaximum.evaluateReal(fieldcache, componentsCount)
        assert result == RESULT_OK
        del fieldNodesetMinimum
        del fieldNodesetMaximum
    return minimums, maximums

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

    arcLength = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
    mag = arcLength/vector.magnitude(d1)
    d1 = [ mag*d for d in d1 ]
    mag = arcLength/vector.magnitude(d2)
    d2 = [ mag*d for d in d2 ]

    xr = 1.0 - xi
    x = interp.interpolateCubicHermite(v1, d1, v2, d2, xi)
    dx_ds = interp.interpolateCubicHermiteDerivative(v1, d1, v2, d2, xi)
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
    d2 = interp.interpolateHermiteLagrangeDerivative(v1, d1, v2, 1.0)
    d2 = [ d*scale2 for d in d2 ]
    return d2

def exnodeStringFromNodeValues(
        nodeValueLabels = [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1 ],
        nodeValues = [
            [ [ 0.0, 0.0, 0.0 ], [ 1.0, 0.0, 0.0 ] ],
            [ [ 1.0, 0.0, 0.0 ], [ 1.0, 0.0, 0.0 ] ] ],
        groupName = 'meshEdits'):
    '''
    Return a string in Zinc EX format defining nodes 1..N with the supplied
    coordinate values and their labels. Works in a private zinc context.
    '''
    # following requires at least one value label and node, assumes consistent values and components counts
    nodeValueLabelsCount = len(nodeValueLabels)
    nodesCount = len(nodeValues)
    componentsCount = len(nodeValues[0][0])
    context = Context('exnodeStringFromNodeValues')
    region = context.getDefaultRegion()
    fieldmodule = region.getFieldmodule()
    with ZincCacheChanges(fieldmodule):
        cache = fieldmodule.createFieldcache()
        coordinates = getOrCreateCoordinateField(fieldmodule, componentsCount = componentsCount)
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        group = fieldmodule.createFieldGroup()
        group.setName(groupName)
        nodesetGroup = group.createFieldNodeGroup(nodes).getNodesetGroup()
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        if not Node.VALUE_LABEL_VALUE in nodeValueLabels:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 0)
        for nodeValueLabel in nodeValueLabels:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, nodeValueLabel, 1)
        # create nodes
        for n in range(nodesCount):
            node = nodesetGroup.createNode(n + 1, nodetemplate)
            cache.setNode(node)
            for v in range(nodeValueLabelsCount):
                coordinates.setNodeParameters(cache, -1, nodeValueLabels[v], 1, nodeValues[n][v])
        # serialise to string
        sir = region.createStreaminformationRegion()
        srm = sir.createStreamresourceMemory()
        sir.setResourceGroupName(srm, groupName)
        region.write(sir)
        result, exString = srm.getBuffer()
    return exString
