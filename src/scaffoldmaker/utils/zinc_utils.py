'''
Utility functions for easing use of Zinc API.
'''

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import Mesh, MeshGroup
from opencmiss.zinc.field import Field, FieldGroup
from opencmiss.zinc.fieldmodule import Fieldmodule
from opencmiss.zinc.node import Node, Nodeset
from opencmiss.zinc.result import RESULT_OK
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import vector


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
    with ChangeManager(fieldmodule):
        cache = fieldmodule.createFieldcache()
        coordinates = findOrCreateFieldCoordinates(fieldmodule, components_count = componentsCount)
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


def createFaceMeshGroupExteriorOnFace(fieldmodule : Fieldmodule, elementFaceType) -> MeshGroup:
    """
    Returns mesh group for the exterior surface on the face described
    by elementFaceType.
    """
    with ChangeManager(fieldmodule):
        isExterior = fieldmodule.createFieldIsExterior()
        isOnFace = fieldmodule.createFieldIsOnFace(elementFaceType)
        mesh2d = fieldmodule.findMeshByDimension(2)
        faceElementGroup = fieldmodule.createFieldElementGroup(mesh2d)
        faceMeshGroup = faceElementGroup.getMeshGroup()
        faceMeshGroup.addElementsConditional(fieldmodule.createFieldAnd(isExterior, isOnFace))
        del isExterior
        del isOnFace
    return faceMeshGroup


def get_highest_dimension_mesh(fieldmodule : Fieldmodule) -> Mesh:
    '''
    Get highest dimension non-empty mesh.
    :return: Zinc Mesh or None if all are empty.
    '''
    for dimension in range(3, 0, -1):
        mesh = fieldmodule.findMeshByDimension(dimension)
        if mesh.getSize() > 0:
            return mesh
    return None


def get_next_unused_node_identifier(nodeset: Nodeset, start_identifier=1) -> int:
    """
    :return: Unused node identifier >= start_identifier.
    """
    identifier = start_identifier
    node = nodeset.findNodeByIdentifier(identifier)
    while node.isValid():
        identifier += 1
        node = nodeset.findNodeByIdentifier(identifier)
    return identifier


def group_add_group_elements(group : FieldGroup, other_group : FieldGroup, only_dimension=None):
    '''
    Add to group elements and/or nodes from other_group.
    :param only_dimension: If set, only add objects of this dimension.
    '''
    fieldmodule = group.getFieldmodule()
    with ChangeManager(fieldmodule):
        for dimension in [ only_dimension ] if only_dimension else range(4):
            if dimension > 0:
                mesh = fieldmodule.findMeshByDimension(dimension)
                element_group = group.getFieldElementGroup(mesh)
                if not element_group.isValid():
                    element_group = group.createFieldElementGroup(mesh)
                mesh_group = element_group.getMeshGroup()
                mesh_group.addElementsConditional(other_group.getFieldElementGroup(mesh))
            elif dimension == 0:
                nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
                node_group = group.getFieldNodeGroup(nodeset)
                if not node_group.isValid():
                    node_group = group.createFieldNodeGroup(nodeset)
                nodeset_group = node_group.getNodesetGroup()
                nodeset_group.addNodesConditional(other_group.getFieldNodeGroup(nodeset))


def group_get_highest_dimension(group : FieldGroup):
    '''
    Get highest dimension of elements or nodes in group.
    :return: Dimensions from 3-0, or -1 if empty.
    '''
    fieldmodule = group.getFieldmodule()
    for dimension in range(3, 0, -1):
        mesh = fieldmodule.findMeshByDimension(dimension)
        element_group = group.getFieldElementGroup(mesh)
        if element_group.isValid() and (element_group.getMeshGroup().getSize() > 0):
            return dimension
    nodeset = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    node_group = group.getFieldNodeGroup(nodeset)
    if node_group.isValid() and (node_group.getNodesetGroup().getSize() > 0):
        return 0
    return -1


def identifier_ranges_fix(identifier_ranges):
    '''
    Sort from lowest to highest identifier and merge adjacent and overlapping
    ranges.
    :param identifier_ranges: List of identifier ranges. Modified in situ.
    '''
    identifier_ranges.sort()
    i = 1
    while i < len(identifier_ranges):
        if identifier_ranges[i][0] <= (identifier_ranges[i - 1][1] + 1):
            if identifier_ranges[i][1] > identifier_ranges[i - 1][1]:
                identifier_ranges[i - 1][1] = identifier_ranges[i][1]
            identifier_ranges.pop(i)
        else:
            i += 1


def identifier_ranges_from_string(identifier_ranges_string):
    '''
    Parse string containing identifiers and identifier ranges.
    Function is suitable for processing manual input with whitespace, trailing non-digits.
    Ranges are sorted so strictly increasing. Overlapping ranges are merged.
    Future: migrate to use .. as separator for compatibility with EX file groups and cmgui.
    :param identifier_ranges_string: Identifier ranges as a string e.g. '1-30,55,66-70'.
    '30-1, 55,66-70s' also produces the same result.
    :return: Ordered list of identifier ranges e.g. [[1,30],[55,55],[66,70]]
    '''
    identifier_ranges = []
    for identifier_range_string in identifier_ranges_string.split(','):
        try:
            identifier_range_ends = identifier_range_string.split('-')
            # after leading whitespace, stop at first non-digit
            for e in range(len(identifier_range_ends)):
                # strip whitespace, trailing non digits
                digits = identifier_range_ends[e].strip()
                for i in range(len(digits)):
                    if not digits[i].isdigit():
                        digits = digits[:i]
                        break;
                identifier_range_ends[e] = digits
            start = int(identifier_range_ends[0])
            if len(identifier_range_ends) == 1:
                stop = start
            else:
                stop = int(identifier_range_ends[1])
                # ensure range is low-high
                if stop < start:
                    start, stop = stop, start
            identifier_ranges.append([start, stop])
        except:
            pass
    identifier_ranges_fix(identifier_ranges)
    return identifier_ranges


def identifier_ranges_to_string(identifier_ranges):
    '''
    Convert ranges to a string, contracting single object ranges.
    Future: migrate to use .. as separator for compatibility with EX file groups and cmgui.
    :param identifier_ranges: Ordered list of identifier ranges e.g. [[1,30],[55,55],[66,70]]
    :return: Identifier ranges as a string e.g. '1-30,55,66-70'
    '''
    identifier_ranges_string = ''
    first = True
    for identifier_range in identifier_ranges:
        if identifier_range[0] == identifier_range[1]:
            identifier_range_string = str(identifier_range[0])
        else:
            identifier_range_string = str(identifier_range[0]) + '-' + str(identifier_range[1])
        if first:
            identifier_ranges_string = identifier_range_string
            first = False
        else:
            identifier_ranges_string += ',' + identifier_range_string
    return identifier_ranges_string


def domain_iterator_to_identifier_ranges(iterator):
    '''
    Extract sorted identifier ranges from iterator.
    Currently requires iterator to be in lowest-highest identifier order.
    Objects must support getIdentifier() method returning unique integer.
    :param iterator: A Zinc Elementiterator or Nodeiterator.
    :return: List of sorted identifier ranges [start,stop] e.g. [[1,30],[55,55],[66,70]]
    '''
    identifier_ranges = []
    obj = iterator.next()
    if obj.isValid():
        stop = start = obj.getIdentifier()
        obj = iterator.next()
        while obj.isValid():
            identifier = obj.getIdentifier()
            if identifier == (stop + 1):
                stop = identifier
            else:
                identifier_ranges.append([ start, stop ])
                stop = start = identifier
            obj = iterator.next()
        identifier_ranges.append([ start, stop ])
    return identifier_ranges


def mesh_group_add_identifier_ranges(mesh_group, identifier_ranges):
    '''
    Add elements with the supplied identifier ranges to mesh_group.
    :param mesh_group: Zinc MeshGroup to modify.
    '''
    mesh = mesh_group.getMasterMesh()
    fieldmodule = mesh.getFieldmodule()
    with ChangeManager(fieldmodule):
        for identifier_range in identifier_ranges:
            for identifier in range(identifier_range[0], identifier_range[1] + 1):
                element = mesh.findElementByIdentifier(identifier)
                mesh_group.addElement(element)


def mesh_group_to_identifier_ranges(mesh_group):
    '''
    :param mesh_group: Zinc MeshGroup.
    :return: Ordered list of element identifier ranges e.g. [[1,30],[55,55],[66,70]]
    '''
    return domain_iterator_to_identifier_ranges(mesh_group.createElementiterator())


def nodeset_group_add_identifier_ranges(nodeset_group, identifier_ranges):
    '''
    Add nodes with the supplied identifier ranges to nodeset_group.
    :param nodeset_group: Zinc NodesetGroup to modify.
    '''
    nodeset = nodeset_group.getMasterNodeset()
    fieldmodule = nodeset.getFieldmodule()
    with ChangeManager(fieldmodule):
        for identifier_range in identifier_ranges:
            for identifier in range(identifier_range[0], identifier_range[1] + 1):
                node = nodeset.findNodeByIdentifier(identifier)
                nodeset_group.addNode(node)


def nodeset_group_to_identifier_ranges(nodeset_group):
    '''
    :param nodeset_group: Zinc NodesetGroup.
    :return: Ordered list of node identifier ranges e.g. [[1,30],[55,55],[66,70]]
    '''
    return domain_iterator_to_identifier_ranges(nodeset_group.createNodeiterator())


def mesh_destroy_elements_and_nodes_by_identifiers(mesh, element_identifiers):
    '''
    Deletes elements and related nodes using element identifiers.
    :param element_identifiers: Element identifiers for elements to be deleted.
    '''
    fm = mesh.getFieldmodule()
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    with ChangeManager(fm):
        # put the elements in a group and use subelement handling to get nodes in use by it
        destroyGroup = fm.createFieldGroup()
        destroyGroup.setSubelementHandlingMode(FieldGroup.SUBELEMENT_HANDLING_MODE_FULL)
        destroyElementGroup = destroyGroup.createFieldElementGroup(mesh)
        destroyMesh = destroyElementGroup.getMeshGroup()
        for elementIdentifier in element_identifiers:
            element = mesh.findElementByIdentifier(elementIdentifier)
            destroyMesh.addElement(element)
        if destroyMesh.getSize() > 0:
            destroyNodeGroup = destroyGroup.getFieldNodeGroup(nodes)
            destroyNodes = destroyNodeGroup.getNodesetGroup()
            # must destroy elements first as Zinc won't destroy nodes that are in use
            mesh.destroyElementsConditional(destroyElementGroup)
            nodes.destroyNodesConditional(destroyNodeGroup)
            # clean up group so no external code hears is notified of its existence
            del destroyNodes
            del destroyNodeGroup
        del destroyMesh
        del destroyElementGroup
        del destroyGroup
    return


def extract_node_field_parameters(nodeset, field, only_value_labels=None):
    '''
    Returns parameters of field from nodes in nodeset in identifier order.
    Assumes all components have the same labels and versions.
    :param onlyValueLabels: Optional list of node value labels to limit extraction from e.g. [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1].
    :return: list of valueLabels returned, list of (node identifier, list over value labels of list of versions of parameters.
    '''
    fieldmodule = nodeset.getFieldmodule()
    componentsCount = field.getNumberOfComponents()
    fieldcache = fieldmodule.createFieldcache()
    valueLabels = only_value_labels if only_value_labels else \
        [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3, Node.VALUE_LABEL_D2_DS2DS3, Node.VALUE_LABEL_D3_DS1DS2DS3 ]
    valueLabelsCount = len(valueLabels)
    fieldParameters = []
    valueLabelParameterCounts = [ 0 for valueLabel in valueLabels ]
    nodeIter = nodeset.createNodeiterator()
    node = nodeIter.next()
    while node.isValid():
        fieldcache.setNode(node)
        nodeParameters = []
        fieldDefinedAtNode = False
        for i in range(valueLabelsCount):
            valueParameters = []
            version = 1
            while True:
                result, parameters = field.getNodeParameters(fieldcache, -1, valueLabels[i], version, componentsCount)
                if result != RESULT_OK:
                    break;
                fieldDefinedAtNode = True
                valueParameters.append(parameters)
                valueLabelParameterCounts[i] += 1
                version += 1
            nodeParameters.append(valueParameters)
        if fieldDefinedAtNode:
            fieldParameters.append( ( node.getIdentifier(), nodeParameters ) )
        node = nodeIter.next()
    for i in range(valueLabelsCount - 1, -1, -1):
        if valueLabelParameterCounts[i] == 0:
            valueLabels.pop(i)
            for nodeParameters in fieldParameters:
                nodeParameters[1].pop(i)
    return valueLabels, fieldParameters


def parameter_lists_to_string(valuesList, format_string):
    '''
    :return: 'None' if values is an empty list, the values in the first item if only one, otherwise the lists of values.
    '''
    if not valuesList:
        return 'None'
    if len(valuesList) == 1:
        return '[' + ','.join(format_string.format(value) for value in valuesList[0]) + ']'
    return '[ ' + ', '.join('[' + (','.join(format_string.format(value) for value in values)) for values in valuesList)  + '] ]'


def print_node_field_parameters(value_labels, node_field_parameters, format_string='{: 11e}'):
    '''
    Print value labels and parameters returned by extract_node_field_parameters ready for
    pasting into python code.
    :param format_string: Default format uses exponential notation with 7 significant digits = single precision.
    '''
    valueLabelsStrings = [ None, 'Node.VALUE_LABEL_VALUE', 'Node.VALUE_LABEL_D_DS1', 'Node.VALUE_LABEL_D_DS2', 'Node.VALUE_LABEL_D2_DS1DS2', 'Node.VALUE_LABEL_D_DS3', 'Node.VALUE_LABEL_D2_DS1DS3', 'Node.VALUE_LABEL_D2_DS2DS3', 'Node.VALUE_LABEL_D3_DS1DS2DS3']
    print('[ ' + ', '.join(valueLabelsStrings[v] for v in value_labels) + ' ]')
    print('[')
    lastNodeIdentifier  = node_field_parameters[-1][0]
    for nodeParameters in node_field_parameters:
        nodeIdentifier = nodeParameters[0]
        print('( ' + str(nodeIdentifier) + ', [ ' + ', '.join(parameter_lists_to_string(valueParameters, format_string) for valueParameters in nodeParameters[1]) + ' ] )' +  (', ' if (nodeIdentifier != lastNodeIdentifier) else ''))
    print(']\n')

def disconnectFieldMeshGroupBoundaryNodes(coordinateFields, meshGroup1, meshGroup2, nextNodeIdentifier):
    """
    Duplicate nodes in use by coordinateField on meshGroup1 and meshGroup2 and use these exclusively in meshGroup2.
    :param coordinateFields: The sequence of field to disconnect, type Zinc FieldFiniteElement. These are expected to
    be define with the same Elementfieldtemplate for all components in any element.
    :param meshGroup1: A Zinc MeshGroup containing elements in first group.
    :param meshGroup2: A Zinc MeshGroup containing elements in second group, which will be changed to use the new copies of nodes.
    :param nextNodeIdentifier: First available node identifier.
    :return: Final nextNodeIdentifier to use after this call, list of created node identifiers.
    """
    elemiter = meshGroup1.createElementiterator()
    element = elemiter.next()
    nodeIdentifiers1 = set()
    while element.isValid():
        eft = element.getElementfieldtemplate(coordinateFields[0], -1)
        nodeCount = eft.getNumberOfLocalNodes()
        for n in range(1, nodeCount + 1):
            node = element.getNode(eft, n)
            nodeIdentifiers1.add(node.getIdentifier())
        element = elemiter.next()
    copyIdentifiersMap = {}
    fieldmodule = coordinateFields[0].getFieldmodule()
    fieldcache = fieldmodule.createFieldcache()
    nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()

    allValueLabels = [
        Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
        Node.VALUE_LABEL_D2_DS1DS3, Node.VALUE_LABEL_D2_DS2DS3, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D3_DS1DS2DS3]

    with ChangeManager(fieldmodule):
        elemiter = meshGroup2.createElementiterator()
        element = elemiter.next()
        while element.isValid():
            eft = element.getElementfieldtemplate(coordinateFields[0], -1)
            nodeCount = eft.getNumberOfLocalNodes()
            for n in range(1, nodeCount + 1):
                existingNode = element.getNode(eft, n)
                existingNodeIdentifier = existingNode.getIdentifier()
                copyNodeIdentifier = copyIdentifiersMap.get(existingNodeIdentifier)
                copyNode = None

                if copyNodeIdentifier:
                    copyNode = nodes.findNodeByIdentifier(copyNodeIdentifier)
                else:
                    if existingNodeIdentifier in nodeIdentifiers1:
                        copyNodeIdentifier = nextNodeIdentifier
                        for coordinateField in coordinateFields:
                            nodetemplate.defineFieldFromNode(coordinateField, existingNode)
                        copyNode = nodes.createNode(nextNodeIdentifier, nodetemplate)
                        copyIdentifiersMap[existingNodeIdentifier] = copyNodeIdentifier
                        nextNodeIdentifier += 1

                        # copy field parameters from existingNode to copyNode
                        for coordinateField in coordinateFields:
                            componentsCount = coordinateField.getNumberOfComponents()

                            for valueLabel in allValueLabels:
                                versionCount = nodetemplate.getValueNumberOfVersions(coordinateField, -1, valueLabel)

                                for version in range(1, versionCount + 1):
                                    fieldcache.setNode(existingNode)
                                    result, values = coordinateField.getNodeParameters(
                                        fieldcache, -1, valueLabel, version, componentsCount)
                                    fieldcache.setNode(copyNode)
                                    coordinateField.setNodeParameters(
                                        fieldcache, -1, valueLabel, version, values)

                if copyNode:
                    result = element.setNode(eft, n, copyNode)
                    assert result == 1

            element = elemiter.next()

    return nextNodeIdentifier, list(copyIdentifiersMap.values())
