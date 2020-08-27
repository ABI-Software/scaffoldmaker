'''
Utility functions for easing use of Zinc API.
'''

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import MeshGroup
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
