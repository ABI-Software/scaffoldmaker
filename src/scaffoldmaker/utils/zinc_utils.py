'''
Utility functions for easing use of Zinc API.
'''

from opencmiss.utils.zinc.field import getOrCreateFieldCoordinates
from opencmiss.utils.zinc.general import ZincCacheChanges
from opencmiss.zinc.context import Context
from opencmiss.zinc.element import MeshGroup
from opencmiss.zinc.field import Field
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
    with ZincCacheChanges(fieldmodule):
        cache = fieldmodule.createFieldcache()
        coordinates = getOrCreateFieldCoordinates(fieldmodule, components_count = componentsCount)
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
    with ZincCacheChanges(fieldmodule):
        isExterior = fieldmodule.createFieldIsExterior()
        isOnFace = fieldmodule.createFieldIsOnFace(elementFaceType)
        mesh2d = fieldmodule.findMeshByDimension(2)
        faceElementGroup = fieldmodule.createFieldElementGroup(mesh2d)
        faceMeshGroup = faceElementGroup.getMeshGroup()
        faceMeshGroup.addElementsConditional(fieldmodule.createFieldAnd(isExterior, isOnFace))
        del isExterior
        del isOnFace
    return faceMeshGroup
