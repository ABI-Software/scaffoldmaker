import os
import csv

from cmlibs.utils.zinc.field import get_group_list, findOrCreateFieldCoordinates, findOrCreateFieldStoredString, \
    findOrCreateFieldGroup
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node


def load_exf_data(data_region):
    """
    Extract data related to vagus from supplied exf file,
    separate out data related to vagus trunk, vagus branches, markers (anatomical landmarks),
    does not extract fascicle data at the moment
    """

    fm = data_region.getFieldmodule()
    fc = fm.createFieldcache()

    group_list = get_group_list(fm)
    group_map = {}
    trunk_group_name = None
    branch_group_names = []
    trunk_keywords = ['vagus x nerve trunk', 'left vagus nerve', 'right vagus nerve']
    branch_keywords = ['branch', 'nerve']
    for group in group_list:
        group_name = group.getName()
        group_map[group_name] = group

        if any([keyword in group_name.lower() for keyword in trunk_keywords]) and 'branch' not in group_name.lower():
            trunk_group_name = group_name
            continue
        if any([keyword in group_name.lower() for keyword in branch_keywords]):
            branch_group_names.append(group_name)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)
    radius = fm.findFieldByName("radius").castFiniteElement()

    markers = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    marker_names = fm.findFieldByName("marker_name")

    # extract markers data - name, coordinates, radius (? - not now)
    marker_data = {}
    marker_group = group_map.get("marker")
    if marker_group:
        marker_nodes = marker_group.getNodesetGroup(markers)
        marker_node_iter = marker_nodes.createNodeiterator()
        marker_node = marker_node_iter.next()
        while marker_node.isValid():
            fc.setNode(marker_node)
            _, x = coordinates.evaluateReal(fc, 3)
            marker_name = marker_names.evaluateString(fc)
            marker_data[marker_name] = x
            marker_node = marker_node_iter.next()

    # extract trunk data - coordinates, nodes, radius - assume only one trunk group is used
    trunk_coordinates, trunk_nodes = get_nodeset_fieldgroup_parameters(nodes, coordinates, trunk_group_name, [Node.VALUE_LABEL_VALUE])
    if radius.isValid():
        trunk_radius, _ = get_nodeset_fieldgroup_parameters(nodes, radius, trunk_group_name, [Node.VALUE_LABEL_VALUE])
    else:
        trunk_radius = [2 for i in range(1, len(trunk_coordinates))]

    # TODO:
    # project markers onto trunk if markers are projectable between first and last trunk node
    # if they are higher or lower, leave them where they are


    # extract branch data - name, coordinates, nodes, radius
    branch_nodeslist = {}
    branch_coordinates_data = {}
    branch_radius_data = {}

    for branch_name in branch_group_names:
        branch_parameters, branch_nodes = get_nodeset_fieldgroup_parameters(nodes, coordinates, branch_name, [Node.VALUE_LABEL_VALUE])
        branch_coordinates_data[branch_name] = branch_parameters
        branch_nodeslist[branch_name] = branch_nodes

        if radius.isValid():
            branch_radius, _ = get_nodeset_fieldgroup_parameters(nodes, radius, branch_name, [Node.VALUE_LABEL_VALUE])
        else:
            branch_radius = [1 for i in range(1, len(branch_parameters))]
        branch_radius_data[branch_name] = branch_radius

    # find parent branch where it connects to
    branch_parents = {}
    for branch_name, branch_nodes in branch_nodeslist.items():
        branch_first_node = sorted(branch_nodes)[0]

        #  check if trunk is a parent by searching for a common node
        parent_name = ''
        if branch_first_node in trunk_nodes:
            parent_name = trunk_group_name
        else:
            for parent_branch_name, parent_branch_nodes in branch_nodeslist.items():
                parent_first_node = sorted(parent_branch_nodes)[0]
                if parent_branch_name != branch_name and branch_first_node != parent_first_node and branch_first_node in parent_branch_nodes:
                    parent_name = parent_branch_name
                    break
        branch_parents[branch_name] = parent_name
        # print('  ', branch_name, ' -> ', parent_name)

    # write all data in a file for geometry fitter
    sir = data_region.createStreaminformationRegion()
    srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_data.exf")
    data_region.write(sir)

    return marker_data, trunk_group_name, trunk_coordinates, trunk_radius, branch_coordinates_data, branch_parents, branch_radius_data


def load_exf_data_contours_japanese_dataset(region):
    """
    Extract data from supplied exf datafile, separate out data related to
    vagus trunk, vagus branches, fascicles, markers (anatomical landmarks)
    """

    data_region = region.getParent().findChildByName('data')
    assert data_region.isValid()

    fm = data_region.getFieldmodule()
    fc = fm.createFieldcache()

    group_list = get_group_list(fm)
    group_map = {}
    trunk_group_name = None
    branch_group_names = []
    branch_keywords = ['branch', 'nerve']
    for group in group_list:
        group_name = group.getName()
        group_map[group_name] = group
        if 'trunk' in group_name.lower():
            trunk_group_name = group_name
            continue
        if any(branch_keywords) in group_name.lower():
            branch_group_names.append(group_name)

    # extract marker data
    markers = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
    marker_coordinates = fm.findFieldByName("marker_data_coordinates").castFiniteElement()
    marker_name_field = fm.findFieldByName("marker_data_name")
    assert marker_coordinates.isValid() and (marker_coordinates.getNumberOfComponents() == 3)

    marker_data = {}
    marker_group = group_map.get("marker")
    if marker_group:
        marker_nodes = marker_group.getNodesetGroup(markers)
        marker_node_iter = marker_nodes.createNodeiterator()
        marker_node = marker_node_iter.next()
        while marker_node.isValid():
            fc.setNode(marker_node)
            result, x = marker_coordinates.evaluateReal(fc, 3)
            marker_name = marker_name_field.evaluateString(fc)
            marker_data[marker_name] = x
            marker_node = marker_node_iter.next()

    # extract other vagus data
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    coordinates = fm.findFieldByName("coordinates").castFiniteElement()
    assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)
    radius = fm.findFieldByName("radius").castFiniteElement()

    trunk_data_coordinates, trunk_nodes = get_nodeset_fieldgroup_parameters(nodes, coordinates, trunk_group_name, [Node.VALUE_LABEL_VALUE])
    if radius.isValid():
        trunk_radius, _ = get_nodeset_fieldgroup_parameters(nodes, radius, trunk_group_name, [Node.VALUE_LABEL_VALUE])
    else:
        trunk_radius = [2 for i in range(1, len(trunk_data_coordinates))]

    branch_data = {}
    branch_radius_data = {}
    for branch_name in branch_group_names:
        branch_parameters, branch_nodes = get_nodeset_fieldgroup_parameters(nodes, coordinates, branch_name, [Node.VALUE_LABEL_VALUE])
        branch_data[branch_name] = branch_parameters

        if radius.isValid():
            branch_radius, _ = get_nodeset_fieldgroup_parameters(nodes, radius, branch_name, [Node.VALUE_LABEL_VALUE])
        else:
            branch_radius = [1 for i in range(1, len(branch_parameters))]

        branch_radius_data[branch_name] = branch_radius

    # write all data in a file for geometry fitter
    sir = data_region.createStreaminformationRegion()
    srf = sir.createStreamresourceFile("C:/MAP/output/vagus_scaffold_temp/vagus_data.exf")
    data_region.write(sir)

    return marker_data, trunk_group_name, trunk_data_coordinates, trunk_radius, branch_data, branch_radius_data


def get_nodeset_fieldgroup_parameters(nodeset, field, group_name, valueLabels):
    """

    """

    fieldmodule = nodeset.getFieldmodule()
    finite_element_field = field.castFiniteElement()
    assert finite_element_field.isValid(), "get_nodeset_fieldgroup_parameters:  Field is not finite element type"

    components_count = field.getNumberOfComponents()
    fieldcache = fieldmodule.createFieldcache()

    group = fieldmodule.findFieldByName(group_name).castGroup()

    node_fieldgroup_parameters = []
    nodes_list = []
    if group.isValid():
        group_nodes = group.getNodesetGroup(nodeset)
        node_iterator = group_nodes.createNodeiterator()
        node = node_iterator.next()
        while node.isValid():
            fieldcache.setNode(node)
            nodes_list.append(node.getIdentifier())
            node_parameters = []
            field_defined_at_node = False
            for valueLabel in valueLabels:
                result, parameters = finite_element_field.getNodeParameters(fieldcache, -1, valueLabel, 1, components_count)
                field_defined_at_node = True
                node_parameters.append(parameters)
            node_fieldgroup_parameters.append(node_parameters)
            node = node_iterator.next()

    return node_fieldgroup_parameters, nodes_list