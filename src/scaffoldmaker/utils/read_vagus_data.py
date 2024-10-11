import os
import csv
import tempfile
import pandas as pd

from cmlibs.utils.zinc.field import get_group_list, findOrCreateFieldCoordinates, findOrCreateFieldStoredString, \
    findOrCreateFieldGroup
from cmlibs.zinc.context import Context
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.annotation.vagus_terms import access_vagus_branch_terms, get_vagus_marker_term, marker_name_in_terms


def load_vagus_data(data_region):
    """
    Extract data related to vagus from supplied exf file,
    separate out data related to vagus trunk, vagus branches, markers (anatomical landmarks),
    does not extract fascicle data at the moment
    """

    # TODO: figure out the right location for the file
    termlist_location = "../../../packages/scaffoldmaker/src/scaffoldmaker/annotation/vagus_terms.xlsx"
    marker_terms, vagus_branch_terms = load_case_vagus_terms(termlist_location)

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
            if marker_name_in_terms(marker_name):
                marker_data[marker_name] = x
            else:
                # find marker name in CASE termslist
                marker_ilx = marker_terms[marker_name]
                # find marker name in our system
                if marker_ilx != 'NEW':
                    marker_term_name = get_vagus_marker_term(marker_ilx)
                    marker_data[marker_term_name[0]] = x
            marker_node = marker_node_iter.next()

    # extract trunk data - coordinates, nodes, radius - assume only one trunk group is used
    trunk_coordinates, trunk_nodes = get_nodeset_fieldgroup_parameters(nodes, coordinates, trunk_group_name, [Node.VALUE_LABEL_VALUE])
    if radius.isValid():
        trunk_radius, _ = get_nodeset_fieldgroup_parameters(nodes, radius, trunk_group_name, [Node.VALUE_LABEL_VALUE])
    else:
        trunk_radius = [2 for i in range(1, len(trunk_coordinates))]

    # TODO:
    # project markers onto trunk if markers are projectable between first and last trunk node
    # if they are higher or lower than data, leave them where they are


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
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        datafile_path = temp_file.name
        srf = sir.createStreamresourceFile(datafile_path)
        data_region.write(sir)

    return marker_data, trunk_group_name, trunk_coordinates, trunk_radius, branch_coordinates_data, branch_parents, \
           branch_radius_data, vagus_branch_terms, datafile_path


def load_case_vagus_terms(filename):
    def format_interlex_id(value):

        if isinstance(value, str) and 'http://uri.interlex.org/base/ilx_' in value:
            start = value.find('ilx_') + 4
            # Extract the numeric part after 'ilx_' and format it as 'ILX:XXXXXXX'
            return f"ILX:{value[start:]}"
        # If no URL is found, return the value as is
        return value

    def extract_between_markers(df, start_marker, end_marker):
        # Find the index of the row containing the start_marker
        start_row = df[df.eq(start_marker).any(axis=1)].index[0]

        # Find the index of the first row containing the end_marker
        end_row = df[df.eq(end_marker).any(axis=1)].index[0]

        # Select rows between start_row and end_row
        return df.loc[start_row + 1:end_row]

    df = pd.read_excel(filename, sheet_name='Full Termlist')
    df['ILX UUID'] = df['ILX UUID'].apply(format_interlex_id)
    df = df[['Unnamed: 0', 'ILX UUID']]

    # Select rows for 'REVA Gross Anatomy Landmarks and Levels'
    start_row = df[df.eq('REVA Gross Anatomy Landmarks and Levels').any(axis=1)].index[0]
    end_row = df[df.eq('REVA Bony Landmarks').any(axis=1)].index[0]
    filtered_vagus_levels = df.loc[start_row:end_row - 1]
    vagus_levels = filtered_vagus_levels.set_index('Unnamed: 0')['ILX UUID'].to_dict()

    # Select rows for 'REVA Vagus Structures'
    start_row = df[df.eq('REVA Vagus Structures').any(axis=1)].index[0]
    end_row = df[start_row:].isnull().all(axis=1).idxmax()
    filtered_vagus_structures = df.loc[start_row:end_row - 1]
    vagus_branches = filtered_vagus_structures.set_index('Unnamed: 0')['ILX UUID'].to_dict()

    vagus_branch_terms = access_vagus_branch_terms()
    vagus_branch_terms.extend(list(vagus_branches.items()))

    # merge two dictionaries
    #vagus_terms_dict = vagus_branches | vagus_levels
    #vagus_terms = access_vagus_terms()
    #vagus_terms.extend(list(vagus_terms_dict.items()))

    return vagus_levels, vagus_branch_terms


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