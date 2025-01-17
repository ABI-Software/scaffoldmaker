import os
import csv
import tempfile
import pandas as pd

from cmlibs.utils.zinc.field import get_group_list
from cmlibs.utils.zinc.group import groups_have_same_local_contents
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.annotation.vagus_terms import marker_name_in_terms


class VagusInputData:

    def __init__(self, data_region):

        self._trunk_keywords = ['vagus x nerve trunk', 'left vagus nerve', 'right vagus nerve']
        self._branch_keywords = ['branch', 'nerve']
        self._term_keywords = ['fma:', 'fma_', 'ilx:', 'ilx_', 'uberon:', 'uberon_']
        self._orientation_keywords = ['orientation', 'microfil']

        self._annotation_term_map = {}
        self._branch_coordinates_data = {}
        self._branch_parent_map = {}
        self._branch_radius_data = {}
        self._datafile_path = None
        self._level_markers = {}
        self._orientation_data = {}
        self._trunk_group_name = None
        self._trunk_coordinates = []
        self._trunk_radius = []

        fm = data_region.getFieldmodule()
        fc = fm.createFieldcache()

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = fm.findFieldByName("coordinates").castFiniteElement()
        assert coordinates.isValid() and (coordinates.getNumberOfComponents() == 3)
        radius = fm.findFieldByName("radius").castFiniteElement()
        markers = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        marker_names = fm.findFieldByName("marker_name")

        annotation_names = []
        term_annotation_names = []
        group_list = get_group_list(fm)
        group_map = {}
        for group in group_list:
            group_name = group.getName().strip()
            lower_name = group_name.casefold()
            group_map[group_name] = group
            if any([keyword in lower_name for keyword in self._term_keywords]):
                term_annotation_names.append(group_name)
            else:
                annotation_names.append(group_name)
        for annotation_name in annotation_names:
            annotation_group = fm.findFieldByName(annotation_name).castGroup()
            for term_annotation in term_annotation_names:
                term_group = fm.findFieldByName(term_annotation).castGroup()
                if groups_have_same_local_contents(annotation_group, term_group):
                    self._annotation_term_map[annotation_name] = term_annotation
                    break
            else:
                # no matching term is found for annotation group
                self._annotation_term_map[annotation_name] = ""

        branch_group_names = []
        orientation_group_names = []
        for annotation_name in annotation_names:
            lower_name = annotation_name.casefold()
            if any([keyword in lower_name for keyword in self._trunk_keywords]) and 'branch' not in lower_name:
                self._trunk_group_name = annotation_name
                continue
            if any([keyword in lower_name for keyword in self._branch_keywords]):
                branch_group_names.append(annotation_name)
            if any([keyword in lower_name for keyword in self._orientation_keywords]):
                # rename label used by Feinstein
                if lower_name == 'microfil':
                    annotation_name = 'orientation anterior'
                orientation_group_names.append(annotation_name)

        # extract marker data - name, coordinates (no marker terms are available)
        marker_group = group_map.get("marker")
        if marker_group:
            marker_nodes = marker_group.getNodesetGroup(markers)
            marker_node_iter = marker_nodes.createNodeiterator()
            marker_node = marker_node_iter.next()
            while marker_node.isValid():
                fc.setNode(marker_node)
                _, x = coordinates.evaluateReal(fc, 3)
                marker_name = marker_names.evaluateString(fc).strip()
                if marker_name_in_terms(marker_name):
                    # add anatomical landmark marker if in approved terms
                    self._level_markers[marker_name] = x
                marker_node = marker_node_iter.next()

        # extract orientation data
        for orientation_group_name in orientation_group_names:
            orientation_points, _ = get_nodeset_fieldgroup_parameters(nodes, coordinates, orientation_group_name,
                                                                           [Node.VALUE_LABEL_VALUE])
            self._orientation_data[orientation_group_name] = orientation_points[:]

        # extract trunk data - coordinates, nodes, radius - assume only one trunk group is used
        trunk_coordinates, trunk_nodes = get_nodeset_fieldgroup_parameters(nodes, coordinates, self._trunk_group_name,
                                                                           [Node.VALUE_LABEL_VALUE])
        self._trunk_coordinates = trunk_coordinates[:]

        # not used at the moment
        if radius.isValid():
            trunk_radius, _ = get_nodeset_fieldgroup_parameters(nodes, radius, self._trunk_group_name,
                                                                [Node.VALUE_LABEL_VALUE])
            if not all(value == 0.0 for value in trunk_radius):
                self._trunk_radius = trunk_radius[:]

        # project markers onto trunk if markers are projectable between first and last trunk node ?
        # not necessary to do that since we are using our own marker locations defined in a dictionary

        # extract branch data - name, coordinates, nodes, radius
        branch_nodes_data = {}
        for branch_name in branch_group_names:
            branch_parameters, branch_nodes = get_nodeset_fieldgroup_parameters(nodes, coordinates, branch_name,
                                                                                [Node.VALUE_LABEL_VALUE])
            self._branch_coordinates_data[branch_name] = branch_parameters
            branch_nodes_data[branch_name] = branch_nodes

            # not used at the moment
            if radius.isValid():
                branch_radius, _ = get_nodeset_fieldgroup_parameters(nodes, radius, branch_name,
                                                                     [Node.VALUE_LABEL_VALUE])
                if not all(value == 0.0 for value in branch_radius):
                    self._branch_radius_data[branch_name] = branch_radius

        # find parent branch where it connects to
        for branch_name, branch_nodes in branch_nodes_data.items():
            branch_first_node = sorted(branch_nodes)[0]

            #  first check if trunk is a parent by searching for a common node
            parent_name = ''
            if branch_first_node in trunk_nodes:
                parent_name = self._trunk_group_name
            else:
                # check other branches if a common node exists
                for parent_branch_name, parent_branch_nodes in branch_nodes_data.items():
                    if parent_branch_name != branch_name:
                        parent_first_node = sorted(parent_branch_nodes)[0]
                        if branch_first_node != parent_first_node and branch_first_node in parent_branch_nodes:
                            parent_name = parent_branch_name
                            break
            self._branch_parent_map[branch_name] = parent_name
            # print(branch_name, ' -> ', parent_name)

        # write all data in a file for geometry fitter
        sir = data_region.createStreaminformationRegion()
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            datafile_path = temp_file.name
            srf = sir.createStreamresourceFile(datafile_path)
            data_region.write(sir)
        self._datafile_path = datafile_path

    def get_level_markers(self):
        return self._level_markers

    def get_orientation_data(self):
        return self._orientation_data

    def get_trunk_group_name(self):
        return self._trunk_group_name

    def get_trunk_coordinates(self):
        return self._trunk_coordinates

    def get_branch_data(self):
        return self._branch_coordinates_data

    def get_annotation_term_map(self):
        return self._annotation_term_map

    def get_branch_parent_map(self):
        return self._branch_parent_map

    def get_datafile_path(self):
        return self._datafile_path

    def reset_datafile_path(self):
        self._datafile_path = None



def load_vagus_data(region):
    data_region = region.getParent().findChildByName('data')
    assert data_region.isValid(), "Invalid input data file"
    vagus_data = VagusInputData(data_region)
    return vagus_data


def get_nodeset_fieldgroup_parameters(nodeset, field, group_name, value_labels):
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
            for value_label in value_labels:
                _, parameters = finite_element_field.getNodeParameters(fieldcache, -1, value_label, 1, components_count)
                node_parameters.append(parameters)
            node_fieldgroup_parameters.append(node_parameters)
            node = node_iterator.next()

    return node_fieldgroup_parameters, nodes_list
