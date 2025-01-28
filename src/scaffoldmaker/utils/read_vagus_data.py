import os
import re
import csv
import tempfile
import pandas as pd

from cmlibs.utils.zinc.field import get_group_list
from cmlibs.utils.zinc.group import groups_have_same_local_contents
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.annotation.vagus_terms import marker_name_in_terms
from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters

class VagusInputData:
    """
    Categorising and storing input data from data region for vagus box scaffold
    """

    def __init__(self, data_region):
        """
        :param data_region Zinc data region with input data
        """

        self._trunk_keywords = ['vagus x nerve trunk', 'left vagus nerve', 'right vagus nerve']
        self._branch_keywords = ['branch', 'nerve']
        self._term_keywords = ['fma:', 'fma_', 'ilx:', 'ilx_', 'uberon:', 'uberon_']
        self._orientation_keywords = ['orientation', 'microfil']

        self._annotation_term_map = {}
        self._branch_coordinates_data = {}
        self._branch_parent_map = {}
        self._branch_common_group_map = {}
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
        datapoints = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        marker_names = fm.findFieldByName("marker_name")

        annotation_names = []
        term_annotation_names = []
        group_list = get_group_list(fm)
        for group in group_list:
            group_name = group.getName().strip()
            lower_name = group_name.casefold()
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
        # sort markers here?
        marker_group = fm.findFieldByName("marker").castGroup()
        if marker_group:
            marker_nodes = marker_group.getNodesetGroup(datapoints)
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
            group = fm.findFieldByName(orientation_group_name).castGroup()
            nodeset = group.getNodesetGroup(nodes)
            _, values = get_nodeset_field_parameters(nodeset, coordinates, [Node.VALUE_LABEL_VALUE])
            orientation_points = [value[1][0][0] for value in values]
            self._orientation_data[orientation_group_name] = orientation_points[:]

        # extract trunk data - coordinates, nodes, radius - assume only one trunk group is used
        group = fm.findFieldByName(self._trunk_group_name).castGroup()
        nodeset = group.getNodesetGroup(nodes)
        _, values = get_nodeset_field_parameters(nodeset, coordinates, [Node.VALUE_LABEL_VALUE])
        trunk_nodes = [value[0] for value in values]
        trunk_coordinates = [value[1][0] for value in values]
        self._trunk_coordinates = trunk_coordinates[:]

        # not used at the moment
        if radius.isValid():
            _, values = get_nodeset_field_parameters(nodeset, radius, [Node.VALUE_LABEL_VALUE])
            trunk_radius = [value[1][0][0] for value in values]
            if not all(value == 0.0 for value in trunk_radius):
                self._trunk_radius = trunk_radius[:]

        # project markers onto trunk if markers are projectable between first and last trunk node ?
        # not necessary to do that since we are using our own marker locations defined in a dictionary

        # extract branch data - name, coordinates, nodes, radius
        branch_nodes_data = {}
        for branch_name in branch_group_names:
            group = fm.findFieldByName(branch_name).castGroup()
            nodeset = group.getNodesetGroup(nodes)
            _, values = get_nodeset_field_parameters(nodeset, coordinates, [Node.VALUE_LABEL_VALUE])
            branch_nodes = [value[0] for value in values]
            branch_parameters = [value[1][0] for value in values]
            self._branch_coordinates_data[branch_name] = branch_parameters
            branch_nodes_data[branch_name] = branch_nodes

            # not used at the moment
            if radius.isValid():
                _, values = get_nodeset_field_parameters(nodeset, radius, [Node.VALUE_LABEL_VALUE])
                branch_radius = [value[1][0][0] for value in values]
                if not all(value == 0.0 for value in branch_radius):
                    self._branch_radius_data[branch_name] = branch_radius

        # find parent branch where it connects to
        for branch_name, branch_nodes in branch_nodes_data.items():
            # assumes trunk and branch node identifiers are strictly increasing.
            branch_first_node = branch_nodes[0]

            #  first check if trunk is a parent by searching for a common node
            parent_name = ''
            if branch_first_node in trunk_nodes:
                parent_name = self._trunk_group_name
            else:
                # check other branches if a common node exists
                for parent_branch_name, parent_branch_nodes in branch_nodes_data.items():
                    if parent_branch_name != branch_name:
                        parent_first_node = parent_branch_nodes[0]
                        if branch_first_node != parent_first_node and branch_first_node in parent_branch_nodes:
                            parent_name = parent_branch_name
                            break
            self._branch_parent_map[branch_name] = parent_name
            #print(branch_name, ' -> ', parent_name)

        # group common branches by names
        branch_common_map = group_common_branches(branch_group_names)
        self._branch_common_group_map = branch_common_map

        # write all data in a file for geometry fitter
        sir = data_region.createStreaminformationRegion()
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            datafile_path = temp_file.name
            srf = sir.createStreamresourceFile(datafile_path)
            data_region.write(sir)
        self._datafile_path = datafile_path

    def get_level_markers(self):
        """
        Get all level marker names and coordinates from the data.
        :return: Dict mapping marker name to x,y,z coordinates.
        """
        return self._level_markers

    def get_orientation_data(self):
        """
        Get all orientations and coordinates from the data.
        :return Dict mapping 8 possible orientations (anterior, left, right, etc.) to list of x, y, z coordinates.
        """
        return self._orientation_data

    def get_trunk_group_name(self):
        """
        Get the name used for the trunk group in the data.
        return: String with trunk group name.
        """
        return self._trunk_group_name

    def get_trunk_coordinates(self):
        """
        Get the x, y, z coordinates of the trunk in the data.
        return: List of coordinates for the trunk group.
        """
        return self._trunk_coordinates

    def get_branch_data(self):
        """
        Get all branch names and coordinates from the data.
        return: Dict mapping branch name to x, y, z data.
        """
        return self._branch_coordinates_data

    def get_annotation_term_map(self):
        """
        Get all annotation names and terms.
        return: Dict mapping annotation name to term annotation name.
        """
        return self._annotation_term_map

    def get_branch_common_group_map(self):
        """
        Get common branch names from the data.
        return: Dict mapping common branch name to list of branches with common names.
        """
        return self._branch_common_group_map

    def get_branch_parent_map(self):
        """
        Get all branch names and their parent branch names (for first level branches trunk is the parent, etc.).
        return: Dict mapping branch name to the parent branch or trunk name.
        """
        return self._branch_parent_map

    def get_datafile_path(self):
        """
        Get the path to the temporary file with the data.
        return: directory name of the temporary data file.
        """
        return self._datafile_path

    def reset_datafile_path(self):
        """
        Reset the path to the file with the data.
        """
        self._datafile_path = None


def group_common_branches(branch_names):
    """
    Groups branches with the same annotations and destinations, only different by A, B, C keyword.
    :param branch_names: List with supplied branch names.
    :param branch_keywords: List of branch keywords
    :return branch_common_map: Dictionary mapping common branch name to list of branches with common names.
    """

    grouped_names = {}
    for branch_name in branch_names:
        # remove single letters like A, B, C, etc.
        common_key = re.sub(r'\b[A-Z]\b\s?', '', branch_name).strip()
        if common_key not in grouped_names.keys():
            grouped_names[common_key] = [branch_name]
        else:
            grouped_names[common_key].append(branch_name)

    # make plural names
    branch_common_map = {}
    for name, groups in grouped_names.items():
        if 'branch' in name and len(groups) > 1:
            branch_name = 'all ' + name.replace('branch', 'branches', 1)
            branch_common_map[branch_name] = groups

    return branch_common_map




def load_vagus_data(region):
    """
    :param region: Zinc region for model definition.
    return: Provided the input file is supplied, it returns a data region with input data.
    """
    data_region = region.getParent().findChildByName('data')
    assert data_region.isValid(), "Invalid input data file"
    vagus_data = VagusInputData(data_region)
    return vagus_data
