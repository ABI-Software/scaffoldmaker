import re
import math
import logging
import tempfile

from cmlibs.maths.vectorops import distance
from cmlibs.utils.zinc.field import get_group_list
from cmlibs.utils.zinc.finiteelement import get_element_node_identifiers
from cmlibs.utils.zinc.group import groups_have_same_local_contents
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.annotation.vagus_terms import (
    get_vagus_term, marker_name_in_terms, get_left_vagus_marker_locations_list, get_right_vagus_marker_locations_list)
from scaffoldmaker.utils.zinc_utils import get_nodeset_field_parameters


logger = logging.getLogger(__name__)


class VagusInputData:
    """
    Categorising and storing input data from data region for vagus box scaffold
    """

    def __init__(self, data_region):
        """
        :param data_region Zinc data region with input data
        """

        self._trunk_keywords = ['cervical vagus nerve', 'thoracic vagus nerve',
                                'cervical trunk', 'thoracic trunk', 'vagus x nerve trunk']
        self._branch_keywords = ['branch', 'nerve']
        self._non_branch_keywords = ['perineurium', 'epineurium']
        self._term_keywords = ['fma:', 'fma_', 'ilx:', 'ilx_', 'uberon:', 'uberon_']
        self._orientation_keywords = ['orientation']

        self._annotation_term_map = {}
        self._branch_coordinates_data = {}
        self._branch_parent_map = {}
        self._branch_common_group_map = {}
        self._branch_radius_data = {}
        self._datafile_path = None
        self._level_markers = {}
        self._orientation_data = {}
        self._side_label = ""
        self._trunk_group_name = None
        self._trunk_coordinates = []
        self._trunk_radius = []

        fm = data_region.getFieldmodule()
        fc = fm.createFieldcache()

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = fm.findFieldByName("coordinates").castFiniteElement()
        radius = fm.findFieldByName("radius").castFiniteElement()
        datapoints = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_DATAPOINTS)
        marker_names = fm.findFieldByName("marker_name")
        mesh = fm.findMeshByDimension(1)

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

        found_trunk_group_names = []
        branch_group_names = []
        orientation_group_names = []
        for annotation_name in annotation_names:
            lower_name = annotation_name.casefold()
            if any([keyword in lower_name for keyword in self._trunk_keywords]) and 'branch' not in lower_name:
                found_trunk_group_names.append(annotation_name)
            elif any([keyword in lower_name for keyword in self._branch_keywords]) and \
                    all([keyword not in lower_name for keyword in self._non_branch_keywords]):
                branch_group_names.append(annotation_name)
            elif any([keyword in lower_name for keyword in self._orientation_keywords]):
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

        # extract trunk data - coordinates, nodes, radius
        if len(found_trunk_group_names) > 0:
            if 'left' in found_trunk_group_names[0]:
                self._trunk_group_name = 'left vagus nerve'
                self._annotation_term_map[self._trunk_group_name] = get_vagus_term(self._trunk_group_name)[1]
                self._side_label = 'left'
            elif 'right' in found_trunk_group_names[0]:
                self._trunk_group_name = 'right vagus nerve'
                self._annotation_term_map[self._trunk_group_name] = get_vagus_term(self._trunk_group_name)[1]
                self._side_label = 'right'

        if self._trunk_group_name:
            trunk_group_count = 0
            trunk_elements = []
            for found_trunk_group_name in found_trunk_group_names:
                group = fm.findFieldByName(found_trunk_group_name).castGroup()
                nodeset = group.getNodesetGroup(nodes)
                meshgroup = group.getMeshGroup(mesh)
                _, values = get_nodeset_field_parameters(nodeset, coordinates, [Node.VALUE_LABEL_VALUE])
                if trunk_group_count == 0:
                    trunk_nodes = [value[0] for value in values]
                    trunk_coordinates = [value[1][0] for value in values]
                    if radius.isValid():
                        _, values = get_nodeset_field_parameters(nodeset, radius, [Node.VALUE_LABEL_VALUE])
                        trunk_radius = [value[1][0][0] for value in values]
                else:
                    trunk_nodes.extend([value[0] for value in values])
                    trunk_coordinates.extend([value[1][0] for value in values])
                    if radius.isValid():
                        _, values = get_nodeset_field_parameters(nodeset, radius, [Node.VALUE_LABEL_VALUE])
                        trunk_radius.extend([value[1][0][0] for value in values])

                # get trunk elements
                if meshgroup.getSize() > 0:
                    element_iterator = meshgroup.createElementiterator()
                    element = element_iterator.next()
                    while element.isValid():
                        eft = element.getElementfieldtemplate(coordinates, -1)
                        local_node_identifiers = get_element_node_identifiers(element, eft)
                        trunk_elements.append({'id': element.getIdentifier(),
                                               'nodes': local_node_identifiers})
                        element = element_iterator.next()

                trunk_group_count += 1

            # order trunk coordinates top to bottom in case trunk elements are available
            if len(trunk_elements) > 0:
                # build trunk graph
                nid_coords = {node_id: n_coord[0] for node_id, n_coord in zip(trunk_nodes, trunk_coordinates)}
                trunk_graph = {node_id: [] for node_id in nid_coords}
                for element in trunk_elements:
                    local_node_1, local_node_2 = element['nodes']
                    if local_node_1 in trunk_nodes and local_node_2 in trunk_nodes:
                        trunk_graph[local_node_1].append(local_node_2)
                        trunk_graph[local_node_2].append(local_node_1)
                # add any isolated nodes
                for node in trunk_nodes:
                    if node not in trunk_graph.keys():
                        trunk_graph[node] = []
                unconnected_nodes = []
                for el in trunk_graph.keys():
                    if len(trunk_graph[el]) <= 1:
                        unconnected_nodes.append(el)
                #print(trunk_graph)

                trunk_path_ids = []
                start = unconnected_nodes[0]
                start_index = 0
                # BFS from first to next unconnected, all connected in one long path
                while len(unconnected_nodes) > 0:
                    unconnected_nodes.pop(start_index)
                    if start not in trunk_path_ids:
                        local_trunk_path_ids = bfs_to_furthest(trunk_graph, start, trunk_path_ids)
                        trunk_path_ids.extend(local_trunk_path_ids)

                    # find next closest unconnected node
                    closest_distance = math.inf
                    last_node_id_in_path = trunk_path_ids[-1]
                    for index, node_id in enumerate(unconnected_nodes):
                        dist = distance(nid_coords[node_id], nid_coords[last_node_id_in_path])
                        if dist < closest_distance:
                            closest_distance = dist
                            start = node_id
                            start_index = index

                # get one of the top markers to check if trunk path needs to be reversed
                if self._side_label == 'left':
                    markers_ordered_list = get_left_vagus_marker_locations_list()
                else:
                    markers_ordered_list = get_right_vagus_marker_locations_list()

                for marker_name in markers_ordered_list.keys():
                    if marker_name in self._level_markers.keys():
                        top_marker = self._level_markers[marker_name]
                        break
                start_dist = distance(nid_coords[trunk_path_ids[0]], top_marker)
                end_dist = distance(nid_coords[trunk_path_ids[-1]], top_marker)
                if trunk_path_ids and end_dist < start_dist:
                    trunk_path_ids.reverse()

                ordered_trunk_coordinates = []
                for trunk_path_id in trunk_path_ids:
                    index = trunk_nodes.index(trunk_path_id)
                    ordered_trunk_coordinates.append(trunk_coordinates[index])

            if len(trunk_elements) > 0 and trunk_path_ids:
                self._trunk_coordinates = ordered_trunk_coordinates[:]
            else:
                self._trunk_coordinates = trunk_coordinates[:]

            if radius.isValid() and not all(value == 0.0 for value in trunk_radius):
                if len(trunk_elements) > 0 and trunk_path_ids:
                    ordered_trunk_radius = []
                    for trunk_path_id in trunk_path_ids:
                        index = trunk_nodes.index(trunk_path_id)
                        ordered_trunk_radius.append(trunk_radius[index])
                    self._trunk_radius = ordered_trunk_radius[:]
                else:
                    self._trunk_radius = trunk_radius[:]

        # extract branch data - name, coordinates, nodes, radius
        branch_nodes_data = {}
        for branch_name in branch_group_names:
            group = fm.findFieldByName(branch_name).castGroup()
            nodeset = group.getNodesetGroup(nodes)
            if 'xml.ex' in branch_name or nodeset.getSize() < 2:
                # xml.ex are temporary regions from segmentation stitcher that might have keywords in their names
                # branch should have at least two nodes to be connected to parent
                continue
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
            if parent_name == '':
                # assume trunk is a parent by default, if no other is found
                parent_name = self._trunk_group_name
            self._branch_parent_map[branch_name] = parent_name
            # print(branch_name, ' -> ', parent_name)

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

    def get_trunk_radius(self):
        """
        Get the radius values of the trunk in the data. The values are in the same order as trunk coordinates
        and each radius point is associated with a trunk coordinate.
        return: List of radius values for the trunk group.
        """
        return self._trunk_radius

    def get_branch_data(self):
        """
        Get all branch names and coordinates from the data.
        return: Dict mapping branch name to x, y, z data.
        """
        return self._branch_coordinates_data

    def get_branch_radius_data(self):
        """
        Get the radius values of the trunk in the data. The values are in the same order as trunk coordinates
        and each radius point is associated with a trunk coordinate.
        return: List of radius values for the trunk group.
        """
        return self._branch_radius_data

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

    def get_side_label(self):
        """
        Get label indicating side of the vagus (left or right or '')
        """
        return self._side_label

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
    Groups branches with the same annotations and destinations, only different by A, B, C variant character.
    :param branch_names: List with supplied branch names.
    :return branch_common_map: Dictionary mapping common branch name to list of branches with common names.
    Only contains entries for branches with variant names.
    """

    branch_common_map = {}
    for branch_name in branch_names:
        # remove single letters like A, B, C, etc. surrounded by whitespace
        common_key = re.sub(r'\b[A-Z]\b\s?', '', branch_name).strip()
        if common_key != branch_name:
            # branch variant was found
            variant_branch_list = branch_common_map.get(common_key)
            if variant_branch_list:
                variant_branch_list.append(branch_name)
            else:
                branch_common_map[common_key] = [branch_name]
    return branch_common_map


def load_vagus_data(region):
    """
    :param region: Zinc region for model definition.
    return: Provided the input file is supplied, it returns a data region with input data, otherwise None.
    """
    data_region = region.getParent().findChildByName('data')
    if not data_region.isValid():
        logger.warning("Missing input data.")
        return None

    vagus_data = VagusInputData(data_region)
    return vagus_data


def bfs_to_furthest(graph, start, trunk_path_ids):
    """
    :param graph:
    :param start:
    :param trunk_path_ids:
    return: Returns the furthest node and the path to it using BFS.
    """

    visited = set()
    parent = {start: None}
    queue = [start]
    last = start

    while queue:
        current = queue.pop(0)
        visited.add(current)
        last = current
        for neighbor in graph[current]:
            if neighbor not in visited and neighbor not in trunk_path_ids and neighbor not in queue:
                parent[neighbor] = current
                queue.append(neighbor)

    # Trace path from furthest node back to start
    path = []
    while last is not None:
        path.append(last)
        last = parent[last]
    return list(reversed(path))

