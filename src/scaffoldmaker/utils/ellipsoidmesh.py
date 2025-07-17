"""
Utilities for building solid ellipsoid meshes from hexahedral elements.
"""
import copy
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node


class EllipsoidMesh:
    """
    Generates a solid ellipsoid of hexahedral elements with oblique cross axes suited to describing lung geometry.
    """

    def __init__(self, element_counts, transition_element_count, sizes, axis1_rotation_radians, axis3_rotation_radians):
        assert all((count >= 4) and (count % 2 == 0) for count in element_counts)
        assert 1 <= transition_element_count <= (min(element_counts) // 2 - 1)
        self._element_counts = element_counts
        self._transition_element_count = transition_element_count
        self._sizes = sizes
        self._axis1_rotation_radians = axis1_rotation_radians
        self._axis3_rotation_radians = axis3_rotation_radians
        self._box_element_count1 = self._element_counts[0] - self._transition_element_count * 2
        self._box_element_count2 = self._element_counts[1] - self._transition_element_count * 2
        self._rim_element_count_around12 = (self._element_counts[0] + self._element_counts[1]) * 2
        parameters = [None] * 4  # x, d1, d2, d3
        self._box_coordinates = []  # over n3, n2, n1, d
        self._box_nids = []
        self._rim_coordinates = []  # over n3, n2 (through wall), n1 (around), d
        self._rim_nids = []
        for n3 in range(self._element_counts[2] + 1):
            box_coordinates_layer = [[copy.copy(parameters) for n1 in range(self._box_element_count1 + 1)]
                         for n2 in range(self._box_element_count2 + 1)]
            self._box_coordinates.append(box_coordinates_layer)
            box_nids_layer = [[None] * (self._box_element_count1 + 1) for n2 in range(self._box_element_count2 + 1)]
            self._box_nids.append(box_nids_layer)
            # no rim ring in start/end transition layers
            has_rim = self._transition_element_count <= n3 < (self._element_counts[2] - self._transition_element_count)
            rim_coordinates_layer = [[copy.copy(parameters) for nc in range(self._rim_element_count_around12)]
                                     for nr in range(self._transition_element_count)] if has_rim else None
            self._rim_coordinates.append(rim_coordinates_layer)
            rim_nids_layer = [[None] * self._rim_element_count_around12
                              for nr in range(self._transition_element_count)] if has_rim else None
            self._rim_nids.append(rim_nids_layer)

    def build(self):
        pass

    def generateMesh(self, fieldmodule, coordinates):
        """
        Client is expected to run within ChangeManager(fieldmodule)
        :param fieldmodule: Owning fieldmodule to create mesh in.
        :param coordinates: Coordinate field to define.
        """
        fieldcache = fieldmodule.createFieldcache()

        # create nodes

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        for value_label in [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3]:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, value_label, 1)

        node_identifier = 1
        for n3 in range(self._element_counts[2] + 1):
            for n2 in range(self._box_element_count2 + 1):
                for n1 in range(self._box_element_count1 + 1):
                    x, d1, d2, d3 = self._box_coordinates[n3][n2][n1]
                    if x:
                        node = nodes.createNode(node_identifier, nodetemplate)
                        self._box_nids[n3][n2][n1] = node_identifier
                        fieldcache.setNode(node)
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        if d1:
                            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        if d2:
                            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        if d3:
                            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        node_identifier += 1
            if self._rim_coordinates[n3]:
                for nr in range(self._transition_element_count):
                    for nc in range(self._rim_element_count_around12):
                        x, d1, d2, d3 = self._rim_coordinates[n3][nr][nc]
                        if x:
                            node = nodes.createNode(node_identifier, nodetemplate)
                            self._rim_nids[n3][nr][nc] = node_identifier
                            fieldcache.setNode(node)
                            coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                            if d1:
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                            if d2:
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                            if d3:
                                coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                            node_identifier += 1
