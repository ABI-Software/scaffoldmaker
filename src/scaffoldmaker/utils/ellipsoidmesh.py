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

    def __init__(self, element_counts, transition_element_count, axis_lengths,
                 axis2_x_rotation_radians, axis3_x_rotation_radians):
        """

        :param element_counts:
        :param transition_element_count:
        :param axis_lengths: List of 3 ellipse axis lengths [a, b, c] in x, y, z direction
        :param axis2_x_rotation_radians: Rotation of axis 2 about +x direction
        :param axis3_x_rotation_radians: Rotation of axis 3 about +x direction.
        """
        assert all((count >= 4) and (count % 2 == 0) for count in element_counts)
        assert 1 <= transition_element_count <= (min(element_counts) // 2 - 1)
        self._element_counts = element_counts
        self._transition_element_count = transition_element_count
        self._axis_lengths = axis_lengths
        self._axis2_x_rotation_radians = axis2_x_rotation_radians
        self._axis3_x_rotation_radians = axis3_x_rotation_radians
        none_parameters = [None] * 4  # x, d1, d2, d3
        self._nx = []  # shield mesh with holes over n3, n2, n1, d
        self._nids = []
        half_counts = [count // 2 for count in self._element_counts]
        for n3 in range(self._element_counts[2] + 1):
            # index into transition zone
            trans3 = (self._transition_element_count - n3) if (n3 < half_counts[2]) else \
                (self._transition_element_count + n3 - self._element_counts[2])
            nx_layer = []
            nids_layer = []
            # print(n3, trans3)
            for n2 in range(self._element_counts[1] + 1):
                # index into transition zone
                trans2 = (self._transition_element_count - n2) if (n2 < half_counts[1]) else \
                    (self._transition_element_count + n2 - self._element_counts[1])
                nx_row = []
                nids_row = []
                # s = ""
                for n1 in range(self._element_counts[0] + 1):
                    # index into transition zone
                    trans1 = (self._transition_element_count - n1) if (n1 < half_counts[0]) else \
                        (self._transition_element_count + n1 - self._element_counts[0])
                    if (((trans1 <= 0) and (trans2 <= 0) and (trans3 <= 0)) or
                            (trans1 == trans2 == trans3) or
                            ((trans1 < 0) and ((trans2 == trans3) or (trans2 < 0))) or
                            ((trans2 < 0) and ((trans3 == trans1) or (trans3 < 0))) or
                            ((trans3 < 0) and ((trans1 == trans2) or (trans1 < 0)))):
                        parameters = copy.copy(none_parameters)
                        # s += "[]"
                    else:
                        parameters = None
                        # s += "  "
                    nx_row.append(parameters)
                    nids_row.append(None)
                nx_layer.append(nx_row)
                nids_layer.append(nids_row)
                # print(s)
            self._nx.append(nx_layer)
            self._nids.append(nids_layer)

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
            for n2 in range(self._element_counts[1] + 1):
                for n1 in range(self._element_counts[0] + 1):
                    parameters = self._nx[n3][n2][n1]
                    if not parameters:
                        continue
                    x, d1, d2, d3 = parameters
                    if not x:
                        continue  # while in development
                    node = nodes.createNode(node_identifier, nodetemplate)
                    self._nids[n3][n2][n1] = node_identifier
                    fieldcache.setNode(node)
                    coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    if d1:
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    if d2:
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    if d3:
                        coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                    node_identifier += 1
