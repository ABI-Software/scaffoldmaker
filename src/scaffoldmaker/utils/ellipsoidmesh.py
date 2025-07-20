"""
Utilities for building solid ellipsoid meshes from hexahedral elements.
"""
import copy
import math
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.geometry import getEllipsePointAtTrueAngle, getEllipseTangentAtPoint, sampleCurveOnEllipsoid


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
        """
        Determine coordinates and derivatives over and within ellipsoid.
        """

        # get outside curve in 1-2 plane starting in 1 = x direction, 1st quadrant of 1-2 loop

        elements_count_q12 = (self._element_counts[0] // 2) + (self._element_counts[1] // 2) - 2
        start_x = [self._axis_lengths[0], 0.0, 0.0]
        start_d1 = [0.0, math.cos(self._axis2_x_rotation_radians), math.sin(self._axis2_x_rotation_radians)]
        start_d2 = [0.0, math.cos(self._axis3_x_rotation_radians), math.sin(self._axis3_x_rotation_radians)]
        end_x = [0.0] + getEllipsePointAtTrueAngle(
            self._axis_lengths[1], self._axis_lengths[2], self._axis2_x_rotation_radians)
        end_d1 = [-1.0, 0.0, 0.0]
        end_d2 = [0.0] + getEllipseTangentAtPoint(self._axis_lengths[1], self._axis_lengths[2], end_x[1:])
        px, pd1, pd2 = sampleCurveOnEllipsoid(self._axis_lengths[0], self._axis_lengths[1], self._axis_lengths[2],
                                              start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count_q12)
        self._setCoordinatesAroundRim(
            [px, pd1, pd2], [[0, 1, 2]],
            [self._element_counts[0], self._element_counts[1] // 2, self._element_counts[2] // 2],
            [[0, 1, 0], [-1, 0, 0]])
        # mirror and add 2nd quadrant of 1-2 loop
        for x, d1 in zip(px, pd1):
            x[0] = -x[0]
            d1[1] = -d1[1]
            d1[2] = -d1[2]
        self._setCoordinatesAroundRim(
            [px, pd1, pd2], [[0, 1, 2]],
            [0, self._element_counts[1] // 2, self._element_counts[2] // 2],
            [[0, 1, 0], [1, 0, 0]])
        # mirror and add 3rd quadrant of 1-2 loop
        for x, d1 in zip(px, pd1):
            x[1] = -x[1]
            x[2] = -x[2]
            d1[0] = -d1[0]
        self._setCoordinatesAroundRim(
            [px, pd1, pd2], [[0, 1, 2]],
            [0, self._element_counts[1] // 2, self._element_counts[2] // 2],
            [[0, -1, 0], [1, 0, 0]])
        # mirror and add 4th quadrant of 1-2 loop
        for x, d1 in zip(px, pd1):
            x[0] = -x[0]
            d1[1] = -d1[1]
            d1[2] = -d1[2]
        self._setCoordinatesAroundRim(
            [px, pd1, pd2], [[0, 1, 2]],
            [self._element_counts[0], self._element_counts[1] // 2, self._element_counts[2] // 2],
            [[0, -1, 0], [-1, 0, 0]])

        # get outside curve in 1-3 plane starting in 1 = x direction, 1st quadrant of 1-3 loop

        elements_count_q13 = (self._element_counts[0] // 2) + (self._element_counts[2] // 2) - 2
        start_x = [self._axis_lengths[0], 0.0, 0.0]
        start_d1 = [0.0, math.cos(self._axis2_x_rotation_radians), math.sin(self._axis2_x_rotation_radians)]
        start_d2 = [0.0, math.cos(self._axis3_x_rotation_radians), math.sin(self._axis3_x_rotation_radians)]
        end_x = [0.0] + getEllipsePointAtTrueAngle(
            self._axis_lengths[1], self._axis_lengths[2], self._axis3_x_rotation_radians)
        end_d2 = [-1.0, 0.0, 0.0]
        end_d1 = [0.0] + [-d for d in getEllipseTangentAtPoint(self._axis_lengths[1], self._axis_lengths[2], end_x[1:])]
        px, pd2, pd1 = sampleCurveOnEllipsoid(self._axis_lengths[0], self._axis_lengths[1], self._axis_lengths[2],
                                              start_x, start_d2, start_d1, end_x, end_d2, end_d1, elements_count_q13)
        self._setCoordinatesAroundRim(
            [px, pd1, pd2], [[0, 1, 2]],
            [self._element_counts[0], self._element_counts[1] // 2, self._element_counts[2] // 2],
            [[0, 0, 1], [-1, 0, 0]])
        # mirror and add 2nd quadrant of 1-3 loop
        for x, d1, d2 in zip(px, pd1, pd2):
            x[0] = -x[0]
            d1[0] = -d1[0]
            d2[0] = -d2[0]
        self._setCoordinatesAroundRim(
            [px, pd1, pd2], [[0, 1, 2], [0, 1, -2]],
            [0, self._element_counts[1] // 2, self._element_counts[2] // 2],
            [[0, 0, 1], [1, 0, 0]])
        # mirror and add 3rd quadrant of 1-3 loop
        for x, d1, d2 in zip(px, pd1, pd2):
            x[1] = -x[1]
            x[2] = -x[2]
            # d1[1] = -d1[1]
            # d1[2] = -d1[2]
            d2[0] = -d2[0]
        self._setCoordinatesAroundRim(
            [px, pd1, pd2], [[0, 1, 2], [0, -1, -2]],
            [0, self._element_counts[1] // 2, self._element_counts[2] // 2],
            [[0, 0, -1], [1, 0, 0]])
        # mirror and add 4th quadrant of 1-3 loop
        for x, d1, d2 in zip(px, pd1, pd2):
            x[0] = -x[0]
            d1[0] = -d1[0]
            d2[0] = -d2[0]
        self._setCoordinatesAroundRim(
            [px, pd1, pd2], [[0, 1, 2], [0, -1, -2]],
            [self._element_counts[0], self._element_counts[1] // 2, self._element_counts[2] // 2],
            [[0, 0, -1], [-1, 0, 0]])
        # centre
        n3 = self._element_counts[2] // 2
        n2 = self._element_counts[1] // 2
        n1 = self._element_counts[0] // 2
        self._nx[n3][n2][n1] = [
            [0.0, 0.0, 0.0], [self._axis_lengths[0], 0.0, 0.0], [0.0, self._axis_lengths[1], 0.0], [0.0, 0.0, self._axis_lengths[2]]]

    def _setCoordinatesAroundRim(self, parameters, parameter_indexes, start_indexes, index_increments):
        """
        Insert parameters around the rim into the coordinates array.
        :param parameters: List of lists of N node parameters e.g. [px, pd1, pd2]
        :param parameter_indexes: Lists of parameter indexes where x=0, d1=1, d2=2, d3=3. Starts with first and
        advances after each corner, then cycles back to first. Can be negative to invert vector.
        e.g. [[0, 1, 2], [0, -2, 1]] for [x, d1, d2] then [x, -d2, d1] after first corner.
        :param start_indexes: Index of first point.
        :param index_increments: List of increments in indexes. Starts with first and after at each corner, then
        cycles back to first.
        """
        half_counts = [count // 2 for count in self._element_counts]
        indexes = start_indexes
        parameter_number = 0
        parameter_index = parameter_indexes[0]
        increment_number = 0
        index_increment = index_increments[0]
        trans = [(self._transition_element_count - indexes[c]) if (indexes[c] < half_counts[c]) else
                 (self._transition_element_count + indexes[c] - self._element_counts[c]) for c in range(3)]
        max_trans = max(trans)
        for n in range(len(parameters[0])):
            for parameter, pix in zip(parameters, parameter_index):
                if pix < 0:
                    self._nx[indexes[2]][indexes[1]][indexes[0]][-pix] = [-d for d in parameter[n]]
                else:
                    self._nx[indexes[2]][indexes[1]][indexes[0]][pix] = copy.copy(parameter[n])
            trans = [(self._transition_element_count - indexes[c]) if (indexes[c] < half_counts[c]) else
                     (self._transition_element_count + indexes[c] - self._element_counts[c]) for c in range(3)]
            if trans.count(max_trans) > 1:
                parameter_number += 1
                if parameter_number == len(parameter_indexes):
                    parameter_number = 0
                parameter_index = parameter_indexes[parameter_number]
                increment_number += 1
                if increment_number == len(index_increments):
                    increment_number = 0
                index_increment = index_increments[increment_number]
            while True:
                indexes = [indexes[c] + index_increment[c] for c in range(3)]
                # skip over blank transition coordinates
                if self._nx[indexes[2]][indexes[1]][indexes[0]]:
                    break

    def generateMesh(self, fieldmodule, coordinates):
        """
        After build() has been called, generate nodes and elements of ellipsoid.
        Client is expected to run within ChangeManager(fieldmodule).
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
