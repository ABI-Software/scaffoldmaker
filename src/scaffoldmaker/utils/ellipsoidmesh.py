"""
Utilities for building solid ellipsoid meshes from hexahedral elements.
"""
import copy
import math
from cmlibs.maths.vectorops import add, cross, div, magnitude, mult, set_magnitude, sub
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.utils.eft_utils import determineCubicHermiteSerendipityEft, HermiteNodeLayoutManager
from scaffoldmaker.utils.geometry import (
    getEllipsePointAtTrueAngle, getEllipseTangentAtPoint, moveCoordinatesToEllipsoidSurface,
    moveDerivativeToEllipsoidSurface, sampleCurveOnEllipsoid)
from scaffoldmaker.utils.interpolation import (
    smoothCubicHermiteDerivativesLine)
from scaffoldmaker.utils.quadtrianglemesh import QuadTriangleMesh


class EllipsoidMesh:
    """
    Generates a solid ellipsoid of hexahedral elements with oblique cross axes suited to describing lung geometry.
    """

    def __init__(self, a, b, c, element_counts, transition_element_count,
                 axis2_x_rotation_radians, axis3_x_rotation_radians, surface_only=False):
        """
        :param a: Axis length (radius) in x direction.
        :param b: Axis length (radius) in y direction.
        :param c: Axis length (radius) in z direction.
        :param element_counts: Number of elements across full ellipse in a, b, c.
        :param transition_element_count: Number of transition elements around outside >= 1.
        :param axis2_x_rotation_radians: Rotation of axis 2 about +x direction
        :param axis3_x_rotation_radians: Rotation of axis 3 about +x direction.
        :param surface_only: Set to True to only make nodes and 2-D elements on the surface.
        """
        assert all((count >= 4) and (count % 2 == 0) for count in element_counts)
        assert 1 <= transition_element_count <= (min(element_counts) // 2 - 1)
        self._a = a
        self._b = b
        self._c = c
        self._element_counts = element_counts
        self._transition_element_count = transition_element_count
        self._axis2_x_rotation_radians = axis2_x_rotation_radians
        self._axis3_x_rotation_radians = axis3_x_rotation_radians
        self._surface_only = surface_only
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
        half_counts = [count // 2 for count in self._element_counts]
        elements_count_q12 = half_counts[0] + half_counts[1] - 2 * self._transition_element_count
        elements_count_q13 = half_counts[0] + half_counts[2] - 2 * self._transition_element_count
        elements_count_q23 = half_counts[1] + half_counts[2] - 2 * self._transition_element_count
        opp_element_counts = [half_counts[i] - self._transition_element_count for i in range(3)]

        axis1 = [self._a, 0.0, 0.0]
        axis2 = [0.0] + getEllipsePointAtTrueAngle(self._b, self._c, self._axis2_x_rotation_radians)
        axis3 = [0.0] + getEllipsePointAtTrueAngle(self._b, self._c, self._axis3_x_rotation_radians)
        axis_d1 = div(axis1, half_counts[0])
        axis_d2 = div(axis2, half_counts[1])
        axis_d3 = div(axis3, half_counts[2])

        sample_curve_on_ellipsoid = (
            lambda start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count,
                   start_weight=None, end_weight=None, overweighting=1.0, end_transition=False:
            sampleCurveOnEllipsoid(
                self._a, self._b, self._c, start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count,
                start_weight, end_weight, overweighting, end_transition))
        move_x_to_ellipsoid_surface = lambda x: moveCoordinatesToEllipsoidSurface(self._a, self._b, self._c, x)
        move_d_to_ellipsoid_surface = lambda x, d: moveDerivativeToEllipsoidSurface(self._a, self._b, self._c, x, d)
        d3_mag = min(magnitude(axis_d1), magnitude(axis_d2), magnitude(axis_d3))
        # surface normal:
        # evaluate_ellipsoid_surface_d3 = lambda tx, td1, td2: set_magnitude(
        #     [tx[0] / (self._a * self._a), tx[1] / (self._b * self._b), tx[2] / (self._c * self._c)], d3_mag)
        evaluate_ellipsoid_surface_d3 = lambda tx, td1, td2: set_magnitude(tx, d3_mag)

        octant1 = EllipsoidOctantMesh(self._a, self._b, self._c, half_counts, self._transition_element_count)

        # get outside curve from axis 1 to axis 2
        abx, abd1, abd2 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            axis1, axis_d2, axis_d3,
            axis2, [-1.0, 0.0, 0.0], [0.0] + getEllipseTangentAtPoint(self._b, self._c, axis2[1:]),
            elements_count_q12)
        # get outside curve from axis 1 to axis 3
        acx, acd2, acd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            abx[0], abd2[0], abd1[0], axis3, [-1.0, 0.0, 0.0],
            [0.0] + [-d for d in getEllipseTangentAtPoint(self._b, self._c, axis3[1:])],
            elements_count_q13)
        # get outside curve from axis 2 to axis 3
        bcx, bcd2, bcd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            abx[-1], abd2[-1], abd1[-1],
            acx[-1], [-d for d in acd1[-1]], acd2[-1],
            elements_count_q23)
        # fix first/last derivatives so correct regardless of order set in QuadTriangleMesh
        abd2[0] = acd2[0]
        abd2[-1] = bcd2[0]
        acd1[-1] = [-d for d in bcd2[-1]]

        # make outer surface of octant 1
        outer_surface1 = QuadTriangleMesh(
            opp_element_counts[0], opp_element_counts[1], opp_element_counts[2],
            sample_curve_on_ellipsoid, move_x_to_ellipsoid_surface, move_d_to_ellipsoid_surface)
        outer_surface1.set_edge_parameters12(abx, abd1, abd2)
        outer_surface1.set_edge_parameters13(acx, acd1, acd2)
        outer_surface1.set_edge_parameters23(bcx, bcd1, bcd2)
        outer_surface1.build()
        if not self._surface_only:
            d3_mag = min(magnitude(axis_d1), magnitude(axis_d2), magnitude(axis_d3))
            outer_surface1.calculate_d3(evaluate_ellipsoid_surface_d3)
        octant1.set_outer_surface(outer_surface1)

        # copy octant1 into ellipsoid
        octant_parameters = octant1.get_parameters()
        for n3 in range(half_counts[2] + 1):
            octant_nx_layer = octant_parameters[n3]
            nx_layer = self._nx[half_counts[2] + n3]
            for n2 in range(half_counts[1] + 1):
                octant_nx_row = octant_nx_layer[n2]
                nx_row = nx_layer[half_counts[1] + n2]
                for n1 in range(half_counts[0] + 1):
                    nx_row[half_counts[0] + n1] = copy.deepcopy(octant_nx_row[n1])

        # transfer parameters on bottom plane of octant1 to octant2 for blending
        n3 = half_counts[2]
        for n2 in range(half_counts[1] + 1, self._element_counts[1] + 1):
            for n1 in range(half_counts[0], self._element_counts[0] + 1):
                parameters = self._nx[n3][n2][n1]
                if not parameters or not parameters[2]:
                    continue
                x = parameters[0]
                d1 = parameters[1]
                d2 = parameters[2]
                d3 = parameters[3]
                x = [x[0], -x[1], -x[2]]
                d1 = [-d1[0], d1[1], d1[2]]
                d2 = [-d2[0], d2[1], d2[2]]
                d3 = [d3[0], -d3[1], -d3[2]] if d3 else None  # GRC fix
                self._nx[n3][self._element_counts[1] - n2][n1] = [x, d1, d2, d3]

        octant2 = EllipsoidOctantMesh(
            self._b, self._a, self._c, [half_counts[1], half_counts[0], half_counts[2]], self._transition_element_count)

        # get outside curve from axis -2 to axis 1, mirrored from ab curve
        mbax = [[x[0], -x[1], -x[2]] for x in reversed(abx)]
        mbad1 = [[-d1[0], d1[1], d1[2]] for d1 in reversed(abd1)]
        mbad2 = [[-d2[0], d2[1], d2[2]] for d2 in reversed(abd2)]
        # get outside curve from axis -2 to axis 3
        mbcx, mbcd2, mbcd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            mbax[0], mbad2[0], mbad1[0],
            acx[-1], acd1[-1], [-d for d in acd2[-1]],
            elements_count_q23)

        # make outer surface of octant 2
        outer_surface2 = QuadTriangleMesh(
            opp_element_counts[1], opp_element_counts[0], opp_element_counts[2],
            sample_curve_on_ellipsoid, move_x_to_ellipsoid_surface, move_d_to_ellipsoid_surface)
        outer_surface2.set_edge_parameters12(mbax, mbad1, mbad2)
        outer_surface2.set_edge_parameters13(mbcx, mbcd1, mbcd2)
        outer_surface2.set_edge_parameters23(acx, acd1, acd2)
        outer_surface2.build()
        if not self._surface_only:
            d3_mag = min(magnitude(axis_d1), magnitude(axis_d2), magnitude(axis_d3))
            outer_surface2.calculate_d3(evaluate_ellipsoid_surface_d3)
        octant2.set_outer_surface(outer_surface2)

        # copy octant2 into ellipsoid, blending existing derivatives
        octant_parameters = octant2.get_parameters()
        for n3 in range(half_counts[2] + 1):
            octant_nx_layer = octant_parameters[n3]
            nx_layer = self._nx[half_counts[2] + n3]
            for n2 in range(half_counts[0] + 1):
                octant_nx_row = octant_nx_layer[n2]
                for n1 in range(half_counts[1] + 1):
                    octant_nx = octant_nx_row[n1]
                    if octant_nx:
                        nx = nx_layer[half_counts[1] - n1][half_counts[0] + n2]
                        nx[0] = copy.copy(octant_nx[0])
                        for pix in range(1, 4):
                            if (n3 < half_counts[2]) or (pix == 3):
                                d = octant_nx[pix]
                            else:
                                d = octant_nx[2] if (pix == 1) else [-d for d in octant_nx[1]] if octant_nx[1] else None
                            if d:
                                if nx[pix]:
                                    # expect derivatives to be the same direction and blend magnitude with harmonic mean
                                    mean_mag = 2.0 / ((1.0 / magnitude(nx[pix])) + (1.0 / magnitude(d)))
                                    nx[pix] = set_magnitude(nx[pix], mean_mag)
                                else:
                                    nx[pix] = copy.copy(d)

        # transfer blended d2 derivatives on bottom plane of octant2 back to octant1
        n3 = half_counts[2]
        for n2 in range(half_counts[1] + 1, self._element_counts[1] + 1):
            for n1 in range(half_counts[0], self._element_counts[0] + 1):
                parameters = self._nx[n3][self._element_counts[1] - n2][n1]
                if not parameters or not parameters[2]:
                    continue
                d2 = parameters[2]
                self._nx[n3][n2][n1][2] = [-d2[0], d2[1], d2[2]]

        # mirror octants over x = 0
        for n3 in range(half_counts[2], self._element_counts[2] + 1):
            for n2 in range(0, self._element_counts[1] + 1):
                for n1 in range(half_counts[0] + 1, self._element_counts[0] + 1):
                    parameters = self._nx[n3][n2][n1]
                    if not parameters or not parameters[0]:
                        continue
                    x, d1, d2, d3 = parameters
                    x = [-x[0], x[1], x[2]]
                    d1 = [d1[0], -d1[1], -d1[2]] if d1 else None
                    d2 = [-d2[0], d2[1], d2[2]] if d2 else None
                    d3 = [-d3[0], d3[1], d3[2]] if d3 else None
                    self._nx[n3][n2][self._element_counts[0] - n1] = [x, d1, d2, d3]

        # flip top half about both y = 0 and z = 0
        for n3 in range(half_counts[2], self._element_counts[2] + 1):
            for n2 in range(0, self._element_counts[1] + 1):
                for n1 in range(self._element_counts[0] + 1):
                    parameters = self._nx[n3][n2][n1]
                    if not parameters or not parameters[0]:
                        continue
                    x, d1, d2, d3 = parameters
                    x = [x[0], -x[1], -x[2]]
                    d1 = [-d1[0], d1[1], d1[2]] if d1 else None
                    d2 = [-d2[0], d2[1], d2[2]] if d2 else None
                    d3 = [d3[0], d3[1], -d3[2]] if d3 else None
                    self._nx[self._element_counts[2] - n3][self._element_counts[1] - n2][n1] = [x, d1, d2, d3]

    def _next_increment_out_of_bounds(self, indexes, index_increment):
        for c in range(3):
            index = indexes[c] + index_increment[c]
            if (index < 0) or (index > self._element_counts[c]):
                return True
        return False

    def _set_coordinates_around_rim(self, parameters, parameter_indexes, start_indexes, index_increments,
                                    skip_start=False, skip_end=False,
                                    blend_start=False, blend_middle=False, blend_end=False):
        """
        Insert parameters around the rim into the coordinates array.
        :param parameters: List of lists of N node parameters e.g. [px, pd1, pd2]
        :param parameter_indexes: Lists of parameter indexes where x=0, d1=1, d2=2, d3=3. Starts with first and
        advances after each corner, then cycles back to first. Can be negative to invert vector.
        e.g. [[0, 1, 2], [0, -2, 1]] for [x, d1, d2] then [x, -d2, d1] after first corner.
        :param start_indexes: Index of first point.
        :param index_increments: List of increments in indexes. Starts with first and after at each corner, then
        cycles back to first.
        :param skip_start: Set to True to skip the first value.
        :param skip_end: Set to True to skip the last value.
        :param blend_start: Set to True to blend parameters with any old parameters at start location.
        :param blend_middle: Set to True to blend parameters with any old parameters at middle locations.
        :param blend_end: Set to True to blend parameters with any old parameters at end location.
        """
        indexes = start_indexes
        parameter_number = 0
        parameter_index = parameter_indexes[0]
        increment_number = 0
        index_increment = index_increments[0]
        start_n = 1 if skip_start else 0
        last_n = len(parameters[0]) - 1
        limit_n = len(parameters[0]) - (1 if skip_end else 0)
        for n in range(start_n, limit_n):
            if n > 0:
                while True:
                    indexes = [indexes[c] + index_increment[c] for c in range(3)]
                    # skip over blank transition coordinates
                    if self._nx[indexes[2]][indexes[1]][indexes[0]]:
                        break
            if self._next_increment_out_of_bounds(indexes, index_increment):
                parameter_number += 1
                if parameter_number == len(parameter_indexes):
                    parameter_number = 0
                parameter_index = parameter_indexes[parameter_number]
                increment_number += 1
                if increment_number == len(index_increments):
                    increment_number = 0
                index_increment = index_increments[increment_number]
            nx = self._nx[indexes[2]][indexes[1]][indexes[0]]
            for parameter, spix in zip(parameters, parameter_index):
                new_parameter = [-d for d in parameter[n]] if (spix < 0) else copy.copy(parameter[n])
                pix = abs(spix)
                if nx[pix] and ((blend_start and (n == 0)) or (blend_middle and (0 < n < last_n)) or
                        (blend_end and (n == last_n))):
                    if pix == 0:
                        x = [0.5 * (nx[pix][c] + new_parameter[c]) for c in range(3)]
                        nx[pix] = moveCoordinatesToEllipsoidSurface(self._a, self._b, self._c, x)
                    else:
                        # expect derivatives to be the same direction and blend magnitude with harmonic mean
                        mean_mag = 2.0 / ((1.0 / magnitude(nx[pix])) + (1.0 / magnitude(new_parameter)))
                        nx[pix] = set_magnitude(nx[pix], mean_mag)
                else:
                    nx[pix] = new_parameter

    def _smooth_derivatives_around_rim(self, start_indexes, end_indexes, index_increments,
                                       derivative_indexes, end_derivative_index,
                                       fix_start_direction=True, fix_end_direction=True,
                                       blend_start=False, blend_end=False):
        """
        Smooth derivatives around the rim in the coordinates array.
        :param start_indexes: Indexes of first point.
        :param end_indexes: Indexes of last point.
        :param index_increments: List of increments in indexes. Starts with first and after at each corner, then
        cycles back to first.
        :param derivative_indexes: List of signed derivative parameter index to along where d1=1, d2=2, d3=3.
        Starts with first and advances after each corner, then cycles back to first. Can be negative to invert vector.
        e.g. [1, -2] for d1 then -d2 after first corner.
        :param end_derivative_index: List of signed derivative indexes to apply on the last point
        e.g. [1, -2] gives d1 - d2.
        :param fix_start_direction: Set to True to keep the start direction but scale its magnitude.
        :param fix_end_direction: Set to True to keep the end direction but scale its magnitude.
        :param blend_start: Set to True to 50:50 blend parameters with any old parameters at start location.
        :param blend_end: Set to True to 50:50 blend parameters with any old parameters at end location.
        """
        indexes = start_indexes
        derivative_number = 0
        derivative_index = derivative_indexes[0]
        increment_number = 0
        index_increment = index_increments[0]
        indexes_list = []
        derivative_index_list = []
        px = []
        pd = []
        n = 0
        while True:
            if n > 0:
                if indexes == end_indexes:
                    break
                while True:
                    indexes = [indexes[c] + index_increment[c] for c in range(3)]
                    # skip over blank transition coordinates
                    if self._nx[indexes[2]][indexes[1]][indexes[0]]:
                        break
            if self._next_increment_out_of_bounds(indexes, index_increment):
                derivative_number += 1
                if derivative_number == len(derivative_indexes):
                    derivative_number = 0
                derivative_index = derivative_indexes[derivative_number]
                increment_number += 1
                if increment_number == len(index_increments):
                    increment_number = 0
                index_increment = index_increments[increment_number]
            parameters = self._nx[indexes[2]][indexes[1]][indexes[0]]
            x = parameters[0]
            use_derivative_index = end_derivative_index if (indexes == end_indexes) else [derivative_index]
            indexes_list.append(copy.copy(indexes))
            derivative_index_list.append(copy.copy(use_derivative_index))
            d = [0.0, 0.0, 0.0]
            for i in range(len(use_derivative_index)):
                spix = use_derivative_index[i]
                pix = abs(spix)
                values = parameters[pix]
                if values:
                    if spix < 0:
                        values = [-ad for ad in values]
                    d = add(d, values)
            px.append(x)
            pd.append(d)
            n += 1
        sd = smoothCubicHermiteDerivativesLine(
            px, pd, fixStartDirection=fix_start_direction, fixEndDirection=fix_end_direction,
            fixEndDerivative=len(end_derivative_index) > 1)
        for n in range(len(sd)):
            sd[n] = moveDerivativeToEllipsoidSurface(self._a, self._b, self._c, px[n], sd[n])
        sd = smoothCubicHermiteDerivativesLine(
            px, sd, fixAllDirections=True, fixEndDerivative=len(end_derivative_index) > 1)
        last_n = len(sd) - 1
        for n in range(len(sd)):
            indexes = indexes_list[n]
            derivative_index = derivative_index_list[n]
            parameters = self._nx[indexes[2]][indexes[1]][indexes[0]]
            if len(derivative_index) == 1:
                spix = derivative_index[0]
                new_derivative = [-d for d in sd[n]] if (spix < 0) else sd[n]
                pix = abs(spix)
                if parameters[pix] and (blend_start and (n == 0)) or (blend_end and (n == last_n)):
                    # expect derivatives to be the same direction and blend magnitude with harmonic mean
                    mean_mag = 2.0 / ((1.0 / magnitude(parameters[pix])) + (1.0 / magnitude(new_derivative)))
                    parameters[pix] = set_magnitude(parameters[pix], mean_mag)
                else:
                    parameters[pix] = new_derivative
            # else:
            #     # not putting back values if summed parameters

    def generate_mesh(self, fieldmodule, coordinates, start_node_identifier=1, start_element_identifier=1):
        """
        After build() has been called, generate nodes and elements of ellipsoid.
        Client is expected to run within ChangeManager(fieldmodule).
        :param fieldmodule: Owning fieldmodule to create mesh in.
        :param coordinates: Coordinate field to define.
        :param start_node_identifier: Identifier to define for first node and incremented thereafter.
        :param start_element_identifier: Identifier to define for first element and incremented thereafter.
        Note this is for mesh2d if surface_only, otherwise mesh3d.
        return next node identifier, next element identifier for objects after this.
        """
        fieldcache = fieldmodule.createFieldcache()

        # create nodes

        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        value_labels = [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2]
        if not self._surface_only:
            value_labels.append(Node.VALUE_LABEL_D_DS3)
        for value_label in value_labels:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, value_label, 1)

        node_identifier = start_node_identifier
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

        node_layout_manager = HermiteNodeLayoutManager()
        element_identifier = start_element_identifier

        if self._surface_only:
            mesh2d = fieldmodule.findMeshByDimension(2)
            elementtemplate_regular = mesh2d.createElementtemplate()
            elementtemplate_regular.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            bicubic_hermite_serendipity_basis = (
                fieldmodule.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY))
            eft_regular = mesh2d.createElementfieldtemplate(bicubic_hermite_serendipity_basis)
            elementtemplate_regular.defineField(coordinates, -1, eft_regular)
            elementtemplate_special = mesh2d.createElementtemplate()
            elementtemplate_special.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            half_counts = [count // 2 for count in self._element_counts]
            rim_indexes = [[0] + [self._transition_element_count + 1 + j
                                  for j in range(self._element_counts[i] - 2 * self._transition_element_count - 1)] +
                           [self._element_counts[i]] for i in range(3)]
            # bottom rectangle
            bottom_nids = self._nids[0]
            last_nids_row = None
            for i2, n2 in enumerate(rim_indexes[1]):
                nids_row = []
                for i1, n1 in enumerate(reversed(rim_indexes[0])):
                    nids_row.append(bottom_nids[n2][n1])
                    if (i2 > 0) and (i1 > 0):
                        nids = [last_nids_row[i1 - 1], last_nids_row[i1], nids_row[i1 - 1], nids_row[i1]]
                        if None in nids:
                            continue
                        element = mesh2d.createElement(element_identifier, elementtemplate_regular)
                        element.setNodesByIdentifier(eft_regular, nids)
                        # print("Element", element_identifier, "nids", nids)
                        element_identifier += 1
                last_nids_row = nids_row
            # around sides
            node_layout_permuted = node_layout_manager.getNodeLayoutRegularPermuted(d3Defined=False)
            node_layout_triple_points = node_layout_manager.getNodeLayoutTriplePoint2D()
            index_increments = [[0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 0, 0]]
            increment_number = 0
            index_increment = index_increments[0]
            elements_count_around12 = \
                2 * (self._element_counts[0] + self._element_counts[1] - 4 * self._transition_element_count)
            last_nids_row = None
            last_parameters_row = None
            last_corners_row = None
            for n3 in rim_indexes[2]:
                indexes = [self._element_counts[0], half_counts[1], n3]
                nids_row = []
                parameters_row = []
                corners_row = []
                for n in range(elements_count_around12):
                    if n > 0:
                        while True:
                            indexes = [indexes[c] + index_increment[c] for c in range(3)]
                            # skip over blank transition coordinates
                            if self._nx[indexes[2]][indexes[1]][indexes[0]]:
                                break
                    nids_row.append(self._nids[indexes[2]][indexes[1]][indexes[0]])
                    parameters_row.append(self._nx[indexes[2]][indexes[1]][indexes[0]])
                    if self._next_increment_out_of_bounds(indexes, index_increment):
                        corners_row.append(True)
                        increment_number += 1
                        if increment_number == len(index_increments):
                            increment_number = 0
                        index_increment = index_increments[increment_number]
                    else:
                        corners_row.append(False)
                if n3 > 0:
                    quarter_elements_count_around12 = elements_count_around12 // 4
                    for n in range(elements_count_around12):
                        nids = [last_nids_row[n], last_nids_row[(n + 1) % elements_count_around12],
                                nids_row[n], nids_row[(n + 1) % elements_count_around12]]
                        if None in nids:
                            continue
                        q = n // quarter_elements_count_around12
                        elementtemplate = elementtemplate_regular
                        eft = eft_regular
                        scalefactors = None
                        if n3 in (rim_indexes[2][1], self._element_counts[2]):
                            node_parameters = [
                                last_parameters_row[n], last_parameters_row[(n + 1) % elements_count_around12],
                                parameters_row[n], parameters_row[(n + 1) % elements_count_around12]]
                            if n3 == rim_indexes[2][1]:
                                node_layouts = [node_layout_permuted, node_layout_permuted, None, None]
                                if last_corners_row[n]:
                                    node_layouts[0] = node_layout_triple_points[(5 - q) % 4]
                                elif last_corners_row[(n + 1) % elements_count_around12]:
                                    node_layouts[1] = node_layout_triple_points[(5 - q) % 4]
                            else:
                                node_layouts = [None, None, node_layout_permuted, node_layout_permuted]
                                if corners_row[n]:
                                    node_layouts[2] = node_layout_triple_points[q]
                                elif corners_row[(n + 1) % elements_count_around12]:
                                    node_layouts[3] = node_layout_triple_points[q]
                            eft, scalefactors = \
                                determineCubicHermiteSerendipityEft(mesh2d, node_parameters, node_layouts)
                            elementtemplate_special.defineField(coordinates, -1, eft)
                            elementtemplate = elementtemplate_special
                        element = mesh2d.createElement(element_identifier, elementtemplate)
                        element.setNodesByIdentifier(eft, nids)
                        if scalefactors:
                            element.setScaleFactors(eft, scalefactors)
                        # print("Element", element_identifier, "nids", nids)
                        element_identifier += 1
                last_nids_row = nids_row
                last_parameters_row = parameters_row
                last_corners_row = corners_row
            # top rectangle
            top_nids = self._nids[self._element_counts[2]]
            last_nids_row = None
            for i2, n2 in enumerate(rim_indexes[1]):
                nids_row = []
                for i1, n1 in enumerate(rim_indexes[0]):
                    nids_row.append(top_nids[n2][n1])
                    if (i2 > 0) and (i1 > 0):
                        nids = [last_nids_row[i1 - 1], last_nids_row[i1], nids_row[i1 - 1], nids_row[i1]]
                        if None in nids:
                            continue
                        element = mesh2d.createElement(element_identifier, elementtemplate_regular)
                        element.setNodesByIdentifier(eft_regular, nids)
                        # print("Element", element_identifier, "nids", nids)
                        element_identifier += 1
                last_nids_row = nids_row
        return node_identifier, element_identifier

class EllipsoidOctantMesh:
    """
    Generates one octant of an ellipsoid, 2-D surface or full 3-D volume.
    2-D outer surface is merely added from a QuadTriangleMesh.
    3-D volume requires 3 axis-aligned and outer surfaces to each be set from a QuadTriangleMesh,
    then the interior can be built.
    """

    def __init__(self, a, b, c, element_counts, transition_element_count):
        """
        A 2-D or 3-D Octant of an ellipsoid.
        Coordinates nx are indexed in 1, 2, 3 directions from origin at index 0, 0, 0
        with holes around the corners and 3-way points.
        :param a: Axis length (radius) in x direction.
        :param b: Axis length (radius) in y direction.
        :param c: Axis length (radius) in z direction.
        :param element_counts: Number of elements across octant only in 1, 2, 3 axes.
        :param transition_element_count: Number of transition elements around outside >= 1.
        """
        assert all((count >= 2) for count in element_counts)
        assert 1 <= transition_element_count <= (min(element_counts) - 1)
        self._a = a
        self._b = b
        self._c = c
        self._element_counts = element_counts
        self._transition_element_count = transition_element_count
        self._element_count12 = element_counts[0] + element_counts[1] - 2 * transition_element_count
        self._element_count13 = element_counts[0] + element_counts[2] - 2 * transition_element_count
        self._element_count23 = element_counts[1] + element_counts[2] - 2 * transition_element_count
        # counts of elements to 3-way point opposite to 3 node at axis 1, axis 2, axis 3
        self._opp_counts = [self._element_counts[i] - self._transition_element_count for i in range(3)]
        none_parameters = [None] * 4  # x, d1, d2, d3
        self._nx = []  # shield mesh with holes over n3, n2, n1, d
        for n3 in range(element_counts[2] + 1):
            # index into transition zone
            trans3 = transition_element_count + n3 - element_counts[2]
            nx_layer = []
            for n2 in range(element_counts[1] + 1):
                # index into transition zone
                trans2 = transition_element_count + n2 - element_counts[1]
                nx_row = []
                # s = ""
                for n1 in range(element_counts[0] + 1):
                    # index into transition zone
                    trans1 = transition_element_count + n1 - element_counts[0]
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
                nx_layer.append(nx_row)
                # print(s)
            self._nx.append(nx_layer)

    def get_parameters(self):
        """
        Get parameters array e.g. for copying to ellipsoid.
        :return: Internal parameters array self._nx. Not to be modified.
        """
        return self._nx

    def set_outer_surface(self, trimesh: QuadTriangleMesh):
        """
        Set parameters of the outer surface along 1-2 edge of triangle.
        :param trimesh: Coordinates to set on outer surface.
        """
        assert trimesh.get_element_count12() == self._element_count12
        assert trimesh.get_element_count13() == self._element_count13
        assert trimesh.get_element_count23() == self._element_count23
        start_indexes = [self._element_counts[0], 0, 0]
        for n3 in range(self._opp_counts[2]):
            px, pd1, pd2, pd3 = trimesh.get_parameters12(n3)
            self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, 2, 3]], start_indexes, [[0, 1, 0], [-1, 0, 0]])
            start_indexes[2] += 1
        start_indexes = [0, 0, self._element_counts[2]]
        for n2 in range(self._opp_counts[1]):
            px, pd1, pd2, pd3 = trimesh.get_parameters31(n2, self._opp_counts[0] + 1)
            self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, 2, 3]], start_indexes, [[1, 0, 0]])
            start_indexes[1] += 1
        start_indexes = [0, self._element_counts[1], self._element_counts[2]]
        px, pd1, pd2, pd3 = trimesh.get_parameters_diagonal()
        self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, 2, 3]], start_indexes, [[1, 0, 0]])

    def _get_corner(self, indexes):
        """
        Get integer indicating which corner if any the indexes are on.
        :param indexes: Index across 1, 2, 3 directions.
        :return: 0 if not a corner, 1 if 2-3 corner, 2 if 1-3 corner, 3 if 1-2 corner, 4 if 1-2-3 corner.
        """
        # get transition index outward or negative in the regular core box.
        trans = [self._transition_element_count + indexes[i] - self._element_counts[i] for i in range(3)]
        if (trans[0] >= 0) and (trans[0] == trans[1] == trans[2]):
            return 4
        if (trans[2] < 0) and (trans[0] >= 0) and (trans[0] == trans[1]):
            return 3
        if (trans[1] < 0) and (trans[2] >= 0) and (trans[0] == trans[2]):
            return 2
        if (trans[0] < 0) and (trans[1] >= 0) and (trans[1] == trans[2]):
            return 1
        return 0

    def _set_coordinates_across(self, parameters, parameter_indexes, start_indexes, index_increments,
                                skip_start=False, skip_end=False, blend=False):
        """
        Insert parameters across the coordinates array.
        :param parameters: List of lists of N node parameters e.g. [px, pd1, pd2, pd3]
        :param parameter_indexes: Lists of parameter indexes where x=0, d1=1, d2=2, d3=3. Starts with first and
        advances at each corner, then cycles back to first. Can be negative to invert vector.
        e.g. [[0, 1, 2], [0, -2, 1]] for [x, d1, d2] then [x, -d2, d1] from first corner.
        :param start_indexes: Indexes into nx array for start point.
        :param index_increments: List of increments in indexes. Starts with first and uses next from each corner, then
        cycles back to first.
        :param skip_start: Set to True to skip the first value.
        :param skip_end: Set to True to skip the last value.
        :param blend: Set to True to blend parameters with any old parameters at locations.
        """
        indexes = start_indexes
        parameter_number = 0
        parameter_index = parameter_indexes[0]
        increment_number = 0
        index_increment = index_increments[0]
        start_n = 1 if skip_start else 0
        last_n = len(parameters[0]) - 1
        limit_n = len(parameters[0]) - (1 if skip_end else 0)
        for n in range(start_n, limit_n):
            if n > 0:
                while True:
                    indexes = [indexes[c] + index_increment[c] for c in range(3)]
                    # skip over blank transition coordinates
                    if self._nx[indexes[2]][indexes[1]][indexes[0]]:
                        break
            corner = self._get_corner(indexes)
            next_indexes = [indexes[c] + index_increment[c] for c in range(3)]
            if corner and (corner != self._get_corner(next_indexes)):
                parameter_number += 1
                if parameter_number == len(parameter_indexes):
                    parameter_number = 0
                parameter_index = parameter_indexes[parameter_number]
                increment_number += 1
                if increment_number == len(index_increments):
                    increment_number = 0
                index_increment = index_increments[increment_number]
            nx = self._nx[indexes[2]][indexes[1]][indexes[0]]
            for parameter, spix in zip(parameters, parameter_index):
                new_parameter = [-d for d in parameter[n]] if (spix < 0) else copy.copy(parameter[n])
                pix = abs(spix)
                if blend and nx[pix]:
                    if pix == 0:
                        x = [0.5 * (nx[pix][c] + new_parameter[c]) for c in range(3)]
                        nx[pix] = moveCoordinatesToEllipsoidSurface(self._a, self._b, self._c, x)
                    else:
                        # expect derivatives to be the same direction and blend magnitude with harmonic mean
                        mean_mag = 2.0 / ((1.0 / magnitude(nx[pix])) + (1.0 / magnitude(new_parameter)))
                        nx[pix] = set_magnitude(nx[pix], mean_mag)
                else:
                    nx[pix] = new_parameter
