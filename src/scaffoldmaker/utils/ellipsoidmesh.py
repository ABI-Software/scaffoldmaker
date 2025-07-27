"""
Utilities for building solid ellipsoid meshes from hexahedral elements.
"""
import copy
import math
from cmlibs.maths.vectorops import add, magnitude, mult, set_magnitude, sub
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.eft_utils import determineCubicHermiteSerendipityEft, HermiteNodeLayoutManager
from scaffoldmaker.utils.geometry import (
    getEllipsePointAtTrueAngle, getEllipseTangentAtPoint, moveCoordinatesToEllipsoidSurface,
    moveDerivativeToEllipsoidSurface, sampleCurveOnEllipsoid)
from scaffoldmaker.utils.interpolation import (
    computeCubicHermiteArcLength, interpolateCubicHermite, smoothCubicHermiteDerivativesLine)


class EllipsoidMesh:
    """
    Generates a solid ellipsoid of hexahedral elements with oblique cross axes suited to describing lung geometry.
    """

    def __init__(self, a, b, c, element_counts, transition_element_count,
                 axis2_x_rotation_radians, axis3_x_rotation_radians, surface_only=False):
        """

        :param element_counts:
        :param transition_element_count:
        :param a: Axis length (radius) in x direction.
        :param b: Axis length (radius) in y direction.
        :param c: Axis length (radius) in z direction.
        :param axis2_x_rotation_radians: Rotation of axis 2 about +x direction
        :param axis3_x_rotation_radians: Rotation of axis 3 about +x direction.
        :param surface_only: Set to True to only make nodes and 2-D elements on the surface.
        """
        assert all((count >= 4) and (count % 2 == 0) for count in element_counts)
        assert 1 <= transition_element_count <= (min(element_counts) // 2 - 1)
        self._element_counts = element_counts
        self._transition_element_count = transition_element_count
        self._a = a
        self._b = b
        self._c = c
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

        # get outside curve in 1-2 plane starting in 1 = x direction, 1st quadrant of 1-2 loop
        start_x = [self._a, 0.0, 0.0]
        start_d1 = [0.0, math.cos(self._axis2_x_rotation_radians), math.sin(self._axis2_x_rotation_radians)]
        start_d2 = [0.0, math.cos(self._axis3_x_rotation_radians), math.sin(self._axis3_x_rotation_radians)]
        end_x = [0.0] + getEllipsePointAtTrueAngle(self._b, self._c, self._axis2_x_rotation_radians)
        end_d1 = [-1.0, 0.0, 0.0]
        end_d2 = [0.0] + getEllipseTangentAtPoint(self._b, self._c, end_x[1:])
        px, pd1, pd2 = sampleCurveOnEllipsoid(self._a, self._b, self._c,
                                              start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count_q12)
        self._set_coordinates_around_rim(
            [px, pd1, pd2], [[0, 1, 2]],
            [self._element_counts[0], half_counts[1], half_counts[2]],
            [[0, 1, 0], [-1, 0, 0]])
        # get outside curve in 1-3 plane starting in 1 = x direction, 1st quadrant of 1-3 loop
        start_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        start_x, start_d1, start_d2 = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]][:3]
        end_x = [0.0] + getEllipsePointAtTrueAngle(self._b, self._c, self._axis3_x_rotation_radians)
        end_d2 = [-1.0, 0.0, 0.0]
        end_d1 = [0.0] + [-d for d in getEllipseTangentAtPoint(self._b, self._c, end_x[1:])]
        px, pd2, pd1 = sampleCurveOnEllipsoid(self._a, self._b, self._c,
                                              start_x, start_d2, start_d1, end_x, end_d2, end_d1, elements_count_q13)
        self._set_coordinates_around_rim(
            [px, pd1, pd2], [[0, 1, 2], [0, 2, -1]], start_indexes,
            [[0, 0, 1], [-1, 0, 0]])

        # get outside curve in 2-3 plane starting in 2 direction, 1st quadrant of 2-3 loop
        start_indexes = [half_counts[0], self._element_counts[1], half_counts[2]]
        start_x, start_d1, start_d2 = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]][:3]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        end_x, end_d1, end_d2 = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]][:3]
        end_d1 = [-d for d in end_d1]
        end_d2 = [-d for d in end_d2]
        px, pd2, pd1 = sampleCurveOnEllipsoid(self._a, self._b, self._c,
                                              start_x, start_d2, start_d1, end_x, end_d2, end_d1, elements_count_q23)
        self._set_coordinates_around_rim(
            [px, pd1, pd2], [[0, 1, 2], [0, -1, -2]], start_indexes,
            [[0, 0, 1], [0, -1, 0]])

        self._sample_surface_octant1()

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
                x = [x[0], -x[1], -x[2]]
                d1 = [-d1[0], d1[1], d1[2]]
                d2 = [-d2[0], d2[1], d2[2]]
                self._nx[n3][self._element_counts[1] - n2][n1] = [x, d1, d2, None]

        # sample 2nd quadrant of 1-3 loop and blend end derivatives
        start_indexes = [half_counts[0], 0, half_counts[2]]
        start_x, start_d1, start_d2 = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]][:3]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        end_x, end_d1, end_d2 = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]][:3]
        px, pd2, pd1 = sampleCurveOnEllipsoid(self._a, self._b, self._c,
                                              start_x, start_d2, start_d1, end_x, end_d2, end_d1, elements_count_q23)
        self._set_coordinates_around_rim(
            [px, pd1, pd2], [[0, 1, 2]], start_indexes,
            [[0, 0, 1], [0, 1, 0]], blend_start=True, blend_end=True)

        self._sample_surface_octant2()

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

    def _sample_surface_octant1(self):
        """
        Sample surface octant between positive 1, positive 2 and positive 3 axes
        """
        half_counts = [count // 2 for count in self._element_counts]
        regular_row_counts = [half_counts[i] - self._transition_element_count - 1 for i in range(3)]
        elements_count_q12 = half_counts[0] + half_counts[1] - 2 * self._transition_element_count
        elements_count_q13 = half_counts[0] + half_counts[2] - 2 * self._transition_element_count
        elements_count_q23 = half_counts[1] + half_counts[2] - 2 * self._transition_element_count

        # determine 3-way point location from mean of overweighted curves in line with it
        point1 = self._nx[self._element_counts[2]][self._element_counts[1]][half_counts[0]]
        point2 = self._nx[self._element_counts[2]][half_counts[1]][self._element_counts[0]]
        point3 = self._nx[half_counts[2]][self._element_counts[1]][self._element_counts[0]]
        min_weight = 1
        weight = [regular_row_counts[i] + min_weight for i in range(3)]
        weight[0] /= magnitude(point1[0])
        weight[1] /= magnitude(point2[0])
        weight[2] /= magnitude(point3[0])
        offset = 1
        overweighting = 1.5
        ax = sampleCurveOnEllipsoid(
            self._a, self._b, self._c, point1[0], point1[1], None, point2[0], [-d for d in point2[2]], None,
            elements_count_q12, start_weight=weight[0], end_weight=weight[1],
            overweighting=overweighting)[0][regular_row_counts[0] + offset]
        bx = sampleCurveOnEllipsoid(
            self._a, self._b, self._c, point3[0], point3[2], None, point1[0], [-d for d in point1[1]], None,
            elements_count_q13, start_weight=weight[2], end_weight=weight[0],
            overweighting=overweighting)[0][regular_row_counts[2] + offset]
        cx = sampleCurveOnEllipsoid(
            self._a, self._b, self._c, point3[0], point3[2], None, point2[0], [-d for d in point2[2]], None,
            elements_count_q23, start_weight=weight[2], end_weight=weight[1],
            overweighting=overweighting)[0][regular_row_counts[2] + offset]
        x_3way = moveCoordinatesToEllipsoidSurface(
            self._a, self._b, self._c, [(ax[c] + bx[c] + cx[c]) / 3.0 for c in range(3)])

        # sample with hermite-lagrange interpolation from sides to 3-way point derivatives
        a_start_indexes = [half_counts[0], self._element_counts[1], self._element_counts[2]]
        start = self._nx[a_start_indexes[2]][a_start_indexes[1]][a_start_indexes[0]]
        ax, ad1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            start[0], start[1], None, x_3way, None, None, regular_row_counts[0] + 1)
        ad2 = start[2]
        b_start_indexes = [self._element_counts[0], half_counts[1], self._element_counts[2]]
        start = self._nx[b_start_indexes[2]][b_start_indexes[1]][b_start_indexes[0]]
        bx, bd2 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            start[0], start[2], None, x_3way, None, None, regular_row_counts[1] + 1)
        bd1 = start[1]
        c_start_indexes = [self._element_counts[0], self._element_counts[1], half_counts[2]]
        start = self._nx[c_start_indexes[2]][c_start_indexes[1]][c_start_indexes[0]]
        cx, cd2 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            start[0], start[2], None, x_3way, None, None, regular_row_counts[2] + 1)
        cd1 = start[1]
        # use the minimum magnitude in all 3 directions
        d_factor = 0.6
        d_mag = d_factor * min(magnitude(ad1[-1]), magnitude(bd2[-1]), magnitude(cd2[-1]))
        d1_3way = set_magnitude(ad1[-1], d_mag)
        d2_3way = set_magnitude(bd2[-1], d_mag)

        # smooth sample from sides to 3-way points using end derivatives
        min_weight = 2
        ax, ad1, ad2 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            ax[0], ad1[0], ad2, x_3way, d1_3way, None, regular_row_counts[0] + 1,
            start_weight=regular_row_counts[0] + min_weight, end_weight=min_weight, end_transition=True)
        bx, bd2, bd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            bx[0], bd2[0], bd1, x_3way, d2_3way, None, regular_row_counts[1] + 1,
            start_weight=regular_row_counts[1] + min_weight, end_weight=min_weight, end_transition=True)
        cx, cd2, cd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            cx[0], cd2[0], cd1, x_3way, [-d for d in add(d1_3way, d2_3way)], None, regular_row_counts[2] + 1,
            start_weight=regular_row_counts[2] + min_weight, end_weight=min_weight, end_transition=True)
        ad2[-1] = bd2[-1]
        bd1[-1] = ad1[-1]
        self._set_coordinates_around_rim([ax, ad1, ad2], [[0, 1, 2]], a_start_indexes, [[1, 0, 0]])
        self._set_coordinates_around_rim([bx, bd1, bd2], [[0, 1, 2]], b_start_indexes, [[0, 1, 0]])
        self._set_coordinates_around_rim([cx, cd1, cd2], [[0, 1, 2]], c_start_indexes, [[0, 0, 1]], skip_end=True)

        # 1-2 curve in 3-direction
        min_weight = 2
        start_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        corner_indexes = [self._element_counts[0], self._element_counts[1], half_counts[2]]
        end_indexes = [half_counts[0], self._element_counts[1], half_counts[2]]
        for i in range(regular_row_counts[2]):
            start_indexes[2] += 1
            corner_indexes[2] += 1
            end_indexes[2] += 1
            start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                start[0], start[1], None, corner[0], corner[1], None, regular_row_counts[1] + 1,
                start_weight=regular_row_counts[0] + min_weight, end_weight=min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], start_indexes, [[0, 1, 0]], skip_start=True, skip_end=True)
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                corner[0], corner[1], None, end[0], end[1], None, regular_row_counts[0] + 1,
                start_weight=min_weight, end_weight=regular_row_counts[1] + min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], corner_indexes, [[-1, 0, 0]], skip_start=True, skip_end=True)
        # 1-3 curve in 2-direction
        start_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        corner_indexes = [self._element_counts[0], half_counts[1], self._element_counts[2]]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        for i in range(regular_row_counts[1]):
            start_indexes[1] += 1
            corner_indexes[1] += 1
            end_indexes[1] += 1
            start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                start[0], start[2], None, corner[0], [-d for d in corner[1]], None, regular_row_counts[2] + 1,
                start_weight=regular_row_counts[2] + min_weight, end_weight=min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], start_indexes, [[0, 0, 1]], skip_start=True, skip_end=True, blend_middle=True)
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                corner[0], [-d for d in corner[1]], None, end[0], [-d for d in end[1]], None, regular_row_counts[0] + 1,
                    start_weight=min_weight, end_weight=regular_row_counts[0] + min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], corner_indexes, [[-1, 0, 0]], skip_start=True, skip_end=True, blend_middle=True)
        # 2-3 curve in 1-direction
        start_indexes = [half_counts[0], self._element_counts[1], half_counts[2]]
        corner_indexes = [half_counts[0], self._element_counts[1], self._element_counts[2]]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        for i in range(regular_row_counts[0]):
            start_indexes[0] += 1
            corner_indexes[0] += 1
            end_indexes[0] += 1
            start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                start[0], start[2], None, corner[0], [-d for d in corner[2]], None, regular_row_counts[2] + 1,
                start_weight=regular_row_counts[2] + min_weight, end_weight=min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], start_indexes, [[0, 0, 1]], skip_start=True, skip_end=True, blend_middle=True)
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                corner[0], [-d for d in corner[2]], None, end[0], [-d for d in end[2]], None, regular_row_counts[1] + 1,
                start_weight=min_weight, end_weight=regular_row_counts[1] + min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], corner_indexes, [[0, -1, 0]], skip_start=True, skip_end=True, blend_middle=True)

        # smooth 1-2 curves
        start_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        end_indexes = [half_counts[0], self._element_counts[1], half_counts[2]]
        for i in range(regular_row_counts[2]):
            start_indexes[2] += 1
            end_indexes[2] += 1
            self._smooth_derivatives_around_rim(start_indexes, end_indexes, [[0, 1, 0], [-1, 0, 0]], [1, 1], [1],
                fix_start_direction=True, fix_end_direction=True)
        # smooth 1-3 curves
        start_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        for i in range(regular_row_counts[1]):
            start_indexes[1] += 1
            end_indexes[1] += 1
            self._smooth_derivatives_around_rim(start_indexes, end_indexes, [[0, 0, 1], [-1, 0, 0]], [2, -1], [-1],
                fix_start_direction=True, fix_end_direction=True)
        # smooth 2-3 curves
        start_indexes = [half_counts[0], self._element_counts[1], half_counts[2]]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        for i in range(regular_row_counts[0]):
            start_indexes[0] += 1
            end_indexes[0] += 1
            self._smooth_derivatives_around_rim(start_indexes, end_indexes, [[0, 0, 1], [0, -1, 0]], [2, -2], [-2],
                fix_start_direction=True, fix_end_direction=True)

    def _sample_surface_octant2(self):
        """
        Sample surface octant between positive 1, negative 2 and positive 3 axes
        """
        half_counts = [count // 2 for count in self._element_counts]
        regular_row_counts = [half_counts[i] - self._transition_element_count - 1 for i in range(3)]
        elements_count_q12 = half_counts[0] + half_counts[1] - 2 * self._transition_element_count
        elements_count_q13 = half_counts[0] + half_counts[2] - 2 * self._transition_element_count
        elements_count_q23 = half_counts[1] + half_counts[2] - 2 * self._transition_element_count

        # determine 3-way point location from mean of overweighted curves in line with it
        point1 = self._nx[self._element_counts[2]][0][half_counts[0]]
        point2 = self._nx[self._element_counts[2]][half_counts[1]][self._element_counts[0]]
        point3 = self._nx[half_counts[2]][0][self._element_counts[0]]
        min_weight = 1
        weight = [regular_row_counts[i] + min_weight for i in range(3)]
        weight[0] /= magnitude(point1[0])
        weight[1] /= magnitude(point2[0])
        weight[2] /= magnitude(point3[0])
        offset = 1
        overweighting = 1.5
        ax = sampleCurveOnEllipsoid(
            self._a, self._b, self._c, point1[0], point1[1], None, point2[0], point2[2], None,
            elements_count_q12, start_weight=weight[0], end_weight=weight[1],
            overweighting=overweighting)[0][regular_row_counts[0] + offset]
        bx = sampleCurveOnEllipsoid(
            self._a, self._b, self._c, point3[0], point3[2], None, point1[0], [-d for d in point1[1]], None,
            elements_count_q13, start_weight=weight[2], end_weight=weight[0],
            overweighting=overweighting)[0][regular_row_counts[2] + offset]
        cx = sampleCurveOnEllipsoid(
            self._a, self._b, self._c, point3[0], point3[2], None, point2[0], point2[2], None,
            elements_count_q23, start_weight=weight[2], end_weight=weight[1],
            overweighting=overweighting)[0][regular_row_counts[2] + offset]
        x_3way = moveCoordinatesToEllipsoidSurface(
            self._a, self._b, self._c, [(ax[c] + bx[c] + cx[c]) / 3.0 for c in range(3)])
        # visualise ax, bx, cx:
        # self._nx[0][self._element_counts[1]][0][0] = ax
        # self._nx[0][self._element_counts[1]][self._element_counts[0]][0] = bx
        # self._nx[self._element_counts[2]][0][0][0] = cx

        # sample with hermite-lagrange interpolation from sides to 3-way point derivatives
        a_start_indexes = [half_counts[0], 0, self._element_counts[2]]
        start = self._nx[a_start_indexes[2]][a_start_indexes[1]][a_start_indexes[0]]
        ax, ad1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            start[0], start[1], None, x_3way, None, None, regular_row_counts[0] + 1)
        ad2 = start[2]
        b_start_indexes = [self._element_counts[0], half_counts[1], self._element_counts[2]]
        start = self._nx[b_start_indexes[2]][b_start_indexes[1]][b_start_indexes[0]]
        bx, bd2 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            start[0], [-d for d in start[2]], None, x_3way, None, None, regular_row_counts[1] + 1)
        bd1 = start[1]
        c_start_indexes = [self._element_counts[0], 0, half_counts[2]]
        start = self._nx[c_start_indexes[2]][c_start_indexes[1]][c_start_indexes[0]]
        cx, cd2 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            start[0], start[2], None, x_3way, None, None, regular_row_counts[2] + 1)
        cd1 = start[1]
        # use the minimum magnitude in all 3 directions
        d_factor = 0.6
        d_mag = d_factor * min(magnitude(ad1[-1]), magnitude(bd2[-1]), magnitude(cd2[-1]))
        d1_3way = set_magnitude(ad1[-1], d_mag)
        d2_3way = set_magnitude([-d for d in bd2[-1]], d_mag)

        # smooth sample from sides to 3-way points using end derivatives
        min_weight = 2
        ax, ad1, ad2 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            ax[0], ad1[0], ad2, x_3way, d1_3way, None, regular_row_counts[0] + 1,
            start_weight=regular_row_counts[0] + min_weight, end_weight=min_weight, end_transition=True)
        bx, bd2, bd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            bx[0], bd2[0], bd1, x_3way, [-d for d in d2_3way], None, regular_row_counts[1] + 1,
            start_weight=regular_row_counts[1] + min_weight, end_weight=min_weight, end_transition=True)
        cx, cd2, cd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            cx[0], cd2[0], cd1, x_3way, sub(d2_3way, d1_3way), None, regular_row_counts[2] + 1,
            start_weight=regular_row_counts[2] + min_weight, end_weight=min_weight, end_transition=True)
        ad2[-1] = bd2[-1]
        bd1[-1] = ad1[-1]
        self._set_coordinates_around_rim([ax, ad1, ad2], [[0, 1, 2]], a_start_indexes, [[1, 0, 0]])
        self._set_coordinates_around_rim([bx, bd1, bd2], [[0, 1, -2]], b_start_indexes, [[0, -1, 0]], blend_start=True)
        self._set_coordinates_around_rim([cx, cd1, cd2], [[0, 1, 2]], c_start_indexes, [[0, 0, 1]],
                                         blend_start=True, skip_end=True)

        # 1-2 curve in 3-direction
        min_weight = 2
        start_indexes = [half_counts[0], 0, half_counts[2]]
        corner_indexes = [self._element_counts[0], 0, half_counts[2]]
        end_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        for i in range(regular_row_counts[2]):
            start_indexes[2] += 1
            corner_indexes[2] += 1
            end_indexes[2] += 1
            start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                start[0], start[1], None, corner[0], corner[1], None, regular_row_counts[0] + 1,
                start_weight=regular_row_counts[0] + min_weight, end_weight=min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], start_indexes, [[1, 0, 0]], skip_start=True, skip_end=True)
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                corner[0], corner[1], None, end[0], end[1], None, regular_row_counts[1] + 1,
                start_weight=min_weight, end_weight=regular_row_counts[1] + min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], corner_indexes, [[0, 1, 0]], skip_start=True, skip_end=True)
        # 1-3 curve in 2-direction
        start_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        corner_indexes = [self._element_counts[0], half_counts[1], self._element_counts[2]]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        for i in range(regular_row_counts[1]):
            start_indexes[1] -= 1
            corner_indexes[1] -= 1
            end_indexes[1] -= 1
            start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                start[0], start[2], None, corner[0], [-d for d in corner[1]], None, regular_row_counts[2] + 1,
                start_weight=regular_row_counts[2] + min_weight, end_weight=min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], start_indexes, [[0, 0, 1]], skip_start=True, skip_end=True, blend_middle=True)
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                corner[0], [-d for d in corner[1]], None, end[0], [-d for d in end[1]], None, regular_row_counts[0] + 1,
                    start_weight=min_weight, end_weight=regular_row_counts[0] + min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], corner_indexes, [[-1, 0, 0]], skip_start=True, skip_end=True, blend_middle=True)
        # 2-3 curve in 1-direction
        start_indexes = [half_counts[0], 0, half_counts[2]]
        corner_indexes = [half_counts[0], 0, self._element_counts[2]]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        for i in range(regular_row_counts[0]):
            start_indexes[0] += 1
            corner_indexes[0] += 1
            end_indexes[0] += 1
            start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                start[0], start[2], None, corner[0], corner[2], None, regular_row_counts[2] + 1,
                start_weight=regular_row_counts[2] + min_weight, end_weight=min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], start_indexes, [[0, 0, 1]], skip_start=True, skip_end=True, blend_middle=True)
            px, _ = sampleCurveOnEllipsoid(
                self._a, self._b, self._c,
                corner[0], corner[2], None, end[0], end[2], None, regular_row_counts[1] + 1,
                start_weight=min_weight, end_weight=regular_row_counts[1] + min_weight)
            self._set_coordinates_around_rim(
                [px], [[0]], corner_indexes, [[0, 1, 0]], skip_start=True, skip_end=True, blend_middle=True)

        # smooth 1-2 curves
        start_indexes = [half_counts[0], 0, half_counts[2]]
        end_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        for i in range(regular_row_counts[2]):
            start_indexes[2] += 1
            end_indexes[2] += 1
            self._smooth_derivatives_around_rim(start_indexes, end_indexes, [[1, 0, 0], [0, 1, 0]], [1, 1], [1],
                fix_start_direction=True, fix_end_direction=True, blend_start=True, blend_end=True)
        # smooth 1-3 curves
        start_indexes = [self._element_counts[0], half_counts[1], half_counts[2]]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        for i in range(regular_row_counts[1]):
            start_indexes[1] -= 1
            end_indexes[1] -= 1
            self._smooth_derivatives_around_rim(start_indexes, end_indexes, [[0, 0, 1], [-1, 0, 0]], [2, -1], [-1],
                fix_start_direction=True, fix_end_direction=True, blend_start=True, blend_end=True)
        # smooth 2-3 curves
        start_indexes = [half_counts[0], 0, half_counts[2]]
        end_indexes = [half_counts[0], half_counts[1], self._element_counts[2]]
        for i in range(regular_row_counts[0]):
            start_indexes[0] += 1
            end_indexes[0] += 1
            self._smooth_derivatives_around_rim(start_indexes, end_indexes, [[0, 0, 1], [0, 1, 0]], [2], [2],
                fix_start_direction=True, fix_end_direction=True, blend_start=True, blend_end=True)

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
        half_counts = [count // 2 for count in self._element_counts]
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

    def generate_mesh(self, fieldmodule, coordinates):
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

        node_layout_manager = HermiteNodeLayoutManager()

        if self._surface_only:
            mesh2d = fieldmodule.findMeshByDimension(2)
            element_identifier = 1
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
