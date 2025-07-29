"""
Utilities for building 3-D triangle-topology meshes out of quad elements
"""
import copy
from cmlibs.maths.vectorops import add, magnitude, set_magnitude
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine


class QuadTriangleMesh:
    """
    Generates a triangular mesh for interpolating with c1-continuous quad elements, with a 3-way point in the middle.
    3 corner points and 3 side points are on the triangle as follows:
             3
            / \
          13  23
         /      \
        1---12---2
    Derivatives at point 1 are d1 from 1-2, d2 from 1-3.
    Derivatives at point 2 are d1 from 1-2, d2 from 2-3.
    Derivatives at point 3 are d1 from 3-1, d2 from 3-2.
    Side point12 is elements_count2 from point1, elements_count1 from point2.
    Side point13 is elements_count3 from point1, elements_count1 from point3.
    Side point23 is elements_count3 from point2, elements_count2 from point3.
    Parameters on 12, 13 rows to 3-way point use derivatives as for point3.
    """

    def __init__(self, element_count1, element_count2, element_count3, sample_curve,
                 move_x_to_surface=None, move_d_to_surface=None):
        """
        :param element_count1: Number of elements from opposite points to side points by point1.
        :param element_count2: Number of elements from opposite points to side points by point2.
        :param element_count3: Number of elements from opposite points to side points by point3.
        :param sample_curve: Callable sampling a curve between 2 edge points.
        :param move_x_to_surface: Optional callable taking (x) and moving it to a surface.
        :param move_d_to_surface: Optional callable taking (x, d) and adjusting d to be tangential to a surface.
        """
        self._element_count1 = element_count1
        self._element_count2 = element_count2
        self._element_count3 = element_count3
        self._sample_curve = sample_curve
        self._move_x_to_surface = move_x_to_surface
        self._move_d_to_surface = move_d_to_surface
        self._element_count12 = element_count1 + element_count2
        self._element_count13 = element_count1 + element_count3
        self._element_count23 = element_count2 + element_count3
        self._node_count12 = self._element_count12 + 1
        self._node_count13 = self._element_count13 + 1
        self._node_count23 = self._element_count23 + 1
        none_parameters = [None] * 4  # x, d1, d2
        self._nx = []
        for n13 in range(self._node_count13):
            nx_row = []
            # s = ""
            for n12 in range(self._node_count12):
                if ((n12 < element_count2) or (n13 < element_count3) or
                        ((n12 - element_count2) == (n13 - element_count3))):
                    parameters = copy.copy(none_parameters)
                    # s += "[]"
                else:
                    parameters = None
                    # s += "  "
                nx_row.append(parameters)
            # print(s)
            self._nx.append(nx_row)

    def get_element_count1(self):
        return self._element_count1

    def get_element_count2(self):
        return self._element_count2

    def get_element_count3(self):
        return self._element_count3

    def get_element_count12(self):
        return self._element_count12

    def get_element_count13(self):
        return self._element_count13

    def get_element_count23(self):
        return self._element_count23

    def set_edge_parameters12(self, px, pd1, pd2, pd3=None):
        """
        Set parameters along 1-2 edge of triangle.
        :param px: Coordinates x.
        :param pd1: Derivatives in direction from point1 towards point2.
        :param pd2: Derivatives in direction towards point3.
        :param pd3: Optional derivatives in 3rd direction.
        """
        assert len(px) == self._node_count12
        for n12 in range(self._node_count12):
            d3 = copy.copy(pd3[n12]) if pd3 else None
            self._nx[0][n12] = [copy.copy(px[n12]), copy.copy(pd1[n12]), copy.copy(pd2[n12]), d3]

    def set_edge_parameters13(self, px, pd1, pd2, pd3=None):
        """
        Set parameters along 1-3 edge of triangle.
        :param px: Coordinates x.
        :param pd1: Derivatives in direction towards point2.
        :param pd2: Derivatives in direction from point1 towards point3.
        :param pd3: Optional derivatives in 3rd direction.
        """
        assert len(px) == self._node_count13
        for n13 in range(self._node_count13):
            d3 = copy.copy(pd3[n13]) if pd3 else None
            if n13 < self._element_count3:
                parameters = [copy.copy(px[n13]), copy.copy(pd1[n13]), copy.copy(pd2[n13]), d3]
            else:
                parameters = [copy.copy(px[n13]), [-d for d in pd2[n13]], copy.copy(pd1[n13]), d3]
            self._nx[n13][0] = parameters

    def set_edge_parameters23(self, px, pd1, pd2, pd3=None):
        """
        Set parameters along 123 edge of triangle.
        :param px: Coordinates x.
        :param pd1: Derivatives in direction from point2 towards point3.
        :param pd2: Derivatives in direction away from point1.
        :param pd3: Optional derivatives in 3rd direction.
        """
        assert len(px) == self._node_count23
        for n23 in range(self._node_count23):
            d3 = copy.copy(pd3[n23]) if pd3 else None
            if n23 < self._element_count3:
                parameters = [copy.copy(px[n23]), copy.copy(pd1[n23]), copy.copy(pd2[n23]), d3]
                n13 = n23
                n12 = self._element_count12
            else:
                parameters = [copy.copy(px[n23]), [-d for d in pd1[n23]], [-d for d in pd2[n23]], d3]
                n13 = self._element_count13
                n12 = (self._element_count12 if n23 == self._element_count3 else
                       self._element_count2 - (n23 - self._element_count3))
            self._nx[n13][n12] = parameters

    def get_parameters12(self, n13, count=None):
        """
        Get parameters along 1-2 direction of triangle.
        :param n13: Index from 0 at point1 toward point3, up to element_count13.
        :param count: Optional lower number to get, or None for all in row.
        :return: px, pd1, pd2, pd3
        """
        assert 0 <= n13 < self._node_count13
        if count == None:
            count = self._node_count12
        else:
            assert 1 <= count <= self._node_count12
        px = []
        pd1 = []
        pd2 = []
        pd3 = []
        nx_row = self._nx[n13]
        n12 = 0
        for n in range(count):
            while not nx_row[n12]:
                n12 += 1
            nx = nx_row[n12]
            px.append(nx[0])
            pd1.append(nx[1])
            pd2.append(nx[2])
            pd3.append(nx[3])
            n12 += 1
        return px, pd1, pd2, pd3

    def get_parameters31(self, n12, count=None):
        """
        Get parameters along 3-1 direction of triangle.
        :param n12: Index from 0 at point3 toward point2, up to element_count12.
        :param count: Optional lower number to get, or None for all in row.
        :return: px, pd1, pd2, pd3
        """
        assert 0 <= n12 < self._node_count12
        if count == None:
            count = self._node_count13
        else:
            assert 1 <= count <= self._node_count13
        px = []
        pd1 = []
        pd2 = []
        pd3 = []
        n13 = self._element_count13
        for n in range(count):
            nx_row = self._nx[n13]
            while not nx_row[n12]:
                n13 -= 1
            nx = nx_row[n12]
            px.append(nx[0])
            pd1.append(nx[1])
            pd2.append(nx[2])
            pd3.append(nx[3])
            n13 -= 1
        return px, pd1, pd2, pd3

    def get_parameters_diagonal(self):
        """
        Get parameters along diagonal from point23 to 3-way point
        :return: px, pd1, pd2, pd3
        """
        px = []
        pd1 = []
        pd2 = []
        pd3 = []
        for n in range(self._element_count1 + 1):
            n13 = self._element_count13 - n
            n12 = self._element_count12 - n
            nx = self._nx[n13][n12]
            px.append(nx[0])
            pd1.append(nx[1])
            pd2.append(nx[2])
            pd3.append(nx[3])
        return px, pd1, pd2, pd3


    def _get_corner(self, indexes):
        """
        Get integer indicating which corner if any the indexes are on.
        :param indexes: Index across 1-2 direction, index across 1-3 direction.
        :return: 0 if not a corner, 1 if 2-3 corner, 2 if 1-3 corner, 3 if 1-2 corner, 4 if 1-2-3 corner.
        """
        n12, n13 = indexes
        if n12 < self._element_count2:
            if n13 == self._element_count3:
                return 2
        elif n12 == self._element_count2:
            if n13 < self._element_count3:
                return 3
            elif n13 == self._element_count3:
                return 4
        elif (n12 - self._element_count2) == (n13 - self._element_count3):
            return 1
        return 0

    def _set_coordinates_across(self, parameters, parameter_indexes, start_indexes, index_increments,
                                skip_start=False, skip_end=False, blend=False):
        """
        Insert parameters in a row across triangle into the coordinates array.
        :param parameters: List of lists of N node parameters e.g. [px, pd1, pd2]
        :param parameter_indexes: Lists of parameter indexes where 0=x, 1=d1, 2=d2, 3=d3. Starts with first and
        advances at a 3-way line, then cycles back to first. Can be negative to invert vector.
        e.g. [[0, 1, 2], [0, -2, 1]] for [x, d1, d2] then [x, -d2, d1] from 3-way line.
        :param start_indexes: Index of first point, a list of 2-integers within the range of self._nx.
        :param index_increments: List of increments in indexes. Starts with first and moves to next at each 3-way line,
        then cycles back to first.
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
        limit_n = len(parameters[0]) - (1 if skip_end else 0)
        for n in range(start_n, limit_n):
            if n > 0:
                while True:
                    indexes = [indexes[0] + index_increment[0], indexes[1] + index_increment[1]]
                    # skip over blank coordinates around 3-way line
                    if self._nx[indexes[1]][indexes[0]]:
                        break
            corner = self._get_corner(indexes)
            if corner and (
                    corner != self._get_corner([indexes[0] + index_increment[0],  indexes[1] + index_increment[1]])):
                parameter_number += 1
                if parameter_number == len(parameter_indexes):
                    parameter_number = 0
                parameter_index = parameter_indexes[parameter_number]
                increment_number += 1
                if increment_number == len(index_increments):
                    increment_number = 0
                index_increment = index_increments[increment_number]
            nx = self._nx[indexes[1]][indexes[0]]
            for parameter, spix in zip(parameters, parameter_index):
                new_parameter = [-d for d in parameter[n]] if (spix < 0) else copy.copy(parameter[n])
                pix = abs(spix)
                if blend and nx[pix]:
                    if pix == 0:
                        x = [0.5 * (nx[pix][c] + new_parameter[c]) for c in range(3)]
                        if self._move_x_to_surface:
                            nx[pix] = self._move_x_to_surface(x)
                    else:
                        # expect derivatives to be the same direction and blend magnitude with harmonic mean
                        mean_mag = 2.0 / ((1.0 / magnitude(nx[pix])) + (1.0 / magnitude(new_parameter)))
                        nx[pix] = set_magnitude(nx[pix], mean_mag)
                else:
                    nx[pix] = new_parameter

    def _smooth_derivative_across(self, start_indexes, end_indexes, index_increments, derivative_indexes,
                                  fix_start_direction=True, fix_end_direction=True):
        """
        Smooth derivatives across triangle.
        :param start_indexes: Indexes of first point.
        :param end_indexes: Indexes of last point.
        :param index_increments: List of increments in indexes. Starts with first and after at each corner, then
        cycles back to first.
        :param derivative_indexes: List of signed derivative parameter index to along where 1=d1, 2=d2, 3=d3.
        Starts with first and advances at each corner, then cycles back to first. Can be negative to invert vector.
        e.g. [1, -2] for d1 then -d2 from first corner.
        :param fix_start_direction: Set to True to keep the start direction but scale its magnitude.
        :param fix_end_direction: Set to True to keep the end direction but scale its magnitude.
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
                    indexes = [indexes[0] + index_increment[0], indexes[1] + index_increment[1]]
                    # skip over blank coordinates around 3-way line
                    if self._nx[indexes[1]][indexes[0]]:
                        break
            corner = self._get_corner(indexes)
            if corner and (
                    corner != self._get_corner([indexes[0] + index_increment[0],  indexes[1] + index_increment[1]])):
                derivative_number += 1
                if derivative_number == len(derivative_indexes):
                    derivative_number = 0
                derivative_index = derivative_indexes[derivative_number]
                increment_number += 1
                if increment_number == len(index_increments):
                    increment_number = 0
                index_increment = index_increments[increment_number]
            parameters = self._nx[indexes[1]][indexes[0]]
            x = parameters[0]
            indexes_list.append(copy.copy(indexes))
            spix = derivative_index
            derivative_index_list.append(spix)
            pix = abs(spix)
            if parameters[pix]:
                d = [-ad for ad in parameters[pix]] if (spix < 0) else parameters[pix]
            else:
                d = [0.0, 0.0, 0.0]
            px.append(x)
            pd.append(d)
            n += 1
        sd = smoothCubicHermiteDerivativesLine(
            px, pd, fixStartDirection=fix_start_direction, fixEndDirection=fix_end_direction)
        if self._move_d_to_surface:
            for n in range(1, len(sd) - 1):
                sd[n] = self._move_d_to_surface(px[n], sd[n])
        sd = smoothCubicHermiteDerivativesLine(px, sd, fixAllDirections=True)
        for n in range(len(sd)):
            indexes = indexes_list[n]
            spix = derivative_index_list[n]
            new_derivative = [-d for d in sd[n]] if (spix < 0) else sd[n]
            pix = abs(spix)
            self._nx[indexes[1]][indexes[0]][pix] = new_derivative

    def build(self):
        """
        Determine interior coordinates from edge coordinates.
        """
        # determine 3-way point location from mean of overweighted curves in line with it
        point12 = self._nx[0][self._element_count2]
        point13 = self._nx[self._element_count3][0]
        point23 = self._nx[self._element_count13][self._element_count12]
        weight12 = self._element_count3
        weight13 = self._element_count2
        weight23 = self._element_count1
        overweighting = 1.5
        ax = self._sample_curve(
            point12[0], point12[2], None, point13[0], [-d for d in point13[2]], None, self._element_count23,
            start_weight=weight12, end_weight=weight13, overweighting=overweighting)[0][self._element_count3]
        bx = self._sample_curve(
            point12[0], point12[2], None, point23[0], [-d for d in point23[1]], None, self._element_count13,
            start_weight=weight12, end_weight=weight23, overweighting=overweighting)[0][self._element_count3]
        cx = self._sample_curve(
            point13[0], point13[2], None, point23[0], [-d for d in point23[1]], None, self._element_count12,
            start_weight=weight13, end_weight=weight23, overweighting=overweighting)[0][self._element_count2]
        x_3way = [(ax[c] + bx[c] + cx[c]) / 3.0 for c in range(3)]
        if self._move_x_to_surface:
            x_3way = self._move_x_to_surface(x_3way)

        # sample with hermite-lagrange interpolation from sides to 3-way point to get derivatives
        ax, ad1 = self._sample_curve(point23[0], point23[1], None, x_3way, None, None, self._element_count1)
        ad2 = point23[2]
        bx, bd2 = self._sample_curve(point13[0], point13[2], None, x_3way, None, None, self._element_count2)
        bd1 = point13[1]
        cx, cd2 = self._sample_curve(point12[0], point12[2], None, x_3way, None, None, self._element_count3)
        cd1 = point12[1]
        # use the minimum magnitude in all 3 directions
        d_factor = 0.6  # GRC revisit - try an exact triangle
        d_mag = d_factor * min(magnitude(ad1[-1]), magnitude(bd2[-1]), magnitude(cd2[-1]))
        d1_3way = set_magnitude(ad1[-1], d_mag)
        d2_3way = set_magnitude(bd2[-1], d_mag)

        # smooth sample from sides to 3-way points using end derivatives
        reg_count1 = self._element_count1 - 1
        reg_count2 = self._element_count2 - 1
        reg_count3 = self._element_count3 - 1
        min_weight = 2  # GRC revisit - is 1 better?
        ax, ad1, ad2 = self._sample_curve(
            ax[0], ad1[0], ad2, x_3way, d1_3way, None, self._element_count1,
            start_weight=reg_count1 + min_weight, end_weight=min_weight, end_transition=True)
        bx, bd2, bd1 = self._sample_curve(
            bx[0], bd2[0], bd1, x_3way, d2_3way, None, self._element_count2,
            start_weight=reg_count2 + min_weight, end_weight=min_weight, end_transition=True)
        cx, cd2, cd1 = self._sample_curve(
            cx[0], cd2[0], cd1, x_3way, [-d for d in add(d1_3way, d2_3way)], None, self._element_count3,
            start_weight=reg_count3 + min_weight, end_weight=min_weight, end_transition=True)
        ad2[-1] = bd2[-1]
        bd1[-1] = ad1[-1]
        self._set_coordinates_across([ax, ad1, ad2], [[0, 1, 2]], [self._element_count12, self._element_count13],
                                     [[-1, -1]])
        self._set_coordinates_across([bx, bd1, bd2], [[0, 1, 2]], [0, self._element_count3], [[1, 0]])
        self._set_coordinates_across([cx, cd1, cd2], [[0, 1, 2]], [self._element_count2, 0], [[0, 1]], skip_end=True)

        # average point coordinates across 2 directions between edges and 3-way lines.
        # 1-2 curves
        min_weight = 2  # GRC revisit
        start_indexes = [0, 0]
        corner_indexes = [self._element_count2, 0]
        end_indexes = [self._element_count12, 0]
        for i in range(1, self._element_count3):
            start_indexes[1] += 1
            corner_indexes[1] += 1
            end_indexes[1] += 1
            start = self._nx[start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[1]][end_indexes[0]]
            px, _ = self._sample_curve(
                start[0], start[1], None, corner[0], corner[1], None, self._element_count2,
                start_weight=reg_count2 + min_weight, end_weight=min_weight)
            self._set_coordinates_across(
                [px], [[0]], start_indexes, [[1, 0]], skip_start=True, skip_end=True)
            px, _ = self._sample_curve(
                corner[0], corner[1], None, end[0], end[1], None, self._element_count1,
                start_weight=min_weight, end_weight=reg_count1 + min_weight)
            self._set_coordinates_across(
                [px], [[0]], corner_indexes, [[1, 0]], skip_start=True, skip_end=True)
        # 1-3 curves
        start_indexes = [0, 0]
        corner_indexes = [0, self._element_count3]
        end_indexes = [0, self._element_count13]
        for i in range(1, self._element_count2):
            start_indexes[0] += 1
            corner_indexes[0] += 1
            end_indexes[0] += 1
            start = self._nx[start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[1]][end_indexes[0]]
            px, _ = self._sample_curve(
                start[0], start[2], None, corner[0], [-d for d in corner[1]], None, self._element_count3,
                start_weight=reg_count3 + min_weight, end_weight=min_weight)
            self._set_coordinates_across(
                [px], [[0]], start_indexes, [[0, 1]], skip_start=True, skip_end=True, blend=True)
            px, _ = self._sample_curve(
                corner[0], [-d for d in corner[1]], None, end[0], [-d for d in end[1]], None, self._element_count1,
                    start_weight=min_weight, end_weight=reg_count1 + min_weight)
            self._set_coordinates_across(
                [px], [[0]], corner_indexes, [[0, 1]], skip_start=True, skip_end=True, blend=True)
        # 2-3 curves
        start_indexes = [self._element_count12, 0]
        corner_indexes = [self._element_count12, self._element_count13]
        end_indexes = [0, self._element_count13]
        for i in range(1, self._element_count1):
            start_indexes[0] -= 1
            corner_indexes[0] -= 1
            corner_indexes[1] -= 1
            end_indexes[1] -= 1
            start = self._nx[start_indexes[1]][start_indexes[0]]
            corner = self._nx[corner_indexes[1]][corner_indexes[0]]
            end = self._nx[end_indexes[1]][end_indexes[0]]
            px, _ = self._sample_curve(
                start[0], start[2], None, corner[0], [-d for d in corner[2]], None, self._element_count3,
                start_weight=reg_count3 + min_weight, end_weight=min_weight)
            self._set_coordinates_across(
                [px], [[0]], start_indexes, [[0, 1]], skip_start=True, skip_end=True, blend=True)
            px, _ = self._sample_curve(
                corner[0], [-d for d in corner[2]], None, end[0], [-d for d in end[2]], None, self._element_count2,
                start_weight=min_weight, end_weight=reg_count2 + min_weight)
            self._set_coordinates_across(
                [px], [[0]], corner_indexes, [[-1, 0]], skip_start=True, skip_end=True, blend=True)

        # smooth 1-2 curves
        start_indexes = [0, 0]
        end_indexes = [self._element_count12, 0]
        for i in range(1, self._element_count3):
            start_indexes[1] += 1
            end_indexes[1] += 1
            self._smooth_derivative_across(start_indexes, end_indexes, [[1, 0]], [1],
                fix_start_direction=True, fix_end_direction=True)
        # smooth 1-3 curves
        start_indexes = [0, 0]
        end_indexes = [0, self._element_count13]
        for i in range(1, self._element_count2):
            start_indexes[0] += 1
            end_indexes[0] += 1
            self._smooth_derivative_across(start_indexes, end_indexes, [[0, 1]], [2, -1],
                fix_start_direction=True, fix_end_direction=True)
        # smooth 2-3 curves
        start_indexes = [self._element_count12, 0]
        end_indexes = [0, self._element_count13]
        for i in range(1, self._element_count1):
            start_indexes[0] -= 1
            end_indexes[1] -= 1
            self._smooth_derivative_across(start_indexes, end_indexes, [[0, 1], [-1, 0]], [2, -2],
                fix_start_direction=True, fix_end_direction=True)

    def calculate_d3(self, evaluate_d3):
        """
        Calls user-supplied evaluate_d3 function to set all d3 derivatives.
        Assumes all coordinates have x, d1 and d2.
        :param evaluate_d3: Callable taking x, d1, d2 and returning d3.
        """
        for n13 in range(self._node_count13):
            nx_row = self._nx[n13]
            for n12 in range(self._node_count12):
                nx = nx_row[n12]
                if nx:
                    nx[3] = evaluate_d3(nx[0], nx[1], nx[2])
