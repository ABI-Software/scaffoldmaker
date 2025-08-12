"""
Utilities for building solid ellipsoid meshes from hexahedral elements.
"""
from cmlibs.maths.vectorops import add, cross, div, magnitude, mult, set_magnitude, sub
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node

from scaffoldmaker.utils.eft_utils import determineCubicHermiteSerendipityEft, HermiteNodeLayoutManager
from scaffoldmaker.utils.geometry import (
    getEllipsePointAtTrueAngle, getEllipseTangentAtPoint, moveCoordinatesToEllipsoidSurface,
    moveDerivativeToEllipsoidSurface, sampleCurveOnEllipsoid)
from scaffoldmaker.utils.interpolation import (
    DerivativeScalingMode, get_nway_point, linearlyInterpolateVectors, sampleHermiteCurve,
    smoothCubicHermiteDerivativesLine)
from scaffoldmaker.utils.quadtrianglemesh import QuadTriangleMesh
import copy
from enum import Enum
import math


class EllipsoidSurfaceD3Mode(Enum):
    SURFACE_NORMAL = 1  # surface D3 are exact surface normals to ellipsoid
    OBLIQUE_DIRECTION = 2  # surface D3 are in direction of surface point on ellipsoid, gives flat oblique planes


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
        self._trans_count = transition_element_count
        self._axis2_x_rotation_radians = axis2_x_rotation_radians
        self._axis3_x_rotation_radians = axis3_x_rotation_radians
        self._surface_only = surface_only
        self._nway_d_factor = 0.6
        self._surface_d3_mode = EllipsoidSurfaceD3Mode.SURFACE_NORMAL
        self._box_group = None
        self._transition_group = None
        self._octant_group_lists = None
        none_parameters = [None] * 4  # x, d1, d2, d3
        self._nx = []  # shield mesh with holes over n3, n2, n1, d
        self._nids = []
        half_counts = [count // 2 for count in self._element_counts]
        for n3 in range(self._element_counts[2] + 1):
            # index into transition zone
            trans3 = (self._trans_count - n3) if (n3 < half_counts[2]) else \
                (self._trans_count + n3 - self._element_counts[2])
            nx_layer = []
            nids_layer = []
            # print(n3, trans3)
            for n2 in range(self._element_counts[1] + 1):
                # index into transition zone
                trans2 = (self._trans_count - n2) if (n2 < half_counts[1]) else \
                    (self._trans_count + n2 - self._element_counts[1])
                nx_row = []
                nids_row = []
                # s = ""
                for n1 in range(self._element_counts[0] + 1):
                    # index into transition zone
                    trans1 = (self._trans_count - n1) if (n1 < half_counts[0]) else \
                        (self._trans_count + n1 - self._element_counts[0])
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

    def set_box_transition_groups(self, box_group, transition_group):
        """
        Set zinc groups to fill with elements in box and transition regions, if not surface_only.
        :param box_group: Group field to add elements from box region to.
        :param transition_group: Group field to add elements from transition region to.
        """
        self._box_group = box_group
        self._transition_group = transition_group

    def set_octant_group_lists(self, octant_group_lists):
        """
        Set lists of zinc groups to add elements in each of the 8 octants to.
        :param octant_group_lists: List of 8 lists of group fields to put elements into. For example, the
        elements of octant 0 (negative 3 axis, negative 2 axis, negative 1 axis) will be put in the mesh group of
        the appropriate dimension for each group in octant_groups_lists[0]. Order of octants for N (negative) or
        P (positive) 321 axes: NNN, NNP, NPN, NPP, PNN, PNP, PPN, PPP.
        """
        assert (octant_group_lists is None) or (len(octant_group_lists) == 8)
        self._octant_group_lists = octant_group_lists

    def set_nway_derivative_factor(self, nway_derivative_factor):
        """
        Set factor controlling the shape of 3-way and 4-way points in ellipsoid mesh.
        :param nway_d_factor: Value, normally from 0.5 to 1.0 giving n-way derivative magnitude as a proportion
        of the minimum regular magnitude sampled to the n-way point. This reflects that distances from the mid-side
        of a triangle to the centre are shorter, so the derivative in the middle must be smaller.
        """
        self._nway_d_factor = nway_derivative_factor

    def set_surface_d3_mode(self, surface_d3_mode: EllipsoidSurfaceD3Mode):
        """
        Set mode controlling how surface d3 values are calculated.
        :param surface_d3_mode: Value from EllipsoidSurfaceD3Mode.
        """
        self._surface_d3_mode = surface_d3_mode

    def build(self):
        """
        Determine coordinates and derivatives over and within ellipsoid.
        """
        half_counts = [count // 2 for count in self._element_counts]
        box_counts = [half_counts[i] - self._trans_count for i in range(3)]

        octant1 = self._build_ellipsoid_octant(half_counts, self._axis2_x_rotation_radians, self._axis3_x_rotation_radians)

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

        octant2 = self._build_ellipsoid_octant([half_counts[0], half_counts[2], half_counts[1]],
                                               self._axis3_x_rotation_radians - math.pi, self._axis2_x_rotation_radians)

        # transfer parameters on bottom plane of octant1 to target location of octant2 for blending
        n3 = half_counts[2]
        for i2 in range(1, half_counts[1] + 1):
            n2 = half_counts[1] + i2
            for i1 in range(half_counts[0] + 1):
                n1 = half_counts[0] + i1
                box = (i1 <= box_counts[0]) and (i2 <= box_counts[1])
                parameters = self._nx[n3][n2][n1]
                if not parameters or not parameters[0]:
                    continue
                x, d1, d2, d3 = parameters
                x = [x[0], -x[1], -x[2]]
                d1 = [d1[0], -d1[1], -d1[2]] if box else [-d1[0], d1[1], d1[2]]
                d2 = [-d2[0], d2[1], d2[2]]
                d3 = ([-d3[0], d3[1], d3[2]] if box else [d3[0], -d3[1], -d3[2]]) if d3 else None
                self._nx[n3][self._element_counts[1] - n2][n1] = [x, d1, d2, d3]

        # copy and mirror in y and z octant2 into ellipsoid, blending existing derivatives
        octant_parameters = octant2.get_parameters()
        for o3 in range(half_counts[1] + 1):
            octant_nx_layer = octant_parameters[o3]
            for o2 in range(half_counts[2] + 1):
                octant_nx_row = octant_nx_layer[o2]
                nx_row = self._nx[half_counts[2] + o2][half_counts[1] - o3]
                box_row = (o2 <= box_counts[2]) and (o3 <= box_counts[1])
                top_transition = o2 > box_counts[2]
                for o1 in range(half_counts[0] + 1):
                    octant_nx = octant_nx_row[o1]
                    if octant_nx and octant_nx[0]:
                        box = box_row and (o1 <= box_counts[0])
                        x, d1, d2, d3 = octant_nx
                        x = [x[0], -x[1], -x[2]]
                        if top_transition:
                            if o3 > box_counts[1]:
                                d1 = [d1[0], -d1[1], -d1[2]]
                                d2 = [d2[0], -d2[1], -d2[2]]
                                # fix 3-way point case on transition 'corner':
                                if (o3 >= box_counts[1]) and (o1 >= box_counts[0]):
                                    d2 = add(d1, d2)
                            else:
                                d1 = [-d1[0], d1[1], d1[2]]
                                d2 = [-d2[0], d2[1], d2[2]]
                            d3 = [d3[0], -d3[1], -d3[2]] if d3 else None
                        elif box:
                            d1 = [d1[0], -d1[1], -d1[2]]
                            d2, d3 = [-d3[0], d3[1], d3[2]], [d2[0], -d2[1], -d2[2]]
                        else:
                            if o3 > box_counts[1]:
                                d1 = [d1[0], -d1[1], -d1[2]]
                                d2 = [d2[0], -d2[1], -d2[2]]
                            else:
                                d1, d2 = [-d2[0], d2[1], d2[2]], [d1[0], -d1[1], -d1[2]]
                            d3 = [d3[0], -d3[1], -d3[2]] if d3 else None
                        new_nx = [x, d1, d2, d3]
                        nx = nx_row[half_counts[0] + o1]
                        # expect coordinates to be at the same location on boundaries
                        nx[0] = copy.copy(x)
                        for pix in range(1, 4):
                            d = new_nx[pix]
                            if nx[pix] and d:
                                # blend derivatives with harmonic mean magnitude; should already be in same direction
                                d = linearlyInterpolateVectors(
                                    nx[pix], d, 0.5, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
                            nx[pix] = copy.copy(d)

        # transfer blended d2 derivatives on bottom plane of octant2 back to octant1
        n3 = half_counts[2]
        for i2 in range(1, half_counts[1] + 1):
            n2 = half_counts[1] + i2
            for i1 in range(half_counts[0] + 1):
                n1 = half_counts[0] + i1
                box = (i1 <= box_counts[0]) and (i2 <= box_counts[1])
                parameters = self._nx[n3][self._element_counts[1] - n2][n1]
                if not parameters or not parameters[0]:
                    continue
                x, d1, d2, d3 = parameters
                x = [x[0], -x[1], -x[2]]
                d1 = [d1[0], -d1[1], -d1[2]] if box else [-d1[0], d1[1], d1[2]]
                d2 = [-d2[0], d2[1], d2[2]]
                d3 = ([-d3[0], d3[1], d3[2]] if box else [d3[0], -d3[1], -d3[2]]) if d3 else None
                self._nx[n3][n2][n1] = [x, d1, d2, d3]

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
        for n3 in range(half_counts[2] + 1, self._element_counts[2] + 1):
            for n2 in range(0, self._element_counts[1] + 1):
                top_nx_row = self._nx[n3][n2]
                box2 = (self._trans_count <= n2 <=
                        (self._element_counts[1] - self._trans_count))
                box_row = (n3 <= (half_counts[2] + box_counts[2])) and box2
                bottom_nx_row = self._nx[self._element_counts[2] - n3][self._element_counts[1] - n2]
                bottom_transition = box2 and (n3 >= (half_counts[2] + box_counts[2]))
                for n1 in range(self._element_counts[0] + 1):
                    parameters = top_nx_row[n1]
                    if parameters and parameters[0]:
                        box = box_row and (self._trans_count <= n1 <=
                                           (self._element_counts[0] - self._trans_count))
                        x, d1, d2, d3 = parameters
                        x = [x[0], -x[1], -x[2]]
                        d2 = [-d2[0], d2[1], d2[2]] if d2 else None
                        if box and not bottom_transition:
                            d1 = [d1[0], -d1[1], -d1[2]]
                            d3 = [-d3[0], d3[1], d3[2]]
                        else:
                            d1 = [-d1[0], d1[1], d1[2]] if d1 else None
                            d3 = [d3[0], -d3[1], -d3[2]] if d3 else None
                        bottom_nx_row[n1] = [x, d1, d2, d3]

    def _build_ellipsoid_octant(self, half_counts, axis2_x_rotation_radians, axis3_x_rotation_radians):
        """
        Get coordinates of top, right, front octant with supplied angles.
        :param half_counts: Numbers of elements across octant 1, 2 and 3 directions.
        :param axis2_x_rotation_radians: Rotation of axis 2 about +x direction
        :param axis3_x_rotation_radians: Rotation of axis 3 about +x direction.
        """
        elements_count_q12 = half_counts[0] + half_counts[1] - 2 * self._trans_count
        elements_count_q13 = half_counts[0] + half_counts[2] - 2 * self._trans_count
        elements_count_q23 = half_counts[1] + half_counts[2] - 2 * self._trans_count
        box_counts = [half_counts[i] - self._trans_count for i in range(3)]

        origin = [0.0, 0.0, 0.0]
        axis1 = [self._a, 0.0, 0.0]
        axis2 = [0.0] + getEllipsePointAtTrueAngle(self._b, self._c, axis2_x_rotation_radians)
        axis3 = [0.0] + getEllipsePointAtTrueAngle(self._b, self._c, axis3_x_rotation_radians)
        axis_d1 = div(axis1, half_counts[0])
        axis_d2 = div(axis2, half_counts[1])
        axis_d3 = div(axis3, half_counts[2])
        axis_md1 = [-d for d in axis_d1]
        axis_md2 = [-d for d in axis_d2]
        axis_md3 = [-d for d in axis_d3]
        # magnitude of most derivatives indicating only direction so magnitude not known
        dir_mag = min(magnitude(axis_d1), magnitude(axis_d2), magnitude(axis_d3))
        axis2_dt = set_magnitude([0.0] + getEllipseTangentAtPoint(self._b, self._c, axis2[1:]), magnitude(axis_d3))
        axis3_dt = set_magnitude([0.0] + getEllipseTangentAtPoint(self._b, self._c, axis3[1:]), magnitude(axis_d2))
        axis3_mdt = [-d for d in axis3_dt]
        axis2_mag = magnitude(axis2)
        axis3_mag = magnitude(axis3)

        sample_curve_on_ellipsoid = (
            lambda start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count,
                   start_weight=None, end_weight=None, overweighting=1.0, end_transition=False:
            sampleCurveOnEllipsoid(
                self._a, self._b, self._c, start_x, start_d1, start_d2, end_x, end_d1, end_d2, elements_count,
                start_weight, end_weight, overweighting, end_transition))
        move_x_to_ellipsoid_surface = lambda x: moveCoordinatesToEllipsoidSurface(self._a, self._b, self._c, x)
        move_d_to_ellipsoid_surface = lambda x, d: moveDerivativeToEllipsoidSurface(self._a, self._b, self._c, x, d)
        if self._surface_d3_mode == EllipsoidSurfaceD3Mode.SURFACE_NORMAL:
            evaluate_surface_d3_ellipsoid = lambda tx, td1, td2: set_magnitude(
                [tx[0] / (self._a * self._a), tx[1] / (self._b * self._b), tx[2] / (self._c * self._c)], dir_mag)
        else:
            evaluate_surface_d3_ellipsoid = lambda tx, td1, td2: set_magnitude(tx, dir_mag)

        octant = EllipsoidOctantMesh(self._a, self._b, self._c, half_counts, self._trans_count,
                                      nway_d_factor=self._nway_d_factor)

        # get outside curve from axis 1 to axis 2
        abx, abd1, abd2 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            axis1, axis_d2, axis_d3,
            axis2, axis_md1, axis2_dt,
            elements_count_q12)
        # get outside curve from axis 1 to axis 3
        acx, acd2, acd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            abx[0], abd2[0], abd1[0],
            axis3, axis_md1, axis3_mdt,
            elements_count_q13)
        # get outside curve from axis 2 to axis 3
        bcx, bcd2, bcd1 = sampleCurveOnEllipsoid(
            self._a, self._b, self._c,
            abx[-1], abd2[-1], abd1[-1],
            acx[-1], [-d for d in acd1[-1]], acd2[-1],
            elements_count_q23)
        # fix first/last derivatives
        abd2[0] = acd2[0]
        abd2[-1] = bcd2[0]
        acd1[-1] = [-d for d in bcd2[-1]]

        # make outer surface triangle of octant 1
        triangle_abc = QuadTriangleMesh(
            box_counts[0], box_counts[1], box_counts[2],
            sample_curve_on_ellipsoid, move_x_to_ellipsoid_surface, move_d_to_ellipsoid_surface, self._nway_d_factor)
        triangle_abc.set_edge_parameters12(abx, abd1, abd2)
        triangle_abc.set_edge_parameters13(acx, acd1, acd2)
        triangle_abc.set_edge_parameters23(bcx, bcd1, bcd2)
        triangle_abc.build()
        if not self._surface_only:
            triangle_abc.assign_d3(evaluate_surface_d3_ellipsoid)
        octant.set_triangle_abc(triangle_abc)

        if not self._surface_only:
            # extract exact derivatives
            abd2 = triangle_abc.get_edge_parameters12()[2]
            acd1 = triangle_abc.get_edge_parameters13()[1]
            bcd1 = triangle_abc.get_edge_parameters23()[1]

            # build interior lines from axis1, axis2, axis3 to origin
            aox, aod2, aod1 = sampleHermiteCurve(
                axis1, axis_md1, abd1[0], origin, axis_md1, axis_d2, elements_count=half_counts[0])
            box, bod2, bod3 = sampleHermiteCurve(
                axis2, axis_md2, bcd2[0], origin, axis_md2, axis_d3, elements_count=half_counts[1])
            bod1 = [abd1[-1]] + [axis_md1] * (len(box) - 1)
            cox, cod2, cod1 = sampleHermiteCurve(
                axis3, axis_md3, [-d for d in acd1[-1]], origin, axis_md3, axis_md2, elements_count=half_counts[2])

            # make inner surface triangle 1-2-origin
            triangle_abo = QuadTriangleMesh(
                box_counts[0], box_counts[1], self._trans_count, sampleHermiteCurve,
                nway_d_factor=self._nway_d_factor)
            abd3 = [[-d for d in evaluate_surface_d3_ellipsoid(x, None, None)] for x in abx]
            triangle_abo.set_edge_parameters12(abx, abd1, abd3, abd2)
            aod3 = [abd2[0]] + [axis_d3] * (len(aox) - 1)
            triangle_abo.set_edge_parameters13(aox, aod1, aod2, aod3)
            triangle_abo.set_edge_parameters23(box, bod1, bod2, bod3)
            triangle_abo.build()
            triangle_abo.assign_d3(lambda tx, td1, td2:
                linearlyInterpolateVectors(axis_d3, axis2_dt, magnitude(tx[1:]) / axis2_mag))
            octant.set_triangle_abo(triangle_abo)
            # extract exact derivatives
            aod1 = triangle_abo.get_edge_parameters13()[1]
            bod1 = triangle_abo.get_edge_parameters23()[1]

            # make inner surface triangle 1-3-origin
            triangle_aco = QuadTriangleMesh(
                box_counts[0], box_counts[2], self._trans_count, sampleHermiteCurve,
                nway_d_factor=self._nway_d_factor)
            acd3 = [[-d for d in evaluate_surface_d3_ellipsoid(x, None, None)] for x in acx]
            acmd1 = [[-d for d in d1] for d1 in acd1]
            triangle_aco.set_edge_parameters12(acx, acd2, acd3, acmd1)
            aomd1 = [[-d for d in d1] for d1 in aod1]
            triangle_aco.set_edge_parameters13(aox, aod3, aod2, aomd1)
            cod3 = [acd2[-1]] + [axis_md1] * (len(cox) - 1)
            triangle_aco.set_edge_parameters23(cox, cod3, cod2, cod1)
            triangle_aco.build()
            triangle_aco.assign_d3(lambda tx, td1, td2:
                linearlyInterpolateVectors(axis_md2, axis3_dt, magnitude(tx[1:]) / axis3_mag))
            octant.set_triangle_aco(triangle_aco)
            # extract exact derivatives
            cod3 = triangle_aco.get_edge_parameters23()[1]

            # make inner surface 2-3-origin
            triangle_bco = QuadTriangleMesh(
                box_counts[1], box_counts[2], self._trans_count, sampleHermiteCurve,
                nway_d_factor=self._nway_d_factor)
            bcd3 = [bod2[0]] + [[-d for d in evaluate_surface_d3_ellipsoid(x, None, None)] for x in bcx[1:-1]] \
                   + [cod2[-1]]
            bcmd1 = [[-d for d in d1] for d1 in bcd1]
            triangle_bco.set_edge_parameters12(bcx, bcd2, bcd3, bcmd1)
            bomd1 = [[-d for d in d1] for d1 in bod1]
            triangle_bco.set_edge_parameters13(box, bod3, bod2, bomd1)
            comd3 = [[-d for d in d3] for d3 in cod3]
            triangle_bco.set_edge_parameters23(cox, cod1, cod2, comd3)
            triangle_bco.build()
            triangle_bco.assign_d3(lambda tx, td1, td2: axis_d1)
            octant.set_triangle_bco(triangle_bco)

            octant.build_interior()

        return octant

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
                        # for fairness, move to surface before blending
                        new_parameter = moveCoordinatesToEllipsoidSurface(self._a, self._b, self._c, new_parameter)
                        new_parameter = [0.5 * (nx[pix][c] + new_parameter[c]) for c in range(3)]
                    else:
                        # harmonic mean to cope with significant element size differences on boundary
                        new_parameter = linearlyInterpolateVectors(
                            nx[pix], new_parameter, 0.5, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
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
                    new_derivative = linearlyInterpolateVectors(parameters[pix], new_derivative, 0.5)
                parameters[pix] = new_derivative
            # else:
            #     # not putting back values if summed parameters


    def _get_nid_to_node_layout_map_3d(self, node_layout_manager):
        """
        Get map from node identifier to special node layout.
        :param node_layout_manager: Manager of node layouts for getting standard layouts from.
        :return: map from node identifier to special node layout. No entry made if default=None layout.
        """
        node_layout_permuted = node_layout_manager.getNodeLayoutRegularPermuted(d3Defined=True)
        node_layout_3way12 = node_layout_manager.getNodeLayout3WayPoints12()
        node_layout_3way13 = node_layout_manager.getNodeLayout3WayPoints13()
        node_layout_3way23 = node_layout_manager.getNodeLayout3WayPoints23()
        node_layout_4way = node_layout_manager.getNodeLayout4WayPoints()
        nid_to_node_layout = {}
        upper_trans_counts = [self._element_counts[i] - self._trans_count for i in range(3)]
        # bottom and top transition side face nodes are fully permuted
        for nt in range(1, self._trans_count + 1):
            # bottom transition
            for n3 in (self._trans_count - nt, upper_trans_counts[2] + nt):
                for n2 in range(self._trans_count + 1, upper_trans_counts[1]):
                    for n1 in (self._trans_count - nt, upper_trans_counts[0] + nt):
                        nid = self._nids[n3][n2][n1]
                        nid_to_node_layout[nid] = node_layout_permuted
                for n1 in range(self._trans_count + 1, upper_trans_counts[0]):
                    for n2 in (self._trans_count - nt, upper_trans_counts[1] + nt):
                        nid = self._nids[n3][n2][n1]
                        nid_to_node_layout[nid] = node_layout_permuted
        for n2 in range(self._element_counts[1] + 1):
            for n1 in range(self._element_counts[0] + 1):
                # nodes on boundary between bottom transition and box are 4-way on corners, 3-way on edges,
                # fully permuted in between
                nid = self._nids[self._trans_count][n2][n1]
                if nid:
                    node_layout = node_layout_permuted
                    if n2 == self._trans_count:
                        node_layout = node_layout_3way23[0]
                        if n1 == self._trans_count:
                            node_layout = node_layout_4way[0]
                        elif n1 == upper_trans_counts[0]:
                            node_layout = node_layout_4way[1]
                    elif n2 == upper_trans_counts[1]:
                        node_layout = node_layout_3way23[1]
                        if n1 == self._trans_count:
                            node_layout = node_layout_4way[2]
                        elif n1 == upper_trans_counts[0]:
                            node_layout = node_layout_4way[3]
                    elif n1 == self._trans_count:
                        node_layout = node_layout_3way13[0]
                    elif n1 == upper_trans_counts[0]:
                        node_layout = node_layout_3way13[1]
                    nid_to_node_layout[nid] = node_layout
                # nodes on boundary between top transition and box are 4-way on corners, 3-way on edges, None in between
                nid = self._nids[upper_trans_counts[2]][n2][n1]
                if nid:
                    node_layout = None
                    if n2 == self._trans_count:
                        node_layout = node_layout_3way23[2]
                        if n1 == self._trans_count:
                            node_layout = node_layout_4way[4]
                        elif n1 == upper_trans_counts[0]:
                            node_layout = node_layout_4way[5]
                    elif n2 == upper_trans_counts[1]:
                        node_layout = node_layout_3way23[3]
                        if n1 == self._trans_count:
                            node_layout = node_layout_4way[6]
                        elif n1 == upper_trans_counts[0]:
                            node_layout = node_layout_4way[7]
                    elif n1 == self._trans_count:
                        node_layout = node_layout_3way13[2]
                    elif n1 == upper_trans_counts[0]:
                        node_layout = node_layout_3way13[3]
                    if node_layout:
                        nid_to_node_layout[nid] = node_layout
        # 3-way points on box edges up middle, fully permuted on faces
        for n3 in range(self._trans_count + 1, upper_trans_counts[2]):
            for n2 in range(self._trans_count + 1, upper_trans_counts[1]):
                for n1 in (self._trans_count, upper_trans_counts[0]):
                    nid = self._nids[n3][n2][n1]
                    nid_to_node_layout[nid] = node_layout_permuted
            for n1 in range(self._trans_count + 1, upper_trans_counts[0]):
                for n2 in (self._trans_count, upper_trans_counts[1]):
                    nid = self._nids[n3][n2][n1]
                    nid_to_node_layout[nid] = node_layout_permuted
            nid = self._nids[n3][self._trans_count][self._trans_count]
            nid_to_node_layout[nid] = node_layout_3way12[0]
            nid = self._nids[n3][self._trans_count][upper_trans_counts[0]]
            nid_to_node_layout[nid] = node_layout_3way12[1]
            nid = self._nids[n3][upper_trans_counts[1]][self._trans_count]
            nid_to_node_layout[nid] = node_layout_3way12[2]
            nid = self._nids[n3][upper_trans_counts[1]][upper_trans_counts[0]]
            nid_to_node_layout[nid] = node_layout_3way12[3]
        # 3-way points on 8 corner transitions out from 4-way points
        for nt in range(self._trans_count):
            nid = self._nids[nt][nt][nt]
            nid_to_node_layout[nid] = node_layout_3way12[1]
            nid = self._nids[nt][nt][self._element_counts[0] - nt]
            nid_to_node_layout[nid] = node_layout_3way12[0]
            nid = self._nids[nt][self._element_counts[1] - nt][nt]
            nid_to_node_layout[nid] = node_layout_3way12[3]
            nid = self._nids[nt][self._element_counts[1] - nt][self._element_counts[0] - nt]
            nid_to_node_layout[nid] = node_layout_3way12[2]
            nid = self._nids[self._element_counts[2] - nt][nt][nt]
            nid_to_node_layout[nid] = node_layout_3way12[0]
            nid = self._nids[self._element_counts[2] - nt][nt][self._element_counts[0] - nt]
            nid_to_node_layout[nid] = node_layout_3way12[1]
            nid = self._nids[self._element_counts[2] - nt][self._element_counts[1] - nt][nt]
            nid_to_node_layout[nid] = node_layout_3way12[2]
            nid = self._nids[self._element_counts[2] - nt][self._element_counts[1] - nt][self._element_counts[0] - nt]
            nid_to_node_layout[nid] = node_layout_3way12[3]
        return nid_to_node_layout

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

        # create elements

        mesh_dimension = 2 if self._surface_only else 3
        mesh = fieldmodule.findMeshByDimension(mesh_dimension)
        element_identifier = start_element_identifier
        half_counts = [count // 2 for count in self._element_counts]
        node_layout_manager = HermiteNodeLayoutManager()
        octant_mesh_group_lists = None
        if self._octant_group_lists:
            octant_mesh_group_lists = []
            for octant_group_list in self._octant_group_lists:
                octant_mesh_group_list = []
                for octant_group in octant_group_list:
                    octant_mesh_group_list.append(octant_group.getOrCreateMeshGroup(mesh))
                octant_mesh_group_lists.append(octant_mesh_group_list)
        box_mesh_group = None
        transition_mesh_group = None
        if not self._surface_only:
            if self._box_group:
                box_mesh_group = self._box_group.getOrCreateMeshGroup(mesh)
            if self._transition_group:
                transition_mesh_group = self._transition_group.getOrCreateMeshGroup(mesh)

        if self._surface_only:
            # 2-D mesh
            elementtemplate_regular = mesh.createElementtemplate()
            elementtemplate_regular.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            bicubic_hermite_serendipity_basis = (
                fieldmodule.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY))
            eft_regular = mesh.createElementfieldtemplate(bicubic_hermite_serendipity_basis)
            elementtemplate_regular.defineField(coordinates, -1, eft_regular)
            elementtemplate_special = mesh.createElementtemplate()
            elementtemplate_special.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
            # get actual indexes used on rim in 1, 2, 3 directions
            rim_indexes = [[0] + [self._trans_count + 1 + j
                                  for j in range(self._element_counts[i] - 2 * self._trans_count - 1)] +
                           [self._element_counts[i]] for i in range(3)]
            # bottom rectangle
            bottom_nids = self._nids[0]
            last_nids_row = None
            octant_n3 = 0
            for i2, n2 in enumerate(rim_indexes[1]):
                octant_n2 = 2 if (n2 > half_counts[1]) else 0
                nids_row = []
                for i1, n1 in enumerate(reversed(rim_indexes[0])):
                    octant_n1 = 1 if (n1 >= half_counts[0]) else 0
                    nids_row.append(bottom_nids[n2][n1])
                    if (i2 > 0) and (i1 > 0):
                        nids = [last_nids_row[i1 - 1], last_nids_row[i1], nids_row[i1 - 1], nids_row[i1]]
                        if None in nids:
                            continue
                        element = mesh.createElement(element_identifier, elementtemplate_regular)
                        element.setNodesByIdentifier(eft_regular, nids)
                        # print("Element", element_identifier, "nids", nids)
                        if octant_mesh_group_lists:
                            octant = octant_n3 + octant_n2 + octant_n1
                            for mesh_group in octant_mesh_group_lists[octant]:
                                mesh_group.addElement(element)
                        element_identifier += 1
                last_nids_row = nids_row
            # around sides
            node_layout_permuted = node_layout_manager.getNodeLayoutRegularPermuted(d3Defined=False)
            node_layout_triple_points = node_layout_manager.getNodeLayoutTriplePoint2D()
            index_increments = [[0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 0, 0]]
            increment_number = 0
            index_increment = index_increments[0]
            elements_count_around12 = \
                2 * (self._element_counts[0] + self._element_counts[1] - 4 * self._trans_count)
            last_nids_row = None
            last_parameters_row = None
            last_corners_row = None
            for n3 in rim_indexes[2]:
                octant_n3 = 4 if (n3 > half_counts[2]) else 0
                indexes = [self._element_counts[0], half_counts[1], n3]
                nids_row = []
                parameters_row = []
                corners_row = []
                for nc in range(elements_count_around12):
                    if nc > 0:
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
                    octant_nc = []
                    for nc in range(elements_count_around12):
                        ncp = (nc + 1) % elements_count_around12
                        nids = [last_nids_row[nc], last_nids_row[ncp],
                                nids_row[nc], nids_row[ncp]]
                        if None in nids:
                            continue
                        q = nc // quarter_elements_count_around12
                        octant_nc.append(3 if (q == 0) else (2 if (q == 1) else (0 if (q == 2) else 1)))
                        elementtemplate = elementtemplate_regular
                        eft = eft_regular
                        scalefactors = None
                        if n3 in (rim_indexes[2][1], self._element_counts[2]):
                            node_parameters = [
                                last_parameters_row[nc], last_parameters_row[ncp],
                                parameters_row[nc], parameters_row[ncp]]
                            if n3 == rim_indexes[2][1]:
                                node_layouts = [node_layout_permuted, node_layout_permuted, None, None]
                                if last_corners_row[nc]:
                                    node_layouts[0] = node_layout_triple_points[(5 - q) % 4]
                                elif last_corners_row[ncp]:
                                    node_layouts[1] = node_layout_triple_points[(5 - q) % 4]
                            else:
                                node_layouts = [None, None, node_layout_permuted, node_layout_permuted]
                                if corners_row[nc]:
                                    node_layouts[2] = node_layout_triple_points[q]
                                elif corners_row[ncp]:
                                    node_layouts[3] = node_layout_triple_points[q]
                            eft, scalefactors = \
                                determineCubicHermiteSerendipityEft(mesh, node_parameters, node_layouts)
                            elementtemplate_special.defineField(coordinates, -1, eft)
                            elementtemplate = elementtemplate_special
                        element = mesh.createElement(element_identifier, elementtemplate)
                        element.setNodesByIdentifier(eft, nids)
                        if scalefactors:
                            element.setScaleFactors(eft, scalefactors)
                        # print("Element", element_identifier, "nids", nids)
                        if octant_mesh_group_lists:
                            octant = octant_n3 + octant_nc[nc]
                            for mesh_group in octant_mesh_group_lists[octant]:
                                mesh_group.addElement(element)
                        element_identifier += 1
                last_nids_row = nids_row
                last_parameters_row = parameters_row
                last_corners_row = corners_row
            # top rectangle
            top_nids = self._nids[self._element_counts[2]]
            last_nids_row = None
            octant_n3 = 4
            for i2, n2 in enumerate(rim_indexes[1]):
                octant_n2 = 2 if (n2 > half_counts[1]) else 0
                nids_row = []
                for i1, n1 in enumerate(rim_indexes[0]):
                    octant_n1 = 1 if (n1 > half_counts[0]) else 0
                    nids_row.append(top_nids[n2][n1])
                    if (i2 > 0) and (i1 > 0):
                        nids = [last_nids_row[i1 - 1], last_nids_row[i1], nids_row[i1 - 1], nids_row[i1]]
                        if None in nids:
                            continue
                        element = mesh.createElement(element_identifier, elementtemplate_regular)
                        element.setNodesByIdentifier(eft_regular, nids)
                        # print("Element", element_identifier, "nids", nids)
                        if octant_mesh_group_lists:
                            octant = octant_n3 + octant_n2 + octant_n1
                            for mesh_group in octant_mesh_group_lists[octant]:
                                mesh_group.addElement(element)
                        element_identifier += 1
                last_nids_row = nids_row
        else:
            # 3-D mesh
            elementtemplate_regular = mesh.createElementtemplate()
            elementtemplate_regular.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            tricubic_hermite_serendipity_basis = (
                fieldmodule.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE_SERENDIPITY))
            eft_regular = mesh.createElementfieldtemplate(tricubic_hermite_serendipity_basis)
            elementtemplate_regular.defineField(coordinates, -1, eft_regular)
            elementtemplate_special = mesh.createElementtemplate()
            elementtemplate_special.setElementShapeType(Element.SHAPE_TYPE_CUBE)
            box_counts = [half_counts[i] - self._trans_count for i in range(3)]
            dbox_counts = [2 * box_counts[i] for i in range(3)]
            nid_to_node_layout = self._get_nid_to_node_layout_map_3d(node_layout_manager)
            # bottom transition
            last_nids_layer = None
            last_nx_layer = None
            for nt in range(self._trans_count + 1):
                n3 = nt
                octant_n3 = 0
                nids_layer = []
                nx_layer = []
                last_nids_row = None
                last_nx_row = None
                for i2 in range(dbox_counts[1] + 1):
                    n2 = (nt if (i2 == 0)
                          else (self._element_counts[1] - nt) if (i2 == dbox_counts[1])
                          else (self._trans_count + i2))
                    octant_n2 = 2 if (n2 > half_counts[1]) else 0
                    nids_row = []
                    nx_row = []
                    for i1 in range(dbox_counts[0] + 1):
                        n1 = (nt if (i1 == 0)
                              else (self._element_counts[0] - nt) if (i1 == dbox_counts[0])
                              else (self._trans_count + i1))
                        octant_n1 = 1 if (n1 > half_counts[0]) else 0
                        nids_row.append(self._nids[n3][n2][n1])
                        nx_row.append(self._nx[n3][n2][n1])
                        if (nt > 0) and (i2 > 0) and (i1 > 0):
                            nids = [last_nids_row[i1], last_nids_row[i1 - 1],
                                    nids_row[i1], nids_row[i1 - 1],
                                    last_nids_layer[i2 - 1][i1], last_nids_layer[i2 - 1][i1 - 1],
                                    last_nids_layer[i2][i1], last_nids_layer[i2][i1 - 1]]
                            if None in nids:
                                continue
                            elementtemplate = elementtemplate_regular
                            eft = eft_regular
                            scalefactors = None
                            node_layouts = [nid_to_node_layout.get(nid) for nid in nids]
                            if any(node_layout is not None for node_layout in node_layouts):
                                node_parameters = [last_nx_row[i1], last_nx_row[i1 - 1],
                                                   nx_row[i1], nx_row[i1 - 1],
                                                   last_nx_layer[i2 - 1][i1], last_nx_layer[i2 - 1][i1 - 1],
                                                   last_nx_layer[i2][i1], last_nx_layer[i2][i1 - 1]]
                                eft, scalefactors = \
                                    determineCubicHermiteSerendipityEft(mesh, node_parameters, node_layouts)
                                elementtemplate_special.defineField(coordinates, -1, eft)
                                elementtemplate = elementtemplate_special
                            element = mesh.createElement(element_identifier, elementtemplate)
                            element.setNodesByIdentifier(eft, nids)
                            # print("Element", element_identifier, "nids", nids)
                            if scalefactors:
                                element.setScaleFactors(eft, scalefactors)
                            if octant_mesh_group_lists:
                                octant = octant_n3 + octant_n2 + octant_n1
                                for mesh_group in octant_mesh_group_lists[octant]:
                                    mesh_group.addElement(element)
                            if transition_mesh_group:
                                transition_mesh_group.addElement(element)
                            element_identifier += 1
                    nids_layer.append(nids_row)
                    nx_layer.append(nx_row)
                    last_nids_row = nids_row
                    last_nx_row = nx_row
                last_nids_layer = nids_layer
                last_nx_layer = nx_layer
            # middle
            upper_trans_counts = [self._element_counts[i] - self._trans_count for i in range(3)]
            last_nids_layer = None
            last_nx_layer = None
            last_rim_nids_layer = None
            last_rim_nx_layer = None
            for i3 in range(dbox_counts[2] + 1):
                n3 = self._trans_count + i3
                octant_n3 = 4 if (n3 > half_counts[2]) else 0
                nids_layer = []
                nx_layer = []
                last_nids_row = None
                last_nx_row = None
                for i2 in range(dbox_counts[1] + 1):
                    n2 = self._trans_count + i2
                    octant_n2 = 2 if (n2 > half_counts[1]) else 0
                    nids_row = []
                    nx_row = []
                    for i1 in range(dbox_counts[0] + 1):
                        n1 = self._trans_count + i1
                        octant_n1 = 1 if (n1 > half_counts[0]) else 0
                        nids_row.append(self._nids[n3][n2][n1])
                        nx_row.append(self._nx[n3][n2][n1])
                        if (i3 > 0) and (i2 > 0) and (i1 > 0):
                            nids = [last_nids_layer[i2 - 1][i1 - 1], last_nids_layer[i2 - 1][i1],
                                    last_nids_layer[i2][i1 - 1], last_nids_layer[i2][i1],
                                    last_nids_row[i1 - 1], last_nids_row[i1],
                                    nids_row[i1 - 1], nids_row[i1]]
                            if None in nids:
                                continue
                            elementtemplate = elementtemplate_regular
                            eft = eft_regular
                            scalefactors = None
                            node_layouts = [nid_to_node_layout.get(nid) for nid in nids]
                            if any(node_layout is not None for node_layout in node_layouts):
                                node_parameters = [last_nx_layer[i2 - 1][i1 - 1], last_nx_layer[i2 - 1][i1],
                                                   last_nx_layer[i2][i1 - 1], last_nx_layer[i2][i1],
                                                   last_nx_row[i1 - 1], last_nx_row[i1],
                                                   nx_row[i1 - 1], nx_row[i1]]
                                eft, scalefactors = \
                                    determineCubicHermiteSerendipityEft(mesh, node_parameters, node_layouts)
                                elementtemplate_special.defineField(coordinates, -1, eft)
                                elementtemplate = elementtemplate_special
                            element = mesh.createElement(element_identifier, elementtemplate)
                            element.setNodesByIdentifier(eft, nids)
                            # print("Element", element_identifier, "nids", nids)
                            if scalefactors:
                                element.setScaleFactors(eft, scalefactors)
                            if octant_mesh_group_lists:
                                octant = octant_n3 + octant_n2 + octant_n1
                                for mesh_group in octant_mesh_group_lists[octant]:
                                    mesh_group.addElement(element)
                            if box_mesh_group:
                                box_mesh_group.addElement(element)
                            element_identifier += 1
                    nids_layer.append(nids_row)
                    nx_layer.append(nx_row)
                    last_nids_row = nids_row
                    last_nx_row = nx_row
                last_nids_layer = nids_layer
                last_nx_layer = nx_layer

                rim_nids_layer = []
                rim_nx_layer = []
                last_rim_nids_row = None
                last_rim_nx_row = None
                for nt in range(self._trans_count + 1):
                    n3 = ((self._trans_count - nt) if (i3 == 0) else (
                        (upper_trans_counts[2] + nt) if (i3 == dbox_counts[2]) else (self._trans_count + i3)))
                    rim_nids_row = []
                    rim_nx_row = []
                    octant_nc = []
                    n2 = self._trans_count - nt
                    for i1 in range(dbox_counts[0]):
                        n1 = ((self._trans_count - nt) if (i1 == 0) else (
                            (upper_trans_counts[0] + nt) if (i1 == dbox_counts[0]) else (self._trans_count + i1)))
                        rim_nids_row.append(self._nids[n3][n2][n1])
                        rim_nx_row.append(self._nx[n3][n2][n1])
                        octant_nc.append(1 if n1 >= half_counts[0] else 0)
                    n1 = upper_trans_counts[0] + nt
                    for i2 in range(dbox_counts[1]):
                        n2 = ((self._trans_count - nt) if (i2 == 0) else (
                            (upper_trans_counts[1] + nt) if (i2 == dbox_counts[1]) else (self._trans_count + i2)))
                        rim_nids_row.append(self._nids[n3][n2][n1])
                        rim_nx_row.append(self._nx[n3][n2][n1])
                        octant_nc.append(3 if n2 >= half_counts[1] else 1)
                    n2 = upper_trans_counts[1] + nt
                    for i1 in range(dbox_counts[0]):
                        n1 = ((upper_trans_counts[0] + nt) if (i1 == 0) else (
                            (self._trans_count - nt) if (i1 == dbox_counts[0]) else (upper_trans_counts[0] - i1)))
                        rim_nids_row.append(self._nids[n3][n2][n1])
                        rim_nx_row.append(self._nx[n3][n2][n1])
                        octant_nc.append(3 if n1 > half_counts[0] else 2)
                    n1 = self._trans_count - nt
                    for i2 in range(dbox_counts[1]):
                        n2 = ((upper_trans_counts[1] + nt) if (i2 == 0) else (
                             (self._trans_count - nt) if (i2 == dbox_counts[1]) else (upper_trans_counts[1] - i2)))
                        rim_nids_row.append(self._nids[n3][n2][n1])
                        rim_nx_row.append(self._nx[n3][n2][n1])
                        octant_nc.append(2 if n2 > half_counts[1] else 0)
                    if (i3 > 0) and (nt > 0):
                        rim_count = len(rim_nids_row)
                        for nc in range(rim_count):
                            ncp = (nc + 1) % rim_count
                            nids = [last_rim_nids_layer[nt - 1][nc], last_rim_nids_layer[nt - 1][ncp],
                                    last_rim_nids_row[nc], last_rim_nids_row[ncp],
                                    last_rim_nids_layer[nt][nc], last_rim_nids_layer[nt][ncp],
                                    rim_nids_row[nc], rim_nids_row[ncp]]
                            if None in nids:
                                continue
                            elementtemplate = elementtemplate_regular
                            eft = eft_regular
                            scalefactors = None
                            node_layouts = [nid_to_node_layout.get(nid) for nid in nids]
                            if any(node_layout is not None for node_layout in node_layouts):
                                node_parameters = [last_rim_nx_layer[nt - 1][nc], last_rim_nx_layer[nt - 1][ncp],
                                                   last_rim_nx_row[nc], last_rim_nx_row[ncp],
                                                   last_rim_nx_layer[nt][nc], last_rim_nx_layer[nt][ncp],
                                                   rim_nx_row[nc], rim_nx_row[ncp]]
                                eft, scalefactors = \
                                    determineCubicHermiteSerendipityEft(mesh, node_parameters, node_layouts)
                                elementtemplate_special.defineField(coordinates, -1, eft)
                                elementtemplate = elementtemplate_special
                            element = mesh.createElement(element_identifier, elementtemplate)
                            element.setNodesByIdentifier(eft, nids)
                            # print("Element", element_identifier, "nids", nids)
                            if scalefactors:
                                element.setScaleFactors(eft, scalefactors)
                            if octant_mesh_group_lists:
                                octant = octant_n3 + octant_nc[nc]
                                for mesh_group in octant_mesh_group_lists[octant]:
                                    mesh_group.addElement(element)
                            if transition_mesh_group:
                                transition_mesh_group.addElement(element)
                            element_identifier += 1
                    rim_nids_layer.append(rim_nids_row)
                    rim_nx_layer.append(rim_nx_row)
                    last_rim_nids_row = rim_nids_row
                    last_rim_nx_row = rim_nx_row
                last_rim_nids_layer = rim_nids_layer
                last_rim_nx_layer = rim_nx_layer
            # top transition
            last_nids_layer = None
            last_nx_layer = None
            for nt in range(self._trans_count, -1, -1):
                n3 = self._element_counts[2] - nt
                octant_n3 = 4
                nids_layer = []
                nx_layer = []
                last_nids_row = None
                last_nx_row = None
                for i2 in range(dbox_counts[1] + 1):
                    n2 = (nt if (i2 == 0)
                          else (self._element_counts[1] - nt) if (i2 == dbox_counts[1])
                          else (self._trans_count + i2))
                    octant_n2 = 2 if (n2 > half_counts[1]) else 0
                    nids_row = []
                    nx_row = []
                    for i1 in range(dbox_counts[0] + 1):
                        n1 = (nt if (i1 == 0)
                              else (self._element_counts[0] - nt) if (i1 == dbox_counts[0])
                              else (self._trans_count + i1))
                        octant_n1 = 1 if (n1 > half_counts[0]) else 0
                        nids_row.append(self._nids[n3][n2][n1])
                        nx_row.append(self._nx[n3][n2][n1])
                        if (nt < self._trans_count) and (i2 > 0) and (i1 > 0):
                            nids = [last_nids_layer[i2 - 1][i1 - 1], last_nids_layer[i2 - 1][i1],
                                    last_nids_layer[i2][i1 - 1], last_nids_layer[i2][i1],
                                    last_nids_row[i1 - 1], last_nids_row[i1],
                                    nids_row[i1 - 1], nids_row[i1]]
                            if None in nids:
                                continue
                            elementtemplate = elementtemplate_regular
                            eft = eft_regular
                            scalefactors = None
                            node_layouts = [nid_to_node_layout.get(nid) for nid in nids]
                            if any(node_layout is not None for node_layout in node_layouts):
                                node_parameters = [last_nx_layer[i2 - 1][i1 - 1], last_nx_layer[i2 - 1][i1],
                                                   last_nx_layer[i2][i1 - 1], last_nx_layer[i2][i1],
                                                   last_nx_row[i1 - 1], last_nx_row[i1],
                                                   nx_row[i1 - 1], nx_row[i1]]
                                eft, scalefactors = \
                                    determineCubicHermiteSerendipityEft(mesh, node_parameters, node_layouts)
                                elementtemplate_special.defineField(coordinates, -1, eft)
                                elementtemplate = elementtemplate_special
                            element = mesh.createElement(element_identifier, elementtemplate)
                            element.setNodesByIdentifier(eft, nids)
                            # print("Element", element_identifier, "nids", nids)
                            if scalefactors:
                                element.setScaleFactors(eft, scalefactors)
                            if octant_mesh_group_lists:
                                octant = octant_n3 + octant_n2 + octant_n1
                                for mesh_group in octant_mesh_group_lists[octant]:
                                    mesh_group.addElement(element)
                            if transition_mesh_group:
                                transition_mesh_group.addElement(element)
                            element_identifier += 1
                    nids_layer.append(nids_row)
                    nx_layer.append(nx_row)
                    last_nids_row = nids_row
                    last_nx_row = nx_row
                last_nids_layer = nids_layer
                last_nx_layer = nx_layer

        return node_identifier, element_identifier


class EllipsoidOctantMesh:
    """
    Generates one octant of an ellipsoid, 2-D surface or full 3-D volume.
    2-D outer surface is merely added from a QuadTriangleMesh.
    3-D volume requires 3 axis-aligned and outer surfaces to each be set from a QuadTriangleMesh,
    then the interior can be built.
    """

    def __init__(self, a, b, c, element_counts, transition_element_count, nway_d_factor=0.6):
        """
        A 2-D or 3-D Octant of an ellipsoid.
        Coordinates nx are indexed in 1, 2, 3 directions from origin at index 0, 0, 0
        with holes around the corners and 3-way points.
        :param a: Axis length (radius) in x direction.
        :param b: Axis length (radius) in y direction.
        :param c: Axis length (radius) in z direction.
        :param element_counts: Number of elements across octant only in 1, 2, 3 axes.
        :param transition_element_count: Number of transition elements around outside >= 1.
        :param nway_d_factor: Value, normally from 0.5 to 1.0 giving n-way derivative magnitude as a proportion
        of the minimum regular magnitude sampled to the n-way point. This reflects that distances from the mid-side
        of a triangle to the centre are shorter, so the derivative in the middle must be smaller.
        """
        assert all((count >= 2) for count in element_counts)
        assert 1 <= transition_element_count <= (min(element_counts) - 1)
        self._a = a
        self._b = b
        self._c = c
        self._element_counts = element_counts
        self._trans_count = transition_element_count
        self._nway_d_factor = nway_d_factor
        self._element_count12 = element_counts[0] + element_counts[1] - 2 * transition_element_count
        self._element_count13 = element_counts[0] + element_counts[2] - 2 * transition_element_count
        self._element_count23 = element_counts[1] + element_counts[2] - 2 * transition_element_count
        # counts of elements to 3-way point opposite to 3 node at axis 1, axis 2, axis 3
        self._box_counts = [self._element_counts[i] - self._trans_count for i in range(3)]
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

    def set_triangle_abc(self, trimesh: QuadTriangleMesh):
        """
        Set parameters on the outer abc surface triangle of octant.
        :param trimesh: Coordinates to set on outer surface.
        """
        assert trimesh.get_element_count12() == self._element_count12
        assert trimesh.get_element_count13() == self._element_count13
        assert trimesh.get_element_count23() == self._element_count23
        start_indexes = [self._element_counts[0], 0, 0]
        for n3 in range(self._box_counts[2]):
            px, pd1, pd2, pd3 = trimesh.get_parameters12(n3)
            self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, 2, 3]], start_indexes, [[0, 1, 0], [-1, 0, 0]])
            start_indexes[2] += 1
        start_indexes = [0, 0, self._element_counts[2]]
        for n2 in range(self._box_counts[1]):
            px, pd1, pd2, pd3 = trimesh.get_parameters31(n2, self._box_counts[0] + 1)
            self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, 2, 3]], start_indexes, [[1, 0, 0]])
            start_indexes[1] += 1
        start_indexes = [0, self._element_counts[1], self._element_counts[2]]
        px, pd1, pd2, pd3 = trimesh.get_parameters_diagonal()
        self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, 2, 3]], start_indexes, [[1, 0, 0]])

    def set_triangle_abo(self, trimesh: QuadTriangleMesh):
        """
        Set parameters on triangle 1-2-origin, an inner surface of octant.
        :param trimesh: Triangle coordinate data with x, d1, d2, optional d3.
        """
        assert trimesh.get_element_count12() == self._element_count12
        assert trimesh.get_element_count13() == self._element_counts[0]
        assert trimesh.get_element_count23() == self._element_counts[1]
        start_indexes = [self._element_counts[0], 0, 0]
        for n0 in range(self._trans_count):
            px, pd1, pd2, pd3 = trimesh.get_parameters12(n0)
            self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, -3, 2]], start_indexes, [[0, 1, 0], [-1, 0, 0]])
            start_indexes[0] -= 1
        start_indexes = [0, 0, 0]
        for n2 in range(self._box_counts[1] + 1):
            px, pd1, pd2, pd3 = (trimesh.get_parameters31(n2, self._box_counts[0] + 1) if (n2 < self._box_counts[1])
                                 else trimesh.get_parameters_diagonal())
            self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, 2, 3]], start_indexes, [[1, 0, 0]])
            start_indexes[1] += 1

    def set_triangle_aco(self, trimesh: QuadTriangleMesh):
        """
        Set parameters on triangle 1-3-origin, an inner surface of octant.
        :param trimesh: Triangle coordinate data with x, d1, d2, optional d3.
        """
        assert trimesh.get_element_count12() == self._element_count13
        assert trimesh.get_element_count13() == self._element_counts[0]
        assert trimesh.get_element_count23() == self._element_counts[2]
        start_indexes = [self._element_counts[0], 0, 0]
        for n0 in range(self._trans_count):
            px, pd1, pd2, pd3 = trimesh.get_parameters12(n0)
            self._set_coordinates_across(
                [px, pd1, pd2, pd3], [[0, 2, -3, -1], [0, -1, -3, -2]], start_indexes, [[0, 0, 1], [-1, 0, 0]])
            start_indexes[0] -= 1
        start_indexes = [0, 0, 0]
        for n3 in range(self._box_counts[2] + 1):
            px, pd1, pd2, pd3 = (trimesh.get_parameters31(n3, self._box_counts[0] + 1) if (n3 < self._box_counts[2])
                                 else trimesh.get_parameters_diagonal())
            self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 1, 3, -2]], start_indexes, [[1, 0, 0]])
            start_indexes[2] += 1

    def set_triangle_bco(self, trimesh: QuadTriangleMesh):
        """
        Set parameters on triangle 2-3-origin, an inner surface of octant.
        :param trimesh: Triangle coordinate data with x, d1, d2, optional d3.
        """
        assert trimesh.get_element_count12() == self._element_count23
        assert trimesh.get_element_count13() == self._element_counts[1]
        assert trimesh.get_element_count23() == self._element_counts[2]
        start_indexes = [0, self._element_counts[1], 0]
        for n0 in range(self._trans_count):
            px, pd1, pd2, pd3 = trimesh.get_parameters12(n0)
            self._set_coordinates_across(
                [px, pd1, pd2, pd3], [[0, 2, -3, -1], [0, -2, -3, 1]], start_indexes, [[0, 0, 1], [0, -1, 0]])
            start_indexes[1] -= 1
        start_indexes = [0, 0, 0]
        for n3 in range(self._box_counts[2] + 1):
            px, pd1, pd2, pd3 = (trimesh.get_parameters31(n3, self._box_counts[1] + 1) if (n3 < self._box_counts[2])
                                 else trimesh.get_parameters_diagonal())
            self._set_coordinates_across([px, pd1, pd2, pd3], [[0, 2, 3, 1]], start_indexes, [[0, 1, 0]])
            start_indexes[2] += 1

    def _get_transitions(self, indexes):
        """
        For each index direction, get False if in core or True if in transition zone.
        :param indexes: Location indexes in 1, 2, 3 directions.
        :return: Transition 1, 2, 3 directions.
        """
        return [(indexes[i] - self._box_counts[i]) > 0 for i in range(3)]

    def _set_coordinates_across(self, parameters, parameter_indexes, start_indexes, index_increments,
                                skip_start=False, skip_end=False, blend=False):
        """
        Insert parameters across the coordinates array.
        :param parameters: List of lists of N node parameters e.g. [px, pd1, pd2, pd3]
        :param parameter_indexes: Lists of parameter indexes where x=0, d1=1, d2=2, d3=3. Starts with first and
        advances at transitions change, then stays on the last. Can be negative to invert vector.
        e.g. [[0, 1, 2], [0, -2, 1]] for [x, d1, d2] then [x, -d2, d1] from first corner.
        :param start_indexes: Indexes into nx array for start point.
        :param index_increments: List of increments in indexes. Starts with first and uses next at each transition
        change, then stays on the last.
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
        last_trans = self._get_transitions(indexes)
        for n in range(start_n, limit_n):
            if n > 0:
                while True:
                    indexes = [indexes[c] + index_increment[c] for c in range(3)]
                    # skip over blank transition coordinates
                    if self._nx[indexes[2]][indexes[1]][indexes[0]]:
                        break
            trans = self._get_transitions(indexes)
            if last_trans and (trans != last_trans):
                if parameter_number < (len(parameter_indexes) - 1):
                    parameter_number += 1
                    parameter_index = parameter_indexes[parameter_number]
                if increment_number < (len(index_increments) - 1):
                    increment_number += 1
                    index_increment = index_increments[increment_number]
            nx = self._nx[indexes[2]][indexes[1]][indexes[0]]
            for parameter, spix in zip(parameters, parameter_index):
                if not parameter[n]:
                    continue
                new_parameter = [-d for d in parameter[n]] if (spix < 0) else copy.copy(parameter[n])
                pix = abs(spix)
                if blend and nx[pix]:
                    if pix == 0:
                        # for fairness, move to surface before blending
                        if any(indexes[i] == self._element_counts[i] for i in range(3)):
                            new_parameter = moveCoordinatesToEllipsoidSurface(self._a, self._b, self._c, new_parameter)
                        new_parameter = [0.5 * (nx[pix][c] + new_parameter[c]) for c in range(3)]
                    else:
                        # harmonic mean to cope with significant element size differences on boundary
                        new_parameter = linearlyInterpolateVectors(
                            nx[pix], new_parameter, 0.5, magnitudeScalingMode=DerivativeScalingMode.HARMONIC_MEAN)
                nx[pix] = new_parameter
            last_trans = trans

    def _smooth_derivative_across(self, start_indexes, end_indexes, index_increments, derivative_indexes,
                                  fix_start_direction=True, fix_end_direction=True, move_d_to_surface=False):
        """
        Smooth derivatives across octant.
        :param start_indexes: Indexes of first point.
        :param end_indexes: Indexes of last point.
        :param index_increments: List of increments in indexes. Starts with first and advances at transitions change,
        then stays on the last.
        :param derivative_indexes: List of signed derivative parameter index to set along where 1=d1, 2=d2, 3=d3.
        Starts with first and advances at transitions change, then stays on the last. Can be negative to invert vector.
        e.g. [1, -2] for d1 then -d2 from first transition change.
        :param fix_start_direction: Set to True to keep the start direction but scale its magnitude.
        :param fix_end_direction: Set to True to keep the end direction but scale its magnitude.
        :param move_d_to_surface: Set to True to force derivatives to surface tangents at their current location.
        Only use for smoothing over the outer surface of ellipsoid.
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
        last_trans = self._get_transitions(indexes)
        while True:
            if n > 0:
                if indexes == end_indexes:
                    break
                while True:
                    indexes = [indexes[i] + index_increment[i] for i in range(3)]
                    # skip over blank coordinates in transition zone
                    if self._nx[indexes[2]][indexes[1]][indexes[0]]:
                        break
            trans = self._get_transitions(indexes)
            if last_trans and (trans != last_trans):
                if derivative_number < (len(derivative_indexes) - 1):
                    derivative_number += 1
                    derivative_index = derivative_indexes[derivative_number]
                if increment_number < (len(index_increments) - 1):
                    increment_number += 1
                    index_increment = index_increments[increment_number]
            parameters = self._nx[indexes[2]][indexes[1]][indexes[0]]
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
            last_trans = trans
        sd = smoothCubicHermiteDerivativesLine(
            px, pd, fixStartDirection=fix_start_direction, fixEndDirection=fix_end_direction)
        if move_d_to_surface:
            for n in range(1, len(sd) - 1):
                sd[n] = self._move_d_to_surface(px[n], sd[n])
            sd = smoothCubicHermiteDerivativesLine(px, sd, fixAllDirections=True)
        for n in range(len(sd)):
            indexes = indexes_list[n]
            spix = derivative_index_list[n]
            new_derivative = [-d for d in sd[n]] if (spix < 0) else sd[n]
            pix = abs(spix)
            self._nx[indexes[2]][indexes[1]][indexes[0]][pix] = new_derivative

    def build_interior(self):
        """
        Determine interior coordinates from surface coordinates.
        """
        # determine 4-way point location from mean curves between side points linking to it
        point12 = self._nx[0][self._box_counts[1]][self._box_counts[0]]
        point13 = self._nx[self._box_counts[2]][0][self._box_counts[0]]
        point23 = self._nx[self._box_counts[2]][self._box_counts[1]][0]
        point123 = self._nx[self._element_counts[2]][self._element_counts[1]][self._element_counts[0]]

        x_4way, d_4way = get_nway_point(
            [point23[0], point13[0], point12[0], point123[0]],
            [point23[1], point13[2], point12[3], [-d for d in point123[3]]],
            [self._box_counts[0], self._box_counts[1], self._box_counts[2], self._trans_count],
            sampleHermiteCurve, nway_d_factor=self._nway_d_factor)

        # smooth sample from sides to 3-way points using end derivatives
        min_weight = 1  # GRC revisit, remove?
        ax, ad1 = sampleHermiteCurve(
            point23[0], point23[1], None, x_4way, d_4way[0], None, self._box_counts[0],
            start_weight=self._box_counts[0] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
        bx, bd2 = sampleHermiteCurve(
            point13[0], point13[2], None, x_4way, d_4way[1], None, self._box_counts[1],
            start_weight=self._box_counts[1] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
        cx, cd3 = sampleHermiteCurve(
            point12[0], point12[3], None, x_4way, d_4way[2], None, self._box_counts[2],
            start_weight=self._box_counts[2] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
        tx, td3 = sampleHermiteCurve(
            point123[0], [-d for d in point123[3]], None, x_4way, d_4way[3], None, self._trans_count,
            start_weight=self._trans_count + min_weight, end_weight=1.0 + min_weight, end_transition=True)

        self._set_coordinates_across([ax, ad1], [[0, 1]], [0, self._box_counts[1], self._box_counts[2]], [[1, 0, 0]])
        self._set_coordinates_across([bx, bd2], [[0, 2]], [self._box_counts[0], 0, self._box_counts[2]], [[0, 1, 0]])
        self._set_coordinates_across([cx, cd3], [[0, 3]], [self._box_counts[0], self._box_counts[1], 0], [[0, 0, 1]])
        self._set_coordinates_across([tx, td3], [[0, -3]],
                                     [self._element_counts[0], self._element_counts[1], self._element_counts[2]],
                                     [[-1, -1, -1]], skip_end=True)

        # sample up to 3-way lines connecting to 4-way point
        for n3 in range(1, self._box_counts[2]):
            point13 = self._nx[n3][0][self._box_counts[0]]
            point23 = self._nx[n3][self._box_counts[1]][0]
            point123 = self._nx[n3][self._element_counts[1]][self._element_counts[0]]
            point_3way = self._nx[n3][self._box_counts[1]][self._box_counts[0]]

            x_3way, d_3way = get_nway_point(
                [point23[0], point13[0], point123[0]],
                [point23[1], point13[2], [-d for d in point123[3]]],
                [self._box_counts[0], self._box_counts[1], self._trans_count],
                sampleHermiteCurve, prescribed_x_nway=point_3way[0], nway_d_factor=self._nway_d_factor)

            ax, ad1, ad2 = sampleHermiteCurve(
                point23[0], point23[1], point23[2], x_3way, d_3way[0], None, self._box_counts[0],
                start_weight=self._box_counts[0] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
            ad1[-1] = d_3way[0]
            ad2[-1] = d_3way[1]
            bx, bd2, bd1 = sampleHermiteCurve(
                point13[0], point13[2], point13[1], x_3way, d_3way[1], None, self._box_counts[1],
                start_weight=self._box_counts[1] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
            bd1[-1] = d_3way[0]
            bd2[-1] = d_3way[1]
            tx, td3, td1 = sampleHermiteCurve(
                point123[0], [-d for d in point123[3]], point123[1], x_3way, d_3way[2], None,
                self._trans_count, start_weight=self._trans_count + min_weight,
                end_weight=1.0 + min_weight, end_transition=True)

            self._set_coordinates_across([ax, ad1, ad2], [[0, 1, 2]], [0, self._box_counts[1], n3], [[1, 0, 0]])
            self._set_coordinates_across([bx, bd1, bd2], [[0, 1, 2]], [self._box_counts[0], 0, n3], [[0, 1, 0]])
            self._set_coordinates_across(
                [tx, td1, td3], [[0, 1, -3]], [self._element_counts[0], self._element_counts[1], n3], [[-1, -1, 0]],
                skip_end=True)

        for n2 in range(1, self._box_counts[1]):
            point12 = self._nx[0][n2][self._box_counts[0]]
            point23 = self._nx[self._box_counts[2]][n2][0]
            point123 = self._nx[self._element_counts[2]][n2][self._element_counts[0]]
            point_3way = self._nx[self._box_counts[2]][n2][self._box_counts[0]]

            x_3way, d_3way = get_nway_point(
                [point23[0], point12[0], point123[0]],
                [point23[1], point12[3], [-d for d in point123[3]]],
                [self._box_counts[0], self._box_counts[2], self._trans_count],
                sampleHermiteCurve, prescribed_x_nway=point_3way[0], nway_d_factor=self._nway_d_factor)

            ax, ad1, ad3 = sampleHermiteCurve(
                point23[0], point23[1], point23[3], x_3way, d_3way[0], None, self._box_counts[0],
                start_weight=self._box_counts[0] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
            ad1[-1] = d_3way[0]
            ad3[-1] = d_3way[1]
            bx, bd3, bd1 = sampleHermiteCurve(
                point12[0], point12[3], point12[1], x_3way, d_3way[1], None, self._box_counts[2],
                start_weight=self._box_counts[2] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
            bd1[-1] = d_3way[0]
            bd3[-1] = d_3way[1]
            tx, td3, td1 = sampleHermiteCurve(
                point123[0], [-d for d in point123[3]], point123[1], x_3way, d_3way[2], None,
                self._trans_count, start_weight=self._trans_count + min_weight,
                end_weight=1.0 + min_weight, end_transition=True)

            self._set_coordinates_across(
                [ax, ad1, ad3], [[0, 1, 3]], [0, n2, self._box_counts[2]], [[1, 0, 0]], blend=True)
            self._set_coordinates_across(
                [bx, bd1, bd3], [[0, 1, 3]], [self._box_counts[0], n2, 0], [[0, 0, 1]], blend=True)
            self._set_coordinates_across(
                [tx, td1, td3], [[0, 1, -3]], [self._element_counts[0], n2, self._element_counts[2]], [[-1, 0, -1]],
                skip_end=True, blend=True)

        for n1 in range(1, self._box_counts[0]):
            point12 = self._nx[0][self._box_counts[1]][n1]
            point13 = self._nx[self._box_counts[2]][0][n1]
            point123 = self._nx[self._element_counts[2]][self._element_counts[1]][n1]
            point_3way = self._nx[self._box_counts[2]][self._box_counts[1]][n1]

            x_3way, d_3way = get_nway_point(
                [point13[0], point12[0], point123[0]],
                [point13[2], point12[3], [-d for d in point123[3]]],
                [self._box_counts[0], self._box_counts[2], self._trans_count],
                sampleHermiteCurve, prescribed_x_nway=point_3way[0], nway_d_factor=self._nway_d_factor)

            ax, ad2, ad3 = sampleHermiteCurve(
                point13[0], point13[2], point13[3], x_3way, d_3way[0], None, self._box_counts[1],
                start_weight=self._box_counts[1] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
            ad2[-1] = d_3way[0]
            ad3[-1] = d_3way[1]
            bx, bd3, bd2 = sampleHermiteCurve(
                point12[0], point12[3], point12[2], x_3way, d_3way[1], None, self._box_counts[2],
                start_weight=self._box_counts[2] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
            bd2[-1] = d_3way[0]
            bd3[-1] = d_3way[1]
            tx, td3, td2 = sampleHermiteCurve(
                point123[0], [-d for d in point123[3]], point123[2], x_3way, d_3way[2], None,
                self._trans_count, start_weight=self._trans_count + min_weight,
                end_weight=1.0 + min_weight, end_transition=True)

            self._set_coordinates_across(
                [ax, ad2, ad3], [[0, 2, 3]], [n1, 0, self._box_counts[2]], [[0, 1, 0]], blend=True)
            self._set_coordinates_across(
                [bx, bd2, bd3], [[0, 2, 3]], [n1, self._box_counts[1], 0], [[0, 0, 1]], blend=True)
            self._set_coordinates_across(
                [tx, td2, td3], [[0, 2, -3]], [n1, self._element_counts[1], self._element_counts[2]],
                [[0, -1, -1]], skip_end=True, blend=True)

        for nt in range(1, self._trans_count):
            point12 = self._nx[0][self._element_counts[1] - nt][self._element_counts[0] - nt]
            point13 = self._nx[self._element_counts[2] - nt][0][self._element_counts[0] - nt]
            point23 = self._nx[self._element_counts[2] - nt][self._element_counts[1] - nt][0]
            point_3way = \
                self._nx[self._element_counts[2] - nt][self._element_counts[1] - nt][self._element_counts[0] - nt]

            x_3way, d_3way = get_nway_point(
                [point23[0], point13[0], point12[0]],
                [point23[1], point13[2], point12[2]],
                [self._box_counts[0], self._box_counts[1], self._box_counts[2]],
                sampleHermiteCurve, prescribed_x_nway=point_3way[0], nway_d_factor=self._nway_d_factor)

            ax, ad1, ad2 = sampleHermiteCurve(
                point23[0], point23[1], point23[2], x_3way, d_3way[0], None, self._box_counts[0],
                start_weight=self._box_counts[0] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
            ad1[-1] = d_3way[0]
            ad2[-1] = d_3way[1]
            bx, bd2, bd1 = sampleHermiteCurve(
                point13[0], point13[2], point13[1], x_3way, d_3way[1], None, self._box_counts[1],
                start_weight=self._box_counts[1] + min_weight, end_weight=1.0 + min_weight, end_transition=True)
            bd1[-1] = d_3way[0]
            bd2[-1] = d_3way[1]
            cx, cd2, cd1 = sampleHermiteCurve(
                point12[0], point12[2], point12[1], x_3way, d_3way[2], None, self._box_counts[2],
                start_weight=self._box_counts[2] + min_weight, end_weight=1.0 + min_weight,
                end_transition=True)

            self._set_coordinates_across(
                [ax, ad1, ad2], [[0, 1, 2]], [0, self._element_counts[1] - nt, self._element_counts[2] - nt],
                [[1, 0, 0]], blend=True)
            self._set_coordinates_across(
                [bx, bd1, bd2], [[0, 1, 2]], [self._element_counts[0] - nt, 0, self._element_counts[2] - nt],
                [[0, 1, 0]], blend=True)
            self._set_coordinates_across(
                [cx, cd1, cd2], [[0, 1, 2]], [self._element_counts[0] - nt, self._element_counts[1] - nt, 0],
                [[0, 0, 1]], skip_end=True, blend=True)

        # average point coordinates across 3 directions between side faces and surfaces to 4 3-way lines.
        min_weight = 1  # GRC revisit, remove?
        # 1-direction
        for n2 in range(1, self._box_counts[1]):
            for n3 in range(1, self._box_counts[2]):
                start_indexes = [0, n2, n3]
                corner_indexes = [self._box_counts[0], n2, n3]
                end_indexes = [self._element_counts[0], n2, n3]
                start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
                corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
                end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
                px, _ = sampleHermiteCurve(
                    start[0], start[1], None, corner[0], corner[1], None, self._box_counts[0],
                    start_weight=self._box_counts[0] + min_weight, end_weight=1.0 + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], start_indexes, [[1, 0, 0]], skip_start=True, skip_end=True, blend=True)
                px, _ = sampleHermiteCurve(
                    corner[0], corner[1], None, end[0], end[3], None, self._trans_count,
                    start_weight=1.0 + min_weight, end_weight=self._trans_count + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], corner_indexes, [[1, 0, 0]], skip_start=True, skip_end=True, blend=True)
            for nt in range(1, self._trans_count):
                start_indexes = [0, n2, self._box_counts[2] + nt]
                corner_indexes = [self._box_counts[0] + nt, n2, self._box_counts[2] + nt]
                end_indexes = [self._box_counts[0] + nt, n2, 0]
                start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
                corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
                end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
                px, _ = sampleHermiteCurve(
                    start[0], start[1], None, corner[0], corner[1], None, self._box_counts[0],
                    start_weight=self._box_counts[0] + min_weight, end_weight=1.0 + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], start_indexes, [[1, 0, 0]], skip_start=True, skip_end=True, blend=True)
                px, _ = sampleHermiteCurve(
                    corner[0], corner[1], None, end[0], [-d for d in end[2]], None, self._box_counts[2],
                    start_weight=1.0 + min_weight, end_weight=self._box_counts[2] + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], corner_indexes, [[0, 0, -1]], skip_start=True, skip_end=True, blend=True)
        # 2-direction
        for n3 in range(1, self._box_counts[2]):
            for n1 in range(1, self._box_counts[0]):
                start_indexes = [n1, 0, n3]
                corner_indexes = [n1, self._box_counts[1], n3]
                end_indexes = [n1, self._element_counts[1], n3]
                start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
                corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
                end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
                px, _ = sampleHermiteCurve(
                    start[0], start[2], None, corner[0], corner[2], None, self._box_counts[1],
                    start_weight=self._box_counts[1] + min_weight, end_weight=1.0 + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], start_indexes, [[0, 1, 0]], skip_start=True, skip_end=True, blend=True)
                px, _ = sampleHermiteCurve(
                    corner[0], corner[2], None, end[0], end[3], None, self._trans_count,
                    start_weight=1.0 + min_weight, end_weight=self._trans_count + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], corner_indexes, [[0, 1, 0]], skip_start=True, skip_end=True, blend=True)
            for nt in range(1, self._trans_count):
                start_indexes = [self._box_counts[0] + nt, 0, n3]
                corner_indexes = [self._box_counts[0] + nt, self._box_counts[1] + nt, n3]
                end_indexes = [0, self._box_counts[1] + nt, n3]
                start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
                corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
                end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
                px, _ = sampleHermiteCurve(
                    start[0], start[1], None, corner[0], corner[1], None, self._box_counts[1],
                    start_weight=self._box_counts[1] + min_weight, end_weight=1.0 + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], start_indexes, [[0, 1, 0]], skip_start=True, skip_end=True, blend=True)
                px, _ = sampleHermiteCurve(
                    corner[0], corner[1], None, end[0], end[1], None, self._box_counts[0],
                    start_weight=1.0 + min_weight, end_weight=self._box_counts[0] + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], corner_indexes, [[-1, 0, 0]], skip_start=True, skip_end=True, blend=True)
        # 3-direction
        for n1 in range(1, self._box_counts[0]):
            for n2 in range(1, self._box_counts[1]):
                start_indexes = [n1, n2, 0]
                corner_indexes = [n1, n2, self._box_counts[2]]
                end_indexes = [n1, n2, self._element_counts[2]]
                start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
                corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
                end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
                px, _ = sampleHermiteCurve(
                    start[0], start[3], None, corner[0], corner[3], None, self._box_counts[2],
                    start_weight=self._box_counts[2] + min_weight, end_weight=1.0 + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], start_indexes, [[0, 0, 1]], skip_start=True, skip_end=True, blend=True)
                px, _ = sampleHermiteCurve(
                    corner[0], corner[3], None, end[0], end[3], None, self._trans_count,
                    start_weight=1.0 + min_weight, end_weight=self._trans_count + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], corner_indexes, [[0, 0, 1]], skip_start=True, skip_end=True, blend=True)
            for nt in range(1, self._trans_count):
                start_indexes = [n1, self._box_counts[1] + nt, 0]
                corner_indexes = [n1, self._box_counts[1] + nt, self._box_counts[2] + nt]
                end_indexes = [n1, 0, self._box_counts[2] + nt]
                start = self._nx[start_indexes[2]][start_indexes[1]][start_indexes[0]]
                corner = self._nx[corner_indexes[2]][corner_indexes[1]][corner_indexes[0]]
                end = self._nx[end_indexes[2]][end_indexes[1]][end_indexes[0]]
                px, _ = sampleHermiteCurve(
                    start[0], start[2], None, corner[0], [-d for d in corner[2]], None, self._box_counts[2],
                    start_weight=self._box_counts[2] + min_weight, end_weight=1.0 + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], start_indexes, [[0, 0, 1]], skip_start=True, skip_end=True, blend=True)
                px, _ = sampleHermiteCurve(
                    corner[0], [-d for d in corner[2]], None, end[0], [-d for d in end[2]], None, self._box_counts[1],
                    start_weight=1.0 + min_weight, end_weight=self._box_counts[1] + min_weight)
                self._set_coordinates_across(
                    [px], [[0]], corner_indexes, [[0, -1, 0]], skip_start=True, skip_end=True, blend=True)

        # smooth 1-direction
        for n2 in range(1, self._box_counts[1]):
            for n3 in range(1, self._box_counts[2]):
                self._smooth_derivative_across(
                    [0, n2, n3], [self._element_counts[0], n2, n3],
                    [[1, 0, 0]], [1, 3], fix_start_direction=True, fix_end_direction=True)
            for nt in range(1, self._trans_count):
                self._smooth_derivative_across(
                    [0, n2, self._box_counts[2] + nt], [self._box_counts[0] + nt, n2, 0],
                    [[1, 0, 0], [0, 0, -1]], [1, 1, -2], fix_start_direction=True, fix_end_direction=True)
        # smooth 2-direction
        for n3 in range(1, self._box_counts[2]):
            for n1 in range(1, self._box_counts[0]):
                self._smooth_derivative_across(
                    [n1, 0, n3], [n1, self._element_counts[1], n3],
                    [[0, 1, 0]], [2, 3], fix_start_direction=True, fix_end_direction=True)
            for nt in range(1, self._trans_count):
                self._smooth_derivative_across(
                    [self._box_counts[0] + nt, 0, n3], [0, self._box_counts[1] + nt, n3],
                    [[0, 1, 0], [-1, 0, 0]], [1], fix_start_direction=True, fix_end_direction=True)
        # smooth 3-direction
        for n1 in range(1, self._box_counts[0]):
            for n2 in range(1, self._box_counts[1]):
                self._smooth_derivative_across(
                    [n1, n2, 0], [n1, n2, self._element_counts[2]],
                    [[0, 0, 1]], [3], fix_start_direction=True, fix_end_direction=True)
            for nt in range(1, self._trans_count):
                self._smooth_derivative_across(
                    [n1, self._box_counts[1] + nt, 0], [n1, 0, self._box_counts[2] + nt],
                    [[0, 0, 1], [0, -1, 0]], [2, -2], fix_start_direction=True, fix_end_direction=True)
