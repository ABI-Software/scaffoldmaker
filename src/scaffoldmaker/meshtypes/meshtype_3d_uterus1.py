"""
Generates a 3-D uterus mesh from a 1-D network layout, with variable
numbers of elements around, along and through wall.
"""

import copy
import math

from cmlibs.maths.vectorops import add, cross, mult, normalize, sub, magnitude, set_magnitude
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm, \
    findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.uterus_terms import get_uterus_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.interpolation import (
    computeCubicHermiteDerivativeScaling, interpolateCubicHermite, interpolateLagrangeHermiteDerivative,
    smoothCubicHermiteDerivativesLine)
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters


def get_bifurcation_triple_point(p1x, p1d, p2x, p2d, p3x, p3d):
    """
    Get coordinates and derivatives of triple point between p1, p2 and p3 with derivatives.
    :param p1x..p3d: Point coordinates and derivatives, numbered anticlockwise around triple point.
    All derivatives point away from triple point.
    Returned d1 points from triple point to p2, d2 points from triple point to p3.
    :return: x, d1, d2
    """
    scaling12 = computeCubicHermiteDerivativeScaling(p1x, [-d for d in p1d], p2x, p2d)
    scaling23 = computeCubicHermiteDerivativeScaling(p2x, [-d for d in p2d], p3x, p3d)
    scaling31 = computeCubicHermiteDerivativeScaling(p3x, [-d for d in p3d], p1x, p1d)
    trx1 = interpolateCubicHermite(p1x, mult(p1d, -scaling12), p2x, mult(p2d, scaling12), 0.5)
    trx2 = interpolateCubicHermite(p2x, mult(p2d, -scaling23), p3x, mult(p3d, scaling23), 0.5)
    trx3 = interpolateCubicHermite(p3x, mult(p3d, -scaling31), p1x, mult(p1d, scaling31), 0.5)
    trx = [(trx1[c] + trx2[c] + trx3[c]) / 3.0 for c in range(3)]
    td1 = interpolateLagrangeHermiteDerivative(trx, p1x, p1d, 0.0)
    td2 = interpolateLagrangeHermiteDerivative(trx, p2x, p2d, 0.0)
    td3 = interpolateLagrangeHermiteDerivative(trx, p3x, p3d, 0.0)
    n12 = cross(td1, td2)
    n23 = cross(td2, td3)
    n31 = cross(td3, td1)
    norm = normalize([(n12[c] + n23[c] + n31[c]) for c in range(3)])
    sd1 = smoothCubicHermiteDerivativesLine([trx, p1x], [normalize(cross(norm, cross(td1, norm))), p1d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    sd2 = smoothCubicHermiteDerivativesLine([trx, p2x], [normalize(cross(norm, cross(td2, norm))), p2d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    sd3 = smoothCubicHermiteDerivativesLine([trx, p3x], [normalize(cross(norm, cross(td3, norm))), p3d],
                                            fixStartDirection=True, fixEndDerivative=True)[0]
    trd1 = mult(sub(sd2, add(sd3, sd1)), 0.5)
    trd2 = mult(sub(sd3, add(sd1, sd2)), 0.5)
    return trx, trd1, trd2


def get_tube_bifurcation_connection_elements_counts(tCounts):
    """
    Get number of elements directly connecting tubes 1, 2 and 3 from the supplied number around.
    :param tCounts: Number of elements around tubes in order.
    :return: List of elements connect tube with its next neighbour, looping back to first.
    """
    assert len(tCounts) == 3
    return [(tCounts[i] + tCounts[i - 2] - tCounts[i - 1]) // 2 for i in range(3)]


class MeshType_3d_uterus1(Scaffold_base):
    """
    Generates a 3-D uterus mesh from a 1-D network layout with variable numbers of elements around, along and through
    wall.
    Magnitude of D2 and D3 are the radii of the uterus in the respective directions.
    """
    parameterSetStructureStrings = {
        'Mouse 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5, 6-7-8-9-5.2, 5.3-10-11, 11-12-13, 13-14"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[24.20, 0.00, 36.56], [-2.53, 0.00, -14.75], [2.42, 0.00, -0.41], [0.13, 0.00, -0.34], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (2, [[20.88, 0.00, 23.34], [-4.10, 0.00, -11.62], [2.47, 0.00, -0.87], [-0.03, 0.00, -0.57], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (3, [[16.28, 0.00, 13.37], [-5.76, 0.00, -8.94], [2.38, 0.00, -1.54], [-0.41, 0.00, -0.56], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (4, [[9.59, 0.00, 5.63], [-8.24, 0.00, -6.86], [1.67, 0.00, -2.00], [-0.75, 0.00, -0.33], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (5, [[0.00, 0.00, 0.00], [[-10.77, 0.00, -4.33], [10.77, 0.00, -4.33], [0.00, 0.00, -4.41]],
                         [[0.88, 0.00, -2.18], [0.83, 0.11, 2.07], [2.67, 0.00, 0.00]], [-0.32, 0.00, 0.38],
                         [[0.00, -2.50, 0.00], [0.05, -2.50, 0.11], [0.00, -2.50, 0.00]], [0.00, 0.00, 0.00]]),
                    (6, [[-24.20, 0.00, 36.56], [2.53, 0.00, -14.75], [2.42, 0.00, 0.41], [0.40, 0.00, 0.96], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (7, [[-20.88, 0.00, 23.34], [4.10, 0.00, -11.62], [2.47, 0.00, 0.87], [-0.03, 0.00, 0.57], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (8, [[-16.28, 0.00, 13.37], [5.76, 0.00, -8.94], [2.38, 0.00, 1.54], [-0.41, 0.00, 0.56], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (9, [[-9.59, 0.00, 5.63], [8.24, 0.00, -6.86], [1.67, 0.00, 2.00], [-0.05, 0.00, -0.59], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (10, [[0.00, 0.00, -4.05], [0.00, 0.00, -3.69], [2.50, 0.00, 0.00], [0.16, 0.00, -0.39], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (11, [[0.00, 0.00, -7.37], [0.00, 0.00, -3.46], [2.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (12, [[0.00, 0.00, -10.98], [0.00, 0.00, -3.31], [2.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (13, [[0.00, 0.00, -14.00], [0.00, 0.00, -3.51], [2.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]]),
                    (14, [[0.00, 0.00, -18.00], [0.00, 0.00, -4.49], [2.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.50, 0.00], [0.00, 0.00, 0.00]])]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_uterus_term('right uterine horn')[0],
                    'ontId': get_uterus_term('right uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_uterus_term('left uterine horn')[0],
                    'ontId': get_uterus_term('left uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9-10',
                    'name': get_uterus_term('body of uterus')[0],
                    'ontId': get_uterus_term('body of uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-12',
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '13',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                }]
        }),
        'Rat 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5, 6-7-8-9-5.2, 5.3-10, 10-11, 11-12-13-14-15-16"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[25.00, 0.00, 43.53], [0.06, 0.00, -11.85], [2.27, 0.09, 0.01], [0.15, -0.06, -0.57], [0.09, -2.33, 0.00], [-0.08, 0.74, -0.01]]),
                    (2, [[23.55, 0.00, 31.29], [-2.96, 0.00, -12.51], [2.27, 0.04, -0.54], [-0.15, -0.04, -0.53], [0.03, -1.91, -0.01], [-0.04, 0.10, 0.00]]),
                    (3, [[19.03, 0.00, 18.76], [-6.24, 0.00, -11.75], [1.96, 0.01, -1.04], [-0.29, 0.01, -0.58], [0.01, -2.16, 0.00], [0.00, -0.20, -0.01]]),
                    (4, [[11.24, 0.00, 8.09], [-9.64, 0.00, -9.57], [1.69, 0.05, -1.70], [-0.44, -0.01, -0.52], [0.03, -2.32, -0.03], [-0.01, -0.14, 0.01]]),
                    (5, [[0.00, 0.00, 0.00], [[-12.65, 0.00, -6.51], [12.65, 0.00, -6.51], [0.00, 0.00, -4.00]],
                         [[1.07, -0.02, -2.07], [1.07, 0.02, 2.07], [4.65, 0.00, 0.00]], [-0.23, -0.03, 0.15],
                         [[-0.01, -2.44, 0.02], [0.01, -2.44, 0.02], [0.00, -2.53, 0.00]], [-0.01, -0.07, 0.04]]),
                    (6, [[-25.00, 0.00, 43.53], [-0.06, 0.00, -11.85], [2.27, 0.09, -0.01], [0.24, -0.02, 0.84], [0.09, -2.33, 0.00], [-0.03, 0.36, 0.00]]),
                    (7, [[-23.55, 0.00, 31.29], [2.96, 0.00, -12.51], [2.27, 0.04, 0.54], [-0.15, -0.04, 0.53], [0.03, -1.91, 0.01], [-0.04, 0.10, 0.00]]),
                    (8, [[-19.03, 0.00, 18.76], [6.24, 0.00, -11.75], [1.96, 0.01, 1.04], [-0.29, 0.01, 0.58], [0.01, -2.16, 0.00], [0.00, -0.20, 0.01]]),
                    (9, [[-11.24, 0.00, 8.09], [9.64, 0.00, -9.57], [1.69, 0.05, 1.70], [1.17, 0.00, -0.39], [0.03, -2.32, 0.03], [0.00, -0.21, 0.00]]),
                    (10, [[0.00, 0.00, -4.00], [0.00, 0.00, -4.00], [4.65, 0.00, 0.00], [0.58, -0.01, -0.33], [0.00, -2.60, 0.00], [-0.01, -0.05, -0.01]]),
                    (11, [[0.00, 0.00, -8.00], [0.00, 0.00, -4.00], [4.65, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.60, 0.00], [0.00, 0.00, 0.00]]),
                    (12, [[0.00, 0.00, -12.00], [0.00, 0.00, -4.00], [4.65, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.60, 0.00], [0.00, 0.00, 0.00]]),
                    (13, [[0.00, 0.00, -16.00], [0.00, 0.00, -4.00], [4.65, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.60, 0.00], [0.00, 0.00, 0.00]]),
                    (14, [[0.00, 0.00, -20.00], [0.00, 0.00, -4.00], [4.65, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.60, 0.00], [0.00, 0.00, 0.00]]),
                    (15, [[0.00, 0.00, -24.00], [0.00, 0.00, -4.00], [4.65, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.60, 0.00], [0.00, 0.00, 0.00]]),
                    (16, [[0.00, 0.00, -28.00], [0.00, 0.00, -4.00], [4.65, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -2.60, 0.00], [0.00, 0.00, 0.00]])]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_uterus_term('right uterine horn')[0],
                    'ontId': get_uterus_term('right uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_uterus_term('left uterine horn')[0],
                    'ontId': get_uterus_term('left uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9',
                    'name': get_uterus_term('body of uterus')[0],
                    'ontId': get_uterus_term('body of uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '10',
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-15',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                }]
        }),
        'Material': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3, 4-5-3.2, 3.3-6, 6-7, 7-8-9"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[2.00, 0.00, 0.00], [-1.00, 0.00, 0.00], [0.00, 0.00, -0.30], [0.00, 0.00, 0.00], [0.00, -0.30, 0.00], [0.00, 0.00, 0.00]]),
                    (2, [[1.00, 0.00, 0.00], [-1.00, 0.00, 0.00], [0.00, 0.00, -0.30], [0.00, 0.00, 0.00], [0.00, -0.30, 0.00], [0.00, 0.00, 0.00]]),
                    (3, [[0.00, 0.00, 0.00], [[-1.00, 0.00, 0.00], [1.00, 0.00, 0.00], [0.00, 0.00, -1.00]],
                         [[0.00, 0.00, -0.30], [0.00, 0.00, 0.30], [0.30, 0.00, 0.00]], [0.00, 0.00, 0.20],
                         [[0.00, -0.30, 0.00], [0.00, -0.30, 0.00], [0.00, -0.30, 0.00]], [0.00, 0.00, 0.00]]),
                    (4, [[-2.00, 0.00, 0.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.30], [0.00, 0.00, 0.20], [0.00, -0.30, 0.00], [0.00, 0.00, 0.00]]),
                    (5, [[-1.00, 0.00, 0.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.30], [0.12, 0.00, -0.12], [0.00, -0.30, 0.00], [0.00, 0.00, 0.00]]),
                    (6, [[0.00, 0.00, -1.00], [0.00, 0.00, -1.00], [0.30, 0.00, 0.00], [0.12, 0.00, -0.12], [0.00, -0.30, 0.00], [0.00, 0.00, 0.00]]),
                    (7, [[0.00, 0.00, -2.00], [0.00, 0.00, -1.00], [0.30, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -0.30, 0.00], [0.00, 0.00, 0.00]]),
                    (8, [[0.00, 0.00, -3.00], [0.00, 0.00, -1.00], [0.30, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -0.30, 0.00], [0.00, 0.00, 0.00]]),
                    (9, [[0.00, 0.00, -4.00], [0.00, 0.00, -1.00], [0.30, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -0.30, 0.00], [0.00, 0.00, 0.00]])]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_uterus_term('right uterine horn')[0],
                    'ontId': get_uterus_term('right uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-4',
                    'name': get_uterus_term('left uterine horn')[0],
                    'ontId': get_uterus_term('left uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5',
                    'name': get_uterus_term('body of uterus')[0],
                    'ontId': get_uterus_term('body of uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7-8',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                }]
        })
    }

    @staticmethod
    def getName():
        return '3D Uterus 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Mouse 1',
            'Rat 1',
            'Material']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Mouse 1' in parameterSetName:
            networkLayoutOption = cls.parameterSetStructureStrings['Mouse 1']
        elif 'Rat 1' in parameterSetName:
            networkLayoutOption = cls.parameterSetStructureStrings['Rat 1']
        elif 'Material' in parameterSetName:
            networkLayoutOption = cls.parameterSetStructureStrings['Material']
        else:
            networkLayoutOption = cls.parameterSetStructureStrings['Mouse 1']
        options = {
            'Network layout': copy.deepcopy(networkLayoutOption),
            'Target element length': 6.0,
            'Number of elements around': 8,
            'Number of elements around horns': 8,
            'Number of elements through wall': 1,
            'Wall thickness': 1.0,
            'Double uterus': False,
            'Use linear through wall': True,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        if 'Rat' in parameterSetName:
            options['Target element length'] = 5.0
            options['Wall thickness'] = 1.2
            options['Double uterus'] = True
        if 'Material' in parameterSetName:
            options['Target element length'] = 1.0
            options['Wall thickness'] = 0.1
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Network layout',
            'Target element length',
            'Number of elements around',
            'Number of elements around horns',
            'Number of elements through wall',
            'Wall thickness',
            'Double uterus',
            'Refine',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall']
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Network layout':
            return [MeshType_1d_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Network layout':
            return list(cls.parameterSetStructureStrings.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        """
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        """
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Network layout':
            if not parameterSetName:
                parameterSetName = list(cls.parameterSetStructureStrings.keys())[1]
            return copy.deepcopy(cls.parameterSetStructureStrings[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Network layout'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Network layout'):
            options['Network layout'] = cls.getOptionScaffoldPackage('Network layout', MeshType_1d_network_layout1)
        for key in [
            'Number of elements around',
            'Number of elements through wall',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall'
        ]:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements around'] % 2 != 0:
            options['Number of elements around'] += 1
        if options['Number of elements around'] < options['Number of elements around horns']:
            options['Number of elements around'] = options['Number of elements around horns']

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        networkLayout = options['Network layout']
        elementsCountAround = options['Number of elements around']
        elementsCountAroundHorn = options['Number of elements around horns']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        doubleUterus = options['Double uterus']
        targetElementLength = options['Target element length']
        useCrossDerivatives = options['Use cross derivatives']

        elementsCountAroundRightHorn = elementsCountAroundHorn
        elementsCountAroundLeftHorn = elementsCountAroundHorn

        materialNetworkLayout = cls.parameterSetStructureStrings['Material']
        materialWallThickness = 0.1
        materialTargetElementLength = 0.5

        # Geometric coordinates
        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        geometricNetworkLayout = UterusNetworkLayout(region, networkLayout, targetElementLength)
        rightHornLength = geometricNetworkLayout.arcLengthOfGroupsAlong[0]
        leftHornLength = geometricNetworkLayout.arcLengthOfGroupsAlong[1]
        bodyLength = geometricNetworkLayout.arcLengthOfGroupsAlong[2]
        cervixLength = geometricNetworkLayout.arcLengthOfGroupsAlong[3]
        vaginaLength = geometricNetworkLayout.arcLengthOfGroupsAlong[4]

        elementsCountInRightHorn = math.ceil(rightHornLength / targetElementLength)
        elementsCountInLeftHorn = math.ceil(leftHornLength / targetElementLength)
        elementsCountInBody = math.ceil(bodyLength / targetElementLength)
        elementsCountInCervix = math.ceil(cervixLength / targetElementLength)
        elementsCountInVagina = math.ceil(vaginaLength / targetElementLength)

        nodeIdentifier, elementIdentifier, annotationGroups = \
            createUterusMesh3D(region, fm, coordinates, geometricNetworkLayout, elementsCountAround,
                               elementsCountAroundRightHorn, elementsCountAroundLeftHorn, elementsCountThroughWall,
                               elementsCountInRightHorn, elementsCountInLeftHorn, elementsCountInBody,
                               elementsCountInCervix, elementsCountInVagina, wallThickness, doubleUterus,
                               useCrossDerivatives)

        # # Material coordinates
        # tmp_region = region.createRegion()
        # tmp_fm = tmp_region.getFieldmodule()
        # with ChangeManager(tmp_fm):
        #     tmp_uterus_coordinates = findOrCreateFieldCoordinates(tmp_fm, name="uterus coordinates")
        #     materialNetworkLayout = UterusNetworkLayout(tmp_region, materialNetworkLayout, materialTargetElementLength)
        #
        #     nodeIdentifier, elementIdentifier, materialAnnotationGroups = \
        #         createUterusMesh3D(tmp_region, tmp_fm, tmp_uterus_coordinates, materialNetworkLayout,
        #                            elementsCountAround, elementsCountAroundRightHorn, elementsCountAroundLeftHorn,
        #                            elementsCountThroughWall, elementsCountInRightHorn, elementsCountInLeftHorn,
        #                            elementsCountInBody, elementsCountInCervix, elementsCountInVagina, wallThickness,
        #                            doubleUterus, useCrossDerivatives)
        #     # Write two coordinates
        #     sir = tmp_region.createStreaminformationRegion()
        #     srm = sir.createStreamresourceMemory()
        #     tmp_region.write(sir)
        #     result, buffer = srm.getBuffer()
        #
        #     sir = region.createStreaminformationRegion()
        #     srm = sir.createStreamresourceMemoryBuffer(buffer)
        #     region.read(sir)
        #
        #     del srm
        #     del sir
        #     del tmp_uterus_coordinates
        # del tmp_fm
        # del tmp_region

        return annotationGroups, None

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return

    @classmethod
    def defineFaceAnnotations(cls, region, options, annotationGroups):
        """
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        """
        # Create 2d surface mesh groups
        fm = region.getFieldmodule()
        cervixGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("uterine cervix"))
        rightHornGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("right uterine horn"))
        leftHornGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("left uterine horn"))
        uterusGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("uterus"))

        mesh2d = fm.findMeshByDimension(2)

        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))

        is_cervix = cervixGroup.getGroup()
        is_cervix_outer = fm.createFieldAnd(is_cervix, is_exterior_face_xi3_1)
        is_cervix_inner = fm.createFieldAnd(is_cervix, is_exterior_face_xi3_0)

        is_rightHorn = rightHornGroup.getGroup()
        is_rightHorn_outer = fm.createFieldAnd(is_rightHorn, is_exterior_face_xi3_1)
        is_rightHorn_inner = fm.createFieldAnd(is_rightHorn, is_exterior_face_xi3_0)

        is_leftHorn = leftHornGroup.getGroup()
        is_leftHorn_outer = fm.createFieldAnd(is_leftHorn, is_exterior_face_xi3_1)
        is_leftHorn_inner = fm.createFieldAnd(is_leftHorn, is_exterior_face_xi3_0)

        is_uterus = uterusGroup.getGroup()
        is_uterus_outer = fm.createFieldAnd(is_uterus, is_exterior_face_xi3_1)
        is_uterus_inner = fm.createFieldAnd(is_uterus, is_exterior_face_xi3_0)

        serosaOfCervix = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of uterine cervix"))
        serosaOfCervix.getMeshGroup(mesh2d).addElementsConditional(is_cervix_outer)

        lumenOfCervix = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of uterine cervix"))
        lumenOfCervix.getMeshGroup(mesh2d).addElementsConditional(is_cervix_inner)

        serosaOfRightHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_uterus_term("serosa of right horn"))
        serosaOfRightHorn.getMeshGroup(mesh2d).addElementsConditional(is_rightHorn_outer)

        lumenOfRightHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_uterus_term("lumen of right horn"))
        lumenOfRightHorn.getMeshGroup(mesh2d).addElementsConditional(is_rightHorn_inner)

        serosaOfLeftHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_uterus_term("serosa of left horn"))
        serosaOfLeftHorn.getMeshGroup(mesh2d).addElementsConditional(is_leftHorn_outer)

        lumenOfLeftHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                             get_uterus_term("lumen of left horn"))
        lumenOfLeftHorn.getMeshGroup(mesh2d).addElementsConditional(is_leftHorn_inner)

        serosaOfUterus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of uterus"))
        serosaOfUterus.getMeshGroup(mesh2d).addElementsConditional(is_uterus_outer)

        lumenOfUterus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of uterus"))
        lumenOfUterus.getMeshGroup(mesh2d).addElementsConditional(is_uterus_inner)


def findDerivativeBetweenPoints(v1, v2):
    """
    Find vector difference between two points and rescale vector difference using cubic hermite arclength
    between the points to derive the derivative between the points.
    :param v1: start vector
    :param v2: end vector
    :return: derivative of between v1 and v2
    """
    d = [v2[c] - v1[c] for c in range(3)]
    arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
    d = [c * arcLengthAround for c in normalize(d)]

    return d


class UterusNetworkLayout:
    """
    Generates network layout for uterus scaffold.
    """
    def __init__(self, region, networkLayout, targetElementLength):
        """
        :param region: Zinc region to define model in.
        :param networkLayout: Network layout subscaffold from MeshType_1d_network_layout1
        """

        layoutRegion = region.createRegion()
        layoutFieldmodule = layoutRegion.getFieldmodule()
        layoutNodes = layoutFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # layoutMesh = layoutFieldmodule.findMeshByDimension(1)
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        # layoutAnnotationGroups = networkLayout.getAnnotationGroups()
        layoutCoordinates = findOrCreateFieldCoordinates(layoutFieldmodule)
        layoutFieldcache = layoutFieldmodule.createFieldcache()

        networkMesh = networkLayout.getConstructionObject()

        networkSegments = networkMesh.getNetworkSegments()

        arcLengthOfGroupsAlong = []
        elementsCountAlongList = []
        cxGroups = []
        sxGroups = []
        for n1 in range(len(networkSegments)):
            networkSegment = networkSegments[n1]
            segmentNodes = networkSegment.getNetworkNodes()
            segmentNodeCount = len(segmentNodes)
            for n in range(segmentNodeCount):
                segmentNode = segmentNodes[n]
                layoutNodeIdentifier = segmentNode.getNodeIdentifier()
                layoutNode = layoutNodes.findNodeByIdentifier(layoutNodeIdentifier)
                layoutFieldcache.setNode(layoutNode)
                cx, cd1, cd2, cd3 = get_nodeset_path_ordered_field_parameters(
                    layoutNodes, layoutCoordinates, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                     Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3],
                    networkSegments[n1].getNodeIdentifiers(), networkSegments[n1].getNodeVersions())
            cxGroup = [cx, cd1, cd2, [], cd3, []]
            arcLength = 0.0
            for e in range(len(cx) - 1):
                arcLength += interp.getCubicHermiteArcLength(cx[e], cd1[e],
                                                             cx[e + 1], cd1[e + 1])
            arcLengthOfGroupsAlong.append(arcLength)
            elementsAlong = math.ceil(arcLength / targetElementLength)
            elementsCountAlongList.append(elementsAlong)
            cxGroups.append(cxGroup)
            # Sample
            sx, sd1, pe, pxi, psf = interp.sampleCubicHermiteCurves(cx, cd1, elementsAlong)
            # sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, pe, pxi, psf)
            # sd3, sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, pe, pxi, psf)
            sxGroup = [sx, sd1]
            sxGroups.append(sxGroup)

        del layoutCoordinates
        del layoutNodes
        del layoutFieldmodule
        del layoutRegion

        self.cxGroups = cxGroups
        self.sxGroups = sxGroups
        self.arcLengthOfGroupsAlong = arcLengthOfGroupsAlong
        self.elementsCountAlongList = elementsCountAlongList


def createUterusMesh3D(region, fm, coordinates, geometricNetworkLayout, elementsCountAround,
                       elementsCountAroundRightHorn, elementsCountAroundLeftHorn, elementsCountThroughWall,
                       elementsCountInRightHorn, elementsCountInLeftHorn, elementsCountInBody, elementsCountInCervix,
                       elementsCountInVagina, wallThickness, doubleUterus, useCrossDerivatives):

    mesh = fm.findMeshByDimension(3)

    firstNodeIdentifier = 1
    firstElementIdentifier = 1

    cx_right_horn_group = geometricNetworkLayout.cxGroups[0]
    cx_left_horn_group = geometricNetworkLayout.cxGroups[1]
    cx_body_group = geometricNetworkLayout.cxGroups[2]
    cx_cervix_group = geometricNetworkLayout.cxGroups[3]
    cx_vagina_group = geometricNetworkLayout.cxGroups[4]

    sx_right_horn_group = geometricNetworkLayout.sxGroups[0]
    sx_left_horn_group = geometricNetworkLayout.sxGroups[1]
    sx_body_group = geometricNetworkLayout.sxGroups[2]
    sx_cervix_group = geometricNetworkLayout.sxGroups[3]

    # Create annotation groups
    rightHornGroup = AnnotationGroup(region, get_uterus_term("right uterine horn"))
    leftHornGroup = AnnotationGroup(region, get_uterus_term("left uterine horn"))
    bodyGroup = AnnotationGroup(region, get_uterus_term("body of uterus"))
    cervixGroup = AnnotationGroup(region, get_uterus_term("uterine cervix"))
    vaginaGroup = AnnotationGroup(region, get_uterus_term("vagina"))
    uterusGroup = AnnotationGroup(region, get_uterus_term("uterus"))
    annotationGroups = [cervixGroup, bodyGroup, leftHornGroup, rightHornGroup, uterusGroup, vaginaGroup]

    rightHornMeshGroup = rightHornGroup.getMeshGroup(mesh)
    leftHornMeshGroup = leftHornGroup.getMeshGroup(mesh)
    bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
    cervixMeshGroup = cervixGroup.getMeshGroup(mesh)
    vaginaMeshGroup = vaginaGroup.getMeshGroup(mesh)
    uterusMeshGroup = uterusGroup.getMeshGroup(mesh)

    elementsCountAcross = elementsCountAroundRightHorn - elementsCountAround // 2

    # Get right horn nodes
    rightStartRadians = -math.pi * (elementsCountAround / (2 * elementsCountAroundRightHorn))
    rightHornCoordinates = getTubeNodes(cx_right_horn_group, elementsCountAroundRightHorn,
                                        elementsCountInRightHorn, elementsCountThroughWall, wallThickness,
                                        startRadian=rightStartRadians)

    rhLastRingNodeCoordinates = getTargetedRingNodesCoordinates(rightHornCoordinates, elementsCountAroundRightHorn,
                                                                elementsCountInRightHorn, elementsCountThroughWall,
                                                                omitStartRows=0, omitEndRows=1)

    rhLastRingNodeId, nodeCount = getTargetedRingNodesIds(firstNodeIdentifier, elementsCountAroundRightHorn,
                                                          elementsCountInRightHorn, elementsCountThroughWall,
                                                          omitStartRows=0, omitEndRows=1)

    # Get left horn nodes
    leftStartRadians = -math.pi * (elementsCountAcross / elementsCountAroundLeftHorn)
    leftHornCoordinates = getTubeNodes(cx_left_horn_group, elementsCountAroundLeftHorn, elementsCountInLeftHorn,
                                       elementsCountThroughWall, wallThickness, startRadian=leftStartRadians)

    lhLastRingNodeCoordinates = getTargetedRingNodesCoordinates(leftHornCoordinates, elementsCountAroundLeftHorn,
                                                                elementsCountInLeftHorn, elementsCountThroughWall,
                                                                omitStartRows=0, omitEndRows=1)

    lhLastRingNodeId, nodeCount = getTargetedRingNodesIds(nodeCount, elementsCountAroundLeftHorn,
                                                          elementsCountInLeftHorn, elementsCountThroughWall,
                                                          omitStartRows=0, omitEndRows=1)

    if doubleUterus:
        # Get body nodes
        bodyInnerRightCoordinates, bodyInnerLeftCoordinates, bodyCoordinatesOuter, septumBodyCoordinates = \
            getDoubleTubeNodes(cx_body_group, elementsCountInBody, elementsCountAround, elementsCountAroundRightHorn,
                               elementsCountAroundLeftHorn, elementsCountThroughWall, wallThickness)

        # Get cervix nodes
        rightInnerCervixthroughWall, leftInnerCervixthroughWall, cervixCoordinatesOuter, septumCervixCoordinates = \
            getDoubleTubeNodes(cx_cervix_group, elementsCountInCervix, elementsCountAround,
                               elementsCountAroundRightHorn, elementsCountAroundLeftHorn, elementsCountThroughWall,
                               wallThickness)

        # Get the first layer of inner right and left body coordinates for sampling with the last layer of horns
        xrList = []
        d1rList = []
        d2rList = []
        xlList = []
        d1lList = []
        d2lList = []
        for n2 in range(elementsCountInBody + 1):
            for n3 in range(elementsCountThroughWall + 1):
                for n1 in range(elementsCountAroundRightHorn):
                    x = bodyInnerRightCoordinates[0][n2][n3][n1]
                    d1 = bodyInnerRightCoordinates[1][n2][n3][n1]
                    d2 = bodyInnerRightCoordinates[2][n2][n3][n1]
                    xrList.append(x)
                    d1rList.append(d1)
                    d2rList.append(d2)
                for n1 in range(elementsCountAroundLeftHorn):
                    x = bodyInnerLeftCoordinates[0][n2][n3][n1]
                    d1 = bodyInnerLeftCoordinates[1][n2][n3][n1]
                    d2 = bodyInnerLeftCoordinates[2][n2][n3][n1]
                    xlList.append(x)
                    d1lList.append(d1)
                    d2lList.append(d2)
        newbodyInnerRightCoordinates = [xrList, d1rList, d2rList]
        newbodyInnerLeftCoordinates = [xlList, d1lList, d2lList]

        brFirstRingNodeCoordinates = \
            getTargetedRingNodesCoordinates(newbodyInnerRightCoordinates, elementsCountAroundRightHorn,
                                            elementsCountInBody, elementsCountThroughWall, omitStartRows=1,
                                            omitEndRows=0)
        blFirstRingNodeCoordinates = \
            getTargetedRingNodesCoordinates(newbodyInnerLeftCoordinates, elementsCountAroundLeftHorn,
                                            elementsCountInBody, elementsCountThroughWall, omitStartRows=1,
                                            omitEndRows=0)

        # Sample the curve between last horn rings and first body rings to get the inner bifurcation nodes
        xLayersRight = []
        d2LayersRight = []
        xLayersLeft = []
        d2LayersLeft = []
        elementsCountAlongBifurcation = 2
        for n3 in range(elementsCountThroughWall):
            xRawRight = []
            d2RawRight = []
            for n1 in range(elementsCountAroundRightHorn):
                xAlong = []
                d2Along = []
                xAlong.append(rhLastRingNodeCoordinates[0][n3][n1])
                xAlong.append(brFirstRingNodeCoordinates[0][n3][n1])
                d2Along.append(rhLastRingNodeCoordinates[2][n3][n1])
                d2Along.append(brFirstRingNodeCoordinates[2][n3][n1])
                xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                                elementsCountAlongBifurcation,
                                                                                arcLengthDerivatives=True)[0:2]
                xRawRight.append(xSampledAlong)
                d2RawRight.append(d2SampledAlong)
            xLayersRight.append(xRawRight)
            d2LayersRight.append(d2RawRight)
            xRawLeft = []
            d2RawLeft = []
            for n1 in range(elementsCountAroundLeftHorn):
                xAlong = []
                d2Along = []
                xAlong.append(lhLastRingNodeCoordinates[0][n3][n1])
                xAlong.append(blFirstRingNodeCoordinates[0][n3][n1])
                d2Along.append(lhLastRingNodeCoordinates[2][n3][n1])
                d2Along.append(blFirstRingNodeCoordinates[2][n3][n1])
                xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                                elementsCountAlongBifurcation,
                                                                                arcLengthDerivatives=True)[0:2]
                xRawLeft.append(xSampledAlong)
                d2RawLeft.append(d2SampledAlong)
            xLayersLeft.append(xRawLeft)
            d2LayersLeft.append(d2RawLeft)

        # Rearrange x and d2
        xSampledBifurcationRight = []
        d1SampledBifurcationRight = []
        d2SampledBifurcationRight = []
        xSampledBifurcationLeft = []
        d1SampledBifurcationLeft = []
        d2SampledBifurcationLeft = []
        for n2 in range(2 + 1):
            xRawRight = []
            d1RawRight = []
            d2RawRight = []
            xRawLeft = []
            d1RawLeft = []
            d2RawLeft = []
            for n3 in range(elementsCountThroughWall):
                xAroundRight = []
                d1AroundRight = []
                d2AroundRight = []
                for n1 in range(elementsCountAroundRightHorn):
                    x = xLayersRight[n3][n1][n2]
                    d2 = d2LayersRight[n3][n1][n2]
                    xAroundRight.append(x)
                    d2AroundRight.append(d2)
                    # Calculate d1
                    v1 = xLayersRight[n3][n1][n2]
                    v2 = xLayersRight[n3][n1 + 1 if n1 < elementsCountAroundRightHorn - 1 else 0][n2]
                    d1 = findDerivativeBetweenPoints(v1, v2)
                    d1AroundRight.append(d1)
                d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAroundRight, d1AroundRight)
                xRawRight.append(xAroundRight)
                d1RawRight.append(d1Smoothed)
                d2RawRight.append(d2AroundRight)
                xAroundLeft = []
                d1AroundLeft = []
                d2AroundLeft = []
                for n1 in range(elementsCountAroundLeftHorn):
                    x = xLayersLeft[n3][n1][n2]
                    d2 = d2LayersLeft[n3][n1][n2]
                    xAroundLeft.append(x)
                    d2AroundLeft.append(d2)
                    # Calculate d1
                    v1 = xLayersLeft[n3][n1][n2]
                    v2 = xLayersLeft[n3][n1 + 1 if n1 < elementsCountAroundLeftHorn - 1 else 0][n2]
                    d1 = findDerivativeBetweenPoints(v1, v2)
                    d1AroundLeft.append(d1)
                d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAroundLeft, d1AroundLeft)
                xRawLeft.append(xAroundLeft)
                d1RawLeft.append(d1Smoothed)
                d2RawLeft.append(d2AroundLeft)
            xSampledBifurcationRight.append(xRawRight)
            d1SampledBifurcationRight.append(d1RawRight)
            d2SampledBifurcationRight.append(d2RawRight)
            xSampledBifurcationLeft.append(xRawLeft)
            d1SampledBifurcationLeft.append(d1RawLeft)
            d2SampledBifurcationLeft.append(d2RawLeft)

        innerBifurcationRight = [xSampledBifurcationRight, d1SampledBifurcationRight, d2SampledBifurcationRight]
        innerBifurcationLeft = [xSampledBifurcationLeft, d1SampledBifurcationLeft, d2SampledBifurcationLeft]

        # Get bifurcation outer nodes
        paCentre = sx_body_group[0][1]
        c1Centre = sx_right_horn_group[0][-2]
        c2Centre = sx_left_horn_group[0][-2]
        paxList = bodyCoordinatesOuter[0][1]
        pad2 = bodyCoordinatesOuter[2][1]
        pad3 = bodyCoordinatesOuter[3][1]
        c1xList = rhLastRingNodeCoordinates[0][-1]
        c1d2 = rhLastRingNodeCoordinates[2][-1]
        c2xList = lhLastRingNodeCoordinates[0][-1]
        c2d2 = lhLastRingNodeCoordinates[2][-1]
        rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex = \
            find_tube_bifurcation_points_converging(paCentre, paxList, pad2, c1Centre, c1xList, c1d2, c2Centre, c2xList,
                                                    c2d2)
        rod3 = pad3

        # Get d3 for outer bifurcation (for cox nodes)
        cod3 = []
        for n in range(len(cox)):
            v1 = septumBodyCoordinates[0][0][n + 1]
            v2 = cox[n]
            d3 = findDerivativeBetweenPoints(v1, v2)
            cod3.append(d3)

        # # Get vagina nodes
        # vaginaInnerRightCoordinates, vaginaInnerLeftCoordinates, vaginaCoordinatesOuter, septumVaginaCoordinates = \
        #     getDoubleTubeNodes(cx_vagina_group, elementsCountInVagina, elementsCountAround,
        #                        elementsCountAroundRightHorn, elementsCountAroundLeftHorn, elementsCountThroughWall,
        #                        wallThickness)


        # # Create nodes through wall for sampled inner right bifurcation nodes
        # xRaw = []
        # d3Raw = []
        # for n1 in range(elementsCountAroundRightHorn):
        #     xAlong = []
        #     d3Along = []
        #     v1 = xSampledBifurcationRight[1][n1]
        #     v1d1 = d3SampledBifurcationRight[n1]
        #     if n1 <= elementsCountAround // 2:
        #         v2 = rox[n1]
        #         if n1 == 0 or elementsCountAround // 2:
        #             v2d1 = d3SampledBifurcationRight[n1]
        #         else:
        #             v2d1 = rod3[n1]
        #     else:
        #         v2 = cox[n1 - elementsCountAround // 2 - 1]
        #         # v2d1 = cod3[n1 - elementsCountAround // 2 - 1]
        #         v2d1 = d3SampledBifurcationRight[n1]
        #     xAlong.append(v1)
        #     xAlong.append(v2)
        #     d3Along.append(v1d1)
        #     d3Along.append(v2d1)
        #     xSampledWall, d3SampledWall = interp.sampleCubicHermiteCurves(xAlong, d3Along, elementsCountThroughWall,
        #                                                                     arcLengthDerivatives=True)[0:2]
        #     xRaw.append(xSampledWall)
        #     d3Raw.append(d3SampledWall)
        # print('len(xRaw)', len(xRaw))
        # print('len(xRaw[0])', len(xRaw[0]))

    else:
        # Get body nodes
        bodyCoordinates = getTubeNodes(cx_body_group, elementsCountAround, elementsCountInBody,
                                       elementsCountThroughWall, wallThickness, startRadian=-math.pi / 2)

        bFirstRingNodeCoordinates = getTargetedRingNodesCoordinates(bodyCoordinates, elementsCountAround,
                                                                    elementsCountInBody, elementsCountThroughWall,
                                                                    omitStartRows=1, omitEndRows=0)
        # Get cervix nodes
        cervixCoordinates = getTubeNodes(cx_cervix_group, elementsCountAround, elementsCountInCervix,
                                         elementsCountThroughWall, wallThickness, startRadian=-math.pi / 2)

        # Get vagina nodes
        vaginaCoordinates = getTubeNodes(cx_vagina_group, elementsCountAround, elementsCountInVagina,
                                         elementsCountThroughWall, wallThickness, startRadian=-math.pi / 2)

    # Create nodes
    # Create right horn nodes
    nodeIdentifier = createTubeNodes(fm, coordinates, firstNodeIdentifier, rightHornCoordinates,
                                     elementsCountInRightHorn, elementsCountAroundRightHorn, elementsCountThroughWall,
                                     omitStartRows=0, omitEndRows=1)

    # Create left horn nodes
    nodeIdentifier = createTubeNodes(fm, coordinates, nodeIdentifier, leftHornCoordinates, elementsCountInLeftHorn,
                                     elementsCountAroundLeftHorn, elementsCountThroughWall, omitStartRows=0,
                                     omitEndRows=1)

    if doubleUterus:
        # Create bifurcation nodes
        nodeIdentifier, roNodeId, coNodeId, birNodeId, bilNodeId = \
            createDoubleTubeBifurcationNodes(fm, nodeIdentifier, rox, rod1, rod2, rod3, cox, cod1, cod2, cod3,
                                             innerBifurcationRight, innerBifurcationLeft, elementsCountAroundRightHorn,
                                             elementsCountAroundLeftHorn, elementsCountThroughWall)

        # Create body nodes
        nodeIdentifier, bricNodeId, blicNodeId, botNodeId, bsNodeId = \
            createDoubleTubeNodes(fm, nodeIdentifier, bodyInnerRightCoordinates, bodyInnerLeftCoordinates,
                                  bodyCoordinatesOuter, septumBodyCoordinates, elementsCountInBody, elementsCountAround,
                                  elementsCountAcross, elementsCountThroughWall, omitStartRows=1, omitEndRows=1)

        # Create cervix nodes
        nodeIdentifier, cricNodeId, clicNodeId, cotNodeId, csNodeId = \
            createDoubleTubeNodes(fm, nodeIdentifier, rightInnerCervixthroughWall, leftInnerCervixthroughWall,
                                  cervixCoordinatesOuter, septumCervixCoordinates, elementsCountInCervix,
                                  elementsCountAround, elementsCountAcross, elementsCountThroughWall, omitStartRows=0,
                                  omitEndRows=0)  # omitEndRows should be 1 when we add vagina

        # # Create vagina nodes
        # nodeIdentifier, vricNodeId, vlicNodeId, votNodeId, vsNodeId = \
        #     createDoubleTubeNodes(fm, nodeIdentifier, vaginaInnerRightCoordinates, vaginaInnerLeftCoordinates,
        #                           vaginaCoordinatesOuter, septumVaginaCoordinates, elementsCountInVagina,
        #                           elementsCountAround, elementsCountAcross, elementsCountThroughWall, omitStartRows=0,
        #                           omitEndRows=0)

    else:
        # Create bifurcation nodes
        paCentre = sx_cervix_group[0][1]
        c1Centre = sx_right_horn_group[0][-2]
        c2Centre = sx_left_horn_group[0][-2]
        paxList = bFirstRingNodeCoordinates[0]
        pad2 = bFirstRingNodeCoordinates[2]
        c1xList = rhLastRingNodeCoordinates[0]
        c1d2 = rhLastRingNodeCoordinates[2]
        c2xList = lhLastRingNodeCoordinates[0]
        c2d2 = lhLastRingNodeCoordinates[2]
        nodeIdentifier, roNodeId, coNodeId = \
            createTubeBifurcationNodes(fm, coordinates, nodeIdentifier, paCentre, paxList, pad2, c1Centre, c1xList,
                                       c1d2, c2Centre, c2xList, c2d2, elementsCountThroughWall)

        # Create body nodes
        nodeCount = nodeIdentifier
        nodeIdentifier = createTubeNodes(fm, coordinates, nodeIdentifier, bodyCoordinates, elementsCountInBody,
                                         elementsCountAround, elementsCountThroughWall, omitStartRows=1, omitEndRows=0)

        # Create cervix nodes
        nodeIdentifier = createTubeNodes(fm, coordinates, nodeIdentifier, cervixCoordinates, elementsCountInCervix,
                                         elementsCountAround, elementsCountThroughWall, omitStartRows=1, omitEndRows=0)

        # Create vagina nodes
        nodeIdentifier = createTubeNodes(fm, coordinates, nodeIdentifier, vaginaCoordinates, elementsCountInVagina,
                                         elementsCountAround, elementsCountThroughWall, omitStartRows=1, omitEndRows=0)

    # Create elements
    # Create right horn elements
    startNodeId = firstNodeIdentifier
    elementIdentifier = \
        make_tube_elements(fm, coordinates, startNodeId, firstElementIdentifier, elementsCountInRightHorn,
                           elementsCountAroundRightHorn, elementsCountThroughWall, useCrossDerivatives, omitStartRows=0,
                           omitEndRows=1, meshGroups=[rightHornMeshGroup, uterusMeshGroup])

    # Create left horn elements
    startNodeId = rhLastRingNodeId[-1][-1] + 1
    elementIdentifier = make_tube_elements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInLeftHorn,
                                           elementsCountAroundLeftHorn, elementsCountThroughWall, useCrossDerivatives,
                                           omitStartRows=0, omitEndRows=1,
                                           meshGroups=[leftHornMeshGroup, uterusMeshGroup])

    if doubleUterus:
        # Create bifurcation elements
        c1NodeId = rhLastRingNodeId
        c2NodeId = lhLastRingNodeId
        if elementsCountInBody < 2:
            bricNodeId = cricNodeId[:elementsCountThroughWall * elementsCountAroundRightHorn]
            blicNodeId = clicNodeId[:elementsCountThroughWall * elementsCountAroundLeftHorn]
            botNodeId = cotNodeId[:elementsCountAround]
            bsNodeId = csNodeId[:elementsCountAcross - 1]

        elementIdentifier = \
            make_double_tube_bifurcation_elements(fm, coordinates, elementIdentifier, elementsCountThroughWall,
                                                  c1NodeId, c2NodeId, roNodeId, coNodeId, birNodeId, bilNodeId,
                                                  bricNodeId, blicNodeId, botNodeId, bsNodeId,
                                                  meshGroups=[bodyMeshGroup, rightHornMeshGroup, leftHornMeshGroup,
                                                              uterusMeshGroup])

        # Create body elements
        if elementsCountInBody >= 2:
            elementIdentifier = \
                make_double_tube_elements(mesh, coordinates, elementIdentifier, elementsCountInBody - 1,
                                          elementsCountAround, elementsCountAroundRightHorn,
                                          elementsCountAroundLeftHorn, elementsCountThroughWall, bricNodeId, blicNodeId,
                                          botNodeId, bsNodeId, useCrossDerivatives,
                                          meshGroups=[bodyMeshGroup, uterusMeshGroup])

        # Create cervix elements
        elementIdentifier = \
            make_double_tube_elements(mesh, coordinates, elementIdentifier, elementsCountInCervix, elementsCountAround,
                                      elementsCountAroundRightHorn, elementsCountAroundLeftHorn,
                                      elementsCountThroughWall, cricNodeId, clicNodeId, cotNodeId, csNodeId,
                                      useCrossDerivatives, meshGroups=[cervixMeshGroup, uterusMeshGroup])

        # # Create vagina elements
        # elementIdentifier = \
        #     make_double_tube_elements(mesh, coordinates, elementIdentifier, elementsCountInVagina, elementsCountAround,
        #                               elementsCountAroundRightHorn, elementsCountAroundLeftHorn,
        #                               elementsCountThroughWall, vricNodeId, vlicNodeId, votNodeId, vsNodeId,
        #                               useCrossDerivatives, meshGroups=[vaginaMeshGroup])
    else:
        # Create bifurcation elements
        bFirstRingNodeId, nodeCount = getTargetedRingNodesIds(nodeCount, elementsCountAround, elementsCountInBody,
                                                              elementsCountThroughWall, omitStartRows=1, omitEndRows=0)
        paNodeId = bFirstRingNodeId
        c1NodeId = rhLastRingNodeId
        c2NodeId = lhLastRingNodeId
        elementIdentifier = make_tube_bifurcation_elements(fm, coordinates, elementIdentifier, elementsCountThroughWall,
                                                           paNodeId, c1NodeId, c2NodeId, roNodeId, coNodeId,
                                                           meshGroups=[bodyMeshGroup, rightHornMeshGroup,
                                                                       leftHornMeshGroup, uterusMeshGroup])

        # Create body elements
        startNodeId = paNodeId[0][0]
        elementIdentifier = make_tube_elements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInBody,
                                               elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                               omitStartRows=1, omitEndRows=0,
                                               meshGroups=[bodyMeshGroup, uterusMeshGroup])

        # Create cervix elements
        startNodeId = paNodeId[0][0] + (elementsCountInBody - 1) * elementsCountAround * (elementsCountThroughWall + 1)
        elementIdentifier = make_tube_elements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInCervix,
                                               elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                               omitStartRows=0, omitEndRows=0,
                                               meshGroups=[cervixMeshGroup, uterusMeshGroup])

        # Create vagina elements
        startNodeId = startNodeId + elementsCountInCervix * elementsCountAround * (elementsCountThroughWall + 1)
        elementIdentifier = make_tube_elements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInVagina,
                                               elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                               omitStartRows=0, omitEndRows=0, meshGroups=[vaginaMeshGroup])

    return nodeIdentifier, elementIdentifier, annotationGroups


def getTubeNodes(cx_group, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall, wallThickness,
                 startRadian):

    # Create ellipses along tube around the central path
    xEllipsesAlong = []
    d1EllipsesAlong = []
    for n in range(len(cx_group[0])):
        px, pd1 = createEllipsePoints(cx_group[0][n], 2 * math.pi, cx_group[2][n], cx_group[4][n], elementsCountAround,
                                      startRadians=startRadian)
        xEllipsesAlong.append(px)
        d1EllipsesAlong.append(pd1)

    # Find d2
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipsesAlong) - 1):
            v1 = xEllipsesAlong[n2][n1]
            v2 = xEllipsesAlong[n2 + 1][n1]
            d2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v1)
            d2Along.append(d2)
        xAlong.append(xEllipsesAlong[-1][n1])
        d2Along.append(d2)
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
        d2Raw.append(d2Smoothed)

    # Rearrange d2
    d2EllipsesAlong = []
    for n2 in range(len(xEllipsesAlong)):
        d2Around = []
        for n1 in range(elementsCountAround):
            d2 = d2Raw[n1][n2]
            d2Around.append(d2)
        d2EllipsesAlong.append(d2Around)

    # Spread out elements along tube
    xRaw = []
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipsesAlong)):
            xAlong.append(xEllipsesAlong[n2][n1])
            d2Along.append(d2EllipsesAlong[n2][n1])
        xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along, elementsCountAlongTube,
                                                                        arcLengthDerivatives=True)[0:2]
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlong, d2SampledAlong)
        xRaw.append(xSampledAlong)
        d2Raw.append(d2SampledAlong)

    # Rearrange x and d2
    xSampledTube = []
    d1SampledTube = []
    d2SampledTube = []
    for n2 in range(elementsCountAlongTube + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAround):
            x = xRaw[n1][n2]
            d2 = d2Raw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)
            # Calculate d1
            v1 = xRaw[n1][n2]
            v2 = xRaw[n1 + 1 if n1 < elementsCountAround - 1 else 0][n2]
            d1 = findDerivativeBetweenPoints(v1, v2)
            d1Around.append(d1)
        d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        xSampledTube.append(xAround)
        d1SampledTube.append(d1Smoothed)
        d2SampledTube.append(d2Around)

    d3Tube = []
    for n2 in range(elementsCountAlongTube + 1):
        d3Around = []
        for n1 in range(elementsCountAround):
            d3Around.append(normalize(
                cross(normalize(d1SampledTube[n2][n1]), normalize(d2SampledTube[n2][n1]))))
        d3Tube.append(d3Around)

    xInner = []
    d1Inner = []
    d2Inner = []
    d3Inner = []
    for n2 in range(elementsCountAlongTube + 1):
        for n1 in range(elementsCountAround):
            n = n2 * elementsCountAround + n1
            xInner.append(xSampledTube[n2][n1])
            d1Inner.append(d1SampledTube[n2][n1])
            d2Inner.append(d2SampledTube[n2][n1])
            d3Inner.append(d3Tube[n2][n1])

    transitElementList = [0] * elementsCountAround
    relativeThicknessList = []
    xList, d1List, d2List, d3List, curvatureList = \
        tubemesh.extrudeSurfaceCoordinates(xInner, d1Inner, d2Inner, d3Inner,
                                           [wallThickness] * (elementsCountAlongTube + 1), relativeThicknessList,
                                           elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                                           transitElementList, outward=False)[0:5]

    coordinatesList = [xList, d1List, d2List, d3List]

    return coordinatesList


def getTargetedRingNodesCoordinates(tubeCoordinates, elementsCountAround, elementsCountAlongTube,
                                    elementsCountThroughWall, omitStartRows, omitEndRows):

    # Create tube nodes
    lastRingsNodeId = []
    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    d3LastRing = []
    firstRingsNodeId = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    d3FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            xLastRingThroughWall = []
            d1LastRingThroughWall = []
            d2LastRingThroughWall = []
            # d3LastRingThroughWall = []
            firstRingNodeIdThroughWall = []
            xFirstRingThroughWall = []
            d1FirstRingThroughWall = []
            d2FirstRingThroughWall = []
            # d3FirstRingThroughWall = []
            for n1 in range(elementsCountAround):
                n = n2 * elementsCountAround * (elementsCountThroughWall + 1) + n3 * elementsCountAround + n1
                x = tubeCoordinates[0][n]
                d1 = tubeCoordinates[1][n]
                d2 = tubeCoordinates[2][n]
                # d3 = tubeCoordinates[3][n]
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        if n2 == elementsCountAlongTube - 1:
                            xLastRingThroughWall.append(x)
                            d1LastRingThroughWall.append(d1)
                            d2LastRingThroughWall.append(d2)
                            # d3LastRingThroughWall.append(d3)
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        if n2 == 1:
                            xFirstRingThroughWall.append(x)
                            d1FirstRingThroughWall.append(d1)
                            d2FirstRingThroughWall.append(d2)
                            # d3FirstRingThroughWall.append(d3)
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
                    xLastRing.append(xLastRingThroughWall)
                    d1LastRing.append(d1LastRingThroughWall)
                    d2LastRing.append(d2LastRingThroughWall)
                    # d3LastRing.append(d3LastRingThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)
                    xFirstRing.append(xFirstRingThroughWall)
                    d1FirstRing.append(d1FirstRingThroughWall)
                    d2FirstRing.append(d2FirstRingThroughWall)
                    # d3FirstRing.append(d3FirstRingThroughWall)

    if omitStartRows == 1:
        targetedRingCoordinates = [xFirstRing, d1FirstRing, d2FirstRing, d3FirstRing]
    elif omitEndRows == 1:
        targetedRingCoordinates = [xLastRing, d1LastRing, d2LastRing, d3LastRing]
    else:
        targetedRingCoordinates = []

    return targetedRingCoordinates


def getTargetedRingNodesIds(nodeCount, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                            omitStartRows, omitEndRows):

    # Create tube nodes
    lastRingsNodeId = []
    firstRingsNodeId = []
    for n2 in range(elementsCountAlongTube + 1):
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            firstRingNodeIdThroughWall = []
            for n1 in range(elementsCountAround):
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        if n2 == elementsCountAlongTube - 1:
                            lastRingNodeIdThroughWall.append(nodeCount)
                        nodeCount += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        if n2 == 1:
                            firstRingNodeIdThroughWall.append(nodeCount)
                        nodeCount += 1
                else:
                    nodeCount += 1
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)

    if omitStartRows == 1:
        targetedRingNodeId = firstRingsNodeId
    elif omitEndRows == 1:
        targetedRingNodeId = lastRingsNodeId
    else:
        targetedRingNodeId = []

    return targetedRingNodeId, nodeCount


def createTubeNodes(fm, coordinates, nodeIdentifier, tubeCoordinates, elementsCountAlongTube, elementsCountAround,
                    elementsCountThroughWall, omitStartRows, omitEndRows):

    cache = fm.createFieldcache()

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    # Create tube nodes
    lastRingsNodeId = []
    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    d3LastRing = []
    firstRingsNodeId = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    d3FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        for n3 in range(elementsCountThroughWall + 1):
            lastRingNodeIdThroughWall = []
            xLastRingThroughWall = []
            d1LastRingThroughWall = []
            d2LastRingThroughWall = []
            d3LastRingThroughWall = []
            firstRingNodeIdThroughWall = []
            xFirstRingThroughWall = []
            d1FirstRingThroughWall = []
            d2FirstRingThroughWall = []
            d3FirstRingThroughWall = []
            for n1 in range(elementsCountAround):
                n = n2 * elementsCountAround * (elementsCountThroughWall + 1) + n3 * elementsCountAround + n1
                x = tubeCoordinates[0][n]
                d1 = tubeCoordinates[1][n]
                d2 = tubeCoordinates[2][n]
                d3 = tubeCoordinates[3][n]
                if omitEndRows == 1:  # merging to the bifurcation
                    if n2 == elementsCountAlongTube:
                        pass
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        if n2 == elementsCountAlongTube - 1:
                            lastRingNodeIdThroughWall.append(nodeIdentifier)
                            xLastRingThroughWall.append(x)
                            d1LastRingThroughWall.append(d1)
                            d2LastRingThroughWall.append(d2)
                            d3LastRingThroughWall.append(d3)
                        nodeIdentifier += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        node = nodes.createNode(nodeIdentifier, nodetemplate)
                        cache.setNode(node)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                        # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                        if n2 == 1:
                            firstRingNodeIdThroughWall.append(nodeIdentifier)
                            xFirstRingThroughWall.append(x)
                            d1FirstRingThroughWall.append(d1)
                            d2FirstRingThroughWall.append(d2)
                            d3FirstRingThroughWall.append(d3)
                        nodeIdentifier += 1
                else:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
                    nodeIdentifier += 1
            if omitEndRows == 1:
                if n2 == elementsCountAlongTube - 1:
                    lastRingsNodeId.append(lastRingNodeIdThroughWall)
                    xLastRing.append(xLastRingThroughWall)
                    d1LastRing.append(d1LastRingThroughWall)
                    d2LastRing.append(d2LastRingThroughWall)
                    d3LastRing.append(d3LastRingThroughWall)
            elif omitStartRows == 1:
                if n2 == 1:
                    firstRingsNodeId.append(firstRingNodeIdThroughWall)
                    xFirstRing.append(xFirstRingThroughWall)
                    d1FirstRing.append(d1FirstRingThroughWall)
                    d2FirstRing.append(d2FirstRingThroughWall)
                    d3FirstRing.append(d3FirstRingThroughWall)

    return nodeIdentifier


def find_tube_bifurcation_points_converging(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):
    """
    Gets first ring of coordinates and derivatives between parent pa and
    children c1, c2, and over the crotch between c1 and c2.
    :return: rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex
    """
    paCount = len(pax)
    c1Count = len(c1x)
    c2Count = len(c2x)
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts([paCount, c1Count, c2Count])
    # convert to number of nodes, includes both 6-way points
    pac1NodeCount = pac1Count + 1
    pac2NodeCount = pac2Count + 1
    c1c2NodeCount = c1c2Count + 1
    paStartIndex = 0
    c1StartIndex = 0
    c2StartIndex = 0
    pac1x = [None] * pac1NodeCount
    pac1d1 = [None] * pac1NodeCount
    pac1d2 = [None] * pac1NodeCount
    for n in range(pac1NodeCount):
        pan = (paStartIndex + n) % paCount
        c1n = (c1StartIndex + n) % c1Count
        x1, d1, x2, d2 = c1x[c1n], mult(c1d2[c1n], 2.0), pax[pan], mult(pad2[pan], 2.0)
        pac1x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac1d1[n] = [0.0, 0.0, 0.0]
        pac1d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    paStartIndex2 = paStartIndex + pac1Count
    c1StartIndex2 = c1StartIndex + pac1Count
    c2StartIndex2 = c2StartIndex + c1c2Count
    pac2x = [None] * pac2NodeCount
    pac2d1 = [None] * pac2NodeCount
    pac2d2 = [None] * pac2NodeCount
    for n in range(pac2NodeCount):
        pan = (paStartIndex2 + n) % paCount
        c2n = (c2StartIndex2 + n) % c2Count
        x1, d1, x2, d2 = c2x[c2n], mult(c2d2[c2n], 2.0), pax[pan], mult(pad2[pan], 2.0)
        pac2x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        pac2d1[n] = [0.0, 0.0, 0.0]
        pac2d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    c1c2x = [None] * c1c2NodeCount
    c1c2d1 = [None] * c1c2NodeCount
    c1c2d2 = [None] * c1c2NodeCount
    for n in range(c1c2NodeCount):
        c1n = (c1StartIndex2 + n) % c1Count
        c2n = (c2StartIndex2 - n) % c2Count  # note: reversed
        x1, d1, x2, d2 = c1x[c1n], mult(c1d2[c1n], 2.0), c2x[c2n], mult(c2d2[c2n], -2.0)
        c1c2x[n] = interp.interpolateCubicHermite(x1, d1, x2, d2, 0.5)
        c1c2d1[n] = [0.0, 0.0, 0.0]
        c1c2d2[n] = mult(interp.interpolateCubicHermiteDerivative(x1, d1, x2, d2, 0.5), 0.5)
    # get hex triple points
    hex1, hex1d1, hex1d2 = get_bifurcation_triple_point(
        c2x[c1StartIndex], mult(c2d2[c2StartIndex], -1.0),
        c1x[c1StartIndex], mult(c1d2[c1StartIndex], -1.0),
        pax[paStartIndex], pad2[paStartIndex])
    hex2, hex2d1, hex2d2 = get_bifurcation_triple_point(
        c1x[c1StartIndex2], mult(c1d2[c1StartIndex2], -1.0),
        c2x[c2StartIndex2], mult(c2d2[c2StartIndex2], -1.0),
        pax[paStartIndex2], pad2[paStartIndex2])
    # smooth around loops through hex points to get d1
    loop1x = [hex2] + pac2x[1:-1] + [hex1]
    loop1d1 = [[-d for d in hex2d2]] + pac2d1[1:-1] + [hex1d1]
    loop2x = [hex1] + pac1x[1:-1] + [hex2]
    loop2d1 = [[-d for d in hex1d2]] + pac1d1[1:-1] + [hex2d1]
    loop1d1 = interp.smoothCubicHermiteDerivativesLine(loop1x, loop1d1, fixStartDirection=True, fixEndDirection=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    loop2d1 = interp.smoothCubicHermiteDerivativesLine(loop2x, loop2d1, fixStartDirection=True, fixEndDirection=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    # smooth over "crotch" between c1 and c2
    crotchx = [hex2] + c1c2x[1:-1] + [hex1]
    crotchd1 = [add(hex2d1, hex2d2)] + c1c2d1[1:-1] + [[(-hex1d1[c] - hex1d2[c]) for c in range(3)]]
    crotchd1 = interp.smoothCubicHermiteDerivativesLine(crotchx, crotchd1, fixStartDerivative=True,
                                                        fixEndDerivative=True,
                                                        magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    rox = [hex1] + pac1x[1:-1] + [hex2] + pac2x[1:-1]
    rod1 = [hex1d1] + loop2d1[1:-1] + [hex2d1] + loop1d1[1:-1]
    rod2 = [hex1d2] + pac1d2[1:-1] + [hex2d2] + pac2d2[1:-1]
    cox = crotchx[1:-1]
    cod1 = crotchd1[1:-1]
    cod2 = c1c2d2[1:-1]
    return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex


def createTubeBifurcationNodes(fm, coordinates, nodeIdentifier, paCentre, paxList, pad2, c1Centre, c1xList, c1d2,
                               c2Centre, c2xList, c2d2, elementsCountThroughWall):

    cache = fm.createFieldcache()

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    roxList = []
    rod1List = []
    rod2List = []
    coxList = []
    cod1List = []
    cod2List = []
    for n3 in range(elementsCountThroughWall + 1):
        roxOuter, rod1Outer, rod2Outer, coxOuter, cod1Outer, cod2Outer, paStartIndex, c1StartIndex, c2StartIndex = \
            find_tube_bifurcation_points_converging(paCentre, paxList[n3], pad2[n3], c1Centre, c1xList[n3], c1d2[n3],
                                                    c2Centre, c2xList[n3], c2d2[n3])
        roxList.append(roxOuter)
        rod1List.append(rod1Outer)
        rod2List.append(rod2Outer)
        coxList.append(coxOuter)
        cod1List.append(cod1Outer)
        cod2List.append(cod2Outer)

    # Create bifurcation nodes
    roNodeId = []
    coNodeId = []
    for n3 in range(elementsCountThroughWall + 1):
        coNodeIdThroughWall = []
        for n in range(len(coxOuter)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, coxList[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1List[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2List[n3][n])
            coNodeIdThroughWall.append(nodeIdentifier)
            nodeIdentifier += 1
        coNodeId.append(coNodeIdThroughWall)
        roNodeIdThroughWall = []
        for n in range(len(roxOuter)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, roxList[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1List[n3][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2List[n3][n])
            roNodeIdThroughWall.append(nodeIdentifier)
            nodeIdentifier += 1
        roNodeId.append(roNodeIdThroughWall)

    return nodeIdentifier, roNodeId, coNodeId


def make_tube_elements(fm, coordinates, startNodeId, elementIdentifier, elementsCountAlongTube, elementsCountAround,
                       elementsCountThroughWall, useCrossDerivatives, omitStartRows, omitEndRows, startNodes=None,
                       meshGroups=None):

    mesh = fm.findMeshByDimension(3)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    # Create tube elements
    now = elementsCountAround * (elementsCountThroughWall + 1)
    for e2 in range(elementsCountAlongTube - 1 if omitStartRows or omitEndRows == 1 else elementsCountAlongTube):
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                element = mesh.createElement(elementIdentifier, elementtemplate)
                bni11 = e2 * now + e3 * elementsCountAround + e1 + startNodeId
                bni12 = e2 * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + startNodeId
                bni21 = e2 * now + (e3 + 1) * elementsCountAround + e1 + startNodeId
                bni22 = e2 * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + startNodeId
                nodeIdentifiers = [bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now]
                # print('nodeIdentifiers', nodeIdentifiers)
                result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1

    return elementIdentifier


def make_tube_bifurcation_elements(fm, coordinates, elementIdentifier, elementsCountThroughWall, paNodeId, c1NodeId,
                                   c2NodeId, roNodeId, coNodeId, meshGroups=None):

    elementsCountAround = len(paNodeId[0])
    elementsCountAroundRightTube = len(c1NodeId[0])
    elementsCountAroundLeftTube = len(c2NodeId[0])
    pac1Count, pac2Count, c1c2Count = \
        get_tube_bifurcation_connection_elements_counts(
            [elementsCountAround, elementsCountAroundRightTube, elementsCountAroundLeftTube])

    mesh = fm.findMeshByDimension(3)
    eftfactory = eftfactory_bicubichermitelinear(mesh, None)
    eftStd = eftfactory.createEftBasic()

    elementtemplateStd = mesh.createElementtemplate()
    elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplateStd.defineField(coordinates, -1, eftStd)

    elementtemplateMod = mesh.createElementtemplate()
    elementtemplateMod.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    # Right tube part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAroundRightTube):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % elementsCountAroundRightTube]
                bni3 = roNodeId[e3][e1]
                bni4 = roNodeId[e3][(e1 + 1) % elementsCountAroundRightTube]
                bni5 = bni1 + elementsCountAroundRightTube
                bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            elif e1 == elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % elementsCountAroundRightTube]
                bni3 = roNodeId[e3][e1]
                bni4 = coNodeId[e3][e1 - elementsCountAround // 2]
                bni5 = bni1 + elementsCountAroundRightTube
                bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % elementsCountAroundRightTube]
                bni3 = coNodeId[e3][e1 - elementsCountAround // 2 - 1]
                bni5 = bni1 + elementsCountAroundRightTube
                if e1 == elementsCountAroundRightTube - 1:
                    bni4 = roNodeId[e3][0]
                    bni6 = bni5 - elementsCountAroundRightTube + 1
                else:
                    bni4 = bni3 + 1
                    bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 in (0, pac1Count - 1, pac1Count, elementsCountAroundRightTube - 1):
                eft = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                elif e1 == pac1Count - 1:
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == pac1Count:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif e1 == elementsCountAroundRightTube - 1:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft, [1], [])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 1:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # Left tube part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAroundLeftTube):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < c1c2Count:
                bni1 = c2NodeId[e3][e1]
                bni2 = c2NodeId[e3][(e1 + 1) % elementsCountAroundLeftTube]
                if e1 == 0:
                    bni3 = roNodeId[e3][e1]
                    bni4 = coNodeId[e3][-1 - e1]
                else:
                    bni3 = coNodeId[e3][-e1]
                    if e1 == c1c2Count - 1:
                        bni4 = roNodeId[e3][0] + pac1Count
                    else:
                        bni4 = bni3 - 1
                bni5 = bni1 + elementsCountAroundLeftTube
                bni6 = bni2 + elementsCountAroundLeftTube
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni1 = c2NodeId[e3][e1]
                bni2 = c2NodeId[e3][(e1 + 1) % elementsCountAroundLeftTube]
                bni3 = roNodeId[e3][e1 - c1c2Count + elementsCountAround // 2]
                bni4 = roNodeId[e3][(e1 - c1c2Count + elementsCountAround // 2 + 1) % elementsCountAround]
                bni5 = bni1 + elementsCountAroundLeftTube
                bni6 = bni2 + elementsCountAroundLeftTube
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 <= c1c2Count:
                eft = eftfactory.createEftBasic()
                scalefactors = [-1.0]
                setEftScaleFactorIds(eft, [1], [])
                if e1 == 0:
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    scaleEftNodeValueLabels(eft, [4, 8], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                elif e1 < (c1c2Count - 1):
                    scaleEftNodeValueLabels(eft, [3, 4, 7, 8], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                elif e1 == (c1c2Count - 1):
                    scaleEftNodeValueLabels(eft, [3, 7], [Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2], [1])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == c1c2Count:
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft, [3, 7], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            elif e1 == elementsCountAroundLeftTube - 1:
                eft = eftfactory.createEftBasic()
                remapEftNodeValueLabel(eft, [4, 8], Node.VALUE_LABEL_D_DS2,
                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 2:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # parent part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAround):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = roNodeId[e3][e1]
            bni2 = roNodeId[e3][(e1 + 1) % elementsCountAround]
            bni3 = paNodeId[e3][e1]
            bni4 = paNodeId[e3][(e1 + 1) % elementsCountAround]
            bni5 = bni1 + len(roNodeId[0]) + len(coNodeId[0])
            bni6 = bni2 + len(roNodeId[0]) + len(coNodeId[0])
            bni7 = bni3 + elementsCountAround
            bni8 = bni4 + elementsCountAround
            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            # print('nodeIdentifiers', nodeIdentifiers)
            if e1 in (0, pac1Count):
                eft = eftfactory.createEftBasic()
                if e1 in (0, pac1Count):
                    remapEftNodeValueLabel(eft, [1, 5], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft)
                elementtemplate = elementtemplateMod
            element = mesh.createElement(elementIdentifier, elementtemplate)
            result2 = element.setNodesByIdentifier(eft, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 0:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    return elementIdentifier


def getDoubleTubeNodes(cx_tube_group, elementsCountAlong, elementsCountAround, elementsCountAroundRightTube,
                       elementsCountAroundLeftTube, elementsCountThroughWall, wallThickness):
    """
    :return: the coordinates and derivatives of the outer tube, inner right and left tubes and the septum nodes
    """
    elementsCountAcross = elementsCountAroundRightTube - elementsCountAround // 2

    # Get tube right and left paths
    distance = wallThickness / 2
    xrList = []
    d2rList = []
    d3rList = []
    xlList = []
    d2lList = []
    d3lList = []
    for n in range(len(cx_tube_group[0])):
        x = cx_tube_group[0][n]
        v2r = normalize(cx_tube_group[2][n])
        v3r = normalize(cx_tube_group[4][n])
        d2Mag = magnitude(cx_tube_group[2][n])
        r2 = (d2Mag - distance - wallThickness) / 2
        # d3Mag = magnitude(cx_tube_group[4][n])
        # r3 = (d3Mag - wallThickness)
        v_trans = set_magnitude(v2r, distance + r2)

        x_right = [x[c] + v_trans[c] for c in range(3)]
        xrList.append(x_right)
        d2rList.append(set_magnitude(v2r, r2))
        d3rList.append(set_magnitude(v3r, r2))

        x_left = [x[c] - v_trans[c] for c in range(3)]
        v2l = [-v2r[c] for c in range(3)]
        v3l = [-v3r[c] for c in range(3)]
        xlList.append(x_left)
        d2lList.append(set_magnitude(v2l, r2))
        d3lList.append(set_magnitude(v3l, r2))
    cx_tube_group_right = [xrList, cx_tube_group[1], d2rList, [], d3rList]
    cx_tube_group_left = [xlList, cx_tube_group[1], d2lList, [], d3lList]

    # Get right inner tube nodes
    rightStartRadians = -math.pi * (elementsCountAround / (2 * elementsCountAroundRightTube))
    tubeInnerRightCoordinates = findNodesAlongTube2D(cx_tube_group_right, elementsCountAroundRightTube,
                                                     elementsCountAlong, startRadian=rightStartRadians)

    # Get left inner tube nodes
    leftStartRadians = math.pi - math.pi * (elementsCountAcross / elementsCountAroundLeftTube)
    tubeInnerLeftCoordinates = findNodesAlongTube2D(cx_tube_group_left, elementsCountAroundLeftTube,
                                                    elementsCountAlong, startRadian=leftStartRadians)

    # Get outer nodes along tube
    startRadian = -math.pi / 2
    xOuter, d1Outer, d2Outer = findNodesAlongTube2D(cx_tube_group, elementsCountAround, elementsCountAlong, startRadian)
    outerCoordinates = [xOuter, d1Outer, d2Outer]

    # Get coordinates across tube septum, between two inner canals
    xAcrossSeptum = []
    d1AcrossSeptum = []
    for n in range(elementsCountAlong + 1):
        oa = 0
        ob = elementsCountAround // 2
        v1 = xOuter[n][ob]
        v2 = xOuter[n][oa]
        v3 = [v1[c] / 2 + v2[c] / 2 for c in range(3)]
        v1v2 = [v2[c] - v1[c] for c in range(3)]
        nx = [xOuter[n][ob], v3, xOuter[n][oa]]
        nd1 = [[d / elementsCountAcross for d in v1v2], [d / elementsCountAcross for d in v1v2],
               [d / elementsCountAcross for d in v1v2]]
        px, pd1, pe, pxi = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAcross)[0:4]
        xAcrossSeptum.append(px)
        d1AcrossSeptum.append(pd1)

    # Find d2 across tube septum
    d2Raw = []
    for n2 in range(elementsCountAcross + 1):
        xAlongSeptum = []
        d2AlongSeptum = []
        for n1 in range(elementsCountAlong):
            v1 = xAcrossSeptum[n1][n2]
            v2 = xAcrossSeptum[n1 + 1][n2]
            d2 = findDerivativeBetweenPoints(v1, v2)
            xAlongSeptum.append(v1)
            d2AlongSeptum.append(d2)
        xAlongSeptum.append(xAcrossSeptum[-1][n2])
        d2AlongSeptum.append(d2)
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlongSeptum, d2AlongSeptum)
        d2Raw.append(d2Smoothed)

    # Rearrange d2
    d2AcrossSeptum = []
    for n2 in range(elementsCountAlong + 1):
        d2Across = []
        for n1 in range(elementsCountAcross + 1):
            d2 = d2Raw[n1][n2]
            d2Across.append(d2)
        d2AcrossSeptum.append(d2Across)

    septumCoordinates = [xAcrossSeptum, d1AcrossSeptum, d2AcrossSeptum]

    # Find d3 for right inner canal nodes
    d3InnerRight = []
    for n2 in range(0, elementsCountAlong + 1):
        d3Raw = []
        for n1 in range(elementsCountAroundRightTube):
            v1 = tubeInnerRightCoordinates[0][n2][n1]
            if n1 <= elementsCountAround // 2:
                v2 = outerCoordinates[0][n2][n1]
            else:
                v2 = septumCoordinates[0][n2][n1 - elementsCountAround // 2]
            d3 = findDerivativeBetweenPoints(v1, v2)
            d3Raw.append(d3)
        d3InnerRight.append(d3Raw)
    tubeInnerRightCoordinates.append(d3InnerRight)

    # Find d3 for Left inner canal nodes
    d3InnerLeft = []
    for n2 in range(0, elementsCountAlong + 1):
        d3Raw = []
        for n1 in range(elementsCountAroundLeftTube):
            v1 = tubeInnerLeftCoordinates[0][n2][n1]
            if n1 == 0:
                v2 = outerCoordinates[0][n2][n1]
            elif 0 < n1 < elementsCountAcross:
                v2 = septumCoordinates[0][n2][elementsCountAcross - n1]
            elif n1 == elementsCountAcross:
                v2 = outerCoordinates[0][n2][elementsCountAround // 2]
            else:
                v2 = outerCoordinates[0][n2][n1]
            d3 = findDerivativeBetweenPoints(v1, v2)
            d3Raw.append(d3)
        d3InnerLeft.append(d3Raw)
    tubeInnerLeftCoordinates.append(d3InnerLeft)

    # Find d3 for outer nodes
    d3Outer = []
    for n2 in range(0, elementsCountAlong + 1):
        d3Raw = []
        for n1 in range(elementsCountAround):
            if n1 == 0:
                d1 = septumCoordinates[1][n2][n1]
                d3Raw.append(d1)
            elif 0 < n1 < elementsCountAround // 2:
                v1 = tubeInnerRightCoordinates[0][n2][n1]
                v2 = outerCoordinates[0][n2][n1]
                v1v2 = findDerivativeBetweenPoints(v1, v2)
                d3Raw.append(v1v2)
            elif n1 == elementsCountAround // 2:
                d = [-d1[c] for c in range(3)]
                d3Raw.append(d)
            else:
                v1 = tubeInnerLeftCoordinates[0][n2][elementsCountAcross + n1 - elementsCountAround // 2]
                v2 = outerCoordinates[0][n2][n1]
                v1v2 = findDerivativeBetweenPoints(v1, v2)
                d3Raw.append(v1v2)
            d3Outer.append(d3Raw)

    outerCoordinates = [outerCoordinates[0], outerCoordinates[1], outerCoordinates[2], d3Outer]

    # Get all nodes through wall for right and left inner parts
    rightInnerTubethroughWall, leftInnerTubethroughWall = \
        getInnerDoubleTubeCoordinates(tubeInnerRightCoordinates, tubeInnerLeftCoordinates, outerCoordinates,
                                      septumCoordinates, elementsCountAlong, elementsCountAround,
                                      elementsCountAroundRightTube, elementsCountAroundLeftTube,
                                      elementsCountThroughWall)

    return rightInnerTubethroughWall, leftInnerTubethroughWall, outerCoordinates, septumCoordinates


def findNodesAlongTube2D(sx_group, elementsCountAround, elementsCountAlongTube, startRadian):
    """
    Gets the central path nodes and return the coordinates and derivatives of a 2D tube.
    :return: tube2dCoordinates; 2D tube coordinates and derivatives.
    """
    # Create ellipses along tube around the central path
    xEllipsesAlong = []
    d1EllipsesAlong = []
    for n in range(len(sx_group[0])):
        px, pd1 = createEllipsePoints(sx_group[0][n], 2 * math.pi, sx_group[2][n], sx_group[4][n], elementsCountAround,
                                      startRadians=startRadian)
        xEllipsesAlong.append(px)
        d1EllipsesAlong.append(pd1)

    # Find d2
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipsesAlong) - 1):
            v1 = xEllipsesAlong[n2][n1]
            v2 = xEllipsesAlong[n2 + 1][n1]
            d2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v1)
            d2Along.append(d2)
        xAlong.append(xEllipsesAlong[-1][n1])
        d2Along.append(d2)
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
        d2Raw.append(d2Smoothed)

    # Rearrange d2
    d2EllipsesAlong = []
    for n2 in range(len(xEllipsesAlong)):
        d2Around = []
        for n1 in range(elementsCountAround):
            d2 = d2Raw[n1][n2]
            d2Around.append(d2)
        d2EllipsesAlong.append(d2Around)

    # Spread out elements along tube
    xRaw = []
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipsesAlong)):
            xAlong.append(xEllipsesAlong[n2][n1])
            d2Along.append(d2EllipsesAlong[n2][n1])
        xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along, elementsCountAlongTube,
                                                                        arcLengthDerivatives=False)[0:2]
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlong, d2SampledAlong)
        xRaw.append(xSampledAlong)
        d2Raw.append(d2SampledAlong)

    # Rearrange x and d2
    xSampledTube = []
    d1SampledTube = []
    d2SampledTube = []
    for n2 in range(elementsCountAlongTube + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAround):
            x = xRaw[n1][n2]
            d2 = d2Raw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)
            # Calculate d1
            v1 = xRaw[n1][n2]
            v2 = xRaw[n1 + 1 if n1 < elementsCountAround - 1 else 0][n2]
            d1 = findDerivativeBetweenPoints(v1, v2)
            d1Around.append(d1)
        d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        xSampledTube.append(xAround)
        d1SampledTube.append(d1Smoothed)
        d2SampledTube.append(d2Around)

    tube2dCoordinates = [xSampledTube, d1SampledTube, d2SampledTube]
    return tube2dCoordinates


def getInnerDoubleTubeCoordinates(innerRightCoordinates, innerLeftCoordinates, outerCoordinates, septumCoordinates,
                                  elementsCountAlong, elementsCountAround, elementsCountAroundRightTube,
                                  elementsCountAroundLeftTube, elementsCountThroughWall):
    """
    Gets the outer layer, inner left and right layers and septum coordinates.
    :return: the coordinates and derivatives of all nodes through wall for left and right inner tubes.
    """
    elementsCountAcross = elementsCountAroundRightTube - elementsCountAround // 2

    # Get nodes through wall for right and left inner tubes
    xRawRight = []
    d3RawRight = []
    xRawLeft = []
    d3RawLeft = []
    for n2 in range(elementsCountAlong + 1):
        xtRight = []
        d3tRight = []
        for n1 in range(elementsCountAroundRightTube):
            xAlong = []
            d3Along = []
            v1 = innerRightCoordinates[0][n2][n1]
            if n1 <= elementsCountAround // 2:
                v2 = outerCoordinates[0][n2][n1]
            else:
                v2 = septumCoordinates[0][n2][n1 - elementsCountAround // 2]
            v1v2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v1)
            xAlong.append(v2)
            d3Along.append(v1v2)
            d3Along.append(v1v2)
            xSampledWall, d3SampledWall = interp.sampleCubicHermiteCurves(xAlong, d3Along, elementsCountThroughWall,
                                                                          arcLengthDerivatives=True)[0:2]
            xtRight.append(xSampledWall)
            d3tRight.append(d3SampledWall)
        xRawRight.append(xtRight)
        d3RawRight.append(d3tRight)
        xtLeft = []
        d3tLeft = []
        for n1 in range(elementsCountAroundLeftTube):
            xAlongLeft = []
            d3AlongLeft = []
            v1 = innerLeftCoordinates[0][n2][n1]
            if n1 == 0:
                v2 = outerCoordinates[0][n2][n1]
            elif 0 < n1 < elementsCountAcross:
                v2 = septumCoordinates[0][n2][elementsCountAcross - n1]
            else:
                v2 = outerCoordinates[0][n2][n1 - elementsCountAcross + elementsCountAround // 2]
            v1v2 = findDerivativeBetweenPoints(v1, v2)
            xAlongLeft.append(v1)
            xAlongLeft.append(v2)
            d3AlongLeft.append(v1v2)
            d3AlongLeft.append(v1v2)
            xSampledWall, d3SampledWall = interp.sampleCubicHermiteCurves(xAlongLeft, d3AlongLeft,
                                                                          elementsCountThroughWall,
                                                                          arcLengthDerivatives=True)[0:2]
            xtLeft.append(xSampledWall)
            d3tLeft.append(d3SampledWall)
        xRawLeft.append(xtLeft)
        d3RawLeft.append(d3tLeft)

    # Rearrange the nodes
    xWallRight = []
    xWallLeft = []
    for n2 in range(elementsCountAlong + 1):
        xAroundRight = []
        xAroundLeft = []
        for n3 in range(elementsCountThroughWall + 1):
            v1ListRight = []
            v1ListLeft = []
            for n1 in range(elementsCountAroundRightTube):
                v1 = xRawRight[n2][n1][n3]
                v1ListRight.append(v1)
            xAroundRight.append(v1ListRight)
            for n1 in range(elementsCountAroundLeftTube):
                v1 = xRawLeft[n2][n1][n3]
                v1ListLeft.append(v1)
            xAroundLeft.append(v1ListLeft)
        xWallRight.append(xAroundRight)
        xWallLeft.append(xAroundLeft)

    # Find d1 around tube for all nodes
    d1WallRight = []
    d1WallLeft = []
    for n2 in range(elementsCountAlong + 1):
        xRawRight = []
        d1RawRight = []
        xRawLeft = []
        d1RawLeft = []
        for n3 in range(elementsCountThroughWall + 1):
            xAroundRight = []
            d1AroundRight = []
            xAroundLeft = []
            d1AroundLeft = []
            for n1 in range(elementsCountAroundRightTube):
                v1 = xWallRight[n2][n3][n1]
                v2 = xWallRight[n2][n3][(n1 + 1) % elementsCountAroundRightTube]
                d1 = findDerivativeBetweenPoints(v1, v2)
                xAroundRight.append(v1)
                d1AroundRight.append(d1)
            d1SmoothedRight = interp.smoothCubicHermiteDerivativesLoop(xAroundRight, d1AroundRight)
            xRawRight.append(d1AroundRight)
            d1RawRight.append(d1SmoothedRight)
            for n1 in range(elementsCountAroundRightTube):
                v1 = xWallLeft[n2][n3][n1]
                v2 = xWallLeft[n2][n3][(n1 + 1) % elementsCountAroundLeftTube]
                d1 = findDerivativeBetweenPoints(v1, v2)
                xAroundLeft.append(v1)
                d1AroundLeft.append(d1)
            d1SmoothedLeft = interp.smoothCubicHermiteDerivativesLoop(xAroundLeft, d1AroundLeft)
            xRawLeft.append(xAroundLeft)
            d1RawLeft.append(d1SmoothedLeft)
        d1WallRight.append(d1RawRight)
        d1WallLeft.append(d1RawLeft)

    # Find d2 along tube for nodes through wall
    d2NewRight = []
    d2NewLeft = []
    for n1 in range(elementsCountAroundRightTube):
        d2RawRight = []
        d2RawLeft = []
        for n3 in range(elementsCountThroughWall + 1):
            xALongRight = []
            d2AlongRight = []
            xALongLeft = []
            d2AlongLeft = []
            for n2 in range(elementsCountAlong):
                v1 = xWallRight[n2][n3][n1]
                v2 = xWallRight[n2 + 1][n3][n1]
                d2 = findDerivativeBetweenPoints(v1, v2)
                xALongRight.append(v1)
                d2AlongRight.append(d2)
            xALongRight.append(v2)
            d2AlongRight.append(d2)
            d2SmoothedRight = interp.smoothCubicHermiteDerivativesLine(xALongRight, d2AlongRight)
            d2RawRight.append(d2SmoothedRight)
            for n2 in range(elementsCountAlong):
                v1 = xWallLeft[n2][n3][n1]
                v2 = xWallLeft[n2 + 1][n3][n1]
                d2 = findDerivativeBetweenPoints(v1, v2)
                xALongLeft.append(v1)
                d2AlongLeft.append(d2)
            xALongLeft.append(v2)
            d2AlongLeft.append(d2)
            d2SmoothedLeft = interp.smoothCubicHermiteDerivativesLine(xALongLeft, d2AlongLeft)
            d2RawLeft.append(d2SmoothedLeft)
        d2NewRight.append(d2RawRight)
        d2NewLeft.append(d2RawLeft)

    # Rearrange d2
    d2WallRight = []
    d2WallLeft = []
    for n2 in range(elementsCountAlong + 1):
        d2tRight = []
        d2tLeft = []
        for n3 in range(elementsCountThroughWall + 1):
            d2AroundRight = []
            d2AroundLeft = []
            for n1 in range(elementsCountAroundRightTube):
                d2 = d2NewRight[n1][n3][n2]
                d2AroundRight.append(d2)
            d2tRight.append(d2AroundRight)
            for n1 in range(elementsCountAroundLeftTube):
                d2 = d2NewLeft[n1][n3][n2]
                d2AroundLeft.append(d2)
            d2tLeft.append(d2AroundLeft)
        d2WallRight.append(d2tRight)
        d2WallLeft.append(d2tLeft)

    rightInnerTubethroughWall = [xWallRight, d1WallRight, d2WallRight]
    leftInnerTubethroughWall = [xWallLeft, d1WallLeft, d2WallLeft]

    return rightInnerTubethroughWall, leftInnerTubethroughWall


def createDoubleTubeNodes(fm, nodeIdentifier, xInnerRigh, xInnerLeft, xOuter, xAcross, elementsCountAlong,
                          elementsCountAround, elementsCountAcross, elementsCountThroughWall, omitStartRows,
                          omitEndRows):

    """
    Gets coordinates and derivatives of inner canals and outer of a tube and generate the nodes.
    :return: nodeIdentifier, tricNodeId, tlicNodeId, toNodeId, tsNodeId
    """
    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    count = elementsCountAround // 2 + elementsCountAcross
    elementsCountAroundRightHorn = count
    elementsCountAroundLeftHorn = count

    tricNodeId = []  # tube right inner canal node ids
    tlicNodeId = []  # tube left inner canal node ids
    toNodeId = []  # tube outer node ids
    tsNodeId = []   # tube septum node ids
    for n2 in range(1 if omitStartRows == 1 else 0, elementsCountAlong if omitEndRows == 1 else elementsCountAlong + 1):
        for n3 in range(elementsCountThroughWall):
            for n1 in range(elementsCountAroundRightHorn):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xInnerRigh[0][n2][n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xInnerRigh[1][n2][n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xInnerRigh[2][n2][n3][n1])
                tricNodeId.append(nodeIdentifier)
                nodeIdentifier += 1
        for n3 in range(elementsCountThroughWall):
            for n1 in range(elementsCountAroundLeftHorn):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xInnerLeft[0][n2][n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xInnerLeft[1][n2][n3][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xInnerLeft[2][n2][n3][n1])
                tlicNodeId.append(nodeIdentifier)
                nodeIdentifier += 1
        for n1 in range(1, elementsCountAcross):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAcross[0][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xAcross[1][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xAcross[2][n2][n1])
            tsNodeId.append(nodeIdentifier)
            nodeIdentifier += 1
        for n1 in range(elementsCountAround):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOuter[0][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xOuter[1][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xOuter[2][n2][n1])
            if n1 in (0, elementsCountAround // 2):
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, xOuter[3][n2][n1])
            toNodeId.append(nodeIdentifier)
            nodeIdentifier += 1

    return nodeIdentifier, tricNodeId, tlicNodeId, toNodeId, tsNodeId


def createDoubleTubeBifurcationNodes(fm, nodeIdentifier, rox, rod1, rod2, rod3, cox, cod1, cod2, cod3,
                                     innerBifurcationRight, innerBifurcationLeft, elementsCountAroundRightHorn,
                                     elementsCountAroundLeftHorn, elementsCountThroughWall):
    """
    Gets coordinates and derivatives of bifurcation inner and outer and generate the nodes.
    :return: nodeIdentifier, roNodeId, coNodeId, birNodeId, bilNodeId
    """
    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    # Create extra right inner nodes in bifurcation
    birNodeId = []  # bifurcation right inner ring node ids
    for n2 in range(2):
        for n3 in range(elementsCountThroughWall):
            for n1 in range(elementsCountAroundRightHorn):
                if n2 == 1:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, innerBifurcationRight[0][n2][n3][n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, innerBifurcationRight[1][n2][n3][n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, innerBifurcationRight[2][n2][n3][n1])
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3SampledBifurcationRight[n1])
                    birNodeId.append(nodeIdentifier)
                    nodeIdentifier += 1

    # Create extra left inner nodes in bifurcation
    bilNodeId = []  # bifurcation left inner ring node ids
    for n2 in range(2):
        for n3 in range(elementsCountThroughWall):
            for n1 in range(elementsCountAroundLeftHorn):
                if n2 == 1:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, innerBifurcationLeft[0][n2][n3][n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, innerBifurcationLeft[1][n2][n3][n1])
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, innerBifurcationLeft[2][n2][n3][n1])
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3SampledBifurcationRight[n1])
                    bilNodeId.append(nodeIdentifier)
                    nodeIdentifier += 1

    # Create bifurcation nodes
    roNodeId = []
    for n in range(len(rox)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rox[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2[n])
        if n in (0, len(rox)//2):
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rod3[n])
        roNodeId.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1
    coNodeId = []
    for n in range(len(cox)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cox[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, cod3[n])
        coNodeId.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1

    return nodeIdentifier, roNodeId, coNodeId, birNodeId, bilNodeId


def make_double_tube_elements(mesh, coordinates, elementIdentifier, elementsCountAlong, elementsCountAround,
                              elementsCountAroundRightTube, elementsCountAroundLeftTube, elementsCountThroughWall,
                              tricNodeId, tlicNodeId, toNodeId, tsNodeId, useCrossDerivatives, meshGroups=None):

    """
    Gets tube two inner canals' node ides and outer node ids and creates the elements.
    tricNodeId, tlicNodeId: Tube's right and left inner canals node Ids.
    toNodeId: Tube's outer node Ids.
    tsNodeId: Tube's septum node Ids.
    :return: elements of tube with double inner canals.
    """
    eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    elementtemplateMod = mesh.createElementtemplate()
    elementtemplateMod.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    elementsCountAcross = elementsCountAroundRightTube - elementsCountAround // 2

    # Tube elements
    for e3 in range(elementsCountThroughWall):
        for e2 in range(elementsCountAlong):
            # Elements around right canal in tube
            for e1 in range(elementsCountAroundRightTube):
                bni1 = tricNodeId[e2 * elementsCountAroundRightTube * elementsCountThroughWall + e3 * elementsCountAroundRightTube + e1]
                bni2 = tricNodeId[e2 * elementsCountAroundRightTube * elementsCountThroughWall + e3 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundRightTube]
                bni3 = bni1 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                bni4 = bni2 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                if e1 < elementsCountAround // 2:
                    if e3 == elementsCountThroughWall - 1:
                        bni5 = toNodeId[e2 * elementsCountAround + e1]
                    else:
                        bni5 = bni1 + elementsCountAroundRightTube
                    bni6 = bni5 + 1
                    bni7 = bni5 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                    bni8 = bni7 + 1
                elif e1 == elementsCountAround // 2:
                    if e3 == elementsCountThroughWall - 1:
                        bni5 = toNodeId[e2 * elementsCountAround + e1]
                        bni6 = tsNodeId[e2 * (elementsCountAcross - 1)]
                    else:
                        bni5 = bni1 + elementsCountAroundRightTube
                        bni6 = bni5 + 1
                    bni7 = bni5 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                    bni8 = bni6 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                elif elementsCountAround // 2 < e1 < elementsCountAroundRightTube - 1:
                    if e3 == elementsCountThroughWall - 1:
                        bni5 = tsNodeId[e2 * (elementsCountAcross - 1) + e1 - elementsCountAround // 2 - 1]
                        bni6 = bni5 + 1
                    else:
                        bni5 = bni1 + elementsCountAroundRightTube
                        bni6 = bni2 + elementsCountAroundRightTube
                    bni7 = bni5 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                    bni8 = bni6 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                elif e1 == elementsCountAroundRightTube - 1:
                    if e3 == elementsCountThroughWall - 1:
                        bni5 = tsNodeId[e2 * (elementsCountAcross - 1) + elementsCountAcross - 2]
                        bni6 = toNodeId[e2 * elementsCountAround]
                    else:
                        bni5 = bni1 + elementsCountAroundRightTube
                        bni6 = bni2 + elementsCountAroundRightTube
                    bni7 = bni5 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                    bni8 = bni6 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
                # Remap derivatives for elements at the beginning and end of the septum
                if e1 == elementsCountAround // 2 and e3 == elementsCountThroughWall - 1:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif e1 == elementsCountAroundRightTube - 1 and e3 == elementsCountThroughWall - 1:
                    eft1 = eftfactory.createEftBasic()
                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                else:
                    scalefactors = None
                    eft1 = eft
                    elementtemplate1 = elementtemplate
                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = '-'
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1

            # Elements around Left canal in tube
            for e1 in range(elementsCountAroundLeftTube):
                if 0 <= e1 < elementsCountAcross:
                    bni1 = tlicNodeId[e2 * elementsCountAroundLeftTube * elementsCountThroughWall + e3 * elementsCountAroundRightTube + e1]
                    bni2 = tlicNodeId[e2 * elementsCountAroundLeftTube * elementsCountThroughWall + e3 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundLeftTube]
                    bni3 = bni1 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                    bni4 = bni2 + 2 * elementsCountAroundRightTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                    if e1 == 0:
                        if e3 == elementsCountThroughWall - 1:
                            bni5 = toNodeId[e2 * elementsCountAround + e1]
                            bni6 = tsNodeId[e2 * (elementsCountAcross - 1) + elementsCountAcross - 2]
                        else:
                            bni5 = bni1 + elementsCountAroundLeftTube
                            bni6 = bni2 + elementsCountAroundLeftTube
                    elif e1 == elementsCountAcross - 1:
                        if e3 == elementsCountThroughWall - 1:
                            bni5 = tsNodeId[e2 * (elementsCountAcross - 1) + elementsCountAcross - 1 - e1]
                            bni6 = toNodeId[e2 * elementsCountAround + elementsCountAround // 2]
                        else:
                            bni5 = bni1 + elementsCountAroundLeftTube
                            bni6 = bni2 + elementsCountAroundLeftTube
                    else:
                        if e3 == elementsCountThroughWall - 1:
                            bni5 = tsNodeId[e2 * (elementsCountAcross - 1) + elementsCountAcross - 1 - e1]
                            bni6 = bni5 - 1
                        else:
                            bni5 = bni1 + elementsCountAroundLeftTube
                            bni6 = bni2 + elementsCountAroundLeftTube
                    bni7 = bni5 + 2 * elementsCountAroundLeftTube * elementsCountThroughWall + elementsCountAround+ elementsCountAcross - 1
                    bni8 = bni6 + 2 * elementsCountAroundLeftTube * elementsCountThroughWall + elementsCountAround+ elementsCountAcross - 1
                    nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                    # print('nodeIdentifiers', nodeIdentifiers)
                else:
                    bni1 = tlicNodeId[e2 * elementsCountAroundLeftTube * elementsCountThroughWall + e3 * elementsCountAroundLeftTube + e1]
                    bni2 = tlicNodeId[e2 * elementsCountAroundLeftTube * elementsCountThroughWall + e3 * elementsCountAroundLeftTube + (e1 + 1) % elementsCountAroundLeftTube]
                    bni3 = bni1 + 2 * elementsCountAroundLeftTube * elementsCountThroughWall + elementsCountAround + elementsCountAcross - 1
                    bni4 = bni2 + 2 * elementsCountAroundLeftTube * elementsCountThroughWall + elementsCountAround+ elementsCountAcross - 1
                    if e3 == elementsCountThroughWall - 1:
                        bni5 = toNodeId[e2 * elementsCountAround + e1 - elementsCountAcross + elementsCountAround // 2]
                        if e1 == elementsCountAroundLeftTube - 1:
                            bni6 = toNodeId[e2 * elementsCountAround]
                        else:
                            bni6 = toNodeId[e2 * elementsCountAround + (
                                        e1 + 1) % elementsCountAroundLeftTube - elementsCountAcross + elementsCountAround // 2]
                    else:
                        bni5 = bni1 + elementsCountAroundLeftTube
                        bni6 = bni2 + elementsCountAroundLeftTube
                    bni7 = bni5 + 2 * elementsCountAroundLeftTube * elementsCountThroughWall + elementsCountAround+ elementsCountAcross - 1
                    bni8 = bni6 + 2 * elementsCountAroundLeftTube * elementsCountThroughWall + elementsCountAround+ elementsCountAcross - 1
                    nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                    # print('nodeIdentifiers', nodeIdentifiers)
                # Remap derivatives for elements across septum
                if e3 == elementsCountThroughWall - 1 and e1 == 0:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif e3 == elementsCountThroughWall - 1 and 0 < e1 < elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif e3 == elementsCountThroughWall - 1 and e1 == elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                else:
                    scalefactors = None
                    eft1 = eft
                    elementtemplate1 = elementtemplate

                element = mesh.createElement(elementIdentifier, elementtemplate1)
                result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                if scalefactors:
                    result3 = element.setScaleFactors(eft1, scalefactors)
                else:
                    result3 = '-'
                for meshGroup in meshGroups:
                    meshGroup.addElement(element)
                elementIdentifier += 1
    return elementIdentifier


def make_double_tube_bifurcation_elements(fm, coordinates, elementIdentifier, elementsCountThroughWall, c1NodeId,
                                          c2NodeId, roNodeId, coNodeId, birNodeId, bilNodeId, pricNodeId, plicNodeId,
                                          poNodeId, psNodeId, meshGroups=None):
    """
    Gets child1, child2, parent, and inner bifurcation node Ids.
    c1NodeId, c2NodeId: Child1 and child2 node Ids.
    roNodeId, coNodeId: 2D bifurcation node Ids which are in outer bifurcation surface.
    birNodeId, bilNodeId: Right and left bifurcation inner node Ids.
    pricNodeId, plicNodeId: Parent's right and left inner canals node Ids.
    poNodeId: Parent's outer node Ids.
    psNodeId: parent's septum node Ids.
    :return: elements of bifurcation with double inner tube.
    """
    elementsCountAround = len(poNodeId)
    elementsCountAroundRightTube = len(c1NodeId[0])
    elementsCountAroundLeftTube = len(c2NodeId[0])
    pac1Count, pac2Count, c1c2Count = \
        get_tube_bifurcation_connection_elements_counts(
            [elementsCountAround, elementsCountAroundRightTube, elementsCountAroundLeftTube])

    elementsCountAcross = elementsCountAroundRightTube - elementsCountAround // 2

    mesh = fm.findMeshByDimension(3)
    eftfactory = eftfactory_bicubichermitelinear(mesh, None)
    eftStd = eftfactory.createEftBasic()

    elementtemplateStd = mesh.createElementtemplate()
    elementtemplateStd.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplateStd.defineField(coordinates, -1, eftStd)

    elementtemplateMod = mesh.createElementtemplate()
    elementtemplateMod.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    # Right tube part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAroundRightTube):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % elementsCountAroundRightTube]
                bni3 = birNodeId[e3 * elementsCountAroundRightTube + e1]
                bni4 = birNodeId[e3 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundRightTube]
                bni5 = bni1 + elementsCountAroundRightTube
                bni6 = bni5 + 1
                if e3 == elementsCountThroughWall - 1:
                    bni7 = roNodeId[e1]
                    bni8 = bni7 + 1
                else:
                    bni7 = bni3 + elementsCountAroundRightTube
                    bni8 = bni4 + elementsCountAroundRightTube
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            elif e1 == elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % elementsCountAroundRightTube]
                bni3 = birNodeId[e3 * elementsCountAroundRightTube + e1]
                bni4 = birNodeId[e3 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundRightTube]
                bni5 = bni1 + elementsCountAroundRightTube
                bni6 = bni5 + 1
                if e3 == elementsCountThroughWall - 1:
                    bni7 = roNodeId[e1]
                    bni8 = coNodeId[0]
                else:
                    bni7 = bni3 + elementsCountAroundRightTube
                    bni8 = bni4 + elementsCountAroundRightTube
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            else:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % elementsCountAroundRightTube]
                bni3 = birNodeId[e3 * elementsCountAroundRightTube + e1]
                bni5 = bni1 + elementsCountAroundRightTube
                if e3 == elementsCountThroughWall - 1:
                    bni7 = coNodeId[e1 - elementsCountAround // 2 - 1]
                else:
                    bni7 = bni3 + elementsCountAroundRightTube
                if e1 == elementsCountAroundRightTube - 1:
                    bni4 = birNodeId[0 + e3 * elementsCountAroundRightTube]
                    bni6 = bni5 - elementsCountAroundRightTube + 1
                    if e3 == elementsCountThroughWall - 1:
                        bni8 = roNodeId[0]
                    else:
                        bni8 = bni4 + elementsCountAroundRightTube
                else:
                    bni4 = bni3 + 1
                    bni6 = bni5 + 1
                    bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
            if e3 == elementsCountThroughWall - 1 and e1 in (0, pac1Count - 1, pac1Count, elementsCountAroundRightTube - 1):
                eft1 = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                elif e1 == pac1Count - 1:
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == pac1Count:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif e1 == elementsCountAroundRightTube - 1:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateMod
            else:
                scalefactors = None
                eft1 = eft
                elementtemplate1 = elementtemplate
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 1:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # Left tube part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAroundLeftTube):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = c2NodeId[e3][e1]
            bni2 = c2NodeId[e3][(e1 + 1) % elementsCountAroundLeftTube]
            bni3 = bilNodeId[e3 * elementsCountAroundLeftTube + e1]
            bni4 = bilNodeId[e3 * elementsCountAroundLeftTube + (e1 + 1) % elementsCountAroundLeftTube]
            bni5 = bni1 + elementsCountAroundLeftTube
            bni6 = bni2 + elementsCountAroundLeftTube
            if 0 <= e1 < elementsCountAcross:
                if e3 != elementsCountThroughWall - 1:
                    bni7 = bni3 + elementsCountAroundLeftTube
                    bni8 = bni4 + elementsCountAroundLeftTube
                else:
                    if e1 == 0:
                        bni7 = roNodeId[0]
                        bni8 = coNodeId[-1]
                    elif e1 == elementsCountAcross - 1:
                        bni7 = coNodeId[elementsCountAcross - 1 - e1]
                        bni8 = roNodeId[elementsCountAround // 2]
                    else:
                        bni7 = coNodeId[elementsCountAcross - 1 - e1]
                        bni8 = bni7 - 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            else:
                if e3 == elementsCountThroughWall - 1:
                    bni7 = roNodeId[elementsCountAround // 2 + e1 - elementsCountAcross]
                    if e1 == elementsCountAroundLeftTube - 1:
                        bni8 = roNodeId[0]
                    else:
                        bni8 = bni7 + 1
                else:
                    bni7 = bni3 + elementsCountAroundLeftTube
                    bni8 = bni4 + elementsCountAroundLeftTube
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
            if e3 == elementsCountThroughWall - 1 and 0 <= e1 < elementsCountAcross + 1:
                eft1 = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif 0 < e1 < elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif e1 == elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == elementsCountAcross:
                    scalefactors = [-1.0]
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                elementtemplateMod.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateMod
            elif e3 == elementsCountThroughWall - 1 and e1 == elementsCountAroundLeftTube - 1:
                eft1 = eftfactory.createEftBasic()
                remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateMod
            else:
                scalefactors = None
                eft1 = eft
                elementtemplate1 = elementtemplate
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 2:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # parent right part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAroundRightTube):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = birNodeId[e3 * elementsCountAroundRightTube + e1]
            bni2 = birNodeId[e3 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundRightTube]
            bni3 = pricNodeId[e3 * elementsCountAroundRightTube + e1]
            bni4 = pricNodeId[e3 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundRightTube]
            if e3 == elementsCountThroughWall - 1:
                if e1 < elementsCountAround // 2:
                        bni5 = roNodeId[e1]
                        bni6 = roNodeId[(e1 + 1) % elementsCountAround]
                        bni7 = poNodeId[e1]
                        bni8 = poNodeId[(e1 + 1) % elementsCountAround]
                elif e1 == elementsCountAround // 2:
                    bni5 = roNodeId[e1]
                    bni6 = coNodeId[e1 - elementsCountAround // 2]
                    bni7 = poNodeId[e1]
                    bni8 = psNodeId[e1 - elementsCountAround // 2]
                elif elementsCountAround // 2 < e1 < elementsCountAroundRightTube - 1:
                    bni5 = coNodeId[e1 - elementsCountAround // 2 - 1]
                    bni6 = bni5 + 1
                    bni7 = psNodeId[e1 - elementsCountAround // 2 - 1]
                    bni8 = bni7 + 1
                elif e1 == elementsCountAroundRightTube - 1:
                    bni5 = coNodeId[-1]
                    bni6 = roNodeId[0]
                    bni8 = poNodeId[0]
                    bni7 = bni8 - 1
            else:
                bni5 = bni1 + elementsCountAroundRightTube
                bni6 = bni2 + elementsCountAroundRightTube
                bni7 = bni3 + elementsCountAroundRightTube
                bni8 = bni4 + elementsCountAroundRightTube
            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
            if e3 == elementsCountThroughWall - 1:
                if e1 == 0:
                    eft1 = eftfactory.createEftBasic()
                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif e1 == elementsCountAround // 2:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif elementsCountAround // 2 < e1 < elementsCountAroundRightTube - 1:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif e1 == elementsCountAroundRightTube - 1:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                else:
                    scalefactors = None
                    eft1 = eft
                    elementtemplate1 = elementtemplate
            else:
                scalefactors = None
                eft1 = eft
                elementtemplate1 = elementtemplate
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 0:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # parent left part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAroundLeftTube):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = bilNodeId[e3 * elementsCountAroundRightTube + e1]
            bni2 = bilNodeId[e3 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundLeftTube]
            bni3 = plicNodeId[e3 * elementsCountAroundRightTube + e1]
            bni4 = plicNodeId[e3 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundLeftTube]
            if e3 == elementsCountThroughWall - 1:
                if e1 == 0:
                    bni5 = roNodeId[0]
                    bni6 = coNodeId[-1]
                    bni7 = poNodeId[e1]
                    bni8 = psNodeId[e1 + elementsCountAcross - 2]
                elif 0 < e1 < elementsCountAcross - 1:
                    bni5 = coNodeId[-e1 + elementsCountAcross - 1]
                    bni6 = bni5 - 1
                    bni7 = psNodeId[-e1 + elementsCountAcross - 1]
                    bni8 = bni7 - 1
                elif e1 == elementsCountAcross - 1:
                    bni5 = coNodeId[0]
                    bni6 = roNodeId[elementsCountAround // 2]
                    bni7 = psNodeId[0]
                    bni8 = poNodeId[elementsCountAround // 2]
                else:
                    bni5 = roNodeId[elementsCountAround // 2 - elementsCountAcross + e1]
                    bni7 = poNodeId[elementsCountAround // 2 - elementsCountAcross + e1]
                    if e1 == elementsCountAroundLeftTube - 1:
                        bni6 = roNodeId[0]
                        bni8 = poNodeId[0]
                    else:
                        bni6 = bni5 + 1
                        bni8 = bni7 + 1
            else:
                bni5 = bni1 + elementsCountAroundRightTube
                bni6 = bni2 + elementsCountAroundRightTube
                bni7 = bni3 + elementsCountAroundRightTube
                bni8 = bni4 + elementsCountAroundRightTube
            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
            if e3 == elementsCountThroughWall - 1:
                if 0 <= e1 < elementsCountAcross + 1:
                    eft1 = eftfactory.createEftBasic()
                    if e1 == 0:
                        scalefactors = [-1.0]
                        setEftScaleFactorIds(eft1, [1], [])
                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                        remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                        remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    elif 0 < e1 < elementsCountAcross - 1:
                        scalefactors = [-1.0]
                        setEftScaleFactorIds(eft1, [1], [])
                        remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                    elif e1 == elementsCountAcross - 1:
                        scalefactors = [-1.0]
                        setEftScaleFactorIds(eft1, [1], [])
                        remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1])])
                        remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                        remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                    elif e1 == elementsCountAcross:
                        remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                else:
                    scalefactors = None
                    eft1 = eft
                    elementtemplate1 = elementtemplate
            else:
                scalefactors = None
                eft1 = eft
                elementtemplate1 = elementtemplate
            element = mesh.createElement(elementIdentifier, elementtemplate1)
            result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
            if scalefactors:
                result3 = element.setScaleFactors(eft1, scalefactors)
            else:
                result3 = '-'
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 0:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    return elementIdentifier