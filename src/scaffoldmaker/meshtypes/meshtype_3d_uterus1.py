"""
Generates a 3-D uterus mesh from a 1-D network layout, with variable
numbers of elements around, along and through wall.
"""

import copy
import math

from cmlibs.maths.vectorops import add, mult
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, scaleEftNodeValueLabels, setEftScaleFactorIds, \
    remapEftLocalNodes
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm, \
    findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.uterus_terms import get_uterus_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.bifurcation import get_tube_bifurcation_connection_elements_counts, \
    get_bifurcation_triple_point
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters


class MeshType_3d_uterus1(Scaffold_base):
    """
    Generates a 3-D uterus mesh from a 1-D network layout with variable numbers of elements around, along and through
    wall.
    Magnitude of D2 and D3 are the radii of the uterus in the respective directions.
    """
    parameterSetStructureStrings = {
        'Mouse 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5, 6-7-8-9-5.2, 5.3-10-11, 11-12-13-14"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[13.84, 2.07, 20.50], [-0.53, -0.68, -6.63], [1.00, 0.07, -0.09], [-0.02, -0.05, -0.24], [0.08, -0.99, 0.10], [-0.02, -0.01, -0.01]]),
                    (2, [[12.56, 1.40, 14.10], [-2.03, -0.65, -6.11], [0.95, 0.03, -0.32], [-0.08, -0.03, -0.22], [0.06, -1.00, 0.09], [-0.01, -0.00, -0.01]]),
                    (3, [[9.85, 0.78, 8.39], [-3.33, -0.59, -5.33], [0.84, 0.01, -0.53], [-0.15, -0.01, -0.22], [0.06, -1.00, 0.07], [-0.00, -0.00, -0.03]]),
                    (4, [[5.97, 0.24, 3.53], [-4.99, -0.41, -4.33], [0.66, 0.01, -0.76], [-0.24, 0.08, -0.20], [0.05, -1.00, 0.03], [0.01, 0.01, -0.11]]),
                    (5, [[0.00, 0.00, 0.00], [[-6.82, -0.07, -2.67], [6.82, -0.07, -2.67], [0.00, 0.00, -2.00]],
                         [[0.66, 0.01, -0.76], [0.66, -0.01, 0.76], [1.00, 0.00, 0.00]], [-0.10, 0.09, 0.10],
                         [[0.05, -1.00, 0.03], [0.05, -1.00, 0.03], [0.05, -1.00, 0.03]], [-0.01, 0.01, -0.10]]),
                    (6, [[-13.84, 2.07, 20.50], [0.53, -0.68, -6.63], [1.00, -0.04, 0.08], [0.10, -0.04, 0.40], [-0.05, -0.99, 0.10], [-0.03, -0.00, 0.04]]),
                    (7, [[-12.56, 1.40, 14.10], [2.03, -0.65, -6.11], [0.96, -0.03, 0.32], [-0.08, 0.01, 0.22], [-0.06, -0.99, 0.09], [-0.01, -0.01, -0.01]]),
                    (8, [[-9.85, 0.78, 8.39], [3.33, -0.59, -5.33], [0.85, -0.02, 0.53], [-0.15, 0.01, 0.22], [-0.07, -1.00, 0.07], [0.00, -0.01, -0.03]]),
                    (9, [[-5.97, 0.24, 3.53], [4.99, -0.41, -4.33], [0.66, -0.01, 0.76], [0.04, 0.01, -0.20], [-0.05, -1.00, 0.03], [0.03, 0.00, -0.03]]),
                    (10, [[0.00, 0.00, -2.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.07, 0.00, -0.15], [0.00, -1.00, 0.00], [0.01, -0.00, -0.01]]),
                    (11, [[0.00, 0.00, -4.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (12, [[0.00, 0.00, -6.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (13, [[0.00, 0.00, -8.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (14, [[0.00, 0.00, -10.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]])]),
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
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-13',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                }]
        }),
        'Rat 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5, 6-7-8-9-5.2, 5.3-10-11, 11-12-13-14"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[13.84, 2.07, 20.50], [-0.53, -0.68, -6.63], [1.00, 0.07, -0.09], [-0.02, -0.05, -0.24],
                         [0.08, -0.99, 0.10], [-0.02, -0.01, -0.01]]),
                    (2, [[12.56, 1.40, 14.10], [-2.03, -0.65, -6.11], [0.95, 0.03, -0.32], [-0.08, -0.03, -0.22],
                         [0.06, -1.00, 0.09], [-0.01, -0.00, -0.01]]),
                    (3, [[9.85, 0.78, 8.39], [-3.33, -0.59, -5.33], [0.84, 0.01, -0.53], [-0.15, -0.01, -0.22],
                         [0.06, -1.00, 0.07], [-0.00, -0.00, -0.03]]),
                    (4, [[5.97, 0.24, 3.53], [-4.99, -0.41, -4.33], [0.66, 0.01, -0.76], [-0.24, 0.08, -0.20],
                         [0.05, -1.00, 0.03], [0.01, 0.01, -0.11]]),
                    (5, [[0.00, 0.00, 0.00], [[-6.82, -0.07, -2.67], [6.82, -0.07, -2.67], [0.00, 0.00, -2.00]],
                         [[0.66, 0.01, -0.76], [0.66, -0.01, 0.76], [1.00, 0.00, 0.00]], [-0.10, 0.09, 0.10],
                         [[0.05, -1.00, 0.03], [0.05, -1.00, 0.03], [0.05, -1.00, 0.03]], [-0.01, 0.01, -0.10]]),
                    (6, [[-13.84, 2.07, 20.50], [0.53, -0.68, -6.63], [1.00, -0.04, 0.08], [0.10, -0.04, 0.40],
                         [-0.05, -0.99, 0.10], [-0.03, -0.00, 0.04]]),
                    (7, [[-12.56, 1.40, 14.10], [2.03, -0.65, -6.11], [0.96, -0.03, 0.32], [-0.08, 0.01, 0.22],
                         [-0.06, -0.99, 0.09], [-0.01, -0.01, -0.01]]),
                    (8, [[-9.85, 0.78, 8.39], [3.33, -0.59, -5.33], [0.85, -0.02, 0.53], [-0.15, 0.01, 0.22],
                         [-0.07, -1.00, 0.07], [0.00, -0.01, -0.03]]),
                    (9, [[-5.97, 0.24, 3.53], [4.99, -0.41, -4.33], [0.66, -0.01, 0.76], [0.04, 0.01, -0.20],
                         [-0.05, -1.00, 0.03], [0.03, 0.00, -0.03]]),
                    (10, [[0.00, 0.00, -4.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.07, 0.00, -0.15],
                          [0.00, -1.00, 0.00], [0.01, -0.00, -0.01]]),
                    (11, [[0.00, 0.00, -8.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00],
                          [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (12, [[0.00, 0.00, -12.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00],
                          [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (13, [[0.00, 0.00, -16.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00],
                          [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]]),
                    (14, [[0.00, 0.00, -20.00], [0.00, 0.00, -2.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.00],
                          [0.00, -1.00, 0.00], [0.00, 0.00, 0.00]])]),
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
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-13',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                }]
        }),
        'Sheep 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5-6, 7-8-9-10-11-6.2, 6.3-12-13, 13-14-15-16"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[52.49, 0.00, -18.84], [-2.73, 0.00, 11.28], [-2.77, 0.00, -0.67], [0.02, 0.00, -0.43], [0.00, -2.43, 0.00], [0.00, -0.00, 0.00]]),
                    (2, [[48.29, 0.00, -7.35], [-5.65, 0.00, 11.61], [-2.57, 0.00, -1.25], [0.38, 0.00, -0.73], [0.00, -2.43, 0.00], [0.00, -0.00, 0.00]]),
                    (3, [[41.18, 0.00, 4.10], [-9.81, 0.00, 9.08], [-2.00, 0.00, -2.15], [0.95, 0.00, -0.92], [0.00, -2.43, 0.00], [0.00, 0.00, 0.00]]),
                    (4, [[29.62, 0.00, 10.08], [-12.84, 0.00, 2.79], [-0.67, 0.00, -3.09], [1.45, 0.00, -0.46], [0.00, -2.43, 0.00], [0.00, 0.00, 0.00]]),
                    (5, [[16.68, 0.00, 9.44], [-15.56, 0.00, -4.52], [0.89, 0.00, -3.08], [1.50, 0.00, 0.09], [0.00, -2.43, -0.00], [0.00, 0.00, 0.00]]),
                    (6, [[0.00, 0.00, 0.00], [[-17.41, 0.00, -14.04], [17.41, 0.00, -14.04], [0.00, 0.00, -5.00]],
                         [[2.31, 0.00, -2.86], [2.87, 0.00, 3.56], [4.89, 0.00, 0.00]], [-0.24, 0.00, 1.06],
                         [[0.00, -2.43, 0.00], [0.00, -2.43, 0.00], [0.00, -3.50, 0.00]], [0.00, 0.00, 0.00]]),
                    (7, [[-52.49, 0.00, -18.84], [2.73, 0.00, 11.28], [-2.77, 0.00, 0.67], [-0.74, 0.00, 1.11], [0.00, -2.43, 0.00], [0.00, 0.00, 0.00]]),
                    (8, [[-48.29, 0.00, -7.35], [5.65, 0.00, 11.61], [-2.57, 0.00, 1.25], [0.38, 0.00, 0.73], [0.00, -2.43, 0.00], [0.00, 0.00, 0.00]]),
                    (9, [[-41.18, 0.00, 4.10], [9.81, 0.00, 9.08], [-2.00, 0.00, 2.15], [0.95, 0.00, 0.92], [0.00, -2.43, 0.00], [0.00, 0.00, 0.00]]),
                    (10, [[-29.62, 0.00, 10.08], [12.84, 0.00, 2.79], [-0.67, 0.00, 3.09], [1.45, 0.00, 0.46], [0.00, -2.43, 0.00], [0.00, 0.00, 0.00]]),
                    (11, [[-16.68, 0.00, 9.44], [15.56, 0.00, -4.52], [0.89, 0.00, 3.08], [1.91, 0.00, -1.03], [0.00, -2.43, 0.00], [0.00, -0.36, 0.00]]),
                    (12, [[0.00, 0.00, -10.00], [0.00, 0.00, -15.00], [3.50, 0.00, 0.00], [1.12, 0.00, -1.33], [0.00, -3.50, 0.00], [0.00, -0.46, 0.00]]),
                    (13, [[0.00, 0.00, -30.00], [0.00, 0.00, -15.00], [3.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -3.50, 0.00], [0.00, -0.00, 0.00]]),
                    (14, [[0.00, 0.00, -40.00], [0.00, 0.00, -10.00], [3.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -3.50, 0.00], [0.00, -0.00, 0.00]]),
                    (15, [[0.00, 0.00, -50.00], [0.00, 0.00, -10.00], [3.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -3.50, 0.00], [0.00, 0.00, 0.00]]),
                    (16, [[0.00, 0.00, -60.00], [0.00, 0.00, -10.00], [3.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -3.50, 0.00], [0.00, 0.00, 0.00]])]),
                    # (1, [[48.29, 0.00, -7.35], [-4.19, 0.00, 13.21], [-2.72, 0.00, -0.86], [0.42, 0.00, -1.48], [0.00, -2.43, 0.00], [0.00, -0.83, 0.00]]),
                    # (2, [[41.18, 0.00, 4.10], [-9.84, 0.00, 9.10], [-1.99, 0.00, -2.16], [1.03, 0.00, -1.11], [0.00, -2.43, 0.00], [0.00, -0.37, 0.00]]),
                    # (3, [[29.62, 0.00, 10.08], [-12.84, 0.00, 2.79], [-0.67, 0.00, -3.09], [1.44, 0.00, -0.46], [0.00, -2.43, 0.00], [0.00, -0.27, 0.00]]),
                    # (4, [[16.68, 0.00, 9.44], [-15.56, 0.00, -4.52], [0.89, 0.00, -3.08], [1.50, 0.00, 0.09], [0.00, -2.43, -0.00], [0.00, -0.14, 0.00]]),
                    # (5, [[0.00, 0.00, 0.00], [[-17.41, 0.00, -14.04], [17.41, 0.00, -14.04], [0.00, 0.00, -5.00]],
                    #      [[2.31, 0.00, -2.86], [2.87, 0.0, 3.56], [4.89, 0.00, 0.00]], [-0.36, -0.09, 1.19],
                    #      [[0.00, -2.43, 0.00], [0.00, -2.43, 0.00], [0.00, -3.50, 0.00]], [0.07, 0.43, -0.02]]),
                    # (6, [[-48.29, 0.00, -7.35], [4.19, 0.00, 13.21], [-2.72, -0.31, 0.86], [-0.55, 0.17, 1.83], [0.00, -2.43, 0.00], [-0.14, -0.29, 0.05]]),
                    # (7, [[-41.18, 0.00, 4.10], [9.84, 0.00, 9.10], [-1.99, 0.00, 2.16], [1.00, 0.15, 1.12], [0.00, -2.43, 0.00], [-0.12, -0.37, 0.04]]),
                    # (8, [[-29.62, 0.00, 10.08], [12.84, 0.00, 2.79], [-0.67, 0.00, 3.09], [1.44, 0.00, 0.46], [0.00, -2.43, 0.00], [0.00, -0.27, 0.00]]),
                    # (9, [[-16.68, 0.00, 9.44], [15.56, 0.00, -4.52], [0.89, 0.00, 3.08], [1.91, 0.00, -1.03], [0.00, -2.43, 0.00], [0.00, -0.24, 0.00]]),
                    # (10, [[0.00, 0.00, -10.00], [0.00, 0.00, -15.00], [3.50, 0.00, 0.00], [1.12, 0.00, -1.33], [0.00, -3.50, 0.00], [0.00, 0.03, 0.00]]),
                    # (11, [[0.00, 0.00, -30.00], [0.00, 0.00, -15.00], [3.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -3.50, 0.00], [0.00, -0.00, 0.00]]),
                    # (12, [[0.00, 0.00, -40.00], [0.00, 0.00, -10.00], [3.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -3.50, 0.00], [0.00, -0.00, 0.00]]),
                    # (13, [[0.00, 0.00, -50.00], [0.00, 0.00, -10.00], [3.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -3.50, 0.00], [0.00, 0.00, 0.00]]),
                    # (14, [[0.00, 0.00, -60.00], [0.00, 0.00, -10.00], [3.50, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -3.50, 0.00], [0.00, 0.00, 0.00]])]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-5',
                    'name': get_uterus_term('right uterine horn')[0],
                    'ontId': get_uterus_term('right uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-10',
                    'name': get_uterus_term('left uterine horn')[0],
                    'ontId': get_uterus_term('left uterine horn')[1]
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
                    'identifierRanges': '13-15',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                }]
        }),
        'Material': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3, 4-5-3.2, 3.3-6, 6-7-8-9"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    (1, [[2.00, 0.00, 0.00], [-1.00, 0.00, 0.00], [0.00, 0.00, -0.10], [0.02, -0.05, -0.24], [0.00, -0.10, 0.00], [-0.02, -0.01, -0.01]]),
                    (2, [[1.00, 0.00, 0.00], [-1.00, 0.00, 0.00], [0.00, 0.00, -0.10], [-0.08, -0.03, -0.22], [0.00, -0.10, 0.00], [-0.01, -0.00, -0.01]]),
                    (3, [[0.00, 0.00, 0.00], [[-1.00, 0.00, 0.00], [1.00, 0.00, 0.00], [0.00, 0.00, -1.00]],
                         [[0.00, 0.00, -0.10], [0.00, 0.00, 0.10], [0.10, 0.00, 0.00]], [-0.10, 0.09, 0.10],
                         [[0.00, -0.10, 0.00], [0.00, -0.10, 0.00], [0.00, -0.10, 0.00]], [-0.01, 0.01, -0.10]]),
                    (4, [[-2.00, 0.00, 0.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.10], [0.10, -0.04, 0.40], [0.00, -0.10, 0.00], [-0.03, -0.00, 0.04]]),
                    (5, [[-1.00, 0.00, 0.00], [1.00, 0.00, 0.00], [0.00, 0.00, 0.10], [-0.08, 0.01, 0.22], [0.00, -0.10, 0.00], [-0.01, -0.01, -0.01]]),
                    (6, [[0.00, 0.00, -1.00], [0.00, 0.00, -1.00], [0.10, 0.00, 0.00], [-0.15, 0.01, 0.22], [0.00, -0.10, 0.00], [0.00, -0.01, -0.03]]),
                    (7, [[0.00, 0.00, -2.00], [0.00, 0.00, -1.00], [0.10, 0.00, 0.00], [0.04, 0.01, -0.20], [0.00, -0.10, 0.00], [0.03, 0.00, -0.03]]),
                    (8, [[0.00, 0.00, -3.00], [0.00, 0.00, -1.00], [0.10, 0.00, 0.00], [0.07, 0.00, -0.15], [0.00, -0.10, 0.00], [0.01, -0.00, -0.01]]),
                    (9, [[0.00, 0.00, -4.00], [0.00, 0.00, -1.00], [0.10, 0.00, 0.00], [0.00, 0.00, 0.00], [0.00, -0.10, 0.00], [0.00, 0.00, 0.00]])]),
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
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-8',
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
            'Sheep 1',
            'Material']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Mouse 1' in parameterSetName:
            networkLayoutOption = cls.parameterSetStructureStrings['Mouse 1']
        elif 'Rat 1' in parameterSetName:
            networkLayoutOption = cls.parameterSetStructureStrings['Rat 1']
        elif 'Sheep 1' in parameterSetName:
            networkLayoutOption = cls.parameterSetStructureStrings['Sheep 1']
        elif 'Material' in parameterSetName:
            networkLayoutOption = cls.parameterSetStructureStrings['Material']
        else:
            networkLayoutOption = cls.parameterSetStructureStrings['Rat 1']
        options = {
            'Network layout': copy.deepcopy(networkLayoutOption),
            'Target element length': 6.0,
            'Number of elements around': 8,
            # 'Number of elements across': 3,
            'Wall thickness': 2.0,
            'Number of elements through wall': 1,
            'Double uterus': False,
            'Use linear through wall': True,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        if 'Rat' in parameterSetName:
            options['Number of elements through wall'] = 1  # only works for 1
            options['Target element length'] = 6.0
            options['Wall thickness'] = 1.5
            options['Double uterus'] = True
        if 'Sheep' in parameterSetName:
            options['Target element length'] = 11.0
            options['Wall thickness'] = 3.0
        if 'Material' in parameterSetName:
            options['Target element length'] = 0.5
            options['Wall thickness'] = 0.1
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Network layout',
            'Target element length',
            'Number of elements around',
            # 'Number of elements across',
            'Wall thickness',
            'Number of elements through wall',
            'Double uterus',
            'Use linear through wall',
            'Use cross derivatives',
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
        # elementsCountAcross = options['Number of elements across']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        doubleUterus = options['Double uterus']
        targetElementLength = options['Target element length']
        useCrossDerivatives = options['Use cross derivatives']

        elementsCountAcross = elementsCountAround // 2

        materialNetworkLayout = cls.parameterSetStructureStrings['Material']
        materialWallThickness = 0.1
        materialTargetElementLength = 0.5

        # Geometric coordinates
        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        geometricNetworkLayout = UterusNetworkLayout(region, networkLayout, targetElementLength)
        rightHornLength = geometricNetworkLayout.arcLengthOfGroupsAlong[0]
        leftHornLength = geometricNetworkLayout.arcLengthOfGroupsAlong[1]
        cervixLength = geometricNetworkLayout.arcLengthOfGroupsAlong[2]
        vaginaLength = geometricNetworkLayout.arcLengthOfGroupsAlong[3]

        elementsCountInRightHorn = math.ceil(rightHornLength / targetElementLength)
        elementsCountInLeftHorn = math.ceil(leftHornLength / targetElementLength)
        elementsCountInCervix = math.ceil(cervixLength / targetElementLength)
        elementsCountInVagina = math.ceil(vaginaLength / targetElementLength)



        if doubleUterus:
            nodeIdentifier, elementIdentifier, annotationGroups = \
                createUterusMesh3DRat(region, fm, coordinates, geometricNetworkLayout, elementsCountAround, elementsCountAcross,
                                   elementsCountThroughWall, elementsCountInRightHorn, elementsCountInLeftHorn,
                                   elementsCountInCervix, elementsCountInVagina, wallThickness, useCrossDerivatives)
        else:
            nodeIdentifier, elementIdentifier, annotationGroups = \
                createUterusMesh3D(region, fm, coordinates, geometricNetworkLayout, elementsCountAround,
                                   elementsCountThroughWall, elementsCountInRightHorn, elementsCountInLeftHorn,
                                   elementsCountInCervix, elementsCountInVagina, wallThickness, useCrossDerivatives)

        # # Material coordinates
        # tmp_region = region.createRegion()
        # tmp_fm = tmp_region.getFieldmodule()
        # with ChangeManager(tmp_fm):
        #     tmp_uterus_coordinates = findOrCreateFieldCoordinates(tmp_fm, name="uterus coordinates")
        #     materialNetworkLayout = UterusNetworkLayout(tmp_region, materialNetworkLayout, materialTargetElementLength)
        #
        #     nodeIdentifier, elementIdentifier, materialAnnotationGroups = \
        #         createUterusMesh3D(tmp_region, tmp_fm, tmp_uterus_coordinates, materialNetworkLayout,
        #                            elementsCountAround, elementsCountThroughWall, elementsCountInRightHorn,
        #                            elementsCountInLeftHorn, elementsCountInCervix, elementsCountInVagina,
        #                            materialWallThickness, useCrossDerivatives)
        #
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
                                                            get_uterus_term("serosa of uerine cervix"))
        serosaOfCervix.getMeshGroup(mesh2d).addElementsConditional(is_cervix_outer)

        lumenOfCervix = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of uerine cervix"))
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
    d = [c * arcLengthAround for c in vector.normalise(d)]

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


def getCoordinatesAlongTube3D(cx_group, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                              wallThickness, startRadian):

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
            d3Around.append(vector.normalise(
                vector.crossproduct3(vector.normalise(d1SampledTube[n2][n1]), vector.normalise(d2SampledTube[n2][n1]))))
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
        tubemesh.getCoordinatesFromInner(xInner, d1Inner, d2Inner, d3Inner, [wallThickness]*(elementsCountAlongTube+1),
                                         relativeThicknessList, elementsCountAround, elementsCountAlongTube,
                                         elementsCountThroughWall, transitElementList)

    coordinatesList = [xList, d1List, d2List, d3List]

    return coordinatesList


def generateTubeNodes(fm, coordinates, nodeIdentifier, tubeCoordinates, elementsCountAlongTube, elementsCountAround,
                      elementsCountThroughWall, omitStartRows, omitEndRows, startNodes=None):

    cache = fm.createFieldcache()
    # coordinates = findOrCreateFieldCoordinates(fm)

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
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
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
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
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
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)
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


def generateTubeElements(fm, coordinates, startNodeId, elementIdentifier, elementsCountAlongTube, elementsCountAround,
                         elementsCountThroughWall, useCrossDerivatives, omitStartRows, omitEndRows, startNodes=None,
                         meshGroups=None):

    mesh = fm.findMeshByDimension(3)
    # coordinates = findOrCreateFieldCoordinates(fm)

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
    for e3 in range(elementsCountThroughWall):
        for e2 in range(elementsCountAlongTube - 1 if omitStartRows or omitEndRows == 1 else elementsCountAlongTube):
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

    # if omitStartRows == 1 or omitEndRows == 1:
    #     lastNodeId = elementsCountAlongTube * elementsCountAround * (elementsCountThroughWall + 1) + startNodeId
    # else:
    #     lastNodeId = (elementsCountAlongTube + 1) * elementsCountAround * (elementsCountThroughWall + 1) + startNodeId
    return elementIdentifier


def make_tube_bifurcation_elements_3d(fm, coordinates, elementIdentifier, elementsCountAround,
                                      elementsCountThroughWall, paNodeId, c1NodeId, c2NodeId, roNodeId, coNodeId,
                                      meshGroups=None):

    paCount = len(paNodeId[0])
    c1Count = len(c1NodeId[0])
    c2Count = len(c2NodeId[0])
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)

    # fm = region.getFieldmodule()
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
        for e1 in range(c1Count):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = roNodeId[e3][(e1 + 1) % c1Count]
                bni5 = bni1 + elementsCountAround
                bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            elif e1 == elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = coNodeId[e3][e1 - elementsCountAround // 2]
                bni5 = bni1 + elementsCountAround
                bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = coNodeId[e3][e1 - elementsCountAround // 2 - 1]
                bni5 = bni1 + elementsCountAround
                if e1 == c1Count - 1:
                    bni4 = roNodeId[e3][0]
                    bni6 = bni5 - elementsCountAround + 1
                else:
                    bni4 = bni3 + 1
                    bni6 = bni5 + 1
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            if e1 in (0, pac1Count - 1, pac1Count, c1Count - 1):
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
                elif e1 == c1Count - 1:
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
        for e1 in range(c2Count):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround//2:
                bni1 = c2NodeId[e3][e1]
                bni2 = c2NodeId[e3][(e1 + 1) % c1Count]
                if e1 == 0:
                    bni3 = roNodeId[e3][e1]
                    bni4 = coNodeId[e3][-1 - e1]
                else:
                    bni3 = coNodeId[e3][-e1]
                    if e1 == elementsCountAround//2 - 1:
                        bni4 = roNodeId[e3][0] + pac1Count
                    else:
                        bni4 = bni3 - 1
                bni5 = bni1 + elementsCountAround
                bni6 = bni2 + elementsCountAround
                bni7 = bni3 + len(roNodeId[0]) + len(coNodeId[0])
                bni8 = bni4 + len(roNodeId[0]) + len(coNodeId[0])
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            else:
                bni1 = c2NodeId[e3][e1]
                bni2 = c2NodeId[e3][(e1 + 1) % c1Count]
                bni3 = roNodeId[e3][e1]
                bni4 = roNodeId[e3][(e1 + 1) % c2Count]
                bni5 = bni1 + elementsCountAround
                bni6 = bni2 + elementsCountAround
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
            elif e1 == c2Count - 1:
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
        for e1 in range(paCount):
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
    # lastNodeId = startNodeIndex + (elementsCountThroughWall + 1) * (len(roNodeId) + len(coNodeId)) + 1
    # lastNodeId = paNodeId[0][0]
    return elementIdentifier


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
        # if n2 == 0:
        #     startNodeId = nodeIdentifier - 1
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
                        if n2 == elementsCountAlongTube - 1:
                            # lastRingNodeIdThroughWall.append(nodeCount)
                            xLastRingThroughWall.append(x)
                            d1LastRingThroughWall.append(d1)
                            d2LastRingThroughWall.append(d2)
                            d3LastRingThroughWall.append(d3)
                        # nodeCount += 1
                elif omitStartRows == 1:  # diverging from bifurcation
                    if n2 == 0:
                        pass
                    else:
                        if n2 == 1:
                            # firstRingNodeIdThroughWall.append(nodeCount)
                            xFirstRingThroughWall.append(x)
                            d1FirstRingThroughWall.append(d1)
                            d2FirstRingThroughWall.append(d2)
                            d3FirstRingThroughWall.append(d3)
                        # nodeCount += 1
                # else:
                #     nodeCount += 1
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

    if omitStartRows == 1:
        targetedRingCoordinates = [xFirstRing, d1FirstRing, d2FirstRing, d3FirstRing]
    elif omitEndRows == 1:
        targetedRingCoordinates = [xLastRing, d1LastRing, d2LastRing, d3LastRing]
    else:
        targetedRingCoordinates = []

    return targetedRingCoordinates


def getTargetedRingNodesId(nodeCount, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                           omitStartRows, omitEndRows):

    # Create tube nodes
    lastRingsNodeId = []
    firstRingsNodeId = []
    for n2 in range(elementsCountAlongTube + 1):
        # if n2 == 0:
        #     startNodeId = nodeIdentifier - 1
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


def create3dBifurcationNodes(fm, coordinates, nodeIdentifier, paCentre, paxList, pad2, c1Centre, c1xList, c1d2,
                             c2Centre, c2xList, c2d2, elementsCountThroughWall):

    cache = fm.createFieldcache()
    # coordinates = findOrCreateFieldCoordinates(fm)

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
            make_tube_bifurcation_points_converging(paCentre, paxList[n3], pad2[n3], c1Centre, c1xList[n3], c1d2[n3],
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


def make_tube_bifurcation_points_converging(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):
    """
    Gets first ring of coordinates and derivatives between parent pa and
    children c1, c2, and over the crotch between c1 and c2.
    :return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex
    """
    paCount = len(pax)
    c1Count = len(c1x)
    c2Count = len(c2x)
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)
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


def createUterusMesh3D(region, fm, coordinates, geometricNetworkLayout, elementsCountAround, elementsCountThroughWall,
                       elementsCountInRightHorn, elementsCountInLeftHorn, elementsCountInCervix, elementsCountInVagina,
                       wallThickness, useCrossDerivatives):

    mesh = fm.findMeshByDimension(3)

    firstNodeIdentifier = 1
    firstElementIdentifier = 1

    cx_right_horn_group = geometricNetworkLayout.cxGroups[0]
    cx_left_horn_group = geometricNetworkLayout.cxGroups[1]
    cx_cervix_group = geometricNetworkLayout.cxGroups[2]
    cx_vagina_group = geometricNetworkLayout.cxGroups[3]

    sx_right_horn_group = geometricNetworkLayout.sxGroups[0]
    sx_left_horn_group = geometricNetworkLayout.sxGroups[1]
    sx_cervix_group = geometricNetworkLayout.sxGroups[2]
    # sx_vagina_group = geometricNetworkLayout.sxGroups[3]

    # Create annotation groups
    rightHornGroup = AnnotationGroup(region, get_uterus_term("right uterine horn"))
    leftHornGroup = AnnotationGroup(region, get_uterus_term("left uterine horn"))
    cervixGroup = AnnotationGroup(region, get_uterus_term("uterine cervix"))
    vaginaGroup = AnnotationGroup(region, get_uterus_term("vagina"))
    uterusGroup = AnnotationGroup(region, get_uterus_term("uterus"))
    annotationGroups = [cervixGroup, vaginaGroup, leftHornGroup, rightHornGroup, uterusGroup]

    rightHornMeshGroup = rightHornGroup.getMeshGroup(mesh)
    leftHornMeshGroup = leftHornGroup.getMeshGroup(mesh)
    cervixMeshGroup = cervixGroup.getMeshGroup(mesh)
    vaginaMeshGroup = vaginaGroup.getMeshGroup(mesh)
    uterusMeshGroup = uterusGroup.getMeshGroup(mesh)

    # Get right horn nodes
    rightHornCoordinates = getCoordinatesAlongTube3D(cx_right_horn_group, elementsCountAround,
                                                     elementsCountInRightHorn, elementsCountThroughWall,
                                                     wallThickness, startRadian=-math.pi / 2)

    rhLastRingNodeCoordinates = getTargetedRingNodesCoordinates(rightHornCoordinates, elementsCountAround,
                                                                elementsCountInRightHorn, elementsCountThroughWall,
                                                                omitStartRows=0, omitEndRows=1)

    rhLastRingNodeId, nodeCount = getTargetedRingNodesId(firstNodeIdentifier, elementsCountAround,
                                                         elementsCountInRightHorn, elementsCountThroughWall,
                                                         omitStartRows=0, omitEndRows=1)

    # Get left horn nodes
    leftHornCoordinates = getCoordinatesAlongTube3D(cx_left_horn_group, elementsCountAround,
                                                    elementsCountInLeftHorn, elementsCountThroughWall,
                                                    wallThickness, startRadian=-math.pi / 2)

    lhLastRingNodeCoordinates = getTargetedRingNodesCoordinates(leftHornCoordinates, elementsCountAround,
                                                                elementsCountInLeftHorn, elementsCountThroughWall,
                                                                omitStartRows=0, omitEndRows=1)

    lhLastRingNodeId, nodeCount = getTargetedRingNodesId(nodeCount, elementsCountAround, elementsCountInLeftHorn,
                                                         elementsCountThroughWall, omitStartRows=0, omitEndRows=1)

    # Get cervix nodes
    cervixCoordinates = getCoordinatesAlongTube3D(cx_cervix_group, elementsCountAround, elementsCountInCervix,
                                                  elementsCountThroughWall, wallThickness, startRadian=-math.pi / 2)

    cFirstRingNodeCoordinates = getTargetedRingNodesCoordinates(cervixCoordinates, elementsCountAround,
                                                                elementsCountInCervix, elementsCountThroughWall,
                                                                omitStartRows=1, omitEndRows=0)

    # Get vagina nodes
    vaginaCoordinates = getCoordinatesAlongTube3D(cx_vagina_group, elementsCountAround, elementsCountInVagina,
                                                  elementsCountThroughWall, wallThickness, startRadian=-math.pi / 2)

    # vFirstRingNodeCoordinates = getTargetedRingNodesCoordinates(vaginaCoordinates, elementsCountAround,
    #                                                             elementsCountInVagina, elementsCountThroughWall,
    #                                                             omitStartRows=0, omitEndRows=0)

    # Create nodes
    # Create right horn nodes
    nodeIdentifier = generateTubeNodes(fm, coordinates, firstNodeIdentifier, rightHornCoordinates,
                                       elementsCountInRightHorn, elementsCountAround, elementsCountThroughWall,
                                       omitStartRows=0, omitEndRows=1, startNodes=None)

    # Create left horn nodes
    nodeIdentifier = generateTubeNodes(fm, coordinates, nodeIdentifier, leftHornCoordinates, elementsCountInLeftHorn,
                                       elementsCountAround, elementsCountThroughWall, omitStartRows=0,
                                       omitEndRows=1, startNodes=None)

    # Create bifurcation nodes
    paCentre = sx_cervix_group[0][1]
    c1Centre = sx_right_horn_group[0][-2]
    c2Centre = sx_left_horn_group[0][-2]
    paxList = cFirstRingNodeCoordinates[0]
    # pad1List = cFirstRingNodeCoordinates[1]
    # pad2List = cFirstRingNodeCoordinates[2]
    pad2 = cFirstRingNodeCoordinates[2]
    c1xList = rhLastRingNodeCoordinates[0]
    c1d2 = rhLastRingNodeCoordinates[2]
    c2xList = lhLastRingNodeCoordinates[0]
    c2d2 = lhLastRingNodeCoordinates[2]
    nodeIdentifier, roNodeId, coNodeId = create3dBifurcationNodes(fm, coordinates, nodeIdentifier, paCentre, paxList,
                                                                  pad2, c1Centre, c1xList, c1d2, c2Centre, c2xList,
                                                                  c2d2, elementsCountThroughWall)

    # Create cervix nodes
    nodeCount = nodeIdentifier
    nodeIdentifier = generateTubeNodes(fm, coordinates, nodeIdentifier, cervixCoordinates, elementsCountInCervix,
                                       elementsCountAround, elementsCountThroughWall, omitStartRows=1,
                                       omitEndRows=0, startNodes=None)

    # Create vagina nodes
    nodeIdentifier = generateTubeNodes(fm, coordinates, nodeIdentifier, vaginaCoordinates, elementsCountInVagina,
                                       elementsCountAround, elementsCountThroughWall, omitStartRows=1,
                                       omitEndRows=0, startNodes=None)

    # Create elements
    # Create right horn elements
    startNodeId = firstNodeIdentifier
    elementIdentifier = \
        generateTubeElements(fm, coordinates, startNodeId, firstElementIdentifier, elementsCountInRightHorn,
                             elementsCountAround, elementsCountThroughWall, useCrossDerivatives, omitStartRows=0,
                             omitEndRows=1, meshGroups=[rightHornMeshGroup, uterusMeshGroup])

    # Create left horn elements
    startNodeId = rhLastRingNodeId[-1][-1] + 1
    elementIdentifier = generateTubeElements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInLeftHorn,
                                             elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                             omitStartRows=0, omitEndRows=1,
                                             meshGroups=[leftHornMeshGroup, uterusMeshGroup])

    # Create bifurcation elements
    cFirstRingNodeId, nodeCount = getTargetedRingNodesId(nodeCount, elementsCountAround, elementsCountInCervix,
                                                         elementsCountThroughWall, omitStartRows=1, omitEndRows=0)
    paNodeId = cFirstRingNodeId
    c1NodeId = rhLastRingNodeId
    c2NodeId = lhLastRingNodeId
    elementIdentifier = make_tube_bifurcation_elements_3d(fm, coordinates, elementIdentifier,
                                                          elementsCountAround, elementsCountThroughWall, paNodeId,
                                                          c1NodeId, c2NodeId, roNodeId, coNodeId,
                                                          meshGroups=[cervixMeshGroup, rightHornMeshGroup,
                                                                      leftHornMeshGroup, uterusMeshGroup])

    # Create cervix elements
    startNodeId = paNodeId[0][0]
    elementIdentifier = generateTubeElements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInCervix,
                                             elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                             omitStartRows=1, omitEndRows=0,
                                             meshGroups=[cervixMeshGroup, uterusMeshGroup])

    # Create vagina elements
    startNodeId = paNodeId[0][0] + (elementsCountInCervix - 1) * elementsCountAround * (elementsCountThroughWall + 1)
    elementIdentifier = generateTubeElements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInVagina,
                                             elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
                                             omitStartRows=0, omitEndRows=0,
                                             meshGroups=[vaginaMeshGroup])

    return nodeIdentifier, elementIdentifier, annotationGroups


def createUterusMesh3DRat(region, fm, coordinates, geometricNetworkLayout, elementsCountAround, elementsCountAcross,
                          elementsCountThroughWall, elementsCountInRightHorn, elementsCountInLeftHorn,
                          elementsCountInCervix, elementsCountInVagina, wallThickness, useCrossDerivatives):

    mesh = fm.findMeshByDimension(3)

    firstNodeIdentifier = 1
    firstElementIdentifier = 1

    cx_right_horn_group = geometricNetworkLayout.cxGroups[0]
    cx_left_horn_group = geometricNetworkLayout.cxGroups[1]
    cx_cervix_group = geometricNetworkLayout.cxGroups[2]
    cx_vagina_group = geometricNetworkLayout.cxGroups[3]

    sx_right_horn_group = geometricNetworkLayout.sxGroups[0]
    sx_left_horn_group = geometricNetworkLayout.sxGroups[1]
    sx_cervix_group = geometricNetworkLayout.sxGroups[2]
    # sx_vagina_group = geometricNetworkLayout.sxGroups[3]

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

    # elementsCountsAroundVessels, elementsCountAroundMid = \
    #     getOstiumElementsCountsAroundVessels(elementsCountAroundOstium, elementsCountAcross, vesselsCount)

    # elementsCountAroundMid = ((elementsCountAround + 1) // 6)
    # countInner = 2 * (elementsCountAroundMid + elementsCountAcross)
    # countOuter = (elementsCountAround - 2 * elementsCountAroundMid) // 2 + elementsCountAcross

    count = elementsCountAround // 2 + elementsCountAcross
    elementsCountAroundRightHorn = count
    elementsCountAroundLeftHorn = count
    # print('count', count)
    # print('countInner', countInner)
    # print('countOuter', countOuter)

    # Get right horn nodes
    rightHornCoordinates = getCoordinatesAlongTube3D(cx_right_horn_group, elementsCountAroundRightHorn,
                                                     elementsCountInRightHorn, elementsCountThroughWall,
                                                     wallThickness, startRadian=-math.pi / 2)

    rhLastRingNodeCoordinates = getTargetedRingNodesCoordinates(rightHornCoordinates, elementsCountAroundRightHorn,
                                                                elementsCountInRightHorn, elementsCountThroughWall,
                                                                omitStartRows=0, omitEndRows=1)

    rhLastRingNodeId, nodeCount = getTargetedRingNodesId(firstNodeIdentifier, elementsCountAroundRightHorn,
                                                         elementsCountInRightHorn, elementsCountThroughWall,
                                                         omitStartRows=0, omitEndRows=1)

    # Get left horn nodes
    leftHornCoordinates = getCoordinatesAlongTube3D(cx_left_horn_group, elementsCountAroundLeftHorn,
                                                    elementsCountInLeftHorn, elementsCountThroughWall,
                                                    wallThickness, startRadian=-math.pi / 2)

    lhLastRingNodeCoordinates = getTargetedRingNodesCoordinates(leftHornCoordinates, elementsCountAroundLeftHorn,
                                                                elementsCountInLeftHorn, elementsCountThroughWall,
                                                                omitStartRows=0, omitEndRows=1)

    lhLastRingNodeId, nodeCount = getTargetedRingNodesId(nodeCount, elementsCountAroundLeftHorn, elementsCountInLeftHorn,
                                                         elementsCountThroughWall, omitStartRows=0, omitEndRows=1)


    # Get cervix right and left path
    cx_cervix_group_right = cx_cervix_group[1:]
    cx_cervix_group_left = cx_cervix_group[1:]
    distance = wallThickness
    xrList = []
    xlList = []
    for n in range(len(cx_cervix_group[0])):
        x = cx_cervix_group[0][n]
        v = vector.normalise(cx_cervix_group[2][n])
        v_trans = vector.setMagnitude(v, distance)
        x_right = [x[c] + v_trans[c] for  c in range(3)]
        x_left = [x[c] - v_trans[c] for  c in range(3)]
        xrList.append(x_right)
        xlList.append(x_left)
    cx_cervix_group_right.insert(0, xrList)
    cx_cervix_group_left.insert(0, xlList)

    # Get right inner cervix nodes
    cervixInnerRightCoordinates = getCoordinatesAlongTube2D(cx_cervix_group_right, elementsCountAroundRightHorn,
                                                            elementsCountInCervix, startRadian=-math.pi / 2)

    # Get left inner cervix nodes
    cervixInnerLeftCoordinates = getCoordinatesAlongTube2D(cx_cervix_group_left, elementsCountAroundLeftHorn,
                                                           elementsCountInCervix, startRadian=-math.pi / 2)

    # cFirstRingNodeCoordinates = getTargetedRingNodesCoordinates(cervixCoordinates, elementsCountAround,
    #                                                             elementsCountInCervix, elementsCountThroughWall,
    #                                                             omitStartRows=1, omitEndRows=0)

    # Find outer nodes along cervix
    # sx_cervix = []
    # for n in range(len(sx_right_tube_group_cervix[0])):
    #     x = [(sx_right_tube_group_cervix[0][n][c] + sx_left_tube_group_cervix[0][n][c]) / 2 for c in range(3)]
    #     sx_cervix.append(x)
    cervix_radius1 = []
    cervix_radius2 = []
    sx_cervix = cx_cervix_group[0]
    for n in range(len(cx_cervix_group[0])):
        v2 = cx_cervix_group_right[0][n]
        v1 = cx_cervix_group_left[0][n]
        v1v2 = [v2[c] - v1[c] for c in range(3)]
        d1Dir = vector.normalise(v1v2)
        d1Mag = vector.magnitude(v1v2)
        sd2RMag = vector.magnitude(cx_cervix_group_right[2][n])
        sd3RMag = vector.magnitude(cx_cervix_group_right[4][n])
        # sd2LMag = sx_left_tube_group[2][n]
        radius1 = d1Mag / 2 + sd2RMag + wallThickness
        radius2 = sd3RMag + wallThickness
        radius1_v = vector.setMagnitude(vector.normalise(v1v2), radius1)
        radius2_v = vector.setMagnitude(vector.normalise(cx_cervix_group_right[4][n]), radius2)
        # radius2_v = vector.setMagnitude(vector.normalise(sx_right_tube_group_cervix[4][n]), radius1)
        cervix_radius1.append(radius1_v)
        cervix_radius2.append(radius2_v)
    sx_cervix_group = [sx_cervix, [], cervix_radius1, [], cervix_radius2]
    # sx_cervix_group = [sx_cervix, [], cervix_radius1, [], cervix_radius1]

    cervixLength = geometricNetworkLayout.arcLengthOfGroupsAlong[2]
    startRadian = -math.pi / 2
    xCervix, d1Cervix, d2Cervix, _ = \
        findNodesAlongTubes2D(sx_cervix_group, elementsCountAround, elementsCountInCervix,
                              elementsCountThroughWall, wallThickness, cervixLength, startRadian)
    # cervixCoordinatesOuter = [xCervix, d1Cervix, d2Cervix, d3Cervix]

    # # Find d3 for cervix outer nodes
    # d3Cervix = []
    # for n2 in range(0, elementsCountInCervix + 1):
    #     d3CervixRaw = []
    #     for n1 in range(elementsCountAround):
    #         if n1 == 0:
    #             d3CervixRaw.append([0.0, 0.0, 0.0])
    #         elif 0 < n1 <= elementsCountAround // 2:
    #             v1 = cervixInnerRightCoordinates[n2][n1]
    #             v2 = xCervix[n2][n1]
    #             v1v2 = findDerivativeBetweenPoints(v1, v2)
    #             d3CervixRaw.append(v1v2)
    #         else:
    #             v1 = cervixInnerLeftCoordinates[n2][n1]
    #             v2 = xCervix[n2][n1]
    #             v1v2 = findDerivativeBetweenPoints(v1, v2)
    #             d3CervixRaw.append(v1v2)
    #             # d3CervixRaw.append([1.0, 0.0, 0.0])
    #         d3Cervix.append(d3CervixRaw)
    #
    # cervixCoordinatesOuter = [xCervix, d1Cervix, d2Cervix, d3Cervix]


    # Sample/create a layer of nodes in bifurcation means between end of inner horns and begining of inner cervix canal
    # Right part
    firstRingRightInnerCervix = [cervixInnerRightCoordinates[0][elementsCountAroundRightHorn:2*elementsCountAroundRightHorn], cervixInnerRightCoordinates[1][elementsCountAroundRightHorn:2*elementsCountAroundRightHorn],
                                 cervixInnerRightCoordinates[2][elementsCountAroundRightHorn:2*elementsCountAroundRightHorn]]
    lastRingRightInnerHorn = [rhLastRingNodeCoordinates[0][0], rhLastRingNodeCoordinates[1][0],
                              rhLastRingNodeCoordinates[2][0]]

    xRaw = []
    d2Raw = []
    elementsCountAlongBifurcation = 2
    for n1 in range(elementsCountAroundRightHorn):
        xAlong = []
        d2Along = []
        xAlong.append(lastRingRightInnerHorn[0][n1])
        xAlong.append(firstRingRightInnerCervix[0][n1])
        d2Along.append(lastRingRightInnerHorn[2][n1])
        d2Along.append(firstRingRightInnerCervix[2][n1])
        xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along, elementsCountAlongBifurcation,
                                                                        arcLengthDerivatives=True)[0:2]
        xRaw.append(xSampledAlong)
        d2Raw.append(d2SampledAlong)

    # Rearrange x and d2
    xSampledBifurcationRight = []
    d1SampledBifurcationRight = []
    d2SampledBifurcationRight = []
    for n2 in range(2 + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAroundRightHorn):
            x = xRaw[n1][n2]
            d2 = d2Raw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)
            # Calculate d1
            v1 = xRaw[n1][n2]
            v2 = xRaw[n1 + 1 if n1 < elementsCountAroundRightHorn - 1 else 0][n2]
            d1 = findDerivativeBetweenPoints(v1, v2)
            d1Around.append(d1)
        d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        xSampledBifurcationRight.append(xAround)
        d1SampledBifurcationRight.append(d1Smoothed)
        d2SampledBifurcationRight.append(d2Around)

    # # should be fixed
    # d3SampledBifurcationRight = []
    # for n2 in range(2 + 1):
    #     d3Around = []
    #     for n1 in range(elementsCountAround):
    #         v1 = d1SampledBifurcationRight[n2][n1]
    #         v2 = d2SampledBifurcationRight[n2][n1]
    #         v3 = vector.crossproduct3(v1, v2)
    #         d3Around.append(v3)
    #     d3SampledBifurcationRight.append(d3Around)

    innerBifurcationRight = [xSampledBifurcationRight, d1SampledBifurcationRight, d2SampledBifurcationRight]

    # Left part
    firstRingLeftInnerCervix = [cervixInnerLeftCoordinates[0][elementsCountAroundLeftHorn:2 * elementsCountAroundLeftHorn],
                                 cervixInnerLeftCoordinates[1][elementsCountAroundLeftHorn:2 * elementsCountAroundLeftHorn],
                                 cervixInnerLeftCoordinates[2][elementsCountAroundLeftHorn:2 * elementsCountAroundLeftHorn]]
    lastRingLeftInnerHorn = [lhLastRingNodeCoordinates[0][0], lhLastRingNodeCoordinates[1][0],
                              lhLastRingNodeCoordinates[2][0]]

    xRaw = []
    d2Raw = []
    elementsCountAlongBifurcation = 2
    for n1 in range(elementsCountAroundLeftHorn):
        xAlong = []
        d2Along = []
        xAlong.append(lastRingLeftInnerHorn[0][n1])
        xAlong.append(firstRingLeftInnerCervix[0][n1])
        d2Along.append(lastRingLeftInnerHorn[2][n1])
        d2Along.append(firstRingLeftInnerCervix[2][n1])
        xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                        elementsCountAlongBifurcation,
                                                                        arcLengthDerivatives=True)[0:2]
        xRaw.append(xSampledAlong)
        d2Raw.append(d2SampledAlong)

    # Rearrange x and d2
    xSampledBifurcationLeft = []
    d1SampledBifurcationLeft = []
    d2SampledBifurcationLeft = []
    for n2 in range(2 + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAroundLeftHorn):
            x = xRaw[n1][n2]
            d2 = d2Raw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)
            # Calculate d1
            v1 = xRaw[n1][n2]
            v2 = xRaw[n1 + 1 if n1 < elementsCountAroundLeftHorn - 1 else 0][n2]
            d1 = findDerivativeBetweenPoints(v1, v2)
            d1Around.append(d1)
        d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        xSampledBifurcationLeft.append(xAround)
        d1SampledBifurcationLeft.append(d1Smoothed)
        d2SampledBifurcationLeft.append(d2Around)

    innerBifurcationLeft = [xSampledBifurcationLeft, d1SampledBifurcationLeft, d2SampledBifurcationLeft]

    # Get bifurcation outer nodes
    paCentre = sx_cervix_group[0][1]
    c1Centre = sx_right_horn_group[0][-2]
    c2Centre = sx_left_horn_group[0][-2]
    paxList = xCervix[1]
    pad2 = d2Cervix[1]
    # pad3 = d3Cervix[1]
    # paxList = cFirstRingNodeCoordinates[0]
    # pad2 = cFirstRingNodeCoordinates[2]
    c1xList = rhLastRingNodeCoordinates[0][1]
    c1d2 = rhLastRingNodeCoordinates[2][1]
    c2xList = lhLastRingNodeCoordinates[0][1]
    c2d2 = lhLastRingNodeCoordinates[2][1]
    rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex = \
        make_tube_bifurcation_points_converging_2d(paCentre, paxList, pad2, c1Centre, c1xList, c1d2, c2Centre, c2xList, c2d2)
    # rod3 = pad3

    # Get coordinates across cervix septum, between two inner canals
    # elementsCountAcross = len(cox) + 1
    xAcrossSeptum = []
    d1AcrossSeptum = []
    for n in range(elementsCountInCervix + 1):
        oa = 0
        ob = elementsCountAround // 2
        v1 = xCervix[n][ob]
        v2 = xCervix[n][oa]
        v3 = [v1[c] / 2 + v2[c] / 2 for c in range(3)]
        v1v2 = [v2[c] - v1[c] for c in range(3)]
        nx = [xCervix[n][ob], v3, xCervix[n][oa]]
        nd1 = [[d / elementsCountAcross for d in v1v2], [d / elementsCountAcross for d in v1v2],
               [d / elementsCountAcross for d in v1v2]]
        px, pd1, pe, pxi = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAcross)[0:4]
        xAcrossSeptum.append(px)
        d1AcrossSeptum.append(pd1)

    # Find d2 across cervix septum
    d2Raw = []
    for n2 in range(elementsCountAcross + 1):
        xAlongSeptum = []
        d2AlongSeptum = []
        for n1 in range(elementsCountInCervix):
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
    for n2 in range(elementsCountInCervix + 1):
        d2Across = []
        for n1 in range(elementsCountAcross + 1):
            d2 = d2Raw[n1][n2]
            d2Across.append(d2)
        d2AcrossSeptum.append(d2Across)

    septumCervixCoordinates = [xAcrossSeptum, d1AcrossSeptum, d2AcrossSeptum]


    # # Add nodes between the two nodes of bifurcation in rox (along septum)
    # # Get coordinates across bifurcation septum, between two inner canals
    # # elementsCountAcross = len(cox) + 1
    # v1 = rox[elementsCountAround // 2]
    # v2 = rox[0]
    # v3 = [v1[c] / 2 + v2[c] / 2 for c in range(3)]
    # v1v2 = [v2[c] - v1[c] for c in range(3)]
    # nx = [v1, v3, v2]
    # nd1 = [[d / elementsCountAcross for d in v1v2], [d / elementsCountAcross for d in v1v2],
    #        [d / elementsCountAcross for d in v1v2]]
    # xSeptumBifurcation, d1SeptumBifurcation, pe, pxi = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAcross)[
    #                                                    0:4]
    # xSeptumBifurcation = xSeptumBifurcation[1:-1]
    # d1SeptumBifurcation = d1SeptumBifurcation[1:-1]
    #
    # # Find d2 across bifurcation septum
    # d2SeptumBifurcation = []
    # for n2 in range(elementsCountAcross - 1):
    #     v1 = xSeptumBifurcation[n2]
    #     v2 = xAcrossSeptum[1][n2]
    #     d2 = findDerivativeBetweenPoints(v1, v2)
    #     d2SeptumBifurcation.append(d2)
    #
    # # Find d3 across bifurcation septum
    # d3SeptumBifurcation = []
    # for n2 in range(elementsCountAcross - 1):
    #     v1 = d1SeptumBifurcation[n2]
    #     v2 = d2AcrossSeptum[1][n2]
    #     v3 = vector.crossproduct3(v1, v2)
    #     v3Norm = vector.normalise(v3)
    #     d1mag = vector.magnitude(v1)
    #     d3 = vector.setMagnitude(v3Norm, d1mag)
    #     d3SeptumBifurcation.append(d3)
    # # print('len(d3SeptumBifurcation)', len(d3SeptumBifurcation))
    #
    # septumBifurCoordinates = [xSeptumBifurcation, d1SeptumBifurcation, d2SeptumBifurcation, d3SeptumBifurcation]


    # Get d3 for inner right bifurcation tube
    d3SampledBifurcationRight = []
    for n1 in range(elementsCountAroundRightHorn):
        v1 = xSampledBifurcationRight[1][n1]
        if n1 <= elementsCountAround // 2:
            v2 = rox[n1]
        else:
            # v2 = xSeptumBifurcation[n1 - elementsCountAround // 2 - 1]
            v2 = cox[n1 - elementsCountAround // 2 - 1]
        d3 = findDerivativeBetweenPoints(v1, v2)
        d3SampledBifurcationRight.append(d3)

    # Get d3 for inner left bifurcation tube
    d3SampledBifurcationLeft = []
    for n1 in range(elementsCountAroundRightHorn):
        v1 = xSampledBifurcationLeft[1][n1]
        if n1 == 0:
            v2 = rox[n1]
        elif 0 < n1 < elementsCountAcross:
            v2 = cox[elementsCountAcross - n1 - 1]
        elif n1 == elementsCountAcross:
            v2 = rox[elementsCountAround // 2]
        else:
            v2 = rox[n1]
        d3 = findDerivativeBetweenPoints(v1, v2)
        d3SampledBifurcationLeft.append(d3)

    # Find d3 for cervix right inner canal nodes
    d3CervixInnerRight = []
    for n2 in range(0, elementsCountInCervix + 1):
        for n1 in range(elementsCountAroundRightHorn):
            v1 = cervixInnerRightCoordinates[0][n2 * elementsCountAroundRightHorn + n1]
            if n1 <= elementsCountAround // 2:
                v2 = xCervix[n2][n1]
            else:
                v2 = xAcrossSeptum[n2][n1 - elementsCountAround // 2]
            d3 = findDerivativeBetweenPoints(v1, v2)
            d3CervixInnerRight.append(d3)
    cervixInnerRightCoordinates.append(d3CervixInnerRight)

    # Find d3 for cervix Left inner canal nodes
    d3CervixInnerLeft = []
    for n2 in range(0, elementsCountInCervix + 1):
        for n1 in range(elementsCountAroundLeftHorn):
            v1 = cervixInnerLeftCoordinates[0][n2 * elementsCountAroundLeftHorn + n1]
            if n1 == 0:
                v2 = xCervix[n2][n1]
            elif 0 < n1 < elementsCountAcross:
                v2 = xAcrossSeptum[n2][elementsCountAcross - n1]
            elif n1 == elementsCountAcross:
                v2 = xCervix[n2][elementsCountAround // 2]
            else:
                v2 = xCervix[n2][n1]
            d3 = findDerivativeBetweenPoints(v1, v2)
            d3CervixInnerLeft.append(d3)
    cervixInnerLeftCoordinates.append(d3CervixInnerLeft)

    # Find d3 for cervix outer nodes
    d3Cervix = []
    for n2 in range(0, elementsCountInCervix + 1):
        d3CervixRaw = []
        for n1 in range(elementsCountAround):
            if n1 == 0:
                d1 = septumCervixCoordinates[1][n2][n1]
                d3CervixRaw.append(d1)
            elif 0 < n1 < elementsCountAround // 2:
                v1 = cervixInnerRightCoordinates[n2][n1]
                v2 = xCervix[n2][n1]
                v1v2 = findDerivativeBetweenPoints(v1, v2)
                d3CervixRaw.append(v1v2)
            elif n1 == elementsCountAround // 2:
                d = [-d1[c] for c in range(3)]
                d3CervixRaw.append(d)
            else:
                v1 = cervixInnerLeftCoordinates[n2][n1]
                v2 = xCervix[n2][n1]
                v1v2 = findDerivativeBetweenPoints(v1, v2)
                d3CervixRaw.append(v1v2)
                # d3CervixRaw.append([1.0, 0.0, 0.0])
            d3Cervix.append(d3CervixRaw)

    cervixCoordinatesOuter = [xCervix, d1Cervix, d2Cervix, d3Cervix]

    # Get d3 for outer bifurcation (for rox nodes)
    pad3 = d3Cervix[1] # parent d3 which is d3 for first row of cervix
    rod3 = []
    for n in range(elementsCountAround):
        if n == 0:
            rod3.append(pad3[n])
        elif 0 < n < elementsCountAround // 2:
            rod3.append(d3SampledBifurcationRight[n])
        elif n == elementsCountAround // 2:
            rod3.append(pad3[n])
        else:
            # should be d3 of the d3SampledBifurcationLeft
            rod3.append(d3SampledBifurcationLeft[n])
            # rod3.append([0.0, 0.0, 1.0])

    # Get d3 for outer bifurcation (for cox nodes)
    cod3 = []
    for n in range(len(cox)):
        # v1 = xSeptumBifurcation[n]
        v1 = septumCervixCoordinates[0][0][n + 1]
        v2 = cox[n]
        d3 = findDerivativeBetweenPoints(v1, v2)
        cod3.append(d3)

    # Create nodes
    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)


    # Create right horn nodes
    nodeIdentifier = generateTubeNodes(fm, coordinates, firstNodeIdentifier, rightHornCoordinates,
                                       elementsCountInRightHorn, elementsCountAroundRightHorn, elementsCountThroughWall,
                                       omitStartRows=0, omitEndRows=1, startNodes=None)

    # Create left horn nodes
    nodeIdentifier = generateTubeNodes(fm, coordinates, nodeIdentifier, leftHornCoordinates, elementsCountInLeftHorn,
                                       elementsCountAroundLeftHorn, elementsCountThroughWall, omitStartRows=0,
                                       omitEndRows=1, startNodes=None)

    # Create extra right inner nodes in bifurcation
    birNodeId = []  # bifurcation right inner ring node id
    for n2 in range(2):
        for n1 in range(elementsCountAroundRightHorn):
            if n2 == 1:
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledBifurcationRight[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledBifurcationRight[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledBifurcationRight[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3SampledBifurcationRight[n1])
                birNodeId.append(nodeIdentifier)
                nodeIdentifier += 1

    # Create extra left inner nodes in bifurcation
    bilNodeId = []  # bifurcation left inner ring node id
    for n2 in range(2):
        for n1 in range(elementsCountAroundLeftHorn):
            if n2 == 1:
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledBifurcationLeft[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledBifurcationLeft[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledBifurcationLeft[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3SampledBifurcationLeft[n1])
                bilNodeId.append(nodeIdentifier)
                nodeIdentifier += 1

    # Create bifurcation nodes
    septumBifurCoordinates = []
    nodeIdentifier, rox, cox, roNodeId, coNodeId, sbNodeId, nextNodeId = \
        create2DBifurcationNodes_mod(fm, nodeIdentifier, rox, rod1, rod2, rod3, cox, cod1, cod2, cod3, septumBifurCoordinates)
    #old
    # paCentre = sx_cervix_group[0][1]
    # c1Centre = sx_right_horn_group[0][-2]
    # c2Centre = sx_left_horn_group[0][-2]
    # paxList = xCervix[1]
    # pad2 = d2Cervix[1]
    # pad3 = d3Cervix[1]
    # # paxList = cFirstRingNodeCoordinates[0]
    # # pad2 = cFirstRingNodeCoordinates[2]
    # c1xList = rhLastRingNodeCoordinates[0][1]
    # c1d2 = rhLastRingNodeCoordinates[2][1]
    # c2xList = lhLastRingNodeCoordinates[0][1]
    # c2d2 = lhLastRingNodeCoordinates[2][1]
    # nodeIdentifier, rox, cox, roNodeId, coNodeId, nextNodeId, paStartIndex, c1StartIndex, c2StartIndex = \
    #     create2DBifurcationNodes(fm, nodeIdentifier, paCentre, paxList, pad2, pad3, c1Centre, c1xList, c1d2, c2Centre,
    #                              c2xList, c2d2)

    # Create cervix nodes
    nodeIdentifier, cricNodeId, clicNodeId, cotNodeId, csNodeId = generateCervixNodes(fm, nodeIdentifier, cervixInnerRightCoordinates, cervixInnerLeftCoordinates,
                                         cervixCoordinatesOuter, septumCervixCoordinates, elementsCountInCervix, elementsCountAround,
                                         elementsCountAcross, omitStartRows=1, omitEndRows=0, startNodes=None)


    # Create elements
    # Create right horn elements
    startNodeId = firstNodeIdentifier
    elementIdentifier = \
        generateTubeElements(fm, coordinates, startNodeId, firstElementIdentifier, elementsCountInRightHorn,
                             elementsCountAroundRightHorn, elementsCountThroughWall, useCrossDerivatives, omitStartRows=0,
                             omitEndRows=1, meshGroups=[rightHornMeshGroup, uterusMeshGroup])

    # Create left horn elements
    startNodeId = rhLastRingNodeId[-1][-1] + 1
    elementIdentifier = generateTubeElements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInLeftHorn,
                                             elementsCountAroundLeftHorn, elementsCountThroughWall, useCrossDerivatives,
                                             omitStartRows=0, omitEndRows=1,
                                             meshGroups=[leftHornMeshGroup, uterusMeshGroup])


    # Create bifurcation elements
    cFirstRingNodeId, nodeCount = getTargetedRingNodesId(nodeCount, elementsCountAround, elementsCountInCervix,
                                                         elementsCountThroughWall, omitStartRows=1, omitEndRows=0)
    paNodeId = cFirstRingNodeId
    c1NodeId = rhLastRingNodeId
    c2NodeId = lhLastRingNodeId
    # elementIdentifier = make_rat_uterus_bifurcation_elements(fm, coordinates, elementIdentifier,
    #                                                       elementsCountAround, elementsCountThroughWall, paNodeId,
    #                                                       c1NodeId, c2NodeId, roNodeId, coNodeId, birNodeId, bilNodeId,
    #                                                       cricNodeId, clicNodeId, cotNodeId, sbNodeId, csNodeId,
    #                                                        meshGroups=[cervixMeshGroup, rightHornMeshGroup, leftHornMeshGroup, uterusMeshGroup])
    elementIdentifier = make_rat_uterus_bifurcation_elements_modified(fm, coordinates, elementIdentifier,
                                                          elementsCountAround, elementsCountAcross, elementsCountThroughWall, paNodeId,
                                                          c1NodeId, c2NodeId, roNodeId, coNodeId, birNodeId, bilNodeId,
                                                          cricNodeId, clicNodeId, cotNodeId, sbNodeId, csNodeId,
                                                           meshGroups=[bodyMeshGroup, rightHornMeshGroup, leftHornMeshGroup, uterusMeshGroup])

    # Create cervix elements
    elementIdentifier = make_cervix_elements(mesh, coordinates, elementIdentifier, elementsCountInCervix,
                                             elementsCountAround, elementsCountAcross, elementsCountAroundRightHorn,
                                             elementsCountAroundLeftHorn, cricNodeId, clicNodeId, cotNodeId, csNodeId, useCrossDerivatives,
                                             meshGroups=[bodyMeshGroup, cervixMeshGroup, uterusMeshGroup])

    # # # Create vagina elements
    # # startNodeId = paNodeId[0][0] + (elementsCountInCervix - 1) * elementsCountAround * (elementsCountThroughWall + 1)
    # # elementIdentifier = generateTubeElements(fm, coordinates, startNodeId, elementIdentifier, elementsCountInVagina,
    # #                                          elementsCountAround, elementsCountThroughWall, useCrossDerivatives,
    # #                                          omitStartRows=0, omitEndRows=0,
    # #                                          meshGroups=[vaginaMeshGroup])
    # elementIdentifier = firstElementIdentifier

    return nodeIdentifier, elementIdentifier, annotationGroups

def getCoordinatesAlongTube2D(cx_group, elementsCountAround, elementsCountAlongTube, startRadian):

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

    xList = []
    d1List = []
    d2List = []
    for n2 in range(elementsCountAlongTube + 1):
        for n1 in range(elementsCountAround):
            xList.append(xSampledTube[n2][n1])
            d1List.append(d1SampledTube[n2][n1])
            d2List.append(d2SampledTube[n2][n1])

    coordinatesList = [xList, d1List, d2List]

    return coordinatesList

def generateTubeNodes2D(fm, nodeIdentifier, tubeCoordinates, elementsCountAlongTube, elementsCountAround,
                        omitStartRows, omitEndRows, startNodes=None):

    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)

    # Create tube nodes
    lastRingsNodeId = []
    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    firstRingsNodeId = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        lastRingNodeIdThroughWall = []
        xLastRingThroughWall = []
        d1LastRingThroughWall = []
        d2LastRingThroughWall = []
        firstRingNodeIdThroughWall = []
        xFirstRingThroughWall = []
        d1FirstRingThroughWall = []
        d2FirstRingThroughWall = []
        for n1 in range(elementsCountAround):
            n = n2 * elementsCountAround + n1
            x = tubeCoordinates[0][n]
            d1 = tubeCoordinates[1][n]
            d2 = tubeCoordinates[2][n]
            if omitEndRows == 1:  # merging to the bifurcation
                if n2 == elementsCountAlongTube:
                    pass
                else:
                    node = nodes.createNode(nodeIdentifier, nodetemplate)
                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                    if n2 == elementsCountAlongTube - 1:
                        lastRingNodeIdThroughWall.append(nodeIdentifier)
                        xLastRingThroughWall.append(x)
                        d1LastRingThroughWall.append(d1)
                        d2LastRingThroughWall.append(d2)
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
                    if n2 == 1:
                        firstRingNodeIdThroughWall.append(nodeIdentifier)
                        xFirstRingThroughWall.append(x)
                        d1FirstRingThroughWall.append(d1)
                        d2FirstRingThroughWall.append(d2)
                    nodeIdentifier += 1
            else:
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
                nodeIdentifier += 1
        if omitEndRows == 1:
            if n2 == elementsCountAlongTube - 1:
                lastRingsNodeId.append(lastRingNodeIdThroughWall)
                xLastRing.append(xLastRingThroughWall)
                d1LastRing.append(d1LastRingThroughWall)
                d2LastRing.append(d2LastRingThroughWall)
        elif omitStartRows == 1:
            if n2 == 1:
                firstRingsNodeId.append(firstRingNodeIdThroughWall)
                xFirstRing.append(xFirstRingThroughWall)
                d1FirstRing.append(d1FirstRingThroughWall)
                d2FirstRing.append(d2FirstRingThroughWall)

    return nodeIdentifier

def findNodesAlongTubes2D(sx_group, elementsCountAround, elementsCountAlongTube, elementsCountThroughWall,
                          wallThickness, tubeLength, startRadian):

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

    d3Tube = []
    for n2 in range(elementsCountAlongTube + 1):
        d3Around = []
        for n1 in range(elementsCountAround):
            d3Around.append(vector.normalise(
                vector.crossproduct3(vector.normalise(d1SampledTube[n2][n1]), vector.normalise(d2SampledTube[n2][n1]))))
        d3Tube.append(d3Around)

    return xSampledTube, d1SampledTube, d2SampledTube, d3Tube


def generateCervixNodes(fm, nodeIdentifier, xInnerRigh, xInnerLeft, xOuter, xAcross, elementsCountInCervix,
                        elementsCountAround, elementsCountAcross, omitStartRows, omitEndRows, startNodes=None):

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

    cricNodeId = []  # cervix right inner canal node ids
    clicNodeId = []  # cervix left inner canal node ids
    cotNodeId = []  # cervix outer tube node id
    csNodeId = []   # cervix septum node Id
    for n2 in range(1 if omitStartRows == 1 else 0, elementsCountInCervix + 1):
        for n1 in range(elementsCountAroundRightHorn):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            n = n2 * elementsCountAroundRightHorn + n1
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xInnerRigh[0][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xInnerRigh[1][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xInnerRigh[2][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, xInnerRigh[3][n])
            # if n2 == 1:
            cricNodeId.append(nodeIdentifier)
            nodeIdentifier += 1
        for n1 in range(elementsCountAroundLeftHorn):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            n = n2 * elementsCountAroundLeftHorn + n1
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xInnerLeft[0][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xInnerLeft[1][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xInnerLeft[2][n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, xInnerLeft[3][n])
            clicNodeId.append(nodeIdentifier)
            nodeIdentifier += 1
        for n1 in range(1, elementsCountAcross):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAcross[0][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xAcross[1][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xAcross[2][n2][n1])
            csNodeId.append(nodeIdentifier)
            nodeIdentifier += 1
        for n1 in range(elementsCountAround):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOuter[0][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, xOuter[1][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, xOuter[2][n2][n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, xOuter[3][n2][n1])
            cotNodeId.append(nodeIdentifier)
            nodeIdentifier += 1

    # # Create tube nodes
    # lastRingsNodeId = []
    # xLastRing = []
    # d1LastRing = []
    # d2LastRing = []
    # firstRingsNodeId = []
    # xFirstRing = []
    # d1FirstRing = []
    # d2FirstRing = []
    # for n2 in range(elementsCountAlongTube + 1):
    #     lastRingNodeIdThroughWall = []
    #     xLastRingThroughWall = []
    #     d1LastRingThroughWall = []
    #     d2LastRingThroughWall = []
    #     firstRingNodeIdThroughWall = []
    #     xFirstRingThroughWall = []
    #     d1FirstRingThroughWall = []
    #     d2FirstRingThroughWall = []
    #     for n1 in range(elementsCountAround):
    #         n = n2 * elementsCountAround + n1
    #         x = tubeCoordinates[0][n]
    #         d1 = tubeCoordinates[1][n]
    #         d2 = tubeCoordinates[2][n]
    #         if omitEndRows == 1:  # merging to the bifurcation
    #             if n2 == elementsCountAlongTube:
    #                 pass
    #             else:
    #                 node = nodes.createNode(nodeIdentifier, nodetemplate)
    #                 cache.setNode(node)
    #                 coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    #                 coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
    #                 coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
    #                 if n2 == elementsCountAlongTube - 1:
    #                     lastRingNodeIdThroughWall.append(nodeIdentifier)
    #                     xLastRingThroughWall.append(x)
    #                     d1LastRingThroughWall.append(d1)
    #                     d2LastRingThroughWall.append(d2)
    #                 nodeIdentifier += 1
    #         elif omitStartRows == 1:  # diverging from bifurcation
    #             if n2 == 0:
    #                 pass
    #             else:
    #                 node = nodes.createNode(nodeIdentifier, nodetemplate)
    #                 cache.setNode(node)
    #                 coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    #                 coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
    #                 coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
    #                 if n2 == 1:
    #                     firstRingNodeIdThroughWall.append(nodeIdentifier)
    #                     xFirstRingThroughWall.append(x)
    #                     d1FirstRingThroughWall.append(d1)
    #                     d2FirstRingThroughWall.append(d2)
    #                 nodeIdentifier += 1
    #         else:
    #             node = nodes.createNode(nodeIdentifier, nodetemplate)
    #             cache.setNode(node)
    #             coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    #             coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1)
    #             coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
    #             nodeIdentifier += 1
    #     if omitEndRows == 1:
    #         if n2 == elementsCountAlongTube - 1:
    #             lastRingsNodeId.append(lastRingNodeIdThroughWall)
    #             xLastRing.append(xLastRingThroughWall)
    #             d1LastRing.append(d1LastRingThroughWall)
    #             d2LastRing.append(d2LastRingThroughWall)
    #     elif omitStartRows == 1:
    #         if n2 == 1:
    #             firstRingsNodeId.append(firstRingNodeIdThroughWall)
    #             xFirstRing.append(xFirstRingThroughWall)
    #             d1FirstRing.append(d1FirstRingThroughWall)
    #             d2FirstRing.append(d2FirstRingThroughWall)

    return nodeIdentifier, cricNodeId, clicNodeId, cotNodeId, csNodeId


# def create2DBifurcationNodes(fm, nodeIdentifier, paCentre, pax, pad2, pad3, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):
#
#     cache = fm.createFieldcache()
#     coordinates = findOrCreateFieldCoordinates(fm)
#
#     nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
#     nodetemplate = nodes.createNodetemplate()
#     nodetemplate.defineField(coordinates)
#     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
#     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
#     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
#     nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
#
#     rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex = \
#         make_tube_bifurcation_points_converging_2d(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2)
#
#     rod3 = pad3
#     # Create bifurcation nodes
#     roNodeId = []
#     for n in range(len(rox)):
#         node = nodes.createNode(nodeIdentifier, nodetemplate)
#         cache.setNode(node)
#         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rox[n])
#         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1[n])
#         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2[n])
#         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, rod3[n])
#         roNodeId.append(nodeIdentifier)
#         nodeIdentifier = nodeIdentifier + 1
#     coNodeId = []
#     for n in range(len(cox)):
#         node = nodes.createNode(nodeIdentifier, nodetemplate)
#         cache.setNode(node)
#         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cox[n])
#         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cod1[n])
#         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cod2[n])
#         coNodeId.append(nodeIdentifier)
#         nodeIdentifier = nodeIdentifier + 1
#
#     nextNodeId = nodeIdentifier
#     return nodeIdentifier, rox, cox, roNodeId, coNodeId, nextNodeId, paStartIndex, c1StartIndex, c2StartIndex

def create2DBifurcationNodes_mod(fm, nodeIdentifier, rox, rod1, rod2, rod3, cox, cod1, cod2, cod3, septumBifurCoordinates):

    cache = fm.createFieldcache()
    coordinates = findOrCreateFieldCoordinates(fm)

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplate = nodes.createNodetemplate()
    nodetemplate.defineField(coordinates)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

    # rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex = \
    #     make_tube_bifurcation_points_converging_2d(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2)
    #
    # rod3 = pad3
    # Create bifurcation nodes
    roNodeId = []
    for n in range(len(rox)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, rox[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, rod1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, rod2[n])
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

    # Create septum nodes in bifurcation
    sbNodeId = []
    # for n in range(len(septumBifurCoordinates[0])):
    #     node = nodes.createNode(nodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, septumBifurCoordinates[0][n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, septumBifurCoordinates[1][n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, septumBifurCoordinates[2][n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, septumBifurCoordinates[3][n])
    #     sbNodeId.append(nodeIdentifier)
    #     nodeIdentifier += 1

    nextNodeId = nodeIdentifier
    return nodeIdentifier, rox, cox, roNodeId, coNodeId, sbNodeId, nextNodeId

def getTargetedRingNodesCoordinates2D(tubeCoordinates, elementsCountAround, elementsCountAlongTube, omitStartRows,
                                      omitEndRows):

    xLastRing = []
    d1LastRing = []
    d2LastRing = []
    xFirstRing = []
    d1FirstRing = []
    d2FirstRing = []
    for n2 in range(elementsCountAlongTube + 1):
        for n1 in range(elementsCountAround):
            n = n2 * elementsCountAround + n1
            x = tubeCoordinates[0][n]
            d1 = tubeCoordinates[1][n]
            d2 = tubeCoordinates[2][n]
            if omitEndRows == 1:  # merging to the bifurcation
                if n2 == elementsCountAlongTube:
                    pass
                else:
                    if n2 == elementsCountAlongTube - 1:
                        xLastRing.append(x)
                        d1LastRing.append(d1)
                        d2LastRing.append(d2)
            elif omitStartRows == 1:  # diverging from bifurcation
                if n2 == 0:
                    pass
                else:
                    if n2 == 1:
                        xFirstRing.append(x)
                        d1FirstRing.append(d1)
                        d2FirstRing.append(d2)

    if omitStartRows == 1:
        targetedRingCoordinates = [xFirstRing, d1FirstRing, d2FirstRing]
    elif omitEndRows == 1:
        targetedRingCoordinates = [xLastRing, d1LastRing, d2LastRing]
    else:
        targetedRingCoordinates = []

    return targetedRingCoordinates


def make_tube_bifurcation_points_converging_2d(paCentre, pax, pad2, c1Centre, c1x, c1d2, c2Centre, c2x, c2d2):
    """
    Gets first ring of coordinates and derivatives between parent pa and
    children c1, c2, and over the crotch between c1 and c2.
    :return rox, rod1, rod2, cox, cod1, cod2, paStartIndex, c1StartIndex, c2StartIndex
    """
    paCount = len(pax)
    c1Count = len(c1x)
    c2Count = len(c2x)
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)
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


def make_rat_uterus_bifurcation_elements(fm, coordinates, elementIdentifier, elementsCountAround,
                                      elementsCountThroughWall, paNodeId, c1NodeId, c2NodeId, roNodeId, coNodeId,
                                      birNodeId, bilNodeId, cricNodeId, clicNodeId, cotNodeId, sbNodeId, csNodeId, meshGroups=None):

    paCount = len(paNodeId[0])
    c1Count = len(c1NodeId[0])
    c2Count = len(c2NodeId[0])
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)

    elementsCountAcross = len((sbNodeId)) + 1

    # fm = region.getFieldmodule()
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
        for e1 in range(elementsCountAround):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 < elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = birNodeId[e1]
                bni4 = birNodeId[(e1 + 1) % c1Count]
                bni5 = bni1 + elementsCountAround
                bni6 = bni5 + 1
                bni7 = bni3 + 2 * elementsCountAround
                bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            elif e1 == elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = birNodeId[e1]
                bni4 = birNodeId[(e1 + 1) % c1Count]
                bni5 = bni1 + elementsCountAround
                bni6 = bni5 + 1
                bni7 = bni3 + 2 * elementsCountAround
                bni8 = coNodeId[0]
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            else:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = birNodeId[e1]
                bni5 = bni1 + elementsCountAround
                bni7 = coNodeId[e1 - elementsCountAround // 2 - 1]
                if e1 == c1Count - 1:
                    bni4 = birNodeId[0]
                    bni6 = bni5 - elementsCountAround + 1
                    bni8 = roNodeId[0]
                else:
                    bni4 = bni3 + 1
                    bni6 = bni5 + 1
                    bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
            if e1 in (0, pac1Count - 1, pac1Count, c1Count - 1):
                eft1 = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                elif e1 == pac1Count - 1:
                    # eft = eftfactory.createEftBasic()
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == pac1Count:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif e1 == c1Count - 1:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
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

    # Fundus right tube part
    elementsCountAcross = len((sbNodeId)) + 1
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAcross):
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            if e1 == 0:
                bni1 = coNodeId[e1]
                bni2 = birNodeId[e1 + elementsCountAround // 2]
                bni3 = bni2 + 1
                bni4 = roNodeId[elementsCountAround // 2]
                bni5 = sbNodeId[e1]
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5]
                # print('nodeIdentifiers', nodeIdentifiers)
            elif 0 < e1 < elementsCountAcross - 1:
                bni1 = coNodeId[e1 - 1]
                bni2 = bni1 + 1
                bni3 = birNodeId[e1 + elementsCountAround // 2]
                bni4 = bni3 + 1
                bni5 = sbNodeId[e1 - 1]
                bni6 = bni5 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6]
                # print('nodeIdentifiers', nodeIdentifiers)
            elif e1 == elementsCountAcross - 1:
                bni1 = coNodeId[-1]
                bni2 = birNodeId[e1 + elementsCountAround // 2]
                bni3 = birNodeId[0]
                bni4 = sbNodeId[-1]
                bni5 = roNodeId[0]
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5]
                # print('nodeIdentifiers', nodeIdentifiers)

            # Remap derivatives for elements at the beginning and end of the septum
            if e1 == 0:
                va = e1
                vb = e1 + 1
                # set general linear map coefficients
                radiansPerElementAround = math.pi / elementsCountAcross
                radiansAround = va * radiansPerElementAround
                radiansAroundNext = vb * radiansPerElementAround
                scalefactors = [-1.0,
                                math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                                math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround]
                eft1 = eftfactory.createEftPyramidBottomSimple(va * 10000, vb * 10000)
                # remapEftNodeValueLabel(eft1, [1], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1]), (Node.VALUE_LABEL_D_DS3, [])])
                # remapEftNodeValueLabel(eft1, [4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                # remapEftNodeValueLabel(eft1, [2], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [1]), (Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                elementtemplateMod.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateMod
            elif 0 < e1 < elementsCountAcross - 1:
                eft1 = eftfactory.createEftWedgeCollapseXi2Quadrant([5, 6])
                # remapEftNodeValueLabel(eft1, [1, 2 ], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                elementtemplateMod.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateMod

                # eft1 = eftfactory.createEftBasic()
                # nodes = [1, 2, 3, 4]
                # collapseNodes = [1, 2]
                # remapEftNodeValueLabel(eft1, nodes, Node.VALUE_LABEL_D_DS2, [])
                # # remapEftNodeValueLabel(eft1, collapseNodes, Node.VALUE_LABEL_D_DS3,
                # #                        [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS3, [])])
                # remapEftNodeValueLabel(eft1, [1, 2, 5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                # ln_map = [1, 2, 1, 2, 3, 4, 5, 6]
                # remapEftLocalNodes(eft1, 6, ln_map)
                # elementtemplateMod.defineField(coordinates, -1, eft1)
                # elementtemplate1 = elementtemplateMod
            # elif e1 == elementsCountAcross - 1:
            #     va = e1
            #     vb = e1 + 1
            #     # set general linear map coefficients
            #     radiansPerElementAround = math.pi / elementsCountAcross
            #     radiansAround = va * radiansPerElementAround
            #     radiansAroundNext = vb * radiansPerElementAround
            #     scalefactors = [-1.0,
            #                     math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
            #                     math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround]
            #     eft1 = eftfactory.createEftPyramidBottomSimple(va * 10000, vb * 10000)
            #     elementtemplateMod.defineField(coordinates, -1, eft1)
            #     elementtemplate1 = elementtemplateMod
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


    # # Left tube part
    # for e3 in range(elementsCountThroughWall):
    #     for e1 in range(elementsCountAround):
    #         eft = eftStd
    #         elementtemplate = elementtemplateStd
    #         scalefactors = None
    #         bni1 = c2NodeId[e3][e1]
    #         bni2 = c2NodeId[e3][(e1 + 1) % elementsCountAround]
    #         bni3 = bilNodeId[e1]
    #         bni4 = bilNodeId[(e1 + 1) % elementsCountAround]
    #         bni5 = bni1 + elementsCountAround
    #         bni6 = bni2 + elementsCountAround
    #         if 0 <= e1 < elementsCountAcross:
    #             if e1 == 0:
    #                 # bni7 = roNodeId[0]
    #                 # bni8 = sbNodeId[-1]
    #                 bni7 = roNodeId[0]
    #                 bni8 = coNodeId[-1]
    #             elif e1 == elementsCountAcross - 1:
    #                 # bni7 = sbNodeId[elementsCountAcross - 1 - e1]
    #                 # bni8 = roNodeId[elementsCountAround // 2]
    #                 bni7 = coNodeId[elementsCountAcross - 1 - e1]
    #                 bni8 = roNodeId[elementsCountAround // 2]
    #             else:
    #                 # bni7 = sbNodeId[elementsCountAcross - 1 - e1]
    #                 # bni8 = bni7 - 1
    #                 bni7 = coNodeId[elementsCountAcross - 1 - e1]
    #                 bni8 = bni7 - 1
    #             nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
    #             # print('nodeIdentifiers', nodeIdentifiers)
    #         else:
    #             bni7 = roNodeId[e1]
    #             bni8 = roNodeId[(e1+1) % elementsCountAround]
    #             nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
    #             # print('nodeIdentifiers', nodeIdentifiers)
    #         # Remap the derivatives
    #         if 0 <= e1 < elementsCountAcross + 1:
    #             eft1 = eftfactory.createEftBasic()
    #             if e1 == 0:
    #                 scalefactors = [-1.0]
    #                 # eft = eftfactory.createEftBasic()
    #                 setEftScaleFactorIds(eft1, [1], [])
    #                 # remapEftNodeValueLabel(eft, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
    #                 # remapEftNodeValueLabel(eft, [7], Node.VALUE_LABEL_D_DS2,
    #                 #                        [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                 # remapEftNodeValueLabel(eft, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
    #                 remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
    #                                        [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                 remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
    #             elif 0 < e1 < elementsCountAcross - 1:
    #                 scalefactors = [-1.0]
    #                 # eft = eftfactory.createEftBasic()
    #                 setEftScaleFactorIds(eft1, [1], [])
    #                 # remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
    #             elif e1 == elementsCountAcross - 1:
    #                 scalefactors = [-1.0]
    #                 # eft = eftfactory.createEftBasic()
    #                 setEftScaleFactorIds(eft1, [1], [])
    #                 # remapEftNodeValueLabel(eft, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 # remapEftNodeValueLabel(eft, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
    #                 # remapEftNodeValueLabel(eft, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
    #                 remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
    #             elif e1 == elementsCountAcross:
    #                 scalefactors = [-1.0]
    #                 # eft = eftfactory.createEftBasic()
    #                 setEftScaleFactorIds(eft1, [1], [])
    #                 remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                 remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
    #             elementtemplateMod.defineField(coordinates, -1, eft1)
    #             elementtemplate1 = elementtemplateMod
    #         elif e1 == elementsCountAround - 1:
    #             eft1 = eftfactory.createEftBasic()
    #             remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #             elementtemplateMod.defineField(coordinates, -1, eft1)
    #             elementtemplate1 = elementtemplateMod
    #         else:
    #             scalefactors = None
    #             eft1 = eft
    #             elementtemplate1 = elementtemplate
    #         element = mesh.createElement(elementIdentifier, elementtemplate1)
    #         result2 = element.setNodesByIdentifier(eft1, nodeIdentifiers)
    #         if scalefactors:
    #             result3 = element.setScaleFactors(eft1, scalefactors)
    #         else:
    #             result3 = '-'
    #         elementIdentifier += 1
    #         for meshGroup in meshGroups:
    #             if meshGroups.index(meshGroup) == 2:
    #                 meshGroup.addElement(element)
    #             elif meshGroups.index(meshGroup) == 3:
    #                 meshGroup.addElement(element)
    #
    # # parent right part
    # for e3 in range(elementsCountThroughWall):
    #     for e1 in range(elementsCountAround): # was paCount
    #         eft = eftStd
    #         elementtemplate = elementtemplateStd
    #         scalefactors = None
    #         bni1 = birNodeId[e1] # was [e3][e1]
    #         bni2 = birNodeId[(e1 + 1) % elementsCountAround]
    #         bni3 = cricNodeId[e1]
    #         bni4 = cricNodeId[(e1 + 1) % elementsCountAround]
    #         if e1 < elementsCountAround // 2:
    #             bni5 = roNodeId[e1]
    #             bni6 = roNodeId[(e1 + 1) % elementsCountAround]
    #             bni7 = cotNodeId[e1]
    #             bni8 = cotNodeId[(e1 + 1) % elementsCountAround]
    #         elif e1 == elementsCountAround // 2:
    #             bni5 = roNodeId[e1]
    #             bni6 = sbNodeId [e1 - elementsCountAround // 2]
    #             bni7 = cotNodeId[e1]
    #             bni8 = csNodeId[e1 - elementsCountAround // 2]
    #         elif elementsCountAround // 2 < e1 < elementsCountAround - 1:
    #             bni5 = sbNodeId[e1 - elementsCountAround // 2 - 1]
    #             bni6 = bni5 + 1
    #             bni7 = csNodeId[e1 - elementsCountAround // 2 - 1]
    #             bni8 = bni7 + 1
    #         elif e1 == elementsCountAround - 1:
    #             bni5 = sbNodeId[-1]
    #             bni6 = roNodeId[0]
    #             bni8 = cotNodeId[0]
    #             bni7 = bni8 - 1
    #         nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
    #         # print('nodeIdentifiers', nodeIdentifiers)
    #         # Remap the derivatives
    #         if e1 == 0:
    #             eft1 = eftfactory.createEftBasic()
    #             remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #             elementtemplateMod.defineField(coordinates, -1, eft1)
    #             elementtemplate1 = elementtemplateMod
    #         elif e1 == elementsCountAround // 2:
    #             eft1 = eftfactory.createEftBasic()
    #             remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
    #             elementtemplateMod.defineField(coordinates, -1, eft1)
    #             elementtemplate1 = elementtemplateMod
    #         elif e1 == elementsCountAround - 1:
    #             eft1 = eftfactory.createEftBasic()
    #             remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
    #             elementtemplateMod.defineField(coordinates, -1, eft1)
    #             elementtemplate1 = elementtemplateMod
    #         else:
    #             scalefactors = None
    #             eft1 = eft
    #             elementtemplate1 = elementtemplate
    #         element = mesh.createElement(elementIdentifier, elementtemplate1)
    #         result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
    #         if scalefactors:
    #             result3 = element.setScaleFactors(eft1, scalefactors)
    #         else:
    #             result3 = '-'
    #         elementIdentifier += 1
    #         for meshGroup in meshGroups:
    #             if meshGroups.index(meshGroup) == 0:
    #                 meshGroup.addElement(element)
    #             elif meshGroups.index(meshGroup) == 3:
    #                 meshGroup.addElement(element)
    #
    # # parent left part
    # for e3 in range(elementsCountThroughWall):
    #     for e1 in range(elementsCountAround):  # was paCount
    #         eft = eftStd
    #         elementtemplate = elementtemplateStd
    #         scalefactors = None
    #         bni1 = bilNodeId[e1]  # was [e3][e1]
    #         bni2 = bilNodeId[(e1 + 1) % elementsCountAround]
    #         bni3 = clicNodeId[e1]
    #         bni4 = clicNodeId[(e1 + 1) % elementsCountAround]
    #         if e1 == 0:
    #             # [169, 170, 199, 200, 177, 190, 210, 209]
    #             bni5 = roNodeId[0]
    #             bni6 = sbNodeId[-1]
    #             bni7 = cotNodeId[e1]
    #             bni8 = csNodeId[e1 + elementsCountAcross -2]
    #         elif 0 < e1 < elementsCountAcross - 1:
    #             # [170, 171, 200, 201, 190, 189, 209, 208]
    #             bni5 = sbNodeId[-e1 + elementsCountAcross - 1]
    #             bni6 = bni5 - 1
    #             bni7 = csNodeId[-e1 + elementsCountAcross - 1]
    #             bni8 = bni7 - 1
    #         elif e1 == elementsCountAcross - 1:
    #             bni5 = sbNodeId[0]
    #             bni6 = roNodeId[elementsCountAround // 2]
    #             bni7 = csNodeId[0]
    #             bni8 = cotNodeId[elementsCountAround // 2]
    #         elif e1 == elementsCountAround - 1:
    #             bni5 = roNodeId[e1 - elementsCountAround // 2 + elementsCountAcross]
    #             bni6 = roNodeId[0]
    #             bni7 = cotNodeId[e1 - elementsCountAround // 2 + elementsCountAcross ]
    #             bni8 = cotNodeId[0]
    #         else:
    #             bni5 = roNodeId[e1 - elementsCountAround // 2 + elementsCountAcross]
    #             bni6 = bni5 + 1
    #             bni7 = cotNodeId[e1 - elementsCountAround // 2 + elementsCountAcross ]
    #             bni8 = bni7 + 1
    #         nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
    #         # print('nodeIdentifiers', nodeIdentifiers)
    #         # Remap the derivatives
    #         if 0 <= e1 < elementsCountAcross + 1:
    #             eft1 = eftfactory.createEftBasic()
    #             if e1 == 0:
    #                 scalefactors = [-1.0]
    #                 # eft1 = eftfactory.createEftBasic()
    #                 setEftScaleFactorIds(eft1, [1], [])
    #                 remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
    #                 remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 # elementtemplateMod.defineField(coordinates, -1, eft1)
    #                 # elementtemplate1 = elementtemplateMod
    #             elif 0 < e1 < elementsCountAcross - 1:
    #                 scalefactors = [-1.0]
    #                 # eft1 = eftfactory.createEftBasic()
    #                 setEftScaleFactorIds(eft1, [1], [])
    #                 remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 # elementtemplateMod.defineField(coordinates, -1, eft1)
    #                 # elementtemplate1 = elementtemplateMod
    #             elif e1 == elementsCountAcross - 1:
    #                 scalefactors = [-1.0]
    #                 # eft1 = eftfactory.createEftBasic()
    #                 setEftScaleFactorIds(eft1, [1], [])
    #                 remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
    #                 remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
    #                 # elementtemplateMod.defineField(coordinates, -1, eft1)
    #                 # elementtemplate1 = elementtemplateMod
    #             elif e1 == elementsCountAcross:
    #                 # eft1 = eftfactory.createEftBasic()
    #                 remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                 # elementtemplateMod.defineField(coordinates, -1, eft1)
    #                 # elementtemplate1 = elementtemplateMod
    #             elementtemplateMod.defineField(coordinates, -1, eft1)
    #             elementtemplate1 = elementtemplateMod
    #         else:
    #             scalefactors = None
    #             eft1 = eft
    #             elementtemplate1 = elementtemplate
    #         element = mesh.createElement(elementIdentifier, elementtemplate1)
    #         result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
    #         if scalefactors:
    #             result3 = element.setScaleFactors(eft1, scalefactors)
    #         else:
    #             result3 = '-'
    #         elementIdentifier += 1
    #         for meshGroup in meshGroups:
    #             if meshGroups.index(meshGroup) == 0:
    #                 meshGroup.addElement(element)
    #             elif meshGroups.index(meshGroup) == 3:
    #                 meshGroup.addElement(element)

    return elementIdentifier


def make_cervix_elements(mesh, coordinates, elementIdentifier, elementsCountAlongCervix, elementsCountAround,
                         elementsCountAcross, elementsCountAroundRightTube, elementsCountAroundLeftTube,
                         cricNodeId, clicNodeId, cotNodeId, csNodeId, useCrossDerivatives, meshGroups=None):

    eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    elementtemplateMod = mesh.createElementtemplate()
    elementtemplateMod.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    newl1 = elementsCountAround * (elementsCountAlongCervix + 1)
    newl2 = (elementsCountAround + elementsCountAroundRightTube) * (elementsCountAlongCervix + 1)
    newl0 = (elementsCountAround + 2 * elementsCountAroundRightTube) * (elementsCountAlongCervix + 1)

    # Cervix elements
    for e3 in range(1):
        for e2 in range(elementsCountAlongCervix - 1):
            # Elements around right canal in cervix
            for e1 in range(elementsCountAroundRightTube):
                bni1 = cricNodeId[e2 * elementsCountAroundRightTube + e1]
                bni2 = cricNodeId[e2 * elementsCountAroundRightTube + (e1 + 1) % elementsCountAroundRightTube]
                # bni3 = bni1 + elementsCountAround * 3 + elementsCountAcross - 1
                # bni4 = bni2 + elementsCountAround * 3 + elementsCountAcross - 1
                bni3 = bni1 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                bni4 = bni2 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                if e1 < elementsCountAround // 2:
                    bni5 = cotNodeId[e2 * elementsCountAround + e1]
                    bni6 = bni5 + 1
                    # bni7 = bni5 + elementsCountAround * 3 + elementsCountAcross - 1
                    bni7 = bni5 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                    bni8 = bni7 + 1
                elif e1 == elementsCountAround // 2:
                    bni5 = cotNodeId[e2 * elementsCountAround + e1]
                    bni6 = csNodeId[e2 * (elementsCountAcross - 1)]
                    # bni7 = bni5 + elementsCountAround * 3 + elementsCountAcross - 1
                    # bni8 = bni6 + elementsCountAround * 3 + elementsCountAcross - 1
                    bni7 = bni5 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                    bni8 = bni6 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                elif elementsCountAround // 2 < e1 < elementsCountAroundLeftTube - 1:
                    bni5 = csNodeId[e2 * (elementsCountAcross - 1) + e1 - elementsCountAround // 2 - 1]
                    bni6 = bni5 + 1
                    # bni7 = bni5 + elementsCountAround * 3 + elementsCountAcross - 1
                    # bni8 = bni6 + elementsCountAround * 3 + elementsCountAcross - 1
                    bni7 = bni5 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                    bni8 = bni6 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                elif e1 == elementsCountAroundLeftTube - 1:
                    bni5 = csNodeId[e2 * (elementsCountAcross - 1) + elementsCountAcross - 2]
                    bni6 = cotNodeId[e2 * elementsCountAround]
                    # bni7 = bni5 + elementsCountAround * 3 + elementsCountAcross - 1
                    # bni8 = bni6 + elementsCountAround * 3 + elementsCountAcross - 1
                    bni7 = bni5 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                    bni8 = bni6 + elementsCountAroundRightTube * 2 + elementsCountAround + elementsCountAcross - 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
                # Remap derivatives for elements at the beginning and end of the septum
                if e1 == elementsCountAround // 2:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif e1 == elementsCountAroundRightTube - 1:
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
                if elementsCountAlongCervix == 2:
                    for meshGroup in meshGroups:
                        if meshGroups.index(meshGroup) == 1:
                                meshGroup.addElement(element)
                        elif meshGroups.index(meshGroup) == 2:
                            meshGroup.addElement(element)
                elif elementsCountAlongCervix > 2:
                    if e2 < elementsCountAlongCervix - 2:
                        for meshGroup in meshGroups:
                            if meshGroups.index(meshGroup) == 0:
                                meshGroup.addElement(element)
                            elif meshGroups.index(meshGroup) == 2:
                                meshGroup.addElement(element)
                    else:
                        for meshGroup in meshGroups:
                            if meshGroups.index(meshGroup) == 1:
                                    meshGroup.addElement(element)
                            elif meshGroups.index(meshGroup) == 2:
                                meshGroup.addElement(element)
                    # else:
                    #     if meshGroups.index(meshGroup) == 1:
                    #             meshGroup.addElement(element)
                    #     elif meshGroups.index(meshGroup) == 2:
                    #         meshGroup.addElement(element)

                    # meshGroup.addElement(element)
                # for meshGroup in meshGroups:
                #     if meshGroups.index(meshGroup) == 0:
                #         meshGroup.addElement(element)
                #     elif meshGroups.index(meshGroup) == 3:
                #         meshGroup.addElement(element)
                elementIdentifier += 1

            # Elements around Left canal in cervix
            # print('clicNodeId', clicNodeId)
            for e1 in range(elementsCountAroundLeftTube):
                if 0 <= e1 < elementsCountAcross:
                    bni1 = clicNodeId[e2 * elementsCountAroundLeftTube + e1]
                    bni2 = clicNodeId[e2 * elementsCountAroundLeftTube + (e1 + 1) % elementsCountAroundLeftTube]
                    # bni3 = bni1 + elementsCountAround * 3 + elementsCountAcross - 1
                    # bni4 = bni2 + elementsCountAround * 3 + elementsCountAcross - 1
                    bni3 = bni1 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    bni4 = bni2 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    if e1 == 0:
                        bni5 = cotNodeId[e2 * elementsCountAround + e1]
                        bni6 = csNodeId[e2 * (elementsCountAcross - 1) + elementsCountAcross - 2]
                        # bni6 = bni5 - 1
                        bni7 = bni5 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                        bni8 = bni6 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    elif e1 == elementsCountAcross - 1:
                        bni5 = csNodeId[e2 * (elementsCountAcross - 1) + elementsCountAcross - 1 - e1]
                        bni6 = cotNodeId[e2 * elementsCountAround + elementsCountAround // 2]
                        bni7 = bni5 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                        bni8 = bni6 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    else:
                        bni5 = csNodeId[e2 * (elementsCountAcross - 1) + elementsCountAcross - 1 - e1]
                        bni6 = bni5 - 1
                        bni7 = bni5 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                        bni8 = bni6 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                    print('nodeIdentifiers', nodeIdentifiers)
                else:
                    bni1 = clicNodeId[e2 * elementsCountAroundLeftTube + e1]
                    bni2 = clicNodeId[e2 * elementsCountAroundLeftTube + (e1 + 1) % elementsCountAroundLeftTube]
                    bni3 = bni1 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    bni4 = bni2 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    bni5 = cotNodeId[e2 * elementsCountAround + e1 - elementsCountAcross + elementsCountAround // 2]
                    if e1 == elementsCountAroundLeftTube - 1:
                       bni6 = cotNodeId[e2 * elementsCountAround]
                    else:
                        bni6 = cotNodeId[e2 * elementsCountAround + (e1 + 1) % elementsCountAroundLeftTube - elementsCountAcross + elementsCountAround // 2]
                    bni7 = bni5 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    bni8 = bni6 + elementsCountAroundLeftTube * 2 + elementsCountAround+ elementsCountAcross - 1
                    nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                    print('nodeIdentifiers', nodeIdentifiers)
                # Remap derivatives for elements across septum
                if e1 == 0:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif 0 < e1 < elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    elementtemplateMod.defineField(coordinates, -1, eft1)
                    elementtemplate1 = elementtemplateMod
                elif e1 == elementsCountAcross - 1:
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
                # for meshGroup in meshGroups:
                #     meshGroup.addElement(element)
                if elementsCountAlongCervix == 2:
                    for meshGroup in meshGroups:
                        if meshGroups.index(meshGroup) == 1:
                                meshGroup.addElement(element)
                        elif meshGroups.index(meshGroup) == 2:
                            meshGroup.addElement(element)
                elif elementsCountAlongCervix > 2:
                    if e2 < elementsCountAlongCervix - 2:
                        for meshGroup in meshGroups:
                            if meshGroups.index(meshGroup) == 0:
                                meshGroup.addElement(element)
                            elif meshGroups.index(meshGroup) == 2:
                                meshGroup.addElement(element)
                    else:
                        for meshGroup in meshGroups:
                            if meshGroups.index(meshGroup) == 1:
                                    meshGroup.addElement(element)
                            elif meshGroups.index(meshGroup) == 2:
                                meshGroup.addElement(element)
                elementIdentifier += 1
    return elementIdentifier


def make_rat_uterus_bifurcation_elements_modified(fm, coordinates, elementIdentifier, elementsCountAround, elementsCountAcross,
                                      elementsCountThroughWall, paNodeId, c1NodeId, c2NodeId, roNodeId, coNodeId,
                                      birNodeId, bilNodeId, cricNodeId, clicNodeId, cotNodeId, sbNodeId, csNodeId, meshGroups=None):

    paCount = len(paNodeId[0])
    c1Count = len(c1NodeId[0])
    c2Count = len(c2NodeId[0])
    pac1Count, pac2Count, c1c2Count = get_tube_bifurcation_connection_elements_counts(paCount, c1Count, c2Count)

    # elementsCountAcross = len((sbNodeId)) + 1

    # elementsCountAcross = 3
    count = elementsCountAround // 2 + elementsCountAcross
    elementsCountAroundRightTube = count
    elementsCountAroundLeftTube = count

    # fm = region.getFieldmodule()
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
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = birNodeId[e1]
                bni4 = birNodeId[(e1 + 1) % c1Count]
                bni5 = bni1 + elementsCountAroundRightTube
                bni6 = bni5 + 1
                # bni7 = bni3 + 2 * elementsCountAround
                bni7 = roNodeId[e1]
                bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            elif e1 == elementsCountAround // 2:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = birNodeId[e1]
                bni4 = birNodeId[(e1 + 1) % c1Count]
                bni5 = bni1 + elementsCountAroundRightTube
                bni6 = bni5 + 1
                # bni7 = bni3 + 2 * elementsCountAround
                # bni8 = coNodeId[0]
                bni7 = bni3 + 2 * elementsCountAroundRightTube
                bni8 = coNodeId[0]
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            else:
                bni1 = c1NodeId[e3][e1]
                bni2 = c1NodeId[e3][(e1 + 1) % c1Count]
                bni3 = birNodeId[e1]
                bni5 = bni1 + elementsCountAroundRightTube
                bni7 = coNodeId[e1 - elementsCountAround // 2 - 1]
                if e1 == elementsCountAroundRightTube - 1:
                    bni4 = birNodeId[0]
                    bni6 = bni5 - elementsCountAroundRightTube + 1
                    bni8 = roNodeId[0]
                else:
                    bni4 = bni3 + 1
                    bni6 = bni5 + 1
                    bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
            if e1 in (0, pac1Count - 1, pac1Count, c1Count - 1):
                eft1 = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                elif e1 == pac1Count - 1:
                    # eft = eftfactory.createEftBasic()
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == pac1Count:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif e1 == c1Count - 1:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
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
            bni3 = bilNodeId[e1]
            bni4 = bilNodeId[(e1 + 1) % elementsCountAroundLeftTube]
            bni5 = bni1 + elementsCountAroundLeftTube
            bni6 = bni2 + elementsCountAroundLeftTube
            if 0 <= e1 < elementsCountAcross:
                if e1 == 0:
                    # bni7 = roNodeId[0]
                    # bni8 = sbNodeId[-1]
                    bni7 = roNodeId[0]
                    bni8 = coNodeId[-1]
                elif e1 == elementsCountAcross - 1:
                    # bni7 = sbNodeId[elementsCountAcross - 1 - e1]
                    # bni8 = roNodeId[elementsCountAround // 2]
                    bni7 = coNodeId[elementsCountAcross - 1 - e1]
                    bni8 = roNodeId[elementsCountAround // 2]
                else:
                    # bni7 = sbNodeId[elementsCountAcross - 1 - e1]
                    # bni8 = bni7 - 1
                    bni7 = coNodeId[elementsCountAcross - 1 - e1]
                    bni8 = bni7 - 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            else:
                # if e1 == elementsCountAcross:
                #     bni7 = roNodeId[elementsCountAround // 2]
                #     bni8 = bni7 + 1
                # else:
                bni7 = roNodeId[elementsCountAround // 2 + e1 - elementsCountAcross]
                if e1 == elementsCountAroundLeftTube - 1:
                    bni8 = roNodeId[0]
                else:
                    bni8 = bni7 + 1
                nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
            if 0 <= e1 < elementsCountAcross + 1:
                eft1 = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    # remapEftNodeValueLabel(eft, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                    # remapEftNodeValueLabel(eft, [7], Node.VALUE_LABEL_D_DS2,
                    #                        [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    # remapEftNodeValueLabel(eft, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif 0 < e1 < elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    # remapEftNodeValueLabel(eft, [7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [7, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                elif e1 == elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    # remapEftNodeValueLabel(eft, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    # remapEftNodeValueLabel(eft, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                    # remapEftNodeValueLabel(eft, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                elif e1 == elementsCountAcross:
                    scalefactors = [-1.0]
                    # eft = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                elementtemplateMod.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateMod
            elif e1 == elementsCountAroundLeftTube - 1:
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
        for e1 in range(elementsCountAroundRightTube): # was paCount
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = birNodeId[e1] # was [e3][e1]
            bni2 = birNodeId[(e1 + 1) % elementsCountAroundRightTube]
            bni3 = cricNodeId[e1]
            bni4 = cricNodeId[(e1 + 1) % elementsCountAroundRightTube]
            if e1 < elementsCountAround // 2:
                bni5 = roNodeId[e1]
                bni6 = roNodeId[(e1 + 1) % elementsCountAround]
                bni7 = cotNodeId[e1]
                bni8 = cotNodeId[(e1 + 1) % elementsCountAround]
            elif e1 == elementsCountAround // 2:
                bni5 = roNodeId[e1]
                bni6 = coNodeId [e1 - elementsCountAround // 2]
                bni7 = cotNodeId[e1]
                bni8 = csNodeId[e1 - elementsCountAround // 2]
            elif elementsCountAround // 2 < e1 < elementsCountAroundRightTube - 1:
                bni5 = coNodeId[e1 - elementsCountAround // 2 - 1]
                bni6 = bni5 + 1
                bni7 = csNodeId[e1 - elementsCountAround // 2 - 1]
                bni8 = bni7 + 1
            elif e1 == elementsCountAroundRightTube - 1:
                bni5 = coNodeId[-1]
                bni6 = roNodeId[0]
                bni8 = cotNodeId[0]
                bni7 = bni8 - 1
            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
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
                remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                elementtemplateMod.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateMod
            elif elementsCountAround // 2 < e1 < elementsCountAroundRightTube - 1:
                # scalefactors = [-1.0]
                eft1 = eftfactory.createEftBasic()
                # setEftScaleFactorIds(eft1, [1], [])
                # remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                remapEftNodeValueLabel(eft1, [5, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                # remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                elementtemplateMod.defineField(coordinates, -1, eft1)
                elementtemplate1 = elementtemplateMod
            elif e1 == elementsCountAroundRightTube - 1:
                eft1 = eftfactory.createEftBasic()
                remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
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
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 0:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    # parent left part
    for e3 in range(elementsCountThroughWall):
        for e1 in range(elementsCountAroundLeftTube):  # was paCount
            eft = eftStd
            elementtemplate = elementtemplateStd
            scalefactors = None
            bni1 = bilNodeId[e1]  # was [e3][e1]
            bni2 = bilNodeId[(e1 + 1) % elementsCountAroundLeftTube]
            bni3 = clicNodeId[e1]
            bni4 = clicNodeId[(e1 + 1) % elementsCountAroundLeftTube]
            if e1 == 0:
                # [137, 138, 164, 165, 146, 155, 175, 174]
                bni5 = roNodeId[0]
                bni6 = coNodeId[-1]
                bni7 = cotNodeId[e1]
                bni8 = csNodeId[e1 + elementsCountAcross -2]
            elif 0 < e1 < elementsCountAcross - 1:
                # [170, 171, 200, 201, 190, 189, 209, 208]
                bni5 = coNodeId[-e1 + elementsCountAcross - 1]
                bni6 = bni5 - 1
                bni7 = csNodeId[-e1 + elementsCountAcross - 1]
                bni8 = bni7 - 1
            elif e1 == elementsCountAcross - 1:
                bni5 = coNodeId[0]
                bni6 = roNodeId[elementsCountAround // 2]
                bni7 = csNodeId[0]
                bni8 = cotNodeId[elementsCountAround // 2]
            # elif e1 == elementsCountAroundLeftTube - 1:
            #     bni5 = roNodeId[e1 - elementsCountAround // 2 + elementsCountAcross]
            #     bni6 = roNodeId[0]
            #     bni7 = cotNodeId[e1 - elementsCountAround // 2 + elementsCountAcross ]
            #     bni8 = cotNodeId[0]
            else:
                bni5 = roNodeId[elementsCountAround // 2 - elementsCountAcross + e1]
                bni7 = cotNodeId[elementsCountAround // 2 - elementsCountAcross + e1]
                if e1 == elementsCountAroundLeftTube - 1:
                    bni6 = roNodeId[0]
                    bni8 = cotNodeId[0]
                else:
                    bni6 = bni5 + 1
                    bni8 = bni7 + 1
            nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
            # print('nodeIdentifiers', nodeIdentifiers)
            # Remap the derivatives
            if 0 <= e1 < elementsCountAcross + 1:
                eft1 = eftfactory.createEftBasic()
                if e1 == 0:
                    scalefactors = [-1.0]
                    # eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft1, [7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [1])])
                    remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    # elementtemplateMod.defineField(coordinates, -1, eft1)
                    # elementtemplate1 = elementtemplateMod
                elif 0 < e1 < elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    # eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5, 6, 7, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS3, [])])
                    # elementtemplateMod.defineField(coordinates, -1, eft1)
                    # elementtemplate1 = elementtemplateMod
                elif e1 == elementsCountAcross - 1:
                    scalefactors = [-1.0]
                    # eft1 = eftfactory.createEftBasic()
                    setEftScaleFactorIds(eft1, [1], [])
                    remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    remapEftNodeValueLabel(eft1, [6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                    remapEftNodeValueLabel(eft1, [8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                    # remapEftNodeValueLabel(eft1, [5, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [1])])
                    # remapEftNodeValueLabel(eft1, [6, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS3, [])])
                    # elementtemplateMod.defineField(coordinates, -1, eft1)
                    # elementtemplate1 = elementtemplateMod
                elif e1 == elementsCountAcross:
                    # eft1 = eftfactory.createEftBasic()
                    remapEftNodeValueLabel(eft1, [5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                    # elementtemplateMod.defineField(coordinates, -1, eft1)
                    # elementtemplate1 = elementtemplateMod
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
            elementIdentifier += 1
            for meshGroup in meshGroups:
                if meshGroups.index(meshGroup) == 0:
                    meshGroup.addElement(element)
                elif meshGroups.index(meshGroup) == 3:
                    meshGroup.addElement(element)

    return elementIdentifier