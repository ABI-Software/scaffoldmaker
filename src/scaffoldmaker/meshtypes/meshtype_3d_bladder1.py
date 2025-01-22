"""
Generates a 3-D bladder mesh along the central line, with variable numbers of elements around , along and through wall.
"""

import copy
import math

from cmlibs.maths.vectorops import angle, set_magnitude, normalize, magnitude, cross, axis_angle_to_rotation_matrix
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.bladder_terms import get_bladder_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, get_nodeset_path_field_parameters
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters


class MeshType_3d_bladder1(Scaffold_base):
    """
    Generates a 3-D bladder mesh with variable numbers of elements around, along the central line, and through wall.
    The bladder is created using a central path as the longitudinal axis of the bladder. Magnitude of D2 and D3 are
    the radii of the bladder in the respective directions.
    """
    centralPathDefaultScaffoldPackages_Bladder = {
        'Cat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                 Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3],
                [[
                    (1, [[1.98482, -0.23213, 0.00000], [-0.27986, 0.04637, 0.00000], [-0.03267, -0.19717, 0.00000],
                         [-0.02887, -0.18604, 0.00000], [0.00000, 0.00000, 0.31783], [0.00000, 0.00000, 0.22083]]),
                    (2, [[1.69774, -0.18588, 0.00000], [-0.29430, 0.04613, 0.00000], [-0.05231, -0.33376, 0.00000],
                         [-0.01041, -0.08714, 0.00000], [0.00000, 0.00000, 0.46417], [0.00000, 0.00000, 0.07185]]),
                    (3, [[1.39622, -0.13993, 0.00000], [-0.29537, 0.04245, 0.00000], [-0.05304, -0.36902, 0.00000],
                         [0.00337, -0.00098, 0.00000], [0.00000, 0.00000, 0.45788], [0.00000, 0.00000, -0.04945]]),
                    (4, [[1.10705, -0.10086, 0.00000], [-0.29137, 0.03953, 0.00000], [-0.04575, -0.33720, 0.00000],
                         [0.01169, 0.05485, 0.00000], [0.00000, 0.00000, 0.36713], [0.00000, 0.00000, -0.09543]]),
                    (5, [[0.81350, -0.06086, 0.00000], [-0.29780, 0.03404, 0.00000], [-0.02960, -0.25897, 0.00000],
                         [0.01563, 0.09020, 0.00000], [0.00000, 0.00000, 0.26695], [0.00000, 0.00000, -0.10351]]),
                    (6, [[0.51160, -0.03295, 0.00000], [-0.28489, 0.02641, 0.00000], [-0.01451, -0.15651, 0.00000],
                         [0.01211, 0.08865, 0.00000], [0.00000, 0.00000, 0.16004], [0.00000, 0.00000, -0.08389]]),
                    (7, [[0.24374, -0.00806, 0.00000], [-0.25596, 0.01611, 0.00000], [-0.00504, -0.08011, 0.00000],
                         [0.00715, 0.05112, 0.00000], [0.00000, 0.00000, 0.09658], [0.00000, 0.00000, -0.04254]]),
                    (8, [[0.00000, 0.00000, 0.00000], [-0.23137, 0.00000, 0.00000], [-0.00000, -0.05191, 0.00000],
                         [0.00294, 0.00527, 0.00000], [0.00000, 0.00000, 0.07300], [0.00000, 0.00000, -0.00462]])
                ]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('dome of the bladder')[0],
                    'ontId': get_bladder_term('dome of the bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-7',
                    'name': get_bladder_term('neck of urinary bladder')[0],
                    'ontId': get_bladder_term('neck of urinary bladder')[1]
                }]
        }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                 Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3],
                [[
                    (1, [[2.00000, 0.00000, 0.00000], [-0.29240, 0.00000, 0.00000], [0.00000, -0.27079, 0.00000],
                         [0.00000, -0.26496, 0.00000], [0.00000, 0.00000, 0.53835], [0.00000, 0.00000, 0.14344]]),
                    (2, [[1.70113, 0.00000, 0.00000], [-0.30534, 0.00000, 0.00000], [0.00000, -0.46797, 0.00000],
                         [0.00000, -0.12939, 0.00000], [0.00000, 0.00000, 0.67711], [0.00000, 0.00000, 0.13406]]),
                    (3, [[1.38931, 0.00000, 0.00000], [-0.31410, 0.00000, 0.00000], [0.00000, -0.52663, 0.00000],
                         [0.00000, -0.00331, 0.00000], [0.00000, 0.00000, 0.80628], [0.00000, 0.00000, 0.00115]]),
                    (4, [[1.07292, 0.00000, 0.00000], [-0.30565, 0.00000, 0.00000], [0.00000, -0.47377, 0.00000],
                         [0.00000, 0.08253, -0.00000], [0.00000, 0.00000, 0.67753], [0.00000, -0.00000, -0.13332]]),
                    (5, [[0.77801, 0.00000, 0.00000], [-0.29228, 0.00000, 0.00000], [0.00000, -0.36358, 0.00000],
                         [0.00000, 0.10935, -0.00000], [0.00000, 0.00000, 0.53994], [0.00000, -0.00000, -0.15490]]),
                    (6, [[0.48836, 0.00000, 0.00000], [-0.26692, 0.00000, 0.00000], [0.00000, -0.25504, 0.00000],
                         [0.00000, 0.09572, -0.00000], [0.00000, 0.00000, 0.36804], [0.00000, -0.00000, -0.16276]]),
                    (7, [[0.24418, 0.00000, 0.00000], [-0.24418, 0.00000, 0.00000], [0.00000, -0.17013, 0.00000],
                         [0.00000, 0.07176, -0.00000], [0.00000, 0.00000, 0.21298], [0.00000, -0.00000, -0.13616]]),
                    (8, [[0.00000, 0.00000, 0.00000], [-0.24418, 0.00000, 0.00000], [0.00000, -0.11153, 0.00000],
                         [0.00000, 0.04545, -0.00000], [0.00000, 0.00000, 0.09572], [0.00000, -0.00000, -0.09837]])
                ]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('dome of the bladder')[0],
                    'ontId': get_bladder_term('dome of the bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-7',
                    'name': get_bladder_term('neck of urinary bladder')[0],
                    'ontId': get_bladder_term('neck of urinary bladder')[1]
                }]
        }),
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                 Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3],
                [[
                    (1, [[1.94475, -0.28705, 0.00000], [-0.34493, 0.16731, 0.00000], [-0.08535, -0.17595, 0.00000],
                         [-0.04402, -0.10886, 0.00000], [0.00000, 0.00000, 0.18665], [0.00000, 0.00000, 0.24107]]),
                    (2, [[1.61917, -0.14378, 0.00000], [-0.30588, 0.11906, 0.00000], [-0.10871, -0.27931, 0.00000],
                         [-0.00271, -0.09786, 0.00000], [0.00000, 0.00000, 0.36886], [0.00000, 0.00000, 0.12335]]),
                    (3, [[1.33473, -0.04705, 0.00000], [-0.30858, 0.07784, 0.00000], [-0.09397, -0.37253, 0.00000],
                         [0.02932, -0.05147, 0.00000], [0.00000, 0.00000, 0.44247], [0.00000, 0.00000, 0.03241]]),
                    (4, [[1.00375, 0.00854, 0.00000], [-0.32781, 0.04202, 0.00000], [-0.04837, -0.37737, 0.00000],
                         [0.03932, 0.04536, 0.00000], [0.00000, 0.00000, 0.42886], [0.00000, 0.00000, -0.06441]]),
                    (5, [[0.67991, 0.03730, 0.00000], [-0.26339, 0.01406, 0.00000], [-0.01513, -0.28340, 0.00000],
                         [0.02873, 0.09787, 0.00000], [0.00000, 0.00000, 0.31526], [0.00000, 0.00000, -0.10213]]),
                    (6, [[0.47773, 0.04096, 0.00000], [-0.21282, -0.01255, 0.00000], [0.01079, -0.18309, 0.00000],
                         [0.01208, 0.09936, 0.00000], [0.00000, 0.00000, 0.22026], [0.00000, 0.00000, -0.09737]]),
                    (7, [[0.25586, 0.01063, 0.00000], [-0.23915, -0.02127, 0.00000], [0.00754, -0.08478, 0.00000],
                         [-0.00525, 0.07201, 0.00000], [0.00000, 0.00000, 0.12027], [0.00000, 0.00000, -0.08509]]),
                    (8, [[0.00000, 0.00000, 0.00000], [-0.27223, 0.00000, 0.00000], [-0.00000, -0.04285, 0.00000],
                         [-0.00982, 0.01186, 0.00000], [0.00000, 0.00000, 0.05222], [0.00000, 0.00000, -0.05101]])
                ]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('dome of the bladder')[0],
                    'ontId': get_bladder_term('dome of the bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-7',
                    'name': get_bladder_term('neck of urinary bladder')[0],
                    'ontId': get_bladder_term('neck of urinary bladder')[1]
                }]
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3],
                [[
                    (1, [[1.88718, 0.56852, 0.00000], [-0.20674, -0.15833, 0.00000], [0.09002, -0.11754, 0.00000],
                         [0.13867, -0.25647, 0.00000], [0.00000, 0.00000, 0.11203], [0.00000, 0.00000, 0.47891]]),
                    (2, [[1.66209, 0.42003, 0.00000], [-0.24279, -0.13815, 0.00000], [0.18168, -0.31927, 0.00000],
                         [0.04465, -0.14700, 0.00000], [0.00000, 0.00000, 0.46092], [0.00000, 0.00000, 0.21887]]),
                    (3, [[1.40235, 0.29407, 0.00000], [-0.26492, -0.11438, 0.00000], [0.17602, -0.40769, 0.00000],
                         [-0.01921, -0.05867, 0.00000], [0.00000, 0.00000, 0.54065], [0.00000, 0.00000, 0.02648]]),
                    (4, [[1.13299, 0.19156, 0.00000], [-0.28314, -0.09291, 0.00000], [0.14328, -0.43666, 0.00000],
                         [-0.03999, -0.00350, 0.00000], [0.00000, 0.00000, 0.51397], [0.00000, 0.00000, -0.04924]]),
                    (5, [[0.83659, 0.10937, 0.00000], [-0.29480, -0.06821, 0.00000], [0.09555, -0.41298, 0.00000],
                         [-0.04673, 0.07253, 0.00000], [0.00000, 0.00000, 0.44065], [0.00000, 0.00000, -0.08938]]),
                    (6, [[0.54418, 0.05488, 0.00000], [-0.28025, -0.04758, 0.00000], [0.04979, -0.29321, 0.00000],
                         [-0.04042, 0.13710, 0.00000], [0.00000, 0.00000, 0.33575], [0.00000, 0.00000, -0.12832]]),
                    (7, [[0.27635, 0.01382, 0.00000], [-0.27252, -0.02763, 0.00000], [0.01423, -0.14032, 0.00000],
                         [-0.02501, 0.12449, 0.00000], [0.00000, 0.00000, 0.18610], [0.00000, 0.00000, -0.14240]]),
                    (8, [[0.00000, 0.00000, 0.00000], [-0.27970, -0.00000, 0.00000], [0.00000, -0.04485, 0.00000],
                         [-0.00345, 0.06647, 0.00000], [0.00000, 0.00000, 0.05111], [0.00000, 0.00000, -0.12759]])
                ]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('dome of the bladder')[0],
                    'ontId': get_bladder_term('dome of the bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-7',
                    'name': get_bladder_term('neck of urinary bladder')[0],
                    'ontId': get_bladder_term('neck of urinary bladder')[1]
                }]
        }),
        'Rat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3],
                [[
                    (1, [[1.98605, -0.21791, 0.00000], [-0.20755, 0.03159, 0.00000], [-0.00760, -0.04989, 0.00000],
                         [-0.02445, -0.15319, 0.00000], [0.00000, 0.00000, 0.12825], [0.00000, 0.00000, 0.29326]]),
                    (2, [[1.75787, -0.18270, 0.00000], [-0.24882, 0.03881, 0.00000], [-0.02816, -0.18050, 0.00000],
                         [-0.01667, -0.10802, 0.00000], [0.00000, 0.00000, 0.36645], [0.00000, 0.00000, 0.18315]]),
                    (3, [[1.48843, -0.14021, 0.00000], [-0.29410, 0.04520, 0.00000], [-0.04024, -0.26183, 0.00000],
                         [-0.00647, -0.06475, 0.00000], [0.00000, 0.00000, 0.48456], [0.00000, 0.00000, 0.09179]]),
                    (4, [[1.16964, -0.09250, 0.00000], [-0.32163, 0.04199, 0.00000], [-0.04008, -0.30698, 0.00000],
                         [0.00483, -0.01528, 0.00000], [0.00000, 0.00000, 0.54524], [0.00000, 0.00000, -0.00502]]),
                    (5, [[0.84529, -0.05632, 0.00000], [-0.31603, 0.03303, 0.00000], [-0.03052, -0.29202, 0.00000],
                         [0.01097, 0.05004, 0.00000], [0.00000, 0.00000, 0.47369], [0.00000, 0.00000, -0.10225]]),
                    (6, [[0.53762, -0.02632, 0.00000], [-0.28085, 0.02450, 0.00000], [-0.01821, -0.20876, 0.00000],
                         [0.01263, 0.09417, 0.00000], [0.00000, 0.00000, 0.34236], [0.00000, 0.00000, -0.13848]]),
                    (7, [[0.28368, -0.00678, 0.00000], [-0.26891, 0.01355, 0.00000], [-0.00532, -0.10560, 0.00000],
                         [0.00931, 0.08494, 0.00000], [0.00000, 0.00000, 0.19799], [0.00000, 0.00000, -0.12929]]),
                    (8, [[0.00000, 0.00000, 0.00000], [-0.29834, 0.00000, 0.00000], [0.00000, -0.04096, 0.00000],
                         [0.00133, 0.04434, 0.00000], [0.00000, 0.00000, 0.08550], [0.00000, 0.00000, -0.09569]])
                ]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('dome of the bladder')[0],
                    'ontId': get_bladder_term('dome of the bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-7',
                    'name': get_bladder_term('neck of urinary bladder')[0],
                    'ontId': get_bladder_term('neck of urinary bladder')[1]
                }]
        }),
        'Material': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 10
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
                 Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3],
                [[
                    (1, [[2.00000, 0.00000, 0.00000], [-0.04772, 0.00000, 0.00000], [0.00000, -0.21530, 0.00000],
                         [0.00000, 0.04601, 0.00000], [0.00000, 0.00000, 0.21530], [0.00000, 0.00000, -0.04601]]),
                    (2, [[1.90255, 0.00000, 0.00000], [-0.14718, 0.00000, 0.00000], [0.00000, -0.21530, 0.00000],
                         [0.00000, -0.04601, 0.00000], [0.00000, 0.00000, 0.21530], [0.00000, 0.00000, 0.04601]]),
                    (3, [[1.70564, 0.00000, 0.00000], [-0.21176, 0.00000, 0.00000], [0.00000, -0.35428, 0.00000],
                         [0.00000, -0.11371, 0.00000], [0.00000, 0.00000, 0.35428], [0.00000, 0.00000, 0.11371]]),
                    (4, [[1.47903, 0.00000, 0.00000], [-0.23206, 0.00000, 0.00000], [0.00000, -0.43890, 0.00000],
                         [0.00000, -0.06591, 0.00000], [0.00000, 0.00000, 0.43890], [0.00000, 0.00000, 0.06591]]),
                    (5, [[1.24152, 0.00000, 0.00000], [-0.23952, 0.00000, 0.00000], [0.00000, -0.48520, 0.00000],
                         [0.00000, -0.03068, 0.00000], [0.00000, 0.00000, 0.48520], [0.00000, 0.00000, 0.03068]]),
                    (6, [[1.00000, 0.00000, 0.00000], [-0.24157, 0.00000, 0.00000], [0.00000, -0.50000, 0.00000],
                         [0.00000, 0.00001, 0.00000], [0.00000, 0.00000, 0.50000], [0.00000, 0.00000, -0.00001]]),
                    (7, [[0.75837, 0.00000, 0.00000], [-0.23957, 0.00000, 0.00000], [0.00000, -0.48518, 0.00000],
                         [0.00000, 0.03070, 0.00000], [0.00000, 0.00000, 0.48518], [0.00000, 0.00000, -0.03070]]),
                    (8, [[0.52087, 0.00000, 0.00000], [-0.23198, 0.00000, 0.00000], [0.00000, -0.43887, 0.00000],
                         [0.00000, 0.06589, 0.00000], [0.00000, 0.00000, 0.43887], [0.00000, 0.00000, -0.06589]]),
                    (9, [[0.29442, 0.00000, 0.00000], [-0.21169, 0.00000, 0.00000], [0.00000, -0.35431, 0.00000],
                         [0.00000, 0.11366, 0.00000], [0.00000, 0.00000, 0.35431], [0.00000, 0.00000, -0.11366]]),
                    (10, [[0.09750, 0.00000, 0.00000], [-0.14721, 0.00000, 0.00000], [0.00000, -0.21534, 0.00000],
                          [0.00000, 0.11291, 0.00000], [0.00000, 0.00000, 0.21534], [0.00000, 0.00000, -0.11291]]),
                    (11, [[0.00000, 0.00000, 0.00000], [-0.04779, 0.00000, 0.00000], [0.00000, -0.11534, 0.00000],
                          [0.00000, 0.08709, 0.00000], [0.00000, 0.00000, 0.11534], [0.00000, 0.00000, -0.08709]])
                ]]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-6',
                    'name': get_bladder_term('dome of the bladder')[0],
                    'ontId': get_bladder_term('dome of the bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7-10',
                    'name': get_bladder_term('neck of urinary bladder')[0],
                    'ontId': get_bladder_term('neck of urinary bladder')[1]
                }]
        }),
    }

    @staticmethod
    def getName():
        return '3D Bladder 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Cat 1',
            'Human 1',
            'Mouse 1',
            'Pig 1',
            'Rat 1',
            'Material']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages_Bladder['Human 1']
        elif 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages_Bladder['Mouse 1']
        elif 'Pig 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages_Bladder['Pig 1']
        elif 'Rat 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages_Bladder['Rat 1']
        elif 'Material' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages_Bladder['Material']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages_Bladder['Cat 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements along dome': 8,
            'Number of elements along neck': 4,
            'Number of elements around': 8,
            'Number of elements through wall': 1,
            'Wall thickness': 0.007,
            'Ureter position around': 0.67,  # should be on the dorsal part (> 0.5). It's a fixed material coordinate and thus not to be edited.
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        if 'Cat 1' in parameterSetName:
            options['Wall thickness'] = 0.007  # was 1.5 * 2.0 / 178.08271473110773
        if 'Human 1' in parameterSetName:
            options['Wall thickness'] = 0.043  # was 3.0 * 2.0 / 327.63
        if 'Mouse 1' in parameterSetName:
            options['Wall thickness'] = 0.008  # was 0.5 * 2.0 / 72.004931756029
        if 'Pig 1' in parameterSetName:
            options['Wall thickness'] = 0.016  # was 2.5 * 2.0 / 505.2045520616655
        if 'Rat 1' in parameterSetName:
            options['Wall thickness'] = 0.006  # was 0.3 * 2.0 / 110.40996641878101
        if 'Material' in parameterSetName:
            options['Number of elements along dome'] = 8
            options['Number of elements along neck'] = 4
            options['Number of elements around'] = 8
            options['Wall thickness'] = 0.02  # an average of the wall thicknesses for above species, rounded.
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Central path',
            'Number of elements along dome',
            'Number of elements along neck',
            'Number of elements around',
            'Number of elements through wall',
            'Wall thickness',
            # 'Ureter position around',  # it's a fixed material coordinate and thus not to be edited.
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages_Bladder.keys())
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
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages_Bladder.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages_Bladder[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        for key in [
            'Number of elements along dome',
            'Number of elements along neck',
            'Number of elements around',
            'Number of elements through wall',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall'
        ]:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements around'] % 2 != 0:
            if options['Number of elements around'] % 4 > 1:
                options['Number of elements around'] += 1
            else:
                options['Number of elements around'] -= 1
        else:
            if options['Number of elements around'] % 4 != 0:
                options['Number of elements around'] += 2
        if options['Ureter position around'] < 0.5:
            options['Ureter position around'] = 0.5  # ureters are on the dorsal part of the bladder

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        centralPath = options['Central path']
        elementsCountAlongDome = options['Number of elements along dome']
        elementsCountAlongNeck = options['Number of elements along neck']
        elementsCountAround = options['Number of elements around']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = False
        useCubicHermiteThroughWall = not (options['Use linear through wall'])

        materialCentralPath = cls.centralPathDefaultScaffoldPackages_Bladder['Material']
        materialWallThickness = 0.02

        elementsCountAlongBladder = elementsCountAlongDome + elementsCountAlongNeck

        fm = region.getFieldmodule()
        fm.beginChange()
        mesh = fm.findMeshByDimension(3)

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Geometric coordinates
        geometricCentralPath = BladderCentralPath(region, centralPath, elementsCountAlongDome, elementsCountAlongNeck)
        xFinal, d1Final, d2Final, d3Final = \
            getBladderCoordinates(elementsCountAlongDome, elementsCountAlongNeck, elementsCountAround,
                                  elementsCountThroughWall, wallThickness, geometricCentralPath)

        # # find biunit scale and rotation data
        # pathLength = geometricCentralPath.bladderCentralPathLength
        # print("path arcLength =", pathLength)
        # print("biunit scale =", 2.0 / pathLength)
        # urethra_d1 = geometricCentralPath.sx_neck_group[1][-1]
        # y_angle = math.degrees(math.atan2(urethra_d1[1], -urethra_d1[0]))
        # print("path y_angle (", urethra_d1[1], ",", -urethra_d1[0], " ) =", y_angle)

        sx_dome_group = geometricCentralPath.sx_dome_group
        sx_neck_group = geometricCentralPath.sx_neck_group
        bladderCentralPathLength = geometricCentralPath.bladderCentralPathLength
        domeLength = geometricCentralPath.domeLength

        # Material coordinates
        tmp_region = region.createRegion()
        materialCentralPath = BladderCentralPath(tmp_region, materialCentralPath, elementsCountAlongDome,
                                                 elementsCountAlongNeck)
        xOrgan, d1Organ, d2Organ, d3Organ = \
            getBladderCoordinates(elementsCountAlongDome, elementsCountAlongNeck, elementsCountAround,
                                  elementsCountThroughWall, materialWallThickness, materialCentralPath)

        # # Obtain the central path coordinates and derivatives for an ellipsoid to put for the "Material" central path
        # bladderLength = 2.0
        # diameter = 1.0
        # centre = [bladderLength / 2, 0.0, 0.0]
        # height = bladderLength
        # poleAxis = [bladderLength / 2, 0.0, 0.0]
        # sideAxis = [0.0, diameter / 2, 0.0]
        # xEllipsoid, d1Ellipsoid, d2Ellipsoid = createEllipsoidPoints(centre, poleAxis, sideAxis, elementsCountAround,
        #                                                              elementsCountAlongBladder, height)
        #
        # pathNodesList = []
        # for n1 in range(elementsCountAlongDome + elementsCountAlongNeck + 1):
        #     n = n1 * elementsCountAround
        #     pathNode = [xEllipsoid[n][0], 0.0, 0.0]
        #     pathNodesList.append(pathNode)
        # print('pathNodesList', pathNodesList)
        #
        # d2path = []
        # d3path = []
        # for n1 in range(1, elementsCountAlongDome + elementsCountAlongNeck):
        #     n = n1 * elementsCountAround
        #     v1 = xEllipsoid[n]
        #     v2 = xEllipsoid[n + elementsCountAround // 2]
        #     v1v2 = [v2[c] - v1[c] for c in range(3)]
        #     mag = magnitude(v1v2)
        #     halfMag = mag / 2
        #     d2path.append([0.0, -halfMag, 0.0])
        #     d3path.append([0.0, 0.0, halfMag])
        # print('d2path', d2path)
        # print('d3path', d3path)

        del tmp_region

        # Create annotation groups for bladder
        bodyGroup = AnnotationGroup(region, get_bladder_term("dome of the bladder"))
        neckGroup = AnnotationGroup(region, get_bladder_term("neck of urinary bladder"))
        bladderGroup = AnnotationGroup(region, get_bladder_term("urinary bladder"))
        bladderDorsalGroup = AnnotationGroup(region, get_bladder_term("dorsal part of bladder"))
        bladderVentralGroup = AnnotationGroup(region, get_bladder_term("ventral part of bladder"))

        elementsCountAlongGroups = [elementsCountAlongDome, elementsCountAlongNeck]
        annotationGroupAlong = [[bladderGroup, bodyGroup], [bladderGroup, neckGroup]]
        annotationGroupsAlong = []
        for i in range(len(elementsCountAlongGroups)):
            elementsCount = elementsCountAlongGroups[i]
            for n in range(elementsCount):
                annotationGroupsAlong.append(annotationGroupAlong[i])

        elementsCountAroundGroups = [elementsCountAround // 4, elementsCountAround // 2, elementsCountAround // 4]
        annotationGroupAround = [[bladderGroup, bladderVentralGroup], [bladderGroup, bladderDorsalGroup],
                                 [bladderGroup, bladderVentralGroup]]
        annotationGroupsAround = []
        for i in range(len(elementsCountAroundGroups)):
            elementsCount = elementsCountAroundGroups[i]
            for n in range(elementsCount):
                annotationGroupsAround.append(annotationGroupAround[i])

        annotationGroupsThroughWall = []
        for i in range(elementsCountThroughWall):
            annotationGroupsThroughWall.append([])

        # Flat coordinates
        xFlat = d1Flat = d2Flat = []
        # urethraOpeningRadius = magnitude(sx_neck_group[2][-1])
        # xFlat, d1Flat, d2Flat = obtainBladderFlatNodes(elementsCountAlongBladder, elementsCountAround,
        #                                                elementsCountThroughWall, xFinal, d1Final, d2Final,
        #                                                bladderCentralPathLength, urethraOpeningRadius, wallThickness)

        # Create nodes and elements
        bladderCoordinatesFieldName = "bladder coordinates"
        nodeIdentifier, elementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
            region, xFinal, d1Final, d2Final, d3Final, xFlat, d1Flat, d2Flat, xOrgan, d1Organ, d2Organ,
            bladderCoordinatesFieldName, elementsCountAround, elementsCountAlongBladder, elementsCountThroughWall,
            annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier,
            useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd=True)[0:3]

        bladderCoordinates = fm.findFieldByName(bladderCoordinatesFieldName)

        # Define trackSurface to put the ureter markers on
        xTrackSurface = []
        d1TrackSurface = []
        d2TrackSurface = []
        for n2 in range(elementsCountAlongBladder + 1):
            for n1 in range(elementsCountAround):
                idx = n2 * elementsCountAround + elementsCountAround + n1
                xTrackSurface.append(xFinal[idx])
                d1TrackSurface.append(d1Final[idx])
                d2TrackSurface.append(d2Final[idx])
        trackSurfaceBladder = TrackSurface(elementsCountAround, elementsCountAlongBladder,
                                           xTrackSurface, d1TrackSurface, d2TrackSurface, loop1=True)

        # Ureter position
        neckStartPositionAlongFactor = domeLength / bladderCentralPathLength
        ureterPositionAroundFactor = options['Ureter position around'] / 2
        ureter1Position = trackSurfaceBladder.createPositionProportion(ureterPositionAroundFactor,
                                                                       neckStartPositionAlongFactor)
        ureterElementPositionAround = ureter1Position.e1

        # Define markers for apex and ureter junctions with bladder
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

        apexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                       get_bladder_term("apex of urinary bladder"))
        leftUreterGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                             get_bladder_term("left ureter junction with bladder"))
        rightUreterGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_bladder_term("right ureter junction with bladder"))

        idx1 = elementsCountAround * elementsCountThroughWall
        xi1 = [0.0, 0.0, 1.0]
        markerList = []
        markerList.append({"group": apexGroup, "elementId": idx1, "xi": xi1})
        idx2 = elementsCountAlongDome * elementsCountAround * elementsCountThroughWall + ureterElementPositionAround + 1
        xi2 = [ureter1Position.xi1, 0.0, 0.0]
        markerList.append({"group": rightUreterGroup, "elementId": idx2, "xi": xi2})
        idx3 = elementsCountAlongDome * elementsCountAround * elementsCountThroughWall + elementsCountAround - ureterElementPositionAround
        xi3 = [1 - ureter1Position.xi1, 0.0, 0.0]
        markerList.append({"group": leftUreterGroup, "elementId": idx3, "xi": xi3})

        bladderNodesetGroup = bladderGroup.getNodesetGroup(nodes)
        for marker in markerList:
            annotationGroup = marker["group"]
            element = mesh.findElementByIdentifier(marker["elementId"])
            # create marker points with element locations
            markerNode = annotationGroup.createMarkerNode(nodeIdentifier, element=element, xi=marker["xi"])
            # calculate bladder coordinates automatically from element:xi
            annotationGroup.setMarkerMaterialCoordinates(bladderCoordinates)
            bladderNodesetGroup.addNode(markerNode)
            nodeIdentifier += 1

        fm.endChange()
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
        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("dome of the bladder"))
        neckGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("neck of urinary bladder"))
        urinaryBladderGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("urinary bladder"))
        bladderVentralGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("ventral part of bladder"))
        bladderDorsalGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("dorsal part of bladder"))

        mesh2d = fm.findMeshByDimension(2)

        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))

        is_body = bodyGroup.getGroup()
        is_body_serosa = fm.createFieldAnd(is_body, is_exterior_face_xi3_1)
        is_body_lumen = fm.createFieldAnd(is_body, is_exterior_face_xi3_0)

        is_neck = neckGroup.getGroup()
        is_neck_serosa = fm.createFieldAnd(is_neck, is_exterior_face_xi3_1)
        is_neck_lumen = fm.createFieldAnd(is_neck, is_exterior_face_xi3_0)

        is_urinaryBladder = urinaryBladderGroup.getGroup()
        is_urinaryBladder_serosa = fm.createFieldAnd(is_urinaryBladder, is_exterior_face_xi3_1)
        is_urinaryBladder_lumen = fm.createFieldAnd(is_urinaryBladder, is_exterior_face_xi3_0)

        is_dorsal_bladder = bladderDorsalGroup.getGroup()
        is_ventral_bladder = bladderVentralGroup.getGroup()

        serosaOfUrinaryBladder = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("serosa of urinary bladder"))
        serosaOfUrinaryBladder.getMeshGroup(mesh2d).addElementsConditional(is_urinaryBladder_serosa)
        lumenOfUrinaryBladder = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("bladder lumen"))
        lumenOfUrinaryBladder.getMeshGroup(mesh2d).addElementsConditional(is_urinaryBladder_lumen)

        serosaOfBody = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("serosa of body of urinary bladder"))
        serosaOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_serosa)
        lumenOfBody = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("lumen of body of urinary bladder"))
        lumenOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_lumen)

        serosaOfNeck = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                          get_bladder_term("serosa of neck of urinary bladder"))
        serosaOfNeck.getMeshGroup(mesh2d).addElementsConditional(is_neck_serosa)
        lumenOfNeck = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_bladder_term("lumen of neck of urinary bladder"))
        lumenOfNeck.getMeshGroup(mesh2d).addElementsConditional(is_neck_lumen)

        is_bladder_serosa_dorsal = fm.createFieldAnd(is_urinaryBladder_serosa, is_dorsal_bladder)
        serosaOfBladder_dorsal = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("dorsal part of serosa of urinary bladder"))
        serosaOfBladder_dorsal.getMeshGroup(mesh2d).addElementsConditional(is_bladder_serosa_dorsal)

        is_bladder_serosa_ventral = fm.createFieldAnd(is_urinaryBladder_serosa, is_ventral_bladder)
        serosaOfBladder_ventral = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("ventral part of serosa of urinary bladder"))
        serosaOfBladder_ventral.getMeshGroup(mesh2d).addElementsConditional(is_bladder_serosa_ventral)

        is_bladder_lumen_dorsal = fm.createFieldAnd(is_urinaryBladder_lumen, is_dorsal_bladder)
        lumenOfBladder_dorsal = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("dorsal part of urinary bladder lumen"))
        lumenOfBladder_dorsal.getMeshGroup(mesh2d).addElementsConditional(is_bladder_lumen_dorsal)

        is_bladder_lumen_ventral = fm.createFieldAnd(is_urinaryBladder_lumen, is_ventral_bladder)
        lumenOfBladder_ventral = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("ventral part of urinary bladder lumen"))
        lumenOfBladder_ventral.getMeshGroup(mesh2d).addElementsConditional(is_bladder_lumen_ventral)

        is_body_serosa_dorsal = fm.createFieldAnd(is_body_serosa, is_dorsal_bladder)
        serosaOfBody_dorsal = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("dorsal part of serosa of body of urinary bladder"))
        serosaOfBody_dorsal.getMeshGroup(mesh2d).addElementsConditional(is_body_serosa_dorsal)

        is_body_serosa_ventral = fm.createFieldAnd(is_body_serosa, is_ventral_bladder)
        serosaOfBody_ventral = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("ventral part of serosa of body of urinary bladder"))
        serosaOfBody_ventral.getMeshGroup(mesh2d).addElementsConditional(is_body_serosa_ventral)

        is_body_lumen_dorsal = fm.createFieldAnd(is_body_lumen, is_dorsal_bladder)
        lumenOfBody_dorsal = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("dorsal part of lumen of body of urinary bladder"))
        lumenOfBody_dorsal.getMeshGroup(mesh2d).addElementsConditional(is_body_lumen_dorsal)

        is_body_lumen_ventral = fm.createFieldAnd(is_body_lumen, is_ventral_bladder)
        lumenOfBody_ventral = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("ventral part of lumen of body of urinary bladder"))
        lumenOfBody_ventral.getMeshGroup(mesh2d).addElementsConditional(is_body_lumen_ventral)

        is_neck_serosa_dorsal = fm.createFieldAnd(is_neck_serosa, is_dorsal_bladder)
        serosaOfNeck_dorsal = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("dorsal part of serosa of neck of urinary bladder"))
        serosaOfNeck_dorsal.getMeshGroup(mesh2d).addElementsConditional(is_neck_serosa_dorsal)

        is_neck_serosa_ventral = fm.createFieldAnd(is_neck_serosa, is_ventral_bladder)
        serosaOfNeck_ventral = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("ventral part of serosa of neck of urinary bladder"))
        serosaOfNeck_ventral.getMeshGroup(mesh2d).addElementsConditional(is_neck_serosa_ventral)

        is_neck_lumen_dorsal = fm.createFieldAnd(is_neck_lumen, is_dorsal_bladder)
        lumenOfNeck_dorsal = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("dorsal part of lumen of neck of urinary bladder"))
        lumenOfNeck_dorsal.getMeshGroup(mesh2d).addElementsConditional(is_neck_lumen_dorsal)

        is_neck_lumen_ventral = fm.createFieldAnd(is_neck_lumen, is_ventral_bladder)
        lumenOfNeck_ventral = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("ventral part of lumen of neck of urinary bladder"))
        lumenOfNeck_ventral.getMeshGroup(mesh2d).addElementsConditional(is_neck_lumen_ventral)


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


def findNodesAlongBladderDome(sx_dome_group, elementsCountAround, elementsCountAlongDome):

    # Find apex nodes d2
    d2Apex = []
    d2 = sx_dome_group[2][0]
    for n1 in range(elementsCountAround):
        rotAngle = n1 * 2.0 * math.pi / elementsCountAround
        rotAxis = normalize(cross(normalize(sx_dome_group[2][0]),
                                                        normalize(sx_dome_group[4][0])))
        rotFrame = axis_angle_to_rotation_matrix(rotAxis, rotAngle)
        d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
        d2Apex.append(d2Rot)

    # Create ellipses along dome around the central path
    xEllipses_dome = []
    d1Ellipses_dome = []
    for n in range(1, len(sx_dome_group[0])):
        px, pd1 = createEllipsePoints(sx_dome_group[0][n], 2 * math.pi, sx_dome_group[2][n], sx_dome_group[4][n],
                                      elementsCountAround, startRadians=0.0)
        xEllipses_dome.append(px)
        d1Ellipses_dome.append(pd1)

    # Find d2
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipses_dome)):
            if n2 == 0:
                v1 = sx_dome_group[0][0]
                v2 = xEllipses_dome[n2][n1]
            else:
                v1 = xEllipses_dome[n2 - 1][n1]
                v2 = xEllipses_dome[n2][n1]
            d2 = findDerivativeBetweenPoints(v1, v2)
            if n2 != 0:
                xAlong.append(v1)
            d2Along.append(d2)
        xAlong.append(xEllipses_dome[-1][n1])
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
        d2Raw.append(d2Smoothed)

    # Rearrange d2
    d2Ellipses_dome = []
    for n2 in range(len(xEllipses_dome)):
        d2Around = []
        for n1 in range(elementsCountAround):
            d2 = d2Raw[n1][n2]
            d2Around.append(d2)
        d2Ellipses_dome.append(d2Around)

    # Merge apex and dome
    domeNodes_x = [[sx_dome_group[0][0]] * elementsCountAround] + xEllipses_dome
    domeNodes_d2 = [d2Apex] + d2Ellipses_dome

    # Spread out elements along dome of the bladder
    xRaw = []
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(domeNodes_x)):
            xAlong.append(domeNodes_x[n2][n1])
            d2Along.append(domeNodes_d2[n2][n1])
        xSampledAlongDome, d2SampledAlongDome = interp.sampleCubicHermiteCurves(xAlong, d2Along, elementsCountAlongDome,
                                                                                arcLengthDerivatives=True)[0:2]
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlongDome, d2SampledAlongDome)
        xRaw.append(xSampledAlongDome)
        d2Raw.append(d2Smoothed)

    # Rearrange x and d2
    xSampledDome = []
    d1SampledDome = []
    d2SampledDome = []
    for n2 in range(elementsCountAlongDome + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAround):
            x = xRaw[n1][n2]
            d2 = d2Raw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)
            # Calculate d1
            if n2 > 0:
                v1 = xRaw[n1][n2]
                v2 = xRaw[n1 + 1 if n1 < elementsCountAround - 2 else 0][n2]
                d1 = findDerivativeBetweenPoints(v1, v2)
                d1Around.append(d1)
            else:
                d1Around.append(d2Raw[n2][0])
        if n2 > 0:
            d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        else:
            d1Smoothed = d1Around
        xSampledDome.append(xAround)
        d1SampledDome.append(d1Smoothed)
        d2SampledDome.append(d2Around)

    return xSampledDome, d1SampledDome, d2SampledDome


def findNodesAlongBladderNeck(sx_dome_group, sx_neck_group, d2SampledDome, domeSegmentLength, neckSegmentLength,
                              elementsCountAround, elementsCountAlongNeck):

    # Transition
    transitLength = (domeSegmentLength + neckSegmentLength) / 2
    addingLength = transitLength - neckSegmentLength
    # if transitLength <= neckSegmentLength:
    #     e = 1
    #     xi = transitLength / neckSegmentLength
    # else:
    #     e = int(transitLength / neckSegmentLength) + 1
    #     xi = transitLength / (e * neckSegmentLength)
    # xTransition = interp.interpolateCubicHermite(sx_dome_group[0][-1], sx_dome_group[1][-1], sx_neck_group[0][e], sx_neck_group[1][e], xi)
    # # d1Transition = interp.interpolateCubicHermiteDerivative(sx_dome_group[0][-1], sx_dome_group[1][-1], sx_neck_group[0][e], sx_neck_group[1][e], xi)
    # d2Transition = interp.interpolateCubicHermite(sx_dome_group[2][-1], sx_dome_group[3][-1], sx_neck_group[2][e], sx_neck_group[3][e], xi)
    # d3Transition = interp.interpolateCubicHermite(sx_dome_group[4][-1], sx_dome_group[5][-1], sx_neck_group[4][e], sx_neck_group[5][e], xi)
    # px_transit, pd1_transit = createEllipsePoints(xTransition, 2 * math.pi, d2Transition, d3Transition, elementsCountAround,
    #                                               startRadians=0.0)

    # Create ellipses along neck around the central path
    xEllipses_neck = []
    d1Ellipses_neck = []
    for n in range(0, len(sx_neck_group[0])):
        px, pd1 = createEllipsePoints(sx_neck_group[0][n], 2 * math.pi, sx_neck_group[2][n], sx_neck_group[4][n],
                                      elementsCountAround, startRadians=0.0)
        xEllipses_neck.append(px)
        d1Ellipses_neck.append(pd1)

    # Find d2
    d2Raw = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipses_neck)):
            if n2 == 0:
                v2 = xEllipses_neck[n2][n1]
                d2 = d2SampledDome[-1][n1]
            elif (n2 == len(xEllipses_neck) - 1):
                v1 = xEllipses_neck[n2 - 1][n1]
                v2 = xEllipses_neck[n2][n1]
                d2vec = findDerivativeBetweenPoints(v1, v2)
                d2mag = magnitude(d2vec)
                d2dir = normalize(sx_neck_group[1][-1])
                d2 = set_magnitude(d2dir, d2mag)
            else:
                v1 = xEllipses_neck[n2 - 1][n1]
                v2 = xEllipses_neck[n2][n1]
                d2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v2)
            d2Along.append(d2)
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along, fixEndDirection=True)
        d2Raw.append(d2Smoothed)

    # Rearrange d2
    d2Ellipses_neck = []
    for n2 in range(len(xEllipses_neck)):
        d2Around = []
        for n1 in range(elementsCountAround):
            d2 = d2Raw[n1][n2]
            d2Around.append(d2)
        d2Ellipses_neck.append(d2Around)

    # Spread out elements along neck of the bladder
    xRawNeck = []
    d2RawNeck = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipses_neck)):
            xAlong.append(xEllipses_neck[n2][n1])
            d2Along.append(d2Ellipses_neck[n2][n1])
        xSampledAlongNeck, d2SampledAlongNeck = interp.sampleCubicHermiteCurves(xAlong, d2Along, elementsCountAlongNeck,
                                                                                addLengthStart=addingLength,
                                                                                arcLengthDerivatives=True)[0:2]
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlongNeck, d2SampledAlongNeck)
        xRawNeck.append(xSampledAlongNeck)
        d2RawNeck.append(d2SampledAlongNeck)

    # Rearrange x and d2
    xSampledNeck = []
    d1SampledNeck = []
    d2SampledNeck = []
    for n2 in range(elementsCountAlongNeck + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAround):
            x = xRawNeck[n1][n2]
            d2 = d2RawNeck[n1][n2]
            xAround.append(x)
            d2Around.append(d2)
            # Calculate d1
            v1 = xRawNeck[n1][n2]
            v2 = xRawNeck[n1 + 1 if n1 < elementsCountAround - 1 else 0][n2]
            d1 = findDerivativeBetweenPoints(v1, v2)
            d1Around.append(d1)
        d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        if n2 != 0:
            xSampledNeck.append(xAround)
            d1SampledNeck.append(d1Smoothed)
            d2SampledNeck.append(d2Around)

    px_transit = xSampledNeck[0]

    return xSampledNeck, d1SampledNeck, d2SampledNeck, px_transit


def obtainBladderFlatNodes(elementsCountAlongBladder, elementsCountAround, elementsCountThroughWall, xFinal, d1Final,
                           d2Final, bladderLength, urethraOpeningRadius, bladderWallThickness):
    """
    Calculates flat coordinates for the bladder when it is opened into a flat preparation.
    :param elementsCountAlongBladder: Number of elements along bladder.
    :param elementsCountAround: Number of elements around bladder.
    :param elementsCountThroughWall: Number of elements through wall.
    :param xFinal, d1Final, d2Final: Coordinates and derivatives of all nodes of coordinates field.
    :param bladderLength: Bladder length along the central path.
    :param urethraOpeningRadius: Major diameter of the bladder opening to urethra.
    :param bladderWallThickness: Thickness of wall.
    :return: Coordinates and derivatives of flat coordinates field.
    """
    # Calculate ellipses circumference along the bladder inner layer(lumen)
    circumferenceList = []
    minorNodeAlong_x = []
    minorNodeAlong_d2 = []
    for n1 in range(1, elementsCountAlongBladder + 1):
        xNodesOnEllipse = []
        d1NodesOnEllipse = []
        idMinor = (n1 - 1) * elementsCountAround * (elementsCountThroughWall + 1) + elementsCountThroughWall + 1
        minorNodeAlong_x.append(xFinal[idMinor])
        minorNodeAlong_d2.append(d2Final[idMinor])
        for n2 in range(elementsCountAround):
            id = n2 + idMinor
            x = xFinal[id]
            d1 = d1Final[id]
            xNodesOnEllipse.append(x)
            d1NodesOnEllipse.append(d1)
        xNodesOnEllipse.append(xFinal[idMinor])
        d1NodesOnEllipse.append(d1Final[idMinor])
        circumference = interp.getCubicHermiteCurvesLength(xNodesOnEllipse, d1NodesOnEllipse)
        circumferenceList.append(circumference)
    maxCircumference = max(circumferenceList)

    # Find the angle at the bottom of the bladder neck
    v1 = [0.0, 0.0, bladderLength]
    v2 = [urethraOpeningRadius, 0.0, bladderLength]
    alpha = angle(v1, v2)

    # Find apex to urethra arcLength in minor radius
    xApexInner = xFinal[0]
    d2ApexInner = d2Final[0]
    minorNodeAlong_x.insert(0, xApexInner)
    minorNodeAlong_d2.insert(0, d2ApexInner)
    minorarcLength = interp.getCubicHermiteCurvesLength(minorNodeAlong_x, minorNodeAlong_d2)

    # calculate nodes coordinates based on Hammer projection formulation
    xfnList = []
    angleAlongUnit = (math.pi - alpha) / elementsCountAlongBladder
    angleAroundUnit = 2 * math.pi / elementsCountAround
    for n2 in range(1, elementsCountAlongBladder + 1):
        for n1 in range(elementsCountAround + 1):
            phi = -math.pi / 2 - alpha + n2 * angleAlongUnit
            if n1 < elementsCountAround // 2 + 1:
                theta = n1 * angleAroundUnit
            else:
                theta = math.pi - n1 * angleAroundUnit
            t = math.sqrt(1 + math.cos(phi) * math.cos(theta / 2))
            xScale = maxCircumference / 2
            yScale = bladderLength / 2
            x = [xScale * math.cos(phi) * math.sin(theta / 2) / t,
                 yScale * (math.sin(phi) / t + 1),
                 0.0]
            xfnList.append(x)

    # Rearrange the nodes
    xfnListRearranged = []
    for n2 in range(elementsCountAlongBladder):
        for n1 in range(elementsCountAround + 1):
            if n1 < elementsCountAround // 2 + 1:
                nodeIndex = n2 * (elementsCountAround + 1) + elementsCountAround // 2 - n1
            else:
                nodeIndex = n2 * (elementsCountAround + 1) + n1
            xfnListRearranged.append(xfnList[nodeIndex])

    # Define d1 for flat nodes
    d1fnListRearranged = []
    for n2 in range(elementsCountAlongBladder):
        for n1 in range(elementsCountAround + 1):
            if n1 == elementsCountAround:
                id = (elementsCountAround + 1) * n2
                d1 = [d1fnListRearranged[id][0], -d1fnListRearranged[id][1], d1fnListRearranged[id][2]]
            else:
                v1 = xfnListRearranged[(elementsCountAround + 1) * n2 + n1]
                v2 = xfnListRearranged[(elementsCountAround + 1) * n2 + n1 + 1]
                d1 = [v2[c] - v1[c] for c in range(3)]
            d1fnListRearranged.append(d1)
    # Smoothing the derivatives
    smoothed_d1 = []
    for n2 in range(elementsCountAlongBladder):
        xLineSmoothing = []
        d1LineSmoothing = []
        for n1 in range(elementsCountAround + 1):
            xLineSmoothing.append(xfnListRearranged[n2 * (elementsCountAround + 1) + n1])
            d1LineSmoothing.append(d1fnListRearranged[n2 * (elementsCountAround + 1) + n1])
        smd1 = interp.smoothCubicHermiteDerivativesLine(xLineSmoothing, d1LineSmoothing, fixAllDirections=False,
                                                        fixStartDerivative=True, fixEndDerivative=True,
                                                        fixStartDirection=False, fixEndDirection=False)
        smoothed_d1 += smd1

    # Define d2 for flat nodes
    # Create lines from top go down the bladder
    xNodesToDown = []
    d2NodesToDown = []
    for n2 in range(elementsCountAround + 1):
        x = [xfnListRearranged[n1 * (elementsCountAround + 1) + n2] for n1 in range(0, elementsCountAlongBladder)]
        xNodesToDown += x
        for n1 in range(1, elementsCountAlongBladder + 1):
            if n1 == elementsCountAlongBladder:
                d2 = [0.0, minorarcLength / elementsCountAlongBladder, 0.0]
            else:
                v1 = xNodesToDown[elementsCountAlongBladder * n2 + n1 - 1]
                v2 = xNodesToDown[elementsCountAlongBladder * n2 + n1]
                d2 = [v2[c] - v1[c] for c in range(3)]
            d2NodesToDown.append(d2)
    # Smoothing the derivatives
    smoothed_d2 = []
    for n1 in range(elementsCountAround + 1):
        lineSmoothingNodes = []
        lineSmoothingNodes_d2 = []
        for n2 in range(elementsCountAlongBladder):
            lineSmoothingNodes.append(xNodesToDown[n1 * elementsCountAlongBladder + n2])
            lineSmoothingNodes_d2.append(d2NodesToDown[n1 * elementsCountAlongBladder + n2])
        smd2 = interp.smoothCubicHermiteDerivativesLine(lineSmoothingNodes, lineSmoothingNodes_d2,
                                                        fixAllDirections=False, fixStartDerivative=False,
                                                        fixEndDerivative=True, fixStartDirection=False,
                                                        fixEndDirection=False)
        smoothed_d2 += smd2
    # Re-arrange the derivatives order
    d2fnListRearranged = []
    for n2 in range(elementsCountAlongBladder):
        for n1 in range(elementsCountAround + 1):
            rd2 = smoothed_d2[n1 * elementsCountAlongBladder + n2]
            d2fnListRearranged.append(rd2)

    # Create the nodes through the wall
    xflatList = []
    d1flatList = []
    d2flatList = []
    for n2 in range(elementsCountThroughWall + 1):
        for n1 in range(len(xfnListRearranged)):
            x = [xfnListRearranged[n1][0], xfnListRearranged[n1][1], n2 * bladderWallThickness / elementsCountThroughWall]
            d1 = smoothed_d1[n1]
            d2 = d2fnListRearranged[n1]
            xflatList.append(x)
            d1flatList.append(d1)
            d2flatList.append(d2)

    # Apex derivatives
    v1 = [0.0, 0.0, 0.0]
    v2 = xfnListRearranged[0]
    v3 = xfnListRearranged[elementsCountAround // 2]
    v21 = [v2[c] - v1[c] for c in range(3)]
    v31 = [v3[c] - v1[c] for c in range(3)]
    d1Mag = magnitude(v21)
    d2Mag = magnitude(v31)

    # Add apex nodes to the list
    xFlat = []
    d1Flat = []
    d2Flat = []
    for n1 in range(elementsCountThroughWall + 1):
        xApex = [0.0, 0.0, 0.0 + n1 * bladderWallThickness / elementsCountThroughWall]
        xFlat.append(xApex)
        d1Flat.append([-d1Mag, 0.0, 0.0])
        d2Flat.append([0.0, d2Mag, 0.0])

    # Re-arrange the nodes
    for n3 in range(1, elementsCountAlongBladder + 1):
        for n2 in range(elementsCountThroughWall + 1):
            for n1 in range(elementsCountAround + 1):
                i = (n3 - 1) * (elementsCountAround + 1) + n2 * (elementsCountAround + 1) * elementsCountAlongBladder + n1
                xFlat.append(xflatList[i])
                d1Flat.append(d1flatList[i])
                d2Flat.append(d2flatList[i])

    return xFlat, d1Flat, d2Flat


class BladderCentralPath:
    """
    Generates sampled central path for bladder scaffold.
    """
    def __init__(self, region, centralPath, elementsCountAlongDome, elementsCountAlongNeck):
        """
        :param region: Zinc region to define model in.
        :param centralPath: Central path subscaffold from meshtype_1d_path1
        :param elementsCountAlongDome, elementsCountAlongNeck: Nummber of elements along bladder dome and neck.
        """

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        tmpFieldmodule = tmpRegion.getFieldmodule()
        tmpCoordinates = tmpFieldmodule.findFieldByName('coordinates')
        tmpNodes = tmpFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        tmpDomeGroup = tmpFieldmodule.findFieldByName('dome of the bladder').castGroup()
        tmpDomeNodesetGroup = tmpDomeGroup.getNodesetGroup(tmpNodes)
        cx_dome, cd1_dome, cd2_dome, cd3_dome, cd12_dome, cd13_dome = get_nodeset_path_field_parameters(
            tmpDomeNodesetGroup, tmpCoordinates, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                  Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                  Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3])
        tmpNeckGroup = tmpFieldmodule.findFieldByName('neck of urinary bladder').castGroup()
        tmpNeckNodesetGroup = tmpNeckGroup.getNodesetGroup(tmpNodes)
        cx_neck, cd1_neck, cd2_neck, cd3_neck, cd12_neck, cd13_neck = get_nodeset_path_field_parameters(
            tmpNeckNodesetGroup, tmpCoordinates, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                  Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                  Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3])
        # Find arcLength
        # Dome
        self.domeLength = 0.0
        elementsCountInDome = len(cx_dome) - 1
        for e in range(elementsCountInDome):
            arcLength = interp.getCubicHermiteArcLength(cx_dome[e], cd1_dome[e],
                                                        cx_dome[e + 1], cd1_dome[e + 1])
            self.domeLength += arcLength
        self.domeSegmentLength = self.domeLength / elementsCountAlongDome
        # Neck
        self.neckLength = 0.0
        elementsCountInNeck = len(cx_neck) - 1
        for e in range(elementsCountInNeck):
            arcLength = interp.getCubicHermiteArcLength(cx_neck[e], cd1_neck[e],
                                                        cx_neck[e + 1], cd1_neck[e + 1])
            self.neckLength += arcLength
        self.neckSegmentLength = self.neckLength / elementsCountAlongNeck
        self.bladderCentralPathLength = self.domeLength + self.neckLength
        del tmpRegion

        # Sample central path
        # Dome
        sx_dome, sd1_dome, se_dome, sxi_dome, ssf_dome = interp.sampleCubicHermiteCurves(cx_dome, cd1_dome, len(cx_dome))
        sd2_dome, sd12_dome = interp.interpolateSampleCubicHermite(cd2_dome, cd12_dome, se_dome, sxi_dome, ssf_dome)
        sd3_dome, sd13_dome = interp.interpolateSampleCubicHermite(cd3_dome, cd13_dome, se_dome, sxi_dome, ssf_dome)
        # Neck
        sx_neck, sd1_neck, se_neck, sxi_neck, ssf_neck = interp.sampleCubicHermiteCurves(cx_neck, cd1_neck, elementsCountAlongNeck)
        sd2_neck, sd12_neck = interp.interpolateSampleCubicHermite(cd2_neck, cd12_neck, se_neck, sxi_neck, ssf_neck)
        sd3_neck, sd13_neck = interp.interpolateSampleCubicHermite(cd3_neck, cd13_neck, se_neck, sxi_neck, ssf_neck)

        # Find nodes of bladder dome and neck
        # self.sx_dome_group = [sx_dome, sd1_dome, sd2_dome, sd12_dome, sd3_dome, sd13_dome]
        # self.sx_neck_group = [sx_neck, sd1_neck, sd2_neck, sd12_neck, sd3_neck, sd13_neck]
        self.sx_dome_group = [cx_dome, cd1_dome, cd2_dome, cd12_dome, cd3_dome, cd13_dome]
        self.sx_neck_group = [cx_neck, cd1_neck, cd2_neck, cd12_neck, cd3_neck, cd13_neck]


def getBladderCoordinates(elementsCountAlongDome, elementsCountAlongNeck, elementsCountAround, elementsCountThroughWall,
                          wallThickness, centralPath):

    sx_dome_group = centralPath.sx_dome_group
    sx_neck_group = centralPath.sx_neck_group
    domeSegmentLength = centralPath.domeSegmentLength
    neckSegmentLength = centralPath.neckSegmentLength

    elementsCountAlongBladder = elementsCountAlongDome + elementsCountAlongNeck

    xDome, d1Dome, d2Dome = findNodesAlongBladderDome(sx_dome_group, elementsCountAround, elementsCountAlongDome)

    xNeck, d1Neck, d2Neck, px_transit = \
        findNodesAlongBladderNeck(sx_dome_group, sx_neck_group, d2Dome, domeSegmentLength, neckSegmentLength,
                                  elementsCountAround, elementsCountAlongNeck)

    # Replace d2 for the last layer of the dome based on the transition nodes
    d2Around = []
    for n1 in range(elementsCountAround):
        v1 = xDome[-1][n1]
        v2 = px_transit[n1]
        d2 = findDerivativeBetweenPoints(v1, v2)
        # xAlong.append(v1)
        d2Around.append(d2)
    d2Dome = d2Dome[:-1] + [d2Around]

    xSampledAll = xDome + xNeck
    d1SampledAll = d1Dome + d1Neck
    d2SampledAll = d2Dome + d2Neck

    # # Smoothing d2 from apex to down the neck
    # d2Raw = []
    # for n1 in range(elementsCountAround):
    #     xAlong = []
    #     d2Along = []
    #     for n2 in range(elementsCountAlongBladder + 1):
    #         x = xSampledAll[n2][n1]
    #         d2 = d2SampledAll[n2][n1]
    #         xAlong.append(x)
    #         d2Along.append(d2)
    #     d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along, fixAllDirections=True)
    #     d2Raw.append(d2Smoothed)
    #
    # # Rearrange d2
    # d2SampledAll = []
    # for n2 in range(elementsCountAlongBladder + 1):
    #     d2Around = []
    #     for n1 in range(elementsCountAround):
    #         d2 = d2Raw[n1][n2]
    #         d2Around.append(d2)
    #     d2SampledAll.append(d2Around)

    d3UnitOuter = []
    for n2 in range(1, elementsCountAlongBladder + 1):
        d3Around = []
        for n1 in range(elementsCountAround):
            d3Around.append(normalize(
                cross(normalize(d1SampledAll[n2][n1]), normalize(d2SampledAll[n2][n1]))))
        d3UnitOuter.append(d3Around)

    # Inner nodes
    xInner = []
    d1Inner = []
    d2Inner = []
    d3Inner = []
    for n2 in range(elementsCountAlongBladder + 1):
        for n1 in range(elementsCountAround):
            x = xSampledAll[n2][n1]
            d1 = d1SampledAll[n2][n1]
            d2 = d2SampledAll[n2][n1]
            d3 = d3UnitOuter[n2 - 1][n1]
            xInner.append(x)
            d1Inner.append(d1)
            d2Inner.append(d2)
            d3Inner.append(d3)

    # Create inner layers from the outer nodes
    wallThicknessList = [wallThickness] * (elementsCountAlongBladder + 1)
    relativeThicknessList = []
    transitElementList = [0] * elementsCountAround
    xList, d1List, d2List, d3List, curvatureList = \
        tubemesh.extrudeSurfaceCoordinates(xInner, d1Inner, d2Inner, d3Inner, wallThicknessList, relativeThicknessList,
                                         elementsCountAround, elementsCountAlongBladder, elementsCountThroughWall,
                                         transitElementList, outward=False)[0:5]

    # Deal with multiple nodes at the start point for closed proximal end
    n = elementsCountAround * (elementsCountThroughWall + 1)
    xApexOuter = xList[n - 1]
    # Arclength between apex point and corresponding point on next face
    mag = interp.getCubicHermiteArcLength(xList[n - 1], d2List[n - 1], xList[2 * n - 1],
                                          d2List[2 * n - 1])
    d2ApexOuter = set_magnitude(sx_dome_group[2][0], mag)

    d1ApexOuter = cross(sx_dome_group[1][0], d2ApexOuter)
    d1ApexOuter = set_magnitude(d1ApexOuter, mag)
    d3ApexUnit = normalize(
        cross(normalize(d1ApexOuter), normalize(d2ApexOuter)))
    d3ApexOuter = [d3ApexUnit[c] * wallThickness / elementsCountThroughWall for c in range(3)]

    # Final nodes on the bladder
    xFinal = []
    d1Final = []
    d2Final = []
    d3Final = []
    for n3 in range(elementsCountThroughWall + 1):
        xApex = [xApexOuter[c] -
                 d3ApexUnit[c] * wallThickness / elementsCountThroughWall * (elementsCountThroughWall - n3) for c in range(3)]
        xFinal.append(xApex)
        d1Final.append(d1ApexOuter)
        d2Final.append(d2ApexOuter)
        d3Final.append(d3ApexOuter)

    xFinal += xList[(elementsCountThroughWall + 1) * elementsCountAround:]
    d1Final += d1List[(elementsCountThroughWall + 1) * elementsCountAround:]
    d2Final += d2List[(elementsCountThroughWall + 1) * elementsCountAround:]
    d3Final += d3List[(elementsCountThroughWall + 1) * elementsCountAround:]

    return xFinal, d1Final, d2Final, d3Final
