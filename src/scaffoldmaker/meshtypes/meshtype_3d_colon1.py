"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy

from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.colon_terms import get_colon_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1, ColonSegmentTubeMeshInnerPoints, \
    getTeniaColi, createFlatCoordinatesTeniaColi, createColonCoordinatesTeniaColi, createNodesAndElementsTeniaColi
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, \
    get_nodeset_path_field_parameters


class MeshType_3d_colon1(Scaffold_base):
    '''
    Generates a 3-D colon mesh with variable numbers
    of elements around, along the central line, and through wall.
    The colon is created by a function that generates a colon
    segment and uses tubemesh to map the segment along a central
    line profile.
    '''

    centralPathDefaultScaffoldPackages = {
        'Cattle 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 52
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    (-1, [[ -245.3, 444.6, -49.1 ], [ -267.7,  -53.1, -20.2 ], [   0.0,   0.0,  15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -380.3, 484.8, -45.0 ], [   24.5,  102.7,  15.7 ], [   0.0,   0.0,  15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -298.1, 510.4, -36.8 ], [   73.6,    9.9, -16.4 ], [   0.0,   0.0,  15.0 ], [ -10.8,  -2.7,  -6.5 ]]),
                    (-1, [[ -213.1, 527.9, -22.5 ], [   -1.0,  -10.8, 125.6 ], [ -22.6,  -5.6,   1.4 ], [  -5.3, -12.9, -22.2 ]]),
                    (-1, [[ -315.5, 570.2,  18.9 ], [ -107.9,    9.3,  21.9 ], [  -3.9, -18.1, -31.1 ], [  14.1,   6.3, -14.4 ]]),
                    (-1, [[ -417.4, 555.0,  14.6 ], [  -83.0,  -41.3,  -0.8 ], [   6.1,   1.7, -31.0 ], [   1.4,  10.6,   7.1 ]]),
                    (-1, [[ -497.3, 488.9,  13.6 ], [  -44.6,  -81.6,  10.0 ], [  -1.7,  -2.6, -21.1 ], [  -3.0,  -0.8,   8.0 ]]),
                    (-1, [[ -527.0, 392.5,   2.7 ], [   47.4,  -82.0,  -7.9 ], [   0.0,   0.0, -15.0 ], [   1.6,   2.9,   1.1 ]]),
                    (-1, [[ -461.2, 345.9,  -0.8 ], [   56.9,  -44.5,   2.4 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -415.6, 293.8,   3.9 ], [   93.2,  -62.6,   3.1 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -232.2, 264.9,   0.2 ], [  140.1,   58.2,  -1.0 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -168.4, 357.2,   1.3 ], [   10.1,   78.6,  -3.2 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -185.3, 419.1,  -0.7 ], [  -45.1,   57.1,  -0.9 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -253.2, 466.7,  -0.3 ], [  -63.4,   24.7,   0.2 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -323.8, 482.5,   0.1 ], [  -68.2,    2.9,  -1.2 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -387.5, 485.4,  -0.2 ], [  -44.2,  -17.1,  -1.0 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -435.6, 433.5,   3.3 ], [    3.4, -109.5,   1.4 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -370.6, 376.3,  -1.1 ], [   66.9,  -29.2,  -0.9 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -313.0, 357.9,  -0.1 ], [   40.0,  -33.5,   9.6 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -259.2, 340.7,   2.1 ], [   48.9,    6.4,   1.4 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -246.5, 380.3,  -0.8 ], [  -29.7,   33.6,  -0.7 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -297.3, 387.1,   0.6 ], [  -59.7,   12.6,  -0.0 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -340.2, 415.6,  -1.0 ], [  -86.2,   28.9,  -2.9 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -398.3, 443.1,  -0.1 ], [   10.6,   82.1,  -2.6 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -329.8, 449.1,  -2.1 ], [   53.2,   14.0,  -0.5 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -251.3, 425.9,  -0.3 ], [   43.9,  -19.3,   0.0 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -209.1, 390.6,   0.0 ], [   26.0,  -38.8,   0.9 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -207.8, 350.8,   1.4 ], [   -9.4,  -43.6,   1.8 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -245.8, 299.4,   7.6 ], [  -70.3,  -36.0,   1.4 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -345.3, 304.1,   3.1 ], [ -100.2,   27.9,  -1.9 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -418.4, 361.1,  -0.2 ], [  -57.8,   55.8,  -1.7 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -479.2, 415.6,   2.2 ], [   -8.8,   73.1,  -1.6 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -439.6, 495.7,  -2.1 ], [   61.1,   57.1,  -1.3 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -361.6, 522.6,  -3.0 ], [   78.6,    9.9,   0.2 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -270.1, 506.5,  -3.8 ], [  103.6,  -33.3,   1.0 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -148.9, 441.4,  -2.1 ], [   79.7,  -91.5,   2.8 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -130.9, 313.3,   4.0 ], [   -4.0, -107.2,   3.1 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -183.9, 251.0,   3.8 ], [  -65.5,  -60.2,   3.6 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -280.3, 213.0,   3.4 ], [ -165.1,  -18.6,   0.1 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -400.8, 247.5,   6.8 ], [ -127.1,   36.8,   1.3 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -530.5, 290.7,   5.2 ], [  -89.0,   86.5,   0.3 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -568.8, 392.3,   6.9 ], [  -77.4,   67.7,  -5.5 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.5 ]]),
                    (-1, [[ -511.2, 535.1,   2.5 ], [   86.2,  111.4,  -1.0 ], [   0.0,   0.0, -15.0 ], [   0.0,   0.0,   0.0 ]]),
                    (-1, [[ -405.0, 601.7,   6.4 ], [  143.6,   52.2,   2.6 ], [   0.0,   0.0, -15.0 ], [   3.2,  -3.1,  -1.4 ]]),
                    (-1, [[ -238.8, 615.9,  16.6 ], [   63.3,   -9.1,  19.1 ], [   8.0,  -7.7, -18.4 ], [ -11.9,   7.9,   0.2 ]]),
                    (-1, [[ -146.2, 605.9,  36.5 ], [   49.3,   -9.9, -50.6 ], [ -23.9,  -3.7, -21.6 ], [  -9.9,   8.8,  14.5 ]]),
                    (-1, [[ -218.4, 585.3,  -2.0 ], [ -124.0,    0.4, -37.5 ], [  -9.2,  13.6,  21.7 ], [   2.9,  -0.2,  23.1 ]]),
                    (-1, [[ -376.3, 579.6, -40.8 ], [ -189.2,  -50.7,  -8.8 ], [  -5.0,   3.3,  23.8 ], [   9.6,  -7.2,  10.6 ]]),
                    (-1, [[ -557.9, 493.9, -24.9 ], [  -30.3,   24.1, 152.8 ], [  27.2,  29.2,   3.0 ], [   3.8, -10.1, -29.6 ]]),
                    (-1, [[ -484.8, 594.4,   0.7 ], [  132.7,   97.0,   3.5 ], [  12.4,  -4.5, -32.3 ], [ -12.3, -22.3, -22.5 ]]),
                    (-1, [[ -318.1, 641.9,  -8.5 ], [  166.7,   17.6,   5.5 ], [   3.0, -13.0, -39.4 ], [  -8.3,  -3.3,  -0.9 ]]),
                    (-1, [[ -158.3, 634.7,  -1.9 ], [  176.5,  -14.0,  10.8 ], [  -4.3, -11.5, -34.6 ], [  -3.1,   2.7,   5.3 ]]),
                    (-1, [[   32.7, 611.7,  13.6 ], [  205.5,  -32.2,  20.0 ], [  -2.4,  -7.3, -28.7 ], [   6.9,   5.6,   6.4 ]])]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-8',
                    'name': get_colon_term('right colon')[0],
                    'ontId': get_colon_term('right colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9-46',
                    'name': get_colon_term('transverse colon')[0],
                    'ontId': get_colon_term('transverse colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '47-52',
                    'name': get_colon_term('left colon')[0],
                    'ontId': get_colon_term('left colon')[1]
                }]
        }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    (1, [[   0.0,   0.0, 0.0 ], [ -50.7,  178.2, 0.0 ], [ -24.0,  -6.0, -12.0 ], [ -14.0,  -1.0, -12.0 ]]),
                    (2, [[ -47.4, 188.6, 0.0 ], [ -19.3,  177.1, 0.0 ], [ -22.0,  -4.0,  -8.0 ], [  -4.0,  19.0,  22.0 ]]),
                    (3, [[  -4.4, 396.5, 0.0 ], [ 206.0,   40.1, 0.0 ], [ -10.0,  20.0,   8.0 ], [  -6.0,   0.0,  51.0 ]]),
                    (4, [[ 130.0, 384.1, 0.0 ], [ 130.8,  -40.5, 0.0 ], [  -5.0,   4.0,  29.0 ], [   0.0,   1.0,  24.0 ]]),
                    (5, [[ 279.4, 383.0, 0.0 ], [ 118.0,   48.7, 0.0 ], [  -2.0,  10.0,  22.0 ], [   5.0,  25.0, -20.0 ]]),
                    (6, [[ 443.9, 390.8, 0.0 ], [ 111.3,  -97.0, 0.0 ], [  10.0,  17.0,   6.0 ], [   1.0,  -6.0, -35.0 ]]),
                    (7, [[ 475.2, 168.0, 0.0 ], [  -0.8, -112.4, 0.0 ], [  20.0,   0.0, -20.0 ], [  15.0,  -1.0, -10.0 ]]),
                    (8, [[ 432.6, -32.3, 0.0 ], [ -90.5,  -59.0, 0.0 ], [   6.0,  -9.0, -14.0 ], [   8.0, -11.0, -13.0 ]]),
                    (9, [[ 272.4,   7.5, 0.0 ], [ -79.0,   47.4, 0.0 ], [   1.0, -11.0, -18.0 ], [   4.0, -12.0, -12.0 ]])
                ]),

            'userAnnotationGroups': [
            {
                '_AnnotationGroup': True,
                'dimension': 1,
                'identifierRanges': '1-2',
                'name': get_colon_term('ascending colon')[0],
                'ontId': get_colon_term('ascending colon')[1]
            },
            {
                '_AnnotationGroup': True,
                'dimension': 1,
                'identifierRanges': '3-5',
                'name': get_colon_term('transverse colon')[0],
                'ontId': get_colon_term('transverse colon')[1]
            },
            {
                '_AnnotationGroup': True,
                'dimension': 1,
                'identifierRanges': '6-8',
                'name': get_colon_term('descending colon')[0],
                'ontId': get_colon_term('descending colon')[1]
            }]
        }),
        'Human 2': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    (1, [[   0.0,   0.0,    0.0 ], [ -34.7,  104.1, -18.1 ], [ -24.0,  -6.0, -12.0 ], [ -14.0,  -1.0, -12.0 ]]),
                    (2, [[ -34.5,  114.0, -18.1 ], [   1.2,   86.6,  -3.4 ], [ -22.0,  -4.0,  -8.0 ], [  -4.0,  19.0,  22.0 ]]),
                    (3, [[ -19.1,  218.5,   5.5 ], [  78.7,   -7.1,  94.5 ], [ -10.0,  20.0,   8.0 ], [  -6.0,   0.0,  51.0 ]]),
                    (4, [[  82.5,  189.1,  94.2 ], [  84.5,    7.1,  71.6 ], [  -5.0,   4.0,  29.0 ], [   0.0,   1.0,  24.0 ]]),
                    (5, [[ 226.6,  218.7,  85.7 ], [  95.0,   91.3, -58.5 ], [  -2.0,  10.0,  22.0 ], [   5.0,  25.0, -20.0 ]]),
                    (6, [[ 325.5,  381.7, -57.9 ], [ 229.2,  -66.7, -20.4 ], [  10.0,  17.0,   6.0 ], [   1.0,  -6.0, -35.0 ]]),
                    (7, [[ 354.0,  105.3, -24.4 ], [  -6.3, -143.7,  20.3 ], [  20.0,   0.0, -20.0 ], [  15.0,  -1.0, -10.0 ]]),
                    (8, [[ 296.5, -121.2,  -0.6 ], [ -90.5,  -59.0,   0.0 ], [   6.0,  -9.0, -14.0 ], [   8.0, -11.0, -13.0 ]]),
                    (9, [[ 169.8,  -73.4, -33.5 ], [ -72.2,   43.4, -27.4 ], [   1.0, -11.0, -18.0 ], [   4.0, -12.0, -12.0 ]])
                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_colon_term('ascending colon')[0],
                    'ontId': get_colon_term('ascending colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-5',
                    'name': get_colon_term('transverse colon')[0],
                    'ontId': get_colon_term('transverse colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-8',
                    'name': get_colon_term('descending colon')[0],
                    'ontId': get_colon_term('descending colon')[1]
                }]
        }),
        'Human 3': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 48
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[-87.21, -111.06, 890.54], [-4.75, 0.41, 12.39], [21.05, -2.83, 7.13], [2.46, -0.39, -2.95],
                     [2.79, 22.12, 0.41], [1.83, 0.46, -4.31]],
                    [[-89.99, -110.57, 902.65], [-2.05, 0.53, 12.72], [22.40, -2.97, 3.74], [0.07, 0.14, -2.81],
                     [3.10, 22.66, -0.58], [0.43, 0.10, -3.97]],
                    [[-91.25, -110.01, 915.92], [-0.55, 0.53, 13.39], [21.27, -2.54, 0.97], [-1.46, 0.55, -1.92],
                     [2.57, 21.27, -0.82], [-1.42, 0.54, -2.13]],
                    [[-91.08, -109.52, 929.39], [0.13, 0.65, 13.47], [19.48, -1.87, -0.09], [-0.53, 0.22, -1.06],
                     [1.86, 19.46, -1.05], [-0.53, 0.21, -1.17]],
                    [[-91.00, -108.72, 942.85], [0.83, 0.67, 13.37], [20.22, -2.11, -1.15], [-0.88, 0.27, -1.85],
                     [2.05, 20.22, -1.27], [-0.92, 0.27, -2.03]],
                    [[-89.43, -108.18, 956.08], [2.81, 0.32, 13.11], [17.75, -1.34, -3.78], [-1.96, 0.62, -1.35],
                     [1.22, 18.09, -0.76], [-2.00, 0.62, -1.49]],
                    [[-85.39, -108.09, 968.92], [3.10, 0.44, 13.08], [16.31, -0.87, -3.84], [-0.76, 0.32, 0.66],
                     [0.72, 16.68, -0.80], [-0.73, 0.31, 0.71]],
                    [[-83.24, -107.30, 982.14], [2.05, 0.85, 13.23], [16.22, -0.70, -2.47], [-0.24, 0.09, 0.57],
                     [0.54, 16.36, -1.26], [-0.21, 0.08, 0.63]],
                    [[-81.29, -106.40, 995.37], [2.32, 1.34, 13.29], [15.83, -0.69, -2.70], [0.14, 0.24, -0.11],
                     [0.41, 15.95, -1.86], [0.05, 0.31, -0.57]],
                    [[-78.60, -104.62, 1008.65], [2.14, -1.50, 13.27], [16.51, -0.22, -2.69], [-0.91, 1.85, -0.22],
                     [0.74, 16.43, 3.01], [-0.95, 1.98, -0.15]],
                    [[-77.17, -109.18, 1020.82], [4.07, -9.07, 9.68], [14.05, 2.96, -3.13], [-2.15, 4.06, 0.80],
                     [-0.01, 10.40, 10.89], [-1.97, 3.89, 1.43]],
                    [[-71.19, -121.03, 1025.58], [7.48, -11.42, 0.49], [12.23, 7.97, -0.98], [-1.27, 2.69, 3.17],
                     [0.00, 1.72, 15.29], [-1.27, 2.72, 3.23]],
                    [[-63.35, -130.56, 1022.18], [7.99, -9.44, -3.40], [11.46, 8.59, 3.10], [-0.81, 0.70, 2.42],
                     [0.00, -4.74, 14.55], [-0.91, 0.85, 2.43]],
                    [[-55.27, -139.84, 1018.80], [8.88, -8.59, -3.55], [10.61, 9.36, 3.87], [-1.07, 0.87, 0.81],
                     [0.00, -5.35, 14.33], [-1.09, 0.87, 0.81]],
                    [[-45.66, -147.66, 1015.10], [9.90, -7.39, -3.36], [9.31, 10.33, 4.71], [-1.63, 1.38, 0.03],
                     [-0.01, -5.79, 14.11], [-1.62, 1.36, 0.02]],
                    [[-35.52, -154.60, 1012.08], [11.08, -6.08, -1.98], [7.36, 12.12, 3.95], [-1.99, 1.15, 0.25],
                     [-0.01, -4.33, 14.74], [-1.97, 1.13, 0.24]],
                    [[-23.71, -159.68, 1011.20], [11.97, -4.32, -1.79], [5.34, 12.62, 5.23], [-1.74, -0.44, 2.51],
                     [0.01, -5.34, 14.28], [-1.79, -0.44, 2.52]],
                    [[-11.72, -163.20, 1008.53], [12.26, -2.58, -2.06], [3.87, 11.25, 8.97], [-2.16, -0.11, 1.44],
                     [-0.01, -8.84, 12.28], [-2.14, -0.12, 1.43]],
                    [[0.67, -164.82, 1007.10], [12.71, -0.73, -0.48], [1.02, 12.38, 8.14], [-2.67, 0.76, -0.81],
                     [-0.01, -7.75, 13.06], [-2.58, 0.75, -0.81]],
                    [[13.53, -164.63, 1007.59], [12.75, 0.98, 0.85], [-1.47, 12.77, 7.36], [-1.75, 0.42, -1.04],
                     [-0.28, -7.05, 13.49], [-1.71, 0.43, -1.04]],
                    [[26.09, -162.88, 1008.79], [12.63, 1.81, 1.23], [-2.49, 13.21, 6.05], [-2.26, -0.16, -1.14],
                     [-0.39, -5.90, 14.09], [-2.23, -0.14, -1.14]],
                    [[38.73, -161.01, 1010.04], [11.76, 4.53, 2.77], [-6.00, 12.45, 5.07], [-3.08, -0.42, -3.34],
                     [-0.85, -5.62, 14.19], [-3.12, -0.44, -3.35]],
                    [[48.79, -154.12, 1014.15], [9.86, 7.10, 4.16], [-8.65, 12.37, -0.62], [-1.93, -0.52, -2.09],
                     [-4.15, -2.21, 15.11], [-1.97, -0.53, -2.09]],
                    [[58.41, -146.84, 1018.34], [9.34, 7.75, 4.35], [-9.86, 11.41, 0.85], [-1.82, -1.10, 1.29],
                     [-3.19, -3.77, 15.02], [-1.80, -1.06, 1.30]],
                    [[67.43, -138.64, 1022.84], [7.99, 8.77, 4.67], [-12.30, 10.17, 1.95], [-2.50, -1.97, -0.34],
                     [-2.29, -5.48, 15.73], [-2.59, -2.16, -0.41]],
                    [[74.34, -129.42, 1027.61], [5.08, 9.97, 3.83], [-14.87, 7.50, 0.22], [0.03, -2.06, 0.33],
                     [-2.33, -5.10, 16.50], [-1.53, -1.17, 0.49]],
                    [[77.56, -119.25, 1030.43], [5.04, 10.78, -0.51], [-12.53, 5.98, 2.38], [0.10, 1.23, -1.77],
                     [3.79, 1.01, 17.93], [0.40, 1.97, -2.97]],
                    [[83.86, -109.53, 1026.21], [6.93, 6.51, -9.45], [-14.97, 10.33, -3.86], [-2.81, 1.30, -3.77],
                     [6.42, 13.14, 12.02], [-1.24, 0.30, -4.60]],
                    [[89.55, -108.42, 1013.56], [3.91, 0.86, -13.16], [-18.20, 8.23, -4.87], [-1.85, -1.68, 0.87],
                     [7.52, 18.71, 3.84], [-2.14, -1.95, 1.17]],
                    [[91.64, -107.83, 1000.37], [1.78, 0.55, -13.27], [-18.73, 6.95, -2.22], [-0.80, -1.12, 1.76],
                     [6.79, 18.82, 1.85], [-0.86, -1.12, 1.92]],
                    [[93.10, -107.33, 987.04], [1.01, 0.30, -13.42], [-19.80, 6.00, -1.36], [-0.56, -0.65, 1.35],
                     [5.95, 19.84, 0.98], [-0.58, -0.66, 1.49]],
                    [[93.66, -107.24, 973.55], [-0.02, 1.05, -13.41], [-19.85, 5.64, 0.47], [0.93, -2.25, 1.33],
                     [5.66, 19.79, 1.71], [0.92, -2.27, 1.47]],
                    [[93.06, -105.25, 960.32], [-0.77, 2.25, -13.20], [-17.96, 1.50, 1.30], [0.92, -2.34, 0.15],
                     [1.70, 17.72, 3.23], [0.93, -2.36, 0.17]],
                    [[92.12, -102.74, 947.16], [-0.42, 2.82, -13.16], [-18.00, 0.97, 0.78], [0.74, -0.11, -0.85],
                     [1.11, 17.57, 4.15], [0.73, -0.09, -0.94]],
                    [[92.23, -99.62, 934.04], [0.47, 1.78, -13.24], [-16.47, 1.28, -0.41], [1.01, 0.01, -1.14],
                     [1.22, 16.31, 2.49], [1.03, 0.02, -1.27]],
                    [[93.04, -99.17, 920.86], [1.21, -0.74, -13.33], [-15.97, 1.00, -1.51], [0.28, -0.61, 0.54],
                     [1.07, 16.00, -0.87], [0.28, -0.69, -0.02]],
                    [[94.65, -101.12, 907.56], [-0.61, -1.62, -13.47], [-15.91, 0.05, 0.71], [-0.21, -1.55, 3.28],
                     [-0.04, 15.82, -2.00], [-0.49, -1.71, 2.07]],
                    [[91.84, -102.35, 894.35], [-3.47, -1.64, -11.95], [-16.38, -2.10, 5.05], [0.91, -2.20, 3.65],
                     [-2.64, 16.98, -1.59], [-0.73, -2.58, 3.87]],
                    [[87.88, -104.31, 883.75], [-5.21, -2.13, -10.32], [-14.30, -4.34, 8.12], [2.73, -2.87, 2.50],
                     [-5.97, 18.48, -0.99], [0.46, -3.19, 4.53]],
                    [[81.45, -106.59, 873.90], [-6.99, -2.14, -9.27], [-10.90, -7.86, 10.04], [0.33, -3.15, 2.69],
                     [-9.26, 19.09, 1.96], [2.23, -2.67, 3.37]],
                    [[73.98, -108.57, 865.27], [-8.31, -0.22, -8.52], [-13.53, -10.65, 13.47], [-1.28, -1.76, 1.51],
                     [-8.01, 18.53, 7.34], [3.18, -1.46, 2.05]],
                    [[65.06, -106.92, 857.15], [-8.97, 1.81, -7.69], [-13.39, -11.32, 12.96], [1.28, 0.43, -0.09],
                     [-4.01, 17.55, 11.89], [2.48, 0.16, 0.52]],
                    [[56.07, -104.97, 849.91], [-10.00, 1.92, -6.88], [-11.02, -9.85, 13.27], [2.49, 5.00, 1.79],
                     [-2.03, 16.12, 10.97], [2.28, 4.69, 1.10]],
                    [[45.11, -103.12, 843.53], [-11.00, 4.81, -5.26], [-8.40, -0.98, 16.68], [2.56, 4.66, 2.88],
                     [6.12, 16.27, 4.28], [2.42, 4.35, 0.87]],
                    [[34.86, -95.59, 839.81], [-9.59, 8.68, -2.67], [-5.91, -0.68, 19.02], [0.94, -0.41, 0.56],
                     [11.03, 11.58, 3.43], [0.00, 0.16, -0.54]],
                    [[26.14, -85.99, 838.23], [-8.26, 10.16, -2.00], [-6.49, -1.78, 17.83], [-0.32, -0.12, -0.85],
                     [12.91, 9.40, 5.55], [-1.05, -0.47, 0.41]],
                    [[18.39, -75.31, 835.82], [-8.42, 10.08, -2.66], [-6.54, -0.88, 17.32], [-0.56, 0.61, 0.80],
                     [14.15, 8.98, 5.99], [-0.36, -0.28, 1.59]],
                    [[9.35, -65.88, 832.92], [-11.90, 13.24, -4.28], [-7.61, -0.56, 19.43], [-0.85, 0.48, 1.28],
                     [15.89, 10.29, 6.56], [-1.39, -0.78, 1.04]],
                    [[-5.12, -48.68, 827.11], [-17.03, 21.16, -7.34], [-8.07, 0.18, 19.25], [-0.08, 1.01, -1.64],
                     [15.89, 10.29, 6.56], [-1.39, -0.78, 1.04]]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-48',
                    'name': get_colon_term('colon')[0],
                    'ontId': get_colon_term('colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-10',
                    'name': get_colon_term('ascending colon')[0],
                    'ontId': get_colon_term('ascending colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-26',
                    'name': get_colon_term('transverse colon')[0],
                    'ontId': get_colon_term('transverse colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '27-48',
                    'name': get_colon_term('descending colon')[0],
                    'ontId': get_colon_term('descending colon')[1]
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
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    (1, [[   0.0,   0.0,  0.0 ], [  6.0, 12.0,  -2.0 ], [ 2.0,  1.0,  2.0 ], [ 6.0, 0.0, 3.0 ]]),
                    (2, [[  -2.0,  11.0, -3.0 ], [ -8.0,  4.0,   9.0 ], [ 2.0,  2.0,  1.0 ], [ 0.0, 1.0, 2.0 ]]),
                    (3, [[  -3.0,   2.0,  3.0 ], [ -4.0, -8.0,   0.0 ], [ 2.0, -1.0,  2.0 ], [ 1.0, 0.0, 2.0 ]]),
                    (4, [[ -11.0,  -3.0, -4.0 ], [ -8.0, -3.0,  -7.0 ], [ 1.0, -2.0,  1.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (5, [[ -16.0,  -4.0,  0.0 ], [  4.0, -3.0,  14.0 ], [ 1.0, -3.0,  0.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (6, [[  -7.0,  -8.0,  0.0 ], [  5.0, -1.0, -14.0 ], [ 0.0, -3.0,  0.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (7, [[  -1.0,  -6.0, -1.0 ], [  2.0, -2.0,   9.0 ], [ 1.0, -3.0, -1.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (8, [[  -2.0, -14.0,  5.0 ], [ -2.0, -4.0,   2.0 ], [ 1.0, -2.0, -2.0 ], [ 0.0, 0.0, 0.5 ]])
                ]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_colon_term('right colon')[0],
                    'ontId': get_colon_term('right colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-4',
                    'name': get_colon_term('transverse colon')[0],
                    'ontId': get_colon_term('transverse colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-7',
                    'name': get_colon_term('left colon')[0],
                    'ontId': get_colon_term('left colon')[1]
                }]
        }),
        'Mouse 2': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    (1, [[   0.0,  0.0,   0.0 ], [  0.0,  0.0,  13.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (2, [[   0.0,  0.0,  13.0 ], [  0.0,  2.0,  28.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (3, [[ -14.0, -2.0,  13.0 ], [  0.0, -3.0, -19.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (4, [[ -14.0, -1.0, -10.0 ], [  1.0,  1.0, -17.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (5, [[ -14.0,  0.0, -28.0 ], [  0.0,  0.0, -11.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ]])
                ]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_colon_term('right colon')[0],
                    'ontId': get_colon_term('right colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3',
                    'name': get_colon_term('transverse colon')[0],
                    'ontId': get_colon_term('transverse colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4',
                    'name': get_colon_term('left colon')[0],
                    'ontId': get_colon_term('left colon')[1]
                }]
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 39
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    (-1, [[   -7.2,   83.3,  -20.7 ], [  -65.2,   -8.1,   7.6 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -68.5,   52.8,   -9.6 ], [  -40.1,  -36.1,  10.7 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -97.4,  -26.3,    5.7 ], [   18.0,  -93.2,  13.7 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -56.8,  -90.5,   14.1 ], [   65.5,  -41.4,   7.3 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   48.9, -100.8,   24.0 ], [  112.2,   40.1,  19.0 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  114.8,  -12.6,   38.7 ], [    8.2,   96.1,  14.2 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   60.3,   83.5,   43.7 ], [ -108.7,   54.1,  22.4 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -41.2,   90.7,   56.3 ], [  -89.0,  -32.4,  14.4 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[ -107.9,   -9.7,   76.6 ], [   11.1,  -94.4,  11.3 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -57.3,  -91.9,   81.3 ], [   71.2,  -31.2,   5.7 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   51.2,  -89.4,   97.2 ], [   99.1,   55.4,  12.9 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   91.6,    9.3,  103.6 ], [    4.7,   51.2,   3.4 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   61.6,  111.8,  109.6 ], [  -85.2,   46.1,   2.6 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -54.6,   91.9,  129.4 ], [  -92.7,  -55.0,  14.5 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[ -109.0,    5.6,  156.9 ], [   23.6, -108.2,  27.7 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -59.1,  -62.5,  170.8 ], [   74.0,  -20.1,  14.4 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   23.5,  -53.2,  179.7 ], [   84.6,   47.0,   6.9 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   62.3,   30.1,  187.5 ], [  -12.8,   58.0,   0.8 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   22.4,   45.2,  181.1 ], [  -23.6,  -34.5,  -7.4 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   -1.9,    4.9,  180.5 ], [  -41.3,  -30.9,   7.5 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -45.1,  -12.6,  194.4 ], [  -40.5,   -4.6,   6.9 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -71.7,   -2.2,  197.2 ], [  -25.2,   35.8,  -6.8 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -65.8,   42.1,  182.3 ], [   26.6,   37.6, -15.6 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -14.1,   81.2,  163.5 ], [   41.0,   10.3,  -9.5 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   61.7,   86.1,  156.4 ], [   77.9,  -40.7,   8.9 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   92.9,   20.5,  150.3 ], [    0.0,  -73.3,  -5.2 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   48.9,  -65.0,  142.8 ], [  -82.8,  -80.0,  -1.9 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -54.3,  -90.8,  134.0 ], [  -60.1,   26.4,  -8.2 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -89.9,   11.2,  115.0 ], [   34.9,  125.1, -27.9 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -17.4,   74.2,   91.1 ], [   78.8,   19.1, -15.4 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   43.4,   50.2,   73.3 ], [   30.2,  -36.0,  -9.9 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   62.4,   -5.1,   63.5 ], [   10.9,  -54.2,  -2.7 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   32.7,  -51.7,   56.1 ], [  -38.6,  -29.8,  -8.1 ], [   0.0,   0.0,  5.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -38.1,  -28.6,   46.8 ], [  -62.5,   82.6, -19.2 ], [   4.0,   6.8, 13.2 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[    5.7,   40.4,   22.4 ], [  144.8,   18.6, -20.5 ], [   4.3, -13.3, 12.3 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   53.0,  -14.7,   -4.1 ], [   -6.0,  -25.7, -46.7 ], [ -13.5, -19.4,  4.0 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   24.8,   -0.4,  -48.8 ], [  -13.4,   23.9, -30.6 ], [ -17.4, -15.1,  6.4 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -20.9,   15.3,  -77.9 ], [  -51.2,  -30.6,  21.1 ], [   7.1,  -7.7, 12.2 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[  -47.6,   33.9, -112.2 ], [   32.6,   30.7, -27.8 ], [ -12.8,  -0.5, -1.6 ], [ 0.0, 0.0, 0.5 ]]),
                    (-1, [[   19.6,   96.0, -167.5 ], [   19.9,   19.1, -18.4 ], [ -12.6,   1.3, -8.1 ], [ 0.0, 0.0, 0.5 ]])
                ]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-34',
                    'name': get_colon_term('spiral colon')[0],
                    'ontId': get_colon_term('spiral colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '35-36',
                    'name': get_colon_term('transverse colon')[0],
                    'ontId': get_colon_term('transverse colon')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '37-39',
                    'name': get_colon_term('descending colon')[0],
                    'ontId': get_colon_term('descending colon')[1]
                }]
        }),
    }

    @staticmethod
    def getName():
        return '3D Colon 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Cattle 1',
            'Human 1',
            'Human 2',
            'Human 3',
            'Mouse 1',
            'Mouse 2',
            'Pig 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 2']
        elif 'Human 3' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 3']
        elif 'Cattle 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Cattle 1']
        elif 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
        elif 'Mouse 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 2']
        elif 'Pig 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']

        if 'Human 3' in parameterSetName:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName='Human 2')
        elif 'Cattle' in parameterSetName:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName='Cattle 1')
        elif 'Mouse' in parameterSetName:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName='Mouse 1')
        elif 'Pig' in parameterSetName:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName='Pig 1')
        else:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName='Human 1')
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Segment profile': segmentProfileOption,
            'Number of segments': 30,
            'Start phase': 0.0,
            'Proximal tenia coli width': 10.0,
            'Proximal-transverse tenia coli width': 10.0,
            'Transverse-distal tenia coli width': 10.0,
            'Distal tenia coli width': 10.0,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }
        if 'Cattle 1' in parameterSetName:
            options['Number of segments'] = 40
            options['Proximal tenia coli width'] = 12.0
            options['Proximal-transverse tenia coli width'] = 6.0
            options['Transverse-distal tenia coli width'] = 3.0
            options['Distal tenia coli width'] = 3.0
        elif 'Human 3' in parameterSetName:
            options['Number of segments'] = 33
        elif 'Mouse' in parameterSetName:
            options['Number of segments'] = 10
            options['Proximal tenia coli width'] = 0.5
            options['Proximal-transverse tenia coli width'] = 0.5
            options['Transverse-distal tenia coli width'] = 0.5
            options['Distal tenia coli width'] = 0.5
        elif 'Pig 1' in parameterSetName:
            options['Number of segments'] = 120
            options['Proximal tenia coli width'] = 5.0
            options['Proximal-transverse tenia coli width'] = 4.0
            options['Transverse-distal tenia coli width'] = 3.0
            options['Distal tenia coli width'] = 1.5
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Segment profile',
            'Number of segments',
            'Start phase',
            'Proximal tenia coli width',
            'Proximal-transverse tenia coli width',
            'Transverse-distal tenia coli width',
            'Distal tenia coli width',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        if optionName == 'Segment profile':
            return [MeshType_3d_colonsegment1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Segment profile':
            if not parameterSetName:
                parameterSetName = scaffoldType.getParameterSetNames()[0]
            return ScaffoldPackage(scaffoldType, defaultParameterSetName=parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if not options['Segment profile'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Segment profile'):
            options['Segment profile'] = cls.getOptionScaffoldPackage('Segment profile', MeshType_3d_colonsegment1)
        for key in [
            'Number of segments',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Proximal tenia coli width',
            'Proximal-transverse tenia coli width',
            'Transverse-distal tenia coli width',
            'Distal tenia coli width']:
            if options[key] < 0.0:
                options[key] = 0.0

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        centralPath = options['Central path']
        segmentProfile = options['Segment profile']
        segmentCount = options['Number of segments']
        startPhase = options['Start phase'] % 360.0
        proximalTCWidth = options['Proximal tenia coli width']
        proximalTransverseTCWidth = options['Proximal-transverse tenia coli width']
        transverseDistalTCWidth = options['Transverse-distal tenia coli width']
        distalTCWidth = options['Distal tenia coli width']
        segmentSettings = segmentProfile.getScaffoldSettings()

        elementsCountAroundTC = segmentSettings['Number of elements around tenia coli']
        elementsCountAroundHaustrum = segmentSettings['Number of elements around haustrum']
        cornerInnerRadiusFactor = segmentSettings['Corner inner radius factor']
        haustrumInnerRadiusFactor = segmentSettings['Haustrum inner radius factor']
        segmentLengthEndDerivativeFactor = segmentSettings['Segment length end derivative factor']
        segmentLengthMidDerivativeFactor = segmentSettings['Segment length mid derivative factor']
        tcCount = segmentSettings['Number of tenia coli']
        tcThickness = segmentSettings['Tenia coli thickness']
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum) * tcCount

        elementsCountAlongSegment = segmentSettings['Number of elements along segment']
        elementsCountThroughWall = segmentSettings['Number of elements through wall']
        wallThickness = segmentSettings['Wall thickness']
        mucosaRelThickness = segmentSettings['Mucosa relative thickness']
        submucosaRelThickness = segmentSettings['Submucosa relative thickness']
        circularRelThickness = segmentSettings['Circular muscle layer relative thickness']
        longitudinalRelThickness = segmentSettings['Longitudinal muscle layer relative thickness']
        useCrossDerivatives = segmentSettings['Use cross derivatives']
        useCubicHermiteThroughWall = not (segmentSettings['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongSegment * segmentCount)

        # Colon coordinates
        lengthToDiameterRatio = 24
        wallThicknessToDiameterRatio = 0.1
        teniaColiThicknessToDiameterRatio = 0.25 * wallThicknessToDiameterRatio
        relativeThicknessListColonCoordinates = [1.0 / elementsCountThroughWall for n3 in range(elementsCountThroughWall)]

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Central path
        if tcCount == 1:
            colonTermsAlong = [None, 'right colon', 'transverse colon', 'left colon']
        elif tcCount == 2:
            colonTermsAlong = [None, 'spiral colon', 'transverse colon', 'descending colon']
        elif tcCount == 3:
            colonTermsAlong = [None, 'ascending colon', 'transverse colon', 'descending colon']

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        tmpFieldmodule = tmpRegion.getFieldmodule()
        cx, cd1, cd2, cd12 = get_nodeset_path_field_parameters(
            tmpFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES),
            tmpFieldmodule.findFieldByName('coordinates'),
            [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2])
        del tmpFieldmodule
        del tmpRegion

        # find arclength of colon
        length = 0.0
        elementsCountIn = len(cx) - 1
        sd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections=True,
                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
        for e in range(elementsCountIn):
            arcLength = interp.getCubicHermiteArcLength(cx[e], sd1[e], cx[e + 1], sd1[e + 1])
            # print(e+1, arcLength)
            length += arcLength
        segmentLength = length / segmentCount

        # Sample central path
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

        centralPathLength = arcLengthOfGroupsAlong[0]
        elementAlongLength = centralPathLength / elementsCountAlong

        elementsCountAlongGroups = []
        groupLength = 0.0
        e = 0
        elementsCount = 1
        length = elementAlongLength
        for i in range(1, len(colonTermsAlong)):
            groupLength += arcLengthOfGroupsAlong[i]
            if e == elementsCountAlong - 2:
                elementsCount += 1
                elementsCountAlongGroups.append(elementsCount)
            else:
                while length < groupLength:
                    elementsCount += 1
                    e += 1
                    length += elementAlongLength

                # check which end is grouplength closer to
                distToUpperEnd = abs(length - groupLength)
                distToLowerEnd = abs(groupLength - (length - elementsCountAlong))
                if distToLowerEnd < distToUpperEnd:
                    elementsCount -= 1
                    elementsCountAlongGroups.append(elementsCount)
                    e -= 1
                    length -= elementAlongLength
                else:
                    elementsCountAlongGroups.append(elementsCount)
            elementsCount = 0

        # Generate variation of radius & tc width along length
        lengthList = [0.0, arcLengthOfGroupsAlong[1], arcLengthOfGroupsAlong[1] + arcLengthOfGroupsAlong[2],
                      arcLengthOfGroupsAlong[0]]

        innerRadiusListCP = [vector.magnitude(c) for c in cd2]
        dInnerRadiusListCP = []
        for n in range(len(innerRadiusListCP) - 1):
            dInnerRadiusListCP.append(innerRadiusListCP[n + 1] - innerRadiusListCP[n])
        dInnerRadiusListCP.append(innerRadiusListCP[-1] - innerRadiusListCP[-2])
        innerRadiusAlongElementList, dInnerRadiusAlongElementList = interp.interpolateSampleCubicHermite(
            innerRadiusListCP, dInnerRadiusListCP, se, sxi, ssf)

        tcWidthList = [proximalTCWidth, proximalTransverseTCWidth, transverseDistalTCWidth, distalTCWidth]
        tcWidthAlongElementList, dTCWidthAlongElementList = interp.sampleParameterAlongLine(lengthList,
                                                                                            tcWidthList,
                                                                                            elementsCountAlong)

        # Account for reduced haustrum appearance in transverse and distal pig colon
        if tcCount == 2:
            haustrumInnerRadiusFactorList = [haustrumInnerRadiusFactor, haustrumInnerRadiusFactor * 0.75,
                                             haustrumInnerRadiusFactor * 0.5, haustrumInnerRadiusFactor * 0.2]
            haustrumInnerRadiusFactorAlongElementList = \
                interp.sampleParameterAlongLine(lengthList, haustrumInnerRadiusFactorList, elementsCountAlong)[0]
        else:
            haustrumInnerRadiusFactorAlongElementList = [haustrumInnerRadiusFactor] * (elementsCountAlong + 1)

        # Create annotation groups for colon sections
        colonGroup = AnnotationGroup(region, get_colon_term("colon"))

        if tcCount == 1:
            proximalGroup = AnnotationGroup(region, get_colon_term("right colon"))
            transverseGroup = AnnotationGroup(region, get_colon_term("transverse colon"))
            distalGroup = AnnotationGroup(region, get_colon_term("left colon"))
            annotationGroupAlong = [[colonGroup, proximalGroup],
                                    [colonGroup, transverseGroup],
                                    [colonGroup, distalGroup]]

        elif tcCount == 2:
            spiralGroup = AnnotationGroup(region, get_colon_term("spiral colon"))
            transverseGroup = AnnotationGroup(region, get_colon_term("transverse colon"))
            distalGroup = AnnotationGroup(region, get_colon_term("descending colon"))
            annotationGroupAlong = [[colonGroup, spiralGroup],
                                    [colonGroup, transverseGroup],
                                    [colonGroup, distalGroup]]

        elif tcCount == 3:
            ascendingGroup = AnnotationGroup(region, get_colon_term("ascending colon"))
            transverseGroup = AnnotationGroup(region, get_colon_term("transverse colon"))
            descendingGroup = AnnotationGroup(region, get_colon_term("descending colon"))
            annotationGroupAlong = [[colonGroup, ascendingGroup],
                                    [colonGroup, transverseGroup],
                                    [colonGroup, descendingGroup]]

        annotationGroupsAlong = []
        for i in range(len(elementsCountAlongGroups)):
            elementsCount = elementsCountAlongGroups[i]
            for n in range(elementsCount):
                annotationGroupsAlong.append(annotationGroupAlong[i])

        xExtrude = []
        d1Extrude = []
        d2Extrude = []
        d3UnitExtrude = []
        sxRefExtrudeList = []

        if elementsCountThroughWall == 1:
            relativeThicknessList = [1.0]
            annotationGroupsThroughWall = [[]]
        else:
            relativeThicknessList = [mucosaRelThickness, submucosaRelThickness,
                                     circularRelThickness, longitudinalRelThickness]
            mucosaGroup = AnnotationGroup(region, get_colon_term("colonic mucosa"))
            submucosaGroup = AnnotationGroup(region, get_colon_term("submucosa of colon"))
            circularMuscleGroup = AnnotationGroup(region, get_colon_term("circular muscle layer of colon"))
            longitudinalMuscleGroup = AnnotationGroup(region, get_colon_term("longitudinal muscle layer of colon"))
            annotationGroupsThroughWall = [[mucosaGroup], [submucosaGroup],
                                           [circularMuscleGroup], [longitudinalMuscleGroup]]

        # Create object
        colonSegmentTubeMeshInnerPoints = ColonSegmentTubeMeshInnerPoints(
            region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
            tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
            segmentLength, wallThickness, cornerInnerRadiusFactor, haustrumInnerRadiusFactorAlongElementList,
            innerRadiusAlongElementList, dInnerRadiusAlongElementList, tcWidthAlongElementList,
            startPhase)

        for nSegment in range(segmentCount):
            # Create inner points
            xInner, d1Inner, d2Inner, transitElementList, segmentAxis, annotationGroupsAround \
                = colonSegmentTubeMeshInnerPoints.getColonSegmentTubeMeshInnerPoints(nSegment)

            # Project reference point for warping onto central path
            start = nSegment * elementsCountAlongSegment
            end = (nSegment + 1) * elementsCountAlongSegment + 1
            sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
                tubemesh.getPlaneProjectionOnCentralPath(xInner, elementsCountAround, elementsCountAlongSegment,
                                                         segmentLength, sx[start:end], sd1[start:end], sd2[start:end],
                                                         sd12[start:end])

            # Warp segment points
            xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = tubemesh.warpSegmentPoints(
                xInner, d1Inner, d2Inner, segmentAxis, sxRefList, sd1RefList, sd2ProjectedListRef,
                elementsCountAround, elementsCountAlongSegment, zRefList)

            # Store points along length
            xExtrude += xWarpedList if nSegment == 0 else xWarpedList[elementsCountAround:]
            d1Extrude += d1WarpedList if nSegment == 0 else d1WarpedList[elementsCountAround:]
            d2Extrude += d2WarpedList if nSegment == 0 else d2WarpedList[elementsCountAround:]
            d3UnitExtrude += d3WarpedUnitList if nSegment == 0 else d3WarpedUnitList[elementsCountAround:]
            sxRefExtrudeList += sxRefList if nSegment == 0 else sxRefList[elementsCountAround:]

        contractedWallThicknessList = colonSegmentTubeMeshInnerPoints.getContractedWallThicknessList()

        # Create coordinates and derivatives
        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xExtrude, d1Extrude,
                                                                                        d2Extrude, d3UnitExtrude, contractedWallThicknessList, relativeThicknessList,
                                                                                        elementsCountAround, elementsCountAlong, elementsCountThroughWall, transitElementList)

        relaxedLengthList, xiList = colonSegmentTubeMeshInnerPoints.getRelaxedLengthAndXiList()

        closedProximalEnd = False

        if tcThickness > 0:
            tubeTCWidthList = colonSegmentTubeMeshInnerPoints.getTubeTCWidthList()
            xList, d1List, d2List, d3List, annotationArrayAround = getTeniaColi(
                region, xList, d1List, d2List, d3List, curvatureList, tcCount, elementsCountAroundTC,
                elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                tubeTCWidthList, tcThickness, sxRefExtrudeList, annotationGroupsAround,
                closedProximalEnd)

            # Create flat coordinates
            xFlat, d1Flat, d2Flat = createFlatCoordinatesTeniaColi(
                xiList, relaxedLengthList, length, wallThickness, relativeThicknessList, tcCount, tcThickness,
                elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong,
                elementsCountThroughWall, transitElementList, closedProximalEnd)

            # Create colon coordinates
            xColon, d1Colon, d2Colon = createColonCoordinatesTeniaColi(xiList, relativeThicknessListColonCoordinates,
                                                                       lengthToDiameterRatio,
                                                                       wallThicknessToDiameterRatio,
                                                                       teniaColiThicknessToDiameterRatio, tcCount,
                                                                       elementsCountAroundTC,
                                                                       elementsCountAroundHaustrum,
                                                                       elementsCountAlong, elementsCountThroughWall,
                                                                       transitElementList, closedProximalEnd)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = createNodesAndElementsTeniaColi(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xColon, d1Colon, d2Colon,
                "colon coordinates", elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong,
                elementsCountThroughWall, tcCount, annotationGroupsAround, annotationGroupsAlong,
                annotationGroupsThroughWall, firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall,
                useCrossDerivatives, closedProximalEnd)

        else:
            # Create flat coordinates
            xFlat, d1Flat, d2Flat = tubemesh.createFlatCoordinates(
                xiList, relaxedLengthList, length, wallThickness, relativeThicknessList, elementsCountAround,
                elementsCountAlong, elementsCountThroughWall, transitElementList)

            # Create colon coordinates
            xColon, d1Colon, d2Colon = tubemesh.createOrganCoordinates(xiList, relativeThicknessListColonCoordinates,
                                                                       lengthToDiameterRatio,
                                                                       wallThicknessToDiameterRatio,
                                                                       elementsCountAround,
                                                                       elementsCountAlong, elementsCountThroughWall,
                                                                       transitElementList)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xColon, d1Colon, d2Colon,
                "colon coordinates", elementsCountAround, elementsCountAlong, elementsCountThroughWall,
                annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
                firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
                closedProximalEnd)

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
        '''
        Add face annotation groups from the highest dimension mesh.
        Must have defined faces and added subelements for highest dimension groups.
        :param region: Zinc region containing model.
        :param options: Dict containing options. See getDefaultOptions().
        :param annotationGroups: List of annotation groups for top-level elements.
        New face annotation groups are appended to this list.
        '''
        # Create 2d surface mesh groups
        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)

        colonGroup = getAnnotationGroupForTerm(annotationGroups, get_colon_term("colon"))
        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        is_colon = colonGroup.getFieldElementGroup(mesh2d)
        is_serosa = fm.createFieldAnd(is_colon, is_exterior_face_xi3_1)
        serosa = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_colon_term("serosa of colon"))
        serosa.getMeshGroup(mesh2d).addElementsConditional(is_serosa)
        is_mucosaInnerSurface = fm.createFieldAnd(is_colon, is_exterior_face_xi3_0)
        mucosaInnerSurface = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_colon_term("luminal surface of the colonic mucosa"))
        mucosaInnerSurface.getMeshGroup(mesh2d).addElementsConditional(is_mucosaInnerSurface)
