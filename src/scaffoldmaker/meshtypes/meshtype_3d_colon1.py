"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy

from cmlibs.zinc.node import Node
from cmlibs.zinc.element import Element
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.colon_terms import get_colon_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1, ColonSegmentTubeMeshInnerPoints, \
    getTeniaColi, createFlatCoordinatesTeniaColi, createColonCoordinatesTeniaColi, createNodesAndElementsTeniaColi
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues


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
            'meshEdits': exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3 ], [
                [ [ -245.30, 444.60, -49.10 ], [ -267.70,  -53.10, -20.20 ], [   0.00,   0.00,  35.00 ], [   0.00,  0.00,  -4.41 ], [  -6.81,  34.33,  0.00 ], [   0.00,  0.00,  -4.41] ],
                [ [ -380.30, 484.80, -45.00 ], [   24.50,  102.70,  15.70 ], [   0.00,   0.00,  31.47 ], [   0.00,  0.00,  -2.65 ], [  30.61,  -7.30,  0.00 ], [   0.00,  0.00,  -2.65] ],
                [ [ -298.10, 510.40, -36.80 ], [   73.60,    9.90, -16.40 ], [   0.00,   0.00,  29.34 ], [ -12.52, -3.10, -14.34 ], [   3.91, -29.08,  0.00 ], [ -12.52, -3.10, -14.34] ],
                [ [ -213.10, 527.90, -22.50 ], [   -1.00,  -10.80, 125.60 ], [ -26.23,  -6.50,   1.62 ], [  -4.10, -6.11, -25.37 ], [   6.36, -26.22, -2.20 ], [  -4.10, -6.11, -25.37] ],
                [ [ -315.50, 570.20,  18.90 ], [ -107.90,    9.30,  21.90 ], [  -2.62, -12.12, -20.82 ], [  14.35,  4.82, -10.40 ], [   0.65, -20.97, 12.13 ], [  14.35,  4.82, -10.40] ],
                [ [ -417.40, 555.00,  14.60 ], [  -83.00,  -41.30,  -0.80 ], [   4.22,   1.17, -21.45 ], [   0.54,  4.88,   0.73 ], [   9.74, -19.59,  0.85 ], [   0.54,  4.88,   0.73] ],
                [ [ -497.30, 488.90,  13.60 ], [  -44.60,  -81.60,  10.00 ], [  -1.56,  -2.39, -19.34 ], [  -2.12, -0.60,   2.13 ], [  17.14,  -9.40, -0.22 ], [  -2.12, -0.60,   2.13] ],
                [ [ -527.00, 392.50,   2.70 ], [   47.40,  -82.00,  -7.90 ], [   0.00,   0.00, -17.20 ], [   0.68,  1.04,   1.96 ], [  14.89,   8.61,  0.00 ], [   0.68,  1.04,   1.96] ],
                [ [ -461.20, 345.90,  -0.80 ], [   56.90,  -44.50,   2.40 ], [   0.00,   0.00, -15.38 ], [   0.00,  0.00,   1.09 ], [   9.47,  12.11,  0.00 ], [   0.00,  0.00,   1.09] ],
                [ [ -415.60, 293.80,   3.90 ], [   93.20,  -62.60,   3.10 ], [   0.00,   0.00, -14.91 ], [   0.00,  0.00,   0.43 ], [   8.31,  12.38,  0.00 ], [   0.00,  0.00,   0.43] ],
                [ [ -232.20, 264.90,   0.20 ], [  140.10,   58.20,  -1.00 ], [   0.00,   0.00, -14.58 ], [   0.00,  0.00,   0.25 ], [  -5.59,  13.46,  0.00 ], [   0.00,  0.00,   0.25] ],
                [ [ -168.40, 357.20,   1.30 ], [   10.10,   78.60,  -3.20 ], [   0.00,   0.00, -14.38 ], [   0.00,  0.00,   0.15 ], [ -14.26,   1.83,  0.00 ], [   0.00,  0.00,   0.15] ],
                [ [ -185.30, 419.10,  -0.70 ], [  -45.10,   57.10,  -0.90 ], [   0.00,   0.00, -14.27 ], [   0.00,  0.00,   0.13 ], [ -11.20,  -8.84,  0.00 ], [   0.00,  0.00,   0.13] ],
                [ [ -253.20, 466.70,  -0.30 ], [  -63.40,   24.70,   0.20 ], [   0.00,   0.00, -14.13 ], [   0.00,  0.00,   0.13 ], [  -5.13, -13.17,  0.00 ], [   0.00,  0.00,   0.13] ],
                [ [ -323.80, 482.50,   0.10 ], [  -68.20,    2.90,  -1.20 ], [   0.00,   0.00, -14.00 ], [   0.00,  0.00,   0.12 ], [  -0.59, -13.99,  0.00 ], [   0.00,  0.00,   0.12] ],
                [ [ -387.50, 485.40,  -0.20 ], [  -44.20,  -17.10,  -1.00 ], [   0.00,   0.00, -13.89 ], [   0.00,  0.00,   0.12 ], [   5.01, -12.95,  0.00 ], [   0.00,  0.00,   0.12] ],
                [ [ -435.60, 433.50,   3.30 ], [    3.40, -109.50,   1.40 ], [   0.00,   0.00, -13.76 ], [   0.00,  0.00,   0.14 ], [  13.75,   0.43,  0.00 ], [   0.00,  0.00,   0.14] ],
                [ [ -370.60, 376.30,  -1.10 ], [   66.90,  -29.20,  -0.90 ], [   0.00,   0.00, -13.60 ], [   0.00,  0.00,   0.12 ], [   5.44,  12.46,  0.00 ], [   0.00,  0.00,   0.12] ],
                [ [ -313.00, 357.90,  -0.10 ], [   40.00,  -33.50,   9.60 ], [   0.00,   0.00, -13.50 ], [   0.00,  0.00,   0.10 ], [   8.67,  10.35,  0.00 ], [   0.00,  0.00,   0.10] ],
                [ [ -259.20, 340.70,   2.10 ], [   48.90,    6.40,   1.40 ], [   0.00,   0.00, -13.40 ], [   0.00,  0.00,   0.09 ], [  -1.74,  13.29,  0.00 ], [   0.00,  0.00,   0.09] ],
                [ [ -246.50, 380.30,  -0.80 ], [  -29.70,   33.60,  -0.70 ], [   0.00,   0.00, -13.31 ], [   0.00,  0.00,   0.09 ], [  -9.97,  -8.82,  0.00 ], [   0.00,  0.00,   0.09] ],
                [ [ -297.30, 387.10,   0.60 ], [  -59.70,   12.60,  -0.00 ], [   0.00,   0.00, -13.22 ], [   0.00,  0.00,   0.09 ], [  -2.73, -12.94,  0.00 ], [   0.00,  0.00,   0.09] ],
                [ [ -340.20, 415.60,  -1.00 ], [  -86.20,   28.90,  -2.90 ], [   0.00,   0.00, -13.13 ], [   0.00,  0.00,   0.10 ], [  -4.17, -12.45,  0.00 ], [   0.00,  0.00,   0.10] ],
                [ [ -398.30, 443.10,  -0.10 ], [   10.60,   82.10,  -2.60 ], [   0.00,   0.00, -13.01 ], [   0.00,  0.00,   0.12 ], [ -12.90,   1.67,  0.00 ], [   0.00,  0.00,   0.12] ],
                [ [ -329.80, 449.10,  -2.10 ], [   53.20,   14.00,  -0.50 ], [   0.00,   0.00, -12.88 ], [   0.00,  0.00,   0.14 ], [  -3.28,  12.46,  0.00 ], [   0.00,  0.00,   0.14] ],
                [ [ -251.30, 425.90,  -0.30 ], [   43.90,  -19.30,   0.00 ], [   0.00,   0.00, -12.74 ], [   0.00,  0.00,   0.11 ], [   5.13,  11.66,  0.00 ], [   0.00,  0.00,   0.11] ],
                [ [ -209.10, 390.60,   0.00 ], [   26.00,  -38.80,   0.90 ], [   0.00,   0.00, -12.65 ], [   0.00,  0.00,   0.08 ], [  10.51,   7.04,  0.00 ], [   0.00,  0.00,   0.08] ],
                [ [ -207.80, 350.80,   1.40 ], [   -9.40,  -43.60,   1.80 ], [   0.00,   0.00, -12.57 ], [   0.00,  0.00,   0.09 ], [  12.29,  -2.65,  0.00 ], [   0.00,  0.00,   0.09] ],
                [ [ -245.80, 299.40,   7.60 ], [  -70.30,  -36.00,   1.40 ], [   0.00,   0.00, -12.46 ], [   0.00,  0.00,   0.14 ], [   5.68, -11.09,  0.00 ], [   0.00,  0.00,   0.14] ],
                [ [ -345.30, 304.10,   3.10 ], [ -100.20,   27.90,  -1.90 ], [   0.00,   0.00, -12.29 ], [   0.00,  0.00,   0.17 ], [  -3.30, -11.84,  0.00 ], [   0.00,  0.00,   0.17] ],
                [ [ -418.40, 361.10,  -0.20 ], [  -57.80,   55.80,  -1.70 ], [   0.00,   0.00, -12.13 ], [   0.00,  0.00,   0.15 ], [  -8.42,  -8.73,  0.00 ], [   0.00,  0.00,   0.15] ],
                [ [ -479.20, 415.60,   2.20 ], [   -8.80,   73.10,  -1.60 ], [   0.00,   0.00, -11.98 ], [   0.00,  0.00,   0.15 ], [ -11.89,  -1.43,  0.00 ], [   0.00,  0.00,   0.15] ],
                [ [ -439.60, 495.70,  -2.10 ], [   61.10,   57.10,  -1.30 ], [   0.00,   0.00, -11.82 ], [   0.00,  0.00,   0.15 ], [  -8.07,   8.64,  0.00 ], [   0.00,  0.00,   0.15] ],
                [ [ -361.60, 522.60,  -3.00 ], [   78.60,    9.90,   0.20 ], [   0.00,   0.00, -11.68 ], [   0.00,  0.00,   0.15 ], [  -1.46,  11.59,  0.00 ], [   0.00,  0.00,   0.15] ],
                [ [ -270.10, 506.50,  -3.80 ], [  103.60,  -33.30,   1.00 ], [   0.00,   0.00, -11.52 ], [   0.00,  0.00,   0.19 ], [   3.53,  10.97,  0.00 ], [   0.00,  0.00,   0.19] ],
                [ [ -148.90, 441.40,  -2.10 ], [   79.70,  -91.50,   2.80 ], [   0.00,   0.00, -11.28 ], [   0.00,  0.00,   0.23 ], [   8.51,   7.41,  0.00 ], [   0.00,  0.00,   0.23] ],
                [ [ -130.90, 313.30,   4.00 ], [   -4.00, -107.20,   3.10 ], [   0.00,   0.00, -11.05 ], [   0.00,  0.00,   0.18 ], [  11.04,  -0.41,  0.00 ], [   0.00,  0.00,   0.18] ],
                [ [ -183.90, 251.00,   3.80 ], [  -65.50,  -60.20,   3.60 ], [   0.00,   0.00, -10.90 ], [   0.00,  0.00,   0.16 ], [   7.38,  -8.03,  0.00 ], [   0.00,  0.00,   0.16] ],
                [ [ -280.30, 213.00,   3.40 ], [ -165.10,  -18.60,   0.10 ], [   0.00,   0.00, -10.72 ], [   0.00,  0.00,   0.20 ], [   1.20, -10.65,  0.00 ], [   0.00,  0.00,   0.20] ],
                [ [ -400.80, 247.50,   6.80 ], [ -127.10,   36.80,   1.30 ], [   0.00,   0.00, -10.51 ], [   0.00,  0.00,   0.23 ], [  -2.92, -10.10,  0.00 ], [   0.00,  0.00,   0.23] ],
                [ [ -530.50, 290.70,   5.20 ], [  -89.00,   86.50,   0.30 ], [   0.00,   0.00, -10.27 ], [   0.00,  0.00,   0.21 ], [  -7.16,  -7.36,  0.00 ], [   0.00,  0.00,   0.21] ],
                [ [ -568.80, 392.30,   6.90 ], [  -77.40,   67.70,  -5.50 ], [   0.00,   0.00, -10.08 ], [   0.00,  0.00,   0.23 ], [  -6.64,  -7.59,  0.00 ], [   0.00,  0.00,   0.23] ],
                [ [ -511.20, 535.10,   2.50 ], [   86.20,  111.40,  -1.00 ], [   0.00,   0.00,  -9.80 ], [   0.00,  0.00,   0.25 ], [  -7.75,   6.00,  0.00 ], [   0.00,  0.00,   0.25] ],
                [ [ -405.00, 601.70,   6.40 ], [  143.60,   52.20,   2.60 ], [   0.00,   0.00,  -9.58 ], [   1.49, -1.44,   0.82 ], [  -3.27,   9.00,  0.00 ], [   1.49, -1.44,   0.82] ],
                [ [ -238.80, 615.90,  16.60 ], [   63.30,   -9.10,  19.10 ], [   3.45,  -3.33,  -7.95 ], [  -5.13,  0.21,   1.79 ], [   2.06,   8.63, -2.72 ], [  -5.13,  0.21,   1.79] ],
                [ [ -146.20, 605.90,  36.50 ], [   49.30,   -9.90, -50.60 ], [  -6.72,  -1.05,  -6.07 ], [  -3.12,  3.94,   7.67 ], [   0.10,   8.96, -1.66 ], [  -3.12,  3.94,   7.67] ],
                [ [ -218.40, 585.30,  -2.00 ], [ -124.00,    0.40, -37.50 ], [  -3.04,   4.50,   7.19 ], [   2.76,  2.30,   8.98 ], [   1.33,   7.80, -4.32 ], [   2.76,  2.30,   8.98] ],
                [ [ -376.30, 579.60, -40.80 ], [ -189.20,  -50.70,  -8.80 ], [  -1.85,   1.22,   8.78 ], [   4.16,  0.49,  -2.60 ], [  -2.23,   8.62, -1.67 ], [   4.16,  0.49,  -2.60] ],
                [ [ -557.90, 493.90, -24.90 ], [  -30.30,   24.10, 152.80 ], [   6.19,   6.65,   0.69 ], [   1.49, -2.48,  -8.73 ], [  -6.35,   6.14, -2.23 ], [   1.49, -2.48,  -8.73] ],
                [ [ -484.80, 594.40,   0.70 ], [  132.70,   97.00,   3.50 ], [   3.25,  -1.18,  -8.47 ], [  -2.78, -5.06,  -5.12 ], [  -5.07,   7.03, -2.92 ], [  -2.78, -5.06,  -5.12] ],
                [ [ -318.10, 641.90,  -8.50 ], [  166.70,   17.60,   5.50 ], [   0.67,  -2.88,  -8.72 ], [  -2.15, -0.82,  -0.11 ], [  -0.82,   8.69, -2.93 ], [  -2.15, -0.82,  -0.11] ],
                [ [ -158.30, 634.70,  -1.90 ], [  176.50,  -14.00,  10.80 ], [  -1.08,  -2.89,  -8.71 ], [  -0.81,  0.27,  -0.12 ], [   0.87,   8.70, -2.99 ], [  -0.81,  0.27,  -0.12] ],
                [ [   32.70, 611.70,  13.60 ], [  205.50,  -32.20,  20.00 ], [  -0.76,  -2.28,  -8.98 ], [   1.46,  0.96,  -0.42 ], [   1.62,   8.84, -2.38 ], [   1.46,  0.96,  -0.42] ] ] ),
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
            'meshEdits': exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [   0.00,   0.00, 0.00 ], [ -50.70,  178.20, 0.00 ], [ -37.97,  -9.49, -18.98 ], [ -6.86, -11.39,  -2.36 ], [ -18.61,  -3.98, 39.12 ], [ -14.00,  -1.00, -12.00] ],
                [ [ -47.40, 188.60, 0.00 ], [ -19.30,  177.10, 0.00 ], [ -35.79,  -6.51, -13.01 ], [ 11.23,  17.36,  14.31 ], [ -12.66,  -3.99, 36.28 ], [  -4.00,  19.00,  22.00] ],
                [ [  -4.40, 396.50, 0.00 ], [ 206.00,   40.10, 0.00 ], [ -13.89,  27.78,  11.11 ], [ 13.54,  -1.87,  21.51 ], [  -6.05, -12.50, 29.93 ], [  -6.00,   0.00,  51.00] ],
                [ [ 130.00, 384.10, 0.00 ], [ 130.80,  -40.50, 0.00 ], [  -5.35,   4.28,  31.06 ], [  5.83,  -8.41,   8.86 ], [ -15.28, -27.78,  2.51 ], [   0.00,   1.00,  24.00] ],
                [ [ 279.40, 383.00, 0.00 ], [ 118.00,   48.70, 0.00 ], [  -2.51,  12.57,  27.65 ], [  9.30,   9.73, -10.83 ], [  12.83, -24.62, 12.58 ], [   5.00,  25.00, -20.00] ],
                [ [ 443.90, 390.80, 0.00 ], [ 111.30,  -97.00, 0.00 ], [  14.07,  23.92,   8.44 ], [ 12.36,  -3.74, -23.50 ], [  -9.40,  -3.01, 27.28 ], [   1.00,  -6.00, -35.00] ],
                [ [ 475.20, 168.00, 0.00 ], [  -0.80, -112.40, 0.00 ], [  20.78,   0.00, -20.78 ], [ -2.41, -19.32, -15.36 ], [  20.78,  -0.00, 20.78 ], [  15.00,  -1.00, -10.00] ],
                [ [ 432.60, -32.30, 0.00 ], [ -90.50,  -59.00, 0.00 ], [  10.09, -15.13, -23.54 ], [ -9.58,  -7.07,  -2.37 ], [  13.01, -19.62, 18.19 ], [   8.00, -11.00, -13.00] ],
                [ [ 272.40,   7.50, 0.00 ], [ -79.00,   47.40, 0.00 ], [   1.42, -15.65, -25.60 ], [ -7.76,   6.05,  -1.75 ], [  -5.22, -26.72, 12.68 ], [   4.00, -12.00, -12.00] ] ] ),
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
            'meshEdits': exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [   0.00,    0.00,   0.00 ], [  -56.81,  105.14,  -38.05 ], [ -37.97,  -9.49, -18.98 ], [ -3.67, -11.77,  -1.56 ], [ -21.20,   5.92, 37.52 ], [ -14.00,  -1.00, -12.00 ] ],
                [ [ -34.50,  114.00, -18.10 ], [   -9.51,  117.91,    3.65 ], [ -33.74,  -6.13, -12.27 ], [ 12.13,  18.49,  14.99 ], [ -12.58,  -4.62, 33.86 ], [  -4.00,  19.00,  22.00 ] ],
                [ [ -19.10,  218.50,   5.50 ], [   79.23,   66.40,   77.49 ], [ -13.72,  27.44,  10.98 ], [ 14.61,   7.18,  21.72 ], [ -22.92, -17.43, 15.26 ], [  -6.00,   0.00,  51.00 ] ],
                [ [  82.50,  189.10,  94.20 ], [  140.70,   -1.06,   48.14 ], [  -5.34,   4.27,  30.95 ], [  5.75,  -8.29,   9.07 ], [  11.54, -25.97, 14.03 ], [   0.00,   1.00,  24.00 ] ],
                [ [ 226.60,  218.70,  85.70 ], [  164.08,  101.90,  -75.52 ], [  -2.53,  12.64,  27.81 ], [  8.00,   9.47,  -9.25 ], [  19.49, -20.43, 11.94 ], [   5.00,  25.00, -20.00 ] ],
                [ [ 325.50,  381.70, -57.90 ], [  187.36, -116.61, -173.53 ], [  14.07,  23.92,   8.44 ], [ 12.05,  -5.07, -23.99 ], [   5.51, -10.97, 26.28 ], [   1.00,  -6.00, -35.00 ] ],
                [ [ 354.00,  105.30, -24.40 ], [  -20.59, -269.54,   30.48 ], [  20.87,   0.00, -20.87 ], [ -2.92, -19.10, -14.61 ], [  20.81,   5.79, 20.11 ], [  15.00,  -1.00, -10.00 ] ],
                [ [ 296.50, -121.20,  -0.60 ], [ -170.98, -102.19,  -18.39 ], [  10.15, -15.22, -23.67 ], [ -9.48,  -6.06,  -2.33 ], [  13.09, -19.73, 18.29 ], [   8.00, -11.00, -13.00 ] ],
                [ [ 169.80,  -73.40, -33.50 ], [  -42.47,  101.91,  -24.43 ], [   1.43, -15.71, -25.71 ], [ -7.96,   5.07,  -1.75 ], [ -16.72, -21.84, 12.39 ], [   4.00, -12.00, -12.00 ] ] ] ),
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
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,  Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [   0.00,   0.00,  0.00 ], [  6.00, 12.00,  -2.00 ], [ 0.67,  0.33,  0.67 ], [ -0.03,  0.95, -0.68 ], [  0.86,  0.04, -0.51 ], [ 6.00, 0.00, 3.00 ] ],
                [ [  -2.00,  11.00, -3.00 ], [ -8.00,  4.00,   9.00 ], [ 0.64,  0.64,  0.32 ], [ -0.03, -0.34, -0.02 ], [ -0.37,  0.66, -0.58 ], [ 0.00, 1.00, 2.00 ] ],
                [ [  -3.00,   2.00,  3.00 ], [ -4.00, -8.00,   0.00 ], [ 0.61, -0.30,  0.61 ], [ -0.15, -0.65,  0.01 ], [ -0.55,  0.27,  0.68 ], [ 1.00, 0.00, 2.00 ] ],
                [ [ -11.00,  -3.00, -4.00 ], [ -8.00, -3.00,  -7.00 ], [ 0.33, -0.67,  0.33 ], [ -0.18, -0.18, -0.31 ], [ -0.32,  0.10,  0.75 ], [ 0.00, 0.00, 0.50 ] ],
                [ [ -16.00,  -4.00,  0.00 ], [  4.00, -3.00,  14.00 ], [ 0.23, -0.70,  0.00 ], [ -0.16, -0.02, -0.19 ], [  0.71,  0.18,  0.05 ], [ 0.00, 0.00, 0.50 ] ],
                [ [  -7.00,  -8.00,  0.00 ], [  5.00, -1.00, -14.00 ], [ 0.00, -0.70,  0.00 ], [  0.03,  0.04, -0.12 ], [ -0.64, -0.00, -0.28 ], [ 0.00, 0.00, 0.50 ] ],
                [ [  -1.00,  -6.00, -1.00 ], [  2.00, -2.00,   9.00 ], [ 0.21, -0.63, -0.21 ], [  0.12,  0.11, -0.23 ], [  0.64,  0.25, -0.11 ], [ 0.00, 0.00, 0.50 ] ],
                [ [  -2.00, -14.00,  5.00 ], [ -2.00, -4.00,   2.00 ], [ 0.23, -0.47, -0.47 ], [ -0.08,  0.22, -0.28 ], [  0.53, -0.17,  0.42 ], [ 0.00, 0.00, 0.50 ] ] ] ),
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
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,  Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [   0.00,  0.00,   0.00 ], [ 0.00,  0.00,  13.00 ], [ 0.00, -1.00, 0.00 ], [ 0.00,  0.03, 0.00 ], [  1.00, -0.00,  0.00 ], [ 0.00, 0.00, 0.50 ] ],
                [ [   0.00,  0.00,  13.00 ], [ 0.00,  2.00,  28.00 ], [ 0.00, -0.96, 0.00 ], [ 0.00,  0.05, 0.00 ], [  0.95, -0.00, -0.07 ], [ 0.00, 0.00, 0.50 ] ],
                [ [ -14.00, -2.00,  13.00 ], [ 0.00, -3.00, -19.00 ], [ 0.00, -0.88, 0.00 ], [ 0.00,  0.13, 0.00 ], [ -0.87, -0.02, -0.14 ], [ 0.00, 0.00, 0.50 ] ],
                [ [ -14.00, -1.00, -10.00 ], [ 1.00,  1.00, -17.00 ], [ 0.00, -0.70, 0.00 ], [ 0.00,  0.08, 0.00 ], [ -0.70, -0.00, -0.00 ], [ 0.00, 0.00, 0.50 ] ],
                [ [ -14.00,  0.00, -28.00 ], [ 0.00,  0.00, -11.00 ], [ 0.00, -0.70, 0.00 ], [ 0.00, -0.08, 0.00 ], [ -0.70, -0.00,  0.00 ], [ 0.00, 0.00, 0.50 ] ] ] ),
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
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,  Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [   -7.20,   83.30,  -20.70 ], [  -65.20,   -8.10,   7.60 ], [  0.00,  0.00, 36.00 ], [  0.00,  0.00, -0.50 ], [  -8.49,  34.98,  0.48 ], [  0.00,  0.00, -0.50 ] ],
                [ [  -68.50,   52.80,   -9.60 ], [  -40.10,  -36.10,  10.70 ], [  0.00,  0.00, 35.44 ], [  0.00,  0.00, -0.63 ], [ -28.28,  21.31,  1.34 ], [  0.00,  0.00, -0.63 ] ],
                [ [  -97.40,  -26.30,    5.70 ], [   18.00,  -93.20,  13.70 ], [  0.00,  0.00, 34.73 ], [  0.00,  0.00, -0.67 ], [ -32.82, -11.34,  0.71 ], [  0.00,  0.00, -0.67 ] ],
                [ [  -56.80,  -90.50,   14.10 ], [   65.50,  -41.40,   7.30 ], [  0.00,  0.00, 34.10 ], [  0.00,  0.00, -0.74 ], [ -15.44, -30.39,  0.30 ], [  0.00,  0.00, -0.74 ] ],
                [ [   48.90, -100.80,   24.00 ], [  112.20,   40.10,  19.00 ], [  0.00,  0.00, 33.22 ], [  0.00,  0.00, -0.91 ], [  15.90, -29.15,  0.82 ], [  0.00,  0.00, -0.91 ] ],
                [ [  114.80,  -12.60,   38.70 ], [    8.20,   96.10,  14.20 ], [  0.00,  0.00, 32.28 ], [  0.00,  0.00, -0.93 ], [  32.22,   1.92,  0.68 ], [  0.00,  0.00, -0.93 ] ],
                [ [   60.30,   83.50,   43.70 ], [ -108.70,   54.10,  22.40 ], [  0.00,  0.00, 31.35 ], [  0.00,  0.00, -0.89 ], [   8.73,  30.09,  1.03 ], [  0.00,  0.00, -0.89 ] ],
                [ [  -41.20,   90.70,   56.30 ], [  -89.00,  -32.40,  14.40 ], [  0.00,  0.00, 30.50 ], [  0.00,  0.00, -0.93 ], [ -14.57,  26.79,  0.69 ], [  0.00,  0.00, -0.93 ] ],
                [ [ -107.90,   -9.70,   76.60 ], [   11.10,  -94.40,  11.30 ], [  0.00,  0.00, 29.47 ], [  0.00,  0.00, -0.91 ], [ -28.66,  -6.85,  0.41 ], [  0.00,  0.00, -0.91 ] ],
                [ [  -57.30,  -91.90,   81.30 ], [   71.20,  -31.20,   5.70 ], [  0.00,  0.00, 28.65 ], [  0.00,  0.00, -0.86 ], [  -9.56, -27.01,  0.15 ], [  0.00,  0.00, -0.86 ] ],
                [ [   51.20,  -89.40,   97.20 ], [   99.10,   55.40,  12.90 ], [  0.00,  0.00, 27.75 ], [  0.00,  0.00, -0.90 ], [  16.17, -22.55,  0.35 ], [  0.00,  0.00, -0.90 ] ],
                [ [   91.60,    9.30,  103.60 ], [    4.70,   51.20,   3.40 ], [  0.00,  0.00, 26.86 ], [  0.00,  0.00, -0.90 ], [  26.85,  -0.69,  0.12 ], [  0.00,  0.00, -0.90 ] ],
                [ [   61.60,  111.80,  109.60 ], [  -85.20,   46.10,   2.60 ], [  0.00,  0.00, 25.95 ], [  0.00,  0.00, -0.95 ], [  11.73,  23.14,  0.02 ], [  0.00,  0.00, -0.95 ] ],
                [ [  -54.60,   91.90,  129.40 ], [  -92.70,  -55.00,  14.50 ], [  0.00,  0.00, 24.95 ], [  0.00,  0.00, -0.94 ], [ -15.45,  19.58,  0.44 ], [  0.00,  0.00, -0.94 ] ],
                [ [ -109.00,    5.60,  156.90 ], [   23.60, -108.20,  27.70 ], [  0.00,  0.00, 24.05 ], [  0.00,  0.00, -0.80 ], [ -21.59, -10.51,  1.42 ], [  0.00,  0.00, -0.80 ] ],
                [ [  -59.10,  -62.50,  170.80 ], [   74.00,  -20.10,  14.40 ], [  0.00,  0.00, 23.34 ], [  0.00,  0.00, -0.70 ], [  -1.93, -23.24,  0.79 ], [  0.00,  0.00, -0.70 ] ],
                [ [   23.50,  -53.20,  179.70 ], [   84.60,   47.00,   6.90 ], [  0.00,  0.00, 22.65 ], [  0.00,  0.00, -0.73 ], [  12.38, -18.97,  0.11 ], [  0.00,  0.00, -0.73 ] ],
                [ [   62.30,   30.10,  187.50 ], [  -12.80,   58.00,   0.80 ], [  0.00,  0.00, 21.88 ], [  0.00,  0.00, -0.55 ], [  21.30,   5.00,  0.00 ], [  0.00,  0.00, -0.55 ] ],
                [ [   22.40,   45.20,  181.10 ], [  -23.60,  -34.50,  -7.40 ], [  0.00,  0.00, 21.44 ], [  0.00,  0.00, -0.41 ], [ -15.35,  14.96,  0.65 ], [  0.00,  0.00, -0.41 ] ],
                [ [   -1.90,    4.90,  180.50 ], [  -41.30,  -30.90,   7.50 ], [  0.00,  0.00, 21.06 ], [  0.00,  0.00, -0.39 ], [ -14.89,  14.89,  0.44 ], [  0.00,  0.00, -0.39 ] ],
                [ [  -45.10,  -12.60,  194.40 ], [  -40.50,   -4.60,   6.90 ], [  0.00,  0.00, 20.67 ], [  0.00,  0.00, -0.30 ], [  -5.68,  19.87,  0.58 ], [  0.00,  0.00, -0.30 ] ],
                [ [  -71.70,   -2.20,  197.20 ], [  -25.20,   35.80,  -6.80 ], [  0.00,  0.00, 20.42 ], [  0.00,  0.00, -0.30 ], [  18.28,   9.08,  0.48 ], [  0.00,  0.00, -0.30 ] ],
                [ [  -65.80,   42.10,  182.30 ], [   26.60,   37.60, -15.60 ], [  0.00,  0.00, 20.03 ], [  0.00,  0.00, -0.46 ], [  11.97, -15.92,  2.06 ], [  0.00,  0.00, -0.46 ] ],
                [ [  -14.10,   81.20,  163.50 ], [   41.00,   10.30,  -9.50 ], [  0.00,  0.00, 19.48 ], [  0.00,  0.00, -0.59 ], [   0.59, -19.44,  0.94 ], [  0.00,  0.00, -0.59 ] ],
                [ [   61.70,   86.10,  156.40 ], [   77.90,  -40.70,   8.90 ], [  0.00,  0.00, 18.85 ], [  0.00,  0.00, -0.62 ], [  -7.01, -17.50,  0.19 ], [  0.00,  0.00, -0.62 ] ],
                [ [   92.90,   20.50,  150.30 ], [    0.00,  -73.30,  -5.20 ], [  0.00,  0.00, 18.23 ], [  0.00,  0.00, -0.69 ], [ -18.19,   1.29,  0.09 ], [  0.00,  0.00, -0.69 ] ],
                [ [   48.90,  -65.00,  142.80 ], [  -82.80,  -80.00,  -1.90 ], [  0.00,  0.00, 17.44 ], [  0.00,  0.00, -0.84 ], [ -11.91,  12.74,  0.00 ], [  0.00,  0.00, -0.84 ] ],
                [ [  -54.30,  -90.80,  134.00 ], [  -60.10,   26.40,  -8.20 ], [  0.00,  0.00, 16.54 ], [  0.00,  0.00, -0.92 ], [   8.46,  14.21,  0.25 ], [  0.00,  0.00, -0.92 ] ],
                [ [  -89.90,   11.20,  115.00 ], [   34.90,  125.10, -27.90 ], [  0.00,  0.00, 15.59 ], [  0.00,  0.00, -0.88 ], [  13.83,  -7.18,  0.69 ], [  0.00,  0.00, -0.88 ] ],
                [ [  -17.40,   74.20,   91.10 ], [   78.80,   19.10, -15.40 ], [  0.00,  0.00, 14.77 ], [  0.00,  0.00, -0.67 ], [   0.79, -14.74,  0.51 ], [  0.00,  0.00, -0.67 ] ],
                [ [   43.40,   50.20,   73.30 ], [   30.20,  -36.00,  -9.90 ], [  0.00,  0.00, 14.20 ], [  0.00,  0.00, -0.52 ], [ -12.48,  -6.73,  0.60 ], [  0.00,  0.00, -0.52 ] ],
                [ [   62.40,   -5.10,   63.50 ], [   10.90,  -54.20,  -2.70 ], [  0.00,  0.00, 13.72 ], [  0.00,  0.00, -0.48 ], [ -13.56,  -2.05,  0.03 ], [  0.00,  0.00, -0.48 ] ],
                [ [   32.70,  -51.70,   56.10 ], [  -38.60,  -29.80,  -8.10 ], [  0.00,  0.00, 13.24 ], [  1.37,  2.33, -1.30 ], [  -6.29,  11.65,  0.36 ], [  1.59,  2.70, -1.41 ] ],
                [ [  -38.10,  -28.60,   46.80 ], [  -62.50,   82.60, -19.20 ], [  3.28,  5.57, 10.81 ], [  1.59, -2.89, -2.75 ], [   9.22,   6.17, -5.96 ], [  1.85, -3.40, -3.08 ] ],
                [ [    5.70,   40.40,   22.40 ], [  144.80,   18.60, -20.50 ], [  2.68, -8.27,  7.65 ], [ -4.75, -6.37, -4.66 ], [   0.33,  -7.82, -8.54 ], [ -5.62, -7.55, -5.37 ] ],
                [ [   53.00,  -14.70,   -4.10 ], [   -6.00,  -25.70, -46.70 ], [ -5.67, -8.15,  1.68 ], [ -3.78,  1.55, -1.86 ], [  -8.24,   3.65, -4.49 ], [ -4.60,  1.69, -2.12 ] ],
                [ [   24.80,   -0.40,  -48.80 ], [  -13.40,   23.90, -30.60 ], [ -6.57, -5.70,  2.41 ], [  4.46,  2.06,  2.31 ], [  -1.87,   3.95,  7.90 ], [  5.41,  2.39,  2.88 ] ],
                [ [  -20.90,   15.30,  -77.90 ], [  -51.20,  -30.60,  21.10 ], [  3.74, -4.06,  6.43 ], [ -1.12,  2.71, -1.85 ], [  -2.63,   5.93,  5.46 ], [ -1.44,  3.31, -2.27 ] ],
                [ [  -47.60,   33.90, -112.20 ], [   32.60,   30.70, -27.80 ], [ -8.36, -0.32, -1.04 ], [ -7.25,  2.76, -6.03 ], [  -3.75,   2.25,  7.20 ], [ -8.97,  3.42, -7.45 ] ],
                [ [   19.60,   96.00, -167.50 ], [   19.90,   19.10, -18.40 ], [ -6.99,  0.72, -4.50 ], [  9.97, -0.68, -0.89 ], [  -2.97,   5.84,  5.17 ], [ 12.30, -0.82, -1.11 ] ] ] ),
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
        'Pig 2': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'Length': 90.0,
                'Number of elements': 3
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                [ [  0.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [ 30.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [ 60.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [ 90.0, 0.0, 0.0 ], [ 30.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ] ] )
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
            'Mouse 1',
            'Mouse 2',
            'Pig 1',
            'Pig 2']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 2']
        elif 'Cattle 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Cattle 1']
        elif 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
        elif 'Mouse 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 2']
        elif 'Pig 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
        elif 'Pig 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 2']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
        if 'Cattle' in parameterSetName:
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
        elif 'Pig 2' in parameterSetName:
            options['Number of segments'] = 3
            options['Proximal tenia coli width'] = 5.0
            options['Proximal-transverse tenia coli width'] = 5.0
            options['Transverse-distal tenia coli width'] = 5.0
            options['Distal tenia coli width'] = 5.0
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
        :return: annotationGroups
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

        arcLengthOfGroupsAlong = []
        for i in range(len(colonTermsAlong)):
            tmpRegion = region.createRegion()
            centralPath.generate(tmpRegion)
            cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                            Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                                groupName=colonTermsAlong[i])
            arcLength = 0.0
            for e in range(len(cxGroup) - 1):
                arcLength += interp.getCubicHermiteArcLength(cxGroup[e], cd1Group[e],
                                                             cxGroup[e + 1], cd1Group[e + 1])
            arcLengthOfGroupsAlong.append(arcLength)

            if i == 0:
                cx = cxGroup
                cd1 = cd1Group
                cd2 = cd2Group
                cd3 = cd3Group
                cd12 = cd12Group
                cd13 = cd13Group

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
                elementsCountAround, elementsCountAlongSegment, zRefList, innerRadiusAlongElementList[start:end],
                closedProximalEnd=False)

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

        return annotationGroups

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
