"""
Generates a 3-D bladder mesh along the central line, with variable
numbers of elements around , along and through wall.
"""

from __future__ import division

import copy
import math

from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldGroup, \
    findOrCreateFieldStoredString, findOrCreateFieldStoredMeshLocation, findOrCreateFieldNodeGroup
from opencmiss.utils.zinc.finiteelement import get_element_node_identifiers
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, mergeAnnotationGroups, \
    getAnnotationGroupForTerm, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.bladder_terms import get_bladder_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues, mesh_destroy_elements_and_nodes_by_identifiers


class MeshType_3d_bladder1(Scaffold_base):
    """
    Generates a 3-D bladder mesh with variable numbers of elements around, along the central line, and through wall.
    The bladder is created using a central path as the longitudinal axis of the bladder. Magnitude of D2 and D3 are
    the radii of the bladder in the respective direction.
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
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [[70.80, 72.30, 0.00], [-6.62, -16.41, -1.20], [38.60, -15.67, 1.49], [-8.82, -10.04, -1.62], [-2.26, -1.91, 38.59], [4.19, 4.77, 4.06]],
                    [[57.51, 51.14, -0.10], [-19.84, -25.60, 1.02], [30.18, -23.37, 0.26], [-8.02, -5.36, -0.83], [0.57, 1.21, 41.38], [1.49, 1.45, 1.53]],
                    [[30.29, 22.36, 2.55], [-25.83, -23.54, -0.03], [22.79, -25.01, 0.06], [-8.14, 1.29, 0.15], [-0.08, 0.03, 40.90], [-1.33, -0.83, -5.13]],
                    [[6.64, 3.83, 0.48], [-22.06, -14.48, -1.56], [14.08, -21.50, 0.48], [-6.95, 3.05, 0.10], [-1.93, -0.54, 32.24], [-0.40, -0.03, -9.52]],
                    [[-13.27, -6.98, -0.63], [-21.39, -9.64, -1.27], [8.45, -18.80, 0.34], [-4.94, 3.33, 0.35], [-1.24, -0.16, 22.07], [0.24, 0.77, -7.46]],
                    [[-36.03, -15.27, -2.05], [-21.61, -6.31, -1.44], [4.24, -14.80, 1.22], [-3.51, 4.96, 0.69], [-1.47, 1.02, 17.51], [0.02, 0.97, -6.00]],
                    [[-56.31, -19.78, -3.49], [-25.21, -4.17, -2.30], [1.33, -9.02, 1.75], [-2.24, 4.68, 0.72], [-1.24, 1.81, 10.27], [0.21, 0.99, -5.74]],
                    [[-86.41, -23.06, -6.84], [-34.94, -2.39, -4.40], [0.06, -5.95, 2.75], [-0.30, 1.46, 1.27], [-1.06, 3.10, 6.73], [0.15, 1.57, -1.34]]]),

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
            'meshEdits': exnodeStringFromNodeValues(
                    [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                    [ [  70.80,  72.30,  0.00], [  -6.62, -16.41, -1.20], [ 38.60, -15.67, 1.49], [ -8.82, -10.04, -1.62], [ -2.26, -1.91, 38.59], [  4.19,  4.77,  4.06] ],
                    [ [  57.51,  51.14, -0.10], [ -19.84, -25.60,  1.02], [ 30.18, -23.37, 0.26], [ -8.02,  -5.36, -0.83], [  0.57,  1.21, 41.38], [  1.49,  1.45,  1.53] ],
                    [ [  30.29,  22.36,  2.55], [ -25.83, -23.54, -0.03], [ 22.79, -25.01, 0.06], [ -8.14,   1.29,  0.15], [ -0.08,  0.03, 40.90], [ -1.33, -0.83, -5.13] ],
                    [ [   6.64,   3.83,  0.48], [ -22.06, -14.48, -1.56], [ 14.08, -21.50, 0.48], [ -6.95,   3.05,  0.10], [ -1.93, -0.54, 32.24], [ -0.40, -0.03, -9.52] ],
                    [ [ -13.27,  -6.98, -0.63], [ -21.39,  -9.64, -1.27], [  8.45, -18.80, 0.34], [ -4.94,   3.33,  0.35], [ -1.24, -0.16, 22.07], [  0.24,  0.77, -7.46] ],
                    [ [ -36.03, -15.27, -2.05], [ -21.61,  -6.31, -1.44], [  4.24, -14.80, 1.22], [ -3.51,   4.96,  0.69], [ -1.47,  1.02, 17.51], [  0.02,  0.97, -6.00] ],
                    [ [ -56.31, -19.78, -3.49], [ -25.21,  -4.17, -2.30], [  1.33,  -9.02, 1.75], [ -2.24,   4.68,  0.72], [ -1.24,  1.81, 10.27], [  0.21,  0.99, -5.74] ],
                    [ [ -86.41, -23.06, -6.84], [ -34.94,  -2.39, -4.40], [  0.06,  -5.95, 2.75], [ -0.30,   1.46,  1.27], [ -1.06,  3.10,  6.73], [  0.15,  1.57, -1.34] ] ] ),

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
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
                 Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [  0.9,  3.7, 0.0 ], [ -0.8, -3.6, 0.0 ], [  3.2, -0.6, 0.0 ], [ -1.3, -0.5, 0.0 ], [ 0.0, 0.0, 2.6 ], [ 0.0, 0.0,  0.9 ] ],
                [ [ -0.1,  0.7, 0.0 ], [ -1.2, -2.4, 0.0 ], [  2.0, -1.5, 0.0 ], [ -1.1, -1.3, 0.0 ], [ 0.0, 0.0, 3.1 ], [ 0.0, 0.0,  0.1 ] ],
                [ [ -1.4, -1.1, 0.0 ], [ -1.6, -1.1, 0.0 ], [  1.0, -3.0, 0.0 ], [ -1.3, -0.8, 0.0 ], [ 0.0, 0.0, 3.0 ], [ 0.0, 0.0, -0.2 ] ],
                [ [ -2.9, -1.6, 0.0 ], [ -1.6,  0.2, 0.0 ], [ -0.6, -3.3, 0.0 ], [ -1.4,  0.2, 0.0 ], [ 0.0, 0.0, 2.8 ], [ 0.0, 0.0, -0.1 ] ],
                [ [ -4.3, -0.8, 0.0 ], [ -1.2,  1.1, 0.0 ], [ -1.8, -2.5, 0.0 ], [ -0.8,  1.1, 0.0 ], [ 0.0, 0.0, 2.9 ], [ 0.0, 0.0, -0.1 ] ],
                [ [ -5.2,  0.6, 0.0 ], [ -0.8,  1.6, 0.0 ], [ -2.2, -1.1, 0.0 ], [  0.2,  1.1, 0.0 ], [ 0.0, 0.0, 2.5 ], [ 0.0, 0.0, -0.7 ] ],
                [ [ -5.9,  2.3, 0.0 ], [ -0.5,  1.3, 0.0 ], [ -1.3, -0.4, 0.0 ], [  0.6,  0.3, 0.0 ], [ 0.0, 0.0, 1.4 ], [ 0.0, 0.0, -0.7 ] ],
                [ [ -6.2,  3.2, 0.0 ], [ -0.4,  0.9, 0.0 ], [ -0.8, -0.3, 0.0 ], [  0.1, -0.0, 0.0 ], [ 0.0, 0.0, 0.9 ], [ 0.0, 0.0, -0.2 ] ],
                [ [ -6.8,  4.1, 0.0 ], [ -0.7,  0.9, 0.0 ], [ -1.1, -0.5, 0.0 ], [ -0.7, -0.4, 0.0 ], [ 0.0, 0.0, 1.1 ], [ 0.0, 0.0,  0.6 ] ] ] ),

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
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [  73.1,  50.2, 0.0 ], [ -18.5, -35.4, 0.0 ], [  18.5, -10.9, 0.0 ], [  18.6, -8.2, 0.0 ], [ 0.0, 0.0, 27.8 ], [ 0.0, 0.0,  8.6 ] ],
                [ [  57.3,  20.3, 0.0 ], [ -13.1, -24.4, 0.0 ], [  30.1, -17.0, 0.0 ], [   4.6, -4.0, 0.0 ], [ 0.0, 0.0, 33.3 ], [ 0.0, 0.0,  2.4 ] ],
                [ [  47.0,   1.4, 0.0 ], [ -12.6, -19.8, 0.0 ], [  30.2, -19.7, 0.0 ], [  -4.3, -4.5, 0.0 ], [ 0.0, 0.0, 33.7 ], [ 0.0, 0.0, -0.7 ] ],
                [ [  32.0, -18.9, 0.0 ], [ -19.5, -14.4, 0.0 ], [  20.7, -26.4, 0.0 ], [ -13.7, -4.9, 0.0 ], [ 0.0, 0.0, 31.6 ], [ 0.0, 0.0, -3.7 ] ],
                [ [  10.7, -26.3, 0.0 ], [ -24.3,   1.9, 0.0 ], [   3.1, -29.7, 0.0 ], [ -16.7,  4.8, 0.0 ], [ 0.0, 0.0, 26.5 ], [ 0.0, 0.0, -8.8 ] ],
                [ [ -11.3, -14.4, 0.0 ], [ -14.4,  19.6, 0.0 ], [ -12.7, -15.9, 0.0 ], [  -4.1, 13.5, 0.0 ], [ 0.0, 0.0, 13.5 ], [ 0.0, 0.0, -7.8 ] ],
                [ [ -15.8,   7.8, 0.0 ], [  -8.3,  18.3, 0.0 ], [  -6.4,  -2.7, 0.0 ], [   2.8,  4.4, 0.0 ], [ 0.0, 0.0, 10.4 ], [ 0.0, 0.0, -1.7 ] ],
                [ [ -26.2,  21.4, 0.0 ], [ -11.8,   8.4, 0.0 ], [  -6.3,  -4.9, 0.0 ], [  -2.6, -8.8, 0.0 ], [ 0.0, 0.0,  9.8 ], [ 0.0, 0.0,  0.5 ] ] ] ),

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
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [  11.3, 13.4, 0.0 ], [  0.3, -13.4, 0.0 ], [  9.4,  -1.3, 0.0 ], [ -1.0, -5.6, 0.0 ], [ 0.0, 0.0, 7.4 ], [ 0.0, 0.0,  1.5 ] ],
                [ [   9.3,  2.1, 0.0 ], [ -4.4,  -8.7, 0.0 ], [  7.4,  -6.1, 0.0 ], [ -3.0, -3.9, 0.0 ], [ 0.0, 0.0, 8.5 ], [ 0.0, 0.0,  0.7 ] ],
                [ [   4.0, -3.6, 0.0 ], [ -6.8,  -3.8, 0.0 ], [  3.7,  -9.4, 0.0 ], [ -4.9, -2.4, 0.0 ], [ 0.0, 0.0, 9.0 ], [ 0.0, 0.0,  0.1 ] ],
                [ [  -3.4, -5.1, 0.0 ], [ -6.4,   0.6, 0.0 ], [ -2.4, -10.9, 0.0 ], [ -5.0,  0.9, 0.0 ], [ 0.0, 0.0, 8.8 ], [ 0.0, 0.0, -0.5 ] ],
                [ [  -8.1, -3.2, 0.0 ], [ -4.4,   3.8, 0.0 ], [ -6.7,  -8.3, 0.0 ], [ -2.5,  3.4, 0.0 ], [ 0.0, 0.0, 8.1 ], [ 0.0, 0.0, -1.2 ] ],
                [ [ -11.4,  2.3, 0.0 ], [ -1.4,   6.4, 0.0 ], [ -6.9,  -4.0, 0.0 ], [  1.9,  4.2, 0.0 ], [ 0.0, 0.0, 6.2 ], [ 0.0, 0.0, -2.8 ] ],
                [ [ -10.7,  8.9, 0.0 ], [  0.3,   5.0, 0.0 ], [ -2.9,   0.0, 0.0 ], [  0.9,  1.1, 0.0 ], [ 0.0, 0.0, 2.4 ], [ 0.0, 0.0, -0.6 ] ],
                [ [ -10.7, 12.2, 0.0 ], [ -0.3,   3.0, 0.0 ], [ -3.5,  -0.3, 0.0 ], [ -0.3, -0.1, 0.0 ], [ 0.0, 0.0, 3.4 ], [ 0.0, 0.0,  0.4 ] ],
                [ [ -11.3, 14.8, 0.0 ], [ -0.9,   2.2, 0.0 ], [ -3.5,  -0.3, 0.0 ], [  0.3,  0.1, 0.0 ], [ 0.0, 0.0, 3.4 ], [ 0.0, 0.0, -0.4 ] ] ] ),

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
            'Rat 1']

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
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages_Bladder['Cat 1']
        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements along bladder': 12,
            'Number of elements around bladder': 8,
            'Number of elements through wall': 1,
            'Wall thickness': 0.5,
            'Ureter position around': 0.67,
            'Ureter position down': 0.8,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements along bladder': 4,
            'Refine number of elements around bladder': 4,
            'Refine number of elements through wall': 1
        }
        if 'Human 1' in parameterSetName:
            options['Number of elements along bladder'] = 8
            options['Number of elements around bladder'] = 12
            options['Wall thickness'] = 1.0
            options['Ureter position around'] = 0.82  # should be on the dorsal part (> 0.5)
            options['Ureter position down'] = 0.63
        if 'Mouse 1' in parameterSetName:
            options['Number of elements along bladder'] = 10
            options['Number of elements around bladder'] = 4
            options['Wall thickness'] = 0.2
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
            options['Ureter position down'] = 0.865
        if 'Pig 1' in parameterSetName:
            options['Number of elements along bladder'] = 6
            options['Number of elements around bladder'] = 16
            options['Wall thickness'] = 2.0
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
            options['Ureter position down'] = 0.865
        if 'Rat 1' in parameterSetName:
            options['Number of elements along bladder'] = 12
            options['Number of elements around bladder'] = 12
            options['Wall thickness'] = 0.1
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
            options['Ureter position down'] = 0.83
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Central path',
            'Number of elements along bladder',
            'Number of elements around bladder',
            'Number of elements through wall',
            'Wall thickness',
            'Ureter position around',
            'Ureter position down',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around bladder',
            'Refine number of elements along bladder',
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
            'Number of elements along bladder',
            'Number of elements around bladder',
            'Number of elements through wall',
            'Refine number of elements along bladder',
            'Refine number of elements around bladder',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        # if options['Number of elements around bladder'] % 2:
        #     options['Number of elements around bladder'] += 1
        if options['Number of elements around bladder'] % 2 != 0:
            options['Number of elements around bladder'] += 1
        if options['Ureter position around'] < 0.5:
            options['Ureter position around'] = 0.5  # ureters are on the dorsal part of the bladder
        # elif options['Ureter position around'] > 0.9:
        #     options['Ureter position around'] = 0.9
        # if options['Ureter position down'] < 0.15:
        #     options['Ureter position down'] = 0.15
        # elif options['Ureter position down'] > 0.95:
        #     options['Ureter position down'] = 0.95

        # if options['Number of elements around bladder'] < 12:
        #     options['Number of elements around bladder'] = 12
        # if options['Number of elements through wall'] != (1 or 4):
        #     options['Number of elements through wall'] = 4

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        centralPath = options['Central path']
        elementsCountAlongBladder = options['Number of elements along bladder']
        elementsCountAroundBladder = options['Number of elements around bladder']
        elementsCountThroughWall = options['Number of elements through wall']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = False
        useCubicHermiteThroughWall = not (options['Use linear through wall'])

        fm = region.getFieldmodule()
        fm.beginChange()

        coordinates = findOrCreateFieldCoordinates(fm)
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        if useCubicHermiteThroughWall:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
            if useCrossDerivatives:
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
                nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        cache = fm.createFieldcache()
        mesh = fm.findMeshByDimension(3)

        nodeIdentifier = 1
        elementIdentifier = 1

        # Extract length of each group along bladder from central path
        arcLengthOfGroupsAlong = []
        bladderTermsAlong = [None, 'dome of the bladder', 'neck of urinary bladder']
        for i in range(len(bladderTermsAlong)):
            tmpRegion = region.createRegion()
            centralPath.generate(tmpRegion)
            cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                            Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                                groupName=bladderTermsAlong[i])

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

        print('arcLengthOfGroupsAlong', arcLengthOfGroupsAlong)
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, len(cx))
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        sd3, sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, se, sxi, ssf)

        # Create annotation groups for bladder
        bladderCentralPathLength = sum(arcLengthOfGroupsAlong[1:])
        allAnnotationGroups = []
        bodyGroup = AnnotationGroup(region, get_bladder_term("dome of the bladder"))
        neckGroup = AnnotationGroup(region, get_bladder_term("neck of urinary bladder"))
        bladderGroup = AnnotationGroup(region, get_bladder_term("urinary bladder"))

        annotationGroupAlong = [[bladderGroup, bodyGroup], [bladderGroup, neckGroup]]
        annotationGroupsAround = [[]]
        for i in range(elementsCountAroundBladder):
            annotationGroupsAround.append([])
        annotationGroupsThroughWall = []
        for i in range(elementsCountThroughWall):
            annotationGroupsThroughWall.append([])

        # annotation fiducial points
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        # Fundus diameter
        fundusRadius = vector.magnitude(sd2[0])
        lengthElementAlongTrackSurface = bladderCentralPathLength / elementsCountAlongBladder
        elementsAlongFundus = int(fundusRadius / lengthElementAlongTrackSurface)

        d2Apex = []
        d2 = sd2[0]
        # print('sd1[0]', sd1[0])
        # print('sd2[0]', sd2[0])
        # print('vector.normalise(sd2[0]', vector.normalise(sd2[0]))
        # print('sd3[0]', sd3[0])
        # print('vector.normalise(sd3[0]', vector.normalise(sd3[0]))
        for n1 in range(elementsCountAroundBladder):
            rotAngle = n1 * 2.0 * math.pi / elementsCountAroundBladder
            rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(sd2[0]), vector.normalise(sd3[0])))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
            d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
            d2Apex.append(d2Rot)
            # if n1 == 1:
            #     print('rotAngle', rotAngle)
            #     print('rotAxis', rotAxis)
            #     print('rotFrame', rotFrame)
            #     print('d2Rot', d2Rot)

        xEllipses = []
        d1Ellipses = []
        print('len(sx)', len(sx))
        # print('elementsAlongFundus', elementsAlongFundus)
        for n in range(1, len(sx)):
            px, pd1 = createEllipsePoints(sx[n], 2 * math.pi, sd2[n], sd3[n], elementsCountAroundBladder,
                                          startRadians=0.0)
            xEllipses.append(px)
            d1Ellipses.append(pd1)
        # print('xEllipses', xEllipses)
        print('len(xEllipses)', len(xEllipses))
        # print('len(xEllipses[0])', len(xEllipses[0]))
        # print('len(xEllipses[0][0])', len(xEllipses[0][0]))

        # Find d2
        d2Raw = []
        for n1 in range(elementsCountAroundBladder):
            xAlong = []
            d2Along = []
            for n2 in range(len(xEllipses) - 1):
                v1 = xEllipses[n2][n1]
                v2 = xEllipses[n2 + 1][n1]
                d2 = findDerivativeBetweenPoints(v1, v2)
                xAlong.append(v1)
                d2Along.append(d2)
            xAlong.append(xEllipses[-1][n1])
            d2Along.append(d2)
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
            d2Raw.append(d2Smoothed)

        # Rearrange d2
        d2Ellipses = []
        for n2 in range(len(xEllipses)):
            d2Around = []
            for n1 in range(elementsCountAroundBladder):
                d2 = d2Raw[n1][n2]
                d2Around.append(d2)
            d2Ellipses.append(d2Around)

        # Merge fundus apex and body
        xAll = [[sx[0]] * elementsCountAroundBladder] + xEllipses
        d2All = [d2Apex] + d2Ellipses
        print('len(xAll)', len(xAll))


        # Spread out elements
        xRaw = []
        d2Raw = []
        for n1 in range(elementsCountAroundBladder):
            xAlong = []
            d2Along = []
            for n2 in range(len(xAll)):
                xAlong.append(xAll[n2][n1])
                d2Along.append(d2All[n2][n1])
            xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                            elementsCountAlongBladder,
                                                                            arcLengthDerivatives=True)[0:2]
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlong, d2SampledAlong)
            xRaw.append(xSampledAlong)
            d2Raw.append(d2Smoothed)

        # Rearrange x and d2
        xSampledAll = []
        d1SampledAll = []
        d2SampledAll = []
        for n2 in range(elementsCountAlongBladder + 1):
            xAround = []
            d1Around = []
            d2Around = []
            for n1 in range(elementsCountAroundBladder):
                x = xRaw[n1][n2]
                d2 = d2Raw[n1][n2]
                xAround.append(x)
                d2Around.append(d2)

                # Calculate d1
                if n2 > 0:
                    v1 = xRaw[n1][n2]
                    v2 = xRaw[n1 + 1 if n1 < elementsCountAroundBladder - 2 else 0][n2]
                    d1 = findDerivativeBetweenPoints(v1, v2)
                    d1Around.append(d1)
                else:
                    d1Around.append(d2Raw[int(elementsCountAroundBladder * 0.75)][0])

            if n2 > 0:
                d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
            else:
                d1Smoothed = d1Around

            xSampledAll.append(xAround)
            d1SampledAll.append(d1Smoothed)
            d2SampledAll.append(d2Around)
        print('len(xSampledAll)', len(xSampledAll))
        print('', xSampledAll[0])

        for n2 in range(elementsCountAlongBladder + 1):
            # if n2 == 0:
            #     node = nodes.createNode(nodeIdentifier, nodetemplate)
            #     cache.setNode(node)
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledAll[n2][0])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledAll[n2][0])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledAll[n2][0])
            #     nodeIdentifier += 1
            # else:
            for n1 in range(elementsCountAroundBladder):
                node = nodes.createNode(nodeIdentifier, nodetemplate)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledAll[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledAll[n2][n1])
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledAll[n2][n1])
                # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3List[n2])
                # if useCrossDerivatives:
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                nodeIdentifier += 1


        fm.endChange()

        return allAnnotationGroups

    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements surface']
        refineElementsCountAlong = options['Refine number of elements surface']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return

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
