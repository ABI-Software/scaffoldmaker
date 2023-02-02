"""
Generates a 3-D bladder mesh along the central line, with variable
numbers of elements around , along and through wall.
"""

import copy
import math

from opencmiss.utils.zinc.field import findOrCreateFieldGroup, findOrCreateFieldStoredString, \
     findOrCreateFieldStoredMeshLocation, findOrCreateFieldNodeGroup
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, \
    getAnnotationGroupForTerm, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.bladder_terms import get_bladder_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues

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
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [[163.52, -70.16, 0.00], [-22.71, 11.06, 0.00], [-7.79, -16.00, 0.00], [-7.18, -15.15, 0.00], [0.00, 0.00, 28.30], [0.00, 0.00, 19.67]],
                [[140.19, -58.93, 0.00], [-23.95, 11.40, 0.00], [-12.93, -27.16, 0.00], [-3.10, -7.17, 0.00], [0.00, 0.00, 41.33], [0.00, 0.00, 6.40]],
                [[115.62, -47.36, 0.00], [-24.13, 11.12, 0.00], [-13.89, -30.15, 0.00], [0.26, -0.17, 0.00], [0.00, 0.00, 40.77], [0.00, 0.00, -4.40]],
                [[91.93, -36.69, 0.00], [-23.87, 10.77, 0.00], [-12.46, -27.62, 0.00], [2.71, 5.15, 0.00], [0.00, 0.00, 32.69], [0.00, 0.00, -8.50]],
                [[67.89, -25.83, 0.00], [-24.55, 10.46, 0.00], [-8.44, -19.81, 0.00], [3.63, 7.32, 0.00], [0.00, 0.00, 23.77], [0.00, 0.00, -9.21]],
                [[42.83, -15.79, 0.00], [-23.65, 9.48, 0.00], [-5.21, -12.99, 0.00], [2.97, 6.53, 0.00], [0.00, 0.00, 14.25], [0.00, 0.00, -7.47]],
                [[20.60, -6.87, 0.00], [-21.44, 7.87, 0.00], [-2.46, -6.71, 0.00], [1.91, 4.18, 0.00], [0.00, 0.00, 8.60], [0.00, 0.00, -3.79]],
                [[0.00, 0.00, 0.00], [-19.75, 5.87, 0.00], [-1.32, -4.43, 0.00], [0.38, 0.37, 0.00], [0.00, 0.00, 6.50], [0.00, 0.00, -0.41]]]),
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
                [[290.00, -124.00, 0.00], [-44.53, 10.62, 2.07], [-14.50, -63.17, 12.12], [-24.99, -23.48, 7.76], [5.15, 10.11, 58.84], [0.77, 5.08, 20.03]],
                [[245.52, -109.41, 1.05], [-44.22, 18.50, 0.02], [-35.41, -84.66, 21.54], [-16.82, -19.48, 11.08], [7.21, 17.14, 79.23], [3.37, 8.99, 20.76]],
                [[202.00, -87.00, 0.00], [-43.33, 20.19, -0.51], [-47.95, -102.05, 34.35], [-4.79, -11.99, 8.62], [11.94, 28.18, 100.37], [3.03, 7.38, 12.77]],
                [[159.00, -69.00, 0.00], [-42.00, 17.50, 0.00], [-45.37, -108.90, 38.98], [1.87, -1.71, 3.74], [13.36, 32.08, 105.18], [0.74, 1.97, -2.61]],
                [[118.00, -52.00, 0.00], [-39.50, 16.50, 0.00], [-44.17, -105.72, 41.88], [4.75, 13.33, -2.01], [13.46, 32.21, 95.50], [-1.53, -4.24, -19.52]],
                [[80.00, -36.00, 0.00], [-39.01, 16.99, 0.00], [-36.12, -82.95, 35.31], [14.15, 35.25, -13.90], [10.41, 23.92, 66.83], [-3.83, -9.68, -28.35]],
                [[40.00, -18.00, 0.00], [-40.00, 18.00, 0.00], [-15.49, -34.42, 13.61], [16.42, 37.84, -16.38], [5.75, 12.77, 38.83], [-4.67, -10.77, -29.35]],
                [[0.00, 0.00, 0.00], [-40.00, 18.00, 0.00], [-3.28, -7.28, 2.56], [8.01, 16.45, -5.73], [1.07, 2.38, 8.14], [-4.68, -10.01, -32.03]]]),
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
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [[65.00, -28.00, 0.00], [-10.45, 9.02, 0.00], [-4.60, -5.33, 1.42], [-2.54, -3.38, -1.33], [0.88, 1.02, 6.72], [-0.90, -0.72, 8.69]],
                [[55.00, -20.00, 0.00], [-9.54, 6.98, 0.00], [-6.37, -8.71, -0.03], [-1.00, -3.38, -1.57], [-0.02, -0.03, 13.28], [-0.92, -1.39, 4.44]],
                [[46.00, -14.00, 0.00], [-10.01, 5.57, 0.00], [-6.72, -12.09, -1.70], [0.54, -2.06, -0.66], [-0.95, -1.71, 15.93], [-0.29, -0.68, 1.17]],
                [[35.00, -9.00, 0.00], [-11.02, 4.50, 0.00], [-5.18, -12.68, -1.23], [1.79, 1.21, 0.19], [-0.52, -1.28, 15.44], [0.21, 0.06, -1.51]],
                [[24.00, -5.00, 0.00], [-9.03, 2.93, 0.00], [-3.15, -9.72, -1.31], [2.01, 3.63, 0.37], [-0.51, -1.58, 12.93], [0.21, 0.31, -4.05]],
                [[17.00, -3.00, 0.00], [-7.52, 1.53, 0.00], [-1.16, -5.67, -0.67], [1.39, 3.65, 0.47], [-0.18, -0.90, 7.93], [0.22, 0.55, -4.93]],
                [[9.00, -2.00, 0.00], [-8.52, 1.48, 0.00], [-0.43, -2.47, -0.40], [0.47, 2.41, 0.25], [-0.08, -0.49, 3.07], [0.05, 0.30, -3.15]],
                [[0.00, 0.00, 0.00], [-9.47, 2.52, 0.00], [-0.26, -0.98, -0.18], [-0.13, 0.57, 0.21], [-0.08, -0.32, 1.88], [-0.05, 0.04, 0.77]]]),
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
                [[417.00, 272.00, 0.00], [-39.04, -53.00, 0.83], [30.11, -22.17, 0.29], [51.78, -52.34, -1.70], [0.03, 0.42, 28.30], [1.06, -1.92, 120.97]],
                [[373.00, 220.00, 0.00], [-49.17, -50.77, -0.19], [66.66, -64.55, -1.68], [21.33, -32.41, -2.24], [1.30, -1.69, 116.41], [1.47, -2.29, 55.25]],
                [[319.00, 171.00, 0.00], [-55.84, -46.38, -0.19], [71.68, -86.28, -4.21], [-0.40, -15.45, 2.49], [3.00, -4.17, 136.47], [-0.79, 2.99, 6.58]],
                [[261.00, 127.00, 0.00], [-62.26, -42.93, 1.27], [65.94, -95.54, 3.23], [-9.21, -3.67, 4.74], [-0.25, 4.21, 129.76], [-1.86, 5.07, -12.37]],
                [[195.00, 86.00, 0.00], [-67.15, -38.10, 1.57], [52.92, -93.08, 5.01], [-16.55, 14.33, 1.48], [-0.59, 5.64, 111.17], [-0.96, 1.25, -22.72]],
                [[128.00, 52.00, 0.00], [-64.44, -31.46, 0.87], [33.03, -67.48, 6.20], [-19.80, 30.23, -1.80], [-2.14, 6.72, 84.52], [0.49, -1.71, -32.39]],
                [[66.00, 23.00, 0.00], [-63.30, -25.43, 1.60], [13.32, -33.05, 1.69], [-14.90, 28.31, -1.76], [0.19, 2.48, 46.94], [0.94, -1.85, -35.99]],
                [[0.00, 0.00, 0.00], [-65.62, -18.77, 3.11], [3.23, -10.86, 2.68], [-5.27, 16.06, 3.74], [-0.27, 3.02, 12.55], [-1.86, 2.93, -32.79]]]),
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
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [[95.34, -57.99, -0.26], [-7.80, 9.56, 0.66], [-8.00, -6.50, -0.44], [-0.46, -2.11, 2.10], [0.00, -0.23, 3.33], [0.92, 0.77, 5.96]],
                [[86.37, -48.34, 0.16], [-10.12, 9.72, 0.18], [-9.06, -9.45, 1.15], [-1.66, -3.81, 1.09], [0.60, 0.47, 8.58], [0.28, 0.63, 4.54]],
                [[75.08, -38.67, 0.06], [-12.19, 9.74, -0.30], [-11.40, -14.21, 1.67], [-1.86, -4.70, -0.44], [0.52, 1.03, 12.31], [-0.27, -0.15, 2.74]],
                [[61.99, -28.92, -0.46], [-13.08, 8.83, -0.01], [-12.72, -18.84, 0.18], [0.31, -2.64, -1.79], [0.05, 0.10, 13.97], [-0.44, -1.05, 0.29]],
                [[49.00, -21.00, 0.00], [-14.49, 8.02, 0.25], [-10.90, -19.63, -1.90], [2.67, 1.37, -0.99], [-0.36, -1.06, 12.98], [-0.24, -0.53, -1.86]],
                [[33.00, -13.00, 0.00], [-15.03, 6.92, 0.00], [-7.23, -15.72, -1.62], [3.94, 6.12, -0.06], [-0.39, -0.86, 10.09], [-0.08, -0.14, -3.63]],
                [[19.00, -7.10, 0.00], [-16.49, 6.55, 0.00], [-3.07, -7.72, -1.96], [3.29, 6.75, 0.54], [-0.51, -1.28, 5.84], [0.12, 0.22, -3.84]],
                [[0.00, 0.00, 0.00], [-21.51, 7.65, 0.00], [-0.94, -2.64, -0.25], [0.97, 3.41, 2.88], [-0.08, -0.21, 2.55], [0.74, 1.92, -2.74]]]),
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
                'Number of elements': 7
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [[2.10, 0.00, 0.00], [-0.30, 0.00, 0.00], [0.00, -0.12, 0.00], [0.00, -0.21, -0.00], [0.00, 0.00, 0.15], [0.00, -0.01, 0.17]],
                [[1.80, 0.00, 0.00], [-0.30, 0.00, 0.00], [0.00, -0.28, -0.00], [0.00, -0.11, -0.00], [0.00, -0.00, 0.28], [0.00, -0.00, 0.09]],
                [[1.50, 0.00, 0.00], [-0.30, 0.00, 0.00], [0.00, -0.34, -0.01], [0.00, -0.04, 0.00], [0.00, -0.01, 0.33], [0.00, 0.00, 0.03]],
                [[1.20, 0.00, 0.00], [-0.30, 0.00, 0.00], [0.00, -0.35, 0.00], [0.00, 0.02, 0.00], [0.00, 0.00, 0.35], [0.00, 0.00, -0.01]],
                [[0.90, 0.00, 0.00], [-0.30, 0.00, 0.00], [0.00, -0.30, 0.00], [0.00, 0.07, 0.00], [0.00, 0.00, 0.30], [0.00, 0.00, -0.10]],
                [[0.60, 0.00, 0.00], [-0.30, 0.00, 0.00], [0.00, -0.20, 0.00], [0.00, 0.11, 0.00], [0.00, 0.00, 0.16], [0.00, 0.00, -0.11]],
                [[0.30, 0.00, 0.00], [-0.30, 0.00, 0.00], [0.00, -0.08, 0.00], [0.00, 0.07, -0.00], [0.00, 0.00, 0.09], [0.00, -0.00, -0.05]],
                [[0.00, 0.00, 0.00], [-0.30, 0.00, 0.00], [0.00, -0.05, 0.00], [0.00, -0.02, -0.00], [0.00, 0.00, 0.05], [0.00, -0.00, -0.03]]]),
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
            'Wall thickness': 1.5,
            'Ureter position around': 0.67,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        if 'Human 1' in parameterSetName:
            options['Number of elements along dome'] = 10
            options['Number of elements along neck'] = 6
            options['Number of elements around'] = 12
            options['Wall thickness'] = 3.0
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
        if 'Mouse 1' in parameterSetName:
            options['Number of elements along neck'] = 6
            options['Wall thickness'] = 0.5
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
        if 'Pig 1' in parameterSetName:
            options['Number of elements along dome'] = 10
            options['Number of elements along neck'] = 5
            options['Number of elements around'] = 12
            options['Wall thickness'] = 2.5
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
        if 'Rat 1' in parameterSetName:
            options['Wall thickness'] = 0.3
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
        if 'Material' in parameterSetName:
            options['Number of elements along dome'] = 8
            options['Number of elements along neck'] = 4
            options['Number of elements around'] = 8
            options['Wall thickness'] = 0.01
            options['Ureter position around'] = 0.67
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
            'Ureter position around',
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
            'Refine number of elements through wall']:
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
        :return: None
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
        materialWallThickness = 0.01

        elementsCountAlongBladder = elementsCountAlongDome + elementsCountAlongNeck

        fm = region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()
        mesh = fm.findMeshByDimension(3)

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Geometric coordinates
        geometricCentralPath = BladderCentralPath(region, centralPath, elementsCountAlongDome, elementsCountAlongNeck)
        xFinal, d1Final, d2Final, d3Final = findBladderNodes3D(elementsCountAlongDome, elementsCountAlongNeck, elementsCountAround,
                                            elementsCountThroughWall, wallThickness, geometricCentralPath)

        sx_dome_group = geometricCentralPath.sx_dome_group
        sx_neck_group = geometricCentralPath.sx_neck_group
        bladderCentralPathLength = geometricCentralPath.bladderCentralPathLength
        domeLength = geometricCentralPath.domeLength

        # Material coordinates
        tmp_region = region.createRegion()
        materialCentralPath = BladderCentralPath(tmp_region, materialCentralPath, elementsCountAlongDome, elementsCountAlongNeck)
        xOrgan, d1Organ, d2Organ, d3Organ = findBladderNodes3D(elementsCountAlongDome, elementsCountAlongNeck, elementsCountAround,
                                            elementsCountThroughWall, materialWallThickness, materialCentralPath)
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
        annotationGroupAround = [[bladderGroup, bladderVentralGroup], [bladderGroup, bladderDorsalGroup], [bladderGroup, bladderVentralGroup]]
        annotationGroupsAround = []
        for i in range(len(elementsCountAroundGroups)):
            elementsCount = elementsCountAroundGroups[i]
            for n in range(elementsCount):
                annotationGroupsAround.append(annotationGroupAround[i])

        annotationGroupsThroughWall = []
        for i in range(elementsCountThroughWall):
            annotationGroupsThroughWall.append([])

        # Flat coordinates
        urethraOpeningRadius = vector.magnitude(sx_neck_group[2][-1])
        xFlat, d1Flat, d2Flat = obtainBladderFlatNodes(elementsCountAlongBladder, elementsCountAround, elementsCountThroughWall,
                                                       xFinal, d1Final, d2Final, bladderCentralPathLength, urethraOpeningRadius, wallThickness)

        # Create nodes and elements
        nodeIdentifier, elementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
            region, xFinal, d1Final, d2Final, d3Final, xFlat, d1Flat, d2Flat, xOrgan, d1Organ, d2Organ, "bladder coordinates",
            elementsCountAround, elementsCountAlongBladder, elementsCountThroughWall,
            annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier,
            useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd=True)

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
        ureter1Position = trackSurfaceBladder.createPositionProportion(ureterPositionAroundFactor, neckStartPositionAlongFactor)
        ureterElementPositionAround = ureter1Position.e1

        # Annotation fiducial point
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        # Define markers for apex and ureter junctions with bladder
        apexGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("apex of urinary bladder"))
        leftUreterGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("left ureter junction with bladder"))
        rightUreterGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("right ureter junction with bladder"))

        idx1 = 1
        xi1 = [0.0, 0.0, 0.0]
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
            markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
            cache.setNode(markerPoint)
            markerLocation.assignMeshLocation(cache, mesh.findElementByIdentifier(marker["elementId"]), marker["xi"])
            markerName.assignString(cache, annotationGroup.getName())
            annotationGroup.getNodesetGroup(nodes).addNode(markerPoint)
            bladderNodesetGroup.addNode(markerPoint)
            nodeIdentifier += 1

        fm.endChange()
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

        is_body = bodyGroup.getFieldElementGroup(mesh2d)
        is_body_serosa = fm.createFieldAnd(is_body, is_exterior_face_xi3_1)
        is_body_lumen = fm.createFieldAnd(is_body, is_exterior_face_xi3_0)

        is_neck = neckGroup.getFieldElementGroup(mesh2d)
        is_neck_serosa = fm.createFieldAnd(is_neck, is_exterior_face_xi3_1)
        is_neck_lumen = fm.createFieldAnd(is_neck, is_exterior_face_xi3_0)

        is_urinaryBladder = urinaryBladderGroup.getFieldElementGroup(mesh2d)
        is_urinaryBladder_serosa = fm.createFieldAnd(is_urinaryBladder, is_exterior_face_xi3_1)
        is_urinaryBladder_lumen = fm.createFieldAnd(is_urinaryBladder, is_exterior_face_xi3_0)

        is_dorsal_bladder = bladderDorsalGroup.getFieldElementGroup(mesh2d)
        is_ventral_bladder = bladderVentralGroup.getFieldElementGroup(mesh2d)

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
    d = [c * arcLengthAround for c in vector.normalise(d)]

    return d

def findNodesAlongBladderDome(sx_dome_group, elementsCountAround, elementsCountAlongDome):

    # Find apex nodes d2
    d2Apex = []
    d2 = sx_dome_group[2][0]
    for n1 in range(elementsCountAround):
        rotAngle = n1 * 2.0 * math.pi / elementsCountAround
        rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(sx_dome_group[2][0]), vector.normalise(sx_dome_group[4][0])))
        rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
        d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
        d2Apex.append(d2Rot)

    # Create ellipses along dome around the central path
    xEllipses_dome = []
    d1Ellipses_dome = []
    for n in range(1, len(sx_dome_group[0])):
        px, pd1 = createEllipsePoints(sx_dome_group[0][n], 2 * math.pi, sx_dome_group[2][n], sx_dome_group[4][n], elementsCountAround,
                                      startRadians=0.0)
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
        xSampledAlongDome, d2SampledAlongDome = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                                elementsCountAlongDome,
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

def findNodesAlongBladderNeck(sx_dome_group, sx_neck_group, d2SampledDome, domeSegmentLength, neckSegmentLength, elementsCountAround, elementsCountAlongNeck):

    # Transition
    transitLength = (domeSegmentLength + neckSegmentLength) / 2
    addingLength = transitLength - neckSegmentLength
    if transitLength <= neckSegmentLength:
        e = 1
        xi = transitLength / neckSegmentLength
    else:
        e = int(transitLength / neckSegmentLength) + 1
        xi = transitLength / (e * neckSegmentLength)
    xTransition = interp.interpolateCubicHermite(sx_dome_group[0][-1], sx_dome_group[1][-1], sx_neck_group[0][e], sx_neck_group[1][e], xi)
    # d1Transition = interp.interpolateCubicHermiteDerivative(sx_dome_group[0][-1], sx_dome_group[1][-1], sx_neck_group[0][e], sx_neck_group[1][e], xi)
    d2Transition = interp.interpolateCubicHermite(sx_dome_group[2][-1], sx_dome_group[3][-1], sx_neck_group[2][e], sx_neck_group[3][e], xi)
    d3Transition = interp.interpolateCubicHermite(sx_dome_group[4][-1], sx_dome_group[5][-1], sx_neck_group[4][e], sx_neck_group[5][e], xi)
    px_transit, pd1_transit = createEllipsePoints(xTransition, 2 * math.pi, d2Transition, d3Transition, elementsCountAround,
                                                  startRadians=0.0)

    # Create ellipses along neck around the central path
    xEllipses_neck = []
    d1Ellipses_neck = []
    for n in range(0, len(sx_neck_group[0])):
        px, pd1 = createEllipsePoints(sx_neck_group[0][n], 2 * math.pi, sx_neck_group[2][n], sx_neck_group[4][n], elementsCountAround,
                                      startRadians=0.0)
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
            else:
                v1 = xEllipses_neck[n2 - 1][n1]
                v2 = xEllipses_neck[n2][n1]
                d2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v2)
            d2Along.append(d2)
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
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
        xSampledAlongNeck, d2SampledAlongNeck = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                                elementsCountAlongNeck, addLengthStart=addingLength,
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

    return xSampledNeck, d1SampledNeck, d2SampledNeck, px_transit

def obtainBladderFlatNodes(elementsCountAlongBladder, elementsCountAround, elementsCountThroughWall,
                           xFinal, d1Final, d2Final, bladderLength, urethraOpeningRadius, bladderWallThickness):
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
    alpha = vector.angleBetweenVectors(v1, v2)

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
        smd2 = interp.smoothCubicHermiteDerivativesLine(lineSmoothingNodes, lineSmoothingNodes_d2, fixAllDirections=False,
                                                 fixStartDerivative=False, fixEndDerivative=True,
                                                 fixStartDirection=False, fixEndDirection=False)
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
    d1Mag = vector.magnitude(v21)
    d2Mag = vector.magnitude(v31)

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
        cx_dome, cd1_dome, cd2_dome, cd3_dome, cd12_dome, cd13_dome = \
            extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                        Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                        Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                            groupName='dome of the bladder')
        cx_neck, cd1_neck, cd2_neck, cd3_neck, cd12_neck, cd13_neck = \
            extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                        Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                        Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                            groupName='neck of urinary bladder')
        # Find arcLength
        # Dome
        domeLength = 0.0
        elementsCountInDome = len(cx_dome) - 1
        for e in range(elementsCountInDome):
            arcLength = interp.getCubicHermiteArcLength(cx_dome[e], cd1_dome[e],
                                                        cx_dome[e + 1], cd1_dome[e + 1])
            domeLength += arcLength
        domeSegmentLength = domeLength / elementsCountAlongDome
        # Neck
        neckLength = 0.0
        elementsCountInNeck = len(cx_neck) - 1
        for e in range(elementsCountInNeck):
            arcLength = interp.getCubicHermiteArcLength(cx_neck[e], cd1_neck[e],
                                                        cx_neck[e + 1], cd1_neck[e + 1])
            neckLength += arcLength
        neckSegmentLength = neckLength / elementsCountAlongNeck
        bladderCentralPathLength = domeLength + neckLength
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
        self.sx_dome_group = [sx_dome, sd1_dome, sd2_dome, sd12_dome, sd3_dome, sd13_dome]
        self.sx_neck_group = [sx_neck, sd1_neck, sd2_neck, sd12_neck, sd3_neck, sd13_neck]

        self.domeSegmentLength = domeSegmentLength
        self.neckSegmentLength = neckSegmentLength
        self.bladderCentralPathLength = bladderCentralPathLength
        self.domeLength = domeLength

def findBladderNodes3D(elementsCountAlongDome, elementsCountAlongNeck, elementsCountAround, elementsCountThroughWall,
                        wallThickness, centralPath):

    sx_dome_group = centralPath.sx_dome_group
    sx_neck_group = centralPath.sx_neck_group
    domeSegmentLength = centralPath.domeSegmentLength
    neckSegmentLength = centralPath.neckSegmentLength

    elementsCountAlongBladder = elementsCountAlongDome + elementsCountAlongNeck

    xDome, d1Dome, d2Dome = findNodesAlongBladderDome(sx_dome_group, elementsCountAround, elementsCountAlongDome)

    xNeck, d1Neck, d2Neck, px_transit = findNodesAlongBladderNeck(sx_dome_group, sx_neck_group, d2Dome, domeSegmentLength, \
                                                                  neckSegmentLength, elementsCountAround, elementsCountAlongNeck)

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
            d3Around.append(vector.normalise(
                vector.crossproduct3(vector.normalise(d1SampledAll[n2][n1]), vector.normalise(d2SampledAll[n2][n1]))))
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

    # Create outer layers from the inner nodes
    wallThicknessList = [wallThickness] * (elementsCountAlongBladder + 1)
    relativeThicknessList = []
    transitElementList = [0] * elementsCountAround
    xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xInner, d1Inner,
                                                                                    d2Inner, d3Inner,
                                                                                    wallThicknessList,
                                                                                    relativeThicknessList,
                                                                                    elementsCountAround,
                                                                                    elementsCountAlongBladder,
                                                                                    elementsCountThroughWall,
                                                                                    transitElementList)

    # Deal with multiple nodes at the start point for closed proximal end
    xApexInner = xList[0]
    # Arclength between apex point and corresponding point on next face
    mag = interp.getCubicHermiteArcLength(xList[0], d2List[0], xList[2 * elementsCountAround],
                                          d2List[2 * elementsCountAround])
    d2ApexInner = vector.setMagnitude(sx_dome_group[2][0], mag)

    d1ApexInner = vector.crossproduct3(sx_dome_group[1][0], d2ApexInner)
    d1ApexInner = vector.setMagnitude(d1ApexInner, mag)
    d3ApexUnit = vector.normalise(
        vector.crossproduct3(vector.normalise(d1ApexInner), vector.normalise(d2ApexInner)))
    d3ApexInner = [d3ApexUnit[c] * wallThickness / elementsCountThroughWall for c in range(3)]

    # Final nodes on the bladder
    xFinal = []
    d1Final = []
    d2Final = []
    d3Final = []
    for n3 in range(elementsCountThroughWall + 1):
        xApex = [xApexInner[c] +
                 d3ApexUnit[c] * wallThickness / elementsCountThroughWall * n3 for c in range(3)]
        xFinal.append(xApex)
        d1Final.append(d1ApexInner)
        d2Final.append(d2ApexInner)
        d3Final.append(d3ApexInner)

    xFinal += xList[(elementsCountThroughWall + 1) * elementsCountAround:]
    d1Final += d1List[(elementsCountThroughWall + 1) * elementsCountAround:]
    d2Final += d2List[(elementsCountThroughWall + 1) * elementsCountAround:]
    d3Final += d3List[(elementsCountThroughWall + 1) * elementsCountAround:]

    return xFinal, d1Final, d2Final, d3Final
