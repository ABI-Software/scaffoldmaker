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
                [[70.80, 72.30, 0.00], [-20.79, -22.78, 0.52], [18.43, -16.77, 1.98], [1.58, -5.87, -2.59], [-1.33, 1.87, 28.21], [1.21, -2.78, 20.33]],
                [[49.64, 51.56, 0.11], [-21.49, -18.65, -0.31], [19.75, -22.76, 0.28], [1.06, -6.10, -0.81], [-0.59, -0.01, 41.40], [0.27, -0.96, 6.06]],
                [[28.03, 34.99, -0.59], [-23.01, -16.34, -0.47], [20.57, -28.97, 0.28], [-2.00, -1.83, 0.44], [-0.75, -0.13, 40.89], [-0.26, 0.40, -4.19]],
                [[3.63, 18.98, -0.82], [-26.42, -15.76, -0.50], [15.55, -26.11, 1.19], [-6.09, 5.90, -0.06], [-1.12, 0.83, 32.77], [0.03, -0.01, -9.44]],
                [[-24.82, 3.62, -1.62], [-23.44, -11.53, -0.83], [8.28, -16.84, 0.05], [-4.41, 6.01, -0.02], [-0.65, -0.25, 21.87], [0.08, 0.02, -8.87]],
                [[-43.01, -4.49, -2.39], [-19.63, -8.66, -0.84], [5.63, -12.84, 0.73], [-2.80, 5.02, 0.06], [-0.81, 0.45, 14.26], [0.01, 0.10, -6.72]],
                [[-64.08, -13.68, -3.31], [-22.03, -8.81, -1.71], [2.65, -6.64, 0.07], [-2.12, 4.29, -0.34], [-0.60, -0.16, 8.58], [-0.01, -0.33, -3.96]],
                [[-87.03, -22.04, -5.85], [-23.83, -7.91, -3.38], [1.45, -4.39, 0.06], [-0.29, 0.22, 0.33], [-0.85, -0.19, 6.45], [-0.48, 0.26, -0.30]]]),
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
                [[139.76, 57.34, -16.27], [-48.91, -19.25, 7.52], [25.22, -60.05, 10.30], [13.84, -27.31, 3.97], [4.34, 11.87, 58.57], [0.50, 3.10, 20.20]],
                [[93.16, 38.43, -9.34], [-44.30, -18.57, 6.33], [37.80, -85.02, 15.07], [11.31, -22.63, 5.57], [4.61, 16.17, 79.63], [0.05, 5.49, 21.92]],
                [[51.16, 20.26, -3.58], [-42.36, -18.06, 5.86], [47.96, -105.53, 21.36], [5.99, -12.80, 6.16], [4.45, 22.74, 102.34], [-1.40, 5.62, 13.90]],
                [[8.43, 2.30, 2.37], [-42.52, -17.82, 5.26], [49.71, -110.52, 27.39], [0.24, -1.20, 4.10], [1.79, 27.40, 107.31], [-2.40, 1.97, -2.17]],
                [[-33.87, -15.38, 6.95], [-39.81, -16.66, 4.39], [48.45, -107.98, 29.59], [-5.06, 13.34, -1.82], [-0.35, 26.73, 98.10], [-0.90, -4.12, -19.76]],
                [[-71.20, -31.01, 11.15], [-38.71, -16.87, 4.70], [40.03, -85.11, 24.22], [-15.45, 35.80, -10.05], [-0.15, 19.57, 69.02], [0.48, -8.03, -29.10]],
                [[-111.26, -49.16, 16.38], [-39.85, -17.77, 5.34], [16.93, -35.23, 9.08], [-18.25, 38.75, -11.22], [0.63, 10.58, 39.89], [0.21, -8.76, -30.38]],
                [[-150.90, -66.56, 21.84], [-39.43, -17.02, 5.58], [3.46, -7.44, 1.72], [-8.70, 16.83, -3.49], [0.29, 2.05, 8.29], [-0.90, -8.30, -32.83]]]),
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
                [[54.13, 48.48, 12.30], [-7.00, -7.20, -1.74], [2.81, -2.80, 0.27], [8.05, -8.69, 0.83], [-0.56, -0.25, 3.28], [-1.03, -0.38, 5.83]],
                [[45.64, 40.32, 10.27], [-9.99, -9.12, -2.31], [9.36, -10.47, 0.88], [5.04, -6.65, 0.39], [-1.44, -0.57, 8.47], [-0.71, -0.27, 4.55]],
                [[34.10, 30.30, 7.68], [-11.90, -9.54, -2.50], [12.45, -15.79, 0.99], [1.92, -4.02, -0.41], [-1.95, -0.77, 12.19], [-0.23, -0.37, 2.65]],
                [[21.85, 21.25, 5.27], [-12.81, -9.14, -2.63], [13.20, -18.51, 0.06], [-0.36, -1.77, -0.92], [-1.90, -1.31, 13.78], [0.21, -0.36, 0.38]],
                [[8.49, 12.04, 2.41], [-14.44, -8.62, -2.68], [11.66, -19.26, -0.86], [-2.78, 1.37, -0.56], [-1.51, -1.48, 12.86], [0.47, 0.09, -1.85]],
                [[-6.98, 4.11, -0.07], [-14.87, -7.05, -2.20], [7.56, -15.62, -1.03], [-4.79, 7.26, -0.08], [-0.95, -1.12, 10.03], [0.68, 0.11, -3.53]],
                [[-21.21, -2.10, -2.00], [-16.58, -6.58, -1.81], [2.15, -5.14, -1.02], [-3.60, 6.99, 0.36], [-0.15, -1.23, 5.87], [0.48, 0.38, -3.79]],
                [[-40.17, -8.90, -3.57], [-21.33, -7.01, -1.33], [0.89, -2.66, -0.20], [1.07, -2.03, 1.29], [-0.09, -0.22, 2.55], [-0.35, 1.65, -2.85]]]),
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
                [[177.93, 203.27, -6.98], [-39.04, -53.00, 0.83], [30.11, -22.17, 0.29], [51.78, -52.34, -1.70], [0.03, 0.42, 28.30], [1.06, -1.92, 120.97]],
                [[133.76, 151.30, -6.66], [-49.17, -50.77, -0.19], [66.66, -64.55, -1.68], [21.33, -32.41, -2.24], [1.30, -1.69, 116.41], [1.47, -2.29, 55.25]],
                [[79.63, 102.22, -7.39], [-55.84, -46.38, -0.19], [71.68, -86.28, -4.21], [-0.40, -15.45, 2.49], [3.00, -4.17, 136.47], [-0.79, 2.99, 6.58]],
                [[22.26, 58.63, -7.05], [-62.26, -42.93, 1.27], [65.94, -95.54, 3.23], [-9.21, -3.67, 4.74], [-0.25, 4.21, 129.76], [-1.86, 5.07, -12.37]],
                [[-44.89, 16.76, -4.76], [-67.15, -38.10, 1.57], [52.92, -93.08, 5.01], [-16.55, 14.33, 1.48], [-0.59, 5.64, 111.17], [-0.96, 1.25, -22.72]],
                [[-111.81, -17.62, -3.87], [-64.44, -31.46, 0.87], [33.03, -67.48, 6.20], [-19.80, 30.23, -1.80], [-2.14, 6.72, 84.52], [0.49, -1.71, -32.39]],
                [[-173.66, -46.27, -3.03], [-63.30, -25.43, 1.60], [13.32, -33.05, 1.69], [-14.90, 28.31, -1.76], [0.19, 2.48, 46.94], [0.94, -1.85, -35.99]],
                [[-238.18, -68.39, -0.67], [-65.62, -18.77, 3.11], [3.23, -10.86, 2.68], [-5.27, 16.06, 3.74], [-0.27, 3.02, 12.55], [-1.86, 2.93, -32.79]]]),
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
                [[26.51, 16.41, 4.08], [-9.57, -8.90, -0.96], [4.88, -5.26, 0.16], [1.94, -3.51, 0.19], [-0.47, -0.23, 6.83], [-0.95, -0.39, 8.28]],
                [[16.73, 8.32, 2.96], [-9.96, -7.26, -1.27], [6.35, -8.72, -0.01], [1.01, -3.40, -0.53], [-1.09, -0.82, 13.21], [-0.30, -0.78, 4.48]],
                [[6.67, 1.90, 1.56], [-10.64, -6.01, -1.39], [6.93, -12.06, -0.88], [-0.36, -1.94, -0.43], [-1.08, -1.78, 15.91], [-0.04, -0.39, 1.12]],
                [[-4.52, -3.64, 0.18], [-10.61, -4.65, -1.28], [5.60, -12.53, -0.85], [-1.88, 1.30, -0.10], [-1.18, -1.56, 15.38], [0.03, -0.00, -1.59]],
                [[-14.49, -7.47, -1.01], [-8.85, -2.84, -1.10], [3.25, -9.72, -1.06], [-2.15, 3.53, 0.17], [-1.04, -1.75, 12.88], [0.25, 0.35, -3.93]],
                [[-22.12, -9.45, -1.99], [-7.82, -1.62, -0.91], [1.25, -5.66, -0.61], [-1.42, 3.63, 0.34], [-0.71, -1.00, 7.89], [0.41, 0.62, -4.91]],
                [[-30.11, -10.70, -2.82], [-8.26, -1.29, -0.82], [0.42, -2.47, -0.38], [-0.55, 2.36, 0.22], [-0.22, -0.51, 3.06], [0.30, 0.34, -3.06]],
                [[-38.65, -12.04, -3.63], [-8.81, -1.38, -0.79], [0.17, -1.00, -0.16], [0.04, 0.58, 0.21], [-0.12, -0.33, 1.88], [-0.09, 0.02, 0.71]]]),
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
            'Number of elements along dome': 8,
            'Number of elements along neck': 4,
            'Number of elements around': 8,
            'Number of elements through wall': 1,
            'Wall thickness': 1.5,
            'Ureter position around': 0.67,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine Number of elements around': 4,
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
            'Refine Number of elements around',
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
            'Refine Number of elements around',
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

        elementsCountAlongBladder = elementsCountAlongDome + elementsCountAlongNeck

        fm = region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()
        mesh = fm.findMeshByDimension(3)

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

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
        sx_dome_group = [sx_dome, sd1_dome, sd2_dome, sd12_dome, sd3_dome, sd13_dome]
        sx_neck_group = [sx_neck, sd1_neck, sd2_neck, sd12_neck, sd3_neck, sd13_neck]

        xDome, d1Dome, d2Dome = findNodesAlongBladderDome(sx_dome_group, elementsCountAround, elementsCountAlongDome)

        xNeck, d1Neck, d2Neck, px_transit = findNodesAlongBladderNeck(sx_dome_group, sx_neck_group, domeSegmentLength, neckSegmentLength, elementsCountAround, \
                                                                             elementsCountAlongNeck, xDome, d2Dome)

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

        # Smoothing d2 from apex to down the neck
        d2Raw = []
        for n1 in range(elementsCountAround):
            xAlong = []
            d2Along = []
            for n2 in range(elementsCountAlongBladder):
                v1 = xSampledAll[n2][n1]
                v2 = xSampledAll[n2 + 1][n1]
                d2 = findDerivativeBetweenPoints(v1, v2)
                xAlong.append(v1)
                d2Along.append(d2)
            xAlong.append(xSampledAll[-1][n1])
            d2Along.append(xSampledAll[-1][n1])
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
            d2Raw.append(d2Smoothed)

        # Rearrange d2
        d2SampledAll = []
        for n2 in range(elementsCountAlongBladder + 1):
            d2Around = []
            for n1 in range(elementsCountAround):
                d2 = d2Raw[n1][n2]
                d2Around.append(d2)
            d2SampledAll.append(d2Around)

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
                d3 = d3UnitOuter[n2-1][n1]
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
        d2ApexInner = vector.setMagnitude(sd2_dome[0], mag)

        d1ApexInner = vector.crossproduct3(sd1_dome[0], d2ApexInner)
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

        xFlat = d1Flat = d2Flat = []
        xOrgan = d1Organ = d2Organ = []

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

        # Create nodes and elements
        nodeIdentifier, elementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
            region, xFinal, d1Final, d2Final, d3Final, xFlat, d1Flat, d2Flat, xOrgan, d1Organ, d2Organ, None,
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
                xTrackSurface.append(xList[idx])
                d1TrackSurface.append(d1List[idx])
                d2TrackSurface.append(d2List[idx])
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
        markerList.append({"group": leftUreterGroup, "elementId": idx2, "xi": xi2})
        idx3 = elementsCountAlongDome * elementsCountAround * elementsCountThroughWall + elementsCountAround - ureterElementPositionAround
        xi3 = [1 - ureter1Position.xi1, 0.0, 0.0]
        markerList.append({"group": rightUreterGroup, "elementId": idx3, "xi": xi3})

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
        refineElementsCountAround = options['Refine number of elements along']
        refineElementsCountAlong = options['Refine Number of elements around']
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

    # Find Apex nodes d2
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
        for n2 in range(len(xEllipses_dome) - 1):
            v1 = xEllipses_dome[n2][n1]
            v2 = xEllipses_dome[n2 + 1][n1]
            d2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v1)
            d2Along.append(d2)
        xAlong.append(xEllipses_dome[-1][n1])
        d2Along.append(d2)
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

def findNodesAlongBladderNeck(sx_dome_group, sx_neck_group, domeSegmentLength, neckSegmentLength, elementsCountAround, elementsCountAlongNeck, xDome, d2Dome):

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

    # d2Around = []
    # for n1 in range(elementsCountAround):
    #     v1 = xDome[-1][n1]
    #     v2 = px_transit[n1]
    #     d2 = findDerivativeBetweenPoints(v1, v2)
    #     # xAlong.append(v1)
    #     d2Around.append(d2)
    # d2Dome = d2Dome[:-1] + [d2Around]

    # Create ellipses along neck around the central path
    sx_neck_new = [xTransition] + sx_neck_group[0][e:]
    sd2_neck_new = [d2Transition] + sx_neck_group[2][e:]
    sd3_neck_new = [d3Transition] + sx_neck_group[4][e:]
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
        for n2 in range(len(xEllipses_neck) - 1):
            v1 = xEllipses_neck[n2][n1]
            v2 = xEllipses_neck[n2 + 1][n1]
            d2 = findDerivativeBetweenPoints(v1, v2)
            xAlong.append(v1)
            d2Along.append(d2)
        xAlong.append(xEllipses_neck[-1][n1])
        d2Along.append(d2)
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along)
        d2Raw.append(d2Along)

    # Rearrange d2
    d2Ellipses_neck_new = []
    for n2 in range(len(xEllipses_neck)):
        d2Around = []
        for n1 in range(elementsCountAround):
            d2 = d2Raw[n1][n2]
            d2Around.append(d2)
        d2Ellipses_neck_new.append(d2Around)

    # Spread out elements along neck of the bladder
    xRawNeck = []
    d2RawNeck = []
    for n1 in range(elementsCountAround):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipses_neck)):
            xAlong.append(xEllipses_neck[n2][n1])
            d2Along.append(d2Ellipses_neck_new[n2][n1])
        xSampledAlongNeck, d2SampledAlongNeck = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                                elementsCountAlongNeck, addLengthStart=addingLength,
                                                                                arcLengthDerivatives=True)[0:2]
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlongNeck, d2SampledAlongNeck)
        xRawNeck.append(xSampledAlongNeck)
        d2RawNeck.append(d2Smoothed)

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
