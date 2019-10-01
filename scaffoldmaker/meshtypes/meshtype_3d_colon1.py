"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1, TubeMeshSegmentInnerPoints, getTeniaColi
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import zinc_utils
from opencmiss.zinc.node import Node

class MeshType_3d_colon1(Scaffold_base):
    '''
    Generates a 3-D colon mesh with variable numbers
    of elements around, along the central line, and through wall.
    The colon is created by a function that generates a colon
    segment and uses tubemesh to map the segment along a central
    line profile.
    '''

    centralPathDefaultScaffoldPackages = {
        'Human 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 8
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
                [ [  0.0,  0.0,  0.0 ], [ -1.6,  5.8,  0.0 ], [ -2.4, -0.6, -1.2 ], [ -1.4, -0.1, -1.2 ] ],
                [ [ -1.6,  6.1,  0.0 ], [ -0.8,  6.2,  0.0 ], [ -2.2, -0.4, -0.8 ], [ -0.4,  1.9,  2.2 ] ],
                [ [ -0.3, 12.3,  0.0 ], [  4.8,  1.2,  0.0 ], [ -1.0,  2.0,  0.8 ], [ -0.6,  0.0,  5.1 ] ],
                [ [  3.0, 12.2,  0.0 ], [  4.0, -0.3,  0.0 ], [ -0.5,  0.4,  2.9 ], [  0.0,  0.1,  2.4 ] ],
                [ [  7.3, 12.3,  0.0 ], [  4.6,  1.1,  0.0 ], [ -0.2,  1.0,  2.2 ], [  0.5,  2.5, -2.0 ] ],
                [ [ 12.2, 12.9,  0.0 ], [  5.1, -3.6,  0.0 ], [  1.0,  1.7,  0.6 ], [  0.1, -0.6, -3.5 ] ],
                [ [ 13.9,  6.4,  0.0 ], [  0.2, -4.0,  0.0 ], [  2.0,  0.0, -2.0 ], [  1.5, -0.1, -1.0 ] ],
                [ [ 12.9, -0.3,  0.0 ], [ -4.1, -3.0,  0.0 ], [  0.6, -0.9, -1.4 ], [  0.8, -1.1, -1.3 ] ],
                [ [  7.7,  0.3,  0.0 ], [ -3.0,  1.7,  0.0 ], [  0.1, -1.1, -1.8 ], [  0.4, -1.2, -1.2 ] ] ] )
            } ),
        'Human 2' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 8
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [  0.0,  0.0,  0.0 ], [ -1.6,  5.8,  0.0 ], [ -2.4, -0.6, -1.2 ], [ -1.4, -0.1, -1.2 ] ],
                [ [ -1.6,  6.1,  0.0 ], [ -0.8,  6.2,  0.0 ], [ -2.2, -0.4, -0.8 ], [ -0.4,  1.9,  2.2 ] ],
                [ [ -0.3, 12.3,  0.0 ], [  3.8, -0.4,  2.6 ], [ -1.0,  2.0,  0.8 ], [ -0.6,  0.0,  5.1 ] ],
                [ [  3.6, 10.4,  2.3 ], [  4.4,  3.9,  1.4 ], [ -0.5,  0.4,  2.9 ], [  0.0,  0.1,  2.4 ] ],
                [ [  9.5, 11.7,  2.4 ], [  5.9, -0.6,  2.0 ], [ -0.2,  1.0,  2.2 ], [  0.5,  2.5, -2.0 ] ],
                [ [ 12.9, 17.3, -1.1 ], [  1.8, -1.5, -5.5 ], [  1.0,  1.7,  0.6 ], [  0.1, -0.6, -3.5 ] ],
                [ [ 14.6,  7.8, -1.5 ], [ -0.2, -3.8, -0.1 ], [  2.0,  0.0, -2.0 ], [  1.5, -0.1, -1.0 ] ],
                [ [ 12.9, -0.3,  0.0 ], [ -4.1, -3.0,  0.0 ], [  0.6, -0.9, -1.4 ], [  0.8, -1.1, -1.3 ] ],
                [ [  7.7,  0.3,  0.0 ], [ -3.0,  1.7,  0.0 ], [  0.1, -1.1, -1.8 ], [  0.4, -1.2, -1.2 ] ] ] )
            } ),
        'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 7
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [  0.0,  0.0,  0.0 ], [  0.6,  1.2, -0.2 ], [  0.2,  0.1,  0.2 ], [ 0.6, 0.0, 0.30 ] ],
                [ [ -0.2,  1.1, -0.3 ], [ -0.8,  0.4,  0.9 ], [  0.2,  0.2,  0.1 ], [ 0.0, 0.1, 0.20 ] ],
                [ [ -0.3,  0.2,  0.3 ], [ -0.4, -0.8,  0.0 ], [  0.2, -0.1,  0.2 ], [ 0.1, 0.0, 0.20 ] ],
                [ [ -1.1, -0.3, -0.4 ], [ -0.8, -0.3, -0.7 ], [  0.1, -0.2,  0.1 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.6, -0.4,  0.0 ], [  0.4, -0.3,  1.4 ], [  0.1, -0.3,  0.0 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -0.7, -0.8,  0.0 ], [  0.5, -0.1, -1.4 ], [  0.0, -0.3,  0.0 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -0.1, -0.6, -0.1 ], [  0.2, -0.2,  0.9 ], [  0.1, -0.3, -0.1 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -0.2, -1.4,  0.5 ], [ -0.2, -0.4,  0.2 ], [  0.1, -0.2, -0.2 ], [ 0.0, 0.0, 0.05 ] ] ] )
            } ),
        'Mouse 2' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 4
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [  0.0,  0.0,  0.0 ], [  0.0,  0.0,  1.3 ], [  0.0, -1.0,  0.0 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  0.0,  0.0,  1.3 ], [  0.0,  0.2,  2.8 ], [  0.0, -1.0,  0.0 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.4, -0.2,  1.3 ], [  0.0, -0.3, -1.9 ], [  0.0, -1.0,  0.0 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.4, -0.1, -1.0 ], [  0.1,  0.1, -1.7 ], [  0.0, -1.0,  0.0 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.4,  0.0, -2.8 ], [  0.0,  0.0, -1.1 ], [  0.0, -1.0,  0.0 ], [ 0.0, 0.0, 0.05 ] ] ] )
            } ),
        'Pig 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 43
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [  8.1,  2.4,  5.1 ], [ -1.5,  0.6, -3.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  6.2,  3.0,  1.8 ], [ -3.0,  0.6, -2.3 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  2.2,  3.6,  0.5 ], [ -5.5,  1.0, -0.6 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -2.0,  3.5,  0.8 ], [ -3.6, -2.1,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -4.0,  0.0,  1.3 ], [  0.0, -4.2,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -2.0, -3.5,  1.7 ], [  3.6, -2.1,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  2.0, -3.5,  2.1 ], [  3.6,  2.1,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  4.0,  0.0,  2.5 ], [  0.0,  4.2,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  2.3,  3.9,  3.1 ], [ -4.1,  2.4,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -2.2,  3.9,  3.5 ], [ -4.1, -2.4,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -4.5,  0.0,  3.9 ], [  0.0, -4.7,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -2.2, -3.9,  4.4 ], [  4.1, -2.4,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  2.3, -3.9,  4.8 ], [  4.1,  2.4,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  4.5,  0.0,  5.1 ], [  0.3,  3.6,  0.6 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  2.8,  4.8,  5.5 ], [ -5.1,  2.8,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -2.7,  4.8,  6.5 ], [ -5.0, -2.9,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -5.5,  0.0,  6.3 ], [  0.0, -5.8,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -2.8, -4.8,  6.7 ], [  5.0, -2.9,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  2.7, -4.8,  7.1 ], [  5.0,  2.9,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  5.5,  0.0,  7.5 ], [  0.0,  5.8,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  3.0,  4.2,  8.2 ], [ -2.5,  1.4,  0.6 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.3,  4.9,  8.8 ], [ -5.0, -1.3,  0.0 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -4.7,  1.5,  8.8 ], [ -1.4, -5.1,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -3.9, -2.5,  9.5 ], [  4.4, -0.8,  0.8 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.8,  0.3, 10.0 ], [  0.0,  3.6,  0.9 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -0.4,  3.4, 10.5 ], [  4.1,  2.7,  0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  3.1,  3.0, 10.3 ], [  2.2, -1.4, -0.2 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  4.4, -0.1,  9.6 ], [ -0.6, -4.6, -1.0 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  2.0, -3.1,  8.5 ], [ -2.7, -1.4, -0.6 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.6, -3.6,  7.6 ], [ -2.8,  1.5, -0.6 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -3.9, -0.5,  7.0 ], [ -0.2,  3.5, -0.2 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -2.5,  2.6,  7.0 ], [  1.8,  1.7, -0.2 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  1.4,  3.8,  6.8 ], [  5.3, -1.9, -0.8 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  3.6,  1.1,  6.5 ], [  0.2, -2.4, -0.3 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  1.9, -2.8,  6.2 ], [ -4.1, -3.7, -0.8 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.7, -3.5,  5.7 ], [ -2.0,  1.0, -0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -3.7, -0.2,  5.1 ], [  0.2,  4.4, -0.8 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.9,  3.2,  4.6 ], [  2.3,  1.1, -0.6 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  1.3,  3.6,  4.4 ], [  3.7, -1.5, -0.3 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  3.3,  0.9,  4.0 ], [ -0.4, -2.4, -0.3 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [  2.0, -1.9,  3.6 ], [ -2.7, -1.5, -0.2 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.2, -2.7,  3.1 ], [ -3.0,  0.8, -0.4 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ -1.9,  1.5,  2.2 ], [  5.5,  4.5, -1.9 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ],
                [ [ 14.3, -5.9, -3.0 ], [  5.2, -4.1, -0.3 ], [ 0.0, 0.0, 0.5 ], [ 0.0, 0.0, 0.05 ] ] ] )
            } ),
        'Test Line' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 1
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [ -4.0, 1.0, 3.0 ], [ 5.0, 1.0, -3.0 ], [ 0.8, -4.0, 0.0 ], [0.0, 0.0, 0.0 ] ],
                [ [  1.0, 2.0, 0.0 ], [ 5.0, 1.0, -3.0 ], [ 0.8, -4.0, 0.0 ], [0.0, 0.0, 0.0 ] ] ])
            } )
        }

    @staticmethod
    def getName():
        return '3D Colon 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Human 2',
            'Mouse 1',
            'Mouse 2',
            'Pig 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 2']
        elif 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
        elif 'Mouse 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 2']
        elif 'Pig' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
        if 'Mouse' in parameterSetName:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Mouse 1')
        elif 'Pig' in parameterSetName:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Pig 1')
        else:
            segmentProfileOption = ScaffoldPackage(MeshType_3d_colonsegment1, defaultParameterSetName = 'Human 1')
        options = {
            'Central path' : copy.deepcopy(centralPathOption),
            'Segment profile' : segmentProfileOption,
            'Number of segments': 30,
            'Proximal length': 12.0,
            'Transverse length': 13.0,
            'Distal length': 20.0,
            'Proximal inner radius': 0.8,
            'Proximal tenia coli width': 0.2,
            'Proximal-transverse inner radius': 1.1,
            'Proximal-transverse tenia coli width': 0.2,
            'Transverse-distal inner radius': 1.25,
            'Transverse-distal tenia coli width': 0.2,
            'Distal inner radius': 1.5,
            'Distal tenia coli width': 0.2,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
            }
        if 'Mouse' in parameterSetName:
            options['Number of segments'] = 10
            options['Proximal length'] = 3.0
            options['Transverse length'] = 2.0
            options['Distal length'] = 2.5
            options['Proximal inner radius'] = 0.1
            options['Proximal tenia coli width'] = 0.05
            options['Proximal-transverse inner radius'] = 0.09
            options['Proximal-transverse tenia coli width'] = 0.07
            options['Transverse-distal inner radius'] = 0.07
            options['Transverse-distal tenia coli width'] = 0.1
            options['Distal inner radius'] = 0.07
            options['Distal tenia coli width'] = 0.1
        elif 'Pig' in parameterSetName:
            options['Number of segments'] = 120
            options['Proximal length'] = 50.0
            options['Transverse length'] = 80.0
            options['Distal length'] = 70.0
            options['Proximal inner radius'] = 0.9
            options['Proximal tenia coli width'] = 0.2
            options['Proximal-transverse inner radius'] = 0.9
            options['Proximal-transverse tenia coli width'] = 0.2
            options['Transverse-distal inner radius'] = 0.9
            options['Transverse-distal tenia coli width'] = 0.2
            options['Distal inner radius'] = 0.9
            options['Distal tenia coli width'] = 0.2
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Segment profile',
            'Number of segments',
            'Proximal length',
            'Transverse length',
            'Distal length',
            'Proximal inner radius',
            'Proximal tenia coli width',
            'Proximal-transverse inner radius',
            'Proximal-transverse tenia coli width',
            'Transverse-distal inner radius',
            'Transverse-distal tenia coli width',
            'Distal inner radius',
            'Distal tenia coli width',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall' ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [ MeshType_1d_path1 ]
        if optionName == 'Segment profile':
            return [ MeshType_3d_colonsegment1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
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
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Segment profile':
            if not parameterSetName:
                parameterSetName = scaffoldType.getParameterSetNames()[0]
            return ScaffoldPackage(scaffoldType, defaultParameterSetName = parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if not options['Segment profile'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Segment profile'):
            options['Segment profile'] = cls.getOptionScaffoldPackage('Segment profile', MeshType_3d_colonsegmentteniacoli1)
        for key in [
            'Number of segments',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Proximal length',
            'Transverse length',
            'Distal length',
            'Proximal inner radius',
            'Proximal tenia coli width',
            'Proximal-transverse inner radius',
            'Proximal-transverse tenia coli width',
            'Transverse-distal inner radius',
            'Transverse-distal tenia coli width',
            'Distal inner radius',
            'Distal tenia coli width']:
            if options[key] < 0.0:
                options[key] = 0.0

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """
        centralPath = options['Central path']
        segmentProfile = options['Segment profile']
        segmentCount = options['Number of segments']
        proximalLength = options['Proximal length']
        transverseLength = options['Transverse length']
        distalLength = options['Distal length']
        proximalInnerRadius = options['Proximal inner radius']
        proximalTCWidth = options['Proximal tenia coli width']
        proximalTransverseInnerRadius = options['Proximal-transverse inner radius']
        proximalTransverseTCWidth = options['Proximal-transverse tenia coli width']
        transverseDistalInnerRadius = options['Transverse-distal inner radius']
        transverseDistalTCWidth = options['Transverse-distal tenia coli width']
        distalInnerRadius = options['Distal inner radius']
        distalTCWidth = options['Distal tenia coli width']
        segmentScaffoldType = segmentProfile.getScaffoldType()
        segmentSettings = segmentProfile.getScaffoldSettings()

        elementsCountAroundTC = segmentSettings['Number of elements around tenia coli']
        elementsCountAroundHaustrum = segmentSettings['Number of elements around haustrum']
        cornerInnerRadiusFactor = segmentSettings['Corner inner radius factor']
        haustrumInnerRadiusFactor = segmentSettings['Haustrum inner radius factor']
        segmentLengthEndDerivativeFactor = segmentSettings['Segment length end derivative factor']
        segmentLengthMidDerivativeFactor = segmentSettings['Segment length mid derivative factor']
        tcCount = segmentSettings['Number of tenia coli']
        tcThickness = segmentSettings['Tenia coli thickness']
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*tcCount

        elementsCountAlongSegment = segmentSettings['Number of elements along segment']
        elementsCountThroughWall = segmentSettings['Number of elements through wall']
        startRadius = segmentSettings['Start inner radius']
        startRadiusDerivative = segmentSettings['Start inner radius derivative']
        endRadius = segmentSettings['End inner radius']
        endRadiusDerivative = segmentSettings['End inner radius derivative']
        wallThickness = segmentSettings['Wall thickness']
        useCrossDerivatives = segmentSettings['Use cross derivatives']
        useCubicHermiteThroughWall = not(segmentSettings['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongSegment*segmentCount)

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        # for i in range(len(cx)):
            # print('cx = ', i+1, cx[i])
            # print('cd1 = ', i+1, cd1[i])
            # print('cd2 = ', i+1, cd2[i])
            # print('cd12 = ', i+1, cd12[i])
        del tmpRegion

        # find arclength of colon
        length = 0.0
        elementsCountIn = len(cx) - 1
        sd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections = True,
            magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
        for e in range(elementsCountIn):
            arcLength = interp.getCubicHermiteArcLength(cx[e], sd1[e], cx[e + 1], sd1[e + 1])
            length += arcLength
        segmentLength = length / segmentCount
        # print('Length = ', length)

        # Generate variation of radius & tc width along length
        lengthList = [0.0, proximalLength, proximalLength + transverseLength, length]
        innerRadiusList = [proximalInnerRadius, proximalTransverseInnerRadius, transverseDistalInnerRadius, distalInnerRadius]
        innerRadiusSegmentList, dInnerRadiusSegmentList = interp.sampleParameterAlongCenterLine(innerRadiusList, lengthList, segmentCount)

        tcWidthList = [proximalTCWidth, proximalTransverseTCWidth, transverseDistalTCWidth, distalTCWidth]
        tcWidthSegmentList, dTCWidthSegmentList = interp.sampleParameterAlongCenterLine(tcWidthList, lengthList, segmentCount)

        # Generate inner surface of a colon segment
        tubeMeshSegmentInnerPoints = TubeMeshSegmentInnerPoints(region, elementsCountAroundTC, elementsCountAroundHaustrum,
            elementsCountAlongSegment, tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
            segmentLength, wallThickness, cornerInnerRadiusFactor, haustrumInnerRadiusFactor, innerRadiusSegmentList, dInnerRadiusSegmentList,
            tcWidthSegmentList, dTCWidthSegmentList)

        # Generate tube mesh
        annotationGroups, nextNodeIdentifier, nextElementIdentifier, xList, d1List, d2List, d3List, sx, curvatureAlong, factorList, uList, relaxedLengthList, tubeTCWidthList = tubemesh.generatetubemesh(region,
           elementsCountAround, elementsCountAlongSegment, elementsCountThroughWall, segmentCount, cx, cd1, cd2, cd12, tubeMeshSegmentInnerPoints, wallThickness, segmentLength, useCrossDerivatives, useCubicHermiteThroughWall)

        # Generate tenia coli
        if tcCount > 1:
            annotationGroupsTC, nextNodeIdentifier, nextElementIdentifier = getTeniaColi(region, nextNodeIdentifier, nextElementIdentifier,
                useCrossDerivatives, useCubicHermiteThroughWall, xList, d1List, d2List, d3List, segmentCount, elementsCountAroundTC, elementsCountAroundHaustrum,
                elementsCountAlong, elementsCountThroughWall, wallThickness, tubeTCWidthList, tcThickness, sx, curvatureAlong, factorList,
                length, tcCount, uList, relaxedLengthList)

            annotationGroups += annotationGroupsTC

        return annotationGroups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup for mesh.
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)
        return meshrefinement.getAnnotationGroups()
