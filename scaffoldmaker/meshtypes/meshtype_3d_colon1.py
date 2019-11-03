"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import MeshType_3d_colonsegment1, ColonSegmentTubeMeshInnerPoints, getTeniaColi, createFlatAndTextureCoordinatesTeniaColi, createNodesAndElementsTeniaColi
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
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
                [ [   0.0,   0.0, 0.0 ], [ -50.7,  178.2, 0.0 ], [ -24.0,  -6.0, -12.0 ], [ -14.0,  -1.0, -12.0 ] ],
                [ [ -47.4, 188.6, 0.0 ], [ -19.3,  177.1, 0.0 ], [ -22.0,  -4.0,  -8.0 ], [  -4.0,  19.0,  22.0 ] ],
                [ [  -4.4, 396.5, 0.0 ], [ 206.0,   40.1, 0.0 ], [ -10.0,  20.0,   8.0 ], [  -6.0,   0.0,  51.0 ] ],
                [ [ 130.0, 384.1, 0.0 ], [ 130.8,  -40.5, 0.0 ], [  -5.0,   4.0,  29.0 ], [   0.0,   1.0,  24.0 ] ],
                [ [ 279.4, 383.0, 0.0 ], [ 118.0,   48.7, 0.0 ], [  -2.0,  10.0,  22.0 ], [   5.0,  25.0, -20.0 ] ],
                [ [ 443.9, 390.8, 0.0 ], [ 111.3,  -97.0, 0.0 ], [  10.0,  17.0,   6.0 ], [   1.0,  -6.0, -35.0 ] ],
                [ [ 475.2, 168.0, 0.0 ], [  -0.8, -112.4, 0.0 ], [  20.0,   0.0, -20.0 ], [  15.0,  -1.0, -10.0 ] ],
                [ [ 432.6, -32.3, 0.0 ], [ -90.5,  -59.0, 0.0 ], [   6.0,  -9.0, -14.0 ], [   8.0, -11.0, -13.0 ] ],
                [ [ 272.4,   7.5, 0.0 ], [ -79.0,   47.4, 0.0 ], [   1.0, -11.0, -18.0 ], [   4.0, -12.0, -12.0 ] ] ] )
            } ),
        'Human 2' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 8
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   0.0,   0.0,    0.0 ], [ -34.7,  104.1, -18.1 ], [ -24.0,  -6.0, -12.0 ], [ -14.0,  -1.0, -12.0 ] ],
                [ [ -34.5,  114.0, -18.1 ], [   1.2,   86.6,  -3.4 ], [ -22.0,  -4.0,  -8.0 ], [  -4.0,  19.0,  22.0 ] ],
                [ [ -19.1,  218.5,   5.5 ], [  78.7,   -7.1,  94.5 ], [ -10.0,  20.0,   8.0 ], [  -6.0,   0.0,  51.0 ] ],
                [ [  82.5,  189.1,  94.2 ], [  84.5,    7.1,  71.6 ], [  -5.0,   4.0,  29.0 ], [   0.0,   1.0,  24.0 ] ],
                [ [ 226.6,  218.7,  85.7 ], [  95.0,   91.3, -58.5 ], [  -2.0,  10.0,  22.0 ], [   5.0,  25.0, -20.0 ] ],
                [ [ 325.5,  381.7, -57.9 ], [ 229.2,  -66.7, -20.4 ], [  10.0,  17.0,   6.0 ], [   1.0,  -6.0, -35.0 ] ],
                [ [ 354.0,  105.3, -24.4 ], [  -6.3, -143.7,  20.3 ], [  20.0,   0.0, -20.0 ], [  15.0,  -1.0, -10.0 ] ],
                [ [ 296.5, -121.2,  -0.6 ], [ -90.5,  -59.0,   0.0 ], [   6.0,  -9.0, -14.0 ], [   8.0, -11.0, -13.0 ] ],
                [ [ 169.8,  -73.4, -33.5 ], [ -72.2,   43.4, -27.4 ], [   1.0, -11.0, -18.0 ], [   4.0, -12.0, -12.0 ] ] ] )
            } ),
        'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 7
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   0.0,   0.0,  0.0 ], [  6.0, 12.0,  -2.0 ], [ 2.0,  1.0,  2.0 ], [ 6.0, 0.0, 3.0 ] ],
                [ [  -2.0,  11.0, -3.0 ], [ -8.0,  4.0,   9.0 ], [ 2.0,  2.0,  1.0 ], [ 0.0, 1.0, 2.0 ] ],
                [ [  -3.0,   2.0,  3.0 ], [ -4.0, -8.0,   0.0 ], [ 2.0, -1.0,  2.0 ], [ 1.0, 0.0, 2.0 ] ],
                [ [ -11.0,  -3.0, -4.0 ], [ -8.0, -3.0,  -7.0 ], [ 1.0, -2.0,  1.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -16.0,  -4.0,  0.0 ], [  4.0, -3.0,  14.0 ], [ 1.0, -3.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -7.0,  -8.0,  0.0 ], [  5.0, -1.0, -14.0 ], [ 0.0, -3.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -1.0,  -6.0, -1.0 ], [  2.0, -2.0,   9.0 ], [ 1.0, -3.0, -1.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  -2.0, -14.0,  5.0 ], [ -2.0, -4.0,   2.0 ], [ 1.0, -2.0, -2.0 ], [ 0.0, 0.0, 0.5 ] ] ] )
            } ),
        'Mouse 2' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 4
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   0.0,  0.0,   0.0 ], [  0.0,  0.0,  13.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [   0.0,  0.0,  13.0 ], [  0.0,  2.0,  28.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -14.0, -2.0,  13.0 ], [  0.0, -3.0, -19.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -14.0, -1.0, -10.0 ], [  1.0,  1.0, -17.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -14.0,  0.0, -28.0 ], [  0.0,  0.0, -11.0 ], [  0.0, -10.0,  0.0 ], [ 0.0, 0.0, 0.5 ] ] ] )
            } ),
        'Pig 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 43
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [  81.0,  24.0,  51.0 ], [ -15.0,   6.0, -34.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  62.0,  30.0,  18.0 ], [ -30.0,   6.0, -23.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  22.0,  36.0,   5.0 ], [ -55.0,  10.0,  -6.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -20.0,  35.0,   8.0 ], [ -36.0, -21.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -40.0,  00.0,  13.0 ], [   0.0, -42.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -20.0, -35.0,  17.0 ], [  36.0, -21.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  20.0, -35.0,  21.0 ], [  36.0,  21.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  40.0,   0.0,  25.0 ], [   0.0,  42.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  23.0,  39.0,  31.0 ], [ -41.0,  24.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -22.0,  39.0,  35.0 ], [ -41.0, -24.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -45.0,   0.0,  39.0 ], [   0.0, -47.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -22.0, -39.0,  44.0 ], [  41.0, -24.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  23.0, -39.0,  48.0 ], [  41.0,  24.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  45.0,   0.0,  51.0 ], [   3.0,  36.0,   6.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  28.0,  48.0,  55.0 ], [ -51.0,  28.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -27.0,  48.0,  65.0 ], [ -50.0, -29.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -55.0,   0.0,  63.0 ], [   0.0, -58.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -28.0, -48.0,  67.0 ], [  50.0, -29.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  27.0, -48.0,  71.0 ], [  50.0,  29.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  55.0,   0.0,  75.0 ], [   0.0,  58.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  30.0,  42.0,  82.0 ], [ -25.0,  14.0,   6.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -13.0,  49.0,  88.0 ], [ -50.0, -13.0,   0.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -47.0,  15.0,  88.0 ], [ -14.0, -51.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -39.0, -25.0,  95.0 ], [  44.0, -08.0,   8.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -18.0,   3.0, 100.0 ], [   0.0,  36.0,   9.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -04.0,  34.0, 105.0 ], [  41.0,  27.0,   4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  31.0,  30.0, 103.0 ], [  22.0, -14.0,  -2.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  44.0, - 1.0,  96.0 ], [  -6.0, -46.0, -10.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  20.0, -31.0,  85.0 ], [ -27.0, -14.0,  -6.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -16.0, -36.0,  76.0 ], [ -28.0,  15.0,  -6.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -39.0,  -5.0,  70.0 ], [  -2.0,  35.0,  -2.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -25.0,  26.0,  70.0 ], [  18.0,  17.0,  -2.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  14.0,  38.0,  68.0 ], [  53.0, -19.0,  -8.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  36.0,  11.0,  65.0 ], [   2.0, -24.0,  -3.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  19.0, -28.0,  62.0 ], [ -41.0, -37.0,  -8.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -17.0, -35.0,  57.0 ], [ -20.0,  10.0,  -4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -37.0,  -2.0,  51.0 ], [   2.0,  44.0,  -8.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -19.0,  32.0,  46.0 ], [  23.0,  11.0,  -6.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  13.0,  36.0,  44.0 ], [  37.0, -15.0,  -3.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  33.0,   9.0,  40.0 ], [  -4.0, -24.0,  -3.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [  20.0, -19.0,  36.0 ], [ -27.0, -15.0,  -2.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -12.0, -27.0,  31.0 ], [ -30.0,   8.0,  -4.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ -19.0,  15.0,  22.0 ], [  55.0,  45.0, -19.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ],
                [ [ 143.0, -59.0, -30.0 ], [  52.0, -41.0,  -3.0 ], [ 0.0, 0.0, 5.0 ], [ 0.0, 0.0, 0.5 ] ] ] )
            } ),
        'Test Line' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 1
                },
            'meshEdits' : zinc_utils.exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [ -40.0, 10.0, 30.0 ], [ 50.0, 10.0, -30.0 ], [ 8.0, -40.0, 0.0 ], [0.0, 0.0, 0.0 ] ],
                [ [  10.0, 20.0,  0.0 ], [ 50.0, 10.0, -30.0 ], [ 8.0, -40.0, 0.0 ], [0.0, 0.0, 0.0 ] ] ])
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
            'Proximal length': 420.0,
            'Transverse length': 460.0,
            'Distal length': 620.0,
            'Proximal inner radius': 43.5,
            'Proximal tenia coli width': 10.0,
            'Proximal-transverse inner radius': 33.0,
            'Proximal-transverse tenia coli width': 10.0,
            'Transverse-distal inner radius': 29.0,
            'Transverse-distal tenia coli width': 10.0,
            'Distal inner radius': 31.5,
            'Distal tenia coli width': 10.0,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
            }
        if 'Human 2' in parameterSetName:
            options['Proximal length'] = 180.0
            options['Transverse length'] = 620.0
            options['Distal length'] = 700.0
        elif 'Mouse' in parameterSetName:
            options['Number of segments'] = 10
            options['Proximal length'] = 30.0
            options['Transverse length'] = 20.0
            options['Distal length'] = 25.0
            options['Proximal inner radius'] = 1.0
            options['Proximal tenia coli width'] = 0.5
            options['Proximal-transverse inner radius'] = 0.9
            options['Proximal-transverse tenia coli width'] = 0.7
            options['Transverse-distal inner radius'] = 0.7
            options['Transverse-distal tenia coli width'] = 1.0
            options['Distal inner radius'] = 0.7
            options['Distal tenia coli width'] = 1.0
        elif 'Pig' in parameterSetName:
            options['Number of segments'] = 120
            options['Proximal length'] = 500.0
            options['Transverse length'] = 800.0
            options['Distal length'] = 700.0
            options['Proximal inner radius'] = 9.0
            options['Proximal tenia coli width'] = 2.0
            options['Proximal-transverse inner radius'] = 9.0
            options['Proximal-transverse tenia coli width'] = 2.0
            options['Transverse-distal inner radius'] = 9.0
            options['Transverse-distal tenia coli width'] = 2.0
            options['Distal inner radius'] = 9.0
            options['Distal tenia coli width'] = 2.0
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

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        # Central path
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

        # Sample central path
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongSegment*segmentCount)
        sd2 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)[0]

        # Generate variation of radius & tc width along length
        lengthList = [0.0, proximalLength, proximalLength + transverseLength, length]
        innerRadiusList = [proximalInnerRadius, proximalTransverseInnerRadius, transverseDistalInnerRadius, distalInnerRadius]
        innerRadiusSegmentList, dInnerRadiusSegmentList = interp.sampleParameterAlongLine(lengthList, innerRadiusList, segmentCount)

        tcWidthList = [proximalTCWidth, proximalTransverseTCWidth, transverseDistalTCWidth, distalTCWidth]
        tcWidthSegmentList, dTCWidthSegmentList = interp.sampleParameterAlongLine(lengthList, tcWidthList, segmentCount)

        xExtrude = []
        d1Extrude = []
        d2Extrude = []
        d3UnitExtrude = []

        # Create object
        colonSegmentTubeMeshInnerPoints = ColonSegmentTubeMeshInnerPoints(
            region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
            tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
            segmentLength, wallThickness, cornerInnerRadiusFactor, haustrumInnerRadiusFactor,
            innerRadiusSegmentList, dInnerRadiusSegmentList, tcWidthSegmentList, dTCWidthSegmentList)

        for nSegment in range(segmentCount):
            # Create inner points
            xInner, d1Inner, d2Inner, transitElementList, segmentAxis, annotationGroups, annotationArray = \
                colonSegmentTubeMeshInnerPoints.getColonSegmentTubeMeshInnerPoints(nSegment)

            # Warp segment points
            xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = tubemesh.warpSegmentPoints(
                xInner, d1Inner, d2Inner, segmentAxis, segmentLength, sx, sd1, sd2,
                elementsCountAround, elementsCountAlongSegment, nSegment)

            # Store points along length
            xExtrude = xExtrude + (xWarpedList if nSegment == 0 else xWarpedList[elementsCountAround:])
            d1Extrude = d1Extrude + (d1WarpedList if nSegment == 0 else d1WarpedList[elementsCountAround:])

            # Smooth d2 for nodes between segments and recalculate d3
            if nSegment == 0:
                d2Extrude = d2Extrude + (d2WarpedList[:-elementsCountAround])
                d3UnitExtrude = d3UnitExtrude + (d3WarpedUnitList[:-elementsCountAround])
            else:
                xSecondFace = xWarpedList[elementsCountAround:elementsCountAround*2]
                d1SecondFace = d1WarpedList[elementsCountAround:elementsCountAround*2]
                d2SecondFace = d2WarpedList[elementsCountAround:elementsCountAround*2]
                for n1 in range(elementsCountAround):
                    nx = [xLastTwoFaces[n1], xLastTwoFaces[n1 + elementsCountAround], xSecondFace[n1]]
                    nd2 = [d2LastTwoFaces[n1], d2LastTwoFaces[n1 + elementsCountAround], d2SecondFace[n1]]
                    d2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)[1]
                    d2Extrude.append(d2)
                    d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1LastTwoFaces[n1 + elementsCountAround]), vector.normalise(d2)))
                    d3UnitExtrude.append(d3Unit)
                d2Extrude = d2Extrude + (d2WarpedList[elementsCountAround:-elementsCountAround] if nSegment < segmentCount - 1 else d2WarpedList[elementsCountAround:])
                d3UnitExtrude = d3UnitExtrude + (d3WarpedUnitList[elementsCountAround:-elementsCountAround] if nSegment < segmentCount - 1 else d3WarpedUnitList[elementsCountAround:])
            xLastTwoFaces = xWarpedList[-elementsCountAround*2:]
            d1LastTwoFaces = d1WarpedList[-elementsCountAround*2:]
            d2LastTwoFaces = d2WarpedList[-elementsCountAround*2:]

        contractedWallThicknessList = colonSegmentTubeMeshInnerPoints.getContractedWallThicknessList()

        # Create coordinates and derivatives
        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xExtrude, d1Extrude,
            d2Extrude, d3UnitExtrude, sx, contractedWallThicknessList,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall, transitElementList)

        relaxedLengthList, xiList = colonSegmentTubeMeshInnerPoints.getRelaxedLengthAndXiList()

        if tcThickness > 0:
            tubeTCWidthList = colonSegmentTubeMeshInnerPoints.getTubeTCWidthList()
            xList, d1List, d2List, d3List, annotationGroups, annotationArray = getTeniaColi(
                region, xList, d1List, d2List, d3List, curvatureList, tcCount, elementsCountAroundTC,
                elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                tubeTCWidthList, tcThickness, sx, annotationGroups, annotationArray)

            # Create flat and texture coordinates
            xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = createFlatAndTextureCoordinatesTeniaColi(
                xiList, relaxedLengthList, length, wallThickness, tcCount, tcThickness,
                elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong,
                elementsCountThroughWall, transitElementList)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = createNodesAndElementsTeniaColi(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
                elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                tcCount, annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
                useCubicHermiteThroughWall, useCrossDerivatives)

        else:
            # Create flat and texture coordinates
            xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture = tubemesh.createFlatAndTextureCoordinates(
                xiList, relaxedLengthList, length, wallThickness, elementsCountAround,
                elementsCountAlong, elementsCountThroughWall, transitElementList)

            # Create nodes and elements
            nextNodeIdentifier, nextElementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
                region, xList, d1List, d2List, d3List, xFlat, d1Flat, d2Flat, xTexture, d1Texture, d2Texture,
                elementsCountAround, elementsCountAlong, elementsCountThroughWall,
                annotationGroups, annotationArray, firstNodeIdentifier, firstElementIdentifier,
                useCubicHermiteThroughWall, useCrossDerivatives)

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
