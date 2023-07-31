"""
Generates a 3-D stomach mesh along the central line, with variable
numbers of elements around esophagus and duodenum, along and through
wall, with variable radius and thickness along.
"""

from __future__ import division

import copy
import math

from cmlibs.maths.vectorops import cross, sub
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.utils.zinc.finiteelement import get_maximum_element_identifier, get_maximum_node_identifier
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.element import Element,  Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, mergeAnnotationGroups, \
    getAnnotationGroupForTerm, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.esophagus_terms import get_esophagus_term
from scaffoldmaker.annotation.smallintestine_terms import get_smallintestine_term
from scaffoldmaker.annotation.stomach_terms import get_stomach_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import setEftScaleFactorIds, remapEftNodeValueLabel, remapEftNodeValueLabelsVersion
from scaffoldmaker.utils.geometry import sampleEllipsePoints
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, \
    mesh_destroy_elements_and_nodes_by_identifiers, get_nodeset_path_field_parameters


class MeshType_3d_stomach1(Scaffold_base):
    """
    Generates a 3-D stomach mesh with variable numbers of elements around the esophagus and duodenum,
    along the central line, and through wall. The stomach is created using a central path as the longitudinal axis
    of the stomach. D2 of the central path points to the greater curvature of the stomach and magnitude of D2 and D3
    are the radii of the stomach in the respective direction.
    """
    centralPathDefaultScaffoldPackages = {
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 12
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                (1, [[0.49,0.06,0.00], [-0.03,-0.00,-0.00], [0.01,-0.04,0.00], [0.02,-0.14,0.00], [0.00,0.00,0.04], [0.00,0.00,0.13]]),
                (2, [[0.45,0.05,-0.00], [-0.06,-0.01,-0.00], [0.03,-0.17,0.00], [0.02,-0.12,0.00], [0.00,0.00,0.17], [0.00,0.00,0.13]]),
                (3, [[0.38,0.04,0.00], [-0.08,-0.02,0.00], [0.04,-0.27,0.00], [0.01,-0.10,0.00], [0.00,0.00,0.29], [0.00,0.00,0.11]]),
                (4, [[0.28,0.02,0.00], [-0.11,-0.02,0.00], [0.05,-0.36,0.00], [0.00,-0.07,0.00], [0.00,0.00,0.38], [0.00,0.00,0.07]]),
                (5, [[0.16,0.01,0.00], [-0.14,-0.01,0.00], [0.04,-0.40,0.00], [-0.03,-0.02,0.00], [0.00,0.00,0.42], [0.00,0.00,0.03]]),
                (6, [[0.00,0.00,0.00], [-0.20,0.00,0.00], [-0.01,-0.40,0.00], [-0.06,0.01,0.00], [0.00,0.00,0.43], [0.00,0.00,0.01]]),
                (7, [[-0.23,0.02,0.00], [-0.23,0.05,0.00], [-0.09,-0.37,0.00], [-0.08,0.05,0.00], [0.00,0.00,0.44], [0.00,0.00,-0.01]]),
                (8, [[-0.45,0.10,0.00], [-0.20,0.12,0.00], [-0.18,-0.29,0.00], [-0.06,0.11,0.00], [0.00,0.00,0.41], [0.00,0.00,-0.05]]),
                (9, [[-0.61,0.26,0.00], [-0.15,0.20,0.00], [-0.22,-0.16,0.00], [-0.02,0.13,0.00], [0.00,0.00,0.34], [0.00,0.00,-0.09]]),
                (10, [[-0.73,0.50,0.00], [-0.04,0.23,0.00], [-0.21,-0.04,0.00], [0.04,0.09,0.00], [0.00,0.00,0.23], [0.00,0.00,-0.07]]),
                (11, [[-0.71,0.69,0.00], [0.05,0.17,0.00], [-0.14,0.03,0.00], [0.06,0.05,-0.00], [0.00,0.00,0.18], [0.00,0.00,-0.06]]),
                (12, [[-0.64,0.83,0.00], [0.10,0.14,0.00], [-0.08,0.06,0.00], [0.03,0.03,-0.00], [0.00,0.00,0.11], [0.00,0.00,-0.03]]),
                (13, [[-0.51,0.97,0.00], [0.16,0.14,0.00], [-0.08,0.09,0.00], [-0.03,0.03,0.00], [0.00,0.00,0.13], [0.00,0.00,0.07]])
                ]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-5',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-8',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9-10',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '12',
                    'name': get_smallintestine_term('duodenum')[0],
                    'ontId': get_smallintestine_term('duodenum')[1]
                }]
        }),
        'Human 2': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 12
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                (1, [[61.590,-101.050,1152.900], [0.442,-2.863,-3.864], [11.300,-0.230,1.460], [15.701,-3.029,2.832], [-0.850,-7.440,5.420], [-5.756,-22.456,14.946]]),
                (2, [[61.990,-104.310,1148.310], [0.358,-3.656,-5.314], [24.170,-2.490,3.340], [10.039,-1.491,0.928], [-4.920,-25.060,16.910], [-2.384,-12.784,8.034]]),
                (3, [[62.280,-108.340,1142.260], [0.251,-4.603,-7.040], [30.560,-2.990,3.040], [5.720,-0.448,-0.071], [-5.130,-31.610,20.490], [-0.206,-5.098,2.599]]),
                (4, [[62.470,-113.500,1134.220], [0.315,-6.148,-9.733], [35.400,-3.370,3.270], [4.541,-0.488,0.325], [-5.330,-34.800,21.800], [-0.383,-3.012,1.219]]),
                (5, [[62.940,-120.610,1122.780], [0.448,-8.598,-13.993], [39.520,-4.010,3.730], [3.189,-0.587,-0.295], [-5.970,-37.560,22.890], [-0.099,-1.218,0.549]]),
                (6, [[63.320,-130.670,1106.220], [-0.184,-10.489,-16.967], [41.370,-4.520,2.350], [0.495,-2.468,-4.032], [-5.290,-36.560,22.660], [0.257,0.762,0.246]]),
                (7, [[62.540,-141.560,1088.910], [-4.623,-11.401,-18.639], [40.430,-9.060,-4.490], [-4.947,-4.939,-8.801], [-5.480,-36.050,23.410], [0.770,0.606,1.053]]),
                (8, [[53.670,-152.880,1069.820], [-15.955,-9.817,-16.265], [30.765,-14.468,-15.598], [-14.187,-3.408,-7.577], [-3.580,-35.330,24.820], [2.031,4.346,-1.238]]),
                (9, [[32.230,-159.890,1058.560], [-22.295,-4.560,-7.466], [11.809,-15.766,-19.452], [-15.401,0.331,-0.623], [-1.410,-27.160,20.790], [1.197,8.617,-5.083]]),
                (10, [[10.360,-162.050,1054.830], [-21.091,1.023,-0.270], [-0.460,-14.000,-17.230], [-8.753,2.830,4.765], [-1.070,-18.150,14.780], [-1.004,7.491,-3.259]]),
                (11, [[-8.740,-158.280,1057.630], [-14.760,6.591,2.223], [-6.110,-10.230,-10.220], [-2.404,4.122,5.111], [-3.260,-12.000,13.950], [-1.432,6.623,-3.360]]),
                (12, [[-18.830,-150.690,1059.200], [-11.483,11.369,1.057], [-6.410,-5.880,-6.340], [-0.423,3.102,2.580], [-4.200,-5.070,8.950], [-1.380,4.802,-2.401]]),
                (13, [[-30.740,-135.240,1059.320], [-12.232,19.364,-0.810], [-7.020,-4.680,-5.740], [-0.797,-0.702,-1.380], [-6.250,-3.510,10.510], [-2.720,-1.682,5.521]])
                ]),
                
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9-10',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '12',
                    'name': get_smallintestine_term('duodenum')[0],
                    'ontId': get_smallintestine_term('duodenum')[1]
                }]
        }),
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 11
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                (1, [[0.42,0.72,0.00], [0.04,-0.05,0.00], [0.07,0.03,0.00], [0.12,0.04,0.00], [0.00,0.00,0.04], [0.00,0.00,0.12]]),
                (2, [[0.46,0.65,0.00], [0.04,-0.09,0.00], [0.19,0.04,0.00], [0.11,-0.00,0.00], [0.00,0.00,0.16], [0.00,0.00,0.12]]),
                (3, [[0.48,0.53,0.00], [0.01,-0.14,0.00], [0.26,-0.00,0.00], [0.07,-0.07,0.00], [0.00,0.00,0.27], [0.00,0.00,0.10]]),
                (4, [[0.47,0.36,0.00], [-0.04,-0.17,0.00], [0.30,-0.12,0.01], [0.00,-0.13,0.00], [0.00,0.00,0.35], [0.00,0.00,0.06]]),
                (5, [[0.40,0.20,0.00], [-0.10,-0.15,0.00], [0.29,-0.24,0.00], [-0.05,-0.12,-0.01], [0.00,0.00,0.38], [0.00,0.00,0.02]]),
                (6, [[0.26,0.08,0.00], [-0.21,-0.12,0.00], [0.19,-0.34,0.00], [-0.13,-0.08,0.00], [0.00,0.00,0.39], [0.00,0.00,0.00]]),
                (7, [[0.00,0.00,0.00], [-0.29,-0.00,0.00], [0.01,-0.40,0.00], [-0.19,0.01,0.00], [0.00,0.00,0.38], [0.00,0.00,-0.01]]),
                (8, [[-0.29,0.07,0.00], [-0.24,0.14,0.00], [-0.20,-0.32,0.00], [-0.13,0.14,0.01], [0.00,0.00,0.36], [0.01,0.00,-0.03]]),
                (9, [[-0.44,0.25,0.00], [-0.09,0.18,0.00], [-0.26,-0.15,0.01], [0.01,0.14,0.00], [0.01,0.00,0.32], [0.00,0.00,-0.08]]),
                (10, [[-0.48,0.42,0.00], [-0.02,0.16,0.00], [-0.18,-0.05,0.00], [0.08,0.05,0.00], [0.00,0.00,0.21], [0.00,0.00,-0.09]]),
                (11, [[-0.48,0.56,0.00], [-0.03,0.14,0.00], [-0.10,-0.02,0.00], [0.05,-0.02,0.00], [0.00,0.00,0.13], [0.00,0.00,-0.04]]),
                (12, [[-0.53,0.69,0.00], [-0.06,0.11,0.00], [-0.10,-0.06,0.00], [-0.05,-0.06,0.00], [0.00,0.00,0.12], [0.00,0.00,0.02]])
                ]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-6',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7-8',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '10',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11',
                    'name': get_smallintestine_term('duodenum')[0],
                    'ontId': get_smallintestine_term('duodenum')[1]
                }]
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 10
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [                
                (1, [[0.44,0.01,0.00], [-0.01,-0.00,0.00], [0.01,-0.05,0.00], [-0.01,-0.09,0.00], [0.00,0.00,0.08], [0.00,0.00,0.06]]),
                (2, [[0.42,0.01,0.00], [-0.05,-0.01,0.00], [-0.00,-0.15,0.00], [-0.00,-0.10,0.00], [0.00,0.00,0.15], [0.00,0.00,0.08]]),
                (3, [[0.33,0.00,0.00], [-0.11,-0.00,0.00], [0.00,-0.28,0.00], [-0.00,-0.13,0.00], [0.00,0.00,0.28], [0.00,0.00,0.12]]),
                (4, [[0.19,0.00,0.00], [-0.17,-0.01,0.00], [0.01,-0.39,0.00], [-0.02,-0.08,0.00], [0.00,0.00,0.38], [0.00,0.00,0.07]]),
                (5, [[0.00,0.00,0.00], [-0.21,0.00,0.00], [-0.03,-0.43,0.00], [-0.02,-0.02,0.00], [0.00,0.00,0.40], [0.00,0.00,0.01]]),
                (6, [[-0.22,0.01,0.00], [-0.25,0.03,0.00], [-0.03,-0.43,0.00], [-0.05,0.03,0.00], [0.00,0.00,0.40], [0.00,0.00,-0.01]]),
                (7, [[-0.50,0.05,0.00], [-0.26,0.12,0.00], [-0.16,-0.36,0.00], [-0.12,0.14,0.00], [0.00,0.00,0.37], [0.00,0.00,-0.05]]),
                (8, [[-0.70,0.24,0.00], [-0.11,0.27,0.00], [-0.28,-0.16,0.00], [-0.05,0.21,0.00], [0.00,0.00,0.31], [0.00,0.00,-0.09]]),
                (9, [[-0.70,0.52,0.00], [0.12,0.26,0.00], [-0.27,0.05,0.00], [0.10,0.09,0.00], [0.00,0.00,0.18], [0.00,0.00,-0.11]]),
                (10, [[-0.50,0.70,0.00], [0.14,0.19,0.00], [-0.09,0.03,0.00], [0.08,-0.02,0.00], [0.00,0.00,0.09], [0.00,0.00,-0.04]]),
                (11, [[-0.41,0.88,0.00], [0.03,0.16,0.00], [-0.09,0.01,0.00], [-0.08,-0.03,0.00], [0.00,0.00,0.09], [0.00,0.00,0.04]])
                ]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-7',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '8',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '10',
                    'name': get_smallintestine_term('duodenum')[0],
                    'ontId': get_smallintestine_term('duodenum')[1]
                }]
        }),
        'Rat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 12
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                (1, [[0.492,0.623,0.000], [0.003,-0.035,0.000], [0.060,-0.002,0.000], [0.118,0.016,0.000], [0.000,0.000,0.080], [0.000,0.000,0.077]]),
                (2, [[0.490,0.553,0.000], [-0.008,-0.105,0.000], [0.180,-0.006,0.000], [0.122,-0.024,0.000], [0.000,0.000,0.160], [0.000,0.000,0.083]]),
                (3, [[0.475,0.414,0.000], [-0.038,-0.151,0.000], [0.308,-0.071,0.000], [0.077,-0.090,0.000], [0.000,0.000,0.250], [0.000,0.000,0.072]]),
                (4, [[0.409,0.256,0.000], [-0.082,-0.143,0.000], [0.324,-0.191,0.000], [-0.028,-0.116,0.001], [0.000,0.000,0.300], [0.000,0.000,0.040]]),
                (5, [[0.315,0.129,0.000], [-0.120,-0.114,0.000], [0.256,-0.304,0.001], [-0.090,-0.100,0.000], [0.000,0.000,0.330], [0.000,0.000,0.030]]),
                (6, [[0.171,0.034,0.000], [-0.161,-0.066,0.000], [0.144,-0.389,0.000], [-0.125,-0.059,-0.001], [0.000,0.000,0.360], [0.000,0.000,0.015]]),
                (7, [[0.000,0.000,0.000], [-0.200,0.002,0.000], [0.005,-0.420,0.000], [-0.156,0.002,0.000], [0.000,0.000,0.360], [0.000,0.000,-0.009]]),
                (8, [[-0.218,0.048,0.000], [-0.209,0.094,0.000], [-0.173,-0.374,0.000], [-0.152,0.079,0.000], [0.000,0.000,0.340], [0.000,0.000,-0.015]]),
                (9, [[-0.404,0.184,0.000], [-0.146,0.182,0.000], [-0.299,-0.260,0.000], [-0.056,0.122,0.005], [0.000,0.000,0.330], [0.000,0.000,-0.025]]),
                (10, [[-0.496,0.398,0.000], [-0.040,0.210,0.000], [-0.285,-0.130,0.010], [0.080,0.105,-0.001], [0.000,0.000,0.290], [0.000,0.000,-0.112]]),
                (11, [[-0.490,0.587,0.000], [-0.014,0.157,0.000], [-0.152,-0.045,0.000], [0.076,0.039,-0.004], [0.000,0.000,0.120], [0.000,0.000,-0.067]]),
                (12, [[-0.516,0.708,0.000], [-0.031,0.116,0.000], [-0.111,-0.036,0.000], [0.009,-0.000,0.000], [0.000,0.000,0.120], [0.000,0.000,0.015]]),
                (13, [[-0.552,0.820,0.000], [-0.041,0.107,0.000], [-0.132,-0.045,0.000], [-0.050,-0.018,0.000], [0.000,0.000,0.150], [0.000,0.000,0.045]])
                ]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-6',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7-9',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '10',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '12',
                    'name': get_smallintestine_term('duodenum')[0],
                    'ontId': get_smallintestine_term('duodenum')[1]
                }]
        }),
        'Material': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 15
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                (1, [[0.700,0.000,0.000], [-0.020,0.000,0.000], [0.000,-0.050,0.000], [0.000,-0.100,0.000], [0.000,0.000,0.050], [0.000,0.000,0.100]]),
                (2, [[0.680,0.000,0.000], [-0.050,0.000,0.000], [0.000,-0.150,0.000], [0.000,-0.110,0.000], [0.000,0.000,0.150], [0.000,0.000,0.110]]),
                (3, [[0.600,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.290,0.000], [0.000,-0.130,0.000], [0.000,0.000,0.290], [0.000,0.000,0.130]]),
                (4, [[0.490,0.000,0.000], [-0.130,0.000,0.000], [0.000,-0.400,0.000], [0.000,-0.100,0.000], [0.000,0.000,0.400], [0.000,0.000,0.100]]),
                (5, [[0.350,0.000,0.000], [-0.160,0.000,0.000], [0.000,-0.480,0.000], [0.000,-0.050,0.000], [0.000,0.000,0.480], [0.000,0.000,0.050]]),
                (6, [[0.180,0.000,0.000], [-0.170,0.000,0.000], [0.000,-0.500,0.000], [0.000,-0.010,0.000], [0.000,0.000,0.500], [0.000,0.000,0.010]]),
                (7, [[0.000,0.000,0.000], [-0.190,0.000,0.000], [0.000,-0.500,0.000], [0.000,0.000,0.000], [0.000,0.000,0.500], [0.000,0.000,0.000]]),
                (8, [[-0.200,0.000,0.000], [-0.200,0.000,0.000], [0.000,-0.500,0.000], [0.000,0.000,0.000], [0.000,0.000,0.500], [0.000,0.000,0.000]]),
                (9, [[-0.400,0.000,0.000], [-0.200,0.000,0.000], [0.000,-0.500,0.000], [0.000,0.000,0.000], [0.000,0.000,0.500], [0.000,0.000,0.000]]),
                (10, [[-0.600,0.000,0.000], [-0.180,0.000,0.000], [0.000,-0.500,0.000], [0.000,0.000,0.000], [0.000,0.000,0.500], [0.000,0.000,0.000]]),
                (11, [[-0.750,0.000,0.000], [-0.150,0.000,0.000], [0.000,-0.440,0.000], [0.000,0.100,0.000], [0.000,0.000,0.440], [0.000,0.000,-0.100]]),
                (12, [[-0.900,0.000,0.000], [-0.120,0.000,0.000], [0.000,-0.310,0.000], [0.000,0.120,0.000], [0.000,0.000,0.310], [0.000,0.000,-0.120]]),
                (13, [[-1.000,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.200,0.000], [0.000,0.000,0.000], [0.000,0.000,0.200], [0.000,0.000,0.000]]),
                (14, [[-1.100,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.200,0.000], [0.000,0.000,0.000], [0.000,0.000,0.200], [0.000,0.000,0.000]]),
                (15, [[-1.200,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.200,0.000], [0.000,0.000,0.000], [0.000,0.000,0.200], [0.000,0.000,0.000]]),
                (16, [[-1.300,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.200,0.000], [0.000,0.000,0.000], [0.000,0.000,0.200], [0.000,0.000,0.000]])
                ]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-6',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7-9',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '10-11',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '12-13',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '14-15',
                    'name': get_smallintestine_term('duodenum')[0],
                    'ontId': get_smallintestine_term('duodenum')[1]
                }]
        }),
    }

    ostiumDefaultScaffoldPackages = {
        'Human 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 4,
                'Unit scale': 0.0105,
                'Outlet': False,
                'Ostium diameter': 25.0,
                'Ostium length': 15.0,
                'Ostium wall thickness': 5.0,
                'Ostium wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 5.0,
                'Vessel wall thickness': 3.0,
                'Vessel wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Human 2': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # 4, later
                'Unit scale': 0.0105 * 101,
                'Outlet': False,
                'Ostium diameter': 25.0,
                'Ostium length': 15.0,
                'Ostium wall thickness': 5.0,
                'Ostium wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 5.0,
                'Vessel wall thickness': 3.0,
                'Vessel wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Vessel angle 1 degrees': -40.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': -60.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Mouse 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 4,
                'Unit scale': 0.147,
                'Outlet': False,
                'Ostium diameter': 1.5,
                'Ostium length': 1.5,
                'Ostium wall thickness': 0.35,
                'Ostium wall relative thicknesses': [0.75, 0.05, 0.15, 0.05],
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.5,
                'Vessel wall thickness': 0.2,
                'Vessel wall relative thicknesses': [0.75, 0.05, 0.15, 0.05],
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Pig 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 4,
                'Unit scale': 0.0118,
                'Outlet': False,
                'Ostium diameter': 20.0,
                'Ostium length': 10.0,
                'Ostium wall thickness': 5.0,
                'Ostium wall relative thicknesses': [0.47, 0.1, 0.33, 0.1],
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 3.0,
                'Vessel wall thickness': 3.0,
                'Vessel wall relative thicknesses': [0.47, 0.1, 0.33, 0.1],
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Rat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 4,
                'Unit scale': 0.043,
                'Outlet': False,
                'Ostium diameter': 5.0,
                'Ostium length': 5.0,
                'Ostium wall thickness': 0.5,
                'Ostium wall relative thicknesses': [0.65, 0.12, 0.18, 0.05],
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 2.0,
                'Vessel wall thickness': 0.3,
                'Vessel wall relative thicknesses': [0.65, 0.12, 0.18, 0.05],
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Material': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 4,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 0.3,
                'Ostium length': 0.3,
                'Ostium wall thickness': 0.05,
                'Ostium wall relative thicknesses': [0.25, 0.25, 0.25, 0.25],
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.1,
                'Vessel wall thickness': 0.03,
                'Vessel wall relative thicknesses': [0.25, 0.25, 0.25, 0.25],
                'Vessel angle 1 degrees': 0.0,
                'Vessel angle 1 spread degrees': 0.0,
                'Vessel angle 2 degrees': 0.0,
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
    }

    @staticmethod
    def getName():
        return '3D Stomach 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Human 2',
            'Mouse 1',
            'Pig 1',
            'Rat 1',
            'Material']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 2' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 2']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Human 2']
        elif 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Mouse 1']
        elif 'Pig 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Pig 1']
        elif 'Rat 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Rat 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Rat 1']
        elif 'Material' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Material']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Material']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Human 1']

        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements around duodenum': 16,
            'Number of elements between fundus apex and cardia': 3,
            'Number of elements between cardia and duodenum': 6,
            'Target element unit length': 0.1,
            'Number of elements through wall': 4,
            'Wall thickness': 0.0525,
            'Mucosa relative thickness': 0.55,
            'Submucosa relative thickness': 0.15,
            'Circular muscle layer relative thickness': 0.25,
            'Longitudinal muscle layer relative thickness': 0.05,
            'Limiting ridge': False,
            'Gastro-esophageal junction': copy.deepcopy(ostiumOption),
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements surface': 4,
            'Refine number of elements through wall': 1
        }
        if 'Human 2' in parameterSetName:
            options['Number of elements between cardia and duodenum'] = 8
            options['Number of elements through wall'] = 1  # 4 later
            options['Wall thickness'] = 0.0525 * 101
        elif 'Mouse 1' in parameterSetName:
            options['Number of elements around duodenum'] = 16
            options['Number of elements between fundus apex and cardia'] = 5
            options['Number of elements between cardia and duodenum'] = 5
            options['Wall thickness'] = 0.05145
            options['Mucosa relative thickness'] = 0.75
            options['Submucosa relative thickness'] = 0.05
            options['Circular muscle layer relative thickness'] = 0.15
            options['Longitudinal muscle layer relative thickness'] = 0.05
            options['Limiting ridge'] = True
        elif 'Pig 1' in parameterSetName:
            options['Number of elements around duodenum'] = 16
            options['Number of elements between fundus apex and cardia'] = 3
            options['Number of elements between cardia and duodenum'] = 6
            options['Wall thickness'] = 0.059
            options['Mucosa relative thickness'] = 0.47
            options['Submucosa relative thickness'] = 0.1
            options['Circular muscle layer relative thickness'] = 0.33
            options['Longitudinal muscle layer relative thickness'] = 0.1
            options['Limiting ridge'] = False
        elif 'Rat 1' in parameterSetName:
            options['Number of elements around duodenum'] = 16
            options['Number of elements between fundus apex and cardia'] = 5
            options['Number of elements between cardia and duodenum'] = 6
            options['Wall thickness'] = 0.0215
            options['Mucosa relative thickness'] = 0.65
            options['Submucosa relative thickness'] = 0.12
            options['Circular muscle layer relative thickness'] = 0.18
            options['Longitudinal muscle layer relative thickness'] = 0.05
            options['Limiting ridge'] = True
        elif 'Material' in parameterSetName:
            options['Number of elements around duodenum'] = 16
            options['Number of elements between fundus apex and cardia'] = 5
            options['Number of elements between cardia and duodenum'] = 5
            options['Wall thickness'] = 0.05
            options['Mucosa relative thickness'] = 0.25
            options['Submucosa relative thickness'] = 0.25
            options['Circular muscle layer relative thickness'] = 0.25
            options['Longitudinal muscle layer relative thickness'] = 0.25
            options['Limiting ridge'] = False
        cls.updateSubScaffoldOptions(options)

        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements around duodenum',
            'Number of elements between fundus apex and cardia',
            'Number of elements between cardia and duodenum',
            'Target element unit length',
            'Number of elements through wall',
            'Wall thickness',
            'Mucosa relative thickness',
            'Submucosa relative thickness',
            'Circular muscle layer relative thickness',
            'Longitudinal muscle layer relative thickness',
            'Limiting ridge',
            'Gastro-esophageal junction',
            'Use linear through wall',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        if optionName == 'Gastro-esophageal junction':
            return [MeshType_3d_ostium1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        if optionName == 'Gastro-esophageal junction':
            return list(cls.ostiumDefaultScaffoldPackages.keys())
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
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Gastro-esophageal junction':
            if not parameterSetName:
                parameterSetName = list(cls.ostiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.ostiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if not options['Gastro-esophageal junction'].getScaffoldType() in cls.getOptionValidScaffoldTypes(
                'Gastro-esophageal junction'):
            options['Gastro-esophageal junction'] = cls.getOptionScaffoldPackage('Gastro-esophageal junction',
                                                                                MeshType_3d_ostium1)
        if options['Number of elements around duodenum'] < 12:
            options['Number of elements around duodenum'] = 12
        if options['Number of elements around duodenum'] % 4 > 0:
            options['Number of elements around duodenum'] = options['Number of elements around duodenum'] // 4 * 4
        if options['Number of elements between fundus apex and cardia'] < 2:
            options['Number of elements between fundus apex and cardia'] = 2
        if options['Number of elements between cardia and duodenum'] < 4:
            options['Number of elements between cardia and duodenum'] = 4
        if options['Number of elements through wall'] != (1 or 4):
            options['Number of elements through wall'] = 4
        for key in [
            'Refine number of elements surface',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1

        cls.updateSubScaffoldOptions(options)

    @classmethod
    def updateSubScaffoldOptions(cls, options):
        """
        Update ostium sub-scaffold options which depend on parent options.
        """
        ostiumOptions = options['Gastro-esophageal junction']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
        wallThickness = options['Wall thickness'] / ostiumSettings['Unit scale']
        ostiumSettings['Ostium wall thickness'] = wallThickness
        elementsCountThroughWall = options['Number of elements through wall']
        ostiumSettings['Number of elements through wall'] = elementsCountThroughWall
        if elementsCountThroughWall == 1:
            ostiumSettings['Ostium wall relative thicknesses'] = [1.0]
            ostiumSettings['Vessel wall relative thicknesses'] = [1.0]
        else:
            mucosaRelThickness = options['Mucosa relative thickness']
            submucosaRelThickness = options['Submucosa relative thickness']
            circularRelThickness = options['Circular muscle layer relative thickness']
            longRelThickness = options['Longitudinal muscle layer relative thickness']
            relThicknesses = [mucosaRelThickness, submucosaRelThickness, circularRelThickness, longRelThickness]
            ostiumSettings['Ostium wall relative thicknesses'] = relThicknesses
            ostiumSettings['Vessel wall relative thicknesses'] = relThicknesses

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """
        cls.updateSubScaffoldOptions(options)
        geometricCentralPath = options['Central path']
        materialCentralPath = cls.centralPathDefaultScaffoldPackages['Material']
        limitingRidge = options['Limiting ridge']
        elementsCountThroughWall = options['Number of elements through wall']
        elementsAlongFundusApexToCardia = options['Number of elements between fundus apex and cardia']
        elementsAlongCardiaToDuod = options['Number of elements between cardia and duodenum']
        elementsCountAroundEso = 8
        elementsAroundQuarterEso = int(elementsCountAroundEso * 0.25)
        allAnnotationGroups = []

        stomachTermsAlong = [None, 'fundus of stomach', 'body of stomach',
                             'pyloric antrum', 'pyloric canal', 'duodenum']

        # Geometric coordinates
        fm = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fm)

        geometricCentralPath = StomachCentralPath(region, geometricCentralPath, stomachTermsAlong)
        # geometricCentralPath._path_region.writeFile("C:\\Users\\mlin865\\tmp\\km_stomach_path.exf")
        arcLengthOfGroupsAlong = geometricCentralPath.arcLengthOfGroupsAlong

        # Pre-calculate element numbers in each group
        elementsAlongFromBody = elementsAlongCardiaToDuod + elementsAroundQuarterEso - 1
        arcLengthOfGroupsAlongFromBody = arcLengthOfGroupsAlong[0] - arcLengthOfGroupsAlong[1]
        estElementLengthFromBody = arcLengthOfGroupsAlongFromBody / elementsAlongFromBody

        modGroups = [0.0]
        elementsAlongCPFundus = elementsAlongFundusApexToCardia + elementsAroundQuarterEso - 1
        elementCountGroupList = [elementsAlongCPFundus]
        elementsTally = 0
        for i in range(2, len(stomachTermsAlong)):
            numberOfElementsGroup = int(arcLengthOfGroupsAlong[i] // estElementLengthFromBody)
            if numberOfElementsGroup < 1:
                numberOfElementsGroup = 1

            mod = arcLengthOfGroupsAlong[i] - estElementLengthFromBody * numberOfElementsGroup
            modGroups.append(mod)

            elementsTally += numberOfElementsGroup
            elementCountGroupList.append(numberOfElementsGroup)

        excessElements = elementsAlongFromBody - elementsTally
        for i in range(excessElements):
            maxIdx = max(range(len(modGroups)), key=modGroups.__getitem__)
            elementCountGroupList[maxIdx] += 1
            modGroups[maxIdx] = arcLengthOfGroupsAlong[maxIdx] - estElementLengthFromBody * \
                                elementCountGroupList[maxIdx]

        allAnnotationGroups, nextNodeIdentifier, nextElementIdentifier = \
            createStomachMesh3d(region, fm, coordinates, stomachTermsAlong,
                                allAnnotationGroups, elementCountGroupList, centralPath=geometricCentralPath,
                                options=options, nodeIdentifier=1, elementIdentifier=1, splitCoordinates=False,
                                materialCoordinates=False)

        # # Material coordinates
        # stomach_coordinates = find_or_create_field_coordinates(fm, name="stomach coordinates")
        # allAnnotationGroupsMaterial = []
        # tmp_region = region.createRegion()
        # tmp_fm = tmp_region.getFieldmodule()
        # with ChangeManager(tmp_fm):
        #     tmp_stomach_coordinates = find_or_create_field_coordinates(tmp_fm, name="stomach coordinates")
        #
        #     materialCentralPath = StomachCentralPath(tmp_region, materialCentralPath, stomachTermsAlong)
        #     materialCentralPath._path_region.writeFile("C:\\Users\\mlin865\\tmp\\km_stomach_path_material.exf")
        #
        #     allAnnotationGroupsMaterial, nextNodeIdentifier, nextElementIdentifier = \
        #         createStomachMesh3d(tmp_region, tmp_fm, tmp_stomach_coordinates, stomachTermsAlong,
        #                             allAnnotationGroupsMaterial, elementCountGroupList,
        #                             centralPath=materialCentralPath, options=options, nodeIdentifier=1,
        #                             elementIdentifier=1, splitCoordinates=False, materialCoordinates=True)
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
        #     del tmp_stomach_coordinates
        # del tmp_fm
        # del tmp_region

        # # Create markers
        # markerTermNameStomachCoordinatesMap = {
        #     'body-antrum junction along the greater curvature on luminal surface': [0.6020919166990195, -0.45004378032192227, 0.0],
        #     'body-antrum junction along the greater curvature on serosa': [0.6, -0.5, 0.0],
        #     'distal point of lower esophageal sphincter serosa on the greater curvature of stomach': [1.280068326927052, 0.7999733714717442, 4.9153708368810965e-16],
        #     'distal point of lower esophageal sphincter serosa on the lesser curvature of stomach': [1.1200683310311637, 0.8000096111703247, 5.132730495672042e-16],
        #     'esophagogastric junction along the greater curvature on luminal surface': [1.3499430436270714, 0.44878481258293096, -4.212552457001039e-17],
        #     'esophagogastric junction along the greater curvature on serosa': [1.3499935896233386, 0.4987847870339471, -2.4350160282618435e-17],
        #     'esophagogastric junction along the lesser curvature on luminal surface': [1.0489058130975502, 0.4491717442850351, 3.0345621453573164e-16],
        #     'esophagogastric junction along the lesser curvature on serosa': [1.050012637401148, 0.4991433628042418, 2.8296958630895795e-16],
        #     'gastroduodenal junction along the greater curvature on luminal surface': [0.2, -0.15, 0.0],
        #     'gastroduodenal junction along the greater curvature on serosa': [0.2, -0.2, 0.0],
        #     'gastroduodenal junction along the lesser curvature on luminal surface': [0.2, 0.15, 0.0],
        #     'gastroduodenal junction along the lesser curvature on serosa': [0.20, 0.20, 0.00],
        #     'limiting ridge at the greater curvature on the luminal surface' if limitingRidge else
        #     'fundus-body junction along the greater curvature on luminal surface': [1.1997241080276948, -0.4500013598322351, -0.0002446732391805909],
        #     'limiting ridge at the greater curvature on serosa' if limitingRidge else
        #     'fundus-body junction along the greater curvature on serosa': [1.2, -0.5, 0.0]
        # }
        # if elementsCountThroughWall == 4:
        #     markerTermNameStomachCoordinatesCMLMMap = {
        #         'body-antrum junction along the greater curvature on circular-longitudinal muscle interface': [0.6005229791747548, -0.48751094508048054, 0.0],
        #         'esophagogastric junction along the greater curvature on circular-longitudinal muscle interface': [1.349980953124272, 0.4862847934211931, -2.8794001354466424e-17],
        #         'esophagogastric junction along the lesser curvature on circular-longitudinal muscle interface': [1.0497365634804512, 0.4866625412064305, 3.2195156437946623e-16],
        #         'gastroduodenal junction along the greater curvature on circular-longitudinal muscle interface': [0.2, -0.1875, 0.0],
        #         'gastroduodenal junction along the lesser curvature on circular-longitudinal muscle interface': [0.2, 0.1875, 0.0],
        #         'limiting ridge at the greater curvature on the circular-longitudinal muscle interface' if limitingRidge
        #         else 'fundus-body junction along the greater curvature on circular-longitudinal muscle interface': [1.199934138287874, -0.48750032317766967, -6.116839191743296e-05]
        #     }
        #     markerTermNameStomachCoordinatesMap.update(markerTermNameStomachCoordinatesCMLMMap)
        #
        # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        # nodeIdentifier = max(1, getMaximumNodeIdentifier(nodes) + 1)
        #
        # for termName, stomachCoordinatesValues in markerTermNameStomachCoordinatesMap.items():
        #     annotationGroup = findOrCreateAnnotationGroupForTerm(
        #         allAnnotationGroups, region, get_stomach_term(termName), isMarker=True)
        #     annotationGroup.createMarkerNode(nodeIdentifier, stomach_coordinates, stomachCoordinatesValues)
        #     nodeIdentifier += 1

        return allAnnotationGroups, None

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

    # @classmethod
    # def defineFaceAnnotations(cls, region, options, annotationGroups):
    #     """
    #     Add face annotation groups from the highest dimension mesh.
    #     Must have defined faces and added subelements for highest dimension groups.
    #     :param region: Zinc region containing model.
    #     :param options: Dict containing options. See getDefaultOptions().
    #     :param annotationGroups: List of annotation groups for top-level elements.
    #     New face annotation groups are appended to this list.
    #     """
    #
    #     limitingRidge = options['Limiting ridge']
    #     elementsCountThroughWall = options['Number of elements through wall']
    #
    #     stomachGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("stomach"))
    #     dorsalStomachGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("dorsal stomach"))
    #     ventralStomachGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("ventral stomach"))
    #     bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("body of stomach"))
    #     cardiaGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("cardia of stomach"))
    #     duodenumGroup = getAnnotationGroupForTerm(annotationGroups, get_smallintestine_term("duodenum"))
    #     fundusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("fundus of stomach"))
    #     antrumGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pyloric antrum"))
    #     pylorusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pyloric canal"))
    #     esoGroup = getAnnotationGroupForTerm(annotationGroups, get_esophagus_term("esophagus"))
    #     nearLCGroup = getAnnotationGroupForTerm(annotationGroups,
    #                                             ("elements adjacent to lesser curvature", "None"))
    #
    #     # Create new groups
    #     stomachLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                              get_stomach_term("luminal surface of stomach"))
    #     stomachSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                             get_stomach_term("serosa of stomach"))
    #     bodyLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                           get_stomach_term("luminal surface of body of stomach"))
    #     bodySerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                          get_stomach_term("serosa of body of stomach"))
    #     cardiaLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                             get_stomach_term(
    #                                                                 "luminal surface of cardia of stomach"))
    #     cardiaSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                            get_stomach_term("serosa of cardia of stomach"))
    #     duodenumLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                               get_smallintestine_term("luminal surface of duodenum"))
    #     duodenumSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                              get_smallintestine_term("serosa of duodenum"))
    #     esoLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                          get_esophagus_term("luminal surface of esophagus"))
    #     esoSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                         get_esophagus_term("serosa of esophagus"))
    #     fundusLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                             get_stomach_term(
    #                                                                 "luminal surface of fundus of stomach"))
    #     fundusSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                            get_stomach_term("serosa of fundus of stomach"))
    #     antrumLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                             get_stomach_term("luminal surface of pyloric antrum"))
    #     antrumSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                            get_stomach_term("serosa of pyloric antrum"))
    #     pylorusLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                              get_stomach_term("luminal surface of pyloric canal"))
    #     pylorusSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                             get_stomach_term("serosa of pyloric canal"))
    #     gastroduodenalJunctionGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                                      get_stomach_term("gastroduodenal junction"))
    #
    #     fm = region.getFieldmodule()
    #     mesh2d = fm.findMeshByDimension(2)
    #     mesh1d = fm.findMeshByDimension(1)
    #
    #     is_exterior = fm.createFieldIsExterior()
    #     is_exterior_face_outer = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
    #     is_exterior_face_inner = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
    #
    #     is_gastroduod = fm.createFieldAnd(duodenumGroup.getGroup(), pylorusGroup.getGroup())
    #     gastroduodenalJunctionGroup.getMeshGroup(mesh2d).addElementsConditional(is_gastroduod)
    #
    #     is_dorsal = dorsalStomachGroup.getGroup()
    #     is_ventral = ventralStomachGroup.getGroup()
    #     is_curvatures = fm.createFieldAnd(is_dorsal, is_ventral)
    #
    #     is_nearLC = nearLCGroup.getGroup()
    #     is_lesserCurvature = fm.createFieldAnd(is_curvatures, is_nearLC)
    #     lesserCurvatureGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                               get_stomach_term("lesser curvature of stomach"))
    #     lesserCurvatureGroup.getMeshGroup(mesh2d).addElementsConditional(is_lesserCurvature)
    #
    #     is_nearGC = fm.createFieldNot(is_nearLC)
    #     is_greaterCurvature = fm.createFieldAnd(is_curvatures, is_nearGC)
    #     greaterCurvatureGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                                get_stomach_term("greater curvature of stomach"))
    #     greaterCurvatureGroup.getMeshGroup(mesh2d).addElementsConditional(is_greaterCurvature)
    #
    #     if elementsCountThroughWall == 4:
    #         CMLMInterfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #             "circular-longitudinal muscle interface of stomach"))
    #         circularMuscleGroup = getAnnotationGroupForTerm(annotationGroups,
    #                                                         get_stomach_term("circular muscle layer of stomach"))
    #         longitudinalMuscleGroup = \
    #             getAnnotationGroupForTerm(annotationGroups, get_stomach_term("longitudinal muscle layer of stomach"))
    #         is_CM = circularMuscleGroup.getGroup()
    #         is_LM = longitudinalMuscleGroup.getGroup()
    #         is_CMLMInterface = fm.createFieldAnd(is_CM, is_LM)
    #         CMLMInterfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_CMLMInterface)
    #         is_curvatures_CMLM = fm.createFieldAnd(is_curvatures, is_CMLMInterface)
    #         is_greaterCurvature_CMLM = fm.createFieldAnd(is_greaterCurvature, is_CMLMInterface)
    #         is_lesserCurvature_CMLM = fm.createFieldAnd(is_lesserCurvature, is_CMLMInterface)
    #
    #         dorsalStomach_CMLMGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #             "circular-longitudinal muscle interface of dorsal stomach"))
    #         ventralStomach_CMLMGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #             "circular-longitudinal muscle interface of ventral stomach"))
    #         is_dorsal_CMLM = fm.createFieldAnd(is_dorsal, is_CMLMInterface)
    #         dorsalStomach_CMLMGroup.getMeshGroup(mesh2d).addElementsConditional(is_dorsal_CMLM)
    #         is_ventral_CMLM = fm.createFieldAnd(is_ventral, is_CMLMInterface)
    #         ventralStomach_CMLMGroup.getMeshGroup(mesh2d).addElementsConditional(is_ventral_CMLM)
    #
    #         gastroduod_CMLMGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #             "circular-longitudinal muscle interface of gastroduodenal junction"))
    #         is_gastroduod_CMLM = fm.createFieldAnd(is_gastroduod, is_CMLMInterface)
    #         gastroduod_CMLMGroup.getMeshGroup(mesh1d).addElementsConditional(is_gastroduod_CMLM)
    #
    #         bodyCurvaturesCMLMGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #                 "circular-longitudinal muscle interface of body of stomach along the gastric-omentum attachment"))
    #         duodenumCurvaturesCMLMGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                get_smallintestine_term("circular-longitudinal muscle interface of "
    #                                                                        "first segment of the duodenum along the "
    #                                                                        "gastric-omentum attachment"))
    #         esoCurvaturesCMLMGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_esophagus_term(
    #                 "circular-longitudinal muscle interface of esophagus along the cut margin"))
    #         fundusCurvaturesCMLMGroup =\
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #                 "circular-longitudinal muscle interface of fundus of stomach along the greater curvature"))
    #         antrumGreaterCurvatureCMLMGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #                 "circular-longitudinal muscle interface of pyloric antrum along the greater curvature"))
    #         antrumLesserCurvatureCMLMGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #                 "circular-longitudinal muscle interface of pyloric antrum along the lesser curvature"))
    #         pylorusGreaterCurvatureCMLMGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #                 "circular-longitudinal muscle interface of pyloric canal along the greater curvature"))
    #         pylorusLesserCurvatureCMLMGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #                 "circular-longitudinal muscle interface of pyloric canal along the lesser curvature"))
    #
    #         sectionCurvaturesCMLMGroups = [None, bodyCurvaturesCMLMGroup, None, duodenumCurvaturesCMLMGroup,
    #                                        esoCurvaturesCMLMGroup, fundusCurvaturesCMLMGroup,
    #                                        antrumGreaterCurvatureCMLMGroup, pylorusGreaterCurvatureCMLMGroup,
    #                                        antrumLesserCurvatureCMLMGroup, pylorusLesserCurvatureCMLMGroup]
    #
    #     sectionGroups = [stomachGroup, bodyGroup, cardiaGroup, duodenumGroup, esoGroup, fundusGroup, antrumGroup,
    #                      pylorusGroup]
    #     sectionSerosaGroups = [stomachSerosaGroup, bodySerosaGroup, cardiaSerosaGroup, duodenumSerosaGroup,
    #                            esoSerosaGroup, fundusSerosaGroup, antrumSerosaGroup,
    #                            pylorusSerosaGroup]
    #     sectionLuminalGroups = [stomachLuminalGroup, bodyLuminalGroup, cardiaLuminalGroup, duodenumLuminalGroup,
    #                             esoLuminalGroup, fundusLuminalGroup, antrumLuminalGroup,
    #                             pylorusLuminalGroup]
    #
    #     for i in range(len(sectionGroups)):
    #         is_section = sectionGroups[i].getGroup()
    #         is_sectionSerosa = fm.createFieldAnd(is_section, is_exterior_face_outer)
    #         sectionSerosaGroups[i].getMeshGroup(mesh2d).addElementsConditional(is_sectionSerosa)
    #         is_sectionLuminal = fm.createFieldAnd(is_section, is_exterior_face_inner)
    #         sectionLuminalGroups[i].getMeshGroup(mesh2d).addElementsConditional(is_sectionLuminal)
    #
    #         if elementsCountThroughWall == 4:
    #             if sectionGroups[i] is antrumGroup or sectionGroups[i] is pylorusGroup:
    #                 is_sectionGreaterCurvatureCMLM = fm.createFieldAnd(is_section, is_greaterCurvature_CMLM)
    #                 is_sectionLesserCurvatureCMLM = fm.createFieldAnd(is_section, is_lesserCurvature_CMLM)
    #                 if sectionCurvaturesCMLMGroups[i]:
    #                     sectionCurvaturesCMLMGroups[i].getMeshGroup(mesh1d). \
    #                         addElementsConditional(is_sectionGreaterCurvatureCMLM)
    #                     sectionCurvaturesCMLMGroups[i+2].getMeshGroup(mesh1d). \
    #                         addElementsConditional(is_sectionLesserCurvatureCMLM)
    #             else:
    #                 is_sectionCurvaturesCMLM = fm.createFieldAnd(is_section, is_curvatures_CMLM)
    #                 if sectionCurvaturesCMLMGroups[i]:
    #                     sectionCurvaturesCMLMGroups[i].getMeshGroup(mesh1d). \
    #                         addElementsConditional(is_sectionCurvaturesCMLM)
    #
    #     if limitingRidge:
    #         limitingRidgeGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                get_stomach_term("forestomach-glandular stomach junction"))
    #         innerLimitingRidgeGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                get_stomach_term("limiting ridge on luminal surface"))
    #         outerLimitingRidgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                                      get_stomach_term("limiting ridge on serosa"))
    #
    #         is_antrum = antrumGroup.getGroup()
    #         is_body = bodyGroup.getGroup()
    #         is_cardia = cardiaGroup.getGroup()
    #         is_fundus = fundusGroup.getGroup()
    #         is_limitingRidgeBody = fm.createFieldAnd(is_fundus, is_body)
    #         is_limitingRidgeCardia = fm.createFieldAnd(is_body, is_cardia)
    #         is_limitingRidgeAntrum = fm.createFieldAnd(is_antrum, is_cardia)
    #         is_limitingRidgeBodyCardia = fm.createFieldOr(is_limitingRidgeBody, is_limitingRidgeCardia)
    #         is_limitingRidgeAntrumCardia = fm.createFieldOr(is_limitingRidgeAntrum, is_limitingRidgeCardia)
    #         is_limitingRidge = fm.createFieldOr(is_limitingRidgeBodyCardia, is_limitingRidgeAntrumCardia)
    #
    #         if elementsCountThroughWall == 4:
    #             mucosaGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("mucosa of stomach"))
    #             is_mucosa = mucosaGroup.getGroup()
    #             is_bodyMucosa = fm.createFieldAnd(is_body, is_mucosa)
    #             is_antrumMucosa = fm.createFieldAnd(is_antrum, is_mucosa)
    #             is_bodyAntrumMucosa = fm.createFieldOr(is_bodyMucosa, is_antrumMucosa)
    #
    #             is_xi1Interior = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0),
    #                                                fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1))
    #             is_xi1All = fm.createFieldOr(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0),
    #                                          fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1))
    #             is_xi1BodyAntrumMucosaAll = fm.createFieldAnd(is_bodyAntrumMucosa, is_xi1All)
    #             is_limitingRidgeAroundCardia = fm.createFieldAnd(is_xi1BodyAntrumMucosaAll,
    #                                                              fm.createFieldNot(is_xi1Interior))
    #
    #             is_xi2Interior = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0),
    #                                                fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1))
    #             is_xi2ZeroBodyMucosa = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0), is_bodyMucosa)
    #             is_limitingRidgeBodyBoundary = fm.createFieldAnd(is_xi2ZeroBodyMucosa,
    #                                                              fm.createFieldNot(is_xi2Interior))
    #
    #             is_limitingRidgeMucosa = fm.createFieldOr(is_limitingRidgeAroundCardia, is_limitingRidgeBodyBoundary)
    #             is_limitingRidge = fm.createFieldOr(is_limitingRidge, is_limitingRidgeMucosa)
    #
    #         limitingRidgeGroup.getMeshGroup(mesh2d).addElementsConditional(is_limitingRidge)
    #
    #         is_xi3Interior = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0),
    #                                            fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
    #         is_xi3ZeroLimitingRidge = fm.createFieldAnd(is_limitingRidge,
    #                                                     fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
    #         is_xi3OneLimitingRidge = fm.createFieldAnd(is_limitingRidge,
    #                                                    fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
    #
    #         is_limitingRidgeInner = fm.createFieldAnd(is_xi3ZeroLimitingRidge, fm.createFieldNot(is_xi3Interior))
    #         innerLimitingRidgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_limitingRidgeInner)
    #
    #         is_limitingRidgeOuter = fm.createFieldAnd(is_xi3OneLimitingRidge, fm.createFieldNot(is_xi3Interior))
    #         outerLimitingRidgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_limitingRidgeOuter)
    #
    #         if elementsCountThroughWall == 4:
    #             limitingRidge_CMLMGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #                 "limiting ridge on circular-longitudinal muscle interface"))
    #             is_limitingRidgeCMLM = fm.createFieldAnd(is_CMLMInterface, is_limitingRidge)
    #             limitingRidge_CMLMGroup.getMeshGroup(mesh1d).addElementsConditional(is_limitingRidgeCMLM)
    #
    #     annotationGroups.remove(nearLCGroup)

class StomachCentralPath:
    """
    Generates sampled central path for stomach scaffold.
    """
    def __init__(self, region, centralPath, stomachTermsAlong=[None]):
        """
        :param region: Zinc region needed to create path region to define path in.
        :param centralPath: Central path subscaffold from meshtype_1d_path1
        :param stomachTermsAlong: Annotation terms along length of stomach
        """
        # Extract length of each group along stomach from central path
        arcLengthOfGroupsAlong = []
        cxGroups = []
        cd1Groups = []
        cd2Groups = []
        cd3Groups = []
        cd12Groups = []
        cd13Groups = []

        self._path_region = region.createRegion()
        centralPath.generate(self._path_region)
        tmpFieldmodule = self._path_region.getFieldmodule()
        tmpNodes = tmpFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        tmpCoordinates = tmpFieldmodule.findFieldByName('coordinates')

        for termName in stomachTermsAlong:
            tmpGroup = tmpFieldmodule.findFieldByName(termName).castGroup() if termName else None
            tmpNodeset = tmpGroup.getNodesetGroup(tmpNodes) if tmpGroup else tmpNodes

            cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = get_nodeset_path_field_parameters(
                tmpNodeset, tmpCoordinates,
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                 Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                 Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3])

            arcLength = 0.0
            for e in range(len(cxGroup) - 1):
                arcLength += interp.getCubicHermiteArcLength(cxGroup[e], cd1Group[e],
                                                             cxGroup[e + 1], cd1Group[e + 1])
            arcLengthOfGroupsAlong.append(arcLength)
            cxGroups.append(cxGroup)
            cd1Groups.append(cd1Group)
            cd2Groups.append(cd2Group)
            cd3Groups.append(cd3Group)
            cd12Groups.append(cd12Group)
            cd13Groups.append(cd13Group)

            del tmpNodeset
            del tmpGroup

        del tmpCoordinates
        del tmpNodes
        del tmpFieldmodule

        self.arcLengthOfGroupsAlong = arcLengthOfGroupsAlong
        self.cxGroups = cxGroups
        self.cd1Groups = cd1Groups
        self.cd2Groups = cd2Groups
        self.cd3Groups = cd3Groups
        self.cd12Groups = cd12Groups
        self.cd13Groups = cd13Groups


def findD1CurvatureAround(xAround, d1Around, normsAround):
    """
    Rearrange points around so that the group of points starts from left side of the annulus to the greater curvature
    and ends on the right side of the annulus. This put points in consecutive order for calculating curvature.
    The calculated curvatures are then re-arranged such that it starts from the greater curvature and goes to the right
    side of the annulus, followed by the left side of the annulus and closing the loop back at the greater curvature.
    :param xAround: points around a loop joining to the annulus.
    :param d1Around: derivative of points.
    :param normsAround: radial normal at each point.
    :return: curvature in the original order.
    """
    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
    normsLoop = normsAround[int(len(d1Around) * 0.5 + 1):] + normsAround[: int(len(d1Around) * 0.5 + 1)]
    curvature = findCurvatureAlongLine(xLoop, d1Loop, normsLoop)
    # Rearrange to correct order
    d1CurvatureAround = curvature[int(len(xAround) * 0.5):] + curvature[: int(len(xAround) * 0.5):]

    return d1CurvatureAround


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


def findCurvatureAroundLoop(nx, nd, radialVectors):
    """
    Calculate curvature for points lying along a loop.
    :param nx: points on loop
    :param nd: derivative of points on loop
    :param radialVectors: radial direction, assumed normal to curve tangent at point
    :return: curvatures along points on loop
    """
    curvature = []
    for n in range(len(nx)):
        prevIdx = n - 1 if n > 0 else -1
        nextIdx = n + 1 if n < len(nx) - 1 else 0
        kappam = interp.getCubicHermiteCurvature(nx[prevIdx], nd[prevIdx], nx[n], nd[n], radialVectors[n], 1.0)
        kappap = interp.getCubicHermiteCurvature(nx[n], nd[n], nx[nextIdx], nd[nextIdx], radialVectors[n], 0.0)
        curvature.append(0.5 * (kappam + kappap))

    return curvature


def findCurvatureAlongLine(nx, nd, radialVectors):
    """
    Calculate curvature for points lying along a line.
    :param nx: points on line
    :param nd: derivative of points on line
    :param radialVectors: radial direction, assumed normal to curve tangent at point
    :return: curvatures along points on line
    """
    curvature = []
    for n in range(len(nx)):
        if n == 0:
            curvature.append(interp.getCubicHermiteCurvature(nx[n], nd[n], nx[n + 1], nd[n + 1], radialVectors[n], 0.0))
        elif n == len(nx) - 1:
            curvature.append(interp.getCubicHermiteCurvature(nx[n - 1], nd[n - 1], nx[n], nd[n], radialVectors[n], 1.0))
        else:
            curvature.append(0.5 * (
                    interp.getCubicHermiteCurvature(nx[n], nd[n], nx[n + 1], nd[n + 1], radialVectors[n], 0.0) +
                    interp.getCubicHermiteCurvature(nx[n - 1], nd[n - 1], nx[n], nd[n], radialVectors[n], 1.0)))

    return curvature


def createStomachMesh3d(region, fm, coordinates, stomachTermsAlong, allAnnotationGroups,
                        elementCountGroupList, centralPath, options, nodeIdentifier, elementIdentifier,
                        splitCoordinates, materialCoordinates=False):
    """
    Generates a stomach scaffold in the region using a central path and parameter options.
    :param region: Region to create elements in.
    :param fm: Zinc fieldModule to create elements in.
    :param coordinates: Coordinate field to define nodes and elements.
    :param stomachTermsAlong: Annotation terms along length of stomach.
    :param allAnnotationGroups: List of annotation groups.
    :param elementCountGroupList: List contains number of elements in each annotation group along the stomach.
    :param centralPath: Central path through the axis of the stomach scaffold.
    :param options: Parameter options for stomach scaffold.
    :param nodeIdentifier: First node identifier.
    :param elementIdentifier: First element identifier.
    :param splitCoordinates: Create split coordinates if True.
    :param materialCoordinates: Create material coordinates if True.
    :return allAnnotationGroups, elementsCountGroupList
    """
    elementsCountAroundDuod = options['Number of elements around duodenum']
    elementsAlongFundusApexToCardia = options['Number of elements between fundus apex and cardia']
    elementsAlongCardiaToDuod = options['Number of elements between cardia and duodenum']
    targetLength = options['Target element unit length']
    elementsCountThroughWall = options['Number of elements through wall']
    mucosaRelThickness = options['Mucosa relative thickness']
    submucosaRelThickness = options['Submucosa relative thickness']
    circularRelThickness = options['Circular muscle layer relative thickness']
    longitudinalRelThickness = options['Longitudinal muscle layer relative thickness']
    useCrossDerivatives = False
    useCubicHermiteThroughWall = not (options['Use linear through wall'])

    GEJOptions = options['Gastro-esophageal junction']
    GEJSettings = GEJOptions.getScaffoldSettings()
    elementsAlongEsophagus = GEJSettings['Number of elements along']
    elementsThroughEsophagusWall = GEJSettings['Number of elements through wall']
    ostiumDiameter = GEJSettings['Ostium diameter']
    limitingRidge = options['Limiting ridge']
    wallThickness = options['Wall thickness']

    elementsCountAroundEso = 8
    elementsCountAcrossCardia = 1
    cardiaDiameterFactor = 1.4  # scale to ostium diameter
    sf = (cardiaDiameterFactor - 1) * ostiumDiameter * 0.5 * GEJSettings['Unit scale']

    elementsAroundHalfEso = int(elementsCountAroundEso * 0.5)
    elementsAroundQuarterEso = int(elementsCountAroundEso * 0.25)
    elementsAroundHalfDuod = int(elementsCountAroundDuod * 0.5)
    elementsAroundQuarterDuod = int(elementsCountAroundDuod * 0.25)
    zero = [0.0, 0.0, 0.0]

    if materialCoordinates:
        wallThickness = 0.05
        mucosaRelThickness = 0.25
        submucosaRelThickness = 0.25
        circularRelThickness = 0.25
        longitudinalRelThickness = 0.25
        if elementsCountThroughWall == 4:
            relThicknesses = [mucosaRelThickness, submucosaRelThickness, circularRelThickness, longitudinalRelThickness]
        else:
            relThicknesses = [1.0]

        unitScale = GEJSettings['Unit scale']
        ostiumLength = GEJSettings['Ostium length']
        ostiumWallThickness = GEJSettings['Ostium wall thickness']
        ostiumWallRelThicknesses = GEJSettings['Ostium wall relative thicknesses']
        vesselInnerDiameter = GEJSettings['Vessel inner diameter']
        vesselWallThickness = GEJSettings['Vessel wall thickness']
        vesselWallRelThicknesses = GEJSettings['Vessel wall relative thicknesses']
        vesselAngle1 = GEJSettings['Vessel angle 1 degrees']
        vesselAngle2 = GEJSettings['Vessel angle 2 degrees']

        GEJSettings['Unit scale'] = 1.0
        GEJSettings['Ostium diameter'] = 0.3
        GEJSettings['Ostium length'] = 0.3
        GEJSettings['Ostium wall thickness'] = wallThickness
        GEJSettings['Ostium wall relative thicknesses'] = relThicknesses
        GEJSettings['Vessel inner diameter'] = 0.1
        GEJSettings['Vessel wall thickness'] = wallThickness * 0.6
        GEJSettings['Vessel wall relative thicknesses'] = relThicknesses
        sf = (cardiaDiameterFactor - 1) * GEJSettings['Ostium diameter'] * 0.5
        GEJSettings['Vessel angle 1 degrees'] = 0.0
        GEJSettings['Vessel angle 2 degrees'] = 0.0

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

    # Create annotation groups for stomach sections
    arcLengthRatioForGroupsFromFundusApex = []
    arcLengthOfGroupsAlong = centralPath.arcLengthOfGroupsAlong
    stomachCentralPathLength = arcLengthOfGroupsAlong[0]

    for i in range(1, len(stomachTermsAlong)):
        arcLengthRatio = arcLengthOfGroupsAlong[i] / stomachCentralPathLength
        arcLengthRatioForGroupsFromFundusApex.append(arcLengthRatio)

    stomachGroup = AnnotationGroup(region, get_stomach_term("stomach"))
    fundusGroup = AnnotationGroup(region, get_stomach_term("fundus of stomach"))
    bodyGroup = AnnotationGroup(region, get_stomach_term("body of stomach"))
    antrumGroup = AnnotationGroup(region, get_stomach_term("pyloric antrum"))
    pylorusGroup = AnnotationGroup(region, get_stomach_term("pyloric canal"))
    duodenumGroup = AnnotationGroup(region, get_smallintestine_term("duodenum"))

    annotationGroupAlong = [[stomachGroup, fundusGroup],
                            [stomachGroup, bodyGroup],
                            [stomachGroup, antrumGroup],
                            [stomachGroup, pylorusGroup],
                            [stomachGroup, duodenumGroup]]

    longitudinalMuscleGroup = AnnotationGroup(region, get_stomach_term("longitudinal muscle layer of stomach"))
    circularMuscleGroup = AnnotationGroup(region, get_stomach_term("circular muscle layer of stomach"))
    submucosaGroup = AnnotationGroup(region, get_stomach_term("submucosa of stomach"))
    mucosaGroup = AnnotationGroup(region, get_stomach_term("mucosa of stomach"))

    if elementsCountThroughWall == 1:
        annotationGroupsThroughWall = [[]]
    else:
        annotationGroupsThroughWall = [[mucosaGroup],
                                       [submucosaGroup],
                                       [circularMuscleGroup],
                                       [longitudinalMuscleGroup]]

    # Break central path into elements allocation to each group
    cxSections = []
    cd1Sections = []
    cd2Sections = []
    cd3Sections = []
    # +/- deltas to get d2 along by finite difference
    deltaXi = 1.0E-5
    cxPlusSections = []
    cxMinusSections = []
    cd2PlusSections = []
    cd2MinusSections = []
    cd3PlusSections = []
    cd3MinusSections = []

    targetLengthTS = 0.025

    for i in (list(range(1, len(elementCountGroupList))) + [0]):  # start from body, go back to fundus
        cxGroup = centralPath.cxGroups[i + 1]
        cd1Group = centralPath.cd1Groups[i + 1]
        cd2Group = centralPath.cd2Groups[i + 1]
        cd3Group = centralPath.cd3Groups[i + 1]
        cd12Group = centralPath.cd12Groups[i + 1]
        cd13Group = centralPath.cd13Groups[i + 1]

        for n2 in range(len(cxGroup)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxGroup[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd2Group[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd1Group[n2])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, cd3Group[n2])
            nodeIdentifier += 1

        if materialCoordinates and i == len(elementCountGroupList) - 1:
            for n in range(len(cxGroup)):
                cd12Group[n] = zero
                cd13Group[n] = zero

        # Break body into equal sized elements, all others vary smoothly from end derivative of last section
        # Except for fundus which start at zero derivative and ends at start derivative for body
        elementsOutSection = math.ceil(arcLengthRatioForGroupsFromFundusApex[i]/targetLengthTS)
        cxSection, cd1Section, pe, pxi, psf = interp.sampleCubicHermiteCurvesSmooth(
            cxGroup, cd1Group, elementsOutSection,
            derivativeMagnitudeStart=None if (i == 1) else 0.0 if (i == 0) else vector.magnitude(cd1Sections[-1][-1]),
            derivativeMagnitudeEnd=None if (i != 0) else vector.magnitude(cd1Sections[0][0]))
        cd2Section = interp.interpolateSampleCubicHermite(cd2Group, cd12Group, pe, pxi, psf)[0]
        cd3Section = interp.interpolateSampleCubicHermite(cd3Group, cd13Group, pe, pxi, psf)[0]

        pxiPlus = [xi + deltaXi for xi in pxi]
        cxPlusSection = interp.interpolateSampleCubicHermite(cxGroup, cd1Group, pe, pxiPlus, psf)[0]
        cd2PlusSection = interp.interpolateSampleCubicHermite(cd2Group, cd12Group, pe, pxiPlus, psf)[0]
        cd3PlusSection = interp.interpolateSampleCubicHermite(cd3Group, cd13Group, pe, pxiPlus, psf)[0]
        pxiMinus = [xi - deltaXi for xi in pxi]
        cxMinusSection = interp.interpolateSampleCubicHermite(cxGroup, cd1Group, pe, pxiMinus, psf)[0]
        cd2MinusSection = interp.interpolateSampleCubicHermite(cd2Group, cd12Group, pe, pxiMinus, psf)[0]
        cd3MinusSection = interp.interpolateSampleCubicHermite(cd3Group, cd13Group, pe, pxiMinus, psf)[0]

        cxSections.append(cxSection)
        cd1Sections.append(cd1Section)
        cd2Sections.append(cd2Section)
        cd3Sections.append(cd3Section)
        cxPlusSections.append(cxPlusSection)
        cd2PlusSections.append(cd2PlusSection)
        cd3PlusSections.append(cd3PlusSection)
        cxMinusSections.append(cxMinusSection)
        cd2MinusSections.append(cd2MinusSection)
        cd3MinusSections.append(cd3MinusSection)

    # put fundus section first
    for values in [cxSections, cd1Sections, cd2Sections, cd3Sections,
                   cxPlusSections, cd2PlusSections, cd3PlusSections,
                   cxMinusSections, cd2MinusSections, cd3MinusSections]:
        values.insert(0, values.pop())

    nodeStartEndSections = []
    startNode = 0
    for i in range(len(cxSections)):
        endNode = startNode + len(cxSections[i]) - 1
        nodeStartEndSections.append([startNode, endNode])
        startNode = endNode

    # Create ellipses
    cxApex = cxSections[0][0]
    xApex = [cxApex for n1 in range(elementsCountAroundDuod)]
    d1ApexAround = []

    d2Apex = cd2PlusSections[0][0]
    d3Apex = cd3PlusSections[0][0]
    rotAxisApex = vector.normalise(vector.crossproduct3(d3Apex, d2Apex))

    px = sampleEllipsePoints(cxPlusSections[0][0], cd2PlusSections[0][0], cd3PlusSections[0][0],
                             0.0, math.pi * 2.0, elementsCountAroundDuod)[0]
    d2ApexAround = [cross(cross(rotAxisApex, sub(tpx, cxApex)), rotAxisApex) for tpx in px]

    rotAngle = -math.pi * 0.5
    rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxisApex, rotAngle)
    for n in range(len(px)):
        d1ApexAround.append([rotFrame[j][0] * d2ApexAround[n][0] + rotFrame[j][1] * d2ApexAround[n][1] +
                             rotFrame[j][2] * d2ApexAround[n][2] for j in range(3)])

    xEllipseAroundAll = [xApex]
    d1EllipseAroundAll = [d1ApexAround]
    d2EllipseAroundAll = [d2ApexAround]
    d2Curvature = []
    curvature = [1.0 for n in range(elementsCountAroundDuod)]
    d2Curvature.append(curvature)

    for s in range(len(cxSections)):
        for n2 in range(1, len(cxSections[s])):
            px, pd1 = sampleEllipsePoints(cxSections[s][n2], cd2Sections[s][n2], cd3Sections[s][n2],
                                          0.0, math.pi * 2.0, elementsCountAroundDuod)
            px.pop()
            pd1.pop()

            # get d2 from finite difference between plus and minus ellipses. note scale is not right
            pxPlus = sampleEllipsePoints(cxPlusSections[s][n2], cd2PlusSections[s][n2], cd3PlusSections[s][n2],
                                         0.0, math.pi * 2.0, elementsCountAroundDuod)[0]
            pxPlus.pop()

            pxMinus = sampleEllipsePoints(cxMinusSections[s][n2], cd2MinusSections[s][n2], cd3MinusSections[s][n2],
                                         0.0, math.pi * 2.0, elementsCountAroundDuod)[0]
            pxMinus.pop()

            d2Around = [sub(pxPlus[n], pxMinus[n]) for n in range(len(pxPlus))]

            d2CurvatureAround = [0.0 for n in range(len(pd1))]
            xEllipseAroundAll.append(px)
            d1EllipseAroundAll.append(pd1)
            d2EllipseAroundAll.append(d2Around)
            d2Curvature.append(d2CurvatureAround)

            if s == 0 and n2 == len(cxSections[s]) - 1:
                xGEJ = px[elementsAroundHalfDuod]

    # Scale d1 and d2 at apex
    for n in range(len(xEllipseAroundAll[0])):
        d1EllipseAroundAll[0][n] = \
            vector.setMagnitude(d1EllipseAroundAll[0][n],
                                interp.computeCubicHermiteArcLength(xEllipseAroundAll[0][n], d2EllipseAroundAll[0][n],
                                                                    xEllipseAroundAll[1][n], d2EllipseAroundAll[1][n],
                                                                    True))
        d2EllipseAroundAll[0][n] = \
            vector.setMagnitude(d2EllipseAroundAll[0][n],
                                interp.computeCubicHermiteArcLength(xEllipseAroundAll[0][n], d2EllipseAroundAll[0][n],
                                                                    xEllipseAroundAll[1][n], d2EllipseAroundAll[1][n],
                                                                    True))

    # Create track surface
    # Find d2
    d2Raw = []
    for n1 in range(elementsCountAroundDuod):
        xAlong = []
        d2Along = []
        for n2 in range(len(xEllipseAroundAll)):
            xAlong.append(xEllipseAroundAll[n2][n1])
            d2Along.append(d2EllipseAroundAll[n2][n1])
        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xAlong, d2Along, fixAllDirections=True)
        d2Raw.append(d2Smoothed)

    # Rearrange d2
    for n2 in range(len(xEllipseAroundAll)):
        for n1 in range(elementsCountAroundDuod):
            d2EllipseAroundAll[n2][n1] = d2Raw[n1][n2]

    # Copy points on lesser curvature before putting annulus
    xTopCurvature = []
    for n2 in range(len(xEllipseAroundAll)):
        xTopCurvature.append(xEllipseAroundAll[n2][elementsAroundHalfDuod])

    # Create tracksurface
    xTrackSurface = []
    d1TrackSurface = []
    d2TrackSurface = []

    for n2 in range(len(xEllipseAroundAll)):
        for n1 in range(len(xEllipseAroundAll[n2])):
            xTrackSurface.append(xEllipseAroundAll[n2][n1])
            d1TrackSurface.append(d1EllipseAroundAll[n2][n1])
            d2TrackSurface.append(d2EllipseAroundAll[n2][n1])

    trackSurfaceStomach = TrackSurface(elementsCountAroundDuod, len(xEllipseAroundAll) - 1,
                                       xTrackSurface, d1TrackSurface, d2TrackSurface, loop1=True)

    # Visualise track surface
    # tmpRegion = region.createRegion()
    nodeIdentifier, elementIdentifier = trackSurfaceStomach.generateMesh(region) #tmpRegion)

    # nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    # startNodeIdentifier = nodeIdentifier = max(get_maximum_node_identifier(nodes), 0) + 1
    # nodetemplate = nodes.createNodetemplate()
    # nodetemplate.defineField(coordinates)
    # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    # nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    #
    # mesh = fm.findMeshByDimension(2)
    # elementIdentifier = max(get_maximum_element_identifier(mesh), 0) + 1
    # elementtemplate = mesh.createElementtemplate()
    # elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
    # bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    # eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
    # # remove cross derivative terms
    # for n in range(4):
    #     eft.setFunctionNumberOfTerms(n * 4 + 4, 0)
    # elementtemplate.defineField(coordinates, -1, eft)
    #
    # fieldcache = fm.createFieldcache()
    # with ChangeManager(fm):
    #     for n1 in range(len(trackSurfaceStomach.nx)):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         fieldcache.setNode(node)
    #         coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, trackSurfaceStomach.nx[n1])
    #         coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS1, 1, trackSurfaceStomach.nd1[n1])
    #         coordinates.setNodeParameters(fieldcache, -1, Node.VALUE_LABEL_D_DS2, 1, trackSurfaceStomach.nd2[n1])
    #         nodeIdentifier += 1
    #
    #     del n1
    #     nodesCount1 = trackSurfaceStomach.elementsCount1 + 1
    #     for e2 in range(trackSurfaceStomach.elementsCount2):
    #         for e1 in range(trackSurfaceStomach.elementsCount1):
    #             nid1 = startNodeIdentifier + e2 * nodesCount1 + e1
    #             nids = [nid1, nid1 + 1, nid1 + nodesCount1, nid1 + nodesCount1 + 1]
    #             element = mesh.createElement(elementIdentifier, elementtemplate)
    #             element.setNodesByIdentifier(eft, nids)
    #             elementIdentifier += 1

    # tmpRegion.getFieldmodule().defineAllFaces()
    # tmpRegion.writeFile("C:\\Users\\mlin865\\tmp\\km_stomach_tracksurface.exf")
    # del tmpRegion

    # # Set up gastro-esophageal junction with midpoint aligned to fundus-body junction
    # GEJSettings['Number of elements around ostium'] = elementsCountAroundEso
    # GEJPosition = trackSurfaceStomach.findNearestPosition(xGEJ)
    # xCentre, d1Centre, d2Centre = trackSurfaceStomach.evaluateCoordinates(GEJPosition, derivatives=True)
    # axis1 = d1Centre
    #
    # esophagusGroup = AnnotationGroup(region, get_esophagus_term("esophagus"))
    # esophagusMeshGroup = esophagusGroup.getMeshGroup(mesh)
    # abdominalEsoGroup = AnnotationGroup(region, get_esophagus_term("abdominal part of esophagus"))
    # abdominalEsoMeshGroup = abdominalEsoGroup.getMeshGroup(mesh)
    # esophagogastricJunctionGroup = AnnotationGroup(region, get_stomach_term("esophagogastric junction"))
    # esophagogastricJunctionMeshGroup = esophagogastricJunctionGroup.getMeshGroup(mesh)
    # stomachMeshGroup = stomachGroup.getMeshGroup(mesh)
    # allAnnotationGroups += [esophagusGroup, esophagogastricJunctionGroup, abdominalEsoGroup]
    #
    # ostiumWallAnnotationGroups = []
    # if elementsCountThroughWall == 4:
    #     esophagusMucosaGroup = AnnotationGroup(region, get_esophagus_term("esophagus mucosa"))
    #     esophagusSubmucosaGroup = AnnotationGroup(region, get_esophagus_term("submucosa of esophagus"))
    #     esophagusCircularGroup = AnnotationGroup(region, get_esophagus_term("esophagus smooth muscle circular layer"))
    #     esophagusLongitudinalGroup = AnnotationGroup(region,
    #                                                  get_esophagus_term("esophagus smooth muscle longitudinal layer"))
    #
    #     ostiumWallAnnotationGroups = [[esophagusMucosaGroup, mucosaGroup],
    #                                   [esophagusSubmucosaGroup, submucosaGroup],
    #                                   [esophagusCircularGroup, circularMuscleGroup],
    #                                   [esophagusLongitudinalGroup, longitudinalMuscleGroup]]
    #
    #     allAnnotationGroups += [esophagusMucosaGroup, esophagusSubmucosaGroup,
    #                             esophagusCircularGroup, esophagusLongitudinalGroup]
    #
    # nextNodeIdentifier, nextElementIdentifier, (o1_x, o1_d1, o1_d2, o1_d3, o1_NodeId, o1_Positions) = \
    #     generateOstiumMesh(region, GEJSettings, trackSurfaceStomach, GEJPosition, axis1,
    #                        nodeIdentifier, elementIdentifier,
    #                        vesselMeshGroups=[[stomachMeshGroup, esophagusMeshGroup, abdominalEsoMeshGroup]],
    #                        ostiumMeshGroups=[stomachMeshGroup, esophagogastricJunctionMeshGroup],
    #                        wallAnnotationGroups=ostiumWallAnnotationGroups, coordinates=coordinates)

    # stomachStartNode = nextNodeIdentifier
    # nodeIdentifier = nextNodeIdentifier
    # elementIdentifier = nextElementIdentifier

    # Calculate arclength at quarter line between lesser and greater curvature for each region
    xQuarterEllipseAll = []
    d2QuarterEllipseAll = []
    for n2 in range(len(xEllipseAroundAll)):
        xQuarterEllipseAll.append(xEllipseAroundAll[n2][elementsAroundQuarterDuod])
        d2QuarterEllipseAll.append(d2EllipseAroundAll[n2][elementsAroundQuarterDuod])
    totalQuarterLength = interp.getCubicHermiteCurvesLength(xQuarterEllipseAll, d2QuarterEllipseAll)

    ratioArcLengthSections = []
    for i in range(len(nodeStartEndSections)):
        xList = []
        d2List = []
        for n in range(nodeStartEndSections[i][0], nodeStartEndSections[i][1] + 1):
            xList.append(xQuarterEllipseAll[n])
            d2List.append(d2QuarterEllipseAll[n])
        quarterLengthSection = interp.getCubicHermiteCurvesLength(xList, d2List)
        ratioArcLengthSections.append(quarterLengthSection/totalQuarterLength)

    elementsAlongSections = [math.ceil(c / targetLength) for c in ratioArcLengthSections]
    print(elementsAlongSections)

    # if materialCoordinates:
    #     GEJSettings['Unit scale'] = unitScale
    #     GEJSettings['Ostium diameter'] = ostiumDiameter
    #     GEJSettings['Ostium length'] = ostiumLength
    #     GEJSettings['Ostium wall thickness'] = ostiumWallThickness
    #     GEJSettings['Ostium wall relative thicknesses'] = ostiumWallRelThicknesses
    #     GEJSettings['Vessel inner diameter'] = vesselInnerDiameter
    #     GEJSettings['Vessel wall thickness'] = vesselWallThickness
    #     GEJSettings['Vessel wall relative thicknesses'] = vesselWallRelThicknesses
    #     GEJSettings['Vessel angle 1 degrees'] = vesselAngle1
    #     GEJSettings['Vessel angle 2 degrees'] = vesselAngle2
    #
    # # Create location of annulus
    # xAnnulusOuter = [[] for x in range(elementsCountAroundEso)]
    # xAnnulusOuterPosition = [[] for x in range(elementsCountAroundEso)]
    # d1AnnulusNorm = []
    # d1AnnulusOuter = []
    # for n1 in range(elementsCountAroundEso):
    #     normD2 = vector.normalise(o1_d2[-1][n1])
    #     d1AnnulusNorm.append(normD2)
    #     d1AnnulusOuter.append(vector.setMagnitude(o1_d2[-1][n1], sf))
    #     x = [o1_x[-1][n1][c] + sf * normD2[c] for c in range(3)]
    #     nearestPosition = trackSurfaceStomach.findNearestPosition(x)
    #     xAnnulusOuterPosition[n1] = nearestPosition
    #     xAnnulusOuter[n1] = trackSurfaceStomach.evaluateCoordinates(nearestPosition)
    #
    # d2AnnulusOuter = []
    # for n in range(elementsCountAroundEso):
    #     d = findDerivativeBetweenPoints(xAnnulusOuter[n], xAnnulusOuter[(n + 1) % elementsCountAroundEso])
    #     d2AnnulusOuter.append(d)
    # d2AnnulusOuter = interp.smoothCubicHermiteDerivativesLoop(xAnnulusOuter, d2AnnulusOuter)
    # d3Annulus = []
    # for n in range(elementsCountAroundEso):
    #     d3 = vector.normalise(vector.crossproduct3(vector.normalise(d2AnnulusOuter[n]), d1AnnulusNorm[n]))
    #     d3Annulus.append(d3)
    # annulusD2Curvature = findCurvatureAroundLoop(xAnnulusOuter, d2AnnulusOuter, d3Annulus)
    #
    # # Adjust annulus derivatives
    # # Make d2 on second half of eso point towards duodenum
    # for n1 in range(elementsAroundHalfEso - 1):
    #     idx = n1 + elementsAroundHalfEso + 1
    #     rotAxis = d3Annulus[idx]
    #     rotAngle = math.pi
    #     rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
    #     d2 = d2AnnulusOuter[idx]
    #     d2AnnulusOuter[idx] = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] +
    #                            rotFrame[j][2] * d2[2] for j in range(3)]
    #
    # # Make d1 on first half of eso point towards esophagus
    # for n1 in range(1, elementsAroundHalfEso):
    #     rotAxis = d3Annulus[n1]
    #     rotAngle = math.pi
    #     rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
    #     d1 = d1AnnulusOuter[n1]
    #     d1AnnulusOuter[n1] = [rotFrame[j][0] * d1[0] + rotFrame[j][1] * d1[1] +
    #                           rotFrame[j][2] * d1[2] for j in range(3)]
    #
    # # Flip d1 and d2 on xAnnulusOuter[0]
    # d1 = d1AnnulusOuter[0]  # original d1
    # rotAxis = d3Annulus[0]
    # rotAngle = math.pi
    # rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
    # d1Rot = [rotFrame[j][0] * d1[0] + rotFrame[j][1] * d1[1] + rotFrame[j][2] * d1[2] for j in range(3)]
    #
    # d2 = d2AnnulusOuter[0]  # original d2
    # rotAngle = math.pi
    # rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
    # d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
    #
    # d2AnnulusOuter[0] = d1Rot  # d2 is now d1
    # d1AnnulusOuter[0] = d2Rot  # d1 is now d2 - curvature = annulusD2Curvature[0]
    #
    # # Flip d1 and d2 on xAnnulusOuter[halfEso]
    # d1 = d1AnnulusOuter[elementsAroundHalfEso]  # original d1
    # d2 = d2AnnulusOuter[elementsAroundHalfEso]  # original d2
    # d2AnnulusOuter[elementsAroundHalfEso] = d1  # d2 is now d1
    # d1AnnulusOuter[elementsAroundHalfEso] = d2  # d1 is now d2 - curvature = annulusD2Curvature[halfEso]
    #
    # # x along GC and LC row - needs to be in 2 parts
    # # GC
    # xAlongGCHalfDuod = []
    # endGCEsoP2 = trackSurfaceStomach.getProportion(o1_Positions[0])[1]
    #
    # for n2 in range(len(xEllipseAroundAll)):
    #     xPosition = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[n2][elementsAroundHalfDuod if n2 else 0])
    #     p2 = trackSurfaceStomach.getProportion(xPosition)
    #     if p2[1] < endGCEsoP2:
    #         # xAlongGCHalfDuod.append(xAroundAllTransformed[n2][elementsAroundHalfDuod if n2 else 0])
    #         xAlongGCHalfDuod.append(xEllipseAroundAll[n2][elementsAroundHalfDuod if n2 else 0])
    #     else:
    #         break
    #
    # xAlongGCHalfDuod.append(xAnnulusOuter[0])
    #
    # rotAngle = -math.pi * 2.0 / elementsCountAroundDuod * elementsAroundHalfDuod
    # rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxisApex, rotAngle)
    # d2ApexRot = [rotFrame[j][0] * d2Apex[0] + rotFrame[j][1] * d2Apex[1] + rotFrame[j][2] * d2Apex[2] for j in range(3)]
    #
    # d2AlongGCHalfDuod = [d2ApexRot]
    # for n2 in range(1, len(xAlongGCHalfDuod) - 1):
    #     d = findDerivativeBetweenPoints(xAlongGCHalfDuod[n2], xAlongGCHalfDuod[n2 + 1])
    #     d2AlongGCHalfDuod.append(d)
    # d2AlongGCHalfDuod.append(d2AnnulusOuter[0])
    #
    # xSampledAlongGCHalfDuod, d2SampledAlongGCHalfDuod = \
    #     interp.sampleCubicHermiteCurvesSmooth(xAlongGCHalfDuod + [o1_x[-1][0]], d2AlongGCHalfDuod + [d2AnnulusOuter[0]],
    #                                           elementsAlongFundusApexToCardia + 1,
    #                                           derivativeMagnitudeEnd=vector.magnitude(d2AnnulusOuter[0]))[0:2]
    #
    # d2SampledAlongGCHalfDuod = \
    #     interp.smoothCubicHermiteDerivativesLine(xSampledAlongGCHalfDuod, d2SampledAlongGCHalfDuod,
    #                                              fixEndDerivative=True)
    #
    # xAlongGCHalfDuod = xSampledAlongGCHalfDuod[:-1]
    # d2AlongGCHalfDuod = d2SampledAlongGCHalfDuod[:-1]
    # for n2 in range(len(xAlongGCHalfDuod)):
    #     xEllipseAroundAll[n2][elementsAroundHalfDuod] = xAlongGCHalfDuod[n2]
    #
    # # LC part
    # elementsInBody = elementCountGroupList[1]
    # xAnnulusOuterHalfEsoPosition = trackSurfaceStomach.findNearestPosition(xAnnulusOuter[elementsAroundHalfEso])
    # startLCAnnulus = trackSurfaceStomach.getProportion(xAnnulusOuterHalfEsoPosition)
    #
    # xAlongLCHalfDuod = [xAnnulusOuter[elementsAroundHalfEso]]
    # n2IdxEndBody = elementCountGroupList[0] + elementCountGroupList[1]
    #
    # if elementsInBody - elementsAroundQuarterEso + 1 > 1:
    #     xPositionEndBody = \
    #         trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[n2IdxEndBody][elementsAroundHalfDuod])
    #     xProportionEndBody = trackSurfaceStomach.getProportion(xPositionEndBody)
    #
    #     nx, nd1, nd2, nd3, proportions = \
    #         trackSurfaceStomach.createHermiteCurvePoints(
    #             startLCAnnulus[0], startLCAnnulus[1], xProportionEndBody[0], xProportionEndBody[1],
    #             elementsInBody - elementsAroundQuarterEso + 1, derivativeStart=d2AnnulusOuter[elementsAroundHalfEso])
    #
    #     nxNew, nd1New = \
    #         trackSurfaceStomach.resampleHermiteCurvePointsSmooth(
    #             nx, nd1, nd2, nd3, proportions,
    #             derivativeMagnitudeStart=vector.magnitude(d2AnnulusOuter[elementsAroundHalfEso]))[0:2]
    #     xAlongLCHalfDuod += nxNew[1:-1]
    #
    # for n2 in range(sum(elementCountGroupList) - sum(elementCountGroupList[:2]) + 1):
    #     xAlongLCHalfDuod.append(xEllipseAroundAll[n2IdxEndBody + n2][elementsAroundHalfDuod])
    #
    # d2AlongLCHalfDuod = [d2AnnulusOuter[elementsAroundHalfEso]]
    # for n2 in range(1, len(xAlongLCHalfDuod)):
    #     d = findDerivativeBetweenPoints(xAlongLCHalfDuod[n2 - 1], xAlongLCHalfDuod[n2])
    #     d2AlongLCHalfDuod.append(d)
    # d2AlongLCHalfDuod[-1] = cd1Sections[-1][-1]
    #
    # if materialCoordinates:
    #     # Break point along at body-antrum and pylorus-duod intersection and smooth separately to keep derivative at
    #     # body end
    #     d2AlongBody = \
    #         interp.smoothCubicHermiteDerivativesLine(
    #             xAlongLCHalfDuod[:elementsInBody - elementsAroundQuarterEso + 2],
    #             d2AlongLCHalfDuod[:elementsInBody - elementsAroundQuarterEso + 2],
    #             fixStartDirection=True, fixEndDirection=True)
    #
    #     for n in range(elementCountGroupList[-1] + 1):
    #         d2AlongLCHalfDuod[-(n + 1)] = cd1Sections[-1][-1]
    #
    #     d2AlongBodyToPylorus = \
    #         interp.smoothCubicHermiteDerivativesLine(
    #             xAlongLCHalfDuod[elementsInBody - elementsAroundQuarterEso + 1:-n],
    #             d2AlongLCHalfDuod[elementsInBody - elementsAroundQuarterEso + 1:-n],
    #             fixStartDirection=True, fixEndDirection=True)
    #     d2AlongDuod = \
    #         interp.smoothCubicHermiteDerivativesLine(
    #             xAlongLCHalfDuod[-(n + 1):], d2AlongLCHalfDuod[-(n + 1):], fixStartDirection=True, fixEndDirection=True)
    #     d2AlongLCHalfDuod = d2AlongBody[:-1] + d2AlongBodyToPylorus[:-1] + d2AlongDuod
    # else:
    #     d2AlongLCHalfDuod = \
    #         interp.smoothCubicHermiteDerivativesLine(xAlongLCHalfDuod, d2AlongLCHalfDuod,
    #                                                  fixStartDerivative=True, fixEndDirection=True)
    #
    # for n2 in range(len(xAlongLCHalfDuod)):
    #     xEllipseAroundAll[n2 + elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2][elementsAroundHalfDuod] = \
    #         xAlongLCHalfDuod[n2]
    #
    # xAlongCurvatures = xAlongGCHalfDuod + xAlongLCHalfDuod
    # d2AlongCurvatures = d2AlongGCHalfDuod + d2AlongLCHalfDuod
    #
    # esoCount = 1
    # # Find arclengths to points on quarter, GC and 3-quarter fundus as trackSurface search for nearest position doesnt
    # # work well near apex
    # proportionsAlongQuarterMidThree = []
    # count = 0
    # for n1 in [elementsAroundQuarterDuod, elementsAroundHalfDuod, elementsAroundQuarterDuod + elementsAroundHalfDuod]:
    #     xAlong = []
    #     arcLength = 0.0
    #     proportions = []
    #     for n2 in range(elementsAlongFundusApexToCardia):
    #         xAlong.append(xEllipseAroundAll[n2][n1])
    #
    #     rotAngle = -math.pi * 2.0 / elementsCountAroundDuod * n1
    #     rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxisApex, rotAngle)
    #     d2ApexRot = [rotFrame[j][0] * d2Apex[0] + rotFrame[j][1] * d2Apex[1] +
    #                  rotFrame[j][2] * d2Apex[2] for j in range(3)]
    #     d2Along = [d2ApexRot]
    #
    #     for n2 in range(len(xAlong) - 1):
    #         d2Along.append(findDerivativeBetweenPoints(xAlong[n2], xAlong[n2 + 1]))
    #
    #     for n2 in range(len(xAlong) - 1):
    #         arcLength += interp.getCubicHermiteArcLength(xAlong[n2], d2Along[n2], xAlong[n2 + 1], d2Along[n2 + 1])
    #         proportion = arcLength / arcLengthsAlongTS[count]
    #         proportions.append(proportion)
    #     proportionsAlongQuarterMidThree.append(proportions)
    #     count += 1
    #
    # for n2 in range(1, len(xEllipseAroundAll) - 1):
    #     startAnnulus = n2 == elementsAlongFundusApexToCardia
    #     endAnnulus = n2 == elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2
    #     interiorAnnulus = \
    #         elementsAlongFundusApexToCardia < n2 < elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2
    #
    #     if startAnnulus or interiorAnnulus or endAnnulus:
    #         xEllipseAroundAll[n2][elementsAroundHalfDuod - 1] = xAnnulusOuter[esoCount]
    #         d1EllipseAroundAll[n2][elementsAroundHalfDuod - 1] = d1AnnulusOuter[esoCount]
    #         xEllipseAroundAll[n2][elementsAroundHalfDuod + 1] = xAnnulusOuter[-esoCount]
    #         d1EllipseAroundAll[n2][elementsAroundHalfDuod + 1] = d1AnnulusOuter[-esoCount]
    #
    #         if startAnnulus or endAnnulus:
    #             xEllipseAroundAll[n2][elementsAroundHalfDuod] = \
    #                 xAnnulusOuter[0 if startAnnulus else elementsAroundHalfEso]
    #
    #         xPositionA = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[n2][elementsAroundQuarterDuod])
    #         xProportionA = trackSurfaceStomach.getProportion(xPositionA)
    #         xPositionB = trackSurfaceStomach.findNearestPosition(xAnnulusOuter[esoCount])
    #         xProportionB = trackSurfaceStomach.getProportion(xPositionB)
    #
    #         nx, nd1, nd2, nd3, proportions = \
    #             trackSurfaceStomach.createHermiteCurvePoints(
    #                 xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1], elementsAroundQuarterDuod - 1,
    #                 derivativeStart=d1EllipseAroundAll[n2][elementsAroundQuarterDuod])
    #
    #         nx, nd1 = \
    #             trackSurfaceStomach.resampleHermiteCurvePointsSmooth(
    #                 nx, nd1, nd2, nd3, proportions,
    #                 derivativeMagnitudeStart=vector.magnitude(d1EllipseAroundAll[n2][elementsAroundQuarterDuod]))[0:2]
    #
    #         xEllipseAroundAll[n2][elementsAroundQuarterDuod + 1:elementsAroundHalfDuod - 1] = nx[1:-1]
    #         d1EllipseAroundAll[n2][elementsAroundQuarterDuod + 1:elementsAroundHalfDuod - 1] = nd1[1:-1]
    #
    #         # Smooth two halves separately and join together again with d1 at xOuterAnnulus[0] and
    #         # xOuterAnnulus[elementsAroundHalfEso]
    #         d1SmoothedFirstHalf = \
    #             interp.smoothCubicHermiteDerivativesLine(xEllipseAroundAll[n2][:elementsAroundHalfDuod],
    #                                                      d1EllipseAroundAll[n2][:elementsAroundHalfDuod],
    #                                                      fixEndDerivative=True if interiorAnnulus else None)
    #
    #         if interiorAnnulus:
    #             d1EllipseAroundAll[n2][:elementsAroundHalfDuod] = d1SmoothedFirstHalf
    #         else:
    #             d1SmoothedFirstHalf[-1] = vector.setMagnitude(d1SmoothedFirstHalf[-1],
    #                                                           0.5 * (vector.magnitude(d1SmoothedFirstHalf[-2]) +
    #                                                           vector.magnitude(d1AnnulusOuter[1])))
    #
    #         xPositionA = trackSurfaceStomach.findNearestPosition(xAnnulusOuter[-esoCount])
    #         xProportionA = trackSurfaceStomach.getProportion(xPositionA)
    #         xPositionB = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[n2][elementsAroundHalfDuod +
    #                                                                                    elementsAroundQuarterDuod])
    #         xProportionB = trackSurfaceStomach.getProportion(xPositionB)
    #
    #         nx, nd1, nd2, nd3, proportions = \
    #             trackSurfaceStomach.createHermiteCurvePoints(
    #                 xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1],
    #                 elementsAroundQuarterDuod - 1,
    #                 derivativeEnd=d1EllipseAroundAll[n2][elementsAroundHalfDuod + elementsAroundQuarterDuod])
    #
    #         nx, nd1 = \
    #             trackSurfaceStomach.resampleHermiteCurvePointsSmooth(
    #                 nx, nd1, nd2, nd3, proportions,
    #                 derivativeMagnitudeEnd=vector.magnitude(d1EllipseAroundAll[n2][elementsAroundHalfDuod +
    #                                                                                elementsAroundQuarterDuod]))[0:2]
    #
    #         xEllipseAroundAll[n2][elementsAroundHalfDuod + 2:elementsAroundHalfDuod + elementsAroundQuarterDuod] = \
    #             nx[1:-1]
    #         d1EllipseAroundAll[n2][elementsAroundHalfDuod + 2:elementsAroundHalfDuod + elementsAroundQuarterDuod] = \
    #             nd1[1:-1]
    #
    #         esoCount += 1
    #         if startAnnulus or endAnnulus:
    #             d1SmoothedSecondHalf = \
    #                 interp.smoothCubicHermiteDerivativesLine(
    #                     xEllipseAroundAll[n2][elementsAroundHalfDuod + 1:] + [xEllipseAroundAll[n2][0]],
    #                     d1EllipseAroundAll[n2][elementsAroundHalfDuod + 1:] + [d1SmoothedFirstHalf[0]],
    #                     fixEndDerivative=True)
    #             d1SmoothedSecondHalf[0] = vector.setMagnitude(d1SmoothedSecondHalf[0],
    #                                                           0.5 * (vector.magnitude(d1SmoothedFirstHalf[1]) +
    #                                                                  vector.magnitude(d1AnnulusOuter[1])))
    #
    #             d1EllipseAroundAll[n2] = d1SmoothedFirstHalf + \
    #                                      [d1AnnulusOuter[0 if startAnnulus else elementsAroundHalfEso]] + \
    #                                      d1SmoothedSecondHalf[:-1]
    #         else:
    #             d1Smoothed = \
    #                 interp.smoothCubicHermiteDerivativesLine(xEllipseAroundAll[n2][elementsAroundHalfDuod + 1:] +
    #                                                          [xEllipseAroundAll[n2][0]],
    #                                                          d1EllipseAroundAll[n2][elementsAroundHalfDuod + 1:] +
    #                                                          [d1EllipseAroundAll[n2][0]],
    #                                                          fixStartDerivative=True, fixEndDerivative=True)
    #
    #             d1EllipseAroundAll[n2][elementsAroundHalfDuod + 1:] = d1Smoothed[:-1]
    #
    #             del xEllipseAroundAll[n2][elementsAroundHalfDuod]
    #             del d1EllipseAroundAll[n2][elementsAroundHalfDuod]
    #             del d2EllipseAroundAll[n2][elementsAroundHalfDuod]
    #
    #     elif 0 < n2 < elementsAlongFundusApexToCardia or \
    #             elementsAlongFundusApexToCardia + elementsAroundHalfEso - 1 <= n2 < elementsAlongFundusApexToCardia + \
    #             elementsAroundHalfEso + elementsInBody - elementsAroundQuarterEso - 1:
    #             if 0 < n2 < elementsAlongFundusApexToCardia:
    #                 xToSample = [xEllipseAroundAll[n2][elementsAroundQuarterDuod]] + \
    #                             [xAlongGCHalfDuod[n2]] + \
    #                             [xEllipseAroundAll[n2][elementsAroundQuarterDuod + elementsAroundHalfDuod]]
    #                 xToSampleProportion0 = [0.25, 0.5, 0.75]
    #                 xToSampleProportion1 = [proportionsAlongQuarterMidThree[0][n2 - 1],
    #                                         proportionsAlongQuarterMidThree[1][n2 - 1],
    #                                         proportionsAlongQuarterMidThree[2][n2 - 1]]
    #             else:
    #                 n2Idx = n2 - (elementsAlongFundusApexToCardia + elementsAroundHalfEso - 1) + 1
    #                 xToSample = [xEllipseAroundAll[n2][elementsAroundQuarterDuod]] + \
    #                             [xAlongLCHalfDuod[n2Idx]] + \
    #                             [xEllipseAroundAll[n2][elementsAroundQuarterDuod + elementsAroundHalfDuod]]
    #
    #             d1ToSample = [d1EllipseAroundAll[n2][elementsAroundQuarterDuod]] + \
    #                          [d1EllipseAroundAll[n2][elementsAroundHalfDuod]] + \
    #                          [d1EllipseAroundAll[n2][elementsAroundHalfDuod + elementsAroundQuarterDuod]]
    #
    #             if 0 < n2 < elementsAlongFundusApexToCardia:
    #                 xProportionA = [xToSampleProportion0[0], xToSampleProportion1[0]]
    #                 xProportionB = [xToSampleProportion0[1], xToSampleProportion1[1]]
    #                 xProportionC = [xToSampleProportion0[2], xToSampleProportion1[2]]
    #
    #             else:
    #                 xPositionA = trackSurfaceStomach.findNearestPosition(xToSample[0])
    #                 xProportionA = trackSurfaceStomach.getProportion(xPositionA)
    #                 xPositionB = trackSurfaceStomach.findNearestPosition(xToSample[1])
    #                 xProportionB = trackSurfaceStomach.getProportion(xPositionB)
    #                 xPositionC = trackSurfaceStomach.findNearestPosition(xToSample[2])
    #                 xProportionC = trackSurfaceStomach.getProportion(xPositionC)
    #
    #             nx, nd1, nd2, nd3, proportions = \
    #                 trackSurfaceStomach.createHermiteCurvePoints(
    #                     xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1],
    #                     elementsAroundQuarterDuod, derivativeStart=d1ToSample[0], derivativeEnd=d1ToSample[1])
    #             if n2 == 1:
    #                 sampledProportionsOneLeft = proportions
    #             nxLeft, nd1Left = \
    #                 trackSurfaceStomach.resampleHermiteCurvePointsSmooth(
    #                     nx, nd1, nd2, nd3, proportions, derivativeMagnitudeStart=vector.magnitude(d1ToSample[0]))[0:2]
    #
    #             nx, nd1, nd2, nd3, proportions = \
    #                 trackSurfaceStomach.createHermiteCurvePoints(
    #                     xProportionB[0], xProportionB[1],  xProportionC[0], xProportionC[1], elementsAroundQuarterDuod,
    #                     derivativeStart=d1ToSample[1], derivativeEnd=d1ToSample[2])
    #             if n2 == 1:
    #                 sampledProportionsOneRight = proportions
    #             nxRight, nd1Right = \
    #                 trackSurfaceStomach.resampleHermiteCurvePointsSmooth(
    #                     nx, nd1, nd2, nd3, proportions, derivativeMagnitudeEnd=vector.magnitude(d1ToSample[2]))[0:2]
    #
    #             xNew = nxLeft[1:] + nxRight[1:-1]
    #             d1New = nd1Left[1:] + nd1Right[1:-1]
    #
    #             xEllipseAroundAll[n2][elementsAroundQuarterDuod + 1:-elementsAroundQuarterDuod] = xNew
    #             d1EllipseAroundAll[n2][elementsAroundQuarterDuod + 1:-elementsAroundQuarterDuod] = d1New
    #
    #             d1EllipseAroundAll[n2] = interp.smoothCubicHermiteDerivativesLoop(xEllipseAroundAll[n2],
    #                                                                               d1EllipseAroundAll[n2])
    #
    # # Calculate and smooth d2
    # xAlongAll = []
    # for n1 in range(elementsAroundHalfDuod - 1):
    #     xAlong = []
    #     for n2 in range(len(xEllipseAroundAll)):
    #         xAlong.append(xEllipseAroundAll[n2][n1 if n2 else 0])
    #     xAlongAll.append(xAlong)
    #
    # # Row with annulus - left
    # xAlong = []
    # for n2 in range(len(xEllipseAroundAll)):
    #     xAlong.append(xEllipseAroundAll[n2][elementsAroundHalfDuod - 1 if n2 else 0])
    # xAlongAll.append(xAlong)
    #
    # xAlongAll.append(xAlongCurvatures)
    #
    # for n1 in range(elementsAroundHalfDuod - 1):
    #     xAlong = []
    #     for n2 in range(len(xEllipseAroundAll)):
    #         if n2 == 0:
    #             ringIdx = 0
    #         elif elementsAlongFundusApexToCardia < n2 < elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2:
    #             ringIdx = n1 + elementsAroundHalfDuod
    #         else:
    #             ringIdx = n1 + elementsAroundHalfDuod + 1
    #         xAlong.append(xEllipseAroundAll[n2][ringIdx])
    #     xAlongAll.append(xAlong)
    #
    # d2AlongAll = []
    # for n1 in range(elementsCountAroundDuod):
    #     rotAngle = -math.pi * 2.0 / elementsCountAroundDuod * n1
    #     rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxisApex, rotAngle)
    #     d2ApexRot = [rotFrame[j][0] * d2Apex[0] + rotFrame[j][1] * d2Apex[1] +
    #                  rotFrame[j][2] * d2Apex[2] for j in range(3)]
    #
    #     if n1 == elementsAroundHalfDuod:  # Process LC and GC
    #         # GC
    #         # Resample point downstream from apex
    #         nx = \
    #             trackSurfaceStomach.createHermiteCurvePoints(
    #                 sampledProportionsOneLeft[-2][0], sampledProportionsOneLeft[-2][1],
    #                 sampledProportionsOneRight[1][0], sampledProportionsOneRight[1][1], 2,
    #                 derivativeStart=d1EllipseAroundAll[1][elementsAroundHalfDuod - 1],
    #                 derivativeEnd=d1EllipseAroundAll[1][elementsAroundHalfDuod + 1])[0]
    #
    #         xEllipseAroundAll[1][elementsAroundHalfDuod] = nx[1]
    #
    #         d1EllipseAroundAll[1] = interp.smoothCubicHermiteDerivativesLoop(xEllipseAroundAll[1],
    #                                                                          d1EllipseAroundAll[1])
    #
    #         xAlongGCHalfDuodNew = []
    #         for n2 in range(elementsAlongFundusApexToCardia + 1):
    #             xAlongGCHalfDuodNew.append(xEllipseAroundAll[n2][elementsAroundHalfDuod])
    #
    #         d2AlongGCHalfDuodNew = [d2ApexRot]
    #         for n2 in range(1, len(xAlongGCHalfDuodNew) - 1):
    #             d = findDerivativeBetweenPoints(xAlongGCHalfDuodNew[n2], xAlongGCHalfDuodNew[n2 + 1])
    #             d2AlongGCHalfDuodNew.append(d)
    #         d2AlongGCHalfDuodNew.append(d2AnnulusOuter[0])
    #
    #         d2AlongGCHalfDuodNew = \
    #             interp.smoothCubicHermiteDerivativesLine(xAlongGCHalfDuodNew + [o1_x[-1][0]],
    #                                                      d2AlongGCHalfDuodNew + [d2AnnulusOuter[0]],
    #                                                      fixEndDerivative=True)[:-1]
    #         d2AlongCurvatures[:len(d2AlongGCHalfDuodNew)] = d2AlongGCHalfDuodNew
    #         d2AlongAll.append(d2AlongCurvatures)
    #
    #     else:
    #         d2Along = [d2ApexRot]
    #         for n2 in range(1, len(xAlongAll[n1])):
    #             d = findDerivativeBetweenPoints(xAlongAll[n1][n2 - 1], xAlongAll[n1][n2])
    #             d2Along.append(d)
    #         d2Along[-1] = cd1Sections[-1][-1]
    #
    #         if materialCoordinates:
    #             # Break point along at body-antrum and pylorus-duod intersection and smooth separately to keep
    #             # derivative at body end
    #             d2AlongFundusBody = \
    #                 interp.smoothCubicHermiteDerivativesLine(
    #                     xAlongAll[n1][:elementCountGroupList[0] + elementCountGroupList[1] + 1],
    #                     d2Along[:elementCountGroupList[0] + elementCountGroupList[1] + 1],
    #                     fixStartDirection=True, fixEndDirection=True)
    #
    #             for n in range(elementCountGroupList[-1] + 1):
    #                 d2Along[-(n + 1)] = cd1Sections[-1][-1]
    #
    #             d2AlongBodyToPylorus = \
    #                 interp.smoothCubicHermiteDerivativesLine(
    #                     xAlongAll[n1][elementCountGroupList[0] + elementCountGroupList[1]:-n],
    #                     d2Along[elementCountGroupList[0] + elementCountGroupList[1]:-n],
    #                     fixStartDirection=True, fixEndDirection=True)
    #
    #             d2AlongDuod = \
    #                 interp.smoothCubicHermiteDerivativesLine(
    #                     xAlongAll[n1][-(n + 1):], d2Along[-(n + 1):], fixStartDirection=True, fixEndDirection=True)
    #
    #             d2Along = d2AlongFundusBody[:-1] + d2AlongBodyToPylorus[:-1] + d2AlongDuod
    #         else:
    #             d2Along = interp.smoothCubicHermiteDerivativesLine(xAlongAll[n1], d2Along, fixStartDirection=True,
    #                                                                fixEndDirection=True)
    #         if n1 == elementsAroundHalfDuod - 1 or n1 == elementsAroundHalfDuod + 1:
    #             d2Along[elementsAlongFundusApexToCardia] = \
    #                 vector.setMagnitude(d2Along[elementsAlongFundusApexToCardia],
    #                                     0.5 * (vector.magnitude(d2Along[elementsAlongFundusApexToCardia]) +
    #                                            vector.magnitude(d2AnnulusOuter[1])))
    #             d2Along[elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2] = \
    #                 vector.setMagnitude(d2Along[elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2],
    #                                     0.5 * (vector.magnitude(d2Along[elementsAlongFundusApexToCardia +
    #                                                                     elementsAroundHalfEso - 2]) +
    #                                            vector.magnitude(d2AnnulusOuter[1])))
    #         d2AlongAll.append(d2Along)
    #
    # # Replace d2 for points sitting on annulus
    # for n2 in range(2, elementsAroundHalfEso - 1):
    #     n2Idx = n2 + elementsAlongFundusApexToCardia - 1
    #     d2AlongAll[elementsAroundHalfDuod - 1][n2Idx] = d2AnnulusOuter[n2]
    #     d2AlongAll[elementsAroundHalfDuod + 1][n2Idx] = d2AnnulusOuter[-n2]
    #
    # # Re-arrange back to around followed by along
    # for n2 in range(len(xEllipseAroundAll)):
    #     incompleteRingsWithinEso = elementsAlongFundusApexToCardia < n2 < elementsAlongFundusApexToCardia + \
    #                                elementsAroundHalfEso - 2
    #     for n1 in range(len(xEllipseAroundAll[n2]) + (1 if incompleteRingsWithinEso else 0)):
    #         if incompleteRingsWithinEso:
    #             if n1 == elementsAroundHalfDuod:
    #                 pass
    #             else:
    #                 d2EllipseAroundAll[n2][n1 - (1 if n1 > elementsAroundHalfDuod else 0)] = d2AlongAll[n1][n2]
    #         elif n2 >= elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2:
    #             d2EllipseAroundAll[n2][n1] = \
    #                 d2AlongAll[n1][n2 - (elementsAroundHalfEso + 1 - 4 if n1 == elementsAroundHalfDuod else 0)]
    #         else:
    #             d2EllipseAroundAll[n2][n1] = d2AlongAll[n1][n2]
    #
    # xOuter = xEllipseAroundAll
    # d1Outer = d1EllipseAroundAll
    # d2Outer = d2EllipseAroundAll
    #
    # # Calculate d3
    # d3UnitOuter = []
    # for n2 in range(len(xOuter)):
    #     d3Around = []
    #     for n1 in range(len(xOuter[n2])):
    #         if n2 == 0:
    #             d3Around.append(rotAxisApex)
    #         else:
    #             d3Around.append(vector.normalise(
    #                 vector.crossproduct3(vector.normalise(d1Outer[n2][n1]), vector.normalise(d2Outer[n2][n1]))))
    #     d3UnitOuter.append(d3Around)
    #
    # # Calculate curvature around
    # d1Curvature = []
    # d1Curvature.append([1.0 for n in range(len(xOuter[0]))])  # Will be replaced in later step
    # esoCount = 1
    # for n2 in range(1, len(xOuter)):
    #     incompleteRingsWithinEso = elementsAlongFundusApexToCardia < n2 < elementsAlongFundusApexToCardia + \
    #                                elementsAroundHalfEso - 2
    #     completeRingsOnCardia = (n2 == elementsAlongFundusApexToCardia or n2 == elementsAlongFundusApexToCardia +
    #                              elementsAroundHalfEso - 2)
    #     if incompleteRingsWithinEso:
    #         d2 = o1_d2[-1][esoCount]
    #         rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(o1_d1[-1][esoCount]),
    #                                                         vector.normalise(o1_d2[-1][esoCount])))
    #         rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
    #         d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
    #         d1CurvatureFirstHalf = findCurvatureAlongLine(xOuter[n2][:elementsAroundHalfDuod] + [o1_x[-1][esoCount]],
    #                                                       d1Outer[n2][:elementsAroundHalfDuod] + [d2Rot],
    #                                                       d3UnitOuter[n2][:elementsAroundHalfDuod] + [rotAxis])
    #         curvature = interp.getCubicHermiteCurvature(xAnnulusOuter[esoCount], d1AnnulusOuter[esoCount],
    #                                                     o1_x[-1][esoCount], d2Rot,
    #                                                     d3UnitOuter[n2][elementsAroundHalfDuod - 1], 0.0)
    #         d1CurvatureFirstHalf[-2] = curvature
    #
    #         d3 = vector.normalise(vector.crossproduct3(vector.normalise(o1_d1[-1][-esoCount]),
    #                                                    vector.normalise(o1_d2[-1][-esoCount])))
    #         d1CurvatureSecondHalf = findCurvatureAlongLine([o1_x[-1][-esoCount]] +
    #                                                        xOuter[n2][elementsAroundHalfDuod:] + [xOuter[n2][0]],
    #                                                        [o1_d2[-1][-esoCount]] +
    #                                                        d1Outer[n2][elementsAroundHalfDuod:] + [d1Outer[n2][0]],
    #                                                        [d3] + d3UnitOuter[n2][elementsAroundHalfDuod:] +
    #                                                        [d3UnitOuter[n2][0]])[:-1]
    #
    #         curvature = interp.getCubicHermiteCurvature(o1_x[-1][-esoCount], o1_d2[-1][-esoCount],
    #                                                     xAnnulusOuter[-esoCount], d1AnnulusOuter[-esoCount],
    #                                                     d3UnitOuter[n2][elementsAroundHalfDuod], 1.0)
    #         d1CurvatureSecondHalf[1] = curvature
    #         d1Curvature.append(d1CurvatureFirstHalf[:-1] + d1CurvatureSecondHalf[1:])
    #         esoCount += 1
    #
    #     elif completeRingsOnCardia:
    #         d1CurvatureFirstHalf = findCurvatureAlongLine(xOuter[n2][:elementsAroundHalfDuod],
    #                                                       d1Outer[n2][:elementsAroundHalfDuod],
    #                                                       d3UnitOuter[n2][:elementsAroundHalfDuod])
    #
    #         d1CurvatureSecondHalf = findCurvatureAlongLine(xOuter[n2][elementsAroundHalfDuod + 1:] + [xOuter[n2][0]],
    #                                                        d1Outer[n2][elementsAroundHalfDuod + 1:] + [d1Outer[n2][0]],
    #                                                        d3UnitOuter[n2][elementsAroundHalfDuod + 1:] +
    #                                                        [d3UnitOuter[n2][0]])[:-1]
    #
    #         midPtCurvature = annulusD2Curvature[0 if n2 == elementsAlongFundusApexToCardia else elementsAroundHalfEso]
    #         d1Curvature.append(d1CurvatureFirstHalf + [midPtCurvature] + d1CurvatureSecondHalf)
    #         esoCount += 1
    #     else:
    #         d1Curvature.append(findD1CurvatureAround(xOuter[n2], d1Outer[n2], d3UnitOuter[n2]))
    #
    # # Populate d3Along for use to calculate curvature along
    # d3UnitAlongAll = []
    # for n1 in range(elementsAroundHalfDuod - 1):
    #     d3Along = []
    #     for n2 in range(len(d3UnitOuter)):
    #         d3Along.append(d3UnitOuter[n2][n1 if n2 else 0])
    #     d3UnitAlongAll.append(d3Along)
    #
    # # Row wth annulus - left
    # d3Along = []
    # for n2 in range(len(d3UnitOuter)):
    #     d3Along.append(d3UnitOuter[n2][elementsAroundHalfDuod - 1 if n2 else 0])
    # d3UnitAlongAll.append(d3Along)
    #
    # # GC and LC row - needs to be in 2 parts
    # # GC part
    # d3AlongGCHalfDuod = []
    # for n2 in range(elementsAlongFundusApexToCardia + 1):
    #     d3AlongGCHalfDuod.append(d3UnitOuter[n2][elementsAroundHalfDuod if n2 else 0])
    #
    # # LC part
    # d3AlongLCHalfDuod = []
    # for n2 in range(elementsAlongCardiaToDuod + 1):
    #     d3AlongLCHalfDuod.append(d3UnitOuter[n2 + elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2]
    #                             [elementsAroundHalfDuod])
    # d3Along = d3AlongGCHalfDuod + d3AlongLCHalfDuod
    # d3UnitAlongAll.append(d3Along)
    #
    # for n1 in range(elementsAroundHalfDuod - 1):
    #     d3Along = []
    #     for n2 in range(len(d3UnitOuter)):
    #         if n2 == 0:
    #             ringIdx = 0
    #         elif elementsAlongFundusApexToCardia < n2 < elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2:
    #             ringIdx = n1 + elementsAroundHalfDuod
    #         else:
    #             ringIdx = n1 + elementsAroundHalfDuod + 1
    #         d3Along.append(d3UnitOuter[n2][ringIdx])
    #     d3UnitAlongAll.append(d3Along)
    #
    # # Calculate curvature along
    # d2CurvatureAlong = []
    # for n1 in range(len(xAlongAll)):
    #     if n1 == elementsAroundHalfDuod:  # Process LC and GC
    #         # GC
    #         d2CurvatureAlongGCHalfDuod = findCurvatureAlongLine(xAlongGCHalfDuod, d2AlongGCHalfDuod, d3AlongGCHalfDuod)
    #         # LC
    #         d2CurvatureAlongLCHalfDuod = findCurvatureAlongLine(xAlongLCHalfDuod, d2AlongLCHalfDuod, d3AlongLCHalfDuod)
    #         d2CurvaturesAlongCurvature = d2CurvatureAlongGCHalfDuod + d2CurvatureAlongLCHalfDuod
    #         d2CurvatureAlong.append(d2CurvaturesAlongCurvature)
    #     else:
    #         curvature = findCurvatureAlongLine(xAlongAll[n1], d2AlongAll[n1], d3UnitAlongAll[n1])
    #         # replace with curvature from annulus
    #         if n1 == elementsAroundHalfDuod - 1 or n1 == elementsAroundHalfDuod + 1:
    #             for n2 in range(2, elementsAroundHalfEso - 1):
    #                 n2Idx = n2 + elementsAlongFundusApexToCardia + 1
    #                 curvature[n2Idx] = annulusD2Curvature[(n2 + 1 if elementsAroundHalfDuod - 1 else -(n2 + 1))]
    #         d2CurvatureAlong.append(curvature)
    #
    # # Re-arrange back to around followed by along
    # for n2 in range(len(xEllipseAroundAll)):
    #     incompleteRingsWithinEso = elementsAlongFundusApexToCardia < n2 < elementsAlongFundusApexToCardia + \
    #                                elementsAroundHalfEso - 2
    #     for n1 in range(len(xEllipseAroundAll[n2]) + (1 if incompleteRingsWithinEso else 0)):
    #         if incompleteRingsWithinEso:
    #             if n1 == elementsAroundHalfDuod:
    #                 pass
    #             else:
    #                 d2Curvature[n2][n1 - (1 if n1 > elementsAroundHalfDuod else 0)] = d2CurvatureAlong[n1][n2]
    #         elif n2 >= elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2:
    #             d2Curvature[n2][n1] = \
    #                 d2CurvatureAlong[n1][n2 - (elementsAroundHalfEso + 1 - 4 if n1 == elementsAroundHalfDuod else 0)]
    #         else:
    #             d2Curvature[n2][n1] = d2CurvatureAlong[n1][n2]
    #
    # # Replace d1Curvature at apex with d2Curvature
    # for n1 in range(elementsCountAroundDuod):
    #     d1Curvature[0][n1] = d2Curvature[0][n1]
    #
    # # Create inner nodes
    # xList = []
    # d1List = []
    # d2List = []
    # d3List = []
    # nodeIdx = stomachStartNode
    # idxMat = []
    #
    # if elementsCountThroughWall > 1:
    #     thicknessProportionsUI = [0.0, mucosaRelThickness, submucosaRelThickness, circularRelThickness,
    #                               longitudinalRelThickness, longitudinalRelThickness]
    #     thicknessProportions = [thicknessProportion / sum(thicknessProportionsUI[:-1])
    #                             for thicknessProportion in thicknessProportionsUI]
    #
    #     xi3List = []
    #     xi3 = 0.0
    #     for i in range(len(thicknessProportions) - 1):
    #         xi3 += thicknessProportions[i]
    #         xi3List.append(xi3)
    #
    # for n2 in range(len(xOuter)):
    #     idxThroughWall = []
    #     for n3 in range(elementsCountThroughWall + 1):
    #         xi3 = xi3List[n3] if elementsCountThroughWall > 1 else 1.0 / elementsCountThroughWall * n3
    #         idxAround = []
    #         for n1 in range(1 if n2 == 0 else len(xOuter[n2])):
    #             # Coordinates
    #             norm = d3UnitOuter[n2][n1]
    #             xOut = xOuter[n2][n1]
    #             xIn = [xOut[i] - norm[i] * wallThickness for i in range(3)]
    #             dWall = [wallThickness * c for c in norm]
    #             x = interp.interpolateCubicHermite(xIn, dWall, xOut, dWall, xi3)
    #             xList.append(x)
    #
    #             # d1
    #             factor = 1.0 + wallThickness * (1.0 - xi3) * d1Curvature[n2][n1]
    #             d1 = [factor * c for c in d1Outer[n2][n1]]
    #             d1List.append(d1)
    #
    #             # d2
    #             factor = 1.0 + wallThickness * (1.0 - xi3) * d2Curvature[n2][n1]
    #             d2 = [factor * c for c in d2Outer[n2][n1]]
    #             d2List.append(d2)
    #
    #             # d3
    #             d3 = [c * wallThickness * (thicknessProportions[n3 + 1] if elementsCountThroughWall > 1 else 1.0)
    #                   for c in norm]
    #             d3List.append(d3)
    #
    #             idxAround.append(nodeIdx)
    #             nodeIdx += 1
    #         idxThroughWall.append(idxAround)
    #     idxMat.append(idxThroughWall)
    #
    # nodeIdxGC = []
    # for n2 in range(len(idxMat)):
    #     for n3 in range(len(idxMat[n2])):
    #         if n2 == 0:
    #             nodeIdxGC += idxMat[n2][n3]
    #         else:
    #             nodeIdxGC.append(idxMat[n2][n3][0])
    #
    # for n2 in range(1, elementsAlongFundusApexToCardia + 1):
    #     for n3 in range(len(idxMat[n2])):
    #         nodeIdxGC.append(idxMat[n2][n3][int(0.5 * len(xOuter[n2]))])
    #
    # nodeIdxLC = []
    # for n2 in range(elementsAlongCardiaToDuod + 1):
    #     for n3 in range(len(idxMat[n2])):
    #         nodeIdxLC.append(
    #             idxMat[elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2 + n2][n3][elementsAroundHalfDuod])
    #
    # for n2 in range(len(xList)):
    #     node = nodes.createNode(nodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1List[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2List[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3List[n2])
    #     if useCrossDerivatives:
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
    #     nodeIdentifier += 1
    #
    # annotationGroupsAlong = []
    # for i in range(len(elementCountGroupList)):
    #     elementsCount = elementCountGroupList[i]
    #     for n in range(elementsCount):
    #         annotationGroupsAlong.append(annotationGroupAlong[i])
    #
    # # Create elements
    # fundusMucosaElementIdentifiers = []
    # elementIdxMat = []
    # n = 0
    # for n2 in range(elementsAlongEsophagus):
    #     elementIdxThroughWall = []
    #     for n3 in range(elementsThroughEsophagusWall):
    #         elementIdxAround = []
    #         for n1 in range(elementsCountAroundEso):
    #             n += 1
    #             elementIdxAround.append(n)
    #         elementIdxThroughWall.append(elementIdxAround)
    #     elementIdxMat.append(elementIdxThroughWall)
    #
    # if useCubicHermiteThroughWall:
    #     eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    # else:
    #     eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    # eftStandard = eftfactory.createEftBasic()
    #
    # elementtemplateStandard = mesh.createElementtemplate()
    # elementtemplateStandard.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    # elementtemplateStandard.defineField(coordinates, -1, eftStandard)
    #
    # elementtemplateX = mesh.createElementtemplate()
    # elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    #
    # fundusElements = elementCountGroupList[0]
    # radiansPerElementAroundDuod = math.pi * 2.0 / elementsCountAroundDuod
    #
    # for e2 in range(len(xOuter) - 1):
    #     elementIdxThroughWall = []
    #     if e2 == 0:  # pole
    #         for e3 in range(elementsCountThroughWall):
    #             elementIdxAround = []
    #             for e1 in range(elementsCountAroundDuod):
    #                 va = e1
    #                 vb = (e1 + 1) % elementsCountAroundDuod
    #                 eft1 = eftfactory.createEftShellPoleBottom(va*100, vb*100)
    #                 elementtemplateX.defineField(coordinates, -1, eft1)
    #                 element = mesh.createElement(elementIdentifier, elementtemplateX)
    #                 bni1 = e3 + stomachStartNode
    #                 bni2 = (elementsCountThroughWall + 1) + stomachStartNode + elementsCountAroundDuod * e3 + e1
    #                 bni3 = (elementsCountThroughWall + 1) + stomachStartNode + elementsCountAroundDuod * e3 + \
    #                        (e1 + 1) % elementsCountAroundDuod
    #                 nodeIdentifiers = [bni1, bni2, bni3, bni1 + 1, bni2 + elementsCountAroundDuod,
    #                                    bni3 + elementsCountAroundDuod]
    #
    #                 element.setNodesByIdentifier(eft1, nodeIdentifiers)
    #                 # set general linear map coefficients
    #                 radiansAround = e1 * radiansPerElementAroundDuod
    #                 radiansAroundNext = ((e1 + 1) % elementsCountAroundDuod) * radiansPerElementAroundDuod
    #                 scalefactors = [
    #                     -1.0,
    #                     math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAroundDuod,
    #                     math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAroundDuod,
    #                     math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAroundDuod,
    #                     math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAroundDuod
    #                 ]
    #                 element.setScaleFactors(eft1, scalefactors)
    #                 if e2 < fundusElements and limitingRidge and elementsCountThroughWall > 1 and e3 == 0:
    #                     fundusMucosaElementIdentifiers.append(elementIdentifier)
    #                 elementIdxAround.append(elementIdentifier)
    #                 elementIdentifier += 1
    #                 annotationGroups = annotationGroupsAlong[e2] + annotationGroupsThroughWall[e3]
    #                 if annotationGroups:
    #                     allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
    #                     for annotationGroup in annotationGroups:
    #                         meshGroup = annotationGroup.getMeshGroup(mesh)
    #                         meshGroup.addElement(element)
    #             elementIdxThroughWall.append(elementIdxAround)
    #         elementIdxMat.append(elementIdxThroughWall)
    #
    #     elif 0 < e2 < elementsAlongFundusApexToCardia or \
    #             e2 > elementsAlongFundusApexToCardia + elementsAroundHalfEso - 3:
    #         for e3 in range(elementsCountThroughWall):
    #             elementIdxAround = []
    #             for e1 in range(len(xOuter[e2])):
    #                 eft1 = eftStandard
    #                 scaleFactors = []
    #                 elementtemplate1 = elementtemplateStandard
    #                 bni111 = idxMat[e2][e3][e1]
    #                 bni211 = idxMat[e2][e3][(e1 + 1) % len(idxMat[e2][e3])]
    #                 bni121 = idxMat[e2 + 1][e3][e1]
    #                 bni221 = idxMat[e2 + 1][e3][(e1 + 1) % len(idxMat[e2 + 1][e3])]
    #                 bni112 = idxMat[e2][e3 + 1][e1]
    #                 bni212 = idxMat[e2][e3 + 1][(e1 + 1) % len(idxMat[e2][e3])]
    #                 bni122 = idxMat[e2 + 1][e3 + 1][e1]
    #                 bni222 = idxMat[e2 + 1][e3 + 1][(e1 + 1) % len(idxMat[e2 + 1][e3])]
    #                 nodeIdentifiers = [bni111, bni211, bni121, bni221,
    #                                    bni112, bni212, bni122, bni222]
    #
    #                 if e2 == elementsAlongFundusApexToCardia - 1:
    #                     if e1 == int(0.5 * len(xOuter[e2]) - 2):
    #                         eft1 = eftfactory.createEftNoCrossDerivatives()
    #                         remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
    #                                                [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                         elementtemplateX.defineField(coordinates, -1, eft1)
    #                         elementtemplate1 = elementtemplateX
    #
    #                     elif e1 == int(0.5 * len(xOuter[e2]) - 1):
    #                         scaleFactors = [-1.0]
    #                         eft1 = eftfactory.createEftNoCrossDerivatives()
    #                         setEftScaleFactorIds(eft1, [1], [])
    #                         remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
    #                                                [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
    #                         remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
    #                                                [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                         elementtemplateX.defineField(coordinates, -1, eft1)
    #                         elementtemplate1 = elementtemplateX
    #
    #                     elif e1 == int(0.5 * len(xOuter[e2])):
    #                         scaleFactors = [-1.0]
    #                         eft1 = eftfactory.createEftNoCrossDerivatives()
    #                         setEftScaleFactorIds(eft1, [1], [])
    #                         remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
    #                                                [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
    #                         remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
    #                                                [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #
    #                         elementtemplateX.defineField(coordinates, -1, eft1)
    #                         elementtemplate1 = elementtemplateX
    #
    #                     elif e1 == int(0.5 * len(xOuter[e2]) + 1):
    #                         scaleFactors = [-1.0]
    #                         eft1 = eftfactory.createEftNoCrossDerivatives()
    #                         setEftScaleFactorIds(eft1, [1], [])
    #                         remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
    #                                                [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
    #                         elementtemplateX.defineField(coordinates, -1, eft1)
    #                         elementtemplate1 = elementtemplateX
    #
    #                 if e2 == elementsAlongFundusApexToCardia + elementsAroundHalfEso - 2:
    #                     if e1 == int(0.5 * len(xOuter[e2]) - 2):
    #                         scaleFactors = [-1.0]
    #                         eft1 = eftfactory.createEftNoCrossDerivatives()
    #                         setEftScaleFactorIds(eft1, [1], [])
    #                         remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
    #                                                [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
    #                         elementtemplateX.defineField(coordinates, -1, eft1)
    #                         elementtemplate1 = elementtemplateX
    #
    #                     elif e1 == int(0.5 * len(xOuter[e2]) - 1):
    #                         scaleFactors = [-1.0]
    #                         eft1 = eftfactory.createEftNoCrossDerivatives()
    #                         setEftScaleFactorIds(eft1, [1], [])
    #                         remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
    #                                                [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
    #                         remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
    #                                                [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                         elementtemplateX.defineField(coordinates, -1, eft1)
    #                         elementtemplate1 = elementtemplateX
    #
    #                     elif e1 == int(0.5 * len(xOuter[e2])):
    #                         scaleFactors = [-1.0]
    #                         eft1 = eftfactory.createEftNoCrossDerivatives()
    #                         setEftScaleFactorIds(eft1, [1], [])
    #                         remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
    #                                                [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
    #                         remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
    #                                                [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                         elementtemplateX.defineField(coordinates, -1, eft1)
    #                         elementtemplate1 = elementtemplateX
    #
    #                     elif e1 == int(0.5 * len(xOuter[e2]) + 1):
    #                         eft1 = eftfactory.createEftNoCrossDerivatives()
    #                         remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
    #                                                [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                         elementtemplateX.defineField(coordinates, -1, eft1)
    #                         elementtemplate1 = elementtemplateX
    #
    #                 element = mesh.createElement(elementIdentifier, elementtemplate1)
    #                 element.setNodesByIdentifier(eft1, nodeIdentifiers)
    #                 if scaleFactors:
    #                     element.setScaleFactors(eft1, scaleFactors)
    #                 if e2 < fundusElements and limitingRidge and elementsCountThroughWall > 1 and e3 == 0:
    #                     fundusMucosaElementIdentifiers.append(elementIdentifier)
    #                 elementIdxAround.append(elementIdentifier)
    #                 elementIdentifier += 1
    #                 annotationGroups = annotationGroupsAlong[e2] + annotationGroupsThroughWall[e3]
    #                 if annotationGroups:
    #                     allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
    #                     for annotationGroup in annotationGroups:
    #                         meshGroup = annotationGroup.getMeshGroup(mesh)
    #                         meshGroup.addElement(element)
    #             elementIdxThroughWall.append(elementIdxAround)
    #         elementIdxMat.append(elementIdxThroughWall)
    #
    #     elif elementsAlongFundusApexToCardia <= e2 <= elementsAlongFundusApexToCardia + elementsAroundHalfEso - 3:
    #         for e3 in range(elementsCountThroughWall):
    #             elementIdxAround = []
    #             for e1 in range(elementsCountAroundDuod - 2):
    #                 eft1 = eftStandard
    #                 scaleFactors = []
    #                 elementtemplate1 = elementtemplateStandard
    #                 if e1 > elementsAroundHalfDuod - 2:
    #                     if elementsAlongFundusApexToCardia == e2:  # first ring in eso
    #                         e1IdxBni1 = e1 + 2
    #                         e1IdxBni3 = e1 + 1
    #                     elif e2 == elementsAlongFundusApexToCardia + elementsAroundHalfEso - 3:
    #                         # last ring in esophagus
    #                         e1IdxBni1 = e1 + 1
    #                         e1IdxBni3 = e1 + 2
    #                     elif elementsAlongFundusApexToCardia < e2 < elementsAlongFundusApexToCardia + \
    #                             elementsAroundHalfEso - 3:  # incomplete rings interior of eso
    #                         e1IdxBni1 = e1 + 1
    #                         e1IdxBni3 = e1 + 1
    #                 else:
    #                     e1IdxBni1 = e1
    #                     e1IdxBni3 = e1
    #                 bni1 = idxMat[e2][e3][e1IdxBni1]
    #                 bni2 = idxMat[e2][e3][(e1IdxBni1 + 1) % len(idxMat[e2][e3])]
    #                 bni3 = idxMat[e2 + 1][e3][e1IdxBni3]
    #                 bni4 = idxMat[e2 + 1][e3][(e1IdxBni3 + 1) % len(idxMat[e2 + 1][e3])]
    #                 bni5 = idxMat[e2][e3 + 1][e1IdxBni1]
    #                 bni6 = idxMat[e2][e3 + 1][(e1IdxBni1 + 1) % len(idxMat[e2][e3])]
    #                 bni7 = idxMat[e2 + 1][e3 + 1][e1IdxBni3]
    #                 bni8 = idxMat[e2 + 1][e3 + 1][(e1IdxBni3 + 1) % len(idxMat[e2 + 1][e3])]
    #                 nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
    #
    #                 if e2 == elementsAlongFundusApexToCardia and e1 == elementsAroundHalfDuod - 2:
    #                     scaleFactors = [-1.0]
    #                     eft1 = eftfactory.createEftNoCrossDerivatives()
    #                     setEftScaleFactorIds(eft1, [1], [])
    #                     remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
    #                                            [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
    #                     elementtemplateX.defineField(coordinates, -1, eft1)
    #                     elementtemplate1 = elementtemplateX
    #                 elif e2 == elementsAlongFundusApexToCardia + elementsAroundHalfEso - 3 and \
    #                         e1 == elementsAroundHalfDuod - 2:
    #                     eft1 = eftfactory.createEftNoCrossDerivatives()
    #                     remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
    #                                            [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                     elementtemplateX.defineField(coordinates, -1, eft1)
    #                     elementtemplate1 = elementtemplateX
    #                 elif e2 == elementsAlongFundusApexToCardia + elementsAroundHalfEso - 3 and \
    #                         e1 == elementsAroundHalfDuod - 1:
    #                     scaleFactors = [-1.0]
    #                     eft1 = eftfactory.createEftNoCrossDerivatives()
    #                     setEftScaleFactorIds(eft1, [1], [])
    #                     remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
    #                                            [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
    #                     elementtemplateX.defineField(coordinates, -1, eft1)
    #                     elementtemplate1 = elementtemplateX
    #                 elif e2 == elementsAlongFundusApexToCardia and e1 == elementsAroundHalfDuod - 1:
    #                     eft1 = eftfactory.createEftNoCrossDerivatives()
    #                     remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
    #                                            [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
    #                     elementtemplateX.defineField(coordinates, -1, eft1)
    #                     elementtemplate1 = elementtemplateX
    #
    #                 element = mesh.createElement(elementIdentifier, elementtemplate1)
    #                 element.setNodesByIdentifier(eft1, nodeIdentifiers)
    #                 if e2 < fundusElements and limitingRidge and elementsCountThroughWall > 1 and e3 == 0:
    #                     fundusMucosaElementIdentifiers.append(elementIdentifier)
    #                 elementIdxAround.append(elementIdentifier)
    #                 if scaleFactors:
    #                     element.setScaleFactors(eft1, scaleFactors)
    #                 elementIdentifier += 1
    #                 annotationGroups = annotationGroupsAlong[e2] + annotationGroupsThroughWall[e3]
    #                 if annotationGroups:
    #                     allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
    #                     for annotationGroup in annotationGroups:
    #                         meshGroup = annotationGroup.getMeshGroup(mesh)
    #                         meshGroup.addElement(element)
    #             elementIdxThroughWall.append(elementIdxAround)
    #         elementIdxMat.append(elementIdxThroughWall)
    #
    # # Annulus
    # # Assemble endPoints for annulus
    # endPoints_x = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    # endPoints_d1 = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    # endPoints_d2 = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    # endNode_Id = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    # endDerivativesMap = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    # endProportions = []
    #
    # for n3 in range(elementsCountThroughWall + 1):
    #     n1 = 0
    #     for nAround in range(elementsCountAroundEso):
    #         if nAround == 0:
    #             idx = idxMat[elementsAlongFundusApexToCardia][n3][elementsAroundHalfDuod]
    #         elif 0 < nAround < elementsAroundHalfEso:
    #             idx = idxMat[elementsAlongFundusApexToCardia + n1][n3][elementsAroundHalfDuod - 1]
    #             n1 += 1
    #         elif nAround == elementsAroundHalfEso:
    #             n1 -= 1
    #             idx = idxMat[elementsAlongFundusApexToCardia + n1][n3][elementsAroundHalfDuod]
    #         else:
    #             idx = idxMat[elementsAlongFundusApexToCardia + n1][n3][
    #                 elementsAroundHalfDuod + (1 if (n1 == elementsAroundHalfEso - 2 or n1 == 0) else 0)]
    #             n1 -= 1
    #
    #         endPoints_x[n3][nAround] = xList[idx - stomachStartNode]
    #         endPoints_d1[n3][nAround] = d1List[idx - stomachStartNode]
    #         endPoints_d2[n3][nAround] = d2List[idx - stomachStartNode]
    #         endNode_Id[n3][nAround] = idx
    #
    #         if n3 == elementsCountThroughWall:  # outer layer
    #             endPosition = trackSurfaceStomach.findNearestPosition(endPoints_x[n3][nAround])
    #             endProportions.append(trackSurfaceStomach.getProportion(endPosition))
    #
    # for n3 in range(elementsCountThroughWall + 1):
    #     for nAround in range(elementsCountAroundEso):
    #         if nAround == 0:
    #             endDerivativesMap[n3][nAround] = ((-1, 0, 0), (0, -1, 0), None)
    #         elif nAround == 1:
    #             endDerivativesMap[n3][nAround] = ((-1, 1, 0), (-1, -1, 0), None)
    #         elif 1 < nAround < elementsAroundHalfEso - 1:
    #             endDerivativesMap[n3][nAround] = ((0, 1, 0), (-1, 0, 0), None)
    #         elif nAround == elementsAroundHalfEso - 1:
    #             endDerivativesMap[n3][nAround] = ((1, 1, 0), (-1, 1, 0), None)
    #         elif nAround == elementsAroundHalfEso:
    #             endDerivativesMap[n3][nAround] = (None, None, None)
    #         elif nAround == elementsAroundHalfEso + 1:
    #             endDerivativesMap[n3][nAround] = ((1, -1, 0), (1, 1, 0), None)
    #         elif elementsAroundHalfEso + 1 < nAround < elementsCountAroundEso - 1:
    #             endDerivativesMap[n3][nAround] = ((0, -1, 0), (1, 0, 0), None)
    #         elif nAround == elementsCountAroundEso - 1:
    #             endDerivativesMap[n3][nAround] = ((-1, -1, 0), (1, -1, 0), None)
    #
    # startProportions = []
    # for n in range(elementsCountAroundEso):
    #     startProportions.append(trackSurfaceStomach.getProportion(o1_Positions[n]))
    #
    # cardiaGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region,
    #                                                  get_stomach_term("cardia of stomach"))
    # cardiaMeshGroup = cardiaGroup.getMeshGroup(mesh)
    # if cardiaGroup not in allAnnotationGroups:
    #     allAnnotationGroups.append(cardiaGroup)
    #
    # lastDuodenumElementIdentifier = elementIdentifier
    #
    # stomachWallAnnotationGroups = []
    # if elementsCountThroughWall == 4:
    #     stomachWallAnnotationGroups = [[mucosaGroup], [submucosaGroup], [circularMuscleGroup],
    #                                    [longitudinalMuscleGroup]]
    #
    # # Remove mucosa layer from annulus
    # if elementsCountThroughWall == 4 and limitingRidge:
    #     o1_x = o1_x[1:]
    #     o1_d1 = o1_d1[1:]
    #     o1_d2 = o1_d2[1:]
    #     o1_NodeId = o1_NodeId[1:]
    #     endPoints_x = endPoints_x[1:]
    #     endPoints_d1 = endPoints_d1[1:]
    #     endPoints_d2 = endPoints_d2[1:]
    #     endNode_Id = endNode_Id[1:]
    #     endDerivativesMap = endDerivativesMap[1:]
    #     stomachWallAnnotationGroups = stomachWallAnnotationGroups[1:]
    #
    # nextNodeIdentifier, nextElementIdentifier = createAnnulusMesh3d(
    #     nodes, mesh, nodeIdentifier, elementIdentifier,
    #     o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
    #     endPoints_x, endPoints_d1, endPoints_d2, None, endNode_Id, endDerivativesMap,
    #     elementsCountRadial=elementsCountAcrossCardia, meshGroups=[stomachMeshGroup, cardiaMeshGroup],
    #     wallAnnotationGroups=stomachWallAnnotationGroups,
    #     tracksurface=trackSurfaceStomach,
    #     startProportions=startProportions, endProportions=endProportions,
    #     rescaleStartDerivatives=True, rescaleEndDerivatives=True, sampleBlend=0.0, fixMinimumStart=True,
    #     coordinates=coordinates)
    #
    # elementIdxThroughWall = []
    # n = lastDuodenumElementIdentifier - 1
    # for n3 in range(elementsCountThroughWall):
    #     elementIdxAround = []
    #     for n1 in range(elementsCountAroundEso):
    #         n += 1
    #         elementIdxAround.append(n)
    #     elementIdxThroughWall.append(elementIdxAround)
    # elementIdxMat.append(elementIdxThroughWall)
    #
    # # delete mucosa layer in fundus when there is a limiting ridge
    # mesh_destroy_elements_and_nodes_by_identifiers(mesh, fundusMucosaElementIdentifiers)
    #
    # # Create annotation groups for dorsal and ventral parts of the stomach
    # dorsalGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("dorsal stomach"))
    # ventralGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("ventral stomach"))
    # dorsalMeshGroup = dorsalGroup.getMeshGroup(mesh)
    # ventralMeshGroup = ventralGroup.getMeshGroup(mesh)
    #
    # for e2 in range(len(elementIdxMat)):
    #     for e3 in range(len(elementIdxMat[e2])):
    #         for e1 in range(len(elementIdxMat[e2][e3])):
    #             elementIdx = elementIdxMat[e2][e3][e1]
    #             element = mesh.findElementByIdentifier(elementIdx)
    #             if e1 < 0.5 * len(elementIdxMat[e2][e3]):
    #                 ventralMeshGroup.addElement(element)
    #             else:
    #                 dorsalMeshGroup.addElement(element)
    # if dorsalGroup not in allAnnotationGroups:
    #     allAnnotationGroups.append(dorsalGroup)
    # if ventralGroup not in allAnnotationGroups:
    #     allAnnotationGroups.append(ventralGroup)
    #
    # nodesOnLCMargin = []
    # for n2 in range(elementsAlongEsophagus + 1):
    #     for n3 in range(elementsThroughEsophagusWall + 1):
    #         nodeIdxOnLCMargin = 1 + elementsAroundHalfEso + \
    #                             n2 * (elementsThroughEsophagusWall + 1) * elementsCountAroundEso + \
    #                             n3 * elementsCountAroundEso
    #         nodesOnLCMargin.append(nodeIdxOnLCMargin)
    # allNodesOnLC = nodesOnLCMargin + nodeIdxLC
    #
    # nearLCGroup = AnnotationGroup(region, ("elements adjacent to lesser curvature", "None"))
    #
    # elementIter = mesh.createElementiterator()
    # element = elementIter.next()
    # while element.isValid():
    #     eft = element.getElementfieldtemplate(coordinates, -1)
    #     nodeIdentifiers = get_element_node_identifiers(element, eft)
    #     for n in range(len(nodeIdentifiers)):
    #         if nodeIdentifiers[n] in allNodesOnLC:
    #             nearLCGroup.getMeshGroup(mesh).addElement(element)
    #             break
    #     element = elementIter.next()
    # allAnnotationGroups.append(nearLCGroup)
    #
    # # Create split coordinate field
    # if splitCoordinates:
    #     nodesOnSplitMargin = []
    #
    #     for n2 in range(elementsAlongEsophagus + 1):
    #         for n3 in range(elementsThroughEsophagusWall + 1):
    #             nodeIdxOnGCMargin = 1 + n2 * (elementsThroughEsophagusWall + 1) * elementsCountAroundEso + \
    #                                 n3 * elementsCountAroundEso
    #             nodesOnSplitMargin.append(nodeIdxOnGCMargin)
    #             nodeIdxOnLCMargin = 1 + elementsAroundHalfEso + \
    #                                 n2 * (elementsThroughEsophagusWall + 1) * elementsCountAroundEso + \
    #                                 n3 * elementsCountAroundEso
    #             nodesOnSplitMargin.append(nodeIdxOnLCMargin)
    #
    #     nodesOnSplitMargin += nodeIdxGC + nodeIdxLC
    #
    #     splitCoordinates = findOrCreateFieldCoordinates(fm, name="split coordinates")
    #     splitNodetemplate1 = nodes.createNodetemplate()
    #     splitNodetemplate1.defineField(splitCoordinates)
    #     splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    #     splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    #     splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    #     if useCrossDerivatives:
    #         splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
    #     if useCubicHermiteThroughWall:
    #         splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
    #         if useCrossDerivatives:
    #             splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
    #             splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
    #             splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
    #
    #     splitNodetemplate2 = nodes.createNodetemplate()
    #     splitNodetemplate2.defineField(splitCoordinates)
    #     splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_VALUE, 2)
    #     splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
    #     splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS2, 2)
    #     if useCrossDerivatives:
    #         splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 2)
    #     if useCubicHermiteThroughWall:
    #         splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS3, 2)
    #         if useCrossDerivatives:
    #             splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 2)
    #             splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 2)
    #             splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 2)
    #
    #     nodeIter = nodes.createNodeiterator()
    #     node = nodeIter.next()
    #     while node.isValid():
    #         # if not markerPoints.containsNode(node):
    #         cache.setNode(node)
    #         identifier = node.getIdentifier()
    #         marginNode = identifier in nodesOnSplitMargin
    #         x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)[1]
    #         d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)[1]
    #         d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)[1]
    #         if useCrossDerivatives:
    #             d1d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)[1]
    #         if useCubicHermiteThroughWall:
    #             d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)[1]
    #             if useCrossDerivatives:
    #                 d1d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)[1]
    #                 d2d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, 3)[1]
    #                 d1d2d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, 3)[1]
    #
    #         node.merge(splitNodetemplate2 if marginNode else splitNodetemplate1)
    #         versionCount = 2 if marginNode else 1
    #         for vn in range(versionCount):
    #             splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, vn + 1, x)
    #             splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, vn + 1, d1)
    #             splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, vn + 1, d2)
    #             if useCrossDerivatives:
    #                 splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, vn + 1, d1d2)
    #             if useCubicHermiteThroughWall:
    #                 splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, vn + 1, d3)
    #                 if useCrossDerivatives:
    #                     splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, vn + 1, d1d3)
    #                     splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, vn + 1, d2d3)
    #                     splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, vn + 1, d1d2d3)
    #         node = nodeIter.next()
    #
    #     elementIter = mesh.createElementiterator()
    #     element = elementIter.next()
    #     splitElementtemplate1 = mesh.createElementtemplate()
    #     splitElementtemplate2 = mesh.createElementtemplate()
    #
    #     count = 0
    #     elementsInOstium = elementsCountAroundEso * elementsAlongEsophagus * elementsThroughEsophagusWall
    #     closedLoopElementId = nextElementIdentifier - elementsCountAroundEso * elementsCountAcrossCardia - \
    #                           elementsCountAroundDuod * elementsCountThroughWall * (elementsAlongCardiaToDuod + 1)
    #
    #     allValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
    #                       Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3,
    #                       Node.VALUE_LABEL_D2_DS2DS3, Node.VALUE_LABEL_D3_DS1DS2DS3]
    #
    #     while element.isValid():
    #         eft = element.getElementfieldtemplate(coordinates, -1)
    #         nodeIdentifiers = get_element_node_identifiers(element, eft)
    #         elementId = element.getIdentifier()
    #         marginDorsal = False
    #         for n in range(len(nodeIdentifiers)):
    #             marginElement = nodeIdentifiers[n] in nodesOnSplitMargin
    #             if marginElement:
    #                 count += 1
    #                 if count < 3 and (elementId <= elementsInOstium or elementId > closedLoopElementId):
    #                     marginDorsal = True
    #                 elif count >= 3 and (elementId <= elementsInOstium or elementId > closedLoopElementId):
    #                     if count == 4:
    #                         count = 0
    #                 elif elementsInOstium < elementId < elementsInOstium + len(xOuter[0]) + 1:
    #                     marginDorsal = True
    #                     count = 0
    #                 elif elementsInOstium + len(xOuter[0]) < elementId < elementsInOstium + len(xOuter[0]) * 2 + 1:
    #                     count = 0
    #                 elif count < 2 and elementId > elementsInOstium + 2 * (len(xOuter[0])):
    #                     marginDorsal = True
    #                 elif count >= 2 and elementId > elementsInOstium + 2 * (len(xOuter[0])):
    #                     if count == 2:
    #                         count = 0
    #                 break
    #
    #         if marginDorsal:
    #             # Find nodes on margin to remap with version 2
    #             lnRemapV2 = []
    #             for n in range(len(nodeIdentifiers)):
    #                 if nodeIdentifiers[n] in nodesOnSplitMargin:
    #                     lnRemapV2.append(n + 1)
    #             eft2 = eft
    #             remapEftNodeValueLabelsVersion(eft2, lnRemapV2, allValueLabels, 2)
    #
    #             splitElementtemplate2.defineField(splitCoordinates, -1, eft2)
    #             element.merge(splitElementtemplate2)
    #             element.setNodesByIdentifier(eft2, nodeIdentifiers)
    #             if eft2.getNumberOfLocalScaleFactors() > 0:
    #                 element.setScaleFactor(eft2, 1, -1.0)
    #         else:
    #             splitElementtemplate1.defineField(splitCoordinates, -1, eft)
    #             element.merge(splitElementtemplate1)
    #             element.setNodesByIdentifier(eft, nodeIdentifiers)
    #             if eft.getNumberOfLocalScaleFactors() == 1:
    #                 element.setScaleFactors(eft, [-1.0])
    #
    #         element = elementIter.next()
    #
    #     allAnnotationGroups.append(nearLCGroup)

    nextNodeIdentifier = nodeIdentifier
    nextElementIdentifier = elementIdentifier

    return allAnnotationGroups, nextNodeIdentifier, nextElementIdentifier


def findCentreOnCentralPathFromCrossAxisEndPt(xPoint, xCentralPath, dCentralPath):
    """
    Find point on central path which intersects with cross axis endpoint.
    :param xPoint: cross axis endpoint
    :param xCentralPath: List of points on central path
    :param dCentralPath: List of cross axis derivatives along central path
    :return xPtOnCP: intersection point on central path
    """
    arcLength = 0.0
    for n in range(len(xCentralPath) - 1):
        arcLength += interp.getCubicHermiteArcLength(xCentralPath[n], dCentralPath[n],
                                                     xCentralPath[n + 1], dCentralPath[n + 1])
    arcLow = 0.0
    arcUp = arcLength

    for iter in range(100):
        thetaLow = findThetaFromArcDistance(xPoint, arcLow, xCentralPath, dCentralPath)
        thetaUp = findThetaFromArcDistance(xPoint, arcUp, xCentralPath, dCentralPath)
        arcDistance = (arcLow + arcUp) * 0.5

        if abs(thetaLow - thetaUp) < 1e-06:
            xPtOnCP = interp.getCubicHermiteCurvesPointAtArcDistance(xCentralPath, dCentralPath, arcDistance)[0]
            break
        elif thetaLow > thetaUp:
            arcLow = arcDistance
        elif thetaLow < thetaUp:
            arcUp = arcDistance

    if iter > 99:
        print('Search for findCentreOnCentralPathFromCrossAxisEndPt - Max iters reached:', iter)

    return xPtOnCP


def findThetaFromArcDistance(xPt, arcDistance, xCentralPath, dCentralPath):
    """
    Find angle between a vector between a guess point on the central path to a cross axis endpoint and the cross-axis
    of the central path.
    :param xPt: Cross axis endpoint
    :param arcDistance: Arclength along central path where guess point sits
    :param xCentralPath: List of points along central path
    :param dCentralPath: List of cross axis derivatives along central path
    :return theta: angle between vectors
    """
    xGuess, dGuess = interp.getCubicHermiteCurvesPointAtArcDistance(xCentralPath, dCentralPath, arcDistance)[0:2]
    if xGuess == xPt:
        return 0.0
    else:
        dxPt = findDerivativeBetweenPoints(xGuess, xPt)
        cosTheta = vector.dotproduct(dGuess, dxPt) / (vector.magnitude(dGuess) * vector.magnitude(dxPt))
        if abs(cosTheta) > 1.0:
            cosTheta = 1.0 if cosTheta > 0.0 else -1.0
        theta = abs(math.pi * 0.5 - math.acos(cosTheta))
        return theta
