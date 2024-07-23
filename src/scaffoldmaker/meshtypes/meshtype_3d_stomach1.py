"""
Generates a 3-D stomach mesh along the central line, with variable
numbers of elements around esophagus and duodenum, along and through
wall, with variable radius and thickness along.
"""

from __future__ import division

import copy
import math

from cmlibs.maths.vectorops import cross, sub, set_magnitude, normalize, magnitude, axis_angle_to_rotation_matrix
from cmlibs.utils.zinc.field import find_or_create_field_coordinates
from cmlibs.utils.zinc.finiteelement import get_element_node_identifiers, get_maximum_element_identifier, \
    get_maximum_node_identifier
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, mergeAnnotationGroups, \
    getAnnotationGroupForTerm, findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.esophagus_terms import get_esophagus_term
from scaffoldmaker.annotation.smallintestine_terms import get_smallintestine_term
from scaffoldmaker.annotation.stomach_terms import get_stomach_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.meshtype_3d_ostium2 import generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import setEftScaleFactorIds, remapEftNodeValueLabel
from scaffoldmaker.utils.geometry import sampleEllipsePoints
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, \
    mesh_destroy_elements_and_nodes_by_identifiers, get_nodeset_path_ordered_field_parameters, \
    get_nodeset_path_field_parameters

def getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName):
    assert parameterSetName in cls.getParameterSetNames()  # make sure parameter set is in list of parameters of parent scaffold
    if parameterSetName in ("Default", "Human 1"):
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-5-6-7-8-3-9-10-11-12-13-14-15"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[0.000,0.560,0.000], [0.000,-0.070,0.000], [0.060,0.000,0.000], [0.170,0.120,-0.000], [0.000,0.000,0.060], [0.000,0.000,0.110]]),
                (2, [[0.000,0.390,0.000], [0.000,-0.280,0.000], [0.130,0.000,0.000], [0.020,-0.120,-0.000], [0.000,0.000,0.130], [0.000,0.000,0.170]]),
                (3, [[0.000,0.000,0.000], [[-0.200,0.000,0.000],[0.000,-0.490,0.000]], [[0.000,-0.400,0.000],[0.130,0.000,0.000]], [[-0.060,0.010,0.000],[-0.010,-0.400,0.000]], [[0.000,0.000,0.430],[0.000,0.000,0.130]], [[0.000,0.000,0.010],[0.000,0.000,0.010]]]),
                (4, [[0.490,0.060,0.000], [-0.030,-0.000,-0.000], [0.010,-0.040,0.000], [0.020,-0.140,0.000], [0.000,0.000,0.040], [0.000,0.000,0.130]]),
                (5, [[0.450,0.050,-0.000], [-0.060,-0.010,-0.000], [0.030,-0.170,0.000], [0.020,-0.120,0.000], [0.000,0.000,0.170], [0.000,0.000,0.130]]),
                (6, [[0.380,0.040,0.000], [-0.080,-0.020,0.000], [0.040,-0.270,0.000], [0.010,-0.100,0.000], [0.000,0.000,0.290], [0.000,0.000,0.110]]),
                (7, [[0.280,0.020,0.000], [-0.110,-0.020,0.000], [0.050,-0.360,0.000], [0.000,-0.070,0.000], [0.000,0.000,0.380], [0.000,0.000,0.070]]),
                (8, [[0.160,0.010,0.000], [-0.140,-0.010,0.000], [0.040,-0.400,0.000], [-0.030,-0.020,0.000], [0.000,0.000,0.420], [0.000,0.000,0.030]]),
                (9, [[-0.230,0.020,0.000], [-0.230,0.050,0.000], [-0.090,-0.370,0.000], [-0.080,0.050,0.000], [0.000,0.000,0.440], [0.000,0.000,-0.010]]),
                (10, [[-0.450,0.100,0.000], [-0.200,0.120,0.000], [-0.180,-0.290,0.000], [-0.060,0.110,0.000], [0.000,0.000,0.410], [0.000,0.000,-0.050]]),
                (11, [[-0.610,0.260,0.000], [-0.150,0.200,0.000], [-0.220,-0.160,0.000], [-0.020,0.130,0.000], [0.000,0.000,0.340], [0.000,0.000,-0.090]]),
                (12, [[-0.730,0.500,0.000], [-0.040,0.230,0.000], [-0.210,-0.040,0.000], [0.040,0.090,0.000], [0.000,0.000,0.230], [0.000,0.000,-0.070]]),
                (13, [[-0.710,0.690,0.000], [0.050,0.170,0.000], [-0.140,0.030,0.000], [0.060,0.050,-0.000], [0.000,0.000,0.180], [0.000,0.000,-0.060]]),
                (14, [[-0.640,0.830,0.000], [0.100,0.140,0.000], [-0.080,0.060,0.000], [0.030,0.030,-0.000], [0.000,0.000,0.110], [0.000,0.000,-0.030]]),
                (15, [[-0.510,0.970,0.000], [0.160,0.140,0.000], [-0.080,0.090,0.000], [-0.030,0.030,0.000], [0.000,0.000,0.130], [0.000,0.000,0.070]])
                ]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-14',
                    'name': get_stomach_term('stomach')[0],
                    'ontId': get_stomach_term('stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_stomach_term('esophagus part of stomach')[0],
                    'ontId': get_stomach_term('esophagus part of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-7',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '8-10',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-12',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '13',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '14',
                    'name': get_stomach_term('duodenum part of stomach')[0],
                    'ontId': get_stomach_term('duodenum part of stomach')[1]
                }]
        })
    elif "Human 2" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
        'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-5-6-7-3-8-9-10-11-12-13-14-15"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[11.750,-111.874,1127.887], [7.636,-5.715,-7.930], [5.678,1.265,4.556], [-8.397,13.092,24.878], [-0.708,-3.530,1.862], [-0.807,-7.995,7.596]]),
                (2, [[23.789,-117.922,1120.040], [26.354,-6.724,-6.404], [4.223,6.205,10.864], [10.037,1.800,8.968], [-1.192,-11.215,6.869], [-2.926,-13.889,10.204]]),
                (3, [[63.704,-120.094,1123.374], [[0.500,-9.138,-13.405],[50.106,1.267,11.056]], [[37.742,-3.477,3.778],[-2.509,7.605,10.499]], [[3.190,-0.590,-0.290],[3.190,-0.590,-0.290]], [[-5.452,-34.121,23.056],[-1.379,-10.790,7.486]], [[-0.100,-1.220,0.550],[-0.100,-1.220,0.550]]]),
                (4, [[61.247,-99.931,1152.681], [0.346,-2.728,-3.873], [11.320,-0.365,1.269], [15.701,-3.029,2.832], [-0.653,-5.931,4.119], [-5.756,-22.456,14.946]]),
                (5, [[61.743,-103.510,1147.760], [0.413,-3.592,-5.311], [24.159,-2.387,3.493], [10.039,-1.491,0.928], [-2.982,-15.339,10.142], [-2.384,-12.784,8.034]]),
                (6, [[62.381,-107.527,1141.785], [0.249,-4.737,-7.073], [30.559,-2.973,3.067], [5.720,-0.448,-0.071], [-3.839,-23.420,15.550], [-0.206,-5.098,2.599]]),
                (7, [[62.800,-113.150,1133.665], [0.116,-6.546,-9.651], [35.408,-3.630,2.888], [4.541,-0.488,0.325], [-4.677,-29.659,20.061], [-0.383,-3.012,1.219]]),
                (8, [[64.339,-131.197,1107.233], [0.086,-11.682,-16.915], [39.201,-3.705,2.758], [0.490,-2.470,-4.030], [-5.108,-35.712,24.638], [0.260,0.760,0.250]]),
                (9, [[62.912,-143.954,1088.811], [-5.216,-12.408,-17.967], [34.623,-8.161,-4.415], [-4.950,-4.940,-8.800], [-4.917,-34.532,25.275], [0.770,0.610,1.050]]),
                (10, [[53.361,-155.397,1072.006], [-15.833,-9.775,-15.486], [25.117,-13.916,-16.896], [-14.190,-3.410,-7.580], [-2.355,-30.712,21.794], [2.030,4.350,-1.240]]),
                (11, [[32.110,-162.230,1059.680], [-22.173,-3.681,-8.142], [10.067,-16.126,-20.126], [-15.400,0.330,-0.620], [-2.559,-23.629,17.653], [1.200,8.620,-5.080]]),
                (12, [[10.560,-162.970,1055.650], [-20.956,2.164,-0.908], [-0.696,-13.976,-17.241], [-8.750,2.830,4.760], [-2.502,-18.048,14.732], [-1.000,7.490,-3.260]]),
                (13, [[-8.740,-158.280,1057.630], [-17.269,8.027,2.588], [-6.253,-10.164,-10.200], [-2.630,4.430,6.290], [-3.442,-11.913,13.981], [-1.380,7.080,-2.660]]),
                (14, [[-23.260,-147.620,1060.640], [-11.787,11.369,1.720], [-5.807,-5.321,-4.621], [-0.700,3.020,1.710], [-3.207,-4.764,9.515], [-2.080,5.830,-2.460]]),
                (15, [[-32.281,-136.261,1061.249], [-6.125,11.114,-0.493], [-5.735,-3.374,-4.823], [-2.070,-0.220,-4.520], [-4.630,-2.238,7.071], [-4.400,2.440,0.570]]) 
                ]]),
                
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-14',
                    'name': get_stomach_term('stomach')[0],
                    'ontId': get_stomach_term('stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_stomach_term('esophagus part of stomach')[0],
                    'ontId': get_stomach_term('esophagus part of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-6',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7-10',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-12',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '13',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '14',
                    'name': get_stomach_term('duodenum part of stomach')[0],
                    'ontId': get_stomach_term('duodenum part of stomach')[1]
                }]
        })
    elif "Mouse 1" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-5-6-7-8-9-3-10-11-12-13-14"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[-0.047,0.617,-0.000], [0.043,-0.123,0.000], [0.065,0.011,0.000], [0.098,-0.016,0.000], [0.000,0.000,0.066], [0.000,0.000,0.124]]),
                (2, [[-0.010,0.400,-0.000], [0.030,-0.309,0.000], [0.109,0.019,0.000], [0.082,-0.024,0.000], [0.000,-0.000,0.110], [0.000,0.000,0.116]]),
                (3, [[0.000,0.000,0.000], [[-0.288,-0.009,0.000],[-0.010,-0.490,0.000]], [[0.010,-0.400,0.000],[0.109,0.019,0.000]], [[-0.194,0.007,0.000],[-0.194,0.007,0.000]], [[0.000,0.000,0.380],[0.000,-0.000,0.110]], [[0.000,0.000,-0.015],[0.000,0.000,-0.015]]]),
                (4, [[0.540,0.710,0.000], [-0.005,-0.065,0.000], [0.080,-0.010,0.000], [0.098,-0.016,0.000], [0.000,0.000,0.040], [0.000,0.000,0.124]]),
                (5, [[0.530,0.630,0.000], [-0.015,-0.095,0.000], [0.170,-0.030,0.000], [0.082,-0.024,0.000], [0.000,0.000,0.160], [0.000,0.000,0.116]]),
                (6, [[0.510,0.520,0.000], [-0.029,-0.135,0.000], [0.240,-0.060,0.000], [0.066,-0.042,0.000], [0.000,0.000,0.270], [0.000,0.000,0.098]]),
                (7, [[0.470,0.360,0.000], [-0.055,-0.161,0.000], [0.300,-0.120,0.000], [0.026,-0.089,0.000], [0.000,0.000,0.350], [0.000,0.000,0.056]]),
                (8, [[0.400,0.200,0.000], [-0.107,-0.145,0.000], [0.290,-0.240,0.000], [-0.054,-0.110,0.000], [0.000,0.000,0.380], [0.000,0.000,0.020]]),
                (9, [[0.260,0.080,0.000], [-0.202,-0.111,0.000], [0.190,-0.340,0.000], [-0.132,-0.084,0.000], [0.000,0.000,0.390], [0.000,0.000,0.002]]),
                (10, [[-0.290,0.070,0.000], [-0.237,0.128,0.000], [-0.200,-0.320,0.000], [-0.109,0.097,0.000], [0.000,0.000,0.360], [0.000,0.000,-0.031]]),
                (11, [[-0.460,0.230,0.000], [-0.099,0.191,0.000], [-0.230,-0.210,0.000], [0.020,0.121,0.000], [0.000,0.000,0.320], [0.000,0.000,-0.079]]),
                (12, [[-0.486,0.419,0.000], [-0.008,0.167,0.000], [-0.170,-0.080,0.000], [0.069,0.087,0.000], [0.000,0.000,0.210], [0.000,0.000,-0.109]]),
                (13, [[-0.480,0.560,0.000], [-0.012,0.142,0.000], [-0.094,-0.025,0.000], [0.030,0.020,0.000], [0.000,0.000,0.102], [0.000,0.000,-0.046]]),
                (14, [[-0.510,0.700,0.000], [-0.048,0.137,0.000], [-0.110,-0.040,0.000], [-0.062,-0.050,0.000], [0.000,0.000,0.120], [0.000,0.000,0.082]])
                ]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-13',
                    'name': get_stomach_term('stomach')[0],
                    'ontId': get_stomach_term('stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_stomach_term('esophagus part of stomach')[0],
                    'ontId': get_stomach_term('esophagus part of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-8',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9-10',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '12',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '13',
                    'name': get_stomach_term('duodenum part of stomach')[0],
                    'ontId': get_stomach_term('duodenum part of stomach')[1]
                }]
        })
    elif "Pig 1" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-5-6-7-3-8-9-10-11-12-13"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[0.040,0.550,0.000], [0.000,-0.040,-0.000], [0.055,-0.007,0.000], [-0.010,-0.090,0.000], [0.000,0.000,0.0554], [0.000,0.000,0.060]]),
                (2, [[0.030,0.430,0.000], [-0.020,-0.270,-0.000], [0.118,-0.005,-0.000], [-0.000,-0.100,0.000], [-0.000,0.000,0.118], [0.000,0.000,0.080]]),
                (3, [[0.000,0.000,0.000], [[-0.210,0.000,0.000],[-0.040,-0.590,-0.000]], [[0.000,-0.431,0.000],[0.118,-0.005,-0.000]], [[-0.020,-0.020,0.000],[-0.020,-0.020,0.000]], [[0.000,0.000,0.400],[-0.000,0.000,0.118]], [[0.000,0.000,0.010],[0.000,0.000,0.010]]]),
                (4, [[0.440,0.010,0.000], [-0.010,-0.000,0.000], [0.010,-0.050,0.000], [-0.010,-0.090,0.000], [0.000,0.000,0.080], [0.000,0.000,0.060]]),
                (5, [[0.420,0.010,0.000], [-0.050,-0.010,0.000], [-0.000,-0.150,0.000], [-0.000,-0.100,0.000], [0.000,0.000,0.150], [0.000,0.000,0.080]]),
                (6, [[0.330,0.000,0.000], [-0.110,-0.000,0.000], [0.000,-0.280,0.000], [-0.000,-0.130,0.000], [0.000,0.000,0.280], [0.000,0.000,0.120]]),
                (7, [[0.190,0.000,0.000], [-0.170,-0.010,0.000], [0.010,-0.390,0.000], [-0.020,-0.080,0.000], [0.000,0.000,0.380], [0.000,0.000,0.070]]),
                (8, [[-0.220,0.010,0.000], [-0.250,0.030,0.000], [-0.030,-0.430,0.000], [-0.050,0.030,0.000], [0.000,0.000,0.400], [0.000,0.000,-0.010]]),
                (9, [[-0.500,0.050,0.000], [-0.260,0.120,0.000], [-0.160,-0.360,0.000], [-0.120,0.140,0.000], [0.000,0.000,0.370], [0.000,0.000,-0.050]]),
                (10, [[-0.700,0.240,0.000], [-0.110,0.270,0.000], [-0.280,-0.160,0.000], [-0.050,0.210,0.000], [0.000,0.000,0.310], [0.000,0.000,-0.090]]),
                (11, [[-0.700,0.520,0.000], [0.120,0.260,0.000], [-0.270,0.050,0.000], [0.100,0.090,0.000], [0.000,0.000,0.180], [0.000,0.000,-0.110]]),
                (12, [[-0.500,0.700,0.000], [0.140,0.190,0.000], [-0.090,0.030,0.000], [0.080,-0.020,0.000], [0.000,0.000,0.090], [0.000,0.000,-0.040]]),
                (13, [[-0.410,0.880,0.000], [0.030,0.160,0.000], [-0.090,0.010,0.000], [-0.080,-0.030,0.000], [0.000,0.000,0.090], [0.000,0.000,0.040]])
                ]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-12',
                    'name': get_stomach_term('stomach')[0],
                    'ontId': get_stomach_term('stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_stomach_term('esophagus part of stomach')[0],
                    'ontId': get_stomach_term('esophagus part of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-6',
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
                    'name': get_stomach_term('duodenum part of stomach')[0],
                    'ontId': get_stomach_term('duodenum part of stomach')[1]
                }]
        })
    elif "Rat 1" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                    "Structure": "1-2-3.2, 4-5-6-7-8-9-3-10-11-12-13-14-15"
                },
                'meshEdits': exnode_string_from_nodeset_field_parameters(
                    ['coordinates'],
                    [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[-0.011,0.636,0.000], [0.012,-0.116,0.000], [0.057,0.000,0.000], [0.122,-0.056,0.000], [-0.000,0.000,0.057], [0.000,0.000,0.076]]),
                    (2, [[-0.001,0.419,-0.000], [0.008,-0.318,0.000], [0.108,0.000,0.000], [0.101,-0.066,0.000], [0.000,0.000,0.108], [0.000,0.000,0.084]]),
                    (3, [[0.000,0.000,0.000], [[-0.200,0.002,0.000],[-0.006,-0.520,0.000]], [[0.005,-0.420,0.000],[0.108,0.000,0.000]], [[-0.156,0.003,0.000],[-0.156,0.003,0.000]], [[0.000,0.000,0.360],[0.000,0.000,0.108]], [[0.000,0.000,-0.009],[0.000,0.000,-0.009]]]),
                    (4, [[0.601,0.602,0.000], [-0.027,-0.056,0.000], [0.077,-0.036,0.000], [0.122,-0.056,0.000], [0.000,0.000,0.080], [0.000,0.000,0.076]]),
                    (5, [[0.564,0.524,0.000], [-0.047,-0.100,-0.000], [0.189,-0.097,0.000], [0.101,-0.066,0.000], [0.000,0.000,0.160], [0.000,0.000,0.084]]),
                    (6, [[0.507,0.402,0.000], [-0.071,-0.139,-0.000], [0.273,-0.171,0.000], [0.062,-0.070,-0.000], [0.000,0.000,0.250], [0.000,0.000,0.073]]),
                    (7, [[0.420,0.248,0.000], [-0.097,-0.137,-0.000], [0.307,-0.237,0.000], [-0.011,-0.067,-0.000], [0.000,0.000,0.300], [0.000,0.000,0.039]]),
                    (8, [[0.315,0.129,0.000], [-0.125,-0.109,-0.000], [0.256,-0.304,0.000], [-0.080,-0.076,0.000], [0.000,0.000,0.330], [0.000,0.000,0.030]]),
                    (9, [[0.171,0.034,0.000], [-0.161,-0.066,0.000], [0.144,-0.389,0.000], [-0.125,-0.058,0.000], [0.000,0.000,0.360], [0.000,0.000,0.015]]),
                    (10, [[-0.218,0.048,0.000], [-0.208,0.094,0.000], [-0.173,-0.374,0.000], [-0.152,0.079,0.000], [0.000,0.000,0.340], [0.000,0.000,-0.015]]),
                    (11, [[-0.404,0.184,0.000], [-0.142,0.162,0.000], [-0.299,-0.260,0.000], [-0.044,0.113,0.000], [0.000,0.000,0.330], [0.000,0.000,-0.025]]),
                    (12, [[-0.497,0.356,0.000], [-0.049,0.209,0.000], [-0.255,-0.188,0.000], [0.077,0.107,0.000], [0.000,0.000,0.290], [0.000,0.000,-0.111]]),
                    (13, [[-0.490,0.587,0.000], [-0.018,0.189,-0.000], [-0.152,-0.045,0.000], [0.069,0.049,0.000], [0.000,0.000,0.120], [0.000,0.000,-0.073]]),
                    (14, [[-0.523,0.730,0.000], [-0.032,0.116,0.000], [-0.111,-0.036,0.000], [0.003,-0.002,0.000], [0.000,0.000,0.120], [0.000,0.000,0.018]]),
                    (15, [[-0.552,0.820,0.000], [-0.026,0.063,0.000], [-0.132,-0.045,0.000], [-0.045,-0.016,0.000], [0.000,0.000,0.150], [0.000,0.000,0.042]])
                    ]]),

                'userAnnotationGroups': [
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '1-14',
                        'name': get_stomach_term('stomach')[0],
                        'ontId': get_stomach_term('stomach')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '1-2',
                        'name': get_stomach_term('esophagus part of stomach')[0],
                        'ontId': get_stomach_term('esophagus part of stomach')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '3-8',
                        'name': get_stomach_term('fundus of stomach')[0],
                        'ontId': get_stomach_term('fundus of stomach')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '9-11',
                        'name': get_stomach_term('body of stomach')[0],
                        'ontId': get_stomach_term('body of stomach')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '12',
                        'name': get_stomach_term('pyloric antrum')[0],
                        'ontId': get_stomach_term('pyloric antrum')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '13',
                        'name': get_stomach_term('pyloric canal')[0],
                        'ontId': get_stomach_term('pyloric canal')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '14',
                        'name': get_stomach_term('duodenum part of stomach')[0],
                        'ontId': get_stomach_term('duodenum part of stomach')[1]
                    }]
            })
    elif "Material" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-5-6-7-8-9-3-10-11-12-13-14-15-16-17-18"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[0.000,0.800,0.000], [0.000,-0.200,0.000], [0.080,-0.000,0.000], [0.000,-0.100,0.000], [0.000,0.000,0.080], [0.000,0.000,0.100]]),
                (2, [[0.000,0.500,0.000], [0.000,-0.400,0.000], [0.150,0.000,0.000], [0.000,-0.110,0.000], [0.000,0.000,0.150], [0.000,0.000,0.110]]),
                (3, [[0.000,0.000,0.000], [[-0.190,0.000,0.000],[0.000,-0.600,0.000]], [[0.000,-0.500,0.000], [0.210,0.000,0.000]], [[0.000,0.000,0.000],[0.000,0.000,0.000]], [[0.000,0.000,0.500],[0.000,0.000,0.200]], [[0.000,0.000,0.000],[0.000,0.000,0.000]]]),
                (4, [[0.700,0.000,0.000], [-0.020,0.000,0.000], [0.000,-0.050,0.000], [0.000,-0.100,0.000], [0.000,0.000,0.050], [0.000,0.000,0.100]]),
                (5, [[0.680,0.000,0.000], [-0.050,0.000,0.000], [0.000,-0.150,0.000], [0.000,-0.110,0.000], [0.000,0.000,0.150], [0.000,0.000,0.110]]),
                (6, [[0.600,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.290,0.000], [0.000,-0.130,0.000], [0.000,0.000,0.290], [0.000,0.000,0.130]]),
                (7, [[0.490,0.000,0.000], [-0.130,0.000,0.000], [0.000,-0.400,0.000], [0.000,-0.100,0.000], [0.000,0.000,0.400], [0.000,0.000,0.100]]),
                (8, [[0.350,0.000,0.000], [-0.160,0.000,0.000], [0.000,-0.480,0.000], [0.000,-0.050,0.000], [0.000,0.000,0.480], [0.000,0.000,0.050]]),
                (9, [[0.180,0.000,0.000], [-0.170,0.000,0.000], [0.000,-0.500,0.000], [0.000,-0.010,0.000], [0.000,0.000,0.500], [0.000,0.000,0.010]]),
                (10, [[-0.200,0.000,0.000], [-0.200,0.000,0.000], [0.000,-0.500,0.000], [0.000,0.000,0.000], [0.000,0.000,0.500], [0.000,0.000,0.000]]),
                (11, [[-0.400,0.000,0.000], [-0.200,0.000,0.000], [0.000,-0.500,0.000], [0.000,0.000,0.000], [0.000,0.000,0.500], [0.000,0.000,0.000]]),
                (12, [[-0.600,0.000,0.000], [-0.180,0.000,0.000], [0.000,-0.500,0.000], [0.000,0.000,0.000], [0.000,0.000,0.500], [0.000,0.000,0.000]]),
                (13, [[-0.750,0.000,0.000], [-0.150,0.000,0.000], [0.000,-0.440,0.000], [0.000,0.100,0.000], [0.000,0.000,0.440], [0.000,0.000,-0.100]]),
                (14, [[-0.900,0.000,0.000], [-0.120,0.000,0.000], [0.000,-0.310,0.000], [0.000,0.120,0.000], [0.000,0.000,0.310], [0.000,0.000,-0.120]]),
                (15, [[-1.000,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.200,0.000], [0.000,0.000,0.000], [0.000,0.000,0.200], [0.000,0.000,0.000]]),
                (16, [[-1.100,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.200,0.000], [0.000,0.000,0.000], [0.000,0.000,0.200], [0.000,0.000,0.000]]),
                (17, [[-1.200,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.200,0.000], [0.000,0.000,0.000], [0.000,0.000,0.200], [0.000,0.000,0.000]]),
                (18, [[-1.300,0.000,0.000], [-0.100,0.000,0.000], [0.000,-0.200,0.000], [0.000,0.000,0.000], [0.000,0.000,0.200], [0.000,0.000,0.000]])
                ]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-17',
                    'name': get_stomach_term('stomach')[0],
                    'ontId': get_stomach_term('stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_stomach_term('esophagus part of stomach')[0],
                    'ontId': get_stomach_term('esophagus part of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-8',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '9-11',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '12-13',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '14-15',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '16-17',
                    'name': get_stomach_term('duodenum part of stomach')[0],
                    'ontId': get_stomach_term('duodenum part of stomach')[1]
                }]
        })

def getDefaultOstiumSettings():
    """
    Generate list of default options for ostium.
    """
    options = { 'Number of elements around ostium': 8,
                'Number of elements along': 3,
                'Number of elements through wall': 4,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium wall thickness': 0.0525,
                'Ostium wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Use linear through ostium wall': True,
                'Vessel wall thickness': 0.0315,
                'Vessel wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1}

    return options

def updateOstiumOptions(options, ostiumOptions):
    """
    Update ostium sub-scaffold options which depend on parent options.
    """
    ostiumOptions['Number of elements around ostium'] = options['Number of elements around esophagus']
    ostiumOptions['Ostium wall thickness'] = options['Wall thickness']
    ostiumOptions['Vessel wall thickness'] = options['Esophagus wall thickness']
    elementsCountThroughWall = options['Number of elements through wall']
    ostiumOptions['Number of elements through wall'] = elementsCountThroughWall
    ostiumOptions['Use linear through ostium wall'] = options['Use linear through wall']
    ostiumOptions['Use linear through vessel wall'] = options['Use linear through wall']
    if elementsCountThroughWall == 1:
        ostiumOptions['Ostium wall relative thicknesses'] = [1.0]
        ostiumOptions['Vessel wall relative thicknesses'] = [1.0]
    else:
        mucosaRelThickness = options['Mucosa relative thickness']
        submucosaRelThickness = options['Submucosa relative thickness']
        circularRelThickness = options['Circular muscle layer relative thickness']
        longRelThickness = options['Longitudinal muscle layer relative thickness']
        relThicknesses = [mucosaRelThickness, submucosaRelThickness, circularRelThickness, longRelThickness]
        ostiumOptions['Ostium wall relative thicknesses'] = relThicknesses
        ostiumOptions['Vessel wall relative thicknesses'] = relThicknesses

    return ostiumOptions


class MeshType_3d_stomach1(Scaffold_base):
    """
    Generates a 3-D stomach mesh with variable numbers of elements around the esophagus and duodenum,
    along the central line, and through wall. The stomach is created using a network layout as the longitudinal axis
    of the stomach. D2 of the network layout points to the greater curvature of the stomach and magnitude of D2 and D3
    are the radii of the stomach in the respective direction.
    """

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
        options = {
            'Network layout': getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName),
            'Number of elements around esophagus': 8,
            'Number of elements around duodenum': 16,
            'Number of elements along': 14,
            'Number of elements through wall': 4,
            'Wall thickness': 0.0525,
            'Esophagus wall thickness': 0.0315,
            'Mucosa relative thickness': 0.55,
            'Submucosa relative thickness': 0.15,
            'Circular muscle layer relative thickness': 0.25,
            'Longitudinal muscle layer relative thickness': 0.05,
            'Limiting ridge': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements surface': 4,
            'Refine number of elements cardia surface': 2,
            'Refine number of elements through wall': 1
        }
        if 'Human 2' in parameterSetName:
            options['Number of elements around duodenum'] = 12
            options['Number of elements through wall'] = 1
            options['Wall thickness'] = 3.0
            options['Esophagus wall thickness'] = 3.0
        elif 'Mouse 1' in parameterSetName:
            options['Wall thickness'] = 0.05145
            options['Esophagus wall thickness'] = 0.01029
            options['Mucosa relative thickness'] = 0.75
            options['Submucosa relative thickness'] = 0.05
            options['Circular muscle layer relative thickness'] = 0.15
            options['Longitudinal muscle layer relative thickness'] = 0.05
            options['Limiting ridge'] = True
        elif 'Pig 1' in parameterSetName:
            options['Wall thickness'] = 0.059
            options['Esophagus wall thickness'] = 0.0354
            options['Mucosa relative thickness'] = 0.47
            options['Submucosa relative thickness'] = 0.1
            options['Circular muscle layer relative thickness'] = 0.33
            options['Longitudinal muscle layer relative thickness'] = 0.1
            options['Limiting ridge'] = False
        elif 'Rat 1' in parameterSetName:
            options['Wall thickness'] = 0.0215
            options['Esophagus wall thickness'] = 0.0129
            options['Mucosa relative thickness'] = 0.65
            options['Submucosa relative thickness'] = 0.12
            options['Circular muscle layer relative thickness'] = 0.18
            options['Longitudinal muscle layer relative thickness'] = 0.05
            options['Limiting ridge'] = True
        elif 'Material' in parameterSetName:
            options['Wall thickness'] = 0.05
            options['Esophagus wall thickness'] = 0.03
            options['Mucosa relative thickness'] = 0.25
            options['Submucosa relative thickness'] = 0.25
            options['Circular muscle layer relative thickness'] = 0.25
            options['Longitudinal muscle layer relative thickness'] = 0.25
            options['Limiting ridge'] = False

        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Network layout',
            'Number of elements around esophagus',
            'Number of elements around duodenum',
            'Number of elements along',
            'Number of elements through wall',
            'Wall thickness',
            'Esophagus wall thickness',
            'Mucosa relative thickness',
            'Submucosa relative thickness',
            'Circular muscle layer relative thickness',
            'Longitudinal muscle layer relative thickness',
            'Limiting ridge',
            'Use linear through wall',
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements cardia surface',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Network layout':
            return [MeshType_1d_network_layout1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Network layout':
            return cls.getParameterSetNames()
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
                parameterSetName = "Default"
            return getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName)
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Network layout'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Network layout'):
            options['Network layout'] = cls.getOptionScaffoldPackage('Network layout', MeshType_1d_network_layout1)
        if options['Number of elements around esophagus'] < 8:
            options['Number of elements around esophagus'] = 8
        if options['Number of elements around duodenum'] < 12:
            options['Number of elements around duodenum'] = 12
        for key in ['Number of elements around esophagus',
                    'Number of elements around duodenum']:
            if options[key] % 4 > 0:
                options[key] = options[key] // 4 * 4
        if options['Number of elements along'] < 12:
            options['Number of elements along'] = 12
        if options['Number of elements through wall'] != (1 or 4):
            options['Number of elements through wall'] = 4
        for key in [
            'Refine number of elements surface',
            'Refine number of elements cardia surface',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """

        geometricNetworkLayout = options['Network layout']
        materialNetworkLayout = getDefaultNetworkLayoutScaffoldPackage(cls, 'Material')
        limitingRidge = options['Limiting ridge']
        elementsCountThroughWall = options['Number of elements through wall']

        allAnnotationGroups = []

        stomachTermsAlong = ['stomach', 'fundus of stomach', 'body of stomach',
                             'pyloric antrum', 'pyloric canal', 'duodenum part of stomach', 'esophagus part of stomach']

        # Geometric coordinates
        fm = region.getFieldmodule()
        coordinates = find_or_create_field_coordinates(fm)

        geometricNetworkLayout = StomachNetworkLayout(region, geometricNetworkLayout, stomachTermsAlong)

        allAnnotationGroups, nextNodeIdentifier, nextElementIdentifier, elementsAlongGroups = \
            createStomachMesh3d(region, fm, coordinates, stomachTermsAlong,
                                allAnnotationGroups, networkLayout=geometricNetworkLayout,
                                options=options, nodeIdentifier=1, elementIdentifier=1)[0:4]

        # Material coordinates
        stomach_coordinates = find_or_create_field_coordinates(fm, name="stomach coordinates")
        allAnnotationGroupsMaterial = []
        tmp_region = region.createRegion()
        tmp_fm = tmp_region.getFieldmodule()
        with ChangeManager(tmp_fm):
            tmp_stomach_coordinates = find_or_create_field_coordinates(tmp_fm, name="stomach coordinates")

            materialNetworkLayout = StomachNetworkLayout(tmp_region, materialNetworkLayout, stomachTermsAlong)

            allAnnotationGroupsMaterial, nextNodeIdentifier, nextElementIdentifier = \
                createStomachMesh3d(tmp_region, tmp_fm, tmp_stomach_coordinates, stomachTermsAlong,
                                    allAnnotationGroupsMaterial,
                                    networkLayout=materialNetworkLayout, options=options, nodeIdentifier=1,
                                    elementIdentifier=1, elementsAlongSections=elementsAlongGroups,
                                    materialCoordinates=True)[:3]

            # Write two coordinates
            sir = tmp_region.createStreaminformationRegion()
            srm = sir.createStreamresourceMemory()
            tmp_region.write(sir)
            result, buffer = srm.getBuffer()

            sir = region.createStreaminformationRegion()
            srm = sir.createStreamresourceMemoryBuffer(buffer)
            region.read(sir)

            del srm
            del sir
            del tmp_stomach_coordinates
        del tmp_fm
        del tmp_region

        # Create markers
        markerTermNameStomachCoordinatesMap = {
            'body-antrum junction along the greater curvature on luminal surface': [-0.6, -0.45, 6.34622e-18],
            'body-antrum junction along the greater curvature on serosa': [-0.6, -0.5, 0.0],
            'distal point of lower esophageal sphincter serosa on the greater curvature of stomach': [0.08, 0.8, 4.48345e-16],
            'distal point of lower esophageal sphincter serosa on the lesser curvature of stomach': [-0.08, 0.8, 4.67938e-16],
            'esophagogastric junction along the greater curvature on luminal surface': [0.14885, 0.451205, 3.88484e-14],
            'esophagogastric junction along the greater curvature on serosa': [0.149987, 0.501192, 3.72966e-16],
            'esophagogastric junction along the lesser curvature on luminal surface': [-0.15, 0.45, 3.33066e-16],
            'esophagogastric junction along the lesser curvature on serosa': [-0.15, 0.5, 2.28983e-16],
            'gastroduodenal junction along the greater curvature on luminal surface': [-1.1, -0.15, 7.93284e-18],
            'gastroduodenal junction along the greater curvature on serosa': [-1.1, -0.2, 0],
            'gastroduodenal junction along the lesser curvature on luminal surface': [-1.1, 0.15, -4.73333e-17],
            'gastroduodenal junction along the lesser curvature on serosa': [-1.1, 0.2, -2.77556e-16],
            'limiting ridge at the greater curvature on the luminal surface' if limitingRidge else
            'fundus-body junction along the greater curvature on luminal surface': [-2.60734e-23, -0.450001, -0.00024468],
            'fundus-body junction along the greater curvature on serosa': [2.77556e-17, -0.5, 5.74685e-16]
        }
        if elementsCountThroughWall == 4:
            markerTermNameStomachCoordinatesCMLMMap = {
                'body-antrum junction along the greater curvature on circular-longitudinal muscle interface': [-0.6, -0.4875, -8.32667e-17],
                'esophagogastric junction along the greater curvature on circular-longitudinal muscle interface': [0.149703, 0.488695, 9.99176e-15],
                'esophagogastric junction along the lesser curvature on circular-longitudinal muscle interface': [-0.15, 0.4875, 2.76195e-16],
                'gastroduodenal junction along the greater curvature on circular-longitudinal muscle interface': [-1.1, -0.1875, 1.66533e-16],
                'gastroduodenal junction along the lesser curvature on circular-longitudinal muscle interface': [-1.1, 0.1875, -2.24625e-16],
                'limiting ridge at the greater curvature on the circular-longitudinal muscle interface' if limitingRidge
                else 'fundus-body junction along the greater curvature on circular-longitudinal muscle interface': [3.75751e-17, -0.4875, -6.117e-05]
            }
            markerTermNameStomachCoordinatesMap.update(markerTermNameStomachCoordinatesCMLMMap)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodeIdentifier = max(1, get_maximum_node_identifier(nodes) + 1)

        for termName, stomachCoordinatesValues in markerTermNameStomachCoordinatesMap.items():
            annotationGroup = findOrCreateAnnotationGroupForTerm(
                allAnnotationGroups, region, get_stomach_term(termName), isMarker=True)
            annotationGroup.createMarkerNode(nodeIdentifier, stomach_coordinates, stomachCoordinatesValues)
            nodeIdentifier += 1

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
        refineElementsCountAlongCardia = options['Refine number of elements cardia surface']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        sourceFm = meshrefinement._sourceFm
        annotationGroups = meshrefinement._sourceAnnotationGroups
        cardiaGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("cardia of stomach"))
        cardiaMeshGroup = cardiaGroup.getMeshGroup(meshrefinement._sourceMesh)

        lastElementIdentifier = get_maximum_element_identifier(meshrefinement._sourceMesh)

        cache = sourceFm.createFieldcache()
        element = meshrefinement._sourceElementiterator.next()
        while element.isValid():
            elementIdentifier = element.getIdentifier()
            refineElements1 = refineElementsCountAround
            refineElements2 = refineElementsCountAlong
            refineElements3 = refineElementsCountThroughWall
            cache.setElement(element)
            if cardiaMeshGroup.containsElement(element):
                refineElements2 = refineElementsCountAlongCardia

            meshrefinement.refineElementCubeStandard3d(element, refineElements1, refineElements2, refineElements3)
            if elementIdentifier == lastElementIdentifier:
                return  # finish on last so can continue elsewhere
            element = meshrefinement._sourceElementiterator.next()


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

        limitingRidge = options['Limiting ridge']
        elementsCountThroughWall = options['Number of elements through wall']

        stomachGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("stomach"))
        dorsalStomachGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("dorsal stomach"))
        ventralStomachGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("ventral stomach"))
        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("body of stomach"))
        cardiaGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("cardia of stomach"))
        duodenumGroup = getAnnotationGroupForTerm(annotationGroups, get_smallintestine_term("duodenum"))
        fundusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("fundus of stomach"))
        antrumGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pyloric antrum"))
        pylorusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pyloric canal"))
        esoGroup = getAnnotationGroupForTerm(annotationGroups, get_esophagus_term("esophagus"))
        nearLCGroup = getAnnotationGroupForTerm(annotationGroups,
                                                ("elements adjacent to lesser curvature", "None"))

        # Create new groups
        stomachLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_stomach_term("luminal surface of stomach"))
        stomachSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_stomach_term("serosa of stomach"))
        bodyLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_stomach_term("luminal surface of body of stomach"))
        bodySerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                             get_stomach_term("serosa of body of stomach"))
        cardiaLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_stomach_term(
                                                                    "luminal surface of cardia of stomach"))
        cardiaSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_stomach_term("serosa of cardia of stomach"))
        duodenumLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_smallintestine_term("luminal surface of duodenum"))
        duodenumSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_smallintestine_term("serosa of duodenum"))
        esoLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                             get_esophagus_term("luminal surface of esophagus"))
        esoSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_esophagus_term("serosa of esophagus"))
        fundusLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_stomach_term(
                                                                    "luminal surface of fundus of stomach"))
        fundusSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_stomach_term("serosa of fundus of stomach"))
        antrumLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_stomach_term("luminal surface of pyloric antrum"))
        antrumSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                               get_stomach_term("serosa of pyloric antrum"))
        pylorusLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_stomach_term("luminal surface of pyloric canal"))
        pylorusSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_stomach_term("serosa of pyloric canal"))
        gastroduodenalJunctionGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                         get_stomach_term("gastroduodenal junction"))

        fm = region.getFieldmodule()
        mesh2d = fm.findMeshByDimension(2)
        mesh1d = fm.findMeshByDimension(1)

        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_outer = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_inner = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))

        is_gastroduod = fm.createFieldAnd(duodenumGroup.getGroup(), pylorusGroup.getGroup())
        gastroduodenalJunctionGroup.getMeshGroup(mesh2d).addElementsConditional(is_gastroduod)

        is_dorsal = dorsalStomachGroup.getGroup()
        is_ventral = ventralStomachGroup.getGroup()
        is_curvatures = fm.createFieldAnd(is_dorsal, is_ventral)

        is_nearLC = nearLCGroup.getGroup()
        is_lesserCurvature = fm.createFieldAnd(is_curvatures, is_nearLC)
        lesserCurvatureGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_stomach_term("lesser curvature of stomach"))
        lesserCurvatureGroup.getMeshGroup(mesh2d).addElementsConditional(is_lesserCurvature)

        is_nearGC = fm.createFieldNot(is_nearLC)
        is_greaterCurvature = fm.createFieldAnd(is_curvatures, is_nearGC)
        greaterCurvatureGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_stomach_term("greater curvature of stomach"))
        greaterCurvatureGroup.getMeshGroup(mesh2d).addElementsConditional(is_greaterCurvature)

        if elementsCountThroughWall == 4:
            CMLMInterfaceGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                "circular-longitudinal muscle interface of stomach"))
            circularMuscleGroup = getAnnotationGroupForTerm(annotationGroups,
                                                            get_stomach_term("circular muscle layer of stomach"))
            longitudinalMuscleGroup = \
                getAnnotationGroupForTerm(annotationGroups, get_stomach_term("longitudinal muscle layer of stomach"))
            is_CM = circularMuscleGroup.getGroup()
            is_LM = longitudinalMuscleGroup.getGroup()
            is_CMLMInterface = fm.createFieldAnd(is_CM, is_LM)
            CMLMInterfaceGroup.getMeshGroup(mesh2d).addElementsConditional(is_CMLMInterface)
            is_curvatures_CMLM = fm.createFieldAnd(is_curvatures, is_CMLMInterface)
            is_greaterCurvature_CMLM = fm.createFieldAnd(is_greaterCurvature, is_CMLMInterface)
            is_lesserCurvature_CMLM = fm.createFieldAnd(is_lesserCurvature, is_CMLMInterface)

            dorsalStomach_CMLMGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                "circular-longitudinal muscle interface of dorsal stomach"))
            ventralStomach_CMLMGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                "circular-longitudinal muscle interface of ventral stomach"))
            is_dorsal_CMLM = fm.createFieldAnd(is_dorsal, is_CMLMInterface)
            dorsalStomach_CMLMGroup.getMeshGroup(mesh2d).addElementsConditional(is_dorsal_CMLM)
            is_ventral_CMLM = fm.createFieldAnd(is_ventral, is_CMLMInterface)
            ventralStomach_CMLMGroup.getMeshGroup(mesh2d).addElementsConditional(is_ventral_CMLM)

            gastroduod_CMLMGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                "circular-longitudinal muscle interface of gastroduodenal junction"))
            is_gastroduod_CMLM = fm.createFieldAnd(is_gastroduod, is_CMLMInterface)
            gastroduod_CMLMGroup.getMeshGroup(mesh1d).addElementsConditional(is_gastroduod_CMLM)

            bodyCurvaturesCMLMGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                    "circular-longitudinal muscle interface of body of stomach along the gastric-omentum attachment"))
            duodenumCurvaturesCMLMGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_smallintestine_term("circular-longitudinal muscle interface of "
                                                                           "first segment of the duodenum along the "
                                                                           "gastric-omentum attachment"))
            esoCurvaturesCMLMGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_esophagus_term(
                    "circular-longitudinal muscle interface of esophagus along the cut margin"))
            fundusCurvaturesCMLMGroup =\
                findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                    "circular-longitudinal muscle interface of fundus of stomach along the greater curvature"))
            antrumGreaterCurvatureCMLMGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                    "circular-longitudinal muscle interface of pyloric antrum along the greater curvature"))
            antrumLesserCurvatureCMLMGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                    "circular-longitudinal muscle interface of pyloric antrum along the lesser curvature"))
            pylorusGreaterCurvatureCMLMGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                    "circular-longitudinal muscle interface of pyloric canal along the greater curvature"))
            pylorusLesserCurvatureCMLMGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                    "circular-longitudinal muscle interface of pyloric canal along the lesser curvature"))

            sectionCurvaturesCMLMGroups = [None, bodyCurvaturesCMLMGroup, None, duodenumCurvaturesCMLMGroup,
                                           esoCurvaturesCMLMGroup, fundusCurvaturesCMLMGroup,
                                           antrumGreaterCurvatureCMLMGroup, pylorusGreaterCurvatureCMLMGroup,
                                           antrumLesserCurvatureCMLMGroup, pylorusLesserCurvatureCMLMGroup]

        sectionGroups = [stomachGroup, bodyGroup, cardiaGroup, duodenumGroup, esoGroup, fundusGroup, antrumGroup,
                         pylorusGroup]
        sectionSerosaGroups = [stomachSerosaGroup, bodySerosaGroup, cardiaSerosaGroup, duodenumSerosaGroup,
                               esoSerosaGroup, fundusSerosaGroup, antrumSerosaGroup,
                               pylorusSerosaGroup]
        sectionLuminalGroups = [stomachLuminalGroup, bodyLuminalGroup, cardiaLuminalGroup, duodenumLuminalGroup,
                                esoLuminalGroup, fundusLuminalGroup, antrumLuminalGroup,
                                pylorusLuminalGroup]

        for i in range(len(sectionGroups)):
            is_section = sectionGroups[i].getGroup()
            is_sectionSerosa = fm.createFieldAnd(is_section, is_exterior_face_outer)
            sectionSerosaGroups[i].getMeshGroup(mesh2d).addElementsConditional(is_sectionSerosa)
            is_sectionLuminal = fm.createFieldAnd(is_section, is_exterior_face_inner)
            sectionLuminalGroups[i].getMeshGroup(mesh2d).addElementsConditional(is_sectionLuminal)

            if elementsCountThroughWall == 4:
                if sectionGroups[i] is antrumGroup or sectionGroups[i] is pylorusGroup:
                    is_sectionGreaterCurvatureCMLM = fm.createFieldAnd(is_section, is_greaterCurvature_CMLM)
                    is_sectionLesserCurvatureCMLM = fm.createFieldAnd(is_section, is_lesserCurvature_CMLM)
                    if sectionCurvaturesCMLMGroups[i]:
                        sectionCurvaturesCMLMGroups[i].getMeshGroup(mesh1d). \
                            addElementsConditional(is_sectionGreaterCurvatureCMLM)
                        sectionCurvaturesCMLMGroups[i+2].getMeshGroup(mesh1d). \
                            addElementsConditional(is_sectionLesserCurvatureCMLM)
                else:
                    is_sectionCurvaturesCMLM = fm.createFieldAnd(is_section, is_curvatures_CMLM)
                    if sectionCurvaturesCMLMGroups[i]:
                        sectionCurvaturesCMLMGroups[i].getMeshGroup(mesh1d). \
                            addElementsConditional(is_sectionCurvaturesCMLM)

        if limitingRidge:
            limitingRidgeGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_stomach_term("forestomach-glandular stomach junction"))
            innerLimitingRidgeGroup = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_stomach_term("limiting ridge on luminal surface"))
            outerLimitingRidgeGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                         get_stomach_term("limiting ridge on serosa"))

            is_antrum = antrumGroup.getGroup()
            is_body = bodyGroup.getGroup()
            is_cardia = cardiaGroup.getGroup()
            is_fundus = fundusGroup.getGroup()
            is_limitingRidgeBody = fm.createFieldAnd(is_fundus, is_body)
            is_limitingRidgeCardia = fm.createFieldAnd(is_body, is_cardia)
            is_limitingRidgeAntrum = fm.createFieldAnd(is_antrum, is_cardia)
            is_limitingRidgeBodyCardia = fm.createFieldOr(is_limitingRidgeBody, is_limitingRidgeCardia)
            is_limitingRidgeAntrumCardia = fm.createFieldOr(is_limitingRidgeAntrum, is_limitingRidgeCardia)
            is_limitingRidge = fm.createFieldOr(is_limitingRidgeBodyCardia, is_limitingRidgeAntrumCardia)

            if elementsCountThroughWall == 4:
                mucosaGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("mucosa of stomach"))
                is_mucosa = mucosaGroup.getGroup()
                is_bodyMucosa = fm.createFieldAnd(is_body, is_mucosa)
                is_antrumMucosa = fm.createFieldAnd(is_antrum, is_mucosa)
                is_bodyAntrumMucosa = fm.createFieldOr(is_bodyMucosa, is_antrumMucosa)

                is_xi1Interior = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0),
                                                   fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1))
                is_xi1All = fm.createFieldOr(fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_0),
                                             fm.createFieldIsOnFace(Element.FACE_TYPE_XI1_1))
                is_xi1BodyAntrumMucosaAll = fm.createFieldAnd(is_bodyAntrumMucosa, is_xi1All)
                is_limitingRidgeAroundCardia = fm.createFieldAnd(is_xi1BodyAntrumMucosaAll,
                                                                 fm.createFieldNot(is_xi1Interior))

                is_xi2Interior = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0),
                                                   fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1))
                is_xi2ZeroBodyMucosa = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0), is_bodyMucosa)
                is_limitingRidgeBodyBoundary = fm.createFieldAnd(is_xi2ZeroBodyMucosa,
                                                                 fm.createFieldNot(is_xi2Interior))

                is_limitingRidgeMucosa = fm.createFieldOr(is_limitingRidgeAroundCardia, is_limitingRidgeBodyBoundary)
                is_limitingRidge = fm.createFieldOr(is_limitingRidge, is_limitingRidgeMucosa)

            limitingRidgeGroup.getMeshGroup(mesh2d).addElementsConditional(is_limitingRidge)

            is_xi3Interior = fm.createFieldAnd(fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0),
                                               fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
            is_xi3ZeroLimitingRidge = fm.createFieldAnd(is_limitingRidge,
                                                        fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
            is_xi3OneLimitingRidge = fm.createFieldAnd(is_limitingRidge,
                                                       fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))

            is_limitingRidgeInner = fm.createFieldAnd(is_xi3ZeroLimitingRidge, fm.createFieldNot(is_xi3Interior))
            innerLimitingRidgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_limitingRidgeInner)

            is_limitingRidgeOuter = fm.createFieldAnd(is_xi3OneLimitingRidge, fm.createFieldNot(is_xi3Interior))
            outerLimitingRidgeGroup.getMeshGroup(mesh1d).addElementsConditional(is_limitingRidgeOuter)

            if elementsCountThroughWall == 4:
                limitingRidge_CMLMGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
                    "limiting ridge on circular-longitudinal muscle interface"))
                is_limitingRidgeCMLM = fm.createFieldAnd(is_CMLMInterface, is_limitingRidge)
                limitingRidge_CMLMGroup.getMeshGroup(mesh1d).addElementsConditional(is_limitingRidgeCMLM)

        annotationGroups.remove(nearLCGroup)

class StomachNetworkLayout:
    """
    Generates sampled network layout for stomach scaffold.
    """
    def __init__(self, region, networkLayout, stomachTermsAlong=[None], esoSegmentIdx=0, stomachSegmentIdx=[1,2]):
        """
        :param region: Zinc region needed to create path region to define path in.
        :param networkLayout: Network layout subscaffold from meshtype_1d_network_layout1
        :param stomachTermsAlong: Annotation terms along length of stomach
        :param esoSegmentIdx: Segment index of esophagus branch.
        :param stomachSegmentIdx: Segment index of the body of stomach.
        """
        # Extract length of each group along stomach from network layout
        arcLengthOfGroupsAlong = []
        cxGroups = []
        cd1Groups = []
        cd2Groups = []
        cd3Groups = []
        cd12Groups = []
        cd13Groups = []

        tmpRegion = region.createRegion()
        networkLayout.generate(tmpRegion)
        pathNetworkMesh = networkLayout.getConstructionObject()
        tmpFieldmodule = tmpRegion.getFieldmodule()
        tmpNodes = tmpFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        tmpCoordinates = tmpFieldmodule.findFieldByName('coordinates')
        networkSegments = pathNetworkMesh.getNetworkSegments()

        cxGroup = []
        cd1Group = []
        cd2Group = []
        cd3Group = []
        cd12Group = []
        cd13Group = []
        stomachNodes = []
        esoNodes = []
        esoVersions = []
        lowerStomachNodes = []
        lowerStomachVersions = []

        for termName in stomachTermsAlong:
            tmpGroup = tmpFieldmodule.findFieldByName(termName).castGroup() if termName else None
            tmpNodeset = tmpGroup.getNodesetGroup(tmpNodes) if tmpGroup else tmpNodes

            if termName == "stomach":
                nodeiterator = tmpNodeset.createNodeiterator()
                node = nodeiterator.next()
                while node.isValid():
                    stomachNodes.append(node.getIdentifier())
                    node = nodeiterator.next()

                for i in range(len(networkSegments[stomachSegmentIdx[1]].getNodeIdentifiers())):
                    if networkSegments[stomachSegmentIdx[1]].getNodeIdentifiers()[i] in stomachNodes:
                        lowerStomachNodes.append(networkSegments[stomachSegmentIdx[1]].getNodeIdentifiers()[i])
                        lowerStomachVersions.append(networkSegments[stomachSegmentIdx[1]].getNodeVersions()[i])

                for i in range(2):
                    cx, cd1, cd2, cd3, cd12, cd13 = get_nodeset_path_ordered_field_parameters(
                        tmpNodeset, tmpCoordinates,
                        [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                         Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                         Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                        (lowerStomachNodes if i else networkSegments[stomachSegmentIdx[0]].getNodeIdentifiers()),
                        (lowerStomachVersions if i else networkSegments[stomachSegmentIdx[0]].getNodeVersions()))

                    cxGroup += cx[(1 if i else 0):]
                    cd1Group += cd1[(1 if i else 0):]
                    cd2Group += cd2[(1 if i else 0):]
                    cd3Group += cd3[(1 if i else 0):]
                    cd12Group += cd12[(1 if i else 0):]
                    cd13Group += cd13[(1 if i else 0):]

                    if i == 0:
                        xbranchpt = cx[-1]
                        d2branchpt = cd2[-1]
                        d3branchpt = cd3[-1]
                        arcLengthToBranchPt = 0.0
                        for n in range(len(cx) - 1):
                            arcLengthToBranchPt += interp.getCubicHermiteArcLength(cx[n], cd1[n], cx[n + 1], cd1[n + 1])

            elif termName == "esophagus part of stomach":
                for i in range(len(networkSegments[esoSegmentIdx].getNodeIdentifiers())):
                    if networkSegments[esoSegmentIdx].getNodeIdentifiers()[i] in stomachNodes:
                        esoNodes.append(networkSegments[esoSegmentIdx].getNodeIdentifiers()[i])
                        esoVersions.append(networkSegments[esoSegmentIdx].getNodeVersions()[i])

                cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                    get_nodeset_path_ordered_field_parameters(tmpNodeset, tmpCoordinates,
                                                              [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                               Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                               Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                                              esoNodes, esoVersions)

            elif termName == "fundus of stomach":
                cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                    get_nodeset_path_ordered_field_parameters(tmpNodeset, tmpCoordinates,
                                                              [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                               Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                               Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                                              networkSegments[1].getNodeIdentifiers(),
                                                              networkSegments[1].getNodeVersions())
            else:
                cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                    get_nodeset_path_field_parameters(tmpNodeset, tmpCoordinates,
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
        self.xBranchPt = xbranchpt
        self.d2BranchPt = d2branchpt
        self.d3BranchPt = d3branchpt
        self.arcLengthToBranchPt = arcLengthToBranchPt


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


def createStomachMesh3d(region, fm, coordinates, stomachTermsAlong, allAnnotationGroups, networkLayout, options,
                        nodeIdentifier, elementIdentifier, elementsAlongSections = [],
                        materialCoordinates=False, nodeIdProximalEso=[], xProximalEso=[], d1ProximalEso=[],
                        d2ProximalEso=[], d3ProximalEso=[]):
    """
    Generates a stomach scaffold in the region using a network layout and parameter options.
    :param region: Region to create elements in.
    :param fm: Zinc fieldModule to create elements in.
    :param coordinates: Coordinate field to define nodes and elements.
    :param stomachTermsAlong: Annotation terms along length of stomach.
    :param allAnnotationGroups: List of annotation groups.
    :param networkLayout: Network layout through the axis of the stomach scaffold.
    :param options: Parameter options for stomach scaffold.
    :param nodeIdentifier: First node identifier.
    :param elementIdentifier: First element identifier.
    :param elementsAlongSections: Number of elements along each section.
    :param materialCoordinates: Create material coordinates if True.
    :param nodeIdProximalEso, xProximalEso, d1ProximalEso, d2ProximalEso, d3ProximalEso: Identifier, coordinates and
    derivatives of nodes to use at the start of esophagus section joining to the stomach.
    :return allAnnotationGroups, nextNodeIdentifier, nextElementIdentifier, elementsAlongSections, nodeIdxDistal,
    xDistal, d1Distal, d2Distal, d3Distal, arclengthDuodenumCP, xPrev, d2Prev
    """
    elementsCountAroundEso = options['Number of elements around esophagus']
    elementsCountAroundDuod = options['Number of elements around duodenum']
    elementsCountAlong = options['Number of elements along']
    elementsCountThroughWall = options['Number of elements through wall']
    mucosaRelThickness = options['Mucosa relative thickness']
    submucosaRelThickness = options['Submucosa relative thickness']
    circularRelThickness = options['Circular muscle layer relative thickness']
    longitudinalRelThickness = options['Longitudinal muscle layer relative thickness']
    useCrossDerivatives = False
    useCubicHermiteThroughWall = not (options['Use linear through wall'])

    ostiumOptions = getDefaultOstiumSettings()
    GEJSettings = updateOstiumOptions(options, ostiumOptions)
    elementsAlongEsophagus = GEJSettings['Number of elements along']
    elementsThroughEsophagusWall = GEJSettings['Number of elements through wall']
    ostiumRadius = magnitude(networkLayout.cd2Groups[-1][1])
    limitingRidge = options['Limiting ridge']
    wallThickness = options['Wall thickness']
    GEJSettings['Use linear through ostium wall'] = options['Use linear through wall']
    GEJSettings['Use linear through vessel wall'] = options['Use linear through wall']

    elementsCountAcrossCardia = 1
    cardiaDiameterFactor = 1.4  # scale to ostium diameter
    sf = (cardiaDiameterFactor - 1) * ostiumRadius

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
        ostiumWallThickness = GEJSettings['Ostium wall thickness']
        ostiumWallRelThicknesses = GEJSettings['Ostium wall relative thicknesses']
        vesselWallThickness = GEJSettings['Vessel wall thickness']
        vesselWallRelThicknesses = GEJSettings['Vessel wall relative thicknesses']

        GEJSettings['Unit scale'] = 1.0
        GEJSettings['Ostium wall thickness'] = wallThickness
        GEJSettings['Ostium wall relative thicknesses'] = relThicknesses
        GEJSettings['Vessel wall thickness'] = wallThickness * 0.6
        GEJSettings['Vessel wall relative thicknesses'] = relThicknesses

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
    arcLengthOfGroupsAlong = networkLayout.arcLengthOfGroupsAlong
    stomachNetworkLayoutLength = arcLengthOfGroupsAlong[0]

    for i in range(1, len(stomachTermsAlong) - 1):
        arcLengthRatio = arcLengthOfGroupsAlong[i] / stomachNetworkLayoutLength
        arcLengthRatioForGroupsFromFundusApex.append(arcLengthRatio)

    stomachGroup = AnnotationGroup(region, get_stomach_term("stomach"))
    smallIntestineGroup = AnnotationGroup(region, get_smallintestine_term("small intestine"))
    fundusGroup = AnnotationGroup(region, get_stomach_term("fundus of stomach"))
    bodyGroup = AnnotationGroup(region, get_stomach_term("body of stomach"))
    antrumGroup = AnnotationGroup(region, get_stomach_term("pyloric antrum"))
    pylorusGroup = AnnotationGroup(region, get_stomach_term("pyloric canal"))
    duodenumGroup = AnnotationGroup(region, get_smallintestine_term("duodenum"))

    annotationGroupAlong = [[stomachGroup, fundusGroup],
                            [stomachGroup, bodyGroup],
                            [stomachGroup, antrumGroup],
                            [stomachGroup, pylorusGroup],
                            [smallIntestineGroup, duodenumGroup]]

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

    # Break network layout into elements allocation to each group
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

    for i in (list(range(1, len(stomachTermsAlong) - 2)) + [0]):  # start from body, go back to fundus
        cxGroup = networkLayout.cxGroups[i + 1]
        cd1Group = networkLayout.cd1Groups[i + 1]
        cd2Group = networkLayout.cd2Groups[i + 1]
        cd3Group = networkLayout.cd3Groups[i + 1]
        cd12Group = networkLayout.cd12Groups[i + 1]
        cd13Group = networkLayout.cd13Groups[i + 1]

        # for n2 in range(len(cxGroup)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxGroup[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd2Group[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd1Group[n2])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, cd3Group[n2])
        #     nodeIdentifier += 1

        if materialCoordinates and i == len(stomachTermsAlong) - 2:
            for n in range(len(cxGroup)):
                cd12Group[n] = zero
                cd13Group[n] = zero

        # Break body into equal sized elements, all others vary smoothly from end derivative of last section
        # Except for fundus which start at zero derivative and ends at start derivative for body
        elementsOutSection = math.ceil(arcLengthRatioForGroupsFromFundusApex[i]/targetLengthTS)
        cxSection, cd1Section, pe, pxi, psf = interp.sampleCubicHermiteCurvesSmooth(
            cxGroup, cd1Group, elementsOutSection,
            derivativeMagnitudeStart=None if (i == 1) else 0.0 if (i == 0) else magnitude(cd1Sections[-1][-1]),
            derivativeMagnitudeEnd=None if (i != 0) else magnitude(cd1Sections[0][0]))
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
    rotAxisApex = normalize(cross(d3Apex, d2Apex))

    px = sampleEllipsePoints(cxPlusSections[0][0], cd2PlusSections[0][0], cd3PlusSections[0][0],
                             0.0, math.pi * 2.0, elementsCountAroundDuod)[0]
    d2ApexAround = [cross(cross(rotAxisApex, sub(tpx, cxApex)), rotAxisApex) for tpx in px]

    rotAngle = -math.pi * 0.5
    rotFrame = axis_angle_to_rotation_matrix(rotAxisApex, rotAngle)
    for n in range(len(px)):
        d1ApexAround.append([rotFrame[j][0] * d2ApexAround[n][0] + rotFrame[j][1] * d2ApexAround[n][1] +
                             rotFrame[j][2] * d2ApexAround[n][2] for j in range(3)])

    xEllipseAroundAll = [xApex]
    d1EllipseAroundAll = [d1ApexAround]
    d2EllipseAroundAll = [d2ApexAround]
    d2Curvature = []
    curvature = [1.0 for n in range(elementsCountAroundDuod)]
    d2Curvature.append(curvature)

    count = 1
    sectionIdx = [0]
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
            if n2 == len(cxSections[s]) - 1:
                sectionIdx.append(count)
            count += 1

            if s == 0 and n2 == len(cxSections[s]) - 1:
                xGEJ = px[elementsAroundHalfDuod]

    # Scale d1 and d2 at apex
    for n in range(len(xEllipseAroundAll[0])):
        d1EllipseAroundAll[0][n] = \
            set_magnitude(d1EllipseAroundAll[0][n],
                                interp.computeCubicHermiteArcLength(xEllipseAroundAll[0][n], d2EllipseAroundAll[0][n],
                                                                    xEllipseAroundAll[1][n], d2EllipseAroundAll[1][n],
                                                                    True))
        d2EllipseAroundAll[0][n] = \
            set_magnitude(d2EllipseAroundAll[0][n],
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
            d1TrackSurface.append(d1EllipseAroundAll[n2][n1] if n2 else zero)
            d2TrackSurface.append(d2EllipseAroundAll[n2][n1])

    trackSurfaceStomach = TrackSurface(elementsCountAroundDuod, len(xEllipseAroundAll) - 1,
                                       xTrackSurface, d1TrackSurface, d2TrackSurface, loop1=True)

    # Visualise track surface
    # nodeIdentifier, elementIdentifier = trackSurfaceStomach.generateMesh(region)

    # Set up gastro-esophageal junction with network layout
    cxEso = networkLayout.cxGroups[-1]
    cd1Eso = networkLayout.cd1Groups[-1]
    cd2Eso = networkLayout.cd2Groups[-1]
    cd3Eso = networkLayout.cd3Groups[-1]
    cd12Eso = networkLayout.cd12Groups[-1]

    # Find centre position
    # track along esophagus path and since cxEso[1] could be above or below the track surface, we check both side to
    # determine direction to track. At each point, find the nearest position and take the diff between nearest point
    # to the point in line, keep tracking till diff is close to zero.

    xTol = 1.0E-6
    arcStart = 0.0
    arcEnd = networkLayout.arcLengthOfGroupsAlong[-1]
    nearestPosition = trackSurfaceStomach.findNearestPosition(cxEso[0])
    xNearestStart = trackSurfaceStomach.evaluateCoordinates(nearestPosition, derivatives=False)
    distStart = magnitude([cxEso[0][c] - xNearestStart[c] for c in range(3)])
    nearestPosition = trackSurfaceStomach.findNearestPosition(cxEso[-1])
    xNearestEnd = trackSurfaceStomach.evaluateCoordinates(nearestPosition, derivatives=False)
    distEnd = magnitude([cxEso[-1][c] - xNearestEnd[c] for c in range(3)])

    for iter in range(100):
        arcDistance = (arcStart + arcEnd) * 0.5
        x, d1 = interp.getCubicHermiteCurvesPointAtArcDistance(cxEso, cd1Eso, arcDistance)[0:2]
        nearestPosition = trackSurfaceStomach.findNearestPosition(x)
        xNearest = trackSurfaceStomach.evaluateCoordinates(nearestPosition, derivatives=False)
        dist = magnitude([x[c] - xNearest[c] for c in range(3)])

        if abs(distStart - distEnd) > xTol:
            if distStart < distEnd:
                arcEnd = arcDistance
                distEnd = dist
            else:
                arcStart = arcDistance
                distStart = dist

        else:
            xCentre, d1Centre, d2Centre = trackSurfaceStomach.evaluateCoordinates(nearestPosition, derivatives=True)
            normAxis = normalize([-d for d in d1])
            eIdx = interp.getNearestPointIndex(cxEso, xCentre) - 1
            arcLenghtSum = 0.0
            for e in range(eIdx):
                arcLenghtSum += interp.getCubicHermiteArcLength(cxEso[e], cd1Eso[e],
                                                                cxEso[e + 1], cd1Eso[e + 1])
            xi = (arcDistance - arcLenghtSum) / \
                 interp.getCubicHermiteArcLength(cxEso[eIdx], cd1Eso[eIdx], cxEso[eIdx + 1], cd1Eso[eIdx + 1])
            d2Centre = interp.interpolateCubicHermite(cd2Eso[eIdx], cd12Eso[eIdx], cd2Eso[eIdx + 1],
                                                      cd12Eso[eIdx + 1], xi)
            break
    if iter > 98:
        print('Search for ileum entry centre - Max iters reached:', iter)

    GEJSettings['Number of elements around ostium'] = elementsCountAroundEso

    esophagusGroup = AnnotationGroup(region, get_esophagus_term("esophagus"))
    esophagusMeshGroup = esophagusGroup.getMeshGroup(mesh)
    abdominalEsoGroup = AnnotationGroup(region, get_esophagus_term("abdominal part of esophagus"))
    abdominalEsoMeshGroup = abdominalEsoGroup.getMeshGroup(mesh)
    esophagogastricJunctionGroup = AnnotationGroup(region, get_stomach_term("esophagogastric junction"))
    esophagogastricJunctionMeshGroup = esophagogastricJunctionGroup.getMeshGroup(mesh)
    stomachMeshGroup = stomachGroup.getMeshGroup(mesh)
    allAnnotationGroups += [esophagusGroup, esophagogastricJunctionGroup, abdominalEsoGroup]

    ostiumWallAnnotationGroups = []
    if elementsCountThroughWall == 4:
        esophagusMucosaGroup = AnnotationGroup(region, get_esophagus_term("esophagus mucosa"))
        esophagusSubmucosaGroup = AnnotationGroup(region, get_esophagus_term("submucosa of esophagus"))
        esophagusCircularGroup = AnnotationGroup(region, get_esophagus_term("esophagus smooth muscle circular layer"))
        esophagusLongitudinalGroup = AnnotationGroup(region,
                                                     get_esophagus_term("esophagus smooth muscle longitudinal layer"))

        ostiumWallAnnotationGroups = [[esophagusMucosaGroup, mucosaGroup],
                                      [esophagusSubmucosaGroup, submucosaGroup],
                                      [esophagusCircularGroup, circularMuscleGroup],
                                      [esophagusLongitudinalGroup, longitudinalMuscleGroup]]

        allAnnotationGroups += [esophagusMucosaGroup, esophagusSubmucosaGroup,
                                esophagusCircularGroup, esophagusLongitudinalGroup]

    xPath = [cxEso[0], xCentre]
    d1Path = [cd1Eso[0], [-d for d in normAxis]]
    ostiumLength = interp.computeCubicHermiteArcLength(xPath[0], d1Path[0], xPath[1], d1Path[1],
                                                       rescaleDerivatives=True)
    d1Path[1] = set_magnitude(d1Path[1], ostiumLength*0.1)
    d2Path = [cd2Eso[0], d2Centre]
    d3Path = [cd3Eso[0], set_magnitude([-d for d in d1Centre], magnitude(d2Centre))]
    d12Path = [cd2Eso[0], [0.0, 0.0, 0.0]]
    d13Path = [cd3Eso[0], [0.0, 0.0, 0.0]]

    networkLayoutEso = CustomNetworkLayout(xPath, d1Path, d2Path, d3Path, d12Path, d13Path)

    nextNodeIdentifier, nextElementIdentifier, (o1_x, o1_d1, o1_d2, o1_d3, o1_NodeId, o1_Positions) = \
        generateOstiumMesh(region, GEJSettings, trackSurfaceStomach, networkLayoutEso,
                           nodeIdentifier, elementIdentifier, nodeIdProximal=nodeIdProximalEso,
                           xProximal=xProximalEso, d1Proximal=d1ProximalEso, d2Proximal=d2ProximalEso,
                           d3Proximal=d3ProximalEso,
                           vesselMeshGroups=[[esophagusMeshGroup, abdominalEsoMeshGroup]],
                           ostiumMeshGroups=[stomachMeshGroup, esophagogastricJunctionMeshGroup],
                           wallAnnotationGroups=ostiumWallAnnotationGroups, coordinates=coordinates)

    stomachStartNode = nextNodeIdentifier
    nodeIdentifier = nextNodeIdentifier
    elementIdentifier = nextElementIdentifier

    if materialCoordinates:
        GEJSettings['Unit scale'] = unitScale
        GEJSettings['Ostium wall thickness'] = ostiumWallThickness
        GEJSettings['Ostium wall relative thicknesses'] = ostiumWallRelThicknesses
        GEJSettings['Vessel wall thickness'] = vesselWallThickness
        GEJSettings['Vessel wall relative thicknesses'] = vesselWallRelThicknesses

    # Create location of annulus
    xAnnulusOuter = [[] for x in range(elementsCountAroundEso)]
    xAnnulusOuterPosition = [[] for x in range(elementsCountAroundEso)]
    d2AnnulusNorm = []
    d2AnnulusOuter = []
    for n1 in range(elementsCountAroundEso):
        normD2 = normalize(o1_d2[-1][n1])
        d2AnnulusNorm.append(normD2)
        d2AnnulusOuter.append(set_magnitude(o1_d2[-1][n1], sf))
        x = [o1_x[-1][n1][c] + sf * normD2[c] for c in range(3)]
        nearestPosition = trackSurfaceStomach.findNearestPosition(x, startPosition=o1_Positions[n1])
        xAnnulusOuterPosition[n1] = nearestPosition
        xAnnulusOuter[n1] = trackSurfaceStomach.evaluateCoordinates(nearestPosition)

    d1AnnulusOuter = []
    for n in range(elementsCountAroundEso):
        v1 = xAnnulusOuter[n]
        v2 = xAnnulusOuter[(n + 1) % elementsCountAroundEso]
        d = [v2[c] - v1[c] for c in range(3)]
        arcLengthAround = interp.computeCubicHermiteArcLength(v1, d, v2, d, True)
        d1 = [c * arcLengthAround for c in normalize(d)]
        d1AnnulusOuter.append(d1)

    d1AnnulusOuter = interp.smoothCubicHermiteDerivativesLoop(xAnnulusOuter, d2AnnulusOuter)
    d3Annulus = []
    for n in range(elementsCountAroundEso):
        d3 = normalize(cross(normalize(d1AnnulusOuter[n]), d2AnnulusNorm[n]))
        d3Annulus.append(d3)
    d1AnnulusCurvatureOuter = findCurvatureAroundLoop(xAnnulusOuter, d1AnnulusOuter, d3Annulus)

    # for m in range(len(xAnnulusOuter)):
    #     node = nodes.createNode(nodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAnnulusOuter[m])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1AnnulusOuter[m])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2AnnulusOuter[m])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
    #     nodeIdentifier += 1

    # Calculate arclength at quarter line between lesser and greater curvature for each region
    if not elementsAlongSections:
        xQuarterEllipseAll = []
        d2QuarterEllipseAll = []
        for n2 in range(len(xEllipseAroundAll)):
            xQuarterEllipseAll.append(xEllipseAroundAll[n2][elementsAroundQuarterDuod])
            d2QuarterEllipseAll.append(d2EllipseAroundAll[n2][elementsAroundQuarterDuod])
        totalQuarterLength = interp.getCubicHermiteCurvesLength(xQuarterEllipseAll, d2QuarterEllipseAll)

        quarterLengthSections = []
        for i in range(len(nodeStartEndSections)):
            xList = []
            d2List = []
            for n in range(nodeStartEndSections[i][0], nodeStartEndSections[i][1] + 1):
                xList.append(xQuarterEllipseAll[n])
                d2List.append(d2QuarterEllipseAll[n])
            quarterLengthSections.append(interp.getCubicHermiteCurvesLength(xList, d2List))

        targetLengthPerElement = totalQuarterLength/elementsCountAlong
        minElementsInSections = [4, 3, 2, 2, 1]

        excessElements = elementsCountAlong - sum(minElementsInSections)
        diff = [quarterLengthSections[c] - targetLengthPerElement * minElementsInSections[c]
                for c in range(len(minElementsInSections))]
        for i in range(excessElements):
            maxIdx = max(range(len(diff)), key=diff.__getitem__)
            minElementsInSections[maxIdx] += 1
            diff[maxIdx] -= targetLengthPerElement
        elementsAlongSections = minElementsInSections

    totalElementsAlong = sum(elementsAlongSections)
    arclengthDuodenumCP = arcLengthOfGroupsAlong[5]/elementsAlongSections[-1]

    xSampledAlong = [[] for n1 in range(elementsCountAroundDuod)]
    d1SampledAlong = [[] for n1 in range(elementsCountAroundDuod)]
    d2SampledAlong = [[] for n1 in range(elementsCountAroundDuod)]
    d3SampledAlong = [[] for n1 in range(elementsCountAroundDuod)]

    n1IdxAtBodyStartIdxPlusMinusOne = list(range(elementsAroundHalfDuod - 1,
                                             elementsAroundHalfDuod + 2))
    annulusIdxAtBodyStartIdxMinusOne = list(range(1, -2, -1))
    annulusIdxAtBodyStartIdxPlusOne = list(range(elementsAroundHalfEso - 1,
                                                 elementsAroundHalfEso + 2))

    # Sample from quarterDuod to annulus at the quadrant points.
    # 1st Quadrant
    startGuessPosition = \
        trackSurfaceStomach.createPositionProportion(1.0 / elementsCountAroundDuod * elementsAroundQuarterDuod,
                                                     1.0 / len(xEllipseAroundAll) * sectionIdx[1])
    aPosition = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[sectionIdx[1]][elementsAroundQuarterDuod],
                                                        startGuessPosition)
    aProportion = trackSurfaceStomach.getProportion(aPosition)

    bPosition = xAnnulusOuterPosition[elementsAroundQuarterEso]
    bProportion = trackSurfaceStomach.getProportion(bPosition)

    nx, nd1, nd2, nd3, proportions = trackSurfaceStomach.createHermiteCurvePoints(
        aProportion[0], aProportion[1], bProportion[0], bProportion[1], elementsAroundQuarterDuod - 1,
        derivativeStart=d1EllipseAroundAll[sectionIdx[1]][elementsAroundQuarterDuod],
        curveMode=TrackSurface.HermiteCurveMode.UNIFORM_SIZE)

    nxR, nd1R, nd2R, nd3R = \
        trackSurfaceStomach.resampleHermiteCurvePointsSmooth(
            nx, nd1, nd2, nd3, proportions, derivativeMagnitudeStart=
            magnitude(d1EllipseAroundAll[sectionIdx[1]][elementsAroundQuarterDuod]))[0:-1]

    # Replace the values in xEllipseAroundAll at quadrants
    for n in range(len(nxR)):
        xEllipseAroundAll[sectionIdx[1]][n + elementsAroundQuarterDuod] = nxR[n]
        d1EllipseAroundAll[sectionIdx[1]][n + elementsAroundQuarterDuod] = nd1R[n]
        d2EllipseAroundAll[sectionIdx[1]][n + elementsAroundQuarterDuod] = nd2R[n]

    # 2nd quadrant
    aPosition = xAnnulusOuterPosition[-elementsAroundQuarterEso]
    aProportion = trackSurfaceStomach.getProportion(aPosition)

    startGuessPosition = \
        trackSurfaceStomach.createPositionProportion(1.0 / elementsCountAroundDuod *
                                                     (elementsAroundQuarterDuod + elementsAroundHalfDuod),
                                                     1.0 / len(xEllipseAroundAll) * sectionIdx[1])
    bPosition = \
        trackSurfaceStomach.findNearestPosition(
            xEllipseAroundAll[sectionIdx[1]][elementsAroundQuarterDuod + elementsAroundHalfDuod], startGuessPosition)
    bProportion = trackSurfaceStomach.getProportion(bPosition)

    nx, nd1, nd2, nd3, proportions = trackSurfaceStomach.createHermiteCurvePoints(
        aProportion[0], aProportion[1], bProportion[0], bProportion[1], elementsAroundQuarterDuod - 1,
        derivativeEnd=d1EllipseAroundAll[sectionIdx[1]][elementsAroundQuarterDuod + elementsAroundHalfDuod],
        curveMode=TrackSurface.HermiteCurveMode.UNIFORM_SIZE)

    nxR, nd1R, nd2R, nd3R = \
        trackSurfaceStomach.resampleHermiteCurvePointsSmooth(
            nx, nd1, nd2, nd3, proportions, derivativeMagnitudeEnd=
            magnitude(
                d1EllipseAroundAll[sectionIdx[1]][elementsAroundQuarterDuod + elementsAroundHalfDuod]))[0:-1]

    for n in range(len(nxR)):
        xEllipseAroundAll[sectionIdx[1]][elementsAroundHalfDuod + 1 + n] = nxR[n]
        d1EllipseAroundAll[sectionIdx[1]][elementsAroundHalfDuod + 1 + n] = nd1R[n]

    for i in range(len(sectionIdx) - 1):
        s = sectionIdx[i]
        sNext = sectionIdx[i + 1]
        count = 0
        for n1 in range(len(xEllipseAroundAll[s])):
            #for each pt around, we take the point on the sectionIdx as Pt A and the point on sectionIdx + 1 as Pt B,
            # do a tracksurface sampling to divide the elements into equal sized elements while keeping the start and
            # end derivatives direction at both pts
            elementsOut = elementsAlongSections[i]
            startDerivative = d2EllipseAroundAll[s][n1]
            startDerivativeMag = None
            endDerivative = d2EllipseAroundAll[sNext][n1]
            endDerivativeMag = None

            if i == 1 and n1 in n1IdxAtBodyStartIdxPlusMinusOne:
                # find endDerivative by spacing elements out evenly as though there is no ostium
                startGuessPosition = trackSurfaceStomach.createPositionProportion(1.0 / elementsCountAroundDuod * n1,
                                                                                  1.0 / len(xEllipseAroundAll) * s)
                aPosition = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[s][n1], startGuessPosition)
                aProportion = trackSurfaceStomach.getProportion(aPosition)

                startGuessPosition = trackSurfaceStomach.createPositionProportion(1.0 / elementsCountAroundDuod * n1,
                                                                                  1.0 / len(xEllipseAroundAll) * sNext)
                bPosition = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[sNext][n1], startGuessPosition)
                bProportion = trackSurfaceStomach.getProportion(bPosition)
                nx, nd1, nd2, nd3, proportions = trackSurfaceStomach.createHermiteCurvePoints(
                    aProportion[0], aProportion[1], bProportion[0], bProportion[1], elementsOut,
                    curveMode=TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
                d2Uniform = \
                    trackSurfaceStomach.resampleHermiteCurvePointsSmooth(nx, nd1, nd2, nd3, proportions)[1]
                endDerivative = d2Uniform[-1]

                # Sample from annulus to body
                aPosition = xAnnulusOuterPosition[annulusIdxAtBodyStartIdxPlusOne[count]]
                startDerivative = d2AnnulusOuter[annulusIdxAtBodyStartIdxPlusOne[count]]
                elementsOut = elementsAlongSections[i] - (elementsAroundQuarterEso - 1)
                count += 1
            else:
                startGuessPosition = trackSurfaceStomach.createPositionProportion(1.0/elementsCountAroundDuod * n1,
                                                                                  1.0/len(xEllipseAroundAll) * s)
                aPosition = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[s][n1], startGuessPosition)
            aProportion = trackSurfaceStomach.getProportion(aPosition)

            if i == 0 and n1 in n1IdxAtBodyStartIdxPlusMinusOne:
                # find startDerivative by spacing elements out evenly as though there is no ostium
                startGuessPosition = trackSurfaceStomach.createPositionProportion(1.0 / elementsCountAroundDuod * n1,
                                                                                  1.0 / len(xEllipseAroundAll) * sNext)
                bPosition = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[sNext][n1], startGuessPosition)
                bProportion = trackSurfaceStomach.getProportion(bPosition)
                nx, nd1, nd2, nd3, proportions = trackSurfaceStomach.createHermiteCurvePoints(
                    aProportion[0], aProportion[1], bProportion[0], bProportion[1], elementsOut,
                    curveMode=TrackSurface.HermiteCurveMode.UNIFORM_SIZE)
                d2Uniform = \
                    trackSurfaceStomach.resampleHermiteCurvePointsSmooth(nx, nd1, nd2, nd3, proportions)[1]
                startDerivative = d2Uniform[0]
                startDerivativeMag = magnitude(startDerivative)

                # Sample from apex to annulus
                bPosition = xAnnulusOuterPosition[annulusIdxAtBodyStartIdxMinusOne[count]]
                d = d2AnnulusOuter[annulusIdxAtBodyStartIdxMinusOne[count]]
                rotFrame = axis_angle_to_rotation_matrix(d3Annulus[annulusIdxAtBodyStartIdxMinusOne[count]],
                                                                 math.pi)
                endDerivative = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2]
                                 for j in range(3)]
                elementsOut = elementsAlongSections[i] - (elementsAroundQuarterEso - 1)
                count += 1

            else:
                startGuessPosition = trackSurfaceStomach.createPositionProportion(1.0 / elementsCountAroundDuod * n1,
                                                                                  1.0 / len(xEllipseAroundAll) * sNext)
                bPosition = trackSurfaceStomach.findNearestPosition(xEllipseAroundAll[sNext][n1], startGuessPosition)

            bProportion = trackSurfaceStomach.getProportion(bPosition)
            if n1 == 0:
                aProportion[0] = 1.0
                bProportion[0] = 1.0

            nx, nd1, nd2, nd3, proportions = trackSurfaceStomach.createHermiteCurvePoints(
                aProportion[0], aProportion[1], bProportion[0], bProportion[1], elementsOut,
                derivativeStart=startDerivative, derivativeEnd=endDerivative,
                curveMode=TrackSurface.HermiteCurveMode.UNIFORM_SIZE)

            nx, nd1, nd2, nd3 = \
                trackSurfaceStomach.resampleHermiteCurvePointsSmooth(nx, nd1, nd2, nd3, proportions,
                                                                     derivativeMagnitudeStart=startDerivativeMag,
                                                                     derivativeMagnitudeEnd=endDerivativeMag)[:-1]

            # Rotate nd2
            for m in range(len(nx)):
                rotFrame = axis_angle_to_rotation_matrix(nd3[m], math.pi)
                nd2[m] = [rotFrame[j][0] * nd2[m][0] + rotFrame[j][1] * nd2[m][1] +
                          rotFrame[j][2] * nd2[m][2] for j in range(3)]

            # Deal with annulus
            if i == 0:
                if n1 == elementsAroundHalfDuod:
                    for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
                        nx.append(zero)
                        nd1.append(zero)
                        nd2.append(zero)
                        nd3.append(zero)
                elif n1 == elementsAroundHalfDuod - 1:
                    for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
                        annulusIdx = m + 2
                        rotFrame = axis_angle_to_rotation_matrix(d3Annulus[annulusIdx], math.pi)
                        d2 = d2AnnulusOuter[annulusIdx]
                        d2 = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2]
                              for j in range(3)]
                        nx.append(xAnnulusOuter[annulusIdx])
                        nd1.append(d1AnnulusOuter[annulusIdx])
                        nd2.append(d2)
                        nd3.append(d3Annulus[annulusIdx])

                elif n1 == elementsAroundHalfDuod + 1:
                    for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
                        annulusIdx = -2 - m
                        rotFrame = axis_angle_to_rotation_matrix(d3Annulus[annulusIdx], math.pi)
                        d1 = d1AnnulusOuter[annulusIdx]
                        d1 = [rotFrame[j][0] * d1[0] + rotFrame[j][1] * d1[1] + rotFrame[j][2] * d1[2]
                              for j in range(3)]
                        nx.append(xAnnulusOuter[annulusIdx])
                        nd1.append(d1)
                        nd2.append(d2AnnulusOuter[annulusIdx])
                        nd3.append(d3Annulus[annulusIdx])

            if i == 1 and elementsAroundHalfDuod - 1 <= n1 <= elementsAroundHalfDuod + 1:
                xSampledAlong[n1] += nx
                d1SampledAlong[n1] += nd2
                d2SampledAlong[n1] += nd1
                d3SampledAlong[n1] += nd3

            else:
                xSampledAlong[n1] += nx[1:] if i else nx
                d1SampledAlong[n1] += nd2[1:] if i else nd2
                d2SampledAlong[n1] += nd1[1:] if i else nd1
                d3SampledAlong[n1] += nd3[1:] if i else nd3

    # for n1 in range(len(xSampledAlong)):
    #     for n2 in range(len(xSampledAlong[n1])):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledAlong[n1][n2])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledAlong[n1][n2])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledAlong[n1][n2])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3SampledAlong[n1][n2])
    #         nodeIdentifier += 1

    # Rearrange to around first
    xSampledAroundAlong = []
    d1SampledAroundAlong = []
    d2SampledAroundAlong = []
    d2SmoothB4ChangeAroundAlong = []
    d3SampledAroundAlong = []

    for n2 in range(totalElementsAlong + 1):
        xSampledAround = []
        d1SampledAround = []
        d2SampledAround = []
        d2SmoothB4ChangeAround = []
        d3SampledAround = []
        for n1 in range(elementsCountAroundDuod):
            xSampledAround.append(xSampledAlong[n1][n2])
            d1SampledAround.append(d1SampledAlong[n1][n2])
            d2SampledAround.append(d2SampledAlong[n1][n2])
            d2SmoothB4ChangeAround.append([])
            d3SampledAround.append(d3SampledAlong[n1][n2])
        xSampledAroundAlong.append(xSampledAround)
        d1SampledAroundAlong.append(d1SampledAround)
        d2SampledAroundAlong.append(d2SampledAround)
        d2SmoothB4ChangeAroundAlong.append(d2SmoothB4ChangeAround)
        d3SampledAroundAlong.append(d3SampledAround)

    bodyStartIdx = elementsAlongSections[0]
    annulusFundusOpenRingIdx = bodyStartIdx - (elementsAroundQuarterEso - 2)
    annulusBodyOpenRingIdx = bodyStartIdx + (elementsAroundQuarterEso - 2)

    # Smooth d1 around
    d1SmoothedAroundAlong = [d1EllipseAroundAll[0]]
    for n2 in range(1, len(xSampledAroundAlong)):
        if annulusFundusOpenRingIdx <= n2 <= annulusBodyOpenRingIdx:
            d1SmoothedLeft = \
                interp.smoothCubicHermiteDerivativesLine(xSampledAroundAlong[n2][0:elementsAroundHalfDuod],
                                                         d1SampledAroundAlong[n2][0:elementsAroundHalfDuod],
                                                         fixEndDirection=True)
            d1SmoothedRight = \
                interp.smoothCubicHermiteDerivativesLine(xSampledAroundAlong[n2][elementsAroundHalfDuod + 1:] +
                                                         [xSampledAroundAlong[n2][0]],
                                                         d1SampledAroundAlong[n2][elementsAroundHalfDuod + 1:] +
                                                         [d1SampledAroundAlong[n2][0]],
                                                         fixStartDirection=True)

            d1Smoothed = d1SmoothedLeft + [[1.0, 0.0, 0.0]] + d1SmoothedRight[:-1]

        else:
            d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xSampledAroundAlong[n2], d1SampledAroundAlong[n2])
        d1SmoothedAroundAlong.append(d1Smoothed)

    d1SampledAroundAlong = d1SmoothedAroundAlong

    # Smooth d2 along
    d2AnnulusNew = [[] for n in range(elementsCountAroundEso)]
    for n1 in range(elementsCountAroundDuod):
        nx = []
        nd2 = []
        if n1 == elementsAroundHalfDuod:
            for n2 in range(annulusFundusOpenRingIdx):
                nx.append(xSampledAroundAlong[n2][n1])
                nd2.append(d2SampledAroundAlong[n2][n1])
            d2SmoothedAlongGC = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixAllDirections=True)
            d2SmoothedAlongGCB4Change = copy.deepcopy(d2SmoothedAlongGC)
            d2SmoothedAlongGC[-1] = set_magnitude(d2AnnulusOuter[0], magnitude(d2SmoothedAlongGC[-1]))
            d2AnnulusNew[0] = d2SmoothedAlongGC[-1]

            nx = []
            nd2 = []
            for n2 in range(annulusBodyOpenRingIdx + 1, len(xSampledAroundAlong)):
                nx.append(xSampledAroundAlong[n2][n1])
                nd2.append(d2SampledAroundAlong[n2][n1])
            d2SmoothedAlongLC = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixAllDirections=True)
            d2SmoothedAlongLCB4Change = copy.deepcopy(d2SmoothedAlongLC)
            d2SmoothedAlongLC[0] = set_magnitude(d2AnnulusOuter[elementsAroundHalfEso],
                                                       magnitude(d2SmoothedAlongLC[0]))
            d2AnnulusNew[elementsAroundHalfEso] = d2SmoothedAlongLC[0]
            d2Smoothed = d2SmoothedAlongGC + \
                         [[0.0, 1.0, 0.0] for n in range(2 * (elementsAroundQuarterEso - 2) + 1)] + \
                         d2SmoothedAlongLC
            d2SmoothedB4Change = d2SmoothedAlongGCB4Change + \
                         [[0.0, 1.0, 0.0] for n in range(2 * (elementsAroundQuarterEso - 2) + 1)] + \
                         d2SmoothedAlongLCB4Change

        else:
            for n2 in range(len(xSampledAroundAlong)):
                nx.append(xSampledAroundAlong[n2][n1])
                nd2.append(d2SampledAroundAlong[n2][n1])
            d2Smoothed = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixAllDirections=True)
            d2SmoothedB4Change = copy.deepcopy(d2Smoothed)

            if n1 == elementsAroundHalfDuod - 1:
                d2Smoothed[annulusFundusOpenRingIdx - 1] = \
                    set_magnitude(d2AnnulusOuter[1], magnitude(nd2[annulusFundusOpenRingIdx - 1]))
                d2AnnulusNew[1] = d2Smoothed[annulusFundusOpenRingIdx - 1]

                for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
                    annulusIdx = m + 2
                    d2Smoothed[annulusFundusOpenRingIdx + m] = \
                        set_magnitude(d2AnnulusOuter[annulusIdx],
                                            magnitude(d1SampledAroundAlong[annulusFundusOpenRingIdx + m][n1]))
                    d2AnnulusNew[annulusIdx] = d2Smoothed[annulusFundusOpenRingIdx + m]

                d2Smoothed[annulusBodyOpenRingIdx + 1] = \
                    set_magnitude(d2AnnulusOuter[elementsAroundHalfEso - 1],
                                        magnitude(nd2[annulusBodyOpenRingIdx + 1]))
                d2AnnulusNew[elementsAroundHalfEso - 1] = d2Smoothed[annulusBodyOpenRingIdx + 1]

            if n1 == elementsAroundHalfDuod + 1:
                d2Smoothed[annulusFundusOpenRingIdx - 1] = \
                    set_magnitude(d2AnnulusOuter[-1], magnitude(nd2[annulusFundusOpenRingIdx - 1]))
                d2AnnulusNew[-1] = d2Smoothed[annulusFundusOpenRingIdx - 1]

                for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
                    annulusIdx = -(m + 2)
                    d2Smoothed[annulusFundusOpenRingIdx + m] = \
                        set_magnitude(d2AnnulusOuter[annulusIdx],
                                            magnitude(d1SampledAroundAlong[annulusFundusOpenRingIdx + m][n1]))
                    d2AnnulusNew[annulusIdx] = d2Smoothed[annulusFundusOpenRingIdx + m]

                d2Smoothed[annulusBodyOpenRingIdx + 1] = \
                    set_magnitude(d2AnnulusOuter[elementsAroundHalfEso + 1],
                                        magnitude(nd2[annulusBodyOpenRingIdx + 1]))
                d2AnnulusNew[elementsAroundHalfEso + 1] = d2Smoothed[annulusBodyOpenRingIdx + 1]

        for n2 in range(len(d2Smoothed)):
            d2SampledAroundAlong[n2][n1] = d2Smoothed[n2]
            d2SmoothB4ChangeAroundAlong[n2][n1] = d2SmoothedB4Change[n2]

    # for n2 in range(len(xSampledAroundAlong)):
    #     for n1 in range(len(xSampledAroundAlong[n2])):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3SampledAroundAlong[n2][n1])
    #         nodeIdentifier += 1

    # Replace derivatives around annulus
    for n in range(3):
        d1SampledAroundAlong[annulusFundusOpenRingIdx - 1][n1IdxAtBodyStartIdxPlusMinusOne[n]] = \
            d1AnnulusOuter[annulusIdxAtBodyStartIdxMinusOne[n]]
        d1SampledAroundAlong[annulusBodyOpenRingIdx + 1][n1IdxAtBodyStartIdxPlusMinusOne[n]] = \
            d1AnnulusOuter[annulusIdxAtBodyStartIdxPlusOne[n]]

    for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
        annulusIdx = m + 2
        d1SampledAroundAlong[annulusFundusOpenRingIdx + m][elementsAroundHalfDuod - 1] = d1AnnulusOuter[annulusIdx]
        d1SampledAroundAlong[annulusFundusOpenRingIdx + m][elementsAroundHalfDuod + 1] = d1AnnulusOuter[-annulusIdx]

    # calculate d3
    for n2 in range(len(xSampledAroundAlong)):
        for n1 in range(len(xSampledAroundAlong[n2])):
            d3SampledAroundAlong[n2][n1] = normalize(cross(
                normalize(d1SampledAroundAlong[n2][n1]), normalize(d2SampledAroundAlong[n2][n1])))

    # for n2 in range(len(xSampledAroundAlong)):
    #     for n1 in range(len(xSampledAroundAlong[n2])):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3SampledAroundAlong[n2][n1])
    #         nodeIdentifier += 1

    # Calculate curvature around
    d1CurvatureAroundAlong = [[0.0 for n in range(elementsCountAroundDuod)]]
    for n2 in range(1, len(xSampledAroundAlong)):
        if annulusFundusOpenRingIdx <= n2 <= annulusBodyOpenRingIdx:
            d1CurvatureLeft = findCurvatureAlongLine(xSampledAroundAlong[n2][0:elementsAroundHalfDuod],
                                                     d1SampledAroundAlong[n2][0:elementsAroundHalfDuod],
                                                     d3SampledAroundAlong[n2][0:elementsAroundHalfDuod])
            d1CurvatureRight = \
                findCurvatureAlongLine(
                    xSampledAroundAlong[n2][elementsAroundHalfDuod + 1:] + [xSampledAroundAlong[n2][0]],
                    d1SampledAroundAlong[n2][elementsAroundHalfDuod + 1:] + [d1SampledAroundAlong[n2][0]],
                    d3SampledAroundAlong[n2][elementsAroundHalfDuod + 1:] + [d3SampledAroundAlong[n2][0]])

            d1Curvature = d1CurvatureLeft + [0.0] + d1CurvatureRight[:-1]
        else:
            d1Curvature = findCurvatureAroundLoop(xSampledAroundAlong[n2], d1SampledAroundAlong[n2],
                                                  d3SampledAroundAlong[n2])
        d1CurvatureAroundAlong.append(d1Curvature)

    # Replace curvatures around annulus
    for n in range(3):
        d1CurvatureAroundAlong[annulusFundusOpenRingIdx - 1][n1IdxAtBodyStartIdxPlusMinusOne[n]] = \
            d1AnnulusCurvatureOuter[annulusIdxAtBodyStartIdxMinusOne[n]]
        d1CurvatureAroundAlong[annulusBodyOpenRingIdx + 1][n1IdxAtBodyStartIdxPlusMinusOne[n]] = \
            d1AnnulusCurvatureOuter[annulusIdxAtBodyStartIdxPlusOne[n]]

    for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
        annulusIdx = m + 2
        d1CurvatureAroundAlong[annulusFundusOpenRingIdx + m][elementsAroundHalfDuod - 1] = \
            d1AnnulusCurvatureOuter[annulusIdx]
        d1CurvatureAroundAlong[annulusFundusOpenRingIdx + m][elementsAroundHalfDuod + 1] = \
            d1AnnulusCurvatureOuter[-annulusIdx]

    # Calculate curvature along
    d2AnnulusCurvature = []
    for n in range(elementsCountAroundEso):
        d2AnnulusCurvature.append(interp.getCubicHermiteCurvature(o1_x[-1][n], set_magnitude(o1_d2[-1][n], sf),
                                                                  xAnnulusOuter[n], d2AnnulusNew[n], d3Annulus[n], 1.0))

    d2CurvatureAroundAlong = [[[] for n1 in range(len(xSampledAroundAlong[n2]))]
                              for n2 in range(len(xSampledAroundAlong))]

    for n1 in range(elementsCountAroundDuod):
        nx = []
        nd2 = []
        nd3 = []
        if n1 == elementsAroundHalfDuod:
            for n2 in range(annulusFundusOpenRingIdx):
                nx.append(xSampledAroundAlong[n2][n1])
                nd2.append(d2SmoothB4ChangeAroundAlong[n2][n1])
                nd3.append(d3SampledAroundAlong[n2][n1])
            d2CurvatureAlongGC = findCurvatureAlongLine(nx, nd2, nd3)

            nx = []
            nd2 = []
            nd3 = []
            for n2 in range(annulusBodyOpenRingIdx + 1, len(xSampledAroundAlong)):
                nx.append(xSampledAroundAlong[n2][n1])
                nd2.append(d2SmoothB4ChangeAroundAlong[n2][n1])
                nd3.append(d3SampledAroundAlong[n2][n1])
            d2CurvatureAlongLC = findCurvatureAlongLine(nx, nd2, nd3)
            d2CurvatureAlong = d2CurvatureAlongGC + \
                               [0.0 for n in range(2 * (elementsAroundQuarterEso - 2) + 1)] + \
                               d2CurvatureAlongLC

        else:
            for n2 in range(len(xSampledAroundAlong)):
                nx.append(xSampledAroundAlong[n2][n1])
                nd2.append(d2SmoothB4ChangeAroundAlong[n2][n1])
                nd3.append(d3SampledAroundAlong[n2][n1])
            d2CurvatureAlong = findCurvatureAlongLine(nx, nd2, nd3)

            if n1 == elementsAroundHalfDuod - 1:
                d2CurvatureAlong[annulusFundusOpenRingIdx - 1] = \
                    0.5 * (d2CurvatureAlong[annulusFundusOpenRingIdx - 1] + d2AnnulusCurvature[1])
                d2CurvatureAlong[annulusBodyOpenRingIdx + 1] = \
                    0.5 * (d2AnnulusCurvature[elementsAroundHalfEso - 1] +
                           interp.getCubicHermiteCurvature(xAnnulusOuter[elementsAroundHalfEso - 1],
                                                           d2AnnulusNew[elementsAroundHalfEso - 1],
                                                           xSampledAroundAlong[annulusBodyOpenRingIdx + 2][n1],
                                                           d2SampledAroundAlong[annulusBodyOpenRingIdx + 2][n1],
                                                           d3SampledAroundAlong[annulusBodyOpenRingIdx + 1][n1], 0.0))
                for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
                    annulusIdx = m + 2
                    d2CurvatureAlong[annulusFundusOpenRingIdx + m] = \
                        0.5 * (d2AnnulusCurvature[annulusIdx] +
                               d1CurvatureAroundAlong[annulusFundusOpenRingIdx + m][n1])

            if n1 == elementsAroundHalfDuod + 1:
                d2CurvatureAlong[annulusFundusOpenRingIdx - 1] = \
                    0.5 * (d2CurvatureAlong[annulusFundusOpenRingIdx - 1] + d2AnnulusCurvature[-1])
                d2CurvatureAlong[annulusBodyOpenRingIdx + 1] = \
                    0.5 * (d2AnnulusCurvature[elementsAroundHalfEso + 1] +
                           interp.getCubicHermiteCurvature(xAnnulusOuter[elementsAroundHalfEso + 1],
                                                           d2AnnulusNew[elementsAroundHalfEso + 1],
                                                           xSampledAroundAlong[annulusBodyOpenRingIdx + 2][n1],
                                                           d2SampledAroundAlong[annulusBodyOpenRingIdx + 2][n1],
                                                           d3SampledAroundAlong[annulusBodyOpenRingIdx + 1][n1], 0.0))
                for m in range(2 * (elementsAroundQuarterEso - 2) + 1):
                    annulusIdx = m + 2
                    d2CurvatureAlong[annulusFundusOpenRingIdx + m] = \
                        0.5 * (d2AnnulusCurvature[-annulusIdx] +
                               d1CurvatureAroundAlong[annulusFundusOpenRingIdx + m][n1])

        for n2 in range(len(d2CurvatureAlong)):
            d2CurvatureAroundAlong[n2][n1] = d2CurvatureAlong[n2]

    for i in range(annulusFundusOpenRingIdx, annulusBodyOpenRingIdx + 1):
        del xSampledAroundAlong[i][elementsAroundHalfDuod]
        del d1SampledAroundAlong[i][elementsAroundHalfDuod]
        del d2SampledAroundAlong[i][elementsAroundHalfDuod]
        del d3SampledAroundAlong[i][elementsAroundHalfDuod]

    # Remove multiple nodes at apex
    del xSampledAroundAlong[0][1:], d1SampledAroundAlong[0][1:], d2SampledAroundAlong[0][1:], \
        d3SampledAroundAlong[0][1:], d1CurvatureAroundAlong[0][1:], d2CurvatureAroundAlong[0][1:]

    # Set magnitude for d1 at apex
    arcLength = interp.computeCubicHermiteArcLength(xSampledAroundAlong[0][0], d1SampledAroundAlong[0][0],
                                                    xSampledAroundAlong[1][elementsAroundQuarterDuod],
                                                    d2SampledAroundAlong[1][elementsAroundQuarterDuod],
                                                    rescaleDerivatives=True)
    d1SampledAroundAlong[0][0] = set_magnitude(d1SampledAroundAlong[0][0], arcLength)
    d2SampledAroundAlong[0][0] = set_magnitude(d2SampledAroundAlong[0][0], arcLength)

    # Replace d1Curvature with d2Curvature
    d1CurvatureAroundAlong[0][0] = d2CurvatureAroundAlong[0][0]

    # for n2 in range(len(xSampledAroundAlong)):
    #     for n1 in range(len(xSampledAroundAlong[n2])):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1SampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2SampledAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3SampledAroundAlong[n2][n1])
    #         nodeIdentifier += 1

    # Create inner nodes
    xList = []
    d1List = []
    d2List = []
    d3List = []
    nodeIdx = stomachStartNode
    idxMat = []

    xDistal = []
    d1Distal = []
    d2Distal = []
    d3Distal = []
    xPrev = []
    d2Prev = []
    nodeIdxDistal = []

    if elementsCountThroughWall > 1:
        thicknessProportionsUI = [0.0, mucosaRelThickness, submucosaRelThickness, circularRelThickness,
                                  longitudinalRelThickness, longitudinalRelThickness]
        thicknessProportions = [thicknessProportion / sum(thicknessProportionsUI[:-1])
                                for thicknessProportion in thicknessProportionsUI]

        xi3List = []
        xi3 = 0.0
        for i in range(len(thicknessProportions) - 1):
            xi3 += thicknessProportions[i]
            xi3List.append(xi3)

    for n2 in range(len(xSampledAroundAlong)):
        idxThroughWall = []
        for n3 in range(elementsCountThroughWall + 1):
            xi3 = xi3List[n3] if elementsCountThroughWall > 1 else 1.0 / elementsCountThroughWall * n3
            idxAround = []
            xDistalAround = []
            d1DistalAround = []
            d2DistalAround = []
            d3DistalAround = []
            nodeIdxDistalAround = []

            for n1 in range(len(xSampledAroundAlong[n2])):
                # Coordinates
                norm = normalize(d3SampledAroundAlong[n2][n1])
                xOut = xSampledAroundAlong[n2][n1]
                xIn = [xOut[i] - norm[i] * wallThickness for i in range(3)]
                dWall = [wallThickness * c for c in norm]
                x = interp.interpolateCubicHermite(xIn, dWall, xOut, dWall, xi3)
                xList.append(x)

                # d1
                factor = 1.0 + wallThickness * (1.0 - xi3) * d1CurvatureAroundAlong[n2][n1]
                d1 = [factor * c for c in d1SampledAroundAlong[n2][n1]]
                d1List.append(d1)

                # d2
                factor = 1.0 + wallThickness * (1.0 - xi3) * d2CurvatureAroundAlong[n2][n1]
                d2 = [factor * c for c in d2SampledAroundAlong[n2][n1]]
                d2List.append(d2)

                # d3
                d3 = [c * wallThickness * (thicknessProportions[n3 + 1] if elementsCountThroughWall > 1 else 1.0)
                      for c in norm]
                d3List.append(d3)

                if n2 >= len(xSampledAroundAlong) - 2:
                    xDistalAround.append(x)
                    d1DistalAround.append(d1)
                    d2DistalAround.append(d2)
                    d3DistalAround.append(d3)
                    nodeIdxDistalAround.append(nodeIdx)

                idxAround.append(nodeIdx)
                nodeIdx += 1
            idxThroughWall.append(idxAround)

            if n2 == len(xSampledAroundAlong) - 2:
                xPrev.append(xDistalAround)
                d2Prev.append(d2DistalAround)

            if n2 == len(xSampledAroundAlong) - 1:
                xDistal.append(xDistalAround)
                d1Distal.append(d1DistalAround)
                d2Distal.append(d2DistalAround)
                d3Distal.append(d3DistalAround)
                nodeIdxDistal.append(nodeIdxDistalAround)
        idxMat.append(idxThroughWall)

    nodeIdxGC = []
    for n2 in range(len(idxMat)):
        for n3 in range(len(idxMat[n2])):
            if n2 == 0:
                nodeIdxGC += idxMat[n2][n3]
            else:
                nodeIdxGC.append(idxMat[n2][n3][0])

    for n2 in range(1, annulusFundusOpenRingIdx + 1):
        for n3 in range(len(idxMat[n2])):
            nodeIdxGC.append(idxMat[n2][n3][int(0.5 * len(xSampledAroundAlong[n2]))])

    nodeIdxLC = []
    for n2 in range(annulusBodyOpenRingIdx, len(xSampledAroundAlong)):
        for n3 in range(len(idxMat[n2])):
            nodeIdxLC.append(
                idxMat[n2][n3][elementsAroundHalfDuod])

    for n2 in range(len(xList)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1List[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2List[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3List[n2])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier += 1

    annotationGroupsAlong = []
    for i in range(len(elementsAlongSections)):
        elementsCount = elementsAlongSections[i]
        for n in range(elementsCount):
            annotationGroupsAlong.append(annotationGroupAlong[i])

    # Create elements
    fundusMucosaElementIdentifiers = []
    elementIdxMat = []
    n = 0
    for n2 in range(elementsAlongEsophagus):
        elementIdxThroughWall = []
        for n3 in range(elementsThroughEsophagusWall):
            elementIdxAround = []
            for n1 in range(elementsCountAroundEso):
                n += 1
                elementIdxAround.append(n)
            elementIdxThroughWall.append(elementIdxAround)
        elementIdxMat.append(elementIdxThroughWall)

    if useCubicHermiteThroughWall:
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    else:
        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eftStandard = eftfactory.createEftBasic()

    elementtemplateStandard = mesh.createElementtemplate()
    elementtemplateStandard.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplateStandard.defineField(coordinates, -1, eftStandard)

    elementtemplateX = mesh.createElementtemplate()
    elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    fundusElements = elementsAlongSections[0]
    radiansPerElementAroundDuod = math.pi * 2.0 / elementsCountAroundDuod

    for e2 in range(len(xSampledAroundAlong) - 1):
        elementIdxThroughWall = []
        if e2 == 0:  # pole
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(elementsCountAroundDuod):
                    va = e1
                    vb = (e1 + 1) % elementsCountAroundDuod
                    eft1 = eftfactory.createEftShellPoleBottom(va*100, vb*100)
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplateX)
                    bni1 = e3 + stomachStartNode
                    bni2 = (elementsCountThroughWall + 1) + stomachStartNode + elementsCountAroundDuod * e3 + e1
                    bni3 = (elementsCountThroughWall + 1) + stomachStartNode + elementsCountAroundDuod * e3 + \
                           (e1 + 1) % elementsCountAroundDuod
                    nodeIdentifiers = [bni1, bni2, bni3, bni1 + 1, bni2 + elementsCountAroundDuod,
                                       bni3 + elementsCountAroundDuod]

                    element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    # set general linear map coefficients
                    radiansAround = e1 * radiansPerElementAroundDuod
                    radiansAroundNext = ((e1 + 1) % elementsCountAroundDuod) * radiansPerElementAroundDuod
                    scalefactors = [
                        -1.0,
                        math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAroundDuod,
                        math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAroundDuod,
                        math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAroundDuod,
                        math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAroundDuod
                    ]
                    element.setScaleFactors(eft1, scalefactors)
                    if e2 < fundusElements and limitingRidge and elementsCountThroughWall > 1 and e3 == 0:
                        fundusMucosaElementIdentifiers.append(elementIdentifier)
                    elementIdxAround.append(elementIdentifier)
                    elementIdentifier += 1
                    annotationGroups = annotationGroupsAlong[e2] + annotationGroupsThroughWall[e3]
                    if annotationGroups:
                        allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                        for annotationGroup in annotationGroups:
                            meshGroup = annotationGroup.getMeshGroup(mesh)
                            meshGroup.addElement(element)
                elementIdxThroughWall.append(elementIdxAround)
            elementIdxMat.append(elementIdxThroughWall)

        else:
            e1Range = elementsCountAroundDuod - 2 if (annulusFundusOpenRingIdx - 1 <= e2 <= annulusBodyOpenRingIdx) \
                else len(xSampledAroundAlong[e2])
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(e1Range):
                    e1IdxBni1 = e1
                    e1IdxBni3 = e1
                    if e1 > elementsAroundHalfDuod - 2:
                        if e2 == annulusFundusOpenRingIdx - 1:
                            e1IdxBni1 = e1 + 2
                            e1IdxBni3 = e1 + 1
                        elif annulusFundusOpenRingIdx - 1 < e2 < annulusBodyOpenRingIdx:
                            e1IdxBni1 = e1 + 1
                            e1IdxBni3 = e1 + 1
                        elif e2 == annulusBodyOpenRingIdx:
                            e1IdxBni1 = e1 + 1
                            e1IdxBni3 = e1 + 2

                    eft1 = eftStandard
                    scaleFactors = []
                    elementtemplate1 = elementtemplateStandard
                    bni111 = idxMat[e2][e3][e1IdxBni1]
                    bni211 = idxMat[e2][e3][(e1IdxBni1 + 1) % len(idxMat[e2][e3])]
                    bni121 = idxMat[e2 + 1][e3][e1IdxBni3]
                    bni221 = idxMat[e2 + 1][e3][(e1IdxBni3 + 1) % len(idxMat[e2 + 1][e3])]
                    bni112 = idxMat[e2][e3 + 1][e1IdxBni1]
                    bni212 = idxMat[e2][e3 + 1][(e1IdxBni1 + 1) % len(idxMat[e2][e3])]
                    bni122 = idxMat[e2 + 1][e3 + 1][e1IdxBni3]
                    bni222 = idxMat[e2 + 1][e3 + 1][(e1IdxBni3 + 1) % len(idxMat[e2 + 1][e3])]
                    nodeIdentifiers = [bni111, bni211, bni121, bni221,
                                       bni112, bni212, bni122, bni222]

                    if e2 == annulusFundusOpenRingIdx - 2:
                        if e1 == elementsAroundHalfDuod - 2:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) # 145

                        elif e1 == elementsAroundHalfDuod - 1:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print('1', elementIdentifier) # 146

                        elif e1 == elementsAroundHalfDuod:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print('2', elementIdentifier) #147

                        elif e1 == elementsAroundHalfDuod + 1:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) #148

                    if e2 == annulusFundusOpenRingIdx - 1:
                        if e1 == elementsAroundHalfDuod - 2:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) # 165

                        elif e1 == elementsAroundHalfDuod - 1:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) # 166

                    elif (elementsAroundQuarterEso - 2) > 0 and \
                            annulusFundusOpenRingIdx <= e2 < annulusFundusOpenRingIdx + \
                            2.0 * (elementsAroundQuarterEso - 2):
                        if e1 == elementsAroundHalfDuod - 2:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) # 183, 201

                        elif e1 == elementsAroundHalfDuod - 1:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) # 184, 202

                    if e2 == annulusBodyOpenRingIdx:
                        if e1 == elementsAroundHalfDuod - 2:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) # 219

                        elif e1 == elementsAroundHalfDuod - 1:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) # 220

                    if e2 == annulusBodyOpenRingIdx + 1:
                        if e1 == elementsAroundHalfDuod - 2:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print('1', elementIdentifier) #237

                        elif e1 == elementsAroundHalfDuod + 1:
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                            # print(elementIdentifier) #240

                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    if scaleFactors:
                        element.setScaleFactors(eft1, scaleFactors)
                    if e2 < fundusElements and limitingRidge and elementsCountThroughWall > 1 and e3 == 0:
                        fundusMucosaElementIdentifiers.append(elementIdentifier)
                    elementIdxAround.append(elementIdentifier)
                    elementIdentifier += 1
                    annotationGroups = annotationGroupsAlong[e2] + annotationGroupsThroughWall[e3]
                    if annotationGroups:
                        allAnnotationGroups = mergeAnnotationGroups(allAnnotationGroups, annotationGroups)
                        for annotationGroup in annotationGroups:
                            meshGroup = annotationGroup.getMeshGroup(mesh)
                            meshGroup.addElement(element)
                elementIdxThroughWall.append(elementIdxAround)
            elementIdxMat.append(elementIdxThroughWall)

    # Annulus
    # Assemble endPoints for annulus
    endPoints_x = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    endPoints_d1 = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    endPoints_d2 = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    endNode_Id = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    endDerivativesMap = [[None] * elementsCountAroundEso for n3 in range(elementsCountThroughWall + 1)]
    endProportions = []

    for n3 in range(elementsCountThroughWall + 1):
        n1 = 0
        m = -1
        for nAround in range(elementsCountAroundEso):
            if nAround == 0:
                idx = idxMat[annulusFundusOpenRingIdx - 1][n3][elementsAroundHalfDuod]
            elif 0 < nAround < elementsAroundQuarterEso:
                idx = idxMat[annulusFundusOpenRingIdx - 1 + n1][n3][elementsAroundHalfDuod - 1]
                n1 += 1
            elif nAround == elementsAroundQuarterEso or nAround == elementsAroundQuarterEso + elementsAroundHalfEso:
                idx = idxMat[bodyStartIdx][n3][elementsAroundHalfDuod - (1 if nAround == elementsAroundQuarterEso
                                                                         else 0)]
                n1 = 1
            elif elementsAroundQuarterEso < nAround < elementsAroundHalfEso - 1:
                idx = idxMat[bodyStartIdx + n1][n3][elementsAroundHalfDuod - 1]
                n1 += 1
            elif elementsAroundHalfEso - 1 <= nAround <= elementsAroundHalfEso + 1:
                idx = idxMat[annulusBodyOpenRingIdx + 1][n3][elementsAroundHalfDuod + m]
                m += 1
                n1 = 1
            elif elementsAroundHalfEso + 1 < nAround < elementsAroundHalfEso + elementsAroundQuarterEso:
                idx = idxMat[annulusBodyOpenRingIdx + 1 - n1][n3][elementsAroundHalfDuod]
                n1 += 1
            elif elementsAroundHalfEso + elementsAroundQuarterEso < nAround < elementsCountAroundEso - 1:
                idx = idxMat[bodyStartIdx - n1][n3][elementsAroundHalfDuod]
                n1 += 1
            else:
                idx = idxMat[annulusFundusOpenRingIdx - 1][n3][elementsAroundHalfDuod + 1]

            endPoints_x[n3][nAround] = xList[idx - stomachStartNode]
            endPoints_d1[n3][nAround] = d1List[idx - stomachStartNode]
            endPoints_d2[n3][nAround] = d2List[idx - stomachStartNode]
            endNode_Id[n3][nAround] = idx

            if n3 == elementsCountThroughWall:  # outer layer
                startGuessPosition = trackSurfaceStomach.findNearestPositionParameter(endPoints_x[n3][nAround])[0]
                endPosition = trackSurfaceStomach.findNearestPosition(endPoints_x[n3][nAround], startGuessPosition)
                endProportions.append(trackSurfaceStomach.getProportion(endPosition))

    for n3 in range(elementsCountThroughWall + 1):
        for nAround in range(elementsCountAroundEso):
            endDerivativesMap[n3][nAround] = (None, None, None)

    startProportions = []
    for n in range(elementsCountAroundEso):
        startProportions.append(trackSurfaceStomach.getProportion(o1_Positions[n]))

    cardiaGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region,
                                                     get_stomach_term("cardia of stomach"))
    cardiaMeshGroup = cardiaGroup.getMeshGroup(mesh)
    if cardiaGroup not in allAnnotationGroups:
        allAnnotationGroups.append(cardiaGroup)

    lastDuodenumElementIdentifier = elementIdentifier

    stomachWallAnnotationGroups = []
    if elementsCountThroughWall == 4:
        stomachWallAnnotationGroups = [[mucosaGroup], [submucosaGroup], [circularMuscleGroup],
                                       [longitudinalMuscleGroup]]

    # Remove mucosa layer from annulus
    if elementsCountThroughWall == 4 and limitingRidge:
        o1_x = o1_x[1:]
        o1_d1 = o1_d1[1:]
        o1_d2 = o1_d2[1:]
        o1_NodeId = o1_NodeId[1:]
        endPoints_x = endPoints_x[1:]
        endPoints_d1 = endPoints_d1[1:]
        endPoints_d2 = endPoints_d2[1:]
        endNode_Id = endNode_Id[1:]
        endDerivativesMap = endDerivativesMap[1:]
        stomachWallAnnotationGroups = stomachWallAnnotationGroups[1:]

    nextNodeIdentifier, nextElementIdentifier = createAnnulusMesh3d(
        nodes, mesh, nodeIdentifier, elementIdentifier,
        o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
        endPoints_x, endPoints_d1, endPoints_d2, None, endNode_Id, endDerivativesMap,
        elementsCountRadial=elementsCountAcrossCardia, meshGroups=[stomachMeshGroup, cardiaMeshGroup],
        wallAnnotationGroups=stomachWallAnnotationGroups,
        tracksurface=trackSurfaceStomach,
        startProportions=startProportions, endProportions=endProportions,
        rescaleStartDerivatives=True, rescaleEndDerivatives=True, sampleBlend=0.0, fixMinimumStart=True,
        coordinates=coordinates)

    elementIdxThroughWall = []
    n = lastDuodenumElementIdentifier - 1
    for n3 in range(elementsCountThroughWall):
        elementIdxAround = []
        for n1 in range(elementsCountAroundEso):
            n += 1
            elementIdxAround.append(n)
        elementIdxThroughWall.append(elementIdxAround)
    elementIdxMat.append(elementIdxThroughWall)

    # delete mucosa layer in fundus when there is a limiting ridge
    mesh_destroy_elements_and_nodes_by_identifiers(mesh, fundusMucosaElementIdentifiers)

    # Create annotation groups for dorsal and ventral parts of the stomach
    dorsalGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("dorsal stomach"))
    ventralGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("ventral stomach"))
    dorsalMeshGroup = dorsalGroup.getMeshGroup(mesh)
    ventralMeshGroup = ventralGroup.getMeshGroup(mesh)

    for e2 in range(len(elementIdxMat)):
        for e3 in range(len(elementIdxMat[e2])):
            for e1 in range(len(elementIdxMat[e2][e3])):
                elementIdx = elementIdxMat[e2][e3][e1]
                element = mesh.findElementByIdentifier(elementIdx)
                if e1 < 0.5 * len(elementIdxMat[e2][e3]):
                    ventralMeshGroup.addElement(element)
                else:
                    dorsalMeshGroup.addElement(element)
    if dorsalGroup not in allAnnotationGroups:
        allAnnotationGroups.append(dorsalGroup)
    if ventralGroup not in allAnnotationGroups:
        allAnnotationGroups.append(ventralGroup)

    nodesOnLCMargin = []
    for n2 in range(elementsAlongEsophagus + 1):
        for n3 in range(elementsThroughEsophagusWall + 1):
            nodeIdxOnLCMargin = 1 + elementsAroundHalfEso + \
                                n2 * (elementsThroughEsophagusWall + 1) * elementsCountAroundEso + \
                                n3 * elementsCountAroundEso
            nodesOnLCMargin.append(nodeIdxOnLCMargin)
    allNodesOnLC = nodesOnLCMargin + nodeIdxLC

    nearLCGroup = AnnotationGroup(region, ("elements adjacent to lesser curvature", "None"))

    elementIter = mesh.createElementiterator()
    element = elementIter.next()
    while element.isValid():
        eft = element.getElementfieldtemplate(coordinates, -1)
        nodeIdentifiers = get_element_node_identifiers(element, eft)
        for n in range(len(nodeIdentifiers)):
            if nodeIdentifiers[n] in allNodesOnLC:
                nearLCGroup.getMeshGroup(mesh).addElement(element)
                break
        element = elementIter.next()
    allAnnotationGroups.append(nearLCGroup)

    return allAnnotationGroups, nextNodeIdentifier, nextElementIdentifier, elementsAlongSections, nodeIdxDistal, \
           xDistal, d1Distal, d2Distal, d3Distal, arclengthDuodenumCP, xPrev, d2Prev


class CustomNetworkLayout:
    """
    Generates sampled network layout for part of network layout.
    """
    def __init__(self, cx, cd1, cd2, cd3, cd12, cd13):
        self.cxPath = cx
        self.cd1Path = cd1
        self.cd2Path = cd2
        self.cd3Path = cd3
        self.cd12Path = cd12
        self.cd13Path = cd13
