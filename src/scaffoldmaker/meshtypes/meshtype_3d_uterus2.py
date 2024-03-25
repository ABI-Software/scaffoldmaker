"""
Generates a 3-D uterus mesh from a 1-D network layout, with variable
numbers of elements around, along and through wall.
"""

import math

from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.utils.zinc.finiteelement import get_element_node_identifiers
from cmlibs.utils.zinc.general import ChangeManager
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, getAnnotationGroupForTerm, \
    findOrCreateAnnotationGroupForTerm
from scaffoldmaker.annotation.uterus_terms import get_uterus_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.bifurcation import SegmentTubeData, \
    TubeBifurcationData, generateTube, generateTubeBifurcation, generateCurveMesh
from scaffoldmaker.utils.networkmesh import getPathRawTubeCoordinates, resampleTubeCoordinates
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters
from scaffoldmaker.utils.zinc_utils import get_nodeset_path_ordered_field_parameters

def getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName):
    assert parameterSetName in cls.getParameterSetNames() # make sure parameter set is in list of parameters of parent scaffold
    if parameterSetName in ("Default", "Human 1"):
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3, 3-4-5-11.1, 6-7-8, 8-9-10-11.2, 11.3-12-13-14,14-15-16,16-17-18",
                "Define inner coordinates": True,
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[                
                (1, [[-9.16, 10.27,-3.09], [ 0.89, 0.73,-1.15], [ 0.04,-0.23,-0.11], [ 0.08, 0.02, 0.01], [-0.20, 0.03,-0.14], [ 0.08,-0.06,-0.09]]),
                (2, [[-8.06, 10.64,-4.32], [ 1.27,-0.02,-1.26], [-0.00,-0.25, 0.00], [-0.10, 0.01, 0.14], [-0.17,-0.00,-0.17], [-0.02,-0.00, 0.02]]),
                (3, [[-6.71, 10.19,-5.51], [ 1.73,-0.52,-0.85], [-0.06,-0.25, 0.03], [-0.11,-0.25,-0.05], [-0.12,-0.00,-0.24], [ 0.12, 0.01,-0.37]]),
                (4, [[-4.71, 9.65,-5.89], [ 2.18,-0.21,-0.35], [-0.09,-0.77,-0.10], [ 0.27,-1.14,-0.57], [-0.11, 0.11,-0.78], [-0.31, 0.46,-0.76]]),
                (5, [[-2.42, 9.80,-6.20], [ 2.36, 0.21,-0.17], [ 0.13,-2.36,-1.07], [ 0.04,-0.93,-0.64], [-0.20, 0.80,-1.80], [ 0.22, 0.68,-0.92]]),
                (6, [[ 9.16, 10.27,-3.09], [-0.89, 0.73,-1.15], [-0.04,-0.23,-0.11], [-0.08, 0.02, 0.01], [-0.20,-0.03, 0.14], [ 0.08, 0.06, 0.09]]),
                (7, [[ 8.06, 10.64,-4.32], [-1.27,-0.02,-1.26], [ 0.00,-0.25, 0.00], [ 0.10, 0.01, 0.14], [-0.17, 0.00, 0.17], [-0.02, 0.00,-0.02]]),
                (8, [[ 6.71, 10.19,-5.51], [-1.73,-0.52,-0.85], [ 0.06,-0.25, 0.03], [ 0.11,-0.25,-0.05], [-0.12, 0.00, 0.24], [ 0.12,-0.01, 0.37]]),
                (9, [[ 4.71, 9.65,-5.89], [-2.18,-0.21,-0.35], [ 0.09,-0.77,-0.10], [-0.27,-1.14,-0.57], [-0.11,-0.11, 0.78], [-0.31,-0.46, 0.76]]),
                (10, [[ 2.42, 9.80,-6.20], [-2.36, 0.21,-0.17], [-0.13,-2.36,-1.07], [-0.04,-0.93,-0.64], [-0.20,-0.80, 1.80], [ 0.22,-0.68, 0.92]]),
                (11, [[-0.00, 10.08,-6.23], [[ 2.47, 0.34, 0.11],[-2.47, 0.34, 0.11],[ 0.02,-0.35, 2.81]], [[ 0.42,-2.59,-1.43],[-0.42,-2.59,-1.43],[ 0.03, 2.93, 0.36]], [[ 0.19, 0.45,-0.09],[-0.19, 0.45,-0.09],[ 0.01,-0.22,-0.05]], [[-0.08, 1.42,-2.59],[-0.08,-1.42, 2.59],[-3.52, 0.03, 0.03]], [[-0.16, 0.54,-0.68],[-0.16,-0.54, 0.68],[ 0.47, 0.03,-0.12]]]),
                (12, [[-0.00, 9.50,-3.62], [ 0.00,-0.81, 2.39], [ 0.02, 2.37, 0.80], [-0.02,-0.82, 0.29], [-3.08, 0.02, 0.01], [ 0.41,-0.03, 0.02]]),
                (13, [[-0.00, 8.50,-1.48], [ 0.00,-1.16, 1.77], [-0.01, 1.33, 0.87], [-0.01,-0.82, 0.10], [-2.70,-0.01,-0.01], [ 0.29,-0.01,-0.01]]),
                (14, [[-0.00, 7.27,-0.08], [ 0.00,-1.09, 0.82], [ 0.00, 0.75, 0.99], [ 0.01,-0.49, 0.09], [-2.51, 0.00, 0.00], [ 0.18, 0.01, 0.01]]),
                (15, [[-0.00, 6.50, 0.28], [ 0.00,-0.83, 0.18], [ 0.00, 0.22, 1.00], [ 0.00,-0.30,-0.20], [-2.34, 0.00, 0.00], [ 0.22, 0.00, 0.00]]),
                (16, [[-0.00, 5.67, 0.28], [ 0.00,-1.58,-0.05], [ 0.00,-0.03, 0.81], [ 0.00, 0.11,-0.26], [-2.07,-0.00, 0.00], [ 0.63, 0.00, 0.00]]),
                (17, [[-0.00, 3.35, 0.14], [ 0.00,-2.85,-0.14], [ 0.00,-0.02, 0.46], [ 0.00, 0.02,-0.13], [-1.08,-0.00, 0.00], [ 0.76, 0.00, 0.00]]),
                (18, [[ 0.00,-0.03, 0.00], [ 0.00,-3.91,-0.13], [ 0.00,-0.02, 0.55], [-0.00,-0.01, 0.31], [-0.55,-0.00, 0.00], [ 0.30,-0.00, 0.00]])], [
                (1, [[-9.16,10.27,-3.09], [0.89,0.73,-1.15], [0.02,-0.11,-0.06], [0.04,0.01,0.00], [-0.10,0.02,-0.07], [0.04,-0.03,-0.05]]), 
                (2, [[-8.06,10.64,-4.32], [1.27,-0.02,-1.26], [0.00,-0.12,0.00], [-0.05,0.00,0.07], [-0.08,0.00,-0.09], [-0.01,0.00,0.01]]), 
                (3, [[-6.71,10.19,-5.51], [1.73,-0.52,-0.85], [-0.03,-0.13,0.01], [-0.05,-0.13,-0.03], [-0.06,0.00,-0.12], [0.06,0.00,-0.19]]), 
                (4, [[-4.71,9.65,-5.89], [2.18,-0.21,-0.35], [-0.05,-0.38,-0.05], [0.13,-0.57,-0.28], [-0.06,0.06,-0.39], [-0.15,0.23,-0.38]]), 
                (5, [[-2.42,9.80,-6.20], [2.36,0.21,-0.17], [0.07,-1.18,-0.54], [0.02,-0.53,-0.26], [-0.10,0.40,-0.90], [0.11,0.27,-0.49]]), 
                (6, [[9.16,10.27,-3.09], [-0.89,0.73,-1.15], [-0.02,-0.11,-0.06], [-0.04,0.01,0.00], [-0.10,-0.02,0.07], [0.04,0.03,0.05]]), 
                (7, [[8.06,10.64,-4.32], [-1.27,-0.02,-1.26], [0.00,-0.12,0.00], [0.05,0.00,0.07], [-0.08,0.00,0.09], [-0.01,0.00,-0.01]]), 
                (8, [[6.71,10.19,-5.51], [-1.73,-0.52,-0.85], [0.03,-0.13,0.01], [0.05,-0.13,-0.03], [-0.06,0.00,0.12], [0.06,0.00,0.19]]), 
                (9, [[4.71,9.65,-5.89], [-2.18,-0.21,-0.35], [0.05,-0.38,-0.05], [-0.13,-0.57,-0.28], [-0.06,-0.06,0.39], [-0.15,-0.23,0.38]]), 
                (10, [[2.42,9.80,-6.20], [-2.36,0.21,-0.17], [-0.07,-1.18,-0.54], [-0.02,-0.53,-0.26], [-0.10,-0.40,0.90], [0.11,-0.27,0.49]]),
                (11, [[-0.00,10.08,-6.23], [[2.47,0.34,0.11],[-2.47,0.34,0.11],[0.02,-0.35,2.81]], [[0.21,-1.29,-0.72],[-0.21,-1.29,-0.72],[0.01,1.46,0.18]], [[0.09,0.23,-0.05],[-0.09,0.23,-0.05],[0.00,-0.11,-0.02]], [[-0.04,0.71,-1.30],[-0.04,-0.71,1.30],[-1.76,0.01,0.01]], [[-0.08,0.27,-0.34],[-0.08,-0.27,0.34],[0.23,0.02,-0.06]]]),                (12, [[-0.00,9.50,-3.62], [0.00,-0.81,2.39], [0.01,1.18,0.40], [-0.01,-0.38,0.15], [-1.54,0.01,0.00], [0.20,-0.01,0.01]]),
                (13, [[-0.00,8.50,-1.48], [0.00,-1.16,1.77], [-0.00,0.67,0.44], [0.00,-0.39,0.06], [-1.35,0.00,0.00], [0.59,0.00,0.00]]), 
                (14, [[-0.00,7.27,-0.08], [0.00,-1.09,0.82], [0.02,0.40,0.53], [0.00,-0.26,0.06], [-0.36,0.01,0.01], [0.47,0.00,0.00]]), 
                (15, [[-0.00,6.50,0.28], [0.00,-0.83,0.18], [-0.00,0.11,0.50], [-0.01,-0.15,-0.12], [-0.41,0.00,0.00], [0.00,0.00,0.00]]), 
                (16, [[-0.00,5.67,0.28], [0.00,-1.58,-0.05], [0.00,-0.01,0.41], [-0.00,0.06,-0.13], [-0.36,0.00,0.00], [-0.07,0.00,0.00]]), 
                (17, [[-0.00,3.35,0.14], [0.00,-2.85,-0.14], [0.00,-0.01,0.23], [-0.00,0.01,-0.06], [-0.54,0.00,0.00], [0.04,0.00,0.00]]), 
                (18, [[0.00,-0.03,0.00], [0.00,-3.91,-0.13], [0.00,-0.01,0.28], [0.00,-0.01,0.16], [-0.28,0.00,0.00], [0.48,0.00,0.00]])]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-17',
                    'name': get_uterus_term('uterus')[0],
                    'ontId': get_uterus_term('uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2, 6-7',
                    'name': get_uterus_term('fallopian tube')[0],
                    'ontId': get_uterus_term('fallopian tube')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-5, 8-13',
                    'name': get_uterus_term('body of uterus')[0],
                    'ontId': get_uterus_term('body of uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '14-15',
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '16-17',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-10',
                    'name': 'pre-bifurcation segments',
                    'ontId': 'None'
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-17',
                    'name': 'post-bifurcation segments',
                    'ontId': 'None'
                }
            ]
        })
    elif "Mouse 1" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3-4-5, 5-6.1, 7-8-9-10-11, 11-6.2, 6.3-12-13-14, 14-15, 15-16",
                "Define inner coordinates": True,
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ["coordinates", "inner coordinates"],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[                    
                (1, [[-24.20,54.56,0.00], [2.53,-14.75,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.42,0.41,0.00], [0.14,-0.02,0.00]]),
                (2, [[-20.88,41.34,0.00], [4.10,-11.62,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.47,0.87,0.00], [0.05,0.43,0.00]]),
                (3, [[-16.28,31.37,0.00], [5.76,-8.94,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.38,1.54,0.00], [-0.55,0.82,0.00]]),
                (4, [[-9.59,23.63,0.00], [6.74,-6.07,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [1.74,1.94,0.00], [-0.66,0.45,0.00]]),
                (5, [[-3.14,19.13,0.00], [4.91,-2.86,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [1.31,2.26,0.00], [-0.29,0.21,0.00]]),
                (6, [[0.00,17.63,-0.00], [[1.29,-0.13,-0.00],[-1.29,-0.13,-0.00],[0.00,-3.86,-0.00]], [[0.00,0.00,-2.50],[0.00,0.00,-2.50],[0.00,0.00,2.50]], [[-0.00,-0.00,-0.00],[-0.00,0.00,-0.00],[0.00,0.00,0.00]], [[0.28,2.81,0.00],[0.28,-2.81,0.00],[-2.67,-0.00,0.00]], [[-3.19,1.05,-0.00],[-3.19,-1.05,-0.00],[0.25,0.00,0.00]]]),
                (7, [[24.20,54.56,0.00], [-2.53,-14.75,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.42,-0.41,0.00], [0.14,0.02,0.00]]),
                (8, [[20.88,41.34,0.00], [-4.10,-11.62,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.47,-0.87,0.00], [0.05,-0.43,0.00]]),
                (9, [[16.28,31.37,0.00], [-5.76,-8.94,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [2.38,-1.54,0.00], [-0.55,-0.82,0.00]]),
                (10, [[9.59,23.63,0.00], [-6.74,-6.07,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [1.74,-1.94,0.00], [-0.66,-0.45,0.00]]),
                (11, [[3.14,19.13,0.00], [-4.91,-2.86,-0.00], [0.00,0.00,-2.50], [0.00,0.00,0.00], [1.31,-2.26,0.00], [-0.29,-0.21,0.00]]),
                (12, [[0.00,13.95,-0.00], [0.00,-3.50,-0.00], [0.00,0.00,2.50], [0.00,0.00,-0.00], [-2.50,-0.00,0.00], [0.09,0.00,0.00]]),
                (13, [[0.00,10.63,-0.00], [0.00,-3.47,-0.00], [0.00,0.00,2.50], [0.00,0.00,-0.00], [-2.50,-0.00,0.00], [0.00,0.00,0.00]]),
                (14, [[0.00,7.02,-0.00], [0.00,-3.31,-0.00], [0.00,0.00,2.50], [0.00,0.00,-0.00], [-2.50,-0.00,0.00], [0.00,0.00,0.00]]),
                (15, [[0.00,4.00,-0.00], [0.00,-3.51,-0.00], [0.00,0.00,2.50], [0.00,0.00,-0.00], [-2.50,-0.00,0.00], [0.00,0.00,0.00]]),
                (16, [[0.00,0.00,-0.00], [0.00,-4.49,-0.00], [0.00,0.00,2.50], [-0.00,-0.00,0.00], [-2.50,-0.00,0.00], [-0.00,-0.00,-0.00]])],[
                (1, [[-24.20,54.56,0.00], [2.53,-14.75,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.45,0.25,0.00], [0.08,-0.01,0.00]]),
                (2, [[-20.88,41.34,0.00], [4.10,-11.62,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.48,0.52,0.00], [0.03,0.26,0.00]]),
                (3, [[-16.28,31.37,0.00], [5.76,-8.94,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.43,0.92,0.00], [-0.33,0.50,0.00]]),
                (4, [[-9.59,23.63,0.00], [6.74,-6.07,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.05,1.16,0.00], [-0.40,0.27,0.00]]),
                (5, [[-3.14,19.13,0.00], [4.91,-2.86,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [0.79,1.35,0.00], [-0.18,0.12,0.00]]),
                (6, [[0.00,17.63,-0.00], [[1.29,-0.13,-0.00],[-1.29,-0.13,-0.00],[0.00,-3.86,-0.00]], [[0.00,0.00,-1.50],[0.00,0.00,-1.50],[0.00,0.00,1.50]], [[-0.00,-0.00,-0.00],[-0.00,0.00,-0.00],[0.00,0.00,0.00]], [[0.17,1.68,0.00],[0.17,-1.68,0.00],[-1.60,-0.00,0.00]], [[-1.91,0.62,-0.00],[-1.91,-0.62,-0.00],[0.15,0.00,0.00]]]),
                (7, [[24.20,54.56,0.00], [-2.53,-14.75,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.45,-0.25,0.00], [0.08,0.01,0.00]]),
                (8, [[20.88,41.34,0.00], [-4.10,-11.62,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.48,-0.52,0.00], [0.03,-0.26,0.00]]),
                (9, [[16.28,31.37,0.00], [-5.76,-8.94,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.43,-0.92,0.00], [-0.33,-0.50,0.00]]),
                (10, [[9.59,23.63,0.00], [-6.74,-6.07,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [1.05,-1.16,0.00], [-0.40,-0.27,0.00]]),
                (11, [[3.14,19.13,0.00], [-4.91,-2.86,-0.00], [0.00,0.00,-1.50], [0.00,0.00,0.00], [0.79,-1.35,0.00], [-0.18,-0.12,0.00]]),
                (12, [[0.00,13.95,-0.00], [0.00,-3.50,-0.00], [0.00,0.00,1.50], [0.00,0.00,0.00], [-1.50,-0.00,0.00], [0.05,0.00,0.00]]),
                (13, [[0.00,10.63,-0.00], [0.00,-3.47,-0.00], [0.00,0.00,1.50], [0.00,0.00,0.00], [-1.50,-0.00,0.00], [-0.00,0.00,0.00]]),
                (14, [[0.00,7.02,-0.00], [0.00,-3.31,-0.00], [0.00,0.00,1.50], [0.00,0.00,0.00], [-1.50,-0.00,0.00], [-0.00,0.00,0.00]]),
                (15, [[0.00,4.00,-0.00], [0.00,-3.51,-0.00], [0.00,0.00,1.50], [0.00,0.00,0.00], [-1.50,-0.00,0.00], [0.00,0.00,0.00]]),
                (16, [[0.00,0.00,-0.00], [0.00,-4.49,-0.00], [0.00,0.00,1.50], [-0.00,-0.00,-0.00], [-1.50,-0.00,0.00], [-0.00,-0.00,-0.00]])]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-15',
                    'name': get_uterus_term('uterus')[0],
                    'ontId': get_uterus_term('uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_uterus_term('left uterine horn')[0],
                    'ontId': get_uterus_term('left uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6-9',
                    'name': get_uterus_term('right uterine horn')[0],
                    'ontId': get_uterus_term('right uterine horn')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5,10-13',
                    'name': get_uterus_term('body of uterus')[0],
                    'ontId': get_uterus_term('body of uterus')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '14',
                    'name': get_uterus_term('uterine cervix')[0],
                    'ontId': get_uterus_term('uterine cervix')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '15',
                    'name': get_uterus_term('vagina')[0],
                    'ontId': get_uterus_term('vagina')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-10',
                    'name': 'pre-bifurcation segments',
                    'ontId': 'None'
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '11-15',
                    'name': 'post-bifurcation segments',
                    'ontId': 'None'
                },
            ]
        })


class MeshType_3d_uterus2(Scaffold_base):
    """
    Generates a 3-D uterus mesh from a 1-D network layout with variable numbers of elements around, along and through
    wall.
    Magnitude of D2 and D3 are the radii of the uterus in the respective directions.
    """

    @staticmethod
    def getName():
        return '3D Uterus 2'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):

        options = {
            'Base parameter set': parameterSetName,
            'Network layout': getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName),
            'Number of elements around': 10,
            'Number of elements around horns': 8,
            'Number of elements through wall': 1,
            "Annotation elements counts around": [0],
            'Target element density along longest segment': 5.0,
            "Serendipity": True,
            'Use linear through wall': True,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements along': 4,
            'Refine number of elements around': 4,
            'Refine number of elements through wall': 1
        }
        if 'Mouse' in parameterSetName:
            options['Number of elements around'] = 8
            options['Target element density along longest segment'] = 10

        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Network layout',
            'Number of elements around',
            'Number of elements around horns',
            'Number of elements through wall',
            'Target element density along longest segment',
            "Serendipity",
            'Refine',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall']
        return optionNames

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
        for key in [
            'Number of elements through wall',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall'
        ]:
            if options[key] < 1:
                options[key] = 1

        for key in [
            'Number of elements around',
            'Number of elements around horns']:
            if options[key] < 4:
                options[key] = 4
            if options[key] % 2 > 0:
                options[key] = options[key] // 2 * 2
        if options["Number of elements through wall"] < 1:
            options["Number of elements through wall"] = 1

        if options["Target element density along longest segment"] < 1.0:
            options["Target element density along longest segment"] = 1.0
        dependentChanges = False
        return dependentChanges

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        parameterSetName = options['Base parameter set']
        isHuman = parameterSetName in ("Default", "Human 1")
        networkLayout = options['Network layout']
        elementsCountAroundPostBifurcation = options['Number of elements around']
        elementsCountAroundPreBifurcation = options['Number of elements around horns']
        elementsCountThroughWall = options['Number of elements through wall']
        defaultElementsCountAround = elementsCountAroundPostBifurcation
        annotationElementsCountsAround = [0, 0, elementsCountAroundPostBifurcation, elementsCountAroundPreBifurcation]
        targetElementDensityAlongLongestSegment = options['Target element density along longest segment']
        useCrossDerivatives = options['Use cross derivatives']
        serendipity = options["Serendipity"]
        showIntersectionCurves = False
        showTrimSurfaces = False # options["Show trim surfaces"]

        # Geometric coordinates
        layoutRegion = region.createRegion()
        networkLayout.generate(layoutRegion)  # ask scaffold to generate to get user-edited parameters
        layoutAnnotationGroups = networkLayout.getAnnotationGroups()

        fieldmodule = region.getFieldmodule()
        nodes = fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        coordinates = findOrCreateFieldCoordinates(fieldmodule)
        nodeIdentifier = 1
        elementIdentifier = 1

        networkMesh = networkLayout.getConstructionObject()

        layoutRegion = networkMesh.getRegion()
        layoutFieldmodule = layoutRegion.getFieldmodule()
        layoutNodes = layoutFieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        layoutCoordinates = layoutFieldmodule.findFieldByName("coordinates").castFiniteElement()
        layoutInnerCoordinates = layoutFieldmodule.findFieldByName("inner coordinates").castFiniteElement()
        if not layoutInnerCoordinates.isValid():
            layoutInnerCoordinates = None
        dimension = 3 if layoutInnerCoordinates else 2
        layoutMesh = layoutFieldmodule.findMeshByDimension(1)
        assert (elementsCountThroughWall == 1) or (layoutInnerCoordinates and (elementsCountThroughWall >= 1))

        fieldmodule = region.getFieldmodule()
        mesh = fieldmodule.findMeshByDimension(dimension)
        fieldcache = fieldmodule.createFieldcache()

        # make 2D annotation groups from 1D network layout annotation groups
        annotationGroups = []
        layoutAnnotationMeshGroupMap = []  # List of tuples of layout annotation mesh group to final mesh group
        for layoutAnnotationGroup in layoutAnnotationGroups:
            if layoutAnnotationGroup.getDimension() == 1:
                annotationGroup = AnnotationGroup(region, layoutAnnotationGroup.getTerm())
                annotationGroups.append(annotationGroup)
                layoutAnnotationMeshGroupMap.append(
                    (layoutAnnotationGroup.getMeshGroup(layoutMesh), annotationGroup.getMeshGroup(mesh)))

        valueLabels = [
            Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2,
            Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3]

        networkSegments = networkMesh.getNetworkSegments()

        bodyGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("body of uterus"))
        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)
        cervixGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("uterine cervix"))
        cervixMeshGroup = cervixGroup.getMeshGroup(mesh)
        uterusGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("uterus"))

        # map from NetworkSegment to SegmentTubeData
        outerSegmentTubeData = {}
        innerSegmentTubeData = {} if layoutInnerCoordinates else None

        longestSegmentLength = 0.0
        for networkSegment in networkSegments:
            pathParameters = get_nodeset_path_ordered_field_parameters(
                layoutNodes, layoutCoordinates, valueLabels,
                networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())
            elementsCountAround = defaultElementsCountAround
            i = 0
            for layoutAnnotationGroup in layoutAnnotationGroups:
                if i >= len(annotationElementsCountsAround):
                    break
                annotationElementsCountAround = annotationElementsCountsAround[i]
                if annotationElementsCountAround > 0:
                    if networkSegment.hasLayoutElementsInMeshGroup(layoutAnnotationGroup.getMeshGroup(layoutMesh)):
                        elementsCountAround = annotationElementsCountAround
                        break
                i += 1
            outerSegmentTubeData[networkSegment] = tubeData = SegmentTubeData(pathParameters, elementsCountAround)
            px, pd1, pd2, pd12 = getPathRawTubeCoordinates(pathParameters, elementsCountAround)
            tubeData.setRawTubeCoordinates((px, pd1, pd2, pd12))

            segmentLength = tubeData.getSegmentLength()
            if segmentLength > longestSegmentLength:
                longestSegmentLength = segmentLength

            if layoutInnerCoordinates:
                innerPathParameters = get_nodeset_path_ordered_field_parameters(
                    layoutNodes, layoutInnerCoordinates, valueLabels,
                    networkSegment.getNodeIdentifiers(), networkSegment.getNodeVersions())
                px, pd1, pd2, pd12 = getPathRawTubeCoordinates(innerPathParameters, elementsCountAround)
                innerSegmentTubeData[networkSegment] = innerTubeData = SegmentTubeData(innerPathParameters,
                                                                                       elementsCountAround)
                innerTubeData.setRawTubeCoordinates((px, pd1, pd2, pd12))

            for layoutAnnotationMeshGroup, annotationMeshGroup in layoutAnnotationMeshGroupMap:
                if networkSegment.hasLayoutElementsInMeshGroup(layoutAnnotationMeshGroup):
                    tubeData.addAnnotationMeshGroup(annotationMeshGroup)
                    if layoutInnerCoordinates:
                        innerTubeData.addAnnotationMeshGroup(annotationMeshGroup)

        if longestSegmentLength > 0.0:
            targetElementLength = longestSegmentLength / targetElementDensityAlongLongestSegment
        else:
            targetElementLength = 1.0

        # map from NetworkNodes to bifurcation data, resample tube coordinates to fit bifurcation
        outerNodeTubeBifurcationData = {}
        innerNodeTubeBifurcationData = {} if layoutInnerCoordinates else None

        allSegmentTubeData = [outerSegmentTubeData]
        if layoutInnerCoordinates:
            allSegmentTubeData.append(innerSegmentTubeData)

        for segmentTubeData in allSegmentTubeData:
            nodeTubeBifurcationData = innerNodeTubeBifurcationData if (segmentTubeData is innerSegmentTubeData) else \
                outerNodeTubeBifurcationData
            # count = 0
            for networkSegment in networkSegments:
                tubeData = segmentTubeData[networkSegment]
                rawTubeCoordinates = tubeData.getRawTubeCoordinates()
                segmentNodes = networkSegment.getNetworkNodes()
                startSegmentNode = segmentNodes[0]
                startTubeBifurcationData = nodeTubeBifurcationData.get(startSegmentNode)
                startSurface = None
                newBifurcationData = []
                if not startTubeBifurcationData:
                    startInSegments = startSegmentNode.getInSegments()
                    startOutSegments = startSegmentNode.getOutSegments()
                    if ((len(startInSegments) + len(startOutSegments)) == 3):
                        # print("create start", networkSegment, startSegmentNode)
                        startTubeBifurcationData = TubeBifurcationData(startInSegments, startOutSegments, segmentTubeData)
                        nodeTubeBifurcationData[startSegmentNode] = startTubeBifurcationData
                        newBifurcationData.append(startTubeBifurcationData)

                if startTubeBifurcationData:
                    startSurface = startTubeBifurcationData.getSegmentTrimSurface(networkSegment)
                endSegmentNode = segmentNodes[-1]
                endTubeBifurcationData = nodeTubeBifurcationData.get(endSegmentNode)
                endSurface = None
                createEndBifurcationData = not endTubeBifurcationData
                if createEndBifurcationData:
                    endInSegments = endSegmentNode.getInSegments()
                    endOutSegments = endSegmentNode.getOutSegments()
                    if ((len(endInSegments) + len(endOutSegments)) == 3):
                        # print("create end", networkSegment, endSegmentNode.getNodeIdentifier(),
                        #       len(endInSegments), len(endOutSegments))
                        endTubeBifurcationData = TubeBifurcationData(endInSegments, endOutSegments, segmentTubeData)
                        nodeTubeBifurcationData[endSegmentNode] = endTubeBifurcationData
                        newBifurcationData.append(endTubeBifurcationData)
                if endTubeBifurcationData:
                    endSurface = endTubeBifurcationData.getSegmentTrimSurface(networkSegment)
                if segmentTubeData is outerSegmentTubeData:
                    segmentLength = tubeData.getSegmentLength()
                    # Previous code setting number of elements along to satisfy targetElementAspectRatio
                    # ringCount = len(rawTubeCoordinates[0])
                    # sumRingLength = 0.0
                    # for n in range(ringCount):
                    #     ringLength = getCubicHermiteCurvesLength(rawTubeCoordinates[0][n], rawTubeCoordinates[1][n], loop=True)
                    #     sumRingLength += ringLength
                    # meanElementLengthAround = sumRingLength / (ringCount * elementsCountAround)
                    # targetElementLength = targetElementAspectRatio * meanElementLengthAround
                    elementsCountAlong = max(1, math.ceil(segmentLength / targetElementLength))
                    loop = (len(startSegmentNode.getInSegments()) == 1) and \
                           (startSegmentNode.getInSegments()[0] is networkSegment) and \
                           (networkSegment.getNodeVersions()[0] == networkSegment.getNodeVersions()[-1])
                    if (elementsCountAlong == 1) and (startTubeBifurcationData or endTubeBifurcationData):
                        # at least 2 segments if bifurcating at either end, or loop
                        elementsCountAlong = 1
                    elif (elementsCountAlong < 3) and loop:
                        # at least 3 segments around loop; 2 should work, but zinc currently makes incorrect faces
                        elementsCountAlong = 3
                else:
                    # must match count from outer surface!
                    outerTubeData = outerSegmentTubeData[networkSegment]
                    elementsCountAlong = outerTubeData.getSampledElementsCountAlong()
                # print("Resample startSurface", startSurface is not None, "endSurface", endSurface is not None)
                sx, sd1, sd2, sd12 = resampleTubeCoordinates(
                    rawTubeCoordinates, elementsCountAlong, startSurface=startSurface, endSurface=endSurface)
                tubeData.setSampledTubeCoordinates((sx, sd1, sd2, sd12))

        del segmentTubeData

        completedBifurcations = set()  # record so only done once

        with ChangeManager(fieldmodule):
            for networkSegment in networkSegments:
                segmentNodes = networkSegment.getNetworkNodes()
                startSegmentNode = segmentNodes[0]
                startInSegments = startSegmentNode.getInSegments()
                startOutSegments = startSegmentNode.getOutSegments()
                startSkipCount = 1 if ((len(startInSegments) > 1) or (len(startOutSegments) > 1)) else 0
                endSegmentNode = segmentNodes[-1]
                endInSegments = endSegmentNode.getInSegments()
                endOutSegments = endSegmentNode.getOutSegments()
                endSkipCount = 1 if ((len(endInSegments) > 1) or (len(endOutSegments) > 1)) else 0

                for stage in range(3):
                    if stage == 1:
                        # tube
                        outerTubeData = outerSegmentTubeData[networkSegment]
                        outerTubeCoordinates = outerTubeData.getSampledTubeCoordinates()
                        if bodyMeshGroup in outerTubeData.getAnnotationMeshGroups():
                            elementsAlongBodyDistalSegment = len(outerTubeCoordinates)
                            elementsAroundBodyDistalSegment = len(outerTubeCoordinates[0][0])
                        if cervixMeshGroup in outerTubeData.getAnnotationMeshGroups():
                            elementsAlongCervixSegment = len(outerTubeCoordinates)
                            elementsAroundCervixSegment = len(outerTubeCoordinates[0][0])
                        startNodeIds = outerTubeData.getStartNodeIds(startSkipCount)
                        endNodeIds = outerTubeData.getEndNodeIds(endSkipCount)
                        loop = (len(startInSegments) == 1) and (startInSegments[0] is networkSegment) and \
                               (networkSegment.getNodeVersions()[0] == networkSegment.getNodeVersions()[-1])
                        innerTubeData = innerSegmentTubeData[networkSegment] if layoutInnerCoordinates else None

                        if cervixMeshGroup in outerTubeData.getAnnotationMeshGroups():
                            elementStartIdxCervix = elementIdentifier

                        innerTubeCoordinates = innerTubeData.getSampledTubeCoordinates() if layoutInnerCoordinates else None
                        nodeIdentifier, elementIdentifier, startNodeIds, endNodeIds = generateTube(
                            outerTubeCoordinates, innerTubeCoordinates, elementsCountThroughWall,
                            region, fieldcache, coordinates, nodeIdentifier, elementIdentifier,
                            startSkipCount=startSkipCount, endSkipCount=endSkipCount,
                            startNodeIds=startNodeIds, endNodeIds=endNodeIds,
                            annotationMeshGroups=outerTubeData.getAnnotationMeshGroups(),
                            loop=loop, serendipity=serendipity)
                        outerTubeData.setStartNodeIds(startNodeIds, startSkipCount)
                        outerTubeData.setEndNodeIds(endNodeIds, endSkipCount)

                        if (len(startInSegments) == 1) and (startSkipCount == 0):
                            # copy startNodeIds to end of last segment
                            inTubeData = outerSegmentTubeData[startInSegments[0]]
                            inTubeData.setEndNodeIds(startNodeIds, 0)
                        if (len(endOutSegments) == 1) and (endSkipCount == 0):
                            # copy endNodesIds to start of next segment
                            outTubeData = outerSegmentTubeData[endOutSegments[0]]
                            outTubeData.setStartNodeIds(endNodeIds, 0)
                    else:
                        # start, end bifurcation
                        outerTubeBifurcationData = outerNodeTubeBifurcationData.get(
                            startSegmentNode if (stage == 0) else endSegmentNode)
                        if outerTubeBifurcationData and not outerTubeBifurcationData in completedBifurcations:
                            if showIntersectionCurves:
                                lineIdentifier = None
                                for s in range(3):
                                    curve = outerTubeBifurcationData.getIntersectionCurve(s)
                                    cx, cd1, cProportions, loop = curve
                                    if cx:
                                        nodeIdentifier, lineIdentifier = \
                                            generateCurveMesh(region, cx, cd1, loop, nodeIdentifier, lineIdentifier)
                            if showTrimSurfaces:
                                faceIdentifier = elementIdentifier if (dimension == 2) else None
                                for s in range(3):
                                    trimSurface = outerTubeBifurcationData.getTrimSurface(s)
                                    if trimSurface:
                                        nodeIdentifier, faceIdentifier = \
                                            trimSurface.generateMesh(region, nodeIdentifier, faceIdentifier)
                                if dimension == 2:
                                    elementIdentifier = faceIdentifier
                            innerTubeBifurcationData = None
                            if innerNodeTubeBifurcationData:
                                innerTubeBifurcationData = innerNodeTubeBifurcationData.get(
                                    startSegmentNode if (stage == 0) else endSegmentNode)

                            crossIndexes = outerTubeBifurcationData.getCrossIndexes()  # only get these from outer
                            if not crossIndexes:
                                outerTubeBifurcationData.determineCrossIndexes()
                                outerTubeBifurcationData.determineMidCoordinates()
                                if innerTubeBifurcationData:
                                    innerTubeBifurcationData.determineCrossIndexes()
                                    innerTubeBifurcationData.determineMidCoordinates()
                                crossIndexes = outerTubeBifurcationData.getCrossIndexes()

                            outerTubeCoordinates = outerTubeBifurcationData.getConnectingTubeCoordinates()
                            outerMidCoordinates = outerTubeBifurcationData.getMidCoordinates()
                            inward = outerTubeBifurcationData.getSegmentsIn()
                            outerTubeData = outerTubeBifurcationData.getTubeData()
                            tubeNodeIds = [outerTubeData[s].getEndNodeIds(1) if inward[s] else \
                                               outerTubeData[s].getStartNodeIds(1) for s in range(3)]
                            innerTubeCoordinates = None
                            innerMidCoordinates = None
                            if innerTubeBifurcationData:
                                innerTubeCoordinates = innerTubeBifurcationData.getConnectingTubeCoordinates()
                                innerMidCoordinates = innerTubeBifurcationData.getMidCoordinates()
                            annotationMeshGroups = [outerTubeData[s].getAnnotationMeshGroups() for s in range(3)]

                            nodeIdentifier, elementIdentifier = generateTubeBifurcation(
                                outerTubeCoordinates, innerTubeCoordinates, inward, elementsCountThroughWall,
                                outerMidCoordinates, innerMidCoordinates, crossIndexes,
                                region, fieldcache, coordinates, nodeIdentifier, elementIdentifier, tubeNodeIds,
                                annotationMeshGroups, serendipity=serendipity)

                            for s in range(3):
                                if inward[s]:
                                    if not outerTubeData[s].getEndNodeIds(1):
                                        outerTubeData[s].setEndNodeIds(tubeNodeIds[s], 1)
                                else:
                                    if not outerTubeData[s].getStartNodeIds(1):
                                        outerTubeData[s].setStartNodeIds(tubeNodeIds[s], 1)

                            completedBifurcations.add(outerTubeBifurcationData)

        is_uterus = uterusGroup.getGroup()
        myometriumGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_uterus_term("myometrium"))
        myometriumGroup.getMeshGroup(mesh).addElementsConditional(is_uterus)

        # Add human specific annotations
        if isHuman:
            # find elements on left and right edges
            elementsOnRightMarginVentral = []
            elementsOnRightMarginDorsal = []
            elementsOnLeftMarginVentral = []
            elementsOnLeftMarginDorsal = []

            for i in range(elementsAlongBodyDistalSegment):
                elementsOnRightMarginVentral.append(elementStartIdxCervix - 1 - i * elementsAroundBodyDistalSegment - (math.ceil(elementsAroundBodyDistalSegment * 0.25) - 1))
                elementsOnRightMarginDorsal.append(elementStartIdxCervix - 1 - i * elementsAroundBodyDistalSegment - (math.ceil(elementsAroundBodyDistalSegment * 0.25) - 1) - 1)
                elementsOnLeftMarginVentral.append(elementStartIdxCervix - (i + 1) * elementsAroundBodyDistalSegment + (math.ceil(elementsAroundBodyDistalSegment * 0.25) - 1))
                elementsOnLeftMarginDorsal.append(elementStartIdxCervix - (i + 1) * elementsAroundBodyDistalSegment + (math.ceil(elementsAroundBodyDistalSegment * 0.25) - 1) + 1)

            nearRightMarginDorsalGroup = AnnotationGroup(region, ("elements adjacent to right margin dorsal", "None"))
            nearRightMarginVentralGroup = AnnotationGroup(region, ("elements adjacent to right margin ventral", "None"))
            nearLeftMarginDorsalGroup = AnnotationGroup(region, ("elements adjacent to left margin dorsal", "None"))
            nearLeftMarginVentralGroup = AnnotationGroup(region, ("elements adjacent to left margin ventral", "None"))

            # track backwards to find elements near the margins
            lastElements = [elementsOnLeftMarginVentral[-1], elementsOnLeftMarginDorsal[-1],
                            elementsOnRightMarginVentral[-1], elementsOnRightMarginDorsal[-1]]
            marginGroups = [nearLeftMarginVentralGroup, nearLeftMarginDorsalGroup,
                            nearRightMarginVentralGroup, nearRightMarginDorsalGroup]

            for i in range(4):
                adjacentElement = mesh.findElementByIdentifier(lastElements[i])
                eft = adjacentElement.getElementfieldtemplate(coordinates, -1)
                lastNode = get_element_node_identifiers(adjacentElement, eft)[4]
                addElementsToMargin(mesh, coordinates, marginGroups[i], lastNode)

            elementIter = mesh.createElementiterator()
            element = elementIter.next()
            while element.isValid():
                elementIdentifier = element.getIdentifier()
                if elementIdentifier in elementsOnRightMarginVentral:
                    nearRightMarginVentralGroup.getMeshGroup(mesh).addElement(element)
                elif elementIdentifier in elementsOnRightMarginDorsal:
                    nearRightMarginDorsalGroup.getMeshGroup(mesh).addElement(element)
                elif elementIdentifier in elementsOnLeftMarginVentral:
                    nearLeftMarginVentralGroup.getMeshGroup(mesh).addElement(element)
                elif elementIdentifier in elementsOnLeftMarginDorsal:
                    nearLeftMarginDorsalGroup.getMeshGroup(mesh).addElement(element)
                element = elementIter.next()

            annotationGroups.append(nearRightMarginDorsalGroup)
            annotationGroups.append(nearRightMarginVentralGroup)
            annotationGroups.append(nearLeftMarginDorsalGroup)
            annotationGroups.append(nearLeftMarginVentralGroup)

            elementsOnRightMarginVentral = []
            elementsOnRightMarginDorsal = []
            elementsOnLeftMarginVentral = []
            elementsOnLeftMarginDorsal = []

            for i in range(elementsAlongCervixSegment):
                elementsOnRightMarginDorsal.append(elementStartIdxCervix - 1 + i * elementsAroundCervixSegment - (math.ceil(elementsAroundBodyDistalSegment * 0.25) - 1) - 1)
                elementsOnRightMarginVentral.append(elementStartIdxCervix - 1 + i * elementsAroundCervixSegment - (math.ceil(elementsAroundBodyDistalSegment * 0.25) - 1))
                elementsOnLeftMarginVentral.append(elementStartIdxCervix + (math.ceil(elementsAroundBodyDistalSegment * 0.25) - 1) + i * elementsAroundCervixSegment)
                elementsOnLeftMarginDorsal.append(elementStartIdxCervix + (math.ceil(elementsAroundBodyDistalSegment * 0.25) - 1) + 1 + i * elementsAroundCervixSegment)

            nearRightMarginDorsalCervixGroup = AnnotationGroup(region, ("elements adjacent to right margin dorsal of cervix", "None"))
            nearRightMarginVentralCervixGroup = AnnotationGroup(region, ("elements adjacent to right margin ventral of cervix", "None"))
            nearLeftMarginDorsalCervixGroup = AnnotationGroup(region, ("elements adjacent to left margin dorsal of cervix", "None"))
            nearLeftMarginVentralCervixGroup = AnnotationGroup(region, ("elements adjacent to left margin ventral of cervix", "None"))

            elementIter = mesh.createElementiterator()
            element = elementIter.next()
            while element.isValid():
                elementIdentifier = element.getIdentifier()
                if elementIdentifier in elementsOnRightMarginVentral:
                    nearRightMarginVentralCervixGroup.getMeshGroup(mesh).addElement(element)
                elif elementIdentifier in elementsOnRightMarginDorsal:
                    nearRightMarginDorsalCervixGroup.getMeshGroup(mesh).addElement(element)
                elif elementIdentifier in elementsOnLeftMarginVentral:
                    nearLeftMarginVentralCervixGroup.getMeshGroup(mesh).addElement(element)
                elif elementIdentifier in elementsOnLeftMarginDorsal:
                    nearLeftMarginDorsalCervixGroup.getMeshGroup(mesh).addElement(element)
                element = elementIter.next()

            annotationGroups.append(nearRightMarginDorsalCervixGroup)
            annotationGroups.append(nearRightMarginVentralCervixGroup)
            annotationGroups.append(nearLeftMarginDorsalCervixGroup)
            annotationGroups.append(nearLeftMarginVentralCervixGroup)

            allMarkers = {"left round ligament of uterus": {"x": [-4.13264, 9.97541, -7.1167]},
                          "right round ligament of uterus": {"x": [4.13279, 9.97534, -7.1166]}}

            for key in allMarkers:
                x = allMarkers[key]["x"]
                group = findOrCreateAnnotationGroupForTerm( annotationGroups, region,
                                                            get_uterus_term(key), isMarker=True)
                markerNode = group.createMarkerNode(nodeIdentifier, coordinates, x)
                nodeIdentifier = markerNode.getIdentifier() + 1
                for group in annotationGroups:
                    group.getNodesetGroup(nodes).addNode(markerNode)

        preBifurcationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 ("pre-bifurcation segments", "None"))
        postBifurcationGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 ("post-bifurcation segments", "None"))
        annotationGroups.remove(preBifurcationGroup)
        annotationGroups.remove(postBifurcationGroup)

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
        parameterSetName = options['Base parameter set']
        isHuman = parameterSetName in ("Default", "Human 1")
        isMouse = parameterSetName in "Mouse 1"

        # Create 2d surface mesh groups
        fm = region.getFieldmodule()
        uterusGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("uterus"))
        cervixGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("uterine cervix"))
        bodyGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("body of uterus"))
        vaginaGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("vagina"))

        mesh1d = fm.findMeshByDimension(1)
        mesh2d = fm.findMeshByDimension(2)

        is_exterior = fm.createFieldIsExterior()
        is_exterior_face_xi3_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_1))
        is_exterior_face_xi3_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI3_0))
        is_exterior_face_xi2_0 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_0))
        is_exterior_face_xi2_1 = fm.createFieldAnd(is_exterior, fm.createFieldIsOnFace(Element.FACE_TYPE_XI2_1))

        is_uterus = uterusGroup.getGroup()
        is_uterus_outer = fm.createFieldAnd(is_uterus, is_exterior_face_xi3_1)
        is_uterus_inner = fm.createFieldAnd(is_uterus, is_exterior_face_xi3_0)

        is_body = bodyGroup.getGroup()
        is_body_outer = fm.createFieldAnd(is_body, is_exterior_face_xi3_1)
        is_body_inner = fm.createFieldAnd(is_body, is_exterior_face_xi3_0)

        is_cervix = cervixGroup.getGroup()
        is_cervix_outer = fm.createFieldAnd(is_cervix, is_exterior_face_xi3_1)
        is_cervix_inner = fm.createFieldAnd(is_cervix, is_exterior_face_xi3_0)

        is_vagina = vaginaGroup.getGroup()
        is_vagina_outer = fm.createFieldAnd(is_vagina, is_exterior_face_xi3_1)
        is_vagina_inner = fm.createFieldAnd(is_vagina, is_exterior_face_xi3_0)
        is_vagina_xi2_0 = fm.createFieldAnd(is_vagina, is_exterior_face_xi2_0)
        is_vagina_xi2_1 = fm.createFieldAnd(is_vagina, is_exterior_face_xi2_1)
        is_vagina_xi2_01 = fm.createFieldXor(is_vagina_xi2_0, is_vagina_xi2_1)

        serosaOfUterus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of uterus"))
        serosaOfUterus.getMeshGroup(mesh2d).addElementsConditional(is_uterus_outer)

        lumenOfUterus = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of uterus"))
        lumenOfUterus.getMeshGroup(mesh2d).addElementsConditional(is_uterus_inner)

        serosaOfBody = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                          get_uterus_term("serosa of body of uterus"))
        serosaOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_outer)

        lumenOfBody = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                         get_uterus_term("lumen of body of uterus"))
        lumenOfBody.getMeshGroup(mesh2d).addElementsConditional(is_body_inner)

        serosaOfCervix = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of uterine cervix"))
        serosaOfCervix.getMeshGroup(mesh2d).addElementsConditional(is_cervix_outer)

        lumenOfCervix = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of uterine cervix"))
        lumenOfCervix.getMeshGroup(mesh2d).addElementsConditional(is_cervix_inner)

        serosaOfVagina = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("serosa of vagina"))
        serosaOfVagina.getMeshGroup(mesh2d).addElementsConditional(is_vagina_outer)

        lumenOfVagina = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                           get_uterus_term("lumen of vagina"))
        lumenOfVagina.getMeshGroup(mesh2d).addElementsConditional(is_vagina_inner)

        if isHuman:
            fallopianTubeGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("fallopian tube"))
            is_fallopianTube = fallopianTubeGroup.getGroup()
            is_fallopianTube_outer = fm.createFieldAnd(is_fallopianTube, is_exterior_face_xi3_1)
            is_fallopianTube_inner = fm.createFieldAnd(is_fallopianTube, is_exterior_face_xi3_0)

            serosaOfFallopianTube = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_uterus_term("serosa of fallopian tube"))
            serosaOfFallopianTube.getMeshGroup(mesh2d).addElementsConditional(is_fallopianTube_outer)

            lumenOfFallopianTube = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("lumen of fallopian tube"))
            lumenOfFallopianTube.getMeshGroup(mesh2d).addElementsConditional(is_fallopianTube_inner)

            nearRightMarginDorsalGroup = getAnnotationGroupForTerm(annotationGroups,
                                                                   ("elements adjacent to right margin dorsal",
                                                                    "None"))
            nearRightMarginVentralGroup = getAnnotationGroupForTerm(annotationGroups,
                                                                    ("elements adjacent to right margin ventral",
                                                                     "None"))
            nearLeftMarginDorsalGroup = getAnnotationGroupForTerm(annotationGroups,
                                                                  ("elements adjacent to left margin dorsal",
                                                                   "None"))
            nearLeftMarginVentralGroup = getAnnotationGroupForTerm(annotationGroups,
                                                                   ("elements adjacent to left margin ventral",
                                                                    "None"))

            is_pubocervical = fm.createFieldAnd(is_body_outer, is_cervix_outer)
            pubocervical = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("pubocervical ligament (TA98)"))
            pubocervical.getMeshGroup(mesh1d).addElementsConditional(is_pubocervical)

            is_internal_os = fm.createFieldAnd(is_body_inner, is_cervix_inner)
            internalOs = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_uterus_term("internal cervical os"))
            internalOs.getMeshGroup(mesh1d).addElementsConditional(is_internal_os)

            is_external_os = fm.createFieldAnd(is_vagina_inner, is_cervix_inner)
            externalOs = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("external cervical os"))
            externalOs.getMeshGroup(mesh1d).addElementsConditional(is_external_os)

            is_vagina_orifice = fm.createFieldAnd(is_vagina_xi2_01, is_exterior_face_xi3_0)
            vaginaOrifice = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                            get_uterus_term("vagina orifice"))
            vaginaOrifice.getMeshGroup(mesh1d).addElementsConditional(is_vagina_orifice)

            # Broad ligament of uterus
            is_nearLeftMarginDorsal = fm.createFieldAnd(nearLeftMarginDorsalGroup.getGroup(), is_exterior_face_xi3_1)
            is_nearLeftMarginVentral = fm.createFieldAnd(nearLeftMarginVentralGroup.getGroup(), is_exterior_face_xi3_1)
            is_leftBroadLigament = fm.createFieldAnd(is_nearLeftMarginDorsal, is_nearLeftMarginVentral)
            leftBroadLigament = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                              get_uterus_term("left broad ligament of uterus"))
            leftBroadLigament.getMeshGroup(mesh1d).addElementsConditional(is_leftBroadLigament)

            is_nearRightMarginDorsal = fm.createFieldAnd(nearRightMarginDorsalGroup.getGroup(), is_exterior_face_xi3_1)
            is_nearRightMarginVentral = fm.createFieldAnd(nearRightMarginVentralGroup.getGroup(), is_exterior_face_xi3_1)
            is_rightBroadLigament = fm.createFieldAnd(is_nearRightMarginDorsal, is_nearRightMarginVentral)
            rightBroadLigament = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_uterus_term("right broad ligament of uterus"))
            rightBroadLigament.getMeshGroup(mesh1d).addElementsConditional(is_rightBroadLigament)

            # Transverse cervical ligament
            nearRightMarginDorsalCervixGroup = \
                getAnnotationGroupForTerm(annotationGroups,
                                          ("elements adjacent to right margin dorsal of cervix", "None"))
            nearRightMarginVentralCervixGroup = \
                getAnnotationGroupForTerm(annotationGroups,
                                          ("elements adjacent to right margin ventral of cervix", "None"))
            nearLeftMarginDorsalCervixGroup = \
                getAnnotationGroupForTerm(annotationGroups,
                                          ("elements adjacent to left margin dorsal of cervix", "None"))
            nearLeftMarginVentralCervixGroup = \
                getAnnotationGroupForTerm(annotationGroups,
                                          ("elements adjacent to left margin ventral of cervix", "None"))

            is_nearLeftMarginDorsalCervix = fm.createFieldAnd(nearLeftMarginDorsalCervixGroup.getGroup(),
                                                              is_cervix_outer)
            is_nearLeftMarginVentralCervix = fm.createFieldAnd(nearLeftMarginVentralCervixGroup.getGroup(),
                                                               is_cervix_outer)
            is_leftTransverseCervicalLigament = fm.createFieldAnd(is_nearLeftMarginDorsalCervix,
                                                                  is_nearLeftMarginVentralCervix)
            leftTransverseCervicalLigament = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_uterus_term("left transverse cervical ligament"))
            leftTransverseCervicalLigament.getMeshGroup(mesh1d).addElementsConditional(is_leftTransverseCervicalLigament)

            is_nearRightMarginDorsalCervix = fm.createFieldAnd(nearRightMarginDorsalCervixGroup.getGroup(),
                                                              is_cervix_outer)
            is_nearRightMarginVentralCervix = fm.createFieldAnd(nearRightMarginVentralCervixGroup.getGroup(),
                                                               is_cervix_outer)
            is_rightTransverseCervicalLigament = fm.createFieldAnd(is_nearRightMarginDorsalCervix,
                                                                  is_nearRightMarginVentralCervix)
            rightTransverseCervicalLigament = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_uterus_term("right transverse cervical ligament"))
            rightTransverseCervicalLigament.getMeshGroup(mesh1d).addElementsConditional(
                is_rightTransverseCervicalLigament)

            annotationGroups.remove(nearRightMarginDorsalGroup)
            annotationGroups.remove(nearRightMarginVentralGroup)
            annotationGroups.remove(nearLeftMarginDorsalGroup)
            annotationGroups.remove(nearLeftMarginVentralGroup)

            annotationGroups.remove(nearRightMarginDorsalCervixGroup)
            annotationGroups.remove(nearRightMarginVentralCervixGroup)
            annotationGroups.remove(nearLeftMarginDorsalCervixGroup)
            annotationGroups.remove(nearLeftMarginVentralCervixGroup)

        if isMouse:
            rightHornGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("right uterine horn"))
            leftHornGroup = getAnnotationGroupForTerm(annotationGroups, get_uterus_term("left uterine horn"))

            is_rightHorn = rightHornGroup.getGroup()
            is_rightHorn_outer = fm.createFieldAnd(is_rightHorn, is_exterior_face_xi3_1)
            is_rightHorn_inner = fm.createFieldAnd(is_rightHorn, is_exterior_face_xi3_0)

            is_leftHorn = leftHornGroup.getGroup()
            is_leftHorn_outer = fm.createFieldAnd(is_leftHorn, is_exterior_face_xi3_1)
            is_leftHorn_inner = fm.createFieldAnd(is_leftHorn, is_exterior_face_xi3_0)

            serosaOfRightHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                   get_uterus_term("serosa of right horn"))
            serosaOfRightHorn.getMeshGroup(mesh2d).addElementsConditional(is_rightHorn_outer)

            lumenOfRightHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("lumen of right horn"))
            lumenOfRightHorn.getMeshGroup(mesh2d).addElementsConditional(is_rightHorn_inner)

            serosaOfLeftHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                  get_uterus_term("serosa of left horn"))
            serosaOfLeftHorn.getMeshGroup(mesh2d).addElementsConditional(is_leftHorn_outer)

            lumenOfLeftHorn = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_uterus_term("lumen of left horn"))
            lumenOfLeftHorn.getMeshGroup(mesh2d).addElementsConditional(is_leftHorn_inner)


def addElementsToMargin(mesh, coordinates, group, nodeId):
    """
    Track from the node to find elements Add elements lying on the margin and add them to a group.
    :param coordinates: Coordinates field.
    :param group: Group for adding elements lying on the margin.
    :param nodeId: Node to start tracking from.
    """
    elementIter = mesh.createElementiterator()
    element = elementIter.next()
    while element.isValid():
        eft = element.getElementfieldtemplate(coordinates, -1)
        nodeIdentifiers = get_element_node_identifiers(element, eft)
        if nodeIdentifiers[6] == nodeId:
            group.getMeshGroup(mesh).addElement(element)
            adjacentElement = element
            eft = adjacentElement.getElementfieldtemplate(coordinates, -1)
            nodeId = get_element_node_identifiers(adjacentElement, eft)[4]
            elementIter = mesh.createElementiterator()
        element = elementIter.next()

    return
