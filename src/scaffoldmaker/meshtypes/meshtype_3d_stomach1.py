"""
Generates a 3-D stomach mesh along the central line, with variable
numbers of elements around esophagus and duodenum, along and through
wall, with variable radius and thickness along.
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
from scaffoldmaker.annotation.stomach_terms import get_stomach_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.bifurcation import get_bifurcation_triple_point
from scaffoldmaker.utils.eft_utils import scaleEftNodeValueLabels, setEftScaleFactorIds, remapEftNodeValueLabel
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import scaleEftNodeValueLabels, setEftScaleFactorIds, remapEftNodeValueLabel, \
    remapEftNodeValueLabelsVersion
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues, mesh_destroy_elements_and_nodes_by_identifiers


class MeshType_3d_stomach1(Scaffold_base):
    """
    Generates a 3-D stomach mesh with variable numbers of elements around the esophagus and duodenum,
    along the central line, and through wall. The stomach is created using a central path as the longitudinal axis
    of the stomach. D2 of the central path points to the greater curvature of the stomach and magnitude of D2 and D3
    are the radii of the stomach in the respective direction. D2 on the first node of the central path provide
    the radius of the fundus dome.
    """
    centralPathDefaultScaffoldPackages = {
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
                [ [  70.8, 72.3, 0.0 ], [ -17.0, -47.9, 0.0 ], [  39.3, -13.9, 0.0 ], [  -3.7, -5.9, 0.0 ], [ 0.0, 0.0, 38.7 ], [ 0.0, 0.0,  4.5] ],
                [ [  51.1, 30.7, 0.0 ], [ -22.3, -34.9, 0.0 ], [  32.2, -20.5, 0.0 ], [ -10.6, -7.2, 0.0 ], [ 0.0, 0.0, 41.4 ], [ 0.0, 0.0,  0.9] ],
                [ [  27.8,  2.7, 0.0 ], [ -27.7, -18.4, 0.0 ], [  18.7, -28.2, 0.0 ], [ -15.1, -2.0, 0.0 ], [ 0.0, 0.0, 40.9 ], [ 0.0, 0.0, -5.0] ],
                [ [  -0.3, -5.8, 0.0 ], [ -28.0,  -2.6, 0.0 ], [   2.3, -25.6, 0.0 ], [ -12.9,  4.6, 0.0 ], [ 0.0, 0.0, 32.3 ], [ 0.0, 0.0, -9.4] ],
                [ [ -26.5, -2.9, 0.0 ], [ -20.6,   8.0, 0.0 ], [  -7.5, -19.2, 0.0 ], [  -6.4,  8.0, 0.0 ], [ 0.0, 0.0, 22.1 ], [ 0.0, 0.0, -6.8] ],
                [ [ -40.6,  7.4, 0.0 ], [ -11.0,  12.5, 0.0 ], [ -11.6, -10.2, 0.0 ], [  -0.4,  7.8, 0.0 ], [ 0.0, 0.0, 17.6 ], [ 0.0, 0.0, -5.9] ],
                [ [ -48.1, 21.0, 0.0 ], [  -6.4,  15.7, 0.0 ], [  -8.6,  -3.5, 0.0 ], [   0.4,  4.4, 0.0 ], [ 0.0, 0.0, 10.5 ], [ 0.0, 0.0, -3.1] ],
                [ [ -52.9, 38.6, 0.0 ], [  -3.2,  19.3, 0.0 ], [ -11.1,  -1.9, 0.0 ], [  -5.5, -1.1, 0.0 ], [ 0.0, 0.0, 12.0 ], [ 0.0, 0.0,  6.1] ] ] ),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '2-3',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4-5',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_stomach_term('duodenum')[0],
                    'ontId': get_stomach_term('duodenum')[1]
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
                    'identifierRanges': '1-3',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4-5',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '8',
                    'name': get_stomach_term('duodenum')[0],
                    'ontId': get_stomach_term('duodenum')[1]
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
                    'identifierRanges': '1',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '2-4',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_stomach_term('duodenum')[0],
                    'ontId': get_stomach_term('duodenum')[1]
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
                    'identifierRanges': '1-3',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '4-5',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '8',
                    'name': get_stomach_term('duodenum')[0],
                    'ontId': get_stomach_term('duodenum')[1]
                }]
        }),
        'Stomach 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'D3 derivatives': True,
                'Length': 1.0,
                'Number of elements': 7
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                [ [2.0,0.0,0.0], [-0.6,0.0,0.0], [0.0,-0.5,0.0], [0.0,0.0,0.0], [0.0,0.0,0.5], [0.0,0.0,0.0] ] , 
                [ [1.5,0.0,0.0], [-0.4,0.0,0.0], [0.0,-0.5,0.0], [0.0,0.0,0.0], [0.0,0.0,0.5], [0.0,0.0,0.0] ] , 
                [ [1.2,0.0,0.0], [-0.2,0.0,0.0], [0.0,-0.5,0.0], [0.0,0.0,0.0], [0.0,0.0,0.5], [0.0,0.0,0.0] ] , 
                [ [1.0,0.0,0.0], [-0.3,0.0,0.0], [0.0,-0.5,0.0], [0.0,0.0,0.0], [0.0,0.0,0.5], [0.0,0.0,0.0] ] , 
                [ [0.6,0.0,0.0], [-0.3,0.0,0.0], [0.0,-0.5,0.0], [0.0,0.07,0.0], [0.0,0.0,0.5], [0.0,0.0,-0.1] ],
                [ [0.4,0.0,0.0], [-0.2,0.0,0.0], [0.0,-0.4,0.0], [0.0,0.17,0.0], [0.0,0.0,0.4], [0.0,0.0,-0.17] ],
                [ [0.2,0.0,0.0], [-0.2,0.0,0.0], [0.0,-0.15,0.0], [0.0,0.1,0.0], [0.0,0.0,0.15], [0.0,0.0,-0.1] ],
                [ [0.0,0.0,0.0], [-0.2,0.0,0.0], [0.0,-0.2,0.0], [0.0,-0.2,0.0], [0.0,0.0,0.2], [0.0,0.0,0.2] ]]),
                     
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_stomach_term('fundus of stomach')[0],
                    'ontId': get_stomach_term('fundus of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-4',
                    'name': get_stomach_term('body of stomach')[0],
                    'ontId': get_stomach_term('body of stomach')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5',
                    'name': get_stomach_term('pyloric antrum')[0],
                    'ontId': get_stomach_term('pyloric antrum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '6',
                    'name': get_stomach_term('pyloric canal')[0],
                    'ontId': get_stomach_term('pyloric canal')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '7',
                    'name': get_stomach_term('duodenum')[0],
                    'ontId': get_stomach_term('duodenum')[1]
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
                'Unit scale': 1.0,
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
        'Mouse 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 12,
                'Number of elements along': 2,
                'Number of elements through wall': 4,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 1.5,
                'Ostium length': 1.5,
                'Ostium wall thickness': 0.45,
                'Ostium wall relative thicknesses': [0.75, 0.05, 0.15, 0.05],
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.5,
                'Vessel wall thickness': 0.27,
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
                'Unit scale': 1.0,
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
                'Number of elements around ostium': 12,
                'Number of elements along': 2,
                'Number of elements through wall': 4,
                'Unit scale': 1.0,
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
        'Stomach 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 4,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 0.3, # changed
                'Ostium length': 0.3, # changed
                'Ostium wall thickness': 0.05, # changed
                'Ostium wall relative thicknesses': [0.25, 0.25, 0.25, 0.25], # changed
                'Ostium inter-vessel distance': 0.0,
                'Ostium inter-vessel height': 0.0,
                'Use linear through ostium wall': True,
                'Vessel end length factor': 1.0,
                'Vessel inner diameter': 0.1, # changed
                'Vessel wall thickness': 0.03, # changed
                'Vessel wall relative thicknesses': [0.25, 0.25, 0.25, 0.25], # changed
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
            'Mouse 1',
            'Pig 1',
            'Rat 1',
            'Stomach 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Mouse 1']
        elif 'Pig 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Pig 1']
        elif 'Rat 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Rat 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Rat 1']
        elif 'Stomach 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Stomach 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Stomach 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Human 1']

        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of elements around esophagus': 8,
            'Number of elements around duodenum': 12,
            'Number of elements between fundus apex and cardia': 3,
            'Number of elements between cardia and duodenum': 5,
            'Number of elements across cardia': 1,
            'Number of elements through wall': 4,
            'Wall thickness': 5.0,
            'Mucosa relative thickness': 0.55,
            'Submucosa relative thickness': 0.15,
            'Circular muscle layer relative thickness': 0.25,
            'Longitudinal muscle layer relative thickness': 0.05,
            'Limiting ridge': False,
            'Gastro-esophagal junction': copy.deepcopy(ostiumOption),
            # 'Gastro-esophagal junction position along factor': 0.35,
            'Cardia diameter factor': 1.5,
            'Use linear through wall': True,
            'Track surface': False, # KM
            'Refine': False,
            'Refine number of elements surface': 4,
            'Refine number of elements through wall': 1
        }
        if 'Mouse 1' in parameterSetName:
            options['Number of elements around esophagus'] = 8
            options['Number of elements around duodenum'] = 12
            options['Number of elements between fundus apex and cardia'] = 4
            options['Number of elements between cardia and duodenum'] = 4
            options['Wall thickness'] = 0.45
            options['Mucosa relative thickness'] = 0.75
            options['Submucosa relative thickness'] = 0.05
            options['Circular muscle layer relative thickness'] = 0.15
            options['Longitudinal muscle layer relative thickness'] = 0.05
            # options['Gastro-esophagal junction position along factor'] = 0.53
            options['Limiting ridge'] = True
        elif 'Pig 1' in parameterSetName:
            options['Number of elements around esophagus'] = 8
            options['Number of elements around duodenum'] = 12
            options['Number of elements between fundus apex and cardia'] = 3
            options['Number of elements between cardia and duodenum'] = 5
            options['Wall thickness'] = 5.0
            options['Mucosa relative thickness'] = 0.47
            options['Submucosa relative thickness'] = 0.1
            options['Circular muscle layer relative thickness'] = 0.33
            options['Longitudinal muscle layer relative thickness'] = 0.1
            # options['Gastro-esophagal junction position along factor'] = 0.45
            options['Limiting ridge'] = False
        elif 'Rat 1' in parameterSetName:
            options['Number of elements around esophagus'] = 8
            options['Number of elements around duodenum'] = 12
            options['Number of elements between fundus apex and cardia'] = 4
            options['Number of elements between cardia and duodenum'] = 4
            options['Wall thickness'] = 0.5
            options['Mucosa relative thickness'] = 0.65
            options['Submucosa relative thickness'] = 0.12
            options['Circular muscle layer relative thickness'] = 0.18
            options['Longitudinal muscle layer relative thickness'] = 0.05
            # options['Gastro-esophagal junction position along factor'] = 0.55
            options['Limiting ridge'] = True
        elif 'Stomach 1' in parameterSetName:
            options['Number of elements around esophagus'] = 8
            options['Number of elements around duodenum'] = 12
            options['Number of elements between fundus apex and cardia'] = 4
            options['Number of elements between cardia and duodenum'] = 4
            options['Wall thickness'] = 0.05
            options['Mucosa relative thickness'] = 0.25
            options['Submucosa relative thickness'] = 0.25
            options['Circular muscle layer relative thickness'] = 0.25
            options['Longitudinal muscle layer relative thickness'] = 0.25
            # options['Gastro-esophagal junction position along factor'] = 0.4
            options['Limiting ridge'] = False

        cls.updateSubScaffoldOptions(options)
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of elements around esophagus',
            'Number of elements around duodenum',
            'Number of elements between fundus apex and cardia',
            'Number of elements between cardia and duodenum',
            'Number of elements across cardia',
            'Number of elements through wall',
            'Wall thickness',
            'Mucosa relative thickness',
            'Submucosa relative thickness',
            'Circular muscle layer relative thickness',
            'Longitudinal muscle layer relative thickness',
            'Limiting ridge',
            'Gastro-esophagal junction',
            # 'Gastro-esophagal junction position along factor',
            'Cardia diameter factor',
            'Use linear through wall',
            'Track surface', # KM
            'Refine',
            'Refine number of elements surface',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [MeshType_1d_path1]
        if optionName == 'Gastro-esophagal junction':
            return [MeshType_3d_ostium1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
        if optionName == 'Gastro-esophagal junction':
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
        if optionName == 'Gastro-esophagal junction':
            if not parameterSetName:
                parameterSetName = list(cls.ostiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.ostiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_path1)
        if not options['Gastro-esophagal junction'].getScaffoldType() in cls.getOptionValidScaffoldTypes(
                'Gastro-esophagal junction'):
            options['Gastro-esophagal junction'] = cls.getOptionScaffoldPackage('Gastro-esophagal junction',
                                                                                MeshType_3d_ostium1)
        # Add condition for elements between cardia and duodenum wrt to number of annotations groups
        if options['Number of elements around esophagus'] < 8:
            options['Number of elements around esophagus'] = 8
        if options['Number of elements around esophagus'] % 4 > 0:
            options['Number of elements around esophagus'] = options['Number of elements around esophagus'] // 4 * 4
        if options['Number of elements around duodenum'] < 12:
            options['Number of elements around duodenum'] = 12
        if options['Number of elements around duodenum'] % 4 > 0:
            options['Number of elements around duodenum'] = options['Number of elements around duodenum'] // 4 * 4
        if options['Number of elements between fundus apex and cardia'] < 3:
            options['Number of elements between fundus apex and cardia'] = 3
        if options['Number of elements between fundus apex and cardia'] < int(options['Number of elements around duodenum'] * 0.25):
            options['Number of elements between fundus apex and cardia'] = int(options['Number of elements around duodenum'] * 0.25)
        if options['Number of elements between cardia and duodenum'] < 3:
            options['Number of elements between cardia and duodenum'] = 3
        if options['Number of elements through wall'] != (1 or 4):
            options['Number of elements through wall'] = 4
        if options['Cardia diameter factor'] <= 1.5:
            options['Cardia diameter factor'] = 1.5
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
        wallThickness = options['Wall thickness']
        ostiumOptions = options['Gastro-esophagal junction']
        ostiumSettings = ostiumOptions.getScaffoldSettings()
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
        :return: None
        """
        cls.updateSubScaffoldOptions(options)
        geometricCentralPath = options['Central path']
        materialCentralPath = cls.centralPathDefaultScaffoldPackages['Stomach 1']
        elementsCountAlongTrackSurface = 20
        allAnnotationGroups = []

        stomachTermsAlong = [None, 'fundus of stomach', 'body of stomach',
                             'pyloric antrum', 'pyloric canal', 'duodenum']

        # Geometric coordinates
        fm = region.getFieldmodule()
        coordinates = findOrCreateFieldCoordinates(fm)

        geometricCentralPath = StomachCentralPath(region, geometricCentralPath, stomachTermsAlong)

        elementCountGroupList = []
        allAnnotationGroups, elementCountGroupList = \
            createStomachMesh3d(region, fm, coordinates, stomachTermsAlong, elementsCountAlongTrackSurface,
                                allAnnotationGroups, elementCountGroupList, centralPath=geometricCentralPath,
                                options=options, nodeIdentifier=1, elementIdentifier=1, splitCoordinates=True,
                                materialCoordinates=False)

        # Material coordinates
        print('Making Stomach Coordinates')
        tmp_region = region.createRegion()
        tmp_fm = tmp_region.getFieldmodule()
        tmp_stomach_coordinates = findOrCreateFieldCoordinates(tmp_fm, name="stomach coordinates")

        materialCentralPath = StomachCentralPath(tmp_region, materialCentralPath, stomachTermsAlong)

        allAnnotationGroups, elementCountGroupList = \
            createStomachMesh3d(tmp_region, tmp_fm, tmp_stomach_coordinates, stomachTermsAlong,
                                elementsCountAlongTrackSurface, allAnnotationGroups, elementCountGroupList,
                                centralPath=materialCentralPath, options=options, nodeIdentifier=1, elementIdentifier=1,
                                splitCoordinates=False, materialCoordinates=True)

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
        del tmp_fm
        del tmp_stomach_coordinates
        del tmp_region

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
    #     duodenumGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("duodenum"))
    #     fundusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("fundus of stomach"))
    #     antrumGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pyloric antrum"))
    #     pylorusGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("pyloric canal"))
    #     esoGroup = getAnnotationGroupForTerm(annotationGroups, get_stomach_term("esophagus"))
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
    #                                                               get_stomach_term("luminal surface of duodenum"))
    #     duodenumSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                              get_stomach_term("serosa of duodenum"))
    #     esoLuminalGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                          get_stomach_term("luminal surface of esophagus"))
    #     esoSerosaGroup = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
    #                                                         get_stomach_term("serosa of esophagus"))
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
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
    #                 "circular-longitudinal muscle interface of first segment of the duodenum along the "
    #                 "gastric-omentum attachment"))
    #         esoCurvaturesCMLMGroup = \
    #             findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_stomach_term(
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
    def __init__(self, region, centralPath, # elementsCountAlong,
                 stomachTermsAlong=[None]):
        """
        :param region: Zinc region to define model in.
        :param centralPath: Central path subscaffold from meshtype_1d_path1
        :param elementsCountAlong: Number of sampled elements to be returned in central path.
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

        for i in range(len(stomachTermsAlong)):
            tmpRegion = region.createRegion()
            centralPath.generate(tmpRegion)
            cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                            Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                            Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                                groupName=stomachTermsAlong[i])

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

            # if i == 0:
            #     cx = cxGroup
            #     cd1 = cd1Group
            #     cd2 = cd2Group
            #     cd3 = cd3Group
            #     cd12 = cd12Group
            #     cd13 = cd13Group
            del tmpRegion

        # print([vector.magnitude(cd2[n]) for n in range(len(cd2))])
        # self.sx, self.sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
        # self.sd2, self.sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
        # self.sd3, self.sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, se, sxi, ssf)
        self.arcLengthOfGroupsAlong = arcLengthOfGroupsAlong
        self.cxGroups = cxGroups
        self.cd1Groups = cd1Groups
        self.cd2Groups = cd2Groups
        self.cd3Groups = cd3Groups
        self.cd12Groups = cd12Groups
        self.cd13Groups = cd13Groups

        # print('Total arcLengths:', arcLengthOfGroupsAlong[0])
        # print('arcLengthOfGroupsAlong:', arcLengthOfGroupsAlong)


def findClosestPositionAndDerivativeOnTrackSurface(x, nx, trackSurface, nxProportion1, elementsCountAlongTrackSurface):
    """
    Find the closest position and derivative around the tracksurface of a point sitting near the fundus of stomach.
    Use a startPosition to improve the search for nearest position on the track surface as the fundus has a curved
    and complex track surface.
    :param x: coordinates of point of interest
    :param nx: coordinates of points along curve where point of interest lies.
    :param trackSurface: track surface where point sits
    :param nxProportion1: proportion around track surface of curve
    :param elementsCountAlongTrackSurface: number of elements along track surface
    :return: position and derivative of point around track surface
    """
    closestIdxOnNx = interp.getNearestPointIndex(nx, x)
    closestPositionToPoint = trackSurface.createPositionProportion(nxProportion1,
                                                                   closestIdxOnNx / elementsCountAlongTrackSurface)
    xPosition = trackSurface.findNearestPosition(x, closestPositionToPoint)
    d = trackSurface.evaluateCoordinates(xPosition, derivatives=True)[1]

    return xPosition, d


def getSmoothedSampledPointsOnTrackSurface(trackSurface, startProportion1, startProportion2, endProportion1,
                                           endProportion2, elementsOut, startDerivative=None, endDerivative=None,
                                           startDerivativeMagnitude=None, endDerivativeMagnitude=None,
                                           curveMode=TrackSurface.HermiteCurveMode.SMOOTH):
    """
    Create smoothly spaced out hermite curve points between two points a and b on the surface,
    each defined by their proportions over the surface in directions 1 and 2.
    :param curveMode:
    :param trackSurface: track surface
    :param startProportion1, startProportion2: proportion of start point in direction around and along track surface
    :param endProportion1, endProportion2: proportion of end point in direction around and along track surface
    :param elementsOut: number of elements out
    :param startDerivative, endDerivative: optional derivative vectors in 3-D world coordinates
        to match at the start and end of the curves. If omitted, fits in with other derivative or is
        in a straight line from a to b
    :param startDerivativeMagnitude, endDerivativeMagnitude: optional magnitude of derivatives to match at the start and
        end of the curves
    :return: coordinates and derivative of sampled points
    """

    mx, md2, md1, md3, mProportions = \
        trackSurface.createHermiteCurvePoints(startProportion1, startProportion2, endProportion1, endProportion2,
                                              elementsOut, startDerivative, endDerivative, curveMode)

    xSampled, dSampled = trackSurface.resampleHermiteCurvePointsSmooth(mx, md2, md1, md3, mProportions,
                                                                       startDerivativeMagnitude,
                                                                       endDerivativeMagnitude)[0:2]
    return xSampled, dSampled


def smoothD1Around(xAround, d1Around):
    """
    Rearrange points around so that the group of points starts from left side of the annulus to the greater curvature
    and ends on the right side of the annulus. This put points in consecutive order for derivative smoothing.
    The smoothed derivatives are then re-arranged such that it starts from the greater curvature and goes to the right
    side of the annulus, followed by the left side of the annulus and closing the loop back at the greater curvature.
    :param xAround: points around a loop joining to the annulus.
    :param d1Around: derivative of points.
    :return: smoothed derivatives in the original order.
    """
    xLoop = xAround[int(len(xAround) * 0.5 + 1):] + xAround[: int(len(xAround) * 0.5 + 1)]
    d1Loop = d1Around[int(len(d1Around) * 0.5 + 1):] + d1Around[: int(len(d1Around) * 0.5 + 1)]
    d1LoopSmooth = interp.smoothCubicHermiteDerivativesLine(xLoop, d1Loop, fixStartDerivative=True,
                                                            fixEndDerivative=True)
    # Rearrange to correct order
    d1Around = d1LoopSmooth[int(len(xAround) * 0.5):] + d1LoopSmooth[: int(len(xAround) * 0.5):]

    return d1Around


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

def createStomachMesh3d(region, fm, coordinates, stomachTermsAlong, elementsCountAlongTrackSurface, allAnnotationGroups,
                        elementCountGroupList, centralPath, options, nodeIdentifier, elementIdentifier,
                        splitCoordinates, materialCoordinates=False):
    """
    Generates a stomach scaffold in the region using a central path and parameter options.
    :param region: Region to create elements in.
    :param fm: Zinc fieldModule to create elements in.
    :param coordinates: Coordinate field to define nodes and elements.
    :param stomachTermsAlong: Annotation terms along length of stomach.
    :param elementsCountAlongTrackSurface: Number of elements to use for track surface.
    :param allAnnotationGroups: List of annotation groups.
    :param ADD!!!!!!!!!!!!!
    :param ADD!!!!!!!!!!!!!
    :param centralPath: Central path through the axis of the stomach scaffold.
    :param options: Parameter options for stomach scaffold.
    :param nodeIdentifier: First node identifier.
    :param elementIdentifier: First element identifier.
    :param splitCoordinates: Create split coordinates if True.
    :param ADD!!!!!!!!!!!!!
    :return allAnnotationGroups
    """
    elementsCountAroundEso = options['Number of elements around esophagus']
    elementsCountAroundDuod = options['Number of elements around duodenum']
    elementsAlongFundusApexToCardia = options['Number of elements between fundus apex and cardia']
    print('FA to cardia', elementsAlongFundusApexToCardia) # Check that this is updated on UI
    elementsAlongCardiaToDuod = options['Number of elements between cardia and duodenum']
    elementsCountThroughWall = options['Number of elements through wall']
    wallThickness = options['Wall thickness']
    mucosaRelThickness = options['Mucosa relative thickness']
    submucosaRelThickness = options['Submucosa relative thickness']
    circularRelThickness = options['Circular muscle layer relative thickness']
    longitudinalRelThickness = options['Longitudinal muscle layer relative thickness']
    useCrossDerivatives = False
    useCubicHermiteThroughWall = not (options['Use linear through wall'])

    # GEJPositionAlongFactor = options['Gastro-esophagal junction position along factor']
    GEJOptions = options['Gastro-esophagal junction']
    GEJSettings = GEJOptions.getScaffoldSettings()
    elementsAlongEsophagus = GEJSettings['Number of elements along']
    elementsThroughEsophagusWall = GEJSettings['Number of elements through wall']
    ostiumDiameter = GEJSettings['Ostium diameter']
    limitingRidge = options['Limiting ridge']
    elementsCountAcrossCardia = options['Number of elements across cardia']
    cardiaDiameterFactor = options['Cardia diameter factor'] # scale to ostium diameter

    elementsAroundHalfEso = int(elementsCountAroundEso * 0.5)
    elementsAroundQuarterEso = int(elementsCountAroundEso * 0.25)
    elementsAroundHalfDuod = int(elementsCountAroundDuod * 0.5)
    elementsAroundQuarterDuod = int(elementsCountAroundDuod * 0.25)

    sf = (cardiaDiameterFactor - 1) * ostiumDiameter * 0.5

    zero = [0.0, 0.0, 0.0]
    seeTrackSurface = options['Track surface'] # KM

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
        GEJPositionAlongFactor = 0.4

        ostiumLength = GEJSettings['Ostium length']
        ostiumWallThickness = GEJSettings['Ostium wall thickness']
        ostiumWallRelThicknesses = GEJSettings['Ostium wall relative thicknesses']
        vesselInnerDiameter = GEJSettings['Vessel inner diameter']
        vesselWallThickness = GEJSettings['Vessel wall thickness']
        vesselWallRelThicknesses = GEJSettings['Vessel wall relative thicknesses']

        GEJSettings['Ostium diameter'] = 0.3
        GEJSettings['Ostium length'] = 0.3
        GEJSettings['Ostium wall thickness'] = wallThickness
        GEJSettings['Ostium wall relative thicknesses'] = relThicknesses
        GEJSettings['Vessel inner diameter'] = 0.1
        GEJSettings['Vessel wall thickness'] = wallThickness * 0.6
        GEJSettings['Vessel wall relative thicknesses'] = relThicknesses
        sf = (cardiaDiameterFactor - 1) * GEJSettings['Ostium diameter'] * 0.5

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

    # fundusEndPositionAlongFactor = arcLengthOfGroupsAlong[1] / stomachCentralPathLength
    for i in range(1, len(stomachTermsAlong)):
        arcLengthRatio = arcLengthOfGroupsAlong[i] / stomachCentralPathLength
        arcLengthRatioForGroupsFromFundusApex.append(arcLengthRatio)

    stomachGroup = AnnotationGroup(region, get_stomach_term("stomach"))
    fundusGroup = AnnotationGroup(region, get_stomach_term("fundus of stomach"))
    bodyGroup = AnnotationGroup(region, get_stomach_term("body of stomach"))
    antrumGroup = AnnotationGroup(region, get_stomach_term("pyloric antrum"))
    pylorusGroup = AnnotationGroup(region, get_stomach_term("pyloric canal"))
    duodenumGroup = AnnotationGroup(region, get_stomach_term("duodenum"))

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

    # annotation fiducial points
    markerGroup = findOrCreateFieldGroup(fm, "marker")
    markerName = findOrCreateFieldStoredString(fm, name="marker_name")
    markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
    markerTemplateInternal = nodes.createNodetemplate()
    markerTemplateInternal.defineField(markerName)
    markerTemplateInternal.defineField(markerLocation)

    # Extract boundary planes between groups
    cxGroupBoundaries = [centralPath.cxGroups[0][0]]
    cd1GroupBoundaries = [centralPath.cd1Groups[0][0]]
    cd2GroupBoundaries = [centralPath.cd2Groups[0][0]]
    cd3GroupBoundaries = [centralPath.cd3Groups[0][0]]
    cd12GroupBoundaries = [centralPath.cd12Groups[0][0]]
    cd13GroupBoundaries = [centralPath.cd13Groups[0][0]]

    for i in range(1, len(stomachTermsAlong)):
        cxGroupBoundaries.append(centralPath.cxGroups[i][-1])
        cd1GroupBoundaries.append(centralPath.cd1Groups[i][-1])
        cd2GroupBoundaries.append(centralPath.cd2Groups[i][-1])
        cd3GroupBoundaries.append(centralPath.cd3Groups[i][-1])
        cd12GroupBoundaries.append(centralPath.cd12Groups[i][-1])
        cd13GroupBoundaries.append(centralPath.cd13Groups[i][-1])

    # elementsAlongFromBody cannot be less than (number of annotation groups - fundus)
    elementsAlongFromBody = elementsAlongCardiaToDuod + elementsAroundQuarterEso - 1
    arcLengthOfGroupsAlongFromBody = arcLengthOfGroupsAlong[0] - arcLengthOfGroupsAlong[1]
    estElementLengthFromBody = arcLengthOfGroupsAlongFromBody / elementsAlongFromBody

    if not materialCoordinates:
        modGroups = []
        elementsAlongCPFundus = elementsAlongFundusApexToCardia + elementsAroundQuarterEso - 1
        elementCountGroupList = [elementsAlongCPFundus]
        elementsTally = 0
        for i in range(2, len(stomachTermsAlong)):
            numberOfElementsGroup = int(arcLengthOfGroupsAlong[i] // estElementLengthFromBody)
            mod = arcLengthOfGroupsAlong[i] % estElementLengthFromBody
            modGroups.append(mod)
            if numberOfElementsGroup < 1:
                numberOfElementsGroup = 1
            elementsTally += numberOfElementsGroup
            # print(numberOfElementsGroup, mod)
            elementCountGroupList.append(numberOfElementsGroup)

        excessElements = elementsAlongFromBody - elementsTally
        for i in range(excessElements):
            maxIdx = min(range(len(modGroups)), key=modGroups.__getitem__)
            # print('idx', maxIdx)
            elementCountGroupList[maxIdx - 1] += 1
            modGroups[maxIdx - 1] = 0

    print('GroupList', elementCountGroupList)

    # Break central path into elements allocation to each group
    cxSections = []
    cd1Sections = []
    cd2Sections = []
    cd3Sections = []

    for i in range(len(elementCountGroupList)):
        cxBoundary = [cxGroupBoundaries[i], cxGroupBoundaries[i + 1]]
        cd1Boundary = [cd1GroupBoundaries[i], cd1GroupBoundaries[i + 1]]
        cd2Boundary = [cd2GroupBoundaries[i], cd2GroupBoundaries[i + 1]]
        cd3Boundary = [cd3GroupBoundaries[i], cd3GroupBoundaries[i + 1]]
        cd12Boundary = [cd12GroupBoundaries[i], cd12GroupBoundaries[i + 1]]
        cd13Boundary = [cd13GroupBoundaries[i], cd13GroupBoundaries[i + 1]]

        # Sample
        elementsOutSection = elementCountGroupList[i]
        cxSection, cd1Section, pe, pxi, psf = \
            interp.sampleCubicHermiteCurvesSmooth(cxBoundary, cd1Boundary, elementsOutSection,
                                                  derivativeMagnitudeStart=
                                                  (vector.magnitude(cd1Sections[-1][-1]) if i else None),
                                                  derivativeMagnitudeEnd=None)
        cd2Section = interp.interpolateSampleCubicHermite(cd2Boundary, cd12Boundary, pe, pxi, psf)[0]
        cd3Section = interp.interpolateSampleCubicHermite(cd3Boundary, cd13Boundary, pe, pxi, psf)[0]
        cxSections.append(cxSection)
        cd1Sections.append(cd1Section)
        cd2Sections.append(cd2Section)
        cd3Sections.append(cd3Section)

    # Create track surface
    # Determine radius along stomach
    cx = centralPath.cxGroups[0]
    cd1 = centralPath.cd1Groups[0]
    cd2 = centralPath.cd2Groups[0]
    cd3 = centralPath.cd3Groups[0]
    cd12 = centralPath.cd12Groups[0]
    cd13 = centralPath.cd13Groups[0]

    sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlongTrackSurface)
    sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)
    sd3, sd13 = interp.interpolateSampleCubicHermite(cd3, cd13, se, sxi, ssf)

    # # Visualise central path
    # for n in range(len(sx)):
    #     node = nodes.createNode(nodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sx[n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sd1[n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sd2[n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero) #sd3[n])
    #     nodeIdentifier += 1

    fundusRadius = vector.magnitude(sd2[0])
    lengthElementAlongTrackSurface = stomachCentralPathLength / elementsCountAlongTrackSurface
    elementsAlongFundus = int(fundusRadius / lengthElementAlongTrackSurface)

    d2Apex = []
    d2 = sd2[0]
    for n1 in range(elementsCountAroundDuod):
        rotAngle = n1 * 2.0 * math.pi / elementsCountAroundDuod
        rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(sd2[0]), vector.normalise(sd3[0])))
        rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, rotAngle)
        d2Rot = [rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)]
        d2Apex.append(d2Rot)

    xEllipses = []
    d1Ellipses = []
    for n in range(elementsAlongFundus + 1, len(sx)):
        px, pd1 = createEllipsePoints(sx[n], 2 * math.pi, sd2[n], sd3[n], elementsCountAroundDuod,
                                      startRadians=0.0)
        xEllipses.append(px)
        d1Ellipses.append(pd1)

    # Find d2
    d2Raw = []
    for n1 in range(elementsCountAroundDuod):
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
        for n1 in range(elementsCountAroundDuod):
            d2 = d2Raw[n1][n2]
            d2Around.append(d2)
        d2Ellipses.append(d2Around)

    # Merge fundus apex and body
    xAll = [[sx[0]] * elementsCountAroundDuod] + xEllipses
    d2All = [d2Apex] + d2Ellipses

    # Spread out elements
    xRaw = []
    d2Raw = []

    for n1 in range(elementsCountAroundDuod):
        xAlong = []
        d2Along = []
        for n2 in range(len(xAll)):
            xAlong.append(xAll[n2][n1])
            d2Along.append(d2All[n2][n1])
        xSampledAlong, d2SampledAlong = interp.sampleCubicHermiteCurves(xAlong, d2Along,
                                                                        elementsCountAlongTrackSurface,
                                                                        arcLengthDerivatives=True)[0:2]

        d2Smoothed = interp.smoothCubicHermiteDerivativesLine(xSampledAlong, d2SampledAlong)
        xRaw.append(xSampledAlong)
        d2Raw.append(d2Smoothed)

        if n1 == 0:
            xSampledGC = xSampledAlong
            d2SampledGC = d2Smoothed
            arcLengthGC = 0.0
            for n in range(len(xSampledGC) - 1):
                arcLengthGC += interp.getCubicHermiteArcLength(xSampledGC[n], d2Smoothed[n],
                                                               xSampledGC[n + 1], d2Smoothed[n + 1])
        elif n1 == elementsAroundQuarterDuod:
            xSampledQuarterDuod = xSampledAlong
            d2SampledQuarterDuod = d2Smoothed
            arcLengthQuarterDuod = 0.0
            for n in range(len(xSampledQuarterDuod) - 1):
                arcLengthQuarterDuod += interp.getCubicHermiteArcLength(xSampledQuarterDuod[n], d2Smoothed[n],
                                                                        xSampledQuarterDuod[n + 1], d2Smoothed[n + 1])

    # Rearrange x and d2
    xSampledAll = []
    d1SampledAll = []
    d2SampledAll = []
    for n2 in range(elementsCountAlongTrackSurface + 1):
        xAround = []
        d1Around = []
        d2Around = []
        for n1 in range(elementsCountAroundDuod):
            x = xRaw[n1][n2]
            d2 = d2Raw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)

            # Calculate d1
            if n2 > 0:
                v1 = xRaw[n1][n2]
                v2 = xRaw[n1 + 1 if n1 < elementsCountAroundDuod - 2 else 0][n2]
                d1 = findDerivativeBetweenPoints(v1, v2)
                d1Around.append(d1)
            else:
                d1Around.append(d2Raw[int(elementsCountAroundDuod * 0.75)][0])

        if n2 > 0:
            d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, d1Around)
        else:
            d1Smoothed = d1Around

        xSampledAll.append(xAround)
        d1SampledAll.append(d1Smoothed)
        d2SampledAll.append(d2Around)

    # Create tracksurface
    xTrackSurface = []
    d1TrackSurface = []
    d2TrackSurface = []
    for n2 in range(elementsCountAlongTrackSurface + 1):
        for n1 in range(elementsCountAroundDuod):
            xTrackSurface.append(xSampledAll[n2][n1])
            d1TrackSurface.append(d1SampledAll[n2][n1])
            d2TrackSurface.append(d2SampledAll[n2][n1])

    trackSurfaceStomach = TrackSurface(elementsCountAroundDuod, elementsCountAlongTrackSurface,
                                       xTrackSurface, d1TrackSurface, d2TrackSurface, loop1=True)

    # Visualise track surface
    if seeTrackSurface:
        for n1 in range(len(xTrackSurface)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTrackSurface[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2TrackSurface[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d1TrackSurface[n1])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            nodeIdentifier += 1

    # Find closest point to cd2 and cd3 on track surface for fundus
    arcLengthPerSampledElementGC = arcLengthGC/elementsCountAlongTrackSurface
    arcLengthPerSampledElementQuarterDuod = arcLengthQuarterDuod/elementsCountAlongTrackSurface

    findRescaledCrossAxis(xSampledGC, d2SampledGC, arcLengthPerSampledElementGC, cxSections[0], cd2Sections[0])
    findRescaledCrossAxis(xSampledQuarterDuod, d2SampledQuarterDuod, arcLengthPerSampledElementQuarterDuod,
                          cxSections[0], cd3Sections[0])

    xEllipseAroundAll = []
    d1EllipseAroundAll = []
    d2EllipseAroundAll = []

    for s in range(len(cxSections)):
        for n2 in range(1, len(cxSections[s])):
            px, pd1 = createEllipsePoints(cxSections[s][n2], 2 * math.pi, cd2Sections[s][n2], cd3Sections[s][n2],
                                          elementsCountAroundDuod, startRadians=0.0)
            d2Around = [zero for n in range(len(pd1))]

            xEllipseAroundAll.append(px)
            d1EllipseAroundAll.append(pd1)
            d2EllipseAroundAll.append(d2Around)

            if s == 0 and n2 == len(cxSections[s]) - 1:
                xGEJ = px[elementsAroundHalfDuod]

            # for n1 in range(len(px)):
            #     node = nodes.createNode(nodeIdentifier, nodetemplate)
            #     cache.setNode(node)
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, px[n1])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, pd1[n1])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            #     nodeIdentifier += 1

    # for s in range(len(cxSections)):
    #     for n2 in range(len(cxSections[s])):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxSections[s][n2])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd2Sections[s][n2])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd3Sections[s][n2])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
    #         nodeIdentifier += 1

    # Set up gastro-esophagal junction with midpoint aligned to fundus-body junction
    GEJSettings['Number of elements around ostium'] = elementsCountAroundEso
    GEJPosition = trackSurfaceStomach.findNearestPosition(xGEJ)
    xCentre, d1Centre, d2Centre = trackSurfaceStomach.evaluateCoordinates(GEJPosition, derivatives=True)
    axis1 = d1Centre

    esophagusGroup = AnnotationGroup(region, get_stomach_term("esophagus"))
    esophagusMeshGroup = esophagusGroup.getMeshGroup(mesh)
    esophagogastricJunctionGroup = AnnotationGroup(region, get_stomach_term("esophagogastric junction"))
    esophagogastricJunctionMeshGroup = esophagogastricJunctionGroup.getMeshGroup(mesh)
    stomachMeshGroup = stomachGroup.getMeshGroup(mesh)
    allAnnotationGroups += [esophagusGroup, esophagogastricJunctionGroup]

    ostiumWallAnnotationGroups = []
    if elementsCountThroughWall == 4:
        esophagusMucosaGroup = AnnotationGroup(region, get_stomach_term("esophagus mucosa"))
        esophagusSubmucosaGroup = AnnotationGroup(region, get_stomach_term("submucosa of esophagus"))
        esophagusCircularGroup = AnnotationGroup(region, get_stomach_term("esophagus smooth muscle circular layer"))
        esophagusLongitudinalGroup = AnnotationGroup(region,
                                                     get_stomach_term("esophagus smooth muscle longitudinal layer"))

        ostiumWallAnnotationGroups = [[esophagusMucosaGroup, mucosaGroup],
                                      [esophagusSubmucosaGroup, submucosaGroup],
                                      [esophagusCircularGroup, circularMuscleGroup],
                                      [esophagusLongitudinalGroup, longitudinalMuscleGroup]]

        allAnnotationGroups += [esophagusMucosaGroup, esophagusSubmucosaGroup,
                                esophagusCircularGroup, esophagusLongitudinalGroup]

    nextNodeIdentifier, nextElementIdentifier, (o1_x, o1_d1, o1_d2, o1_d3, o1_NodeId, o1_Positions) = \
        generateOstiumMesh(region, GEJSettings, trackSurfaceStomach, GEJPosition, axis1,
                           nodeIdentifier, elementIdentifier,
                           vesselMeshGroups=[[stomachMeshGroup, esophagusMeshGroup]],
                           ostiumMeshGroups=[stomachMeshGroup, esophagogastricJunctionMeshGroup],
                           wallAnnotationGroups=ostiumWallAnnotationGroups, coordinates=coordinates)

    if materialCoordinates:
        GEJSettings['Ostium diameter'] = ostiumDiameter
        GEJSettings['Ostium length'] = ostiumLength
        GEJSettings['Ostium wall thickness'] = ostiumWallThickness
        GEJSettings['Ostium wall relative thicknesses'] = ostiumWallRelThicknesses
        GEJSettings['Vessel inner diameter'] = vesselInnerDiameter
        GEJSettings['Vessel wall thickness'] = vesselWallThickness
        GEJSettings['Vessel wall relative thicknesses'] = vesselWallRelThicknesses

    stomachStartNode = nextNodeIdentifier
    # stomachStartElement = nextElementIdentifier
    nodeIdentifier = nextNodeIdentifier
    elementIdentifier = nextElementIdentifier

    # Create location of annulus
    xAnnulusOuter = [[] for x in range(elementsCountAroundEso)]
    xAnnulusOuterPosition = [[] for x in range(elementsCountAroundEso)]
    d1AnnulusNorm = []
    for n1 in range(elementsCountAroundEso):
        normD2 = vector.normalise(o1_d2[-1][n1])
        d1AnnulusNorm.append(normD2)
        x = [o1_x[-1][n1][c] + sf * normD2[c] for c in range(3)]
        nearestPosition = trackSurfaceStomach.findNearestPosition(x)
        xAnnulusOuterPosition[n1] = nearestPosition
        xAnnulusOuter[n1] = trackSurfaceStomach.evaluateCoordinates(nearestPosition)

    dAnnulusOuter = []
    for n in range(elementsCountAroundEso):
        d = findDerivativeBetweenPoints(xAnnulusOuter[n], xAnnulusOuter[(n + 1) % elementsCountAroundEso])
        dAnnulusOuter.append(d)
    dAnnulusOuter = interp.smoothCubicHermiteDerivativesLoop(xAnnulusOuter, dAnnulusOuter)
    d3Annulus = []
    for n in range(elementsCountAroundEso):
        d3Annulus.append(vector.normalise(vector.crossproduct3(vector.normalise(dAnnulusOuter[n]), d1AnnulusNorm[n])))
    annulusCurvature = findCurvatureAroundLoop(xAnnulusOuter, dAnnulusOuter, d3Annulus)

    # Check if we need to smooth it with its path?? and points[0] & [halfEso]

    # for n1 in range(len(xAnnulusOuter)):
    #     node = nodes.createNode(nodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAnnulusOuter[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dAnnulusOuter[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
    #     # print('annulus', nodeIdentifier)
    #     nodeIdentifier += 1

    esoCount = 2
    xAroundB = []
    d1AroundB = []
    d2AroundB = []
    xAroundC = []
    d1AroundC = []
    d2AroundC = []
    xAroundD = []
    d1AroundD = []
    d2AroundD = []

    xAroundE = []
    d1AroundE = []
    d2AroundE = []

    xAroundF = []
    d1AroundF = []
    d2AroundF = []

    for n2 in range(len(xEllipseAroundAll)):
        # straight from LC and GC to end of duodenum
        if n2 <= elementsAroundQuarterDuod - 3: # xAroundA - DEAL WITH DUOD > 12
            print('link to duod - ', n2)
            xEllipseNew, d1EllipseNew = createEllipsePoints(cxSections[0][1], 2*math.pi, cd2Sections[0][1],
                                                            cd3Sections[0][1], elementsCountAroundDuod - 4)
            xEllipseAroundAll[n2] = xEllipseNew
            d1EllipseAroundAll[n2] = d1EllipseNew
            d2Around = [zero for n in range(len(d1EllipseNew))]
            d2EllipseAroundAll[n2] = d2Around

            xAroundA = xEllipseNew
            d1AroundA = d1EllipseNew
            d2AroundA = d2Around

        elif elementsAroundQuarterDuod - 3 < n2 < elementsAlongFundusApexToCardia - 1: # AroundB?
            xAroundB.append(xEllipseAroundAll[n2])
            d1AroundB.append(d1EllipseAroundAll[n2])
            d2AroundB.append(d2EllipseAroundAll[n2])
            print('complete ring before eso with no change -', n2) # take complete ellipse

        elif n2 == elementsAlongFundusApexToCardia - 1: # AroundC
            print('first ring joining to annulus[0] and annulus[1] -', n2)
            xEllipseAroundAll[n2][elementsAroundHalfDuod - 1] = xAnnulusOuter[1]
            xEllipseAroundAll[n2][elementsAroundHalfDuod] = xAnnulusOuter[0]
            xEllipseAroundAll[n2][elementsAroundHalfDuod + 1] = xAnnulusOuter[-1]
            xAroundC.append(xEllipseAroundAll[n2])
            d1AroundC.append(d1EllipseAroundAll[n2])
            d2AroundC.append(d2EllipseAroundAll[n2])
        elif elementsAlongFundusApexToCardia - 1 < n2 < elementsAlongFundusApexToCardia - 1 + elementsAroundHalfEso - 2: # AroundD
            print('ring interior of eso -', n2)
            xEllipseAroundAll[n2][elementsAroundHalfDuod - 1] = xAnnulusOuter[esoCount]
            xEllipseAroundAll[n2][elementsAroundHalfDuod + 1] = xAnnulusOuter[-esoCount]
            esoCount += 1
            del xEllipseAroundAll[n2][elementsAroundHalfDuod]
            del d1EllipseAroundAll[n2][elementsAroundHalfDuod]
            del d2EllipseAroundAll[n2][elementsAroundHalfDuod]
            xAroundD.append(xEllipseAroundAll[n2])
            d1AroundD.append(d1EllipseAroundAll[n2])
            d2AroundD.append(d2EllipseAroundAll[n2])
        elif n2 == elementsAlongFundusApexToCardia - 1 + elementsAroundHalfEso - 2: # AroundE
            print('last ring joining to annulus[0] and annulus[1] -', n2)
            xEllipseAroundAll[n2][elementsAroundHalfDuod - 1] = xAnnulusOuter[elementsAroundHalfEso - 1]
            xEllipseAroundAll[n2][elementsAroundHalfDuod] = xAnnulusOuter[elementsAroundHalfEso]
            xEllipseAroundAll[n2][elementsAroundHalfDuod + 1] = xAnnulusOuter[elementsAroundHalfEso + 1]
            xAroundE.append(xEllipseAroundAll[n2])
            d1AroundE.append(d1EllipseAroundAll[n2])
            d2AroundE.append(d2EllipseAroundAll[n2])
        else:
            print('complete rings beyond eso -', n2) # AroundF
            xAroundF.append(xEllipseAroundAll[n2])
            d1AroundF.append(d1EllipseAroundAll[n2])
            d2AroundF.append(d2EllipseAroundAll[n2])

    # for n2 in range(len(xEllipseAroundAll)):
    #     for n1 in range(len(xEllipseAroundAll[n2])):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xEllipseAroundAll[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1EllipseAroundAll[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
    #         nodeIdentifier += 1

    # Along GC
    xAlongGC = []
    annulusD2 = vector.normalise(o1_d2[-1][0])
    d2 = [sf * annulusD2[c] for c in range(3)]
    d2AlongGC = [d2]
    for n2 in range(elementsAlongFundusApexToCardia - 1, -1, -1):
        xAlongGC.append(xEllipseAroundAll[n2][int(0.5 * len(xEllipseAroundAll[n2]))])
    xAlongGC.append(cxSections[0][0])
    for n2 in range(len(xEllipseAroundAll)):
        xAlongGC.append(xEllipseAroundAll[n2][0])
    for n in range(len(xAlongGC) - 1):
        d2AlongGC.append(findDerivativeBetweenPoints(xAlongGC[n], xAlongGC[n + 1]))
    d2AlongGCSmoothed = interp.smoothCubicHermiteDerivativesLine(xAlongGC, d2AlongGC, fixStartDerivative=True)
                                                                 #fixEndDerivative=True)

    xAlongGCReverse = copy.deepcopy(xAlongGC)
    xAlongGCReverse.reverse()
    d2AlongGCReverse = []
    for n in range(len(xAlongGCReverse) - 1):
        d = findDerivativeBetweenPoints(xAlongGCReverse[n], xAlongGCReverse[n + 1])
        d2AlongGCReverse.append(d)
    rotAxis = vector.normalise(vector.crossproduct3(vector.normalise(o1_d1[-1][0]), vector.normalise(o1_d2[-1][0])))
    normD2 = vector.normalise(o1_d2[-1][0])
    d2 = [sf * normD2[c] for c in range(3)]
    rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
    d2AlongGCReverse.append([rotFrame[j][0] * d2[0] + rotFrame[j][1] * d2[1] + rotFrame[j][2] * d2[2] for j in range(3)])
    d2AlongGCReverseSmoothed = interp.smoothCubicHermiteDerivativesLine(xAlongGCReverse, d2AlongGCReverse)
    #fixEndDerivative=True)

    # Along LC
    startIdx = elementsAlongFundusApexToCardia + elementsAroundHalfEso - 3
    xAlongLC = []
    d2 = vector.normalise(o1_d2[-1][elementsAroundHalfEso])
    d2 = [sf * d2[c] for c in range(3)]
    d2AlongLC = [d2]
    for n2 in range(elementsAlongCardiaToDuod + 1):
        idx = startIdx + n2
        xAlongLC.append(xEllipseAroundAll[idx][int(0.5 * len(xEllipseAroundAll[idx]))])
    for n in range(len(xAlongLC) - 1):
        d2AlongLC.append(findDerivativeBetweenPoints(xAlongLC[n], xAlongLC[n + 1]))
    d2AlongLCSmoothed = interp.smoothCubicHermiteDerivativesLine(xAlongLC, d2AlongLC, fixStartDerivative=True)

    # Along quarter line
    xSampledQuarterLine = [cxSections[0][0]]
    d2SampledQuarterLine = [cd3Sections[0][0]]
    for n2 in range(len(xEllipseAroundAll)):
        xSampledQuarterLine.append(xEllipseAroundAll[n2][elementsAroundQuarterDuod if n2 else int(0.25 * len(xEllipseAroundAll[n2]))])
    for n in range(len(xSampledQuarterLine) - 1):
        d2SampledQuarterLine.append(findDerivativeBetweenPoints(xSampledQuarterLine[n], xSampledQuarterLine[n + 1]))
    d2SampledQuarterLineSmoothed = interp.smoothCubicHermiteDerivativesLine(xSampledQuarterLine, d2SampledQuarterLine,
                                                                            fixStartDirection=True)

    # Along three quarter line
    xSampledThreeQuarterLine = [cxSections[0][0]]
    d = cd3Sections[0][0]
    rotAxis = vector.crossproduct3(vector.normalise(cd2Sections[0][0]), vector.normalise(cd3Sections[0][0]))
    rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
    d = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
    d2SampledThreeQuarterLine = [d]
    for n2 in range(len(xEllipseAroundAll)):
        xSampledThreeQuarterLine.append(xEllipseAroundAll[n2][-elementsAroundQuarterDuod if n2 else -int(0.25 * len(xEllipseAroundAll[n2]))])
    for n in range(len(xSampledThreeQuarterLine) - 1):
        d2SampledThreeQuarterLine.append(findDerivativeBetweenPoints(xSampledThreeQuarterLine[n],
                                                                     xSampledThreeQuarterLine[n + 1]))
    d2SampledThreeQuarterLineSmoothed = interp.smoothCubicHermiteDerivativesLine(xSampledThreeQuarterLine,
                                                                                 d2SampledThreeQuarterLine,
                                                                                 fixStartDirection=True)

    # for n2 in range(len(xSampledThreeQuarterLine)):
    #     node = nodes.createNode(nodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xSampledThreeQuarterLine[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2SampledThreeQuarterLineSmoothed[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
    #     print('xSampledThreeQuarterLine', nodeIdentifier)
    #     nodeIdentifier += 1

    # Assemble xOuter
    xOuter = []
    d1Outer = []
    d2Outer = []
    x0 = []
    d0 = []
    # Apex
    for n in range(elementsAroundHalfDuod - 5):
        GCIdx = elementsAlongFundusApexToCardia - (elementsAroundQuarterDuod - 2) + 1 + n
        x0.append(xAlongGC[GCIdx])
        d0.append(findDerivativeBetweenPoints(xAlongGC[GCIdx], xAroundA[int(len(xAroundA) * 0.5) - 2 - n]))
    x0Reverse = copy.deepcopy(x0)
    x0Reverse.reverse()

    # Rearrange x0 - apex first, next level along, GC first
    x0Arranged = []
    d0Arranged = []
    x0Arranged.append(x0[int(len(x0) * 0.5)])
    # d0Arranged.append(d2SampledQuarterLineSmoothed[0])
    d0Arranged.append(findDerivativeBetweenPoints(cxSections[0][0],
                                                  xEllipseAroundAll[0][int(0.25 * len(xEllipseAroundAll[0]))]))
    for n2 in range(1, int(len(x0) * 0.5) + 1):
        x0Arranged.append(x0[int(len(x0) * 0.5) + n2])
        x0Arranged.append(x0[int(len(x0) * 0.5) - n2])
        d0Arranged.append(d0[int(len(d0) * 0.5) + n2])
        d0Arranged.append(d0[int(len(d0) * 0.5) - n2])
    xOuter.append(x0Arranged)
    d1Outer.append(d0Arranged)
    d2Outer.append([zero for c in range(len(d0Arranged))])

    xOuter += [xAroundA] + xAroundB + xAroundC + xAroundD + xAroundE + xAroundF
    d1Outer += [d1AroundA] + d1AroundB + d1AroundC + d1AroundD + d1AroundE + d1AroundF
    d2Outer += [d2AroundA] + d2AroundB + d2AroundC + d2AroundD + d2AroundE + d2AroundF

    # completeRowsBeforeEso = len(xSection1) - elementsAroundQuarterDuod
    completeRowsBeforeEso = elementsAlongFundusApexToCardia + 1 - elementsAroundQuarterDuod
    # elementsAlongStomach = len(xSection1) - elementsAroundQuarterDuod + elementsAroundHalfEso + elementsAlongCardiaToDuod
    elementsAlongStomach = elementsAlongFundusApexToCardia + 1 - elementsAroundQuarterDuod + elementsAroundHalfEso + \
                           elementsAlongCardiaToDuod

    d2Curvature = []
    for n2 in range(len(d2Outer)):
        curvature = []
        for n1 in range(len(d2Outer[n2])):
            curvature.append(0.0)
        d2Curvature.append(curvature)

    xAlongAll = []
    d2AlongAll = []
    count = 0
    for n1 in range(elementsCountAroundDuod):
        xAlongN1 = []
        d2AlongN1 = []
        if n1 == 0: # GC
            xAlongN1 = xAlongGC[elementsAlongFundusApexToCardia:]
            d2AlongN1 = d2AlongGCSmoothed[elementsAlongFundusApexToCardia:]

        elif n1 == 1: # Start from 2nd point on quarterline
            for n in range(int(0.25 * len(xOuter[1])), 0, -1):
                xAlongN1.append(xOuter[1][n])
            for n2 in range(2, elementsAlongStomach + 1): # Template
                xAlongN1.append(xOuter[n2][n1])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1)

            # for n in range(len(d2AlongN1)):
            #     node = nodes.createNode(nodeIdentifier, nodetemplate)
            #     cache.setNode(node)
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAlongN1[n])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2AlongN1[n])
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            #     nodeIdentifier += 1

        elif n1 == 2: # Start from xOuter[1][0]
            xAlongsideGC = [xOuter[1][-1]] + xOuter[1][:2]
            d2AlongsideGC = []
            for n in range(len(xAlongsideGC) - 1):
                d = findDerivativeBetweenPoints(xAlongsideGC[n], xAlongsideGC[n + 1])
                d2AlongsideGC.append(d)
            d2AlongsideGC.append(d)
            d2AlongsideGC = interp.smoothCubicHermiteDerivativesLine(xAlongsideGC, d2AlongsideGC)

            xAlongN1 += xOuter[1][:2]
            for n2 in range(2, elementsAlongStomach + 1): # Template
                xAlongN1.append(xOuter[n2][n1])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1[0] = d2AlongsideGC[1]
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1, fixStartDirection=True)

        elif 2 < n1 < elementsAroundQuarterDuod or elementsAroundQuarterDuod < n1 < elementsAroundHalfDuod - 2:
            # Additional elements for duodenum above and below quarter line
            xAlongsideGC = [xOuter[1][int(len(xOuter[1]) * 0.5) + len(x0) - count + 1]] + [x0Reverse[count]] + \
                           [xOuter[1][2 + count]]
            d2AlongsideGC = []
            for n in range(len(xAlongsideGC) - 1):
                d = findDerivativeBetweenPoints(xAlongsideGC[n], xAlongsideGC[n + 1])
                d2AlongsideGC.append(d)
            d2AlongsideGC.append(d)
            d2AlongsideGC = interp.smoothCubicHermiteDerivativesLine(xAlongsideGC, d2AlongsideGC)
            xAlongN1 = [x0Reverse[count]] + [xOuter[1][2 + count]]
            for n2 in range(2, elementsAlongStomach + 1): # Template
                xAlongN1.append(xOuter[n2][n1])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1[0] = d2AlongsideGC[1]
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1, fixStartDirection=True)
            count += 1

        elif n1 == elementsAroundQuarterDuod:
            xAlongN1 = xSampledQuarterLine
            d2AlongN1 = d2SampledQuarterLineSmoothed
            count += 1

        elif n1 == elementsAroundHalfDuod - 2: # start from xOuter[1][half]
            xAlongsideGC = [xOuter[1][int(len(xOuter[1]) * 0.5) + 1]] + [xOuter[1][int(len(xOuter[1]) * 0.5)]] + \
                           [xOuter[1][int(len(xOuter[1]) * 0.5) - 1]]
            d2AlongsideGC = []
            for n in range(len(xAlongsideGC) - 1):
                d = findDerivativeBetweenPoints(xAlongsideGC[n], xAlongsideGC[n + 1])
                d2AlongsideGC.append(d)
            d2AlongsideGC.append(d)
            d2AlongsideGC = interp.smoothCubicHermiteDerivativesLine(xAlongsideGC, d2AlongsideGC)

            xAlongN1 = [xOuter[1][int(len(xOuter[1]) * 0.5)]] + [xOuter[1][int(len(xOuter[1]) * 0.5 - 1)]]
            for n2 in range(2, elementsAlongStomach + 1): # Template
                xAlongN1.append(xOuter[n2][n1])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1[0] = d2AlongsideGC[1]
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1, fixStartDirection=True)

        elif n1 == elementsAroundHalfDuod - 1:
            # Start from 2nd point on quarterline
            for n in range(int(0.25 * len(xOuter[1])), int(len(xOuter[1]) * 0.5)):
                xAlongN1.append(xOuter[1][n])
            for n2 in range(2, elementsAlongStomach + 1):
                xAlongN1.append(xOuter[n2][n1])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1)

            # Replace points on annulus with annulus d2
            annulusStartIdx = 1
            n2StartIdx = completeRowsBeforeEso + 2
            for n in range(elementsAroundHalfEso - 1):
                annulusIdx = annulusStartIdx + n
                n2Idx = n2StartIdx + n + int(0.5 * (len(xOuter[0]) - 1))
                d2AlongN1[n2Idx] = dAnnulusOuter[annulusIdx]

        elif n1 == elementsAroundHalfDuod: # Need to be processed as two parts!
            xAlongN1 = xAlongGCReverse[-(elementsAlongFundusApexToCardia + 1):] + xAlongLC
            d2AlongN1 = d2AlongGCReverseSmoothed[-(elementsAlongFundusApexToCardia + 1):] + d2AlongLCSmoothed

        elif n1 == elementsAroundHalfDuod + 1:
            # Start from 2nd point on threeQuarterline
            for n in range(int(0.75 * len(xOuter[1])), int(len(xOuter[1]) * 0.5), -1):
                xAlongN1.append(xOuter[1][n])
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                xAlongN1.append(xOuter[n2][n1Idx])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1)

            # Replace points on annulus with annulus d2
            annulusStartIdx = elementsAroundHalfEso + 1
            n2StartIdx = completeRowsBeforeEso + elementsAroundHalfEso
            for n in range(elementsAroundHalfEso - 1):
                annulusIdx = annulusStartIdx + n
                n2Idx = n2StartIdx - n + int(0.5 * (len(xOuter[0]) - 1))
                rotAxis = d3Annulus[annulusIdx]
                d = dAnnulusOuter[annulusIdx]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
                d = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
                d2AlongN1[n2Idx] = d

        elif n1 == elementsAroundHalfDuod + 2:
            # Start from xOuter[1][half]
            xAlongsideGC = [xOuter[1][int(len(xOuter[1]) * 0.5) - 1]] + [xOuter[1][int(len(xOuter[1]) * 0.5)]] + \
                           [xOuter[1][int(len(xOuter[1]) * 0.5) + 1]]
            d2AlongsideGC = []
            for n in range(len(xAlongsideGC) - 1):
                d = findDerivativeBetweenPoints(xAlongsideGC[n], xAlongsideGC[n + 1])
                d2AlongsideGC.append(d)
            d2AlongsideGC.append(d)
            d2AlongsideGC = interp.smoothCubicHermiteDerivativesLine(xAlongsideGC, d2AlongsideGC)

            xAlongN1 += [xOuter[1][int(len(xOuter[1]) * 0.5)]] + [xOuter[1][int(len(xOuter[1]) * 0.5) + 1]]
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                xAlongN1.append(xOuter[n2][n1Idx])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1[0] = d2AlongsideGC[1]
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1, fixStartDerivative=True)

        elif elementsAroundHalfDuod + 2 < n1 < int(elementsAroundQuarterDuod * 3) or \
                int(elementsAroundQuarterDuod * 3) < n1 < elementsCountAroundDuod - 2:
            count -= 1
            xAlongsideGC = [xOuter[1][2 + count]] + [x0Reverse[count]] + \
                           [xOuter[1][int(len(xOuter[1]) * 0.5) + len(x0) - count + 1]]
            d2AlongsideGC = []
            for n in range(len(xAlongsideGC) - 1):
                d = findDerivativeBetweenPoints(xAlongsideGC[n], xAlongsideGC[n + 1])
                d2AlongsideGC.append(d)
            d2AlongsideGC.append(d)
            d2AlongsideGC = interp.smoothCubicHermiteDerivativesLine(xAlongsideGC, d2AlongsideGC)

            # Additional elements for duodenum above and below three quarter line
            xAlongN1 = [x0Reverse[count]] + [xOuter[1][int(len(xOuter[1]) * 0.5) + len(x0) - count + 1]]
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                xAlongN1.append(xOuter[n2][n1Idx])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1[0] = d2AlongsideGC[1]
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1, fixStartDirection=True)

        elif n1 == int(elementsAroundQuarterDuod * 3):
            xAlongN1 = xSampledThreeQuarterLine
            d2AlongN1 = d2SampledThreeQuarterLineSmoothed
            count -= 1

        elif n1 == elementsCountAroundDuod - 2: # Start from xOuter[1][0]
            xAlongsideGC = [xOuter[1][1]] + [xOuter[1][0]] + [xOuter[1][-1]]
            d2AlongsideGC = []
            for n in range(len(xAlongsideGC) - 1):
                d = findDerivativeBetweenPoints(xAlongsideGC[n], xAlongsideGC[n + 1])
                d2AlongsideGC.append(d)
            d2AlongsideGC.append(d)
            d2AlongsideGC = interp.smoothCubicHermiteDerivativesLine(xAlongsideGC, d2AlongsideGC)

            xAlongN1 += [xOuter[1][0]] + [xOuter[1][-1]]
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                xAlongN1.append(xOuter[n2][n1Idx])
            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1[0] = d2AlongsideGC[1]
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1, fixStartDirection=True)

        elif n1 == elementsCountAroundDuod - 1: # Start from 2nd point on quarterline
            for n in range(int(0.25 * len(xOuter[1]))):
                xAlongN1.append(xOuter[1][int(0.75 * len(xOuter[1])) + n])
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                xAlongN1.append(xOuter[n2][n1Idx])

            for n in range(len(xAlongN1) - 1):
                d = findDerivativeBetweenPoints(xAlongN1[n], xAlongN1[n + 1])
                d2AlongN1.append(d)
            d2AlongN1.append(d)
            d2AlongN1 = interp.smoothCubicHermiteDerivativesLine(xAlongN1, d2AlongN1)

        xAlongAll.append(xAlongN1)
        d2AlongAll.append(d2AlongN1)

    # for n1 in range(len(xAlongAll)):
    #     for n2 in range(len(xAlongAll[n1])):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAlongAll[n1][n2])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d2AlongAll[n1][n2])
    #         nodeIdentifier += 1

    # Re-arrange back to xOuter
    for n2 in range(len(xOuter)):
        if n2 == 0:
            n1 = 0
            d2Outer[n2][n1] = d2AlongAll[int(3 * elementsAroundQuarterDuod)][n2]
            d1Outer[n2][n2] = d2AlongAll[n1][n2]
            if len(d2Outer[n2]) > 1:
                for n in range(int(len(xOuter[n2]) * 0.5)):
                    d2Outer[n2][2 * n + 1] = d2AlongAll[int(3 * elementsAroundQuarterDuod) - n - 1][n2]
                    d2Outer[n2][2 * n + 2] = d2AlongAll[int(3 * elementsAroundQuarterDuod) + n + 1][n2]
        elif n2 == 1:
            d2New = []
            count = 0
            for n in range(int(2 * (len(xOuter[0]) + 2))):
                n1AllIdx = elementsAroundQuarterDuod - int(0.5 * len(xOuter[0])) - 1 + n + count
                if n == 0:
                    d2New.append(d2AlongAll[-2][0])
                d2New.append(d2AlongAll[n1AllIdx][1])
                if n == len(xOuter[0]) + 1:
                    d2New.append(d2AlongAll[elementsAroundHalfDuod + 2][0])
                    count = 3
            d2Outer[n2] = d2New

        elif n2 > 1:
            d2New = []
            idxFromEnd = -(elementsAlongStomach + 1 - n2)
            for n1 in range(elementsCountAroundDuod +
                            (-1 if 3 + completeRowsBeforeEso <= n2 <
                                   3 + completeRowsBeforeEso + (elementsAroundHalfEso - 3) else 0)):
                if n1 == elementsAroundHalfDuod and 1 < n2 < 3 + completeRowsBeforeEso:
                    d2New.append(d2AlongAll[n1][int(0.5 * len(xOuter[0])) + n2])
                elif n1 >= elementsAroundHalfDuod \
                        and 3 + completeRowsBeforeEso <= n2 < 3 + completeRowsBeforeEso + (elementsAroundHalfEso - 3):
                    d2New.append(d2AlongAll[n1 + 1][idxFromEnd])
                else:
                    d2New.append(d2AlongAll[n1][idxFromEnd])
            d2Outer[n2] = d2New

    # Set d1 on xOuter[0] (apart from apex) to point along GC towards esophagus
    count = 0
    for n in range(1, len(xOuter[0]), 2):
        d1Outer[0][n] = d2AlongGCSmoothed[elementsAlongFundusApexToCardia + int(0.5 * n) + 1]
        d1Outer[0][n + 1] = d2AlongGCSmoothed[elementsAlongFundusApexToCardia - int(0.5 * n) - 1]

        n1IdxAlong = int((len(xOuter[0]) - 1) * 0.5) - count
        d2Outer[0][n] = d2AlongAll[int(3 * elementsAroundQuarterDuod) - n1IdxAlong][0]
        d2Outer[0][n + 1] = d2AlongAll[int(3 * elementsAroundQuarterDuod) + n1IdxAlong][0]
        count += 1

    # Set d1 on xOuter[1] on GC to point along GC towards esophagus
    if len(xOuter[0]) == 1:
        d1Outer[1][0] = d2AlongGCSmoothed[elementsAlongFundusApexToCardia + 1]
        d1Outer[1][int(0.5 * len(xOuter[1]))] = d2AlongGCSmoothed[elementsAlongFundusApexToCardia - 1]
    else:
        d1Outer[1][0] = d2AlongGCSmoothed[elementsAlongFundusApexToCardia + int(0.5 * n) + 2]
        d1Outer[1][int(0.5 * len(xOuter[1]))] = d2AlongGCSmoothed[elementsAlongFundusApexToCardia - int(0.5 * n) - 2]

    d1OuterForCurvature = copy.deepcopy(d1Outer)

    # Replace d1 at annulus
    annulusStartIdx = 1
    n2StartIdx = completeRowsBeforeEso + 2
    for n in range(elementsAroundHalfEso - 1):
        annulusIdx = annulusStartIdx + n
        n2Idx = n2StartIdx + n
        normD2 = vector.normalise(o1_d2[-1][annulusIdx])
        d = [sf * normD2[c] for c in range(3)]
        rotAxis = d3Annulus[annulusIdx]
        rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
        d = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
        d1Outer[n2Idx][elementsAroundHalfDuod - 1] = d

    annulusStartIdx = elementsAroundHalfEso + 1
    n2StartIdx = completeRowsBeforeEso + elementsAroundHalfEso
    for n in range(elementsAroundHalfEso - 1):
        annulusIdx = annulusStartIdx + n
        n2Idx = n2StartIdx - n
        normD2 = vector.normalise(o1_d2[-1][annulusIdx])
        d = [sf * normD2[c] for c in range(3)]
        d1Outer[n2Idx][elementsAroundHalfDuod + (0 if completeRowsBeforeEso + 2 < n2Idx < n2StartIdx else 1)] = d

    # Deal with curvature for these points!!!

    # Calculate d3
    d3UnitOuter = []
    for n2 in range(len(xOuter)):
        d3Around = []
        for n1 in range(len(xOuter[n2])):
            d3Around.append(vector.normalise(
                vector.crossproduct3(vector.normalise(d1Outer[n2][n1]), vector.normalise(d2Outer[n2][n1]))))
        d3UnitOuter.append(d3Around)

    # for n2 in range(len(xOuter)):
    #     for n1 in range(len(xOuter[n2])):
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xOuter[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Outer[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Outer[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3UnitOuter[n2][n1])
    #         nodeIdentifier += 1

    # Calculate curvature around
    d1Curvature = []
    d1Curvature.append([1.0 for n in range(len(xOuter[0]))]) # To be replaced later!
    d1Curvature[0][0] = d2Curvature[int(0.25 * len(xOuter[1]))][0]

    for n2 in range(1, len(xOuter)):
        if 2 + completeRowsBeforeEso < n2 < 2 + completeRowsBeforeEso + elementsAroundHalfEso - 2:
            curvatureLeft = findCurvatureAlongLine(xOuter[n2][0: int(0.5 * len(xOuter[n2])) + 1],
                                                   d1OuterForCurvature[n2][0: int(0.5 * len(xOuter[n2])) + 1],
                                                   d3UnitOuter[n2][0: int(0.5 * len(xOuter[n2])) + 1])
            xRight = xOuter[n2][int(0.5 * len(xOuter[n2])) + 1:] + [xOuter[n2][0]]
            d1Right = d1OuterForCurvature[n2][int(0.5 * len(xOuter[n2])) + 1:] + [d1OuterForCurvature[n2][0]]
            d3Right = d3UnitOuter[n2][int(0.5 * len(xOuter[n2])) + 1:] + [d3UnitOuter[n2][0]]
            curvatureRight = findCurvatureAlongLine(xRight, d1Right, d3Right)
            curvatureBoth = curvatureLeft + curvatureRight[:-1]
            d1Curvature.append(curvatureBoth)
        else:
            d1Curvature.append(findD1CurvatureAround(xOuter[n2], d1OuterForCurvature[n2], d3UnitOuter[n2]))

    # Calculate curvature along
    # Populate d3AlongAll
    d3AlongAll = []
    count = 0
    for n1 in range(elementsCountAroundDuod):
        d3AlongN1 = []
        if n1 == 0: # GC
            d3AlongN1.append(d3UnitOuter[0][0])
            for n in range(1, len(xOuter[0]), 2):
                d3AlongN1.append(d3UnitOuter[0][n])
            for n2 in range(1, len(xOuter)):
                d3AlongN1.append(d3UnitOuter[n2][0])

        elif n1 == 1: # Start from 2nd point on quarterline
            for n in range(int(0.25 * len(xOuter[1])), 0, -1):
                d3AlongN1.append(d3UnitOuter[1][n])
            for n2 in range(2, elementsAlongStomach + 1): # Template
                d3AlongN1.append(d3UnitOuter[n2][n1])

        elif n1 == 2: # Start from xOuter[1][0]
            d3AlongN1 += d3UnitOuter[1][:2]
            for n2 in range(2, elementsAlongStomach + 1): # Template
                d3AlongN1.append(d3UnitOuter[n2][n1])

        elif 2 < n1 < elementsAroundQuarterDuod or elementsAroundQuarterDuod < n1 < elementsAroundHalfDuod - 2:
            # Additional elements for duodenum above and below quarter line
            if 2 < n1 < elementsAroundQuarterDuod:
                x0Idx = int(2 * int(0.5 * len(xOuter[0]))) - 1 - 2 * count
            else:
                x0Idx = -(int(2 * int(0.5 * len(xOuter[0]))) - 2 * count)
            d3AlongN1 = [d3UnitOuter[0][x0Idx]] + [d3UnitOuter[1][2 + count]]

            for n2 in range(2, elementsAlongStomach + 1): # Template
                d3AlongN1.append(d3UnitOuter[n2][n1])
            count += 1

        elif n1 == elementsAroundQuarterDuod:
            d3AlongN1.append(d3UnitOuter[0][0])
            for n2 in range(1, len(xOuter)):
                d3AlongN1.append(d3UnitOuter[n2]
                                 [int(0.25 * len(xOuter[n2])) +
                                  (1 if 2 + completeRowsBeforeEso < n2 < completeRowsBeforeEso + elementsAroundHalfEso
                                   else 0)])
            count += 1

        elif n1 == elementsAroundHalfDuod - 2: # start from xOuter[1][half]
            d3AlongN1 = [d3UnitOuter[1][int(len(xOuter[1]) * 0.5)]] + [d3UnitOuter[1][int(len(xOuter[1]) * 0.5 - 1)]]
            for n2 in range(2, elementsAlongStomach + 1): # Template
                d3AlongN1.append(d3UnitOuter[n2][n1])

        elif n1 == elementsAroundHalfDuod - 1:
            # Start from 2nd point on quarterline
            for n in range(int(0.25 * len(xOuter[1])), int(len(xOuter[1]) * 0.5)):
                d3AlongN1.append(d3UnitOuter[1][n])
            for n2 in range(2, elementsAlongStomach + 1):
                d3AlongN1.append(d3UnitOuter[n2][n1])

        elif n1 == elementsAroundHalfDuod:
            d3AlongN1.append(d3UnitOuter[0][0])
            for n in range(2, len(xOuter[0]), 2):
                d3AlongN1.append(d3UnitOuter[0][n])
            for n2 in range(1, 3 + completeRowsBeforeEso):
                d3AlongN1.append(d3UnitOuter[n2][int(0.5 * len(xOuter[n2]))])
            # Need to be processed as two parts!
            breakIdx = len(d3AlongN1)
            for n2 in range(completeRowsBeforeEso + elementsAroundHalfEso, len(xOuter)):
                d3AlongN1.append(d3UnitOuter[n2][int(0.5 * len(xOuter[n2]))])

        elif n1 == elementsAroundHalfDuod + 1:
            # Start from 2nd point on threeQuarterline
            for n in range(int(0.75 * len(xOuter[1])), int(len(xOuter[1]) * 0.5), -1):
                d3AlongN1.append(d3UnitOuter[1][n])
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                d3AlongN1.append(d3UnitOuter[n2][n1Idx])

        elif n1 == elementsAroundHalfDuod + 2:
            # Start from xOuter[1][half]
            d3AlongN1 += [d3UnitOuter[1][int(len(xOuter[1]) * 0.5)]] + [d3UnitOuter[1][int(len(xOuter[1]) * 0.5) + 1]]
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                d3AlongN1.append(d3UnitOuter[n2][n1Idx])

        elif elementsAroundHalfDuod + 2 < n1 < int(elementsAroundQuarterDuod * 3) or \
                int(elementsAroundQuarterDuod * 3) < n1 < elementsCountAroundDuod - 2:
            count -= 1
            # Additional elements for duodenum above and below three quarter line
            if elementsAroundHalfDuod + 2 < n1 < int(elementsAroundQuarterDuod * 3):
                x0Idx = -(int(2 * int(0.5 * len(xOuter[0]))) - 2 * count)
            else:
                x0Idx = int(2 * int(0.5 * len(xOuter[0]))) - 1 - 2 * count
            d3AlongN1 = [d3UnitOuter[0][x0Idx]] + [d3UnitOuter[1][int(len(xOuter[1]) * 0.5) + len(x0) - count + 1]]
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                d3AlongN1.append(d3UnitOuter[n2][n1Idx])

        elif n1 == int(elementsAroundQuarterDuod * 3):
            d3AlongN1.append(d3UnitOuter[0][0])
            for n2 in range(1, len(xOuter)):
                d3AlongN1.append(d3UnitOuter[n2][int(0.75 * len(xOuter[n2]))])
            count -= 1

        elif n1 == elementsCountAroundDuod - 2: # Start from xOuter[1][0]
            d3AlongN1 += [d3UnitOuter[1][0]] + [d3UnitOuter[1][-1]]
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                d3AlongN1.append(d3UnitOuter[n2][n1Idx])

        elif n1 == elementsCountAroundDuod - 1: # Start from 2nd point on quarterline
            for n in range(int(0.25 * len(xOuter[1]))):
                d3AlongN1.append(d3UnitOuter[1][int(0.75 * len(xOuter[1])) + n])
            for n2 in range(2, elementsAlongStomach + 1): # new condition
                if 2 + len(xAroundB) < n2 < 2 + len(xAroundB) + elementsAroundHalfEso - 2:
                    n1Idx = n1 - 1
                else:
                    n1Idx = n1
                d3AlongN1.append(d3UnitOuter[n2][n1Idx])

        d3AlongAll.append(d3AlongN1)

#     # for n1 in range(len(xAlongAll)):
#     #     for n2 in range(len(xAlongAll[n1])):
#     #         node = nodes.createNode(nodeIdentifier, nodetemplate)
#     #         cache.setNode(node)
#     #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAlongAll[n1][n2])
#     #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
#     #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2AlongAll[n1][n2])
#     #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
#     #         nodeIdentifier += 1

    # Calculate d2Curvature along
    d2CurvatureAlong = []
    for n1 in range(len(xAlongAll)):
        if n1 == elementsAroundHalfDuod:
            curvatureGC = findCurvatureAlongLine(xAlongAll[n1][:breakIdx], d2AlongAll[n1][:breakIdx],
                                                 d3AlongAll[n1][:breakIdx])
            curvatureLC = findCurvatureAlongLine(xAlongAll[n1][breakIdx:], d2AlongAll[n1][breakIdx:],
                                                 d3AlongAll[n1][breakIdx:])
            curvature = curvatureGC + curvatureLC
        else:
            curvature = findCurvatureAlongLine(xAlongAll[n1], d2AlongAll[n1], d3AlongAll[n1])
        d2CurvatureAlong.append(curvature)

    # Replace curvature on annulus
    annulusStartIdx = 1
    n2StartIdx = completeRowsBeforeEso + 2
    for n in range(elementsAroundHalfEso - 1):
        annulusIdx = annulusStartIdx + n
        n2Idx = n2StartIdx + n + int(0.5 * (len(xOuter[0]) - 1))
        d2CurvatureAlong[elementsAroundHalfDuod - 1][n2Idx] = annulusCurvature[annulusIdx]

    annulusStartIdx = elementsAroundHalfEso + 1
    n2StartIdx = completeRowsBeforeEso + elementsAroundHalfEso
    for n in range(elementsAroundHalfEso - 1):
        annulusIdx = annulusStartIdx + n
        n2Idx = n2StartIdx - n + int(0.5 * (len(xOuter[0]) - 1))
        d2CurvatureAlong[elementsAroundHalfDuod + 1][n2Idx] = annulusCurvature[annulusIdx]

    # Re-arrange d2Curvature along to xOuter
    for n2 in range(len(xOuter)):
        if n2 == 0:
            n1 = 0
            d2Curvature[n2][n1] = d2CurvatureAlong[n1][n2]
            if len(d2Curvature[n2]) > 1:
                for n in range(int(len(xOuter[n2]) * 0.5)):
                    d2Curvature[n2][2 * n + 1] = d2CurvatureAlong[elementsAroundQuarterDuod - n - 1][n2]
                    d2Curvature[n2][2 * n + 2] = d2CurvatureAlong[elementsAroundQuarterDuod + n + 1][n2]
        elif n2 == 1:
            d2New = []
            count = 0
            for n in range(int(2 * (len(xOuter[0]) + 2))):
                n1AllIdx = elementsAroundQuarterDuod - int(0.5 * len(xOuter[0])) - 1 + n + count
                if n == 0:
                    d2New.append(d2CurvatureAlong[n1AllIdx][0])
                d2New.append(d2CurvatureAlong[n1AllIdx][1])
                if n == len(xOuter[0]) + 1:
                    d2New.append(d2CurvatureAlong[n1AllIdx][0])
                    count = 3
            d2Curvature[n2] = d2New

        elif n2 > 1:
            d2New = []
            idxFromEnd = -(elementsAlongStomach + 1 - n2)
            for n1 in range(elementsCountAroundDuod +
                            (-1 if 3 + completeRowsBeforeEso <= n2 <
                                   3 + completeRowsBeforeEso + (elementsAroundHalfEso - 3) else 0)):
                if n1 == elementsAroundHalfDuod and 1 < n2 < 3 + completeRowsBeforeEso:
                    d2New.append(d2CurvatureAlong[n1][int(0.5 * len(xOuter[0])) + n2])
                elif n1 >= elementsAroundHalfDuod \
                        and 3 + completeRowsBeforeEso <= n2 < 3 + completeRowsBeforeEso + (elementsAroundHalfEso - 3):
                    d2New.append(d2CurvatureAlong[n1 + 1][idxFromEnd])
                else:
                    d2New.append(d2CurvatureAlong[n1][idxFromEnd])
            d2Curvature[n2] = d2New

    # Replace curvature with change of derivatives on xOuter[0] (apart from apex)
    d3AlongGC = []
    for n2 in range(2 + completeRowsBeforeEso, 0, -1):
        d3AlongGC.append(d3UnitOuter[n2][int(0.5 * len(xOuter[n2]))])
    if len(xOuter) > 1:
        for n in range(len(xOuter[0]) - 1, 0, -2):
            d3AlongGC.append(d3UnitOuter[0][n])

    for n2 in range(len(xAlongAll[0])):
        d3AlongGC.append(d3AlongAll[0][n2])
    curvatureGC = findCurvatureAlongLine(xAlongGC, d2AlongGCSmoothed, d3AlongGC)

    count = 0
    for n in range(1, len(xOuter[0]), 2):
        d1Curvature[0][n] = curvatureGC[-(elementsAlongFundusApexToCardia + int(0.5 * n) + 2)]
        d1Curvature[0][n + 1] = curvatureGC[-(elementsAlongFundusApexToCardia - int(0.5 * n))]

        n1IdxAlong = int((len(xOuter[0]) - 1) * 0.5) - count
        d2Curvature[0][n] = d2CurvatureAlong[int(3 * elementsAroundQuarterDuod) - n1IdxAlong][0]
        d2Curvature[0][n + 1] = d2CurvatureAlong[int(3 * elementsAroundQuarterDuod) + n1IdxAlong][0]
        count += 1

    # Set d1 on xOuter[1] on GC to point along GC towards esophagus
    if len(xOuter[0]) == 1:
        d1Curvature[1][0] = curvatureGC[-(elementsAlongFundusApexToCardia + 2)]
        d1Curvature[1][int(0.5 * len(xOuter[1]))] = curvatureGC[-(elementsAlongFundusApexToCardia)]
    else:
        d1Curvature[1][0] = curvatureGC[-(elementsAlongFundusApexToCardia + int(0.5 * n) + 3)]
        d1Curvature[1][int(0.5 * len(xOuter[1]))] = \
            curvatureGC[-(elementsAlongFundusApexToCardia - int(0.5 * n) - 1)]

    # Replace d1 at tip of wedges
    n2 = 1
    n1IdxOuterAll = [1, int(0.5 * len(xOuter[1])) - 1, int(0.5 * len(xOuter[1])) + 1, len(xOuter[1]) - 1]
    n1IdxAlongAll = [1, elementsAroundHalfDuod - 1, elementsAroundHalfDuod + 1, elementsCountAroundDuod - 1]
    for n1 in range(4):
        n1IdxOuter = n1IdxOuterAll[n1]
        n1IdxAlong = n1IdxAlongAll[n1]
        d = d2AlongAll[n1IdxAlong][1 + int(0.5 * (len(xOuter[0]) - 1))]
        if n1 == 0 or n1 == 2:
            rotAxis = d3AlongAll[n1IdxAlong][1 + int(0.5 * (len(xOuter[0]) - 1))]
            rotFrame = matrix.getRotationMatrixFromAxisAngle(rotAxis, math.pi)
            d = [rotFrame[j][0] * d[0] + rotFrame[j][1] * d[1] + rotFrame[j][2] * d[2] for j in range(3)]
        d1Outer[n2][n1IdxOuter] = d
        d1Curvature[n2][n1IdxOuter] = d2CurvatureAlong[n1IdxAlong][1 + int(0.5 * (len(xOuter[0]) - 1))]

    # Create inner nodes
    xList = []
    d1List = []
    d2List = []
    d3List = []
    nodeIdx = stomachStartNode
    idxMat = []

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

    for n2 in range(len(xOuter)):
        idxThroughWall = []
        for n3 in range(elementsCountThroughWall + 1):
            xi3 = xi3List[n3] if elementsCountThroughWall > 1 else 1.0 / elementsCountThroughWall * n3
            idxAround = []
            for n1 in range(len(xOuter[n2])):
                # Coordinates
                norm = d3UnitOuter[n2][n1]
                xOut = xOuter[n2][n1]
                xIn = [xOut[i] - norm[i] * wallThickness for i in range(3)]
                dWall = [wallThickness * c for c in norm]
                x = interp.interpolateCubicHermite(xIn, dWall, xOut, dWall, xi3)
                xList.append(x)

                # d1
                factor = 1.0 + wallThickness * (1.0 - xi3) * d1Curvature[n2][n1]
                d1 = [factor * c for c in d1Outer[n2][n1]]
                d1List.append(d1)

                # d2
                factor = 1.0 + wallThickness * (1.0 - xi3) * d2Curvature[n2][n1]
                d2 = [factor * c for c in d2Outer[n2][n1]]
                d2List.append(d2)

                # d3
                d3 = [c * wallThickness * (thicknessProportions[n3 + 1] if elementsCountThroughWall > 1 else 1.0)
                      for c in norm]
                d3List.append(d3)

                idxAround.append(nodeIdx)
                nodeIdx += 1
            idxThroughWall.append(idxAround)
        idxMat.append(idxThroughWall)

    nodeIdxGC = []
    for n2 in range(len(idxMat)):
        for n3 in range(len(idxMat[n2])):
            if n2 == 0:
                nodeIdxGC += idxMat[n2][n3]
            else:
                nodeIdxGC.append(idxMat[n2][n3][0])
    for n2 in range(1, completeRowsBeforeEso + 3):
        for n3 in range(len(idxMat[n2])):
            nodeIdxGC.append(idxMat[n2][n3][int(0.5 * len(xOuter[n2]))])

    nodeIdxLC = []
    for n2 in range(elementsAlongCardiaToDuod + 1):
        for n3 in range(len(idxMat[n2])):
            nodeIdxLC.append(idxMat[completeRowsBeforeEso + elementsAroundHalfEso + n2][n3][elementsAroundHalfDuod])

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

#     # Move along elements on GC and assign to annotation group
#     elementsCountAlongGroups = []
#     groupLength = 0.0
#     e = 0
#     elementsCount = 1
#     length = arcLengthsGCFromFundusApexToDuod[e]
#     for i in range(len(arcLengthRatioForGroupsFromFundusApex)):
#         groupLength += arcLengthRatioForGroupsFromFundusApex[i] * arcLengthGCFromFundusApexToDuodTotal
#         if e == len(arcLengthsGCFromFundusApexToDuod) - 2:
#             elementsCount += 1
#             elementsCountAlongGroups.append(elementsCount)
#         else:
#             while length < groupLength:
#                 elementsCount += 1
#                 e += 1
#                 length += arcLengthsGCFromFundusApexToDuod[e]
#
#             # check which end is grouplength closer to
#             distToUpperEnd = abs(length - groupLength)
#             distToLowerEnd = abs(groupLength - (length - arcLengthsGCFromFundusApexToDuod[e]))
#             if distToLowerEnd < distToUpperEnd:
#                 elementsCount -= 1
#                 elementsCountAlongGroups.append(elementsCount)
#                 e -= 1
#                 length -= arcLengthsGCFromFundusApexToDuod[e]
#             else:
#                 elementsCountAlongGroups.append(elementsCount)
#         elementsCount = 0
#
#     #  Account for multiple elements in n2 == 0:
#     if len(xOuter[0]) > 1:
#         elementsCountAlongGroups[0] -= int((len(xOuter[0]) - 1) * 0.5)
#
    annotationGroupsAlong = []
    for i in range(len(elementCountGroupList)):
        elementsCount = elementCountGroupList[i]
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

    idxBottom = []
    idxTop = []
    for n in range(len(xOuter[0]) - 2, -1, -2):
        idxBottom.append(n)
    for n in range(0, len(xOuter[0]), 2):
        idxTop.append(n)
    outer0 = idxBottom + idxTop

    fundusElements = elementCountGroupList[0]

    for e2 in range(len(xOuter) - 1):
        elementIdxThroughWall = []
        if e2 == 0:
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(len(xOuter[e2 + 1]) - 4):
                    eft1 = eftfactory.createEftNoCrossDerivatives()
                    scaleFactors = [-1.0]
                    setEftScaleFactorIds(eft1, [1], [])

                    if e1 <= len(xOuter[0]):
                        bni1 = (idxMat[e2][e3][outer0[e1]] if e1 < len(xOuter[0]) else idxMat[e2 + 1][e3][3 + e1])
                        bni2 = idxMat[e2 + 1][e3][2 + (0 if e1 == 0 else e1)]
                        bni3 = (idxMat[e2 + 1][e3][0] if e1 == 0 else idxMat[e2][e3][outer0[e1 - 1]])
                        bni4 = idxMat[e2 + 1][e3][1 + (0 if e1 == 0 else e1)]
                        bni5 = (idxMat[e2][e3 + 1][outer0[e1]] if e1 < len(xOuter[0]) else
                                idxMat[e2 + 1][e3 + 1][3 + e1])
                        bni6 = idxMat[e2 + 1][e3 + 1][2 + (0 if e1 == 0 else e1)]
                        bni7 = (idxMat[e2 + 1][e3 + 1][0] if e1 == 0 else idxMat[e2][e3 + 1][outer0[e1 - 1]])
                        bni8 = idxMat[e2 + 1][e3 + 1][1 + (0 if e1 == 0 else e1)]

                        remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS2, [])])

                    else:
                        bni1 = idxMat[e2 + 1][e3][3 + e1]
                        bni2 = (idxMat[e2 + 1][e3][2 + e1] if e1 == len(xOuter[0]) + 1 else
                                idxMat[e2][e3][outer0[len(xOuter[0]) - e1 + 1]])
                        bni3 = idxMat[e2 + 1][e3][4 + e1]
                        bni4 = (idxMat[e2][e3][outer0[len(xOuter[0]) - e1]] if e1 < len(xOuter[e2 + 1]) - 5 else
                                idxMat[e2 + 1][e3][0])
                        bni5 = idxMat[e2 + 1][e3 + 1][3 + e1]
                        bni6 = (idxMat[e2 + 1][e3 + 1][2 + e1] if e1 == len(xOuter[0]) + 1 else
                                idxMat[e2][e3 + 1][outer0[len(xOuter[0]) - e1 + 1]])
                        bni7 = idxMat[e2 + 1][e3 + 1][4 + e1]
                        bni8 = (idxMat[e2][e3 + 1][outer0[len(xOuter[0]) - e1]] if e1 < len(xOuter[e2 + 1]) - 5 else
                                idxMat[e2 + 1][e3 + 1][0])

                        remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS1,
                                               [(Node.VALUE_LABEL_D_DS2, [1])])
                        remapEftNodeValueLabel(eft1, [2, 4, 6, 8], Node.VALUE_LABEL_D_DS2,
                                               [(Node.VALUE_LABEL_D_DS1, [])])

                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS1,
                                           [(Node.VALUE_LABEL_D_DS2, [1])])
                    remapEftNodeValueLabel(eft1, [1, 3, 5, 7], Node.VALUE_LABEL_D_DS2,
                                           [(Node.VALUE_LABEL_D_DS1, [])])
                    nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplateX)
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

        elif e2 == 1:
            for e3 in range(elementsCountThroughWall):
                offset = 0
                elementIdxAround = []
                for e1 in range(len(xOuter[e2 + 1])):
                    scaleFactors = []
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard
                    if e1 in [1, int(0.5 * len(xOuter[e2 + 1]) - 2), int(0.5 * len(xOuter[e2 + 1]) + 1),
                              len(xOuter[e2 + 1]) - 2]:
                        offset -= 1
                    bni1 = idxMat[e2][e3][e1 + offset]
                    bni2 = idxMat[e2][e3][(e1 + 1 + offset) % len(idxMat[e2][e3])]
                    bni3 = idxMat[e2 + 1][e3][e1]
                    bni4 = idxMat[e2 + 1][e3][(e1 + 1) % len(idxMat[e2 + 1][e3])]
                    bni5 = idxMat[e2][e3 + 1][e1 + offset]
                    bni6 = idxMat[e2][e3 + 1][(e1 + 1 + offset) % len(idxMat[e2][e3])]
                    bni7 = idxMat[e2 + 1][e3 + 1][e1]
                    bni8 = idxMat[e2 + 1][e3 + 1][(e1 + 1) % len(idxMat[e2 + 1][e3])]
                    nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]

                    if e1 in [0, int(0.5 * len(xOuter[e2 + 1]) - 1)]:
                        scaleFactors = [-1.0]
                        eft1 = eftfactory.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                        remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

                    elif e1 in [1, int(0.5 * len(xOuter[e2 + 1])) + 1]:
                        nodeIdentifiers = [bni2, bni3, bni4, bni6, bni7, bni8]
                        eft1 = eftfactory.createEftWedgeCollapseXi1Quadrant([1, 5])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX
                        setEftScaleFactorIds(eft1, [1], [])

                    elif e1 in [int(0.5 * len(xOuter[e2 + 1])) - 2, len(xOuter[e2 + 1]) - 2]:
                        nodeIdentifiers = [bni2, bni3, bni4, bni6, bni7, bni8]
                        eft1 = eftfactory.createEftWedgeCollapseXi1Quadrant([2, 6])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

                    elif e1 == int(0.5 * len(xOuter[e2 + 1])):
                        scaleFactors = [-1.0]
                        eft1 = eftfactory.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                        remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1])])
                        remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

                    elif e1 == len(xOuter[e2 + 1]) - 1:
                        scaleFactors = [-1.0]
                        eft1 = eftfactory.createEftNoCrossDerivatives()
                        setEftScaleFactorIds(eft1, [1], [])
                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                        remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                        remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                        remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                        elementtemplateX.defineField(coordinates, -1, eft1)
                        elementtemplate1 = elementtemplateX

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

        elif e2 == completeRowsBeforeEso + 2: # first ring in esophagus
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(len(xOuter[e2]) - 1):
                    scaleFactors = []
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard
                    if e1 == int(0.5 * len(xOuter[e2]) - 1):
                        pass
                    else:
                        if e1 > int(0.5 * len(xOuter[e2]) - 1):
                            e1IdxBni1 = e1 + 1
                        else:
                            e1IdxBni1 = e1
                        bni1 = idxMat[e2][e3][e1IdxBni1]
                        bni2 = idxMat[e2][e3][(e1IdxBni1 + 1) % len(idxMat[e2][e3])]
                        bni3 = idxMat[e2 + 1][e3][e1]
                        bni4 = idxMat[e2 + 1][e3][(e1 + 1) % len(idxMat[e2 + 1][e3])]
                        bni5 = idxMat[e2][e3 + 1][e1IdxBni1]
                        bni6 = idxMat[e2][e3 + 1][(e1IdxBni1 + 1) % len(idxMat[e2][e3])]
                        bni7 = idxMat[e2 + 1][e3 + 1][e1]
                        bni8 = idxMat[e2 + 1][e3 + 1][(e1 + 1) % len(idxMat[e2 + 1][e3])]
                        nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
                        if e1 == int(0.5 * len(xOuter[e2])) - 2:
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            # remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []),
                                                                                          (Node.VALUE_LABEL_D_DS2, [1]) ])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                        elif e1 == int(0.5 * len(xOuter[e2])):
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []),
                                                                                           (Node.VALUE_LABEL_D_DS2, []) ])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
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

        elif completeRowsBeforeEso + 2 < e2 < completeRowsBeforeEso + elementsAroundHalfEso - 1: # within esophagus
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(len(xOuter[e2]) - 1):
                    scaleFactors = []
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard
                    if e1 > int(0.5 * len(xOuter[e2]) - 1):
                        e1IdxBni1 = e1 + 1
                    else:
                        e1IdxBni1 = e1
                    bni1 = idxMat[e2][e3][e1IdxBni1]
                    bni2 = idxMat[e2][e3][(e1IdxBni1 + 1) % len(idxMat[e2][e3])]
                    bni3 = idxMat[e2 + 1][e3][e1IdxBni1]
                    bni4 = idxMat[e2 + 1][e3][(e1IdxBni1 + 1) % len(idxMat[e2 + 1][e3])]
                    bni5 = idxMat[e2][e3 + 1][e1IdxBni1]
                    bni6 = idxMat[e2][e3 + 1][(e1IdxBni1 + 1) % len(idxMat[e2][e3])]
                    bni7 = idxMat[e2 + 1][e3 + 1][e1IdxBni1]
                    bni8 = idxMat[e2 + 1][e3 + 1][(e1IdxBni1 + 1) % len(idxMat[e2 + 1][e3])]
                    nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
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

        elif e2 == completeRowsBeforeEso + elementsAroundHalfEso - 1: # last ring in esophagus
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(len(xOuter[e2])):
                    scaleFactors = []
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard
                    if e1 == int(0.5 * len(xOuter[e2])):
                        pass
                    else:
                        if e1 > int(0.5 * len(xOuter[e2])):
                            e1IdxBni3 = e1 + 1
                        else:
                            e1IdxBni3 = e1
                        bni1 = idxMat[e2][e3][e1]
                        bni2 = idxMat[e2][e3][(e1 + 1) % len(idxMat[e2][e3])]
                        bni3 = idxMat[e2 + 1][e3][e1IdxBni3]
                        bni4 = idxMat[e2 + 1][e3][(e1IdxBni3 + 1) % len(idxMat[e2 + 1][e3])]
                        bni5 = idxMat[e2][e3 + 1][e1]
                        bni6 = idxMat[e2][e3 + 1][(e1 + 1) % len(idxMat[e2][e3])]
                        bni7 = idxMat[e2 + 1][e3 + 1][e1IdxBni3]
                        bni8 = idxMat[e2 + 1][e3 + 1][(e1IdxBni3 + 1) % len(idxMat[e2 + 1][e3])]
                        nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6, bni7, bni8]
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

        elif 1 < e2 < 2 + completeRowsBeforeEso or e2 > completeRowsBeforeEso + elementsAroundHalfEso - 1:
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(len(xOuter[e2])):
                    scaleFactors = []
                    eft1 = eftStandard
                    elementtemplate1 = elementtemplateStandard
                    bni111 = idxMat[e2][e3][e1]
                    bni211 = idxMat[e2][e3][(e1 + 1) % len(idxMat[e2][e3])]
                    bni121 = idxMat[e2 + 1][e3][e1]
                    bni221 = idxMat[e2 + 1][e3][(e1 + 1) % len(idxMat[e2 + 1][e3])]
                    bni112 = idxMat[e2][e3 + 1][e1]
                    bni212 = idxMat[e2][e3 + 1][(e1 + 1) % len(idxMat[e2][e3])]
                    bni122 = idxMat[e2 + 1][e3 + 1][e1]
                    bni222 = idxMat[e2 + 1][e3 + 1][(e1 + 1) % len(idxMat[e2 + 1][e3])]
                    nodeIdentifiers = [bni111, bni211, bni121, bni221,
                                       bni112, bni212, bni122, bni222]

                    if e2 == completeRowsBeforeEso + 1:
                        if e1 == int(0.5 * len(xOuter[e2]) - 2):
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                        elif e1 == int(0.5 * len(xOuter[e2]) - 1):
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                        elif e1 == int(0.5 * len(xOuter[e2])):
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                        elif e1 == int(0.5 * len(xOuter[e2]) + 1):
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, []),
                                                                                          (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                    if e2 == completeRowsBeforeEso + elementsAroundHalfEso:
                        if e1 == int(0.5 * len(xOuter[e2]) - 2):
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2, [(Node.VALUE_LABEL_D_DS1, [1]),
                                                                                         (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e1 == int(0.5 * len(xOuter[e2]) - 1):
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,[(Node.VALUE_LABEL_D_DS1, [1]),
                                                                                         (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                        elif e1 == int(0.5 * len(xOuter[e2])):
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,[(Node.VALUE_LABEL_D_DS1, []),
                                                                                         (Node.VALUE_LABEL_D_DS2, [])])

                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX
                        elif e1 == int(0.5 * len(xOuter[e2]) + 1):
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,[(Node.VALUE_LABEL_D_DS1, []),
                                                                                         (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

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

    n2IdxDownRing = completeRowsBeforeEso + elementsAroundHalfEso
    for n3 in range(elementsCountThroughWall + 1):
        n1 = 0
        for nAround in range(elementsCountAroundEso):
            if nAround == 0:
                idx = idxMat[2 + completeRowsBeforeEso][n3][int(0.5 * len(xOuter[2 + completeRowsBeforeEso]))]
            elif 0 < nAround < elementsAroundHalfEso:
                n2Idx = 2 + completeRowsBeforeEso + n1
                n1Idx = int((len(xOuter[n2Idx]) - 1) * 0.5)
                idx = idxMat[n2Idx][n3][n1Idx]
                n1 += 1
            elif nAround == elementsAroundHalfEso:
                idx = idxMat[n2IdxDownRing][n3][int(0.5 * len(xOuter[n2IdxDownRing]))]
            else:
                n2Idx = n2IdxDownRing + (elementsAroundHalfEso + 1 - nAround)
                n1Idx = int(len(xOuter[n2Idx]) * 0.5 + 1)
                idx = idxMat[n2Idx][n3][n1Idx]

            endPoints_x[n3][nAround] = xList[idx - stomachStartNode]
            endPoints_d1[n3][nAround] = d1List[idx - stomachStartNode]
            endPoints_d2[n3][nAround] = d2List[idx - stomachStartNode]
            endNode_Id[n3][nAround] = idx

            if n3 == elementsCountThroughWall:  # outer layer
                endPosition = trackSurfaceStomach.findNearestPosition(endPoints_x[n3][nAround])
                endProportions.append(trackSurfaceStomach.getProportion(endPosition))

    for n3 in range(elementsCountThroughWall + 1):
        for nAround in range(elementsCountAroundEso):
            if nAround == 0:
                endDerivativesMap[n3][nAround] = ((-1, 0, 0), (0, -1, 0), None)
            elif 0 < nAround <= elementsAroundHalfEso - 1:
                endDerivativesMap[n3][nAround] = ((0, 1, 0), (-1, 0, 0), None)
            elif nAround == elementsAroundHalfEso:
                endDerivativesMap[n3][nAround] = (None, None, None)
            elif elementsAroundHalfEso + 1 <= nAround <= elementsCountAroundEso - 1:
                endDerivativesMap[n3][nAround] = ((0, -1, 0), (1, 0, 0), None)

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
        startProportions=startProportions, endProportions=endProportions, coordinates=coordinates)

    elementIdxThroughWall = []
    n = lastDuodenumElementIdentifier - 1
    for n3 in range(elementsCountThroughWall):
        elementIdxAround = []
        for n1 in range(elementsCountAroundEso):
            n += 1
            elementIdxAround.append(n)
        elementIdxThroughWall.append(elementIdxAround)
    elementIdxMat.append(elementIdxThroughWall)

    nodeIdentifier = nextNodeIdentifier

    # delete mucosa layer in fundus when there is a limiting ridge
    mesh_destroy_elements_and_nodes_by_identifiers(mesh, fundusMucosaElementIdentifiers)

#     # # annotation fiducial points for embedding in whole body
#     # markerNames = [["esophagogastric junction along the greater curvature on luminal surface",
#     #                 "esophagogastric junction along the lesser curvature on luminal surface",
#     #                 "gastroduodenal junction along the greater curvature on luminal surface",
#     #                 "gastroduodenal junction along the lesser curvature on luminal surface",
#     #                 "body-antrum junction along the greater curvature on luminal surface",
#     #                 "limiting ridge at the greater curvature on the luminal surface" if limitingRidge else
#     #                 "fundus-body junction along the greater curvature on luminal surface"]]
#     # if elementsCountThroughWall == 4:
#     #     markerNames.append(
#     #         ["esophagogastric junction along the greater curvature on circular-longitudinal muscle interface",
#     #          "esophagogastric junction along the lesser curvature on circular-longitudinal muscle interface",
#     #          "gastroduodenal junction along the greater curvature on circular-longitudinal muscle interface",
#     #          "gastroduodenal junction along the lesser curvature on circular-longitudinal muscle interface",
#     #          "body-antrum junction along the greater curvature on circular-longitudinal muscle interface",
#     #          "limiting ridge at the greater curvature on the circular-longitudinal muscle interface" if limitingRidge
#     #          else "fundus-body junction along the greater curvature on circular-longitudinal muscle interface"])
#     # markerNames.append(["esophagogastric junction along the greater curvature on serosa",
#     #                     "esophagogastric junction along the lesser curvature on serosa",
#     #                     "gastroduodenal junction along the greater curvature on serosa",
#     #                     "gastroduodenal junction along the lesser curvature on serosa",
#     #                     "body-antrum junction along the greater curvature on serosa",
#     #                     "limiting ridge at the greater curvature on serosa" if limitingRidge else
#     #                     "fundus-body junction along the greater curvature on serosa",
#     #                     "distal point of lower esophageal sphincter serosa on the greater curvature of stomach",
#     #                     "distal point of lower esophageal sphincter serosa on the lesser curvature of stomach"])
#     #
#     # markerInnerElementIdentifiers = [stomachStartElement - elementsCountThroughWall * elementsCountAroundEso,
#     #                                  stomachStartElement - (elementsCountThroughWall - 1) * elementsCountAroundEso -
#     #                                  elementsAroundHalfEso,
#     #                                  lastDuodenumElementIdentifier - elementsCountThroughWall *
#     #                                  elementsCountAroundDuod * (elementsCountAlongGroups[-1] + 1),
#     #                                  lastDuodenumElementIdentifier - elementsCountThroughWall *
#     #                                  elementsCountAroundDuod * (elementsCountAlongGroups[-1] + 1) +
#     #                                  elementsAroundHalfDuod,
#     #                                  lastDuodenumElementIdentifier - elementsCountThroughWall *
#     #                                  elementsCountAroundDuod * (sum(elementsCountAlongGroups[-3:]) + 1),
#     #                                  fundusBodyJunctionInnerElementIdentifier,
#     #                                  elementsCountAroundEso * (elementsCountThroughWall - 1) + 1,
#     #                                  elementsCountAroundEso * elementsCountThroughWall - elementsAroundHalfEso + 1]
#     #
#     # elementsCountAroundLayer = [elementsCountAroundEso, elementsCountAroundEso,
#     #                             elementsCountAroundDuod, elementsCountAroundDuod,
#     #                             elementsCountAroundDuod, elementsCountAroundDuod]
#     #
#     # for n3 in range(len(markerNames)):
#     #     for n in range(len(markerNames[n3])):
#     #         markerGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region,
#     #                                                          get_stomach_term(markerNames[n3][n]))
#     #         if n < 6:
#     #             markerElementIdentifier = \
#     #                 markerInnerElementIdentifiers[n] + \
#     #                 (0 if n3 == 0 or elementsCountThroughWall == 1 else elementsCountAroundLayer[n] *
#     #                                                                     (elementsCountThroughWall - 1))
#     #             markerXi = [0.0, 1.0, 0.0 if n3 != len(markerNames) - 1 else 1.0] if n < len(markerNames[n3]) - 1 \
#     #                 else [0.0, 0.0 if limitingRidge else 1.0, 0.0 if n3 != len(markerNames) - 1 else 1.0]
#     #         else:
#     #             markerElementIdentifier = markerInnerElementIdentifiers[n]
#     #             markerXi = [0.0, 0.0, 1.0]
#     #         markerElement = mesh.findElementByIdentifier(markerElementIdentifier)
#     #
#     #         cache.setMeshLocation(markerElement, markerXi)
#     #         markerPoint = markerPoints.createNode(nodeIdentifier, markerTemplateInternal)
#     #         nodeIdentifier += 1
#     #         cache.setNode(markerPoint)
#     #         markerName.assignString(cache, markerGroup.getName())
#     #         markerLocation.assignMeshLocation(cache, markerElement, markerXi)
#     #         for group in [stomachGroup, markerGroup]:
#     #             group.getNodesetGroup(nodes).addNode(markerPoint)
#     #
#     # Create annotation groups for dorsal and ventral parts of the stomach
#     dorsalGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("dorsal stomach"))
#     ventralGroup = findOrCreateAnnotationGroupForTerm(allAnnotationGroups, region, get_stomach_term("ventral stomach"))
#     dorsalMeshGroup = dorsalGroup.getMeshGroup(mesh)
#     ventralMeshGroup = ventralGroup.getMeshGroup(mesh)
#
#     for e2 in range(len(elementIdxMat)):
#         for e3 in range(len(elementIdxMat[e2])):
#             for e1 in range(len(elementIdxMat[e2][e3])):
#                 elementIdx = elementIdxMat[e2][e3][e1]
#                 element = mesh.findElementByIdentifier(elementIdx)
#                 if e1 < 0.5 * len(elementIdxMat[e2][e3]):
#                     ventralMeshGroup.addElement(element)
#                 else:
#                     dorsalMeshGroup.addElement(element)
#     if dorsalGroup not in allAnnotationGroups:
#         allAnnotationGroups.append(dorsalGroup)
#     if ventralGroup not in allAnnotationGroups:
#         allAnnotationGroups.append(ventralGroup)
#
#     # Create split coordinate field
#     if splitCoordinates:
#         nodesOnSplitMargin = []
#         nodesOnLCMargin = []
#
#         for n2 in range(elementsAlongEsophagus + 1):
#             for n3 in range(elementsThroughEsophagusWall + 1):
#                 nodeIdxOnGCMargin = 1 + n2 * (elementsThroughEsophagusWall + 1) * elementsCountAroundEso + \
#                                     n3 * elementsCountAroundEso
#                 nodesOnSplitMargin.append(nodeIdxOnGCMargin)
#                 nodeIdxOnLCMargin = 1 + elementsAroundHalfEso + \
#                                     n2 * (elementsThroughEsophagusWall + 1) * elementsCountAroundEso + \
#                                     n3 * elementsCountAroundEso
#                 nodesOnSplitMargin.append(nodeIdxOnLCMargin)
#                 nodesOnLCMargin.append(nodeIdxOnLCMargin)
#         nodesOnSplitMargin += nodeIdxGC + nodeIdxLC
#         allNodesOnLC = nodesOnLCMargin + nodeIdxLC
#
#         splitCoordinates = findOrCreateFieldCoordinates(fm, name="split coordinates")
#         splitNodetemplate1 = nodes.createNodetemplate()
#         splitNodetemplate1.defineField(splitCoordinates)
#         splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
#         splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
#         splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
#         if useCrossDerivatives:
#             splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
#         if useCubicHermiteThroughWall:
#             splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
#             if useCrossDerivatives:
#                 splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
#                 splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
#                 splitNodetemplate1.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
#
#         splitNodetemplate2 = nodes.createNodetemplate()
#         splitNodetemplate2.defineField(splitCoordinates)
#         splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_VALUE, 2)
#         splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
#         splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS2, 2)
#         if useCrossDerivatives:
#             splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 2)
#         if useCubicHermiteThroughWall:
#             splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D_DS3, 2)
#             if useCrossDerivatives:
#                 splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 2)
#                 splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 2)
#                 splitNodetemplate2.setValueNumberOfVersions(splitCoordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 2)
#
#         nodeIter = nodes.createNodeiterator()
#         node = nodeIter.next()
#         while node.isValid():
#             if not markerPoints.containsNode(node):
#                 cache.setNode(node)
#                 identifier = node.getIdentifier()
#                 marginNode = identifier in nodesOnSplitMargin
#                 x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, 3)[1]
#                 d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, 3)[1]
#                 d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, 3)[1]
#                 if useCrossDerivatives:
#                     d1d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, 3)[1]
#                 if useCubicHermiteThroughWall:
#                     d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, 3)[1]
#                     if useCrossDerivatives:
#                         d1d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, 3)[1]
#                         d2d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, 3)[1]
#                         d1d2d3 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, 3)[1]
#
#                 node.merge(splitNodetemplate2 if marginNode else splitNodetemplate1)
#                 versionCount = 2 if marginNode else 1
#                 for vn in range(versionCount):
#                     splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, vn + 1, x)
#                     splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, vn + 1, d1)
#                     splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, vn + 1, d2)
#                     if useCrossDerivatives:
#                         splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, vn + 1, d1d2)
#                     if useCubicHermiteThroughWall:
#                         splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, vn + 1, d3)
#                         if useCrossDerivatives:
#                             splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, vn + 1, d1d3)
#                             splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, vn + 1, d2d3)
#                             splitCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, vn + 1, d1d2d3)
#
#             node = nodeIter.next()
#
#         nearLCGroup = AnnotationGroup(region, ("elements adjacent to lesser curvature", "None"))
#
#         elementIter = mesh.createElementiterator()
#         element = elementIter.next()
#         splitElementtemplate1 = mesh.createElementtemplate()
#         splitElementtemplate2 = mesh.createElementtemplate()
#
#         count = 0
#         elementsInOstium = elementsCountAroundEso * elementsAlongEsophagus * elementsThroughEsophagusWall
#         closedLoopElementId = nextElementIdentifier - elementsCountAroundEso * elementsCountAcrossCardia - \
#                               elementsCountAroundDuod * elementsCountThroughWall * (elementsAlongCardiaToDuod + 1)
#
#         allValueLabels = [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2,
#                           Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3,
#                           Node.VALUE_LABEL_D2_DS2DS3, Node.VALUE_LABEL_D3_DS1DS2DS3]
#
#         while element.isValid():
#             eft = element.getElementfieldtemplate(coordinates, -1)
#             nodeIdentifiers = get_element_node_identifiers(element, eft)
#             for n in range(len(nodeIdentifiers)):
#                 if nodeIdentifiers[n] in allNodesOnLC:
#                     nearLCGroup.getMeshGroup(mesh).addElement(element)
#                     break
#             elementId = element.getIdentifier()
#             marginDorsal = False
#             for n in range(len(nodeIdentifiers)):
#                 marginElement = nodeIdentifiers[n] in nodesOnSplitMargin
#                 if marginElement:
#                     count += 1
#                     if count < 3 and (elementId <= elementsInOstium or elementId > closedLoopElementId):
#                         marginDorsal = True
#                     elif count >= 3 and (elementId <= elementsInOstium or elementId > closedLoopElementId):
#                         if count == 4:
#                             count = 0
#                     elif elementsInOstium < elementId < elementsInOstium + len(xOuter[0]) + 1:
#                         marginDorsal = True
#                         count = 0
#                     elif elementsInOstium + len(xOuter[0]) < elementId < elementsInOstium + len(xOuter[0]) * 2 + 1:
#                         count = 0
#                     elif count < 2 and elementId > elementsInOstium + 2 * (len(xOuter[0])):
#                         marginDorsal = True
#                     elif count >= 2 and elementId > elementsInOstium + 2 * (len(xOuter[0])):
#                         if count == 2:
#                             count = 0
#                     break
#
#             if marginDorsal:
#                 # Find nodes on margin to remap with version 2
#                 lnRemapV2 = []
#                 for n in range(len(nodeIdentifiers)):
#                     if nodeIdentifiers[n] in nodesOnSplitMargin:
#                         lnRemapV2.append(n + 1)
#                 eft2 = eft
#                 remapEftNodeValueLabelsVersion(eft2, lnRemapV2, allValueLabels, 2)
#
#                 splitElementtemplate2.defineField(splitCoordinates, -1, eft2)
#                 element.merge(splitElementtemplate2)
#                 element.setNodesByIdentifier(eft2, nodeIdentifiers)
#                 if eft2.getNumberOfLocalScaleFactors() > 0:
#                     element.setScaleFactor(eft2, 1, -1.0)
#             else:
#                 splitElementtemplate1.defineField(splitCoordinates, -1, eft)
#                 element.merge(splitElementtemplate1)
#                 element.setNodesByIdentifier(eft, nodeIdentifiers)
#                 if eft.getNumberOfLocalScaleFactors() == 1:
#                     element.setScaleFactors(eft, [-1.0])
#
#             element = elementIter.next()
#
#         allAnnotationGroups.append(nearLCGroup)

    return allAnnotationGroups, elementCountGroupList

def findRescaledCrossAxis(xSampled, d2Sampled, arcLengthPerSampledElement, xCentralPath, dCentralPath):
    """

    :param xSampled:
    :param d2Sampled:
    :param arcLengthPerSampledElement:
    :param xCentralPath:
    :param dCentralPath:
    :return:
    """
    for n in range(1, len(xCentralPath) - 1):
        cdProjected = [xCentralPath[n][c] + dCentralPath[n][c] for c in range(3)]
        index = interp.getNearestPointIndex(xSampled, cdProjected)
        arcDist = arcLengthPerSampledElement * index
        arcLow = arcDist - arcLengthPerSampledElement * 0.5
        arcUp = arcDist + arcLengthPerSampledElement * 0.5

        for iter in range(100):
            thetaLow = findThetaFromArcDistance(xSampled, d2Sampled, arcLow, xCentralPath[n], dCentralPath[n])
            thetaUp = findThetaFromArcDistance(xSampled, d2Sampled, arcUp, xCentralPath[n], dCentralPath[n])
            arcDistance = (arcLow + arcUp) * 0.5

            if thetaLow > thetaUp:
                arcLow = arcDistance
            elif thetaLow < thetaUp:
                arcUp = arcDistance
            else:
                xFundus = interp.getCubicHermiteCurvesPointAtArcDistance(xSampled, d2Sampled, arcDistance)[0]
                dCentralPath[n] = \
                    vector.setMagnitude(dCentralPath[n],
                                        vector.magnitude([xFundus[c] - xCentralPath[n][c] for c in range(3)]))
        if iter > 99:
            print('Search for fundus element boundaries on greater curvature - Max iters reached:', iter)

    return

def findThetaFromArcDistance(xSampled, d2Sampled, arcDistance, xCentralPath, dCentralPath):
    """

    :param xSampled:
    :param d2Sampled:
    :param arcDistance:
    :param xCentralPath:
    :param dCentralPath:
    :return:
    """
    x = interp.getCubicHermiteCurvesPointAtArcDistance(xSampled, d2Sampled, arcDistance)[0]
    dv = [x[c] - xCentralPath[c] for c in range(3)]
    cosTheta = vector.dotproduct(dv, dCentralPath) / (vector.magnitude(dv) * vector.magnitude(dCentralPath))
    cosTheta = 1 if cosTheta > 1.0 else cosTheta
    theta = math.acos(cosTheta)

    return theta