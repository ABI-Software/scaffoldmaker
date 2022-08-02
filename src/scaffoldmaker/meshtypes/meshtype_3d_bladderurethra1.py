"""
Generates 3D bladder and urethra meshes along the central path,
with variable numbers of elements around, along and through
wall.
"""

import copy
import math

from opencmiss.utils.zinc.field import findOrCreateFieldGroup, findOrCreateFieldNodeGroup, \
    findOrCreateFieldStoredMeshLocation, findOrCreateFieldStoredString
from opencmiss.zinc.element import Element
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, findOrCreateAnnotationGroupForTerm, \
    getAnnotationGroupForTerm
from scaffoldmaker.annotation.bladder_terms import get_bladder_term
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.geometry import createEllipsePoints
from scaffoldmaker.utils.interpolation import smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues, mesh_destroy_elements_and_nodes_by_identifiers


class MeshType_3d_bladderurethra1(Scaffold_base):
    """
    Generates 3D bladder and urethra meshes with variable numbers
    of elements around, along the central path, and through the wall.
    """
    centralPathDefaultScaffoldPackages_LUT = {
        'Cat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
                },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0000, 0.0000, 0.0000], [0.0000, 0.0000, 15.000], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, 0.0000, 15.000], [0.0000, 0.0000, 15.000], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, 0.0000, 30.000], [0.0000, 0.0000, 15.000], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, 0.0000, 45.000], [0.0000, 0.0000, 15.000], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, 0.1000, 60.450], [0.0000, 1.1000, 13.500], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, 1.4200, 75.200], [0.0000, 1.4000, 13.200], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, 3.0500, 89.940], [0.0000, 1.8000, 15.400], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, 5.9800, 107.94], [0.0000, 4.5000, 16.800], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, 12.470, 127.07], [0.0000, 6.6700, 14.800], [0.0000, 0.5000, 0.00000], [0.0000, 0.0000, -0.5000]]
                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('urinary bladder')[0],
                    'ontId': get_bladder_term('urinary bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_bladder_term('urethra')[0],
                    'ontId': get_bladder_term('urethra')[1]
                }]
            }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0000,  0.0000, 0.0000], [0.0000,  0.0000, 20.000], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000,  0.0000, 20.000], [0.0000,  0.0000, 20.000], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000,  0.0000, 40.000], [0.0000,  0.0000, 20.000], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000,  0.0000, 60.000], [0.0000,  0.0000, 20.000], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000,  0.0000, 80.000], [0.0000,  0.3000, 17.700], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, -2.1000, 100.60], [0.0000, -1.9000, 17.200], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, -3.8000, 120.30], [0.0000, -0.6000, 18.800], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000, -3.6000, 140.50], [0.0000,  2.0000, 18.400], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]],
                    [[0.0000,  0.0000, 160.00], [0.0000,  4.0000, 15.500], [0.0000, 0.5000, 0.0000], [0.0000, 0.0000, -0.5000]]
                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('urinary bladder')[0],
                    'ontId': get_bladder_term('urinary bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_bladder_term('urethra')[0],
                    'ontId': get_bladder_term('urethra')[1]
                }]
        }),
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.000,  0.0000, 0.0000], [0.0000,  0.2107, 19.185], [0.0000, 0.5000,  0.0006], [0.0000, -0.0002, -0.0098]],
                    [[0.000,  0.1095, 17.595], [0.0000,  0.0083, 16.005], [0.0000, 0.4999, -0.0074], [0.0000, -0.0000, -0.0062]],
                    [[0.000,  0.0347, 32.009], [0.0000,  0.0042, 13.831], [0.0000, 0.4999, -0.0122], [0.0000,  0.0001,  0.0027]],
                    [[0.000,  0.1115, 45.257], [0.0000, -0.0207, 13.330], [0.0000, 0.5000, -0.0026], [0.0000, -0.0015,  0.0336]],
                    [[0.000, -0.0079, 58.668], [0.0000, -0.0381, 13.242], [0.0000, 0.4996,  0.0196], [0.0000, -0.0040,  0.0453]],
                    [[0.000,  0.0332, 71.739], [0.0000, -0.2969, 13.684], [0.0000, 0.4983,  0.0403], [0.0000,  0.0013, -0.0227]],
                    [[0.000, -0.6320, 86.026], [0.0000, -0.5928, 14.006], [0.0000, 0.4997,  0.0193], [0.0000, -0.0027, -0.1007]],
                    [[0.000, -1.1544, 99.742], [0.0000,  0.7147, 16.377], [0.0000, 0.4961, -0.0620], [0.0000, -0.0226, -0.1114]],
                    [[0.000,  1.2247, 118.63], [0.0000,  4.0316, 21.331], [0.0000, 0.4719, -0.1653], [0.0000, -0.0498, -0.0958]]
                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('urinary bladder')[0],
                    'ontId': get_bladder_term('urinary bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_bladder_term('urethra')[0],
                    'ontId': get_bladder_term('urethra')[1]
                }]
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0000,  0.0000, 0.0000], [0.0000,  0.3967, 17.118], [0.0000, 0.4999, -0.0116], [0.0000,  0.0002,  0.0084]],
                    [[0.0000,  0.2526, 16.482], [0.0000,  0.1084, 15.844], [0.0000, 0.5000, -0.0034], [0.0000,  0.0001,  0.0073]],
                    [[0.0000,  0.2277, 31.687], [0.0000, -0.1005, 16.250], [0.0000, 0.5000,  0.0031], [0.0000,  0.0180,  0.0030]],
                    [[0.0000,  0.0421, 48.982], [0.0000, -0.1616, 17.639], [0.0000, 0.5385,  0.0049], [0.0000,  0.0162,  0.0406]],
                    [[0.0000, -0.0945, 66.965], [0.0000, -0.4986, 17.987], [0.0000, 0.5378,  0.0276], [0.0000, -0.0051,  0.0509]],
                    [[0.0000, -0.9552, 84.946], [0.0000, -0.6462, 18.734], [0.0000, 0.5374,  0.0334], [0.0000, -0.0173, -0.0718]],
                    [[0.0000, -1.3669, 104.42], [0.0000,  0.5462, 17.118], [0.0000, 0.4994, -0.0233], [0.0000, -0.0266, -0.1269]],
                    [[0.0000, -0.1150, 119.12], [0.0000,  2.7663, 14.068], [0.0000, 0.4902, -0.0980], [0.0000, -0.0116, -0.0490]],
                    [[0.0000,  4.0321, 132.37], [0.0000,  5.4802, 12.315], [0.0000, 0.4691, -0.1728], [0.0000,  0.0080,  0.0382]],
                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('urinary bladder')[0],
                    'ontId': get_bladder_term('urinary bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_bladder_term('urethra')[0],
                    'ontId': get_bladder_term('urethra')[1]
                }]
        }),
        'Rat 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'D2 derivatives': True,
                'Length': 1.0,
                'Number of elements': 8
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.000,  0.0000, 0.0000], [0.0000, -0.0219, 18.928], [0.0000, 0.5000,  0.0006], [0.0000,  0.0001, -0.0062]],
                    [[0.000,  0.1095, 17.595], [0.0000,  0.2409, 16.260], [0.0000, 0.4999, -0.0074], [0.0000, -0.0000, -0.0051]],
                    [[0.000,  0.4589, 32.519], [0.0000,  0.3415, 14.001], [0.0000, 0.4999, -0.0122], [0.0000, -0.0000,  0.0018]],
                    [[0.000,  0.7906, 45.597], [0.0000,  0.0678, 13.290], [0.0000, 0.5000, -0.0026], [0.0000, -0.0060,  0.0568]],
                    [[0.000,  0.5862, 59.093], [0.0000, -0.5418, 13.798], [0.0000, 0.4969,  0.0553], [0.0000, -0.0080,  0.0623]],
                    [[0.000, -0.3075, 73.183], [0.0000, -1.1224, 13.895], [0.0000, 0.4921,  0.0884], [0.0000,  0.0061, -0.0646]],
                    [[0.000, -1.6517, 86.875], [0.0000, -0.5654, 14.673], [0.0000, 0.5000,  0.0047], [0.0000, -0.0038, -0.1421]],
                    [[0.000, -1.3262, 102.46], [0.0000,  1.8796, 15.041], [0.0000, 0.4872, -0.1123], [0.0000, -0.0646, -0.1586]],
                    [[0.000,  1.9903, 116.76], [0.0000,  4.7166, 13.462], [0.0000, 0.4510, -0.2159], [0.0000, -0.1378, -0.1893]]
                ]),
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_bladder_term('urinary bladder')[0],
                    'ontId': get_bladder_term('urinary bladder')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '5-8',
                    'name': get_bladder_term('urethra')[0],
                    'ontId': get_bladder_term('urethra')[1]
                }]
        })
    }
    ostiumDefaultScaffoldPackages = {
        'Ureter Cat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements around ostium': 8,  # implemented for 8
                'Number of elements along': 1,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 2.2,
                'Ostium length': 0.5,
                'Ostium wall thickness': 0.5,
                'Ostium wall relative thicknesses': [1.0],
                'Use linear through ostium wall': True,
                'Vessel end length factor': 2.0,
                'Vessel inner diameter': 0.8,
                'Vessel wall thickness': 0.25,
                'Vessel wall relative thicknesses': [1.0],
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
        'Ureter Human 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,  # implemented for 8
                'Number of elements along': 1,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 3.0,
                'Ostium length': 0.5,
                'Ostium wall thickness': 0.5,
                'Ostium wall relative thicknesses': [1.0],
                'Use linear through ostium wall': True,
                'Vessel end length factor': 2.0,
                'Vessel inner diameter': 1.0,
                'Vessel wall thickness': 0.1,
                'Vessel wall relative thicknesses': [1.0],
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
        'Ureter Mouse 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,  # implemented for 8
                'Number of elements along': 1,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 1.3,
                'Ostium length': 0.25,
                'Ostium wall thickness': 0.3,
                'Ostium wall relative thicknesses': [1.0],
                'Use linear through ostium wall': True,
                'Vessel end length factor': 2.0,
                'Vessel inner diameter': 0.45,
                'Vessel wall thickness': 0.1,
                'Vessel wall relative thicknesses': [1.0],
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
        'Ureter Pig 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,  # implemented for 8
                'Number of elements along': 1,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 2.0,
                'Ostium length': 0.2,
                'Ostium wall thickness': 0.4,
                'Ostium wall relative thicknesses': [1.0],
                'Use linear through ostium wall': True,
                'Vessel end length factor': 2.0,
                'Vessel inner diameter': 0.9,
                'Vessel wall thickness': 0.2,
                'Vessel wall relative thicknesses': [1.0],
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
        'Ureter Rat 1': ScaffoldPackage(MeshType_3d_ostium1, {
            'scaffoldSettings': {
                'Number of vessels': 1,
                'Number of elements across common': 2,
                'Number of elements around ostium': 8,  # implemented for 8
                'Number of elements along': 1,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium diameter': 1.0,
                'Ostium length': 0.25,
                'Ostium wall thickness': 0.02,
                'Ostium wall relative thicknesses': [1.0],
                'Use linear through ostium wall': True,
                'Vessel end length factor': 2.0,
                'Vessel inner diameter': 0.3,
                'Vessel wall thickness': 0.1,
                'Vessel wall relative thicknesses': [1.0],
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
        })
    }

    @staticmethod
    def getName():
        return '3D Bladder with Urethra 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Cat 1',
            'Human 1',
            'Mouse 1',
            'Pig 1',
            'Rat 1'
            ]

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 1' in parameterSetName:
            centralPathOption_LUT = cls.centralPathDefaultScaffoldPackages_LUT['Human 1']
            ureterOption = cls.ostiumDefaultScaffoldPackages['Ureter Human 1']
        elif 'Mouse 1' in parameterSetName:
            centralPathOption_LUT = cls.centralPathDefaultScaffoldPackages_LUT['Mouse 1']
            ureterOption = cls.ostiumDefaultScaffoldPackages['Ureter Mouse 1']
        elif 'Pig 1' in parameterSetName:
            centralPathOption_LUT = cls.centralPathDefaultScaffoldPackages_LUT['Pig 1']
            ureterOption = cls.ostiumDefaultScaffoldPackages['Ureter Pig 1']
        elif 'Rat 1' in parameterSetName:
            centralPathOption_LUT = cls.centralPathDefaultScaffoldPackages_LUT['Rat 1']
            ureterOption = cls.ostiumDefaultScaffoldPackages['Ureter Rat 1']
        else:
            centralPathOption_LUT = cls.centralPathDefaultScaffoldPackages_LUT['Cat 1']
            ureterOption = cls.ostiumDefaultScaffoldPackages['Ureter Cat 1']
        options = {
            'Central path LUT': copy.deepcopy(centralPathOption_LUT),
            'Number of elements along bladder': 12,
            'Number of elements around': 12,
            'Number of elements through wall': 1,
            'Major diameter': 30.0,
            'Minor diameter': 25.0,
            'Neck diameter 1': 5.0,
            'Neck diameter 2': 4.0,
            'Wall thickness': 0.5,
            'Neck angle degrees': 45,
            'Include urethra': False,
            'Number of elements along urethra': 8,
            'Urethra diameter 1': 1.5,
            'Urethra diameter 2': 1.0,
            'Urethra wall thickness': 0.5,
            'Include ureter': False,
            'Ureter': copy.deepcopy(ureterOption),
            'Ureter position around': 0.67,  # should be on the dorsal part (> 0.5)
            'Ureter position down': 0.83,
            'Number of elements ureter radial': 2,
            'Use cross derivatives': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 4,
            'Refine number of elements along': 4,
            'Refine number of elements through wall': 1
            }
        if 'Human' in parameterSetName:
            options['Number of elements along bladder'] = 8
            options['Number of elements around'] = 12
            options['Major diameter'] = 75.0
            options['Minor diameter'] = 35.0
            options['Wall thickness'] = 0.2
            options['Neck diameter 1'] = 5.0
            options['Neck diameter 2'] = 4.0
            options['Wall thickness'] = 0.5
            options['Neck angle degrees'] = 115
            options['Ureter position around'] = 0.82  # should be on the dorsal part (> 0.5)
            options['Ureter position down'] = 0.63
            options['Number of elements along urethra'] = 8
            options['Urethra diameter 1'] = 1.5
            options['Urethra diameter 2'] = 1.0
            options['Urethra wall thickness'] = 0.5
        if 'Mouse' in parameterSetName:
            options['Number of elements along bladder'] = 8
            options['Major diameter'] = 20.0
            options['Minor diameter'] = 15.0
            options['Wall thickness'] = 0.5
            options['Neck diameter 1'] = 4.0
            options['Neck diameter 2'] = 3.0
            options['Neck angle degrees'] = 45
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
            options['Number of elements along urethra'] = 6
            options['Ureter position down'] = 0.865
            options['Urethra diameter 1'] = 0.75
            options['Urethra diameter 2'] = 0.65
            options['Urethra wall thickness'] = 0.25
        if 'Pig' in parameterSetName:
            options['Number of elements along bladder'] = 8
            options['Major diameter'] = 30.0
            options['Minor diameter'] = 25.0
            options['Wall thickness'] = 0.5
            options['Neck diameter 1'] = 6.0
            options['Neck diameter 2'] = 5.5
            options['Neck angle degrees'] = 45
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
            options['Ureter position down'] = 0.865
            options['Number of elements along urethra'] = 6
            options['Urethra diameter 1'] = 1.25
            options['Urethra diameter 2'] = 1.0
            options['Urethra wall thickness'] = 0.25
        if 'Rat' in parameterSetName:
            options['Number of elements along bladder'] = 12
            options['Major diameter'] = 20.0
            options['Minor diameter'] = 15.0
            options['Wall thickness'] = 0.2
            options['Neck diameter 1'] = 3.5
            options['Neck diameter 2'] = 3.0
            options['Neck angle degrees'] = 45
            options['Ureter position around'] = 0.67  # should be on the dorsal part (> 0.5)
            options['Ureter position down'] = 0.83
            options['Number of elements along urethra'] = 8
            options['Urethra diameter 1'] = 0.75
            options['Urethra diameter 2'] = 0.65
            options['Urethra wall thickness'] = 0.25
        return options

    @staticmethod
    def getOrderedOptionNames():
        optionNames = [
            'Central path LUT',
            'Number of elements along bladder',
            'Number of elements around',
            'Number of elements through wall',
            'Major diameter',
            'Minor diameter',
            'Neck diameter 1',
            'Neck diameter 2',
            'Wall thickness',
            'Neck angle degrees',
            'Include urethra',
            'Number of elements along urethra',
            'Urethra diameter 1',
            'Urethra diameter 2',
            'Urethra wall thickness',
            'Include ureter',
            'Ureter',
            'Ureter position around',
            'Ureter position down',
            'Number of elements ureter radial',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']
        return optionNames

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path LUT':
            return [MeshType_1d_path1]
        if optionName == 'Ureter':
            return [MeshType_3d_ostium1]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path LUT':
            return list(cls.centralPathDefaultScaffoldPackages_LUT.keys())
        if optionName == 'Ureter':
            return list(cls.ostiumDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + 'Invalid option \'' + \
            optionName + '\' scaffold type ' + scaffoldType.getName()
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
        if optionName == 'Central path LUT':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages_LUT.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages_LUT[parameterSetName])
        if optionName == 'Ureter':
            if not parameterSetName:
                parameterSetName = list(cls.ostiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.ostiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path LUT'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path LUT'):
            options['Central path LUT'] = cls.getOptionScaffoldPackage('Central path LUT', MeshType_1d_path1)
        for key in [
            'Number of elements along bladder',
            'Number of elements around',
            'Number of elements through wall',
            'Number of elements along urethra',
            'Number of elements ureter radial',
            'Refine number of elements along',
            'Refine number of elements around',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements around'] % 2 != 0:
            options['Number of elements around'] += 1
        if options['Ureter position around'] < 0.5:
            options['Ureter position around'] = 0.5  # ureters are on the dorsal part of the bladder
        elif options['Ureter position around'] > 0.9:
            options['Ureter position around'] = 0.9
        if options['Ureter position down'] < 0.15:
            options['Ureter position down'] = 0.15
        elif options['Ureter position down'] > 0.95:
            options['Ureter position down'] = 0.95

    @classmethod
    def updateSubScaffoldOptions(cls, options):
        """
        Update ostium sub-scaffold options which depend on parent options.
        """
        bladderWallThickness = options['Wall thickness']
        elementsCountThroughWall = options['Number of elements through wall']
        ureterOptions = options['Ureter']
        ureterDefaultOptions = ureterOptions.getScaffoldSettings()
        ureterDefaultOptions['Ostium wall thickness'] = bladderWallThickness
        ureterDefaultOptions['Number of elements through wall'] = elementsCountThroughWall

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """
        cls.updateSubScaffoldOptions(options)
        centralPath = options['Central path LUT']
        elementsCountAlongBladder = options['Number of elements along bladder']
        elementsCountAround = options['Number of elements around']
        elementsCountThroughWall = options['Number of elements through wall']
        majorDiameter = options['Major diameter']
        minorDiameter = options['Minor diameter']
        neckDiameter1 = options['Neck diameter 1']
        neckDiameter2 = options['Neck diameter 2']
        bladderWallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        neckAngleRadians = math.radians(options['Neck angle degrees'])

        includeUrethra = options['Include urethra']
        elementsCountAlongUrethra = options['Number of elements along urethra']
        urethraDiameter1 = options['Urethra diameter 1']
        urethraDiameter2 = options['Urethra diameter 2']
        urethraWallThickness = options['Urethra wall thickness']

        includeUreter = options['Include ureter']
        ureterOptions = options['Ureter']
        ureterDefaultOptions = ureterOptions.getScaffoldSettings()
        elementsCountAroundUreter = ureterDefaultOptions['Number of elements around ostium']
        elementsCountUreterRadial = options['Number of elements ureter radial']
        ureterPositionAround = options['Ureter position around']
        ureterPositionDown = options['Ureter position down']

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        fm = region.getFieldmodule()
        fm.beginChange()
        mesh = fm.findMeshByDimension(3)

        # Annotation fiducial point
        markerGroup = findOrCreateFieldGroup(fm, "marker")
        markerName = findOrCreateFieldStoredString(fm, name="marker_name")
        markerLocation = findOrCreateFieldStoredMeshLocation(fm, mesh, name="marker_location")

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        markerPoints = findOrCreateFieldNodeGroup(markerGroup, nodes).getNodesetGroup()
        markerTemplateInternal = nodes.createNodetemplate()
        markerTemplateInternal.defineField(markerName)
        markerTemplateInternal.defineField(markerLocation)

        cache = fm.createFieldcache()

        # Central path
        # Bladder part
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx_bladder, cd1_bladder, cd2_bladder, cd12_bladder = \
            extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                        Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2],
                                            groupName='urinary bladder')
        # for i in range(len(cx_bladder)):
        #     print(i, '[', cx_bladder[i], ',', cd1_bladder[i], ',', cd2_bladder[i], ',', cd12_bladder[i], '],')
        del tmpRegion
        # Urethra part
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx_urethra, cd1_urethra, cd2_urethra, cd12_urethra = \
            extractPathParametersFromRegion(tmpRegion, [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                        Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2],
                                            groupName='urethra')
        # for i in range(len(cx_urethra)):
        #     print(i, '[', cx_urethra[i], ',', cd1_urethra[i], ',', cd2_urethra[i], ',', cd12_urethra[i], '],')
        del tmpRegion

        # Find arcLength
        # Bladder part
        bladderLength = 0.0
        elementsCountInBladder = len(cx_bladder) - 1
        for e in range(elementsCountInBladder):
            arcLength = interp.getCubicHermiteArcLength(cx_bladder[e], cd1_bladder[e],
                                                        cx_bladder[e + 1], cd1_bladder[e + 1])
            bladderLength += arcLength
        bladderSegmentLength = bladderLength / elementsCountAlongBladder
        # Urethra part
        urethraLength = 0.0
        elementsCountInUrethra = len(cx_urethra) - 1
        for e in range(elementsCountInUrethra):
            arcLength = interp.getCubicHermiteArcLength(cx_urethra[e], cd1_urethra[e],
                                                        cx_urethra[e + 1], cd1_urethra[e + 1])
            urethraLength += arcLength
        urethraSegmentLength = urethraLength / elementsCountAlongUrethra

        if includeUrethra:
            length = bladderLength + urethraLength
            elementsCountAlong = elementsCountAlongBladder + elementsCountAlongUrethra
            cx = cx_bladder + cx_urethra[1:]
            cd1 = cd1_bladder + cd1_urethra[1:]
            cd2 = cd2_bladder + cd2_urethra[1:]
            cd12 = cd12_bladder + cd12_urethra[1:]
        else:
            length = bladderLength
            elementsCountAlong = elementsCountAlongBladder
            cx = cx_bladder
            cd1 = cd1_bladder
            cd2 = cd2_bladder
            cd12 = cd12_bladder

        # Sample central path
        sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
        sd2, sd12 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

        # Create bladder
        # Create line down the bladder with major radius
        R1_max = 0.5 * majorDiameter
        R2_max = 0.5 * neckDiameter1
        # Create the body part
        b = R1_max
        a = bladderLength / 2
        xLoop, d2Loop = createEllipsePoints([0.0, 0.0, a], math.pi - neckAngleRadians, [0.0, 0.0, a], [0.0, b, 0.0],
                                            elementsCountAlongBladder, startRadians=math.pi)
        d1Loop = [0, 0, 0] * len(xLoop)
        nodesAlongMax_x = []
        nodesAlongMax_d1 = []
        nodesAlongMax_d2 = []
        for n in range(len(xLoop)):
            nodesAlongMax_x.append(xLoop[n])
            nodesAlongMax_d1.append(d1Loop[n])
            nodesAlongMax_d2.append(d2Loop[n])
        # Create the neck part
        px = [[nodesAlongMax_x[-1][0], nodesAlongMax_x[-1][1], nodesAlongMax_x[-1][2]], [0, (0 - R2_max), 2 * a]]
        pd2 = [nodesAlongMax_d2[-1], [0, 0, bladderSegmentLength]]
        pd2_smoothed = interp.smoothCubicHermiteDerivativesLine(px, pd2, fixAllDirections=False,
                                                                fixStartDerivative=True, fixEndDerivative=True,
                                                                fixStartDirection=False, fixEndDirection=False)
        pd1 = [0, 0, 0] * len(px)
        for n in range(1, len(pd2_smoothed)):
            nodesAlongMax_x.append(px[n])
            nodesAlongMax_d2.append(pd2_smoothed[n])
            nodesAlongMax_d1.append(pd1[n])
        nx_max, nd2_max, ne_max, nxi_max, nnf_max = interp.sampleCubicHermiteCurves(nodesAlongMax_x, nodesAlongMax_d2,
                                                                                    elementsCountAlongBladder)

        # Create line down the bladder with minor radius
        R1_min = 0.5 * minorDiameter
        R2_min = 0.5 * neckDiameter2
        # Create the body part
        a = bladderLength / 2
        b = R1_min
        xLoop, d2Loop = createEllipsePoints([0.0, 0.0, a], math.pi - neckAngleRadians, [0.0, 0.0, a], [0.0, b, 0.0],
                                            elementsCountAlongBladder, startRadians=math.pi)
        d1Loop = [0, 0, 0] * len(xLoop)
        nodesAlongMin_x = []
        nodesAlongMin_d1 = []
        nodesAlongMin_d2 = []
        for n in range(len(xLoop)):
            nodesAlongMin_x.append(xLoop[n])
            nodesAlongMin_d1.append(d1Loop[n])
            nodesAlongMin_d2.append(d2Loop[n])
        # Create the neck part
        px = [[nodesAlongMin_x[-1][0], nodesAlongMin_x[-1][1], nodesAlongMin_x[-1][2]], [0, (0 - R2_min), 2 * a]]
        pd2 = [nodesAlongMin_d2[-1], [0, 0, bladderSegmentLength]]
        pd2_smoothed = interp.smoothCubicHermiteDerivativesLine(px, pd2, fixAllDirections=False,
                                                                fixStartDerivative=True, fixEndDerivative=True,
                                                                fixStartDirection=False, fixEndDirection=False)
        pd1 = [0, 0, 0] * len(px)
        for n in range(1, len(pd2_smoothed)):
            nodesAlongMin_x.append(px[n])
            nodesAlongMin_d2.append(pd2_smoothed[n])
            nodesAlongMin_d1.append(pd1[n])
        nx_min, nd2_min, ne_min, nxi_min, nnf_min = interp.sampleCubicHermiteCurves(nodesAlongMin_x, nodesAlongMin_d2,
                                                                                    elementsCountAlongBladder)

        # Create the whole bladder nodes
        # Create ellipses surrounding the center line
        innerNodes_x = []
        innerNodes_d1 = []
        innerNodes_d2 = []
        nd1 = matrix.rotateAboutZAxis(nd2_max[0], 0.5 * math.pi)
        innerNodes_d1 += [nd1] * elementsCountAround
        for n1 in range(0, len(nx_max)):
            xAround, d1Around = createEllipsePoints([0.0, 0.0, nx_max[n1][2]], 2 * math.pi, [nx_max[n1][1], 0.0, 0.0],
                                                    [0.0, nx_min[n1][1], 0.0], elementsCountAround,
                                                    startRadians=-math.pi/2)
            innerNodes_x += xAround
            if n1 >= 1:
                innerNodes_d1 += d1Around
        # Create lines from apex go down the bladder
        nodesInLineToDown_x = []
        nodesInLineToDown_d2 = []
        for n2 in range(elementsCountAround):
            x = [innerNodes_x[n1 * elementsCountAround + n2] for n1 in range(0, elementsCountAlongBladder + 1)]
            nodesInLineToDown_x += x
            for n1 in range(0, elementsCountAlongBladder + 1):
                if n1 == 0:
                    d2 = [bladderSegmentLength, 0.0, 0.0]
                elif n1 == elementsCountAlongBladder:
                    if not includeUrethra:
                        d2 = [0.0, 0.0, bladderSegmentLength]
                else:
                    v1 = nodesInLineToDown_x[(elementsCountAlongBladder + 1) * n2 + n1 - 1]
                    v2 = nodesInLineToDown_x[(elementsCountAlongBladder + 1) * n2 + n1]
                    d2 = [v2[c] - v1[c] for c in range(3)]
                nodesInLineToDown_d2.append(d2)
        # Smoothing the derivatives
        smoothed_d2 = []
        for n1 in range(elementsCountAround):
            lineSmoothingNodes = []
            lineSmoothingNodes_d2 = []
            for n2 in range(elementsCountAlongBladder + 1):
                lineSmoothingNodes.append(nodesInLineToDown_x[n1 * (elementsCountAlongBladder + 1) + n2])
                lineSmoothingNodes_d2.append(nodesInLineToDown_d2[n1 * (elementsCountAlongBladder + 1) + n2])
            smd2 = smoothCubicHermiteDerivativesLine(lineSmoothingNodes, lineSmoothingNodes_d2, fixAllDirections=True,
                                                    fixStartDerivative=True, fixEndDerivative=True,
                                                    fixStartDirection=False, fixEndDirection=False)
            smoothed_d2 += smd2
        # Re-arrange the derivatives order
        for n2 in range(elementsCountAlongBladder + 1):
            for n1 in range(elementsCountAround):
                rd2 = smoothed_d2[n1 * (elementsCountAlongBladder + 1) + n2]
                innerNodes_d2.append(rd2)

        # Store derivatives for transition elements to replace in d2List
        d2transitionList = []
        for n in range(0, elementsCountAround):
            transit_d2 = innerNodes_d2[elementsCountAround * elementsCountAlongBladder + n]
            d2transitionList.append(transit_d2)

        # Create urethra
        if includeUrethra:
            last_z_bladder = innerNodes_x[-1][2]
            radiansPerElementAround = 2.0 * math.pi / elementsCountAround
            transitLength = (bladderSegmentLength + urethraSegmentLength) / 2
            newSegmentLength = (urethraLength - transitLength) / (elementsCountAlongUrethra - 1)
            for n2 in range(0, elementsCountAlongUrethra + 1):
                for n1 in range(elementsCountAround):
                    radiansAround = n1 * radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    if n2 == 1:
                        x = [
                            -urethraDiameter1 * sinRadiansAround,
                            urethraDiameter2 * cosRadiansAround,
                            last_z_bladder + transitLength
                        ]
                        dx_ds1 = [
                            -radiansPerElementAround * urethraDiameter1 * cosRadiansAround,
                            radiansPerElementAround * urethraDiameter2 * -sinRadiansAround,
                            0.0
                        ]
                        dx_ds2 = [0, 0, newSegmentLength]
                    else:
                        x = [
                            -urethraDiameter1 * sinRadiansAround,
                            urethraDiameter2 * cosRadiansAround,
                            last_z_bladder + transitLength + (n2 - 1) * newSegmentLength
                        ]
                        dx_ds1 = [
                            -radiansPerElementAround * urethraDiameter1 * cosRadiansAround,
                            radiansPerElementAround * urethraDiameter2 * -sinRadiansAround,
                            0.0
                        ]
                        dx_ds2 = [0, 0, newSegmentLength]
                    if n2 == 0:
                        pass
                    else:
                        innerNodes_x.append(x)
                        innerNodes_d1.append(dx_ds1)
                        innerNodes_d2.append(dx_ds2)

        # Project reference point for warping onto central path
        sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
             tubemesh.getPlaneProjectionOnCentralPath(innerNodes_x, elementsCountAround, elementsCountAlong,
                                                      length, sx, sd1, sd2, sd12)

        innerRadiusAlong = []
        for n2 in range(elementsCountAlong + 1):
            firstNodeAlong = innerNodes_x[n2 * elementsCountAround]
            radius = vector.magnitude(firstNodeAlong[c] - [0.0, 0.0, firstNodeAlong[2]][c] for c in range(3))
            innerRadiusAlong.append(radius)

        # Warp points
        segmentAxis = [0.0, 0.0, 1.0]
        xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = \
            tubemesh.warpSegmentPoints(innerNodes_x, innerNodes_d1, innerNodes_d2, segmentAxis, sxRefList, sd1RefList,
                                       sd2ProjectedListRef, elementsCountAround, elementsCountAlong,
                                       zRefList, innerRadiusAlong, closedProximalEnd=True)

        if includeUrethra:
            wallThicknessList = [bladderWallThickness] * (elementsCountAlongBladder + 1) + [urethraWallThickness] * \
                                elementsCountAlongUrethra
        else:
            wallThicknessList = [bladderWallThickness] * (elementsCountAlong + 1)
        transitElementList = [0] * elementsCountAround

        relativeThicknessList = []
        xList, d1List, d2List, d3List, curvatureList = tubemesh.getCoordinatesFromInner(xWarpedList, d1WarpedList,
                                                                                        d2WarpedList, d3WarpedUnitList,
                                                                                        wallThicknessList,
                                                                                        relativeThicknessList,
                                                                                        elementsCountAround,
                                                                                        elementsCountAlong,
                                                                                        elementsCountThroughWall,
                                                                                        transitElementList)
        # Call the derivatives from the transition list to be replaced in the d2List
        idx = elementsCountAlongBladder * elementsCountAround
        for n2 in range(elementsCountThroughWall + 1):
            for n1 in range(0, elementsCountAround):
                new_d2 = d2transitionList[n1]
                d2List[(elementsCountThroughWall + 1) * idx + n2 * elementsCountAround + n1] = new_d2

        # Deal with multiple nodes at end point for closed proximal end
        xApexInner = xList[0]
        # Arclength between apex point and corresponding point on next face
        mag = interp.getCubicHermiteArcLength(xList[0], d2List[0], xList[2 * elementsCountAround],
                                              d2List[2 * elementsCountAround])
        d2ApexInner = vector.setMagnitude(sd2[0], mag)
        d1ApexInner = vector.crossproduct3(sd1[0], d2ApexInner)
        d1ApexInner = vector.setMagnitude(d1ApexInner, mag)
        d3ApexUnit = vector.normalise(
            vector.crossproduct3(vector.normalise(d1ApexInner), vector.normalise(d2ApexInner)))
        d3ApexInner = [d3ApexUnit[c] * bladderWallThickness / elementsCountThroughWall for c in range(3)]

        xFinal = []
        d1Final = []
        d2Final = []
        d3Final = []
        for n3 in range(elementsCountThroughWall + 1):
            xApex = [xApexInner[c] +
                     d3ApexUnit[c] * bladderWallThickness / elementsCountThroughWall * n3 for c in range(3)]
            xFinal.append(xApex)
            d1Final.append(d1ApexInner)
            d2Final.append(d2ApexInner)
            d3Final.append(d3ApexInner)

        xFinal += xList[(elementsCountThroughWall + 1) * elementsCountAround:]
        d1Final += d1List[(elementsCountThroughWall + 1) * elementsCountAround:]
        d2Final += d2List[(elementsCountThroughWall + 1) * elementsCountAround:]
        d3Final += d3List[(elementsCountThroughWall + 1) * elementsCountAround:]

        # Obtain elements count along body and neck of the bladder for defining annotation groups
        elementsCountAlongBody = round(ureterPositionDown * elementsCountAlongBladder - 1)
        elementsCountAlongNeck = elementsCountAlongBladder - elementsCountAlongBody

        # Define outer layer of the bladder to create the trackSurface on it
        outerNodes_x = []
        outerNodes_d1 = []
        outerNodes_d2 = []
        for n2 in range(elementsCountAlongBladder + 1):
            for n1 in range(elementsCountAround):
                idx = n2 * (elementsCountThroughWall + 1) * elementsCountAround + \
                      elementsCountThroughWall * elementsCountAround + n1
                outerNodes_x.append(xList[idx])
                outerNodes_d1.append(d1List[idx])
                outerNodes_d2.append(d2List[idx])

        elementsCount1 = elementsCountAround // 2
        elementsCount2 = elementsCountAlongBladder
        # Create trackSurface at the outer layer of the bladder for ureter 1
        nodesOnTrackSurface_x = []
        nodesOnTrackSurface_d1 = []
        nodesOnTrackSurface_d2 = []
        for n2 in range(elementsCountAlongBladder + 1):
            for n1 in range(elementsCountAround // 2 + 1):
                nodesOnTrackSurface_x.append(outerNodes_x[n2 * elementsCountAround + n1])
                nodesOnTrackSurface_d1.append(outerNodes_d1[n2 * elementsCountAround + n1])
                nodesOnTrackSurface_d2.append(outerNodes_d2[n2 * elementsCountAround + n1])
        trackSurfaceUreter1 = TrackSurface(elementsCount1, elementsCount2, nodesOnTrackSurface_x,
                                           nodesOnTrackSurface_d1, nodesOnTrackSurface_d2)

        ureter1Position = trackSurfaceUreter1.createPositionProportion(ureterPositionAround, ureterPositionDown)
        ureterElementPositionAround = ureter1Position.e1
        ureterElementPositionDown = ureter1Position.e2
        elementToDeleteStartIdx1 = elementsCountThroughWall * elementsCountAround * \
                                   (ureterElementPositionDown - (0 if ureter1Position.xi2 > 0.5 else 1)) + \
                                   ureterElementPositionAround + (1 if ureter1Position.xi1 > 0.5 else 0)

        elementIndex = elementToDeleteStartIdx1 - elementsCountAround * ureterElementPositionDown

        # Create trackSurface at the outer layer of the bladder for ureter 2
        nodesOnTrackSurface2_x = []
        nodesOnTrackSurface2_d1 = []
        nodesOnTrackSurface2_d2 = []
        for n2 in range(elementsCountAlongBladder + 1):
            for n1 in range(elementsCountAround // 2, elementsCountAround):
                nodesOnTrackSurface2_x.append(outerNodes_x[n2 * elementsCountAround + n1])
                nodesOnTrackSurface2_d1.append(outerNodes_d1[n2 * elementsCountAround + n1])
                nodesOnTrackSurface2_d2.append(outerNodes_d2[n2 * elementsCountAround + n1])
            nodesOnTrackSurface2_x.append(outerNodes_x[n2 * elementsCountAround])
            nodesOnTrackSurface2_d1.append(outerNodes_d1[n2 * elementsCountAround])
            nodesOnTrackSurface2_d2.append(outerNodes_d2[n2 * elementsCountAround])

        trackSurfaceUreter2 = TrackSurface(elementsCount1, elementsCount2, nodesOnTrackSurface2_x,
                                           nodesOnTrackSurface2_d1, nodesOnTrackSurface2_d2)
        ureter2Position = TrackSurfacePosition(elementsCountAround // 2 -
                                               ureterElementPositionAround + (-1 if ureter1Position.xi1 > 0 else 0),
                                               ureterElementPositionDown,
                                               (1 - ureter1Position.xi1) if ureter1Position.xi1 > 0
                                               else ureter1Position.xi1,
                                               ureter1Position.xi2)

        if includeUreter:
            if ureterElementPositionAround > (elementsCountAround / 4):
                elementsCountVentral1 = (elementsCountAround // 4)
            else:
                elementsCountVentral1 = elementIndex - 1
        else:
            elementsCountVentral1 = elementsCountAround // 4

        # Create annotation groups for bladder and urethra
        bodyGroup = AnnotationGroup(region, get_bladder_term("dome of the bladder"))
        neckGroup = AnnotationGroup(region, get_bladder_term("neck of urinary bladder"))
        bladderGroup = AnnotationGroup(region, get_bladder_term("urinary bladder"))
        urethraGroup = AnnotationGroup(region, get_bladder_term("urethra"))
        ureterGroup = AnnotationGroup(region, get_bladder_term("ureter"))
        if includeUrethra:
            elementsCountAlongGroups = [elementsCountAlongBody, elementsCountAlongNeck, elementsCountAlongUrethra]
            annotationGroupAlong = [[bladderGroup, bodyGroup], [bladderGroup, neckGroup], [urethraGroup]]
        else:
            elementsCountAlongGroups = [elementsCountAlongBody, elementsCountAlongNeck]
            annotationGroupAlong = [[bladderGroup, bodyGroup], [bladderGroup, neckGroup]]

        annotationGroupsAlong = []
        for i in range(len(elementsCountAlongGroups)):
            elementsCount = elementsCountAlongGroups[i]
            for n in range(elementsCount):
                annotationGroupsAlong.append(annotationGroupAlong[i])

        annotationGroupsAround = []
        for i in range(elementsCountAround):
            annotationGroupsAround.append([])

        annotationGroupsThroughWall = []
        for i in range(elementsCountThroughWall):
            annotationGroupsThroughWall.append([])

        neckMeshGroup = neckGroup. getMeshGroup(mesh)
        urinaryBladderMeshGroup = bladderGroup.getMeshGroup(mesh)
        ureterMeshGroup = ureterGroup. getMeshGroup(mesh)

        xOrgan = d1Organ = d2Organ = []

        # Create ureters and urethra only for coordinates (not for flat coordinates)
        if includeUrethra or includeUreter:
            xFlat = d1Flat = d2Flat = []
        else:
            # Obtain flat nodes coordinates
            xFlat, d1Flat, d2Flat = obtainBladderFlatNodes(elementsCountAlongBladder, elementsCountAround, elementsCountThroughWall,
                                                           xFinal, d1Final, d2Final, bladderLength, neckDiameter1, xApexInner, d2ApexInner,
                                                           bladderWallThickness)

        # Create nodes and elements
        nodeIdentifier, elementIdentifier, annotationGroups = tubemesh.createNodesAndElements(
            region, xFinal, d1Final, d2Final, d3Final, xFlat, d1Flat, d2Flat, xOrgan, d1Organ, d2Organ, None,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall,
            annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier,
            useCubicHermiteThroughWall, useCrossDerivatives, closedProximalEnd=True)

        if includeUreter:
            annotationGroups.append(ureterGroup)

        # Create annotation groups for dorsal and ventral parts of the bladder and urethra
        bladderDorsalGroup = AnnotationGroup(region, get_bladder_term("dorsal part of bladder"))
        bladderVentralGroup = AnnotationGroup(region, get_bladder_term("ventral part of bladder"))
        leftHemibladderGroup = AnnotationGroup(region, get_bladder_term("left hemi-bladder"))
        rightHemibladderGroup = AnnotationGroup(region, get_bladder_term("right hemi-bladder"))
        dorsalBladderMeshGroup = bladderDorsalGroup.getMeshGroup(mesh)
        ventralBladderMeshGroup = bladderVentralGroup.getMeshGroup(mesh)
        leftHemibladderMeshGroup = leftHemibladderGroup.getMeshGroup(mesh)
        rightHemibladderMeshGroup = rightHemibladderGroup.getMeshGroup(mesh)
        if includeUrethra:
            urethraDorsalGroup = AnnotationGroup(region, get_bladder_term("dorsal part of urethra"))
            urethraVentralGroup = AnnotationGroup(region, get_bladder_term("ventral part of urethra"))
            dorsalUrethraMeshGroup = urethraDorsalGroup.getMeshGroup(mesh)
            ventralUrethraMeshGroup = urethraVentralGroup.getMeshGroup(mesh)

        for e2 in range(elementsCountAlong):
            for e3 in range(elementsCountThroughWall):
                for e1 in range(elementsCountAround):
                    elementIdx = e3 * elementsCountAround + e2 * elementsCountThroughWall * elementsCountAround + e1 + 1
                    if e1 < elementsCountVentral1 or e1 > elementsCountAround - elementsCountVentral1 - 1:
                        element = mesh.findElementByIdentifier(elementIdx)
                        if e2 < elementsCountAlongBladder:
                            ventralBladderMeshGroup.addElement(element)
                            if e1 < elementsCountAround // 2:
                                leftHemibladderMeshGroup.addElement(element)
                            else:
                                rightHemibladderMeshGroup.addElement(element)
                        else:
                            ventralUrethraMeshGroup.addElement(element)
                    else:
                        element = mesh.findElementByIdentifier(elementIdx)
                        if e2 < elementsCountAlongBladder:
                            dorsalBladderMeshGroup.addElement(element)
                            if e1 < elementsCountAround // 2:
                                leftHemibladderMeshGroup.addElement(element)
                            else:
                                rightHemibladderMeshGroup.addElement(element)
                        else:
                            dorsalUrethraMeshGroup.addElement(element)

        annotationGroups.append(bladderDorsalGroup)
        annotationGroups.append(bladderVentralGroup)
        annotationGroups.append(leftHemibladderGroup)
        annotationGroups.append(rightHemibladderGroup)
        if includeUrethra:
            annotationGroups.append(urethraDorsalGroup)
            annotationGroups.append(urethraVentralGroup)

        if includeUreter:
            bladderMeshGroup = [neckMeshGroup, urinaryBladderMeshGroup, dorsalBladderMeshGroup]
            nodeIdentifier = generateUreterInlets(region, nodes, mesh, ureterDefaultOptions, elementsCountAround,
                                                  elementsCountThroughWall,
                                                  elementsCountAroundUreter, trackSurfaceUreter1, ureter1Position,
                                                  trackSurfaceUreter2,
                                                  ureter2Position, ureterElementPositionDown,
                                                  ureterElementPositionAround, xFinal, d1Final,
                                                  d2Final, nodeIdentifier, elementIdentifier,
                                                  elementsCountUreterRadial,
                                                  ureterMeshGroup, bladderMeshGroup)

        # Define markers for apex, ureter and urethra junctions with bladder
        apexGroup = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region, get_bladder_term("apex of urinary bladder"))
        leftUreterGroup = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("left ureter junction with bladder"))
        rightUreterGroup = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("right ureter junction with bladder"))
        dorsalUrethraGroup = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("urethra junction of dorsal bladder neck"))
        ventralUrethraGroup = \
            findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                               get_bladder_term("urethra junction of ventral bladder neck"))

        idx1 = 1
        xi1 = [0.0, 0.0, 0.0]
        markerList = []
        markerList.append({"group": apexGroup, "elementId": idx1, "xi": xi1})
        if includeUreter:
            idx2 = elementsCountAlong * elementsCountAround * elementsCountThroughWall + elementsCountAroundUreter
            xi2 = [0.0, 1.0, 0.0]
            markerList.append({"group": leftUreterGroup, "elementId": idx2, "xi": xi2})
            idx3 = elementsCountAlong * elementsCountAround * elementsCountThroughWall + \
                   elementsCountThroughWall * elementsCountAroundUreter + elementsCountAroundUreter
            xi3 = [0.0, 1.0, 0.0]
            markerList.append({"group": rightUreterGroup, "elementId": idx3, "xi": xi3})
        else:
            idx2 = ureterElementPositionDown * elementsCountAround * elementsCountThroughWall + \
                   ureterElementPositionAround + 1
            xi2 = [ureter1Position.xi1, ureter1Position.xi2, 0.0]
            markerList.append({"group": leftUreterGroup, "elementId": idx2, "xi": xi2})
            idx3 = ureterElementPositionDown * elementsCountAround * elementsCountThroughWall + \
                   elementsCountAround - ureterElementPositionAround
            xi3 = [1 - ureter1Position.xi1, ureter1Position.xi2, 0.0]
            markerList.append({"group": rightUreterGroup, "elementId": idx3, "xi": xi3})
        idx4 = (elementsCountAlongBladder - 1) * elementsCountAround * elementsCountThroughWall + \
               elementsCountAround // 2
        xi4 = [1.0, 1.0, 0.0]
        markerList.append({"group": dorsalUrethraGroup, "elementId": idx4, "xi": xi4})
        idx5 = (elementsCountAlongBladder - 1) * elementsCountAround * elementsCountThroughWall + 1
        xi5 = [0.0, 1.0, 0.0]
        markerList.append({"group": ventralUrethraGroup, "elementId": idx5, "xi": xi5})

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
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountAround = options['Refine number of elements around']
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

        if options['Include urethra'] == True:
            urethraGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("urethra"))
            urethraVentralGroup = getAnnotationGroupForTerm(annotationGroups,
                                                            get_bladder_term("ventral part of urethra"))
            urethraDorsalGroup = getAnnotationGroupForTerm(annotationGroups, get_bladder_term("dorsal part of urethra"))

            is_urethra = urethraGroup.getFieldElementGroup(mesh2d)
            is_urethra_serosa = fm.createFieldAnd(is_urethra, is_exterior_face_xi3_1)
            is_urethra_lumen = fm.createFieldAnd(is_urethra, is_exterior_face_xi3_0)

            is_dorsal_urethra = urethraDorsalGroup.getFieldElementGroup(mesh2d)
            is_ventral_urethra = urethraVentralGroup.getFieldElementGroup(mesh2d)

            serosaOfUrethra = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                 get_bladder_term("serosa of urethra"))
            serosaOfUrethra.getMeshGroup(mesh2d).addElementsConditional(is_urethra_serosa)
            lumenOfUrethra = findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                                get_bladder_term("lumen of urethra"))
            lumenOfUrethra.getMeshGroup(mesh2d).addElementsConditional(is_urethra_lumen)

            is_urethra_serosa_dorsal = fm.createFieldAnd(is_urethra_serosa, is_dorsal_urethra)
            serosaOfUrethra_dorsal = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_bladder_term("dorsal part of serosa of urethra"))
            serosaOfUrethra_dorsal.getMeshGroup(mesh2d).addElementsConditional(is_urethra_serosa_dorsal)

            is_urethra_serosa_ventral = fm.createFieldAnd(is_urethra_serosa, is_ventral_urethra)
            serosaOfUrethra_ventral = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_bladder_term("ventral part of serosa of urethra"))
            serosaOfUrethra_ventral.getMeshGroup(mesh2d).addElementsConditional(is_urethra_serosa_ventral)

            is_urethra_lumen_dorsal = fm.createFieldAnd(is_urethra_lumen, is_dorsal_urethra)
            lumenOfUrethra_dorsal = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_bladder_term("dorsal part of lumen of urethra"))
            lumenOfUrethra_dorsal.getMeshGroup(mesh2d).addElementsConditional(is_urethra_lumen_dorsal)

            is_urethra_lumen_ventral = fm.createFieldAnd(is_urethra_lumen, is_ventral_urethra)
            lumenOfUrethra_ventral = \
                findOrCreateAnnotationGroupForTerm(annotationGroups, region,
                                                   get_bladder_term("ventral part of lumen of urethra"))
            lumenOfUrethra_ventral.getMeshGroup(mesh2d).addElementsConditional(is_urethra_lumen_ventral)


def generateUreterInlets(region, nodes, mesh, ureterDefaultOptions, elementsCountAround, elementsCountThroughWall,
                         elementsCountAroundUreter, trackSurfaceUreter1, ureter1Position, trackSurfaceUreter2,
                         ureter2Position, ureterElementPositionDown, ureterElementPositionAround,
                         xBladder, d1Bladder, d2Bladder, nextNodeIdentifier, nextElementIdentifier,
                         elementsCountUreterRadial, ureterMeshGroup, bladderMeshGroup):

    # Update ostium and vessel wall relative thickness with elementsCountThroughWall
    # Set to uniform layer thicknesses now, can be changed to varying thickness later
    ureterDefaultOptions['Ostium wall relative thicknesses'] = [1.0 / elementsCountThroughWall] * \
                                                               elementsCountThroughWall
    ureterDefaultOptions['Vessel wall relative thicknesses'] = [1.0 / elementsCountThroughWall] * \
                                                               elementsCountThroughWall

    # Create ureters on the surface
    # Ureter 1
    centerUreter1_x, centerUreter1_d1, centerUreter1_d2 = \
        trackSurfaceUreter1.evaluateCoordinates(ureter1Position, derivatives=True)
    td1, td2, td3 = calculate_surface_axes(centerUreter1_d1, centerUreter1_d2, [1.0, 0.0, 0.0])
    endPointStartId1 = \
        elementsCountThroughWall + 1 + (elementsCountThroughWall + 1) * elementsCountAround * \
        (ureterElementPositionDown - (1 if ureter1Position.xi2 > 0.5 else 2)) + \
        ureterElementPositionAround + (1 if ureter1Position.xi1 > 0.5 else 0)
    ureter1StartCornerx = xBladder[endPointStartId1 - 1]
    v1 = [(ureter1StartCornerx[c] - centerUreter1_x[c]) for c in range(3)]
    ureter1Direction = vector.crossproduct3(td3, v1)
    nodeIdentifier, elementIdentifier, (o1_x, o1_d1, o1_d2, _, o1_NodeId, o1_Positions) = \
        generateOstiumMesh(region, ureterDefaultOptions, trackSurfaceUreter1, ureter1Position, ureter1Direction,
                           startNodeIdentifier=nextNodeIdentifier, startElementIdentifier=nextElementIdentifier,
                           vesselMeshGroups=[[ureterMeshGroup]], ostiumMeshGroups=bladderMeshGroup)

    # Ureter 2
    centerUreter2_x, centerUreter2_d1, centerUreter2_d2 = trackSurfaceUreter2.evaluateCoordinates(ureter2Position,
                                                                                                  derivatives=True)
    ad1, ad2, ad3 = calculate_surface_axes(centerUreter2_d1, centerUreter2_d2, [1.0, 0.0, 0.0])
    endPointStartId2 = elementsCountThroughWall + 1 + \
                       (elementsCountThroughWall + 1) * elementsCountAround * \
                       (ureterElementPositionDown - (1 if ureter1Position.xi2 > 0.5 else 2)) + \
                       elementsCountAround - ureterElementPositionAround + (-1 if ureter1Position.xi1 > 0.5 else 0)
    ureter2StartCornerx = xBladder[endPointStartId2 - 1]
    v2 = [(ureter2StartCornerx[c] - centerUreter2_x[c]) for c in range(3)]
    ureter2Direction = vector.crossproduct3(ad3, v2)
    nodeIdentifier, elementIdentifier, (o2_x, o2_d1, o2_d2, _, o2_NodeId, o2_Positions) = \
        generateOstiumMesh(region, ureterDefaultOptions, trackSurfaceUreter2, ureter2Position, ureter2Direction,
                           startNodeIdentifier=nodeIdentifier, startElementIdentifier=elementIdentifier,
                           vesselMeshGroups=[[ureterMeshGroup]], ostiumMeshGroups=bladderMeshGroup)

    # Create annulus mesh around ureters
    endPoints1_x = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]
    endPoints1_d1 = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]
    endPoints1_d2 = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]
    endNode1_Id = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]
    endDerivativesMap = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]
    endPoints2_x = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]
    endPoints2_d1 = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]
    endPoints2_d2 = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]
    endNode2_Id = [[None] * elementsCountAroundUreter for n3 in range(elementsCountThroughWall + 1)]

    count = 0
    for n2 in range(3):
        endNode1_Id[0][count] = endPointStartId1 + n2 * (elementsCountThroughWall + 1) * elementsCountAround
        endNode2_Id[0][count] = endPointStartId2 + n2 * (elementsCountThroughWall + 1) * elementsCountAround
        count += 1
    for n1 in range(2):
        endNode1_Id[0][count] = endNode1_Id[0][count - 1] + 1
        endNode2_Id[0][count] = endNode2_Id[0][count - 1] + 1
        count += 1
    for n2 in range(2):
        endNode1_Id[0][count] = endNode1_Id[0][count - 1] - (elementsCountThroughWall + 1) * elementsCountAround
        endNode2_Id[0][count] = endNode2_Id[0][count - 1] - (elementsCountThroughWall + 1) * elementsCountAround
        count += 1
    endNode1_Id[0][count] = endNode1_Id[0][count - 1] - 1
    endNode2_Id[0][count] = endNode2_Id[0][count - 1] - 1

    for n3 in range(1, elementsCountThroughWall + 1):
        for n in range(len(endNode1_Id[0])):
            endNode1_Id[n3][n] = endNode1_Id[0][n] + elementsCountAround * n3
            endNode2_Id[n3][n] = endNode2_Id[0][n] + elementsCountAround * n3

    for n3 in range(elementsCountThroughWall + 1):
        for n1 in range(elementsCountAroundUreter):
            nc1 = endNode1_Id[n3][n1] - 1
            endPoints1_x[n3][n1] = xBladder[nc1]
            endPoints1_d1[n3][n1] = d1Bladder[nc1]
            endPoints1_d2[n3][n1] = d2Bladder[nc1]
            nc2 = endNode2_Id[n3][n1] - 1
            endPoints2_x[n3][n1] = xBladder[nc2]
            endPoints2_d1[n3][n1] = d1Bladder[nc2]
            endPoints2_d2[n3][n1] = d2Bladder[nc2]

    for n3 in range(elementsCountThroughWall + 1):
        for n1 in range(elementsCountAroundUreter):
            if n1 == 0:
                endDerivativesMap[n3][n1] = ((-1, 0, 0), (-1, -1, 0), None, (0, 1, 0))
            elif n1 == 1:
                endDerivativesMap[n3][n1] = ((0, 1, 0), (-1, 0, 0), None)
            elif n1 == 2:
                endDerivativesMap[n3][n1] = ((0, 1, 0), (-1, 1, 0), None, (1, 0, 0))
            elif n1 == 3:
                endDerivativesMap[n3][n1] = ((1, 0, 0), (0, 1, 0), None)
            elif n1 == 4:
                endDerivativesMap[n3][n1] = ((1, 0, 0), (1, 1, 0), None, (0, -1, 0))
            elif n1 == 5:
                endDerivativesMap[n3][n1] = ((0, -1, 0), (1, 0, 0), None)
            elif n1 == 6:
                endDerivativesMap[n3][n1] = ((0, -1, 0), (1, -1, 0), None, (-1, 0, 0))
            else:
                endDerivativesMap[n3][n1] = ((-1, 0, 0), (0, -1, 0), None)

    startProportions1 = []
    for n in range(len(o1_Positions)):
        startProportions1.append(trackSurfaceUreter1.getProportion(o1_Positions[n]))

    startProportions2 = []
    for n in range(len(o2_Positions)):
        startProportions2.append(trackSurfaceUreter2.getProportion(o2_Positions[n]))

    endProportions1 = []
    elementsAroundTrackSurface1 = trackSurfaceUreter1.elementsCount1
    elementsAlongTrackSurface1 = trackSurfaceUreter1.elementsCount2

    endProportions2 = []
    elementsAroundTrackSurface2 = trackSurfaceUreter2.elementsCount1
    elementsAlongTrackSurface2 = trackSurfaceUreter2.elementsCount2

    firstIdxAround1 = ureterElementPositionAround + (0 if ureter1Position.xi1 > 0.5 else -1)
    firstIdxAlong = ureterElementPositionDown - (0 if ureter1Position.xi2 > 0.5 else 1)
    firstIdxAround2 = elementsCountAround//2 - ureterElementPositionAround - (2 if ureter1Position.xi1 > 0.5 else 1)

    for n in range(3):
        endProportions1.append([firstIdxAround1/elementsAroundTrackSurface1,
                                (firstIdxAlong + n)/elementsAlongTrackSurface1])
        endProportions2.append([firstIdxAround2 / elementsAroundTrackSurface2,
                                (firstIdxAlong + n) / elementsAlongTrackSurface2])
    for n in range(2):
        endProportions1.append([(firstIdxAround1 + n + 1) / elementsAroundTrackSurface1,
                                (firstIdxAlong + 2) / elementsAlongTrackSurface1])
        endProportions2.append([(firstIdxAround2 + n + 1) / elementsAroundTrackSurface2,
                                (firstIdxAlong + 2) / elementsAlongTrackSurface2])
    for n in range(2):
        endProportions1.append([(firstIdxAround1 + 2) / elementsAroundTrackSurface1,
                                (firstIdxAlong - n + 1) / elementsAlongTrackSurface1])
        endProportions2.append([(firstIdxAround2 + 2) / elementsAroundTrackSurface2,
                                (firstIdxAlong - n + 1) / elementsAlongTrackSurface2])
    endProportions1.append([(firstIdxAround1 + 1) / elementsAroundTrackSurface1,
                            firstIdxAlong / elementsAlongTrackSurface1])
    endProportions2.append([(firstIdxAround2 + 1) / elementsAroundTrackSurface2,
                            firstIdxAlong / elementsAlongTrackSurface2])

    nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
        nodes, mesh, nodeIdentifier, elementIdentifier,
        o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
        endPoints1_x, endPoints1_d1, endPoints1_d2, None, endNode1_Id, endDerivativesMap,
        elementsCountRadial=elementsCountUreterRadial, meshGroups=bladderMeshGroup,
        tracksurface=trackSurfaceUreter1, startProportions=startProportions1, endProportions=endProportions1,
        rescaleEndDerivatives=True)

    nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
        nodes, mesh, nodeIdentifier, elementIdentifier,
        o2_x, o2_d1, o2_d2, None, o2_NodeId, None,
        endPoints2_x, endPoints2_d1, endPoints2_d2, None, endNode2_Id, endDerivativesMap,
        elementsCountRadial=elementsCountUreterRadial, meshGroups=bladderMeshGroup,
        tracksurface=trackSurfaceUreter2, startProportions=startProportions2, endProportions=endProportions2,
        rescaleEndDerivatives=True)

    # Delete elements under annulus mesh
    element_identifiers = []
    elementToDeleteStartIdx1 = elementsCountThroughWall * elementsCountAround * \
                               (ureterElementPositionDown - (0 if ureter1Position.xi2 > 0.5 else 1)) + \
                               ureterElementPositionAround + (1 if ureter1Position.xi1 > 0.5 else 0)
    elementToDeleteStartIdx2 = elementsCountThroughWall * elementsCountAround * \
                               (ureterElementPositionDown - (0 if ureter1Position.xi2 > 0.5 else 1)) + \
                               elementsCountAround - ureterElementPositionAround + \
                               (-1 if ureter1Position.xi1 > 0.5 else 0)
    elementToDeleteStartIdxList = [elementToDeleteStartIdx1, elementToDeleteStartIdx2]
    for i in range(2):
        elementToDeleteStart = elementToDeleteStartIdxList[i]
        baseElementsToDelete = [elementToDeleteStart,
                                elementToDeleteStart + 1,
                                elementToDeleteStart + elementsCountAround * elementsCountThroughWall,
                                elementToDeleteStart + 1 + elementsCountAround * elementsCountThroughWall]
        for e3 in range(elementsCountThroughWall):
            elementsToDelete = [c + elementsCountAround * e3 for c in baseElementsToDelete]
            element_identifiers += elementsToDelete

    mesh_destroy_elements_and_nodes_by_identifiers(mesh, element_identifiers)

    return nodeIdentifier

def obtainBladderFlatNodes(elementsCountAlongBladder, elementsCountAround, elementsCountThroughWall,
                           xFinal, d1Final, d2Final, bladderLength, neckDiameter1, xApexInner, d2ApexInner,
                           bladderWallThickness):
    """
    Calculates flat coordinates for the bladder when it is opened into a flat preparation.
    :param elementsCountAlongBladder: Number of elements along bladder.
    :param elementsCountAround: Number of elements around bladder.
    :param elementsCountThroughWall: Number of elements through wall.
    :param xFinal, d1Final, d2Final: Coordinates and derivatives of all nodes of coordinates field.
    :param bladderLength: Bladder length along the central path.
    :param neckDiameter1: Major diameter of the bladder neck.
    :param xApexInner, d2ApexInner: Coordinates and d2 of the apex located on the lumen of the bladder.
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
    v2 = [0.5 * neckDiameter1, 0.0, bladderLength]
    alpha = vector.angleBetweenVectors(v1, v2)

    # Find apex to urethra arcLength in minor radius
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
        smd1 = smoothCubicHermiteDerivativesLine(xLineSmoothing, d1LineSmoothing, fixAllDirections=False,
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
        smd2 = smoothCubicHermiteDerivativesLine(lineSmoothingNodes, lineSmoothingNodes_d2, fixAllDirections=False,
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
    v1 = xApexInner
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
        xApex = [xApexInner[0], xApexInner[1], xApexInner[2] + n1 * bladderWallThickness / elementsCountThroughWall]
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
