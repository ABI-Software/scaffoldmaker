"""
Generates a 3-D cecum mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy
import math

from cmlibs.maths.vectorops import set_magnitude, normalize, magnitude, cross, dot, rotate_about_z_axis, axis_angle_to_rotation_matrix
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates
from cmlibs.zinc.element import Element
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, mergeAnnotationGroups
from scaffoldmaker.annotation.cecum_terms import get_cecum_term
from scaffoldmaker.annotation.smallintestine_terms import get_smallintestine_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import ColonSegmentTubeMeshOuterPoints, \
    getFullProfileFromHalfHaustrum, getTeniaColi, createNodesAndElementsTeniaColi
from scaffoldmaker.meshtypes.meshtype_3d_ostium2 import generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import setEftScaleFactorIds, remapEftNodeValueLabel
from scaffoldmaker.utils.tracksurface import TrackSurface
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, \
    mesh_destroy_elements_and_nodes_by_identifiers, get_nodeset_path_ordered_field_parameters


def getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName):
    assert parameterSetName in cls.getParameterSetNames()  # make sure parameter set is in list of parameters of parent scaffold
    if parameterSetName in ("Default", "Human 1"):
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                    "Structure": "1-2-3.2, 4-3-5"
                },
                'meshEdits': exnode_string_from_nodeset_field_parameters(
                    ['coordinates'],
                    [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                    (1, [[7.50,-20.00,0.00], [0.00,8.28,0.00], [-2.50,0.00,0.00], [0.00,0.00,0.00], [0.00,-0.00,2.50], [0.00,0.00,0.00]]),
                    (2, [[7.50,-10.86,0.00], [0.00,10.00,0.00], [-4.50,0.00,0.00], [0.00,0.00,0.00], [0.00,-0.00,4.50], [0.00,0.00,0.00]]),
                    (3, [[7.50,0.00,0.00], [[8.44,0.00,0.04],[0.00,11.72,0.00]], [[0.00,11.60,0.00],[-4.50,0.00,0.00]], [[1.02,6.79,0.00],[0.00,0.00,0.00]], [[0.00,0.00,11.60],[0.00,-0.00,4.50]], [[0.00,0.00,5.77],[0.00,0.00,0.00]]]),
                    (4, [[-1.88,0.00,-0.08], [10.32,0.00,0.12], [0.00,11.60,0.00], [0.00,0.00,0.00], [0.00,0.00,11.60], [0.00,0.00,0.00]]),
                    (5, [[15.00,0.00,0.00], [6.56,0.00,-0.04], [0.00,11.60,0.00], [0.00,0.00,0.00], [0.00,0.00,11.60], [0.00,0.00,0.00]])
                    ]]),

                'userAnnotationGroups': [
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '1-4',
                        'name': get_cecum_term('caecum')[0],
                        'ontId': get_cecum_term('caecum')[1]
                    },
                    {
                        '_AnnotationGroup': True,
                        'dimension': 1,
                        'identifierRanges': '1-2',
                        'name': get_cecum_term('ileum part of cecum')[0],
                        'ontId': get_cecum_term('ileum part of cecum')[1]
                    }]
            })
    elif "Human 2" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-3-5"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[-60.630,-80.530,895.250], [-6.170,-10.100,-5.870], [2.660,0.730,-4.050], [2.380,0.990,-1.140], [3.310,-2.970,1.640], [2.380,0.990,-1.140]]),
                (2, [[-68.640,-93.290,888.060], [-9.850,-15.420,-8.510], [3.880,1.200,-6.660], [2.380,0.990,-1.140], [2.820,-2.460,1.200], [2.380,0.990,-1.140]]),
                (3, [[-80.390,-111.370,878.290], [[-7.790,-0.980,12.360],[-13.650,-20.740,-11.030]], [[20.720,-4.040,12.540],[5.590,1.310,-9.380]], [[2.380,0.990,-1.140],[2.380,0.990,-1.140]], [[2.470,23.830,3.610],[1.510,-1.370,0.710]], [[2.380,0.990,-1.140],[2.380,0.990,-1.140]]]),
                (4, [[-71.690,-109.000,866.040], [-9.550,-3.730,12.060], [17.410,-4.850,11.460], [3.730,0.680,4.200], [0.820,19.940,7.200], [3.740,0.690,4.200]]),
                (5, [[-87.210,-111.060,890.540], [-4.750,0.410,12.390], [23.270,-3.130,7.880], [2.460,-0.390,-2.950], [3.090,24.460,0.450], [1.830,0.460,-4.310]])
                ]]),
                
            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-4',
                    'name': get_cecum_term('caecum')[0],
                    'ontId': get_cecum_term('caecum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_cecum_term('ileum part of cecum')[0],
                    'ontId': get_cecum_term('ileum part of cecum')[1]
                }]
        })
    elif "Pig 1" in parameterSetName:
        return ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-5-6-3-7"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                ['coordinates'],
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [[
                (1, [[164.00,29.75,51.53], [0.00,-9.63,-16.67], [-6.00,0.00,0.00], [0.00,0.00,0.00], [0.00,5.20,-3.00], [0.00,0.00,0.00]]),
                (2, [[164.00,17.50,30.31], [0.00,-14.87,-25.77], [-10.00,0.00,0.00], [0.00,0.00,0.00], [0.00,8.66,-5.00], [0.00,0.00,0.00]]),
                (3, [[164.00,0.00,0.00], [[30.00,0.00,0.00],[0.00,-20.13,-34.85]], [[0.00,37.00,0.00],[-10.00,0.00,0.00]], [[0.00,0.00,0.00],[0.00,0.00,0.00]], [[0.00,0.00,37.00],[0.00,8.66,-5.00]], [[0.00,0.00,0.00],[0.00,0.00,0.00]]]),
                (4, [[0.00,0.00,0.00], [60.00,0.00,0.00], [0.00,37.00,0.00], [0.00,0.00,0.00], [0.00,0.00,37.00], [0.00,0.00,0.00]]),
                (5, [[60.00,0.00,0.00], [60.00,0.00,0.00], [0.00,37.00,0.00], [0.00,0.00,0.00], [0.00,0.00,37.00], [0.00,0.00,0.00]]),
                (6, [[120.00,0.00,0.00], [52.00,0.00,0.00], [0.00,37.00,0.00], [0.00,0.00,0.00], [0.00,0.00,37.00], [0.00,0.00,0.00]]),
                (7, [[180.00,0.00,0.00], [2.00,0.00,0.00], [0.00,37.00,0.00], [0.00,0.00,0.00], [0.00,0.00,37.00], [0.00,0.00,0.00]])
                ]]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-6',
                    'name': get_cecum_term('caecum')[0],
                    'ontId': get_cecum_term('caecum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_cecum_term('ileum part of cecum')[0],
                    'ontId': get_cecum_term('ileum part of cecum')[1]
                }]
        })

def getDefaultOstiumSettings():
    """
    Generate list of default options for ostium.
    """
    options = { 'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium wall thickness': 1.6,
                'Ostium wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Use linear through ostium wall': True,
                'Vessel wall thickness': 1.6,
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
    ostiumOptions['Ostium wall thickness'] = options['Wall thickness']
    ostiumOptions['Vessel wall thickness'] = options['Ileum wall thickness']
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
           

class MeshType_3d_cecum1(Scaffold_base):
    '''
    Generates a 3-D cecum mesh with variable numbers of elements around, along the central line, and through wall. The
    cecum is created by a function that generates a cecum segment and uses tubemesh to map the segment along a network
    layout. The proximal end of the cecum is closed up with an apex plate. An ostium is included to generate the
    ileo-cecal junction.
    '''

    @staticmethod
    def getName():
        return '3D Cecum 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Human 2',
            'Pig 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {
            'Network layout': getDefaultNetworkLayoutScaffoldPackage(cls, parameterSetName),
            'Number of segments': 1,
            'Number of elements around tenia coli': 2,
            'Number of elements around haustrum': 8,
            'Number of elements along segment': 8,
            'Number of elements through wall': 1,
            'Corner outer radius factor': 0.536,
            'Haustrum outer radius factor': 0.464,
            'Segment length end derivative factor': 0.5,
            'Segment length mid derivative factor': 1.0,
            'Number of tenia coli': 3,
            'Start tenia coli width': 2.0,
            'Start tenia coli width derivative': 3.3,
            'End tenia coli width': 5.3,
            'End tenia coli width derivative': 3.3,
            'Tenia coli thickness': 0.5,
            'Wall thickness': 1.6,
            'Ileum wall thickness': 1.6,
            'Mucosa relative thickness': 0.55,
            'Submucosa relative thickness': 0.15,
            'Circular muscle layer relative thickness': 0.25,
            'Longitudinal muscle layer relative thickness': 0.05,
            'Use cross derivatives': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }

        if 'Human 2' in parameterSetName:
            options['Number of elements along segment'] = 12
            options['Segment length mid derivative factor'] = 3.0
            options['Start tenia coli width'] = 10.0
            options['Start tenia coli width derivative'] = 0.0
            options['End tenia coli width'] = 10.0
            options['End tenia coli width derivative'] = 0.0
            options[ 'Tenia coli thickness'] = 1.6

        elif 'Pig 1' in parameterSetName:
            options['Number of segments'] = 5
            options['Number of elements around haustrum'] = 12
            options['Haustrum outer radius factor'] = 0.25
            options['Segment length end derivative factor'] = 0.5
            options['Segment length mid derivative factor'] = 4.0
            options['Start tenia coli width'] = 5.0
            options['Start tenia coli width derivative'] = 0.0
            options['End tenia coli width'] = 5.0
            options['End tenia coli width derivative'] = 0.0
            options['Wall thickness'] = 2.0
            options['Ileum wall thickness'] = 2.0

        options['Base parameter set'] = parameterSetName
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Network layout',
            'Number of segments',
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Number of elements along segment',
            'Number of elements through wall',
            'Corner outer radius factor',
            'Haustrum outer radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Number of tenia coli',
            'Start tenia coli width',
            'Start tenia coli width derivative',
            'End tenia coli width',
            'End tenia coli width derivative',
            'Tenia coli thickness',
            'Wall thickness',
            'Ileum wall thickness',
            'Mucosa relative thickness',
            'Submucosa relative thickness',
            'Circular muscle layer relative thickness',
            'Longitudinal muscle layer relative thickness',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Network layout':
            return [ MeshType_1d_network_layout1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Network layout':
            return cls.getParameterSetNames()
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
            'Number of segments',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Corner outer radius factor',
            'Haustrum outer radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Number of tenia coli',
            'Start tenia coli width',
            'Start tenia coli width derivative',
            'End tenia coli width',
            'End tenia coli width derivative',
            'Tenia coli thickness',
            'Wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
            if options['Number of elements along segment'] < 8:
                options['Number of elements along segment'] = 8
            if options['Number of elements along segment'] % 4 > 0:
               options['Number of elements along segment'] = options['Number of elements along segment'] // 4 * 4
            if options['Number of elements through wall'] != (1 or 4):
                options['Number of elements through wall'] = 4

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup, None
        """

        nextNodeIdentifier = 1
        nextElementIdentifier = 1
        geometricNetworkLayout = options['Network layout']
        cecumTermsAlong = ['caecum', 'ileum part of cecum']
        geometricNetworkLayout = CecumNetworkLayout(region, geometricNetworkLayout, cecumTermsAlong)
        annotationGroups, nextNodeIdentifier, nextElementIdentifier = \
            createCecumMesh3d(region, options, geometricNetworkLayout, nextNodeIdentifier,
                              nextElementIdentifier)[0:3]

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

def getApexSegmentForCecum(xOuter, d1Outer, d2Outer, elementsCountAroundHalfHaustrum,
                           elementsCountAroundTC, elementsCountAround, elementsCountAlongSegment, tcCount):
    """
    Generates the outer coordinates and derivatives for a cecum segment on the closed end.
    The closed end is a single node and segment is created by sampling curves between the
    point on closed end with nodes on the length along second half of a colon segment.
    :param xInner: coordinates of a colon segment.
    :param d1Inner: derivative around colon segment.
    :param d2Inner: derivative along colon segment.
    :param elementsCountAroundHalfHaustrum: half of total number of elements in haustrum and tenia coli.
    :param elementsCountAroundTC: number of elements around a tenia coli.
    :param elementsCountAround: number of elements around a cecum.
    :param elementsCountAlongSegment: number of elements along a segment of cecum.
    :param tcCount: number of tenia coli.
    :return: coordinates and derivatives around and along closed segment of cecum, and directional derivative for node
    on middle of tenia coli.
    """

    # Make apex cap - apex points like multiple points
    xFirstSegment = [[0.0, 0.0, 0.0] for c in range(elementsCountAround)]

    # Compile nodes and d2 for sampling
    xFirstSegment += xOuter[elementsCountAround * int(elementsCountAlongSegment * 0.5):] # second half of first regular segment
    d1FirstDirectionVector = normalize(d1Outer[elementsCountAround]) # Store direction vector of first d1 intra-haustral for later
    d2Vector = xOuter[elementsCountAround * int(elementsCountAlongSegment * 0.5):
                      elementsCountAround * (int(elementsCountAlongSegment * 0.5) + 1)] # half face of segment - apex
    d2FirstSegment = []
    for c in range(elementsCountAround):
        d2 = [d2Vector[c][0], d2Vector[c][1], 0.0 ]  # project onto x-y plane to get d2 pointing vertically
        d2FirstSegment.append(d2)
    d2FirstSegment += d2Outer[elementsCountAround * int(elementsCountAlongSegment*0.5):]

    # Sample along first segment
    xFirstSegmentSampledRaw = []
    d2FirstSegmentSampledRaw = []
    xFirstSegmentSampled = []
    d1FirstSegmentSampled = []
    d2FirstSegmentSampled = []

    for n1 in range(elementsCountAround):
        xForSamplingAlong = []
        d2ForSamplingAlong = []
        for n2 in range(1 + elementsCountAlongSegment - int(elementsCountAlongSegment * 0.5) + 1):
            idx = elementsCountAround * n2 + n1
            xForSamplingAlong.append(xFirstSegment[idx])
            d2ForSamplingAlong.append(d2FirstSegment[idx])
        xResampled, d1Resampled, se, sxi, _ = interp.sampleCubicHermiteCurves(xForSamplingAlong, d2ForSamplingAlong,
                                                                              elementsCountAlongSegment,
                                                                              arcLengthDerivatives = True)

        xFirstSegmentSampledRaw.append(xResampled)
        d2FirstSegmentSampledRaw.append(d1Resampled)

    # Re-arrange sample order
    for n2 in range(elementsCountAlongSegment + 1):
        xAround = []
        for n1 in range(elementsCountAround):
            x = xFirstSegmentSampledRaw[n1][n2]
            d2 = d2FirstSegmentSampledRaw[n1][n2]
            if n1 < elementsCountAroundHalfHaustrum + 1:
                xAround.append(x)
            xFirstSegmentSampled.append(x)
            d2FirstSegmentSampled.append(d2)
            if n2 == 0:
                d1 = rotate_about_z_axis(d2, math.pi*0.5)
                d1FirstSegmentSampled.append(d1)

        if n2 > 0:
            d1Around = []
            for n1 in range(elementsCountAroundHalfHaustrum):
                v1 = xAround[n1]
                v2 = xAround[n1 + 1]
                d1 = d1FirstDirectionVector if n1 == 0 else [v2[c] - v1[c] for c in range(3)]
                d2 = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c*arcLengthAround for c in normalize(d1)]
                d1Around.append(dx_ds1)
            # Account for d1 of node sitting on half haustrum
            d1 = normalize(
                [xAround[elementsCountAroundHalfHaustrum][c] - xAround[elementsCountAroundHalfHaustrum - 1][c]
                 for c in range(3)])
            dx_ds1 = [c * arcLengthAround for c in d1]
            d1Around.append(dx_ds1)

            d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, d1Around, fixStartDerivative=True)
            d1TCEdge = set_magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5)],
                                           magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 - 1)]))
            d1Transition = set_magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 1)],
                                               magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 2)]))
            d1Corrected = []
            d1Corrected = d1Corrected + d1Smoothed[:int(elementsCountAroundTC * 0.5)]
            d1Corrected.append(d1TCEdge)
            d1Corrected.append(d1Transition)
            d1Corrected = d1Corrected + d1Smoothed[int(elementsCountAroundTC * 0.5 + 2):]
            d1Full = getD1ForFullProfileFromHalfHaustrum(d1Corrected, tcCount)
            d1FirstSegmentSampled += d1Full

    return xFirstSegmentSampled, d1FirstSegmentSampled, d2FirstSegmentSampled, d1FirstDirectionVector

def getD1ForFullProfileFromHalfHaustrum(d1HaustrumHalfSet, tcCount):
    """
    Get full profile from half haustrum
    :param d1HaustrumHalfSet: d1 for first half of haustrum
    :param tcCount: Number of tenia coli
    :return: d1 for nodes around the entire haustrum
    """
    d1HaustrumHalfSet2 = []
    d1Haustra = []

    rotAng = 2 * math.pi / tcCount
    for n in range(1, len(d1HaustrumHalfSet)):
        idx = -n + len(d1HaustrumHalfSet) - 1
        d1 = d1HaustrumHalfSet[idx]
        d1Reflect = [d1[0], -d1[1], d1[2]]
        d1Rot = [-(d1Reflect[0] * math.cos(rotAng) - d1Reflect[1] * math.sin(rotAng)),
                 -(d1Reflect[0] * math.sin(rotAng) + d1Reflect[1] * math.cos(rotAng)),
                 -d1Reflect[2]]
        d1HaustrumHalfSet2.append(d1Rot)

    d1Haustrum = d1HaustrumHalfSet + d1HaustrumHalfSet2

    # Rotate to get all 3 sectors
    d1Haustra = d1Haustra + d1Haustrum[:-1]
    ang = [2 / 3 * math.pi, -2 / 3 * math.pi] if tcCount == 3 else [math.pi]
    for i in range(tcCount - 1):
        rotAng = ang[i]
        cosRotAng = math.cos(rotAng)
        sinRotAng = math.sin(rotAng)
        for n in range(len(d1Haustrum) - 1):
            d1 = d1Haustrum[n]
            dx_ds1 = [d1[0] * cosRotAng - d1[1] * sinRotAng, d1[0] * sinRotAng + d1[1] * cosRotAng, d1[2]]
            d1Haustra.append(dx_ds1)

    return d1Haustra

def sampleCecumAlongLength(xToSample, d1ToSample, d2ToSample, d1FirstDirectionVector, elementsCountAroundHalfHaustrum,
                           elementsCountAroundTC, elementsCountAround, elementsCountAlong, tcCount):
    """
    Get systematically spaced points and derivatives over cubic Hermite interpolated curves along the
    length of the cecum.
    :param xToSample: coordinates of nodes.
    :param d1ToSample: derivative around elements.
    :param d2ToSample: derivative along cecum length.
    :param d1FirstDirectionVector: directional vector of derivative around for node on the middle of the tenia coli.
    :param elementsCountAroundHalfHaustrum:half the total number of elements around tenia coli and haustrum.
    :param elementsCountAroundTC: number of elements around tenia coli.
    :param elementsCountAround: number of elements around cecum.
    :param elementsCountAlong: number of elements along cecum length.
    :param tcCount: number of tenia coli.
    :return: nodes and derivatives for equally spaced points.
    """

    xOuterRaw = []
    d2OuterRaw = []
    xSampledAlongLength = []
    d1SampledAlongLength = []
    d2SampledAlongLength = []

    for n1 in range(elementsCountAroundHalfHaustrum + 1):
        xForSamplingAlong = []
        d2ForSamplingAlong = []
        for n2 in range(elementsCountAlong + 1):
            idx = n2 * elementsCountAround + n1
            xForSamplingAlong.append(xToSample[idx])
            d2ForSamplingAlong.append(d2ToSample[idx])
        xSampled, d2Sampled, se, sxi, _ = interp.sampleCubicHermiteCurves(xForSamplingAlong, d2ForSamplingAlong,
                                                                          elementsCountAlong,
                                                                          arcLengthDerivatives=True)
        xOuterRaw.append(xSampled)
        d2OuterRaw.append(d2Sampled)

    # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
    for n2 in range(elementsCountAlong + 1):
        xAround = []
        d2Around = []

        for n1 in range(elementsCountAroundHalfHaustrum + 1):
            x = xOuterRaw[n1][n2]
            d2 = d2OuterRaw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)

        d1OuterAroundList = []
        if n2 == 0:
            d1Corrected = d1ToSample[:elementsCountAroundHalfHaustrum + 1]

        else:
            for n1 in range(elementsCountAroundHalfHaustrum):
                v1 = xAround[n1]
                v2 = xAround[n1 + 1]
                d1 = d1FirstDirectionVector if n1 == 0 else [v2[c] - v1[c] for c in range(3)]
                d2 = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c * arcLengthAround for c in normalize(d1)]
                d1OuterAroundList.append(dx_ds1)
            # Account for d1 of node sitting on half haustrum
            d1 = normalize([xAround[elementsCountAroundHalfHaustrum][c] -
                                   xAround[elementsCountAroundHalfHaustrum - 1][c] for c in range(3)])
            dx_ds1 = [c * arcLengthAround for c in d1]
            d1OuterAroundList.append(dx_ds1)

        if d1OuterAroundList:
            d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, d1OuterAroundList, fixStartDerivative=True)
            d1TCEdge = set_magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5)],
                                           magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 - 1)]))
            d1Transition = set_magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 1)],
                                               magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 2)]))
            d1Corrected = []
            d1Corrected = d1Corrected + d1Smoothed[:int(elementsCountAroundTC * 0.5)]
            d1Corrected.append(d1TCEdge)
            d1Corrected.append(d1Transition)
            d1Corrected = d1Corrected + d1Smoothed[int(elementsCountAroundTC * 0.5 + 2):]

        xAlongList, d1AlongList, d2AlongList = getFullProfileFromHalfHaustrum(xAround, d1Corrected, d2Around, tcCount)

        xSampledAlongLength += xAlongList
        d1SampledAlongLength += d1AlongList
        d2SampledAlongLength += d2AlongList

    return xSampledAlongLength, d1SampledAlongLength, d2SampledAlongLength

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
    d = [c * arcLengthAround for c in normalize(d)]

    return d

def createCecumMesh3d(region, options, networkLayout, nodeIdentifier, elementIdentifier, nodeIdProximalIleum=[],
                      xProximalIleum=[], d1ProximalIleum=[], d2ProximalIleum=[], d3ProximalIleum=[]):
    """
    Generates a cecum scaffold in the region using a network layout and parameter options.
    :param region: Region to create elements in.
    :param options: Parameter options for scaffold.
    :param networkLayout: Network layout through the axis of the cecum scaffold.
    :param nodeIdentifier: First node identifier.
    :param elementIdentifier: First element identifier.
    :param nodeIdProximalIleum: Node identifiers of nodes around starting nodes for ileum.
    :param xProximalIleum, d1ProximalIleum, d2ProximalIleum, d3ProximalIleum: coordinates and derivatives of nodes
    around starting nodes for ileum.
    :return allAnnotationGroups, nextNodeIdentifier, nextElementIdentifier.
    """
    parameterSetName = options['Base parameter set']
    isHuman = 'Human' in parameterSetName

    segmentCount = options['Number of segments']
    startPhase = 0.0
    elementsCountAroundTC = options['Number of elements around tenia coli']
    elementsCountAroundHaustrum = options['Number of elements around haustrum']
    elementsCountAlongSegment = options['Number of elements along segment']
    elementsCountThroughWall = options['Number of elements through wall']
    cornerOuterRadiusFactor = options['Corner outer radius factor']
    haustrumOuterRadiusFactor = options['Haustrum outer radius factor']
    segmentLengthEndDerivativeFactor = options['Segment length end derivative factor']
    segmentLengthMidDerivativeFactor = options['Segment length mid derivative factor']
    tcCount = options['Number of tenia coli']
    startTCWidth = options['Start tenia coli width']
    startTCWidthDerivative = options['Start tenia coli width derivative']
    endTCWidth = options['End tenia coli width']
    endTCWidthDerivative = options['End tenia coli width derivative']
    tcThickness = options['Tenia coli thickness']
    wallThickness = options['Wall thickness']
    mucosaRelThickness = options['Mucosa relative thickness']
    submucosaRelThickness = options['Submucosa relative thickness']
    circularRelThickness = options['Circular muscle layer relative thickness']
    longitudinalRelThickness = options['Longitudinal muscle layer relative thickness']
    useCrossDerivatives = options['Use cross derivatives']
    useCubicHermiteThroughWall = not (options['Use linear through wall'])
    elementsCountAlong = int(elementsCountAlongSegment * segmentCount)
    elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum) * tcCount

    ostiumOptions = getDefaultOstiumSettings()
    ostiumSettings = updateOstiumOptions(options, ostiumOptions)

    zero = [0.0, 0.0, 0.0]
    firstNodeIdentifier = nodeIdentifier
    firstElementIdentifier = elementIdentifier
    startNode = nodeIdentifier
    startElement = elementIdentifier

    fm = region.getFieldmodule()
    coordinates = findOrCreateFieldCoordinates(fm)
    cache = fm.createFieldcache()
    mesh = fm.findMeshByDimension(3)

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

    if elementsCountThroughWall == 1:
        relativeThicknessList = [1.0]
        annotationGroupsThroughWall = [[]]
    else:
        relativeThicknessList = [mucosaRelThickness, submucosaRelThickness,
                                 circularRelThickness, longitudinalRelThickness]
        longitudinalMuscleGroup = AnnotationGroup(region, get_cecum_term("longitudinal muscle layer of cecum"))
        circularMuscleGroup = AnnotationGroup(region, get_cecum_term("circular muscle layer of cecum"))
        submucosaGroup = AnnotationGroup(region, get_cecum_term("submucosa of cecum"))
        mucosaGroup = AnnotationGroup(region, get_cecum_term("cecum mucosa"))
        annotationGroupsThroughWall = [[mucosaGroup], [submucosaGroup], [circularMuscleGroup],
                                       [longitudinalMuscleGroup]]

    # Sample network layout along cecum
    # print(len(networkLayout.cxGroups))
    cecumLength = networkLayout.arcLengthOfGroupsAlong[0]
    cx = networkLayout.cxGroups[0]
    cd1 = networkLayout.cd1Groups[0]
    cd2 = networkLayout.cd2Groups[0]
    cd12 = networkLayout.cd12Groups[0]

    cxIleum = networkLayout.cxGroups[1]
    cd1Ileum = networkLayout.cd1Groups[1]
    cd2Ileum = networkLayout.cd2Groups[1]
    cd3Ileum = networkLayout.cd3Groups[1]
    cd12Ileum = networkLayout.cd12Groups[1]

    d2BranchPt = networkLayout.d2BranchPt
    d3BranchPt = networkLayout.d3BranchPt

    sxCecum, sd1Cecum, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
    sd2Cecum, sd12Cecum = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

    # Calculate segment length
    segmentLength = cecumLength / segmentCount

    # Generate variation of radius & tc width along length
    outerRadiusAlongCecum = []
    dOuterRadiusAlongCecum = []
    tcWidthAlongCecum = []

    closedProximalEnd = True
    outerRadiusListCP = [magnitude(c) for c in cd2]
    dOuterRadiusListCP = []
    for n in range(len(outerRadiusListCP) - 1):
        dOuterRadiusListCP.append(outerRadiusListCP[n + 1] - outerRadiusListCP[n])
    dOuterRadiusListCP.append(outerRadiusListCP[-1] - outerRadiusListCP[-2])
    outerRadiusAlongElementList, dOuterRadiusAlongElementList = \
        interp.interpolateSampleCubicHermite(outerRadiusListCP, dOuterRadiusListCP, se, sxi, ssf)

    for n2 in range(elementsCountAlongSegment * segmentCount + 1):
        xi = 1 / (elementsCountAlongSegment * segmentCount) * n2

        radius = outerRadiusAlongElementList[n2]
        outerRadiusAlongCecum.append(radius)
        dRadius = dOuterRadiusAlongElementList[n2]
        dOuterRadiusAlongCecum.append(dRadius)
        tcWidth = interp.interpolateCubicHermite([startTCWidth], [startTCWidthDerivative],
                                                 [endTCWidth], [endTCWidthDerivative], xi)[0]
        tcWidthAlongCecum.append(tcWidth)

    haustrumOuterRadiusFactorAlongCecum = [haustrumOuterRadiusFactor] * (elementsCountAlong + 1)

    xToSample = []
    d1ToSample = []
    d2ToSample = []

    elementsCountAroundHalfHaustrum = int((elementsCountAroundTC + elementsCountAroundHaustrum) * 0.5)

    # Create object
    colonSegmentTubeMeshOuterPoints = ColonSegmentTubeMeshOuterPoints(
        region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
        tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
        segmentLength, wallThickness, cornerOuterRadiusFactor, haustrumOuterRadiusFactorAlongCecum,
        outerRadiusAlongCecum, dOuterRadiusAlongCecum, tcWidthAlongCecum, startPhase)

    # Create annotation
    cecumGroup = AnnotationGroup(region, get_cecum_term("caecum"))
    annotationGroupsAlong = []
    for i in range(elementsCountAlong):
        annotationGroupsAlong.append([cecumGroup])

    for nSegment in range(segmentCount):
        # Make regular segments
        xOuter, d1Outer, d2Outer, transitElementList, segmentAxis, annotationGroupsAround \
            = colonSegmentTubeMeshOuterPoints.getColonSegmentTubeMeshOuterPoints(nSegment)

        # Replace first half of first segment with apex and sample along apex and second half of segment
        if nSegment == 0:
            xFirstSegmentSampled, d1FirstSegmentSampled, d2FirstSegmentSampled, d1FirstDirectionVector = \
                getApexSegmentForCecum(xOuter, d1Outer, d2Outer, elementsCountAroundHalfHaustrum,
                                       elementsCountAroundTC, elementsCountAround, elementsCountAlongSegment,
                                       tcCount)

            xToSample += xFirstSegmentSampled
            d1ToSample += d1FirstSegmentSampled
            d2ToSample += d2FirstSegmentSampled
        else:
            xOuterExtrude = []
            for n in range(len(xOuter)):
                xOuterExtrude.append([xOuter[n][0], xOuter[n][1], xOuter[n][2] + segmentLength * nSegment])
            xToSample += xOuterExtrude[elementsCountAround:]
            d1ToSample += d1Outer[elementsCountAround:]
            d2ToSample += d2Outer[elementsCountAround:]

    # Sample along length
    xToWarp, d1ToWarp, d2ToWarp = sampleCecumAlongLength(xToSample, d1ToSample, d2ToSample, d1FirstDirectionVector,
                                                         elementsCountAroundHalfHaustrum, elementsCountAroundTC,
                                                         elementsCountAround, elementsCountAlong, tcCount)

    # Ensure cecum starts at z = 0.0
    minZ = xToWarp[0][2]
    for n2 in range(elementsCountAlong + 1):
        zFirstNodeAlong = xToWarp[n2 * elementsCountAround][2]
        if zFirstNodeAlong < minZ:
            minZ = zFirstNodeAlong

    for n in range(len(xToWarp)):
        xToWarp[n][2] = xToWarp[n][2] - minZ

    # Project reference point for warping onto network layout
    sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
        tubemesh.getPlaneProjectionOnCentralPath(xToWarp, elementsCountAround, elementsCountAlong,
                                                 cecumLength, sxCecum, sd1Cecum, sd2Cecum, sd12Cecum)

    # Warp points
    xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = \
        tubemesh.warpSegmentPoints(xToWarp, d1ToWarp, d2ToWarp, segmentAxis, sxRefList, sd1RefList,
                                   sd2ProjectedListRef, elementsCountAround, elementsCountAlong, zRefList)

    # Create coordinates and derivatives
    wallThicknessList = [wallThickness] * (elementsCountAlong + 1)

    xList, d1List, d2List, d3List, curvatureList, localIdxDistal, xDistal, d1Distal, d2Distal, d3Distal = \
        tubemesh.extrudeSurfaceCoordinates(xWarpedList, d1WarpedList,d2WarpedList, d3WarpedUnitList,
                                           wallThicknessList, relativeThicknessList, elementsCountAround,
                                           elementsCountAlong, elementsCountThroughWall, transitElementList,
                                           outward=False)

    # Deal with multiple nodes at end point for closed proximal end
    xApexInner = xList[0]
    # arclength between apex point and corresponding point on next face
    magMin = interp.computeCubicHermiteArcLength(xList[0], d2List[0],
                                              xList[elementsCountAround * (elementsCountThroughWall + 1)],
                                              d2List[elementsCountAround * (elementsCountThroughWall + 1)],
                                              rescaleDerivatives=True)
    magMax = interp.computeCubicHermiteArcLength(xList[int(0.5*(elementsCountAroundTC + elementsCountAroundHaustrum))],
                                                 d2List[int(0.5*(elementsCountAroundTC + elementsCountAroundHaustrum))],
                                                 xList[int(0.5*(elementsCountAroundTC + elementsCountAroundHaustrum)) +
                                                       elementsCountAround * (elementsCountThroughWall + 1)],
                                                 d2List[int(0.5*(elementsCountAroundTC + elementsCountAroundHaustrum)) +
                                                        elementsCountAround * (elementsCountThroughWall + 1)],
                                                 rescaleDerivatives=True)
    mag = 0.5*(magMin + magMax)
    d2ApexInner = set_magnitude(sd2Cecum[0], mag)
    d1ApexInner = cross(sd1Cecum[0], d2ApexInner)
    d1ApexInner = set_magnitude(d1ApexInner, mag)
    d3ApexUnit = normalize(cross(normalize(d1ApexInner),
                                                       normalize(d2ApexInner)))
    d3ApexInner = [d3ApexUnit[c] * wallThickness / elementsCountThroughWall for c in range(3)]

    xCecum = []
    d1Cecum = []
    d2Cecum = []
    d3Cecum = []

    for n3 in range(elementsCountThroughWall + 1):
        xApex = [xApexInner[c] + d3ApexUnit[c] * wallThickness / elementsCountThroughWall * n3 for c in range(3)]
        xCecum.append(xApex)
        d1Cecum.append(d1ApexInner)
        d2Cecum.append(d2ApexInner)
        d3Cecum.append(d3ApexInner)

    xCecum += xList[(elementsCountThroughWall + 1) * elementsCountAround:]
    d1Cecum += d1List[(elementsCountThroughWall + 1) * elementsCountAround:]
    d2Cecum += d2List[(elementsCountThroughWall + 1) * elementsCountAround:]
    d3Cecum += d3List[(elementsCountThroughWall + 1) * elementsCountAround:]

    xFlat = d1Flat = d2Flat = []
    xOrgan = d1Organ = d2Organ = []

    # Create nodes and elements
    if tcThickness > 0:
        tubeTCWidthList = colonSegmentTubeMeshOuterPoints.getTubeTCWidthList()
        xCecum, d1Cecum, d2Cecum, d3Cecum, annotationArrayAround, localIdxDistal, xDistal, d1Distal, d2Distal, \
        d3Distal = \
            getTeniaColi(region, xCecum, d1Cecum, d2Cecum, d3Cecum, curvatureList, tcCount, elementsCountAroundTC,
                         elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
                         tubeTCWidthList, tcThickness, annotationGroupsAround, closedProximalEnd, isHuman)

        nextNodeIdentifier, nextElementIdentifier, allAnnotationGroups, nodesIdDistal = createNodesAndElementsTeniaColi(
            region, xCecum, d1Cecum, d2Cecum, d3Cecum, xFlat, d1Flat, d2Flat, xOrgan, d1Organ, d2Organ, None,
            elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
            tcCount, annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
            closedProximalEnd, localIdxDistal)

    else:
        nextNodeIdentifier, nextElementIdentifier, allAnnotationGroups, nodesIdDistal = tubemesh.createNodesAndElements(
            region, xCecum, d1Cecum, d2Cecum, d3Cecum, xFlat, d1Flat, d2Flat, xOrgan, d1Organ, d2Organ, None,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall,
            annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
            closedProximalEnd)

    # nodeIdentifierCecum = nextNodeIdentifier
    # for n2 in range(len(cx)):
    #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    #     nextNodeIdentifier += 1
    #
    # #################
    # # Create elements
    # #################
    #
    # mesh = fm.findMeshByDimension(1)
    # cubicHermiteBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    # eft = mesh.createElementfieldtemplate(cubicHermiteBasis)
    # elementtemplate = mesh.createElementtemplate()
    # elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
    # result = elementtemplate.defineField(coordinates, -1, eft)
    #
    # elementIdentifier = nextElementIdentifier
    # for e in range(len(cx) - 1):
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     element.setNodesByIdentifier(eft, [nodeIdentifierCecum + e, nodeIdentifierCecum + 1 + e])
    #     elementIdentifier = elementIdentifier + 1
    #
    # cx = networkLayout.cxGroups[1]
    # cd1 = networkLayout.cd1Groups[1]
    # cd2 = networkLayout.cd2Groups[1]
    #
    # nodeIdentifierIleum = nextNodeIdentifier
    # for n2 in range(len(cx)):
    #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n2])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    #     nextNodeIdentifier += 1
    #
    # for e in range(2):
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     element.setNodesByIdentifier(eft, [nodeIdentifierIleum + e, nodeIdentifierIleum + 1 + e])
    #     elementIdentifier = elementIdentifier + 1

    # Add ostium on track surface
    elementsAroundTrackSurface = elementsCountAroundHaustrum
    elementsAlongTrackSurface = elementsCountAlongSegment

    # Find region where ostium sits
    # angle between d2 of branch point and vector between branch point and 1st point on ileum
    dV = [cxIleum[0][c] - cxIleum[-1][c] for c in range(3)]
    ostiumPositionAngleAround = math.acos(dot(dV, d2BranchPt)/
                                          (magnitude(dV) * magnitude(d2BranchPt)))
    if dot(dV,d3BranchPt) != 0:
        sign = dot(dV, d3BranchPt)/abs(dot(dV, d3BranchPt))
        if sign < 0:
            ostiumPositionAngleAround = math.pi * 2.0 - ostiumPositionAngleAround
    sectorIdx = ostiumPositionAngleAround // (2 * math.pi / tcCount)
    sectorStartAngle = sectorIdx * (2 * math.pi / tcCount)

    startIdxElementsAround = int((elementsCountAroundHaustrum + elementsCountAroundTC) * sectorIdx +
                                 elementsCountAroundTC * 0.5)

    segmentIdx = int(networkLayout.arcLengthToBranchPt // segmentLength)

    baseNodesIdx = (elementsCountThroughWall + 1) + \
                   + (elementsCountAround * (elementsCountThroughWall + 1) +
                      ((elementsCountAroundTC - 1) * tcCount if tcThickness > 0.0 else 0)) * \
                   (elementsCountAlongSegment * segmentIdx - 1) + elementsCountAround

    # Elements to delete
    deleteElementIdentifier = []
    for n2 in range(elementsCountAlongSegment):
        for n3 in range(elementsCountThroughWall):
            for n1 in range(elementsCountAroundHaustrum):
                elementIdx = \
                    startIdxElementsAround + int(elementsCountAroundTC * (0.5 if tcThickness > 0.0 else 0)) + n1 + \
                    (elementsCountAround * n3) + n2 * (elementsCountAround * elementsCountThroughWall +
                                                       elementsCountAroundTC * (tcCount if tcThickness > 0.0 else 0))\
                    + startElement - 1 + segmentIdx * (elementsCountAlongSegment *
                                                       (elementsCountAround * elementsCountThroughWall +
                                                        elementsCountAroundTC * (tcCount if tcThickness > 0.0 else 0)))
                deleteElementIdentifier.append(elementIdx)

    xTrackSurface = []
    d1TrackSurface = []
    d2TrackSurface = []
    nodesOnLHS = []
    nodesOnRHS = []
    nodesOnDistal = []
    nodesOnProximal = []

    for n2 in range(elementsCountAlongSegment + 1):
        for n1 in range(elementsCountAroundHaustrum + 1):
            if n2 == 0 and segmentIdx == 0:
                xTrackSurface.append(xApex)
                d1TrackSurface.append(d1ApexInner)
                d2TrackSurface.append(d2ApexInner)
            else:
                idx = baseNodesIdx + \
                      (elementsCountAround * (elementsCountThroughWall + 1) +
                       ((elementsCountAroundTC - 1) * tcCount if tcThickness > 0.0 else 0)) * n2 + \
                      startIdxElementsAround + n1 + (elementsCountAround * (elementsCountThroughWall - 1))
                xTrackSurface.append(xCecum[idx])
                d1TrackSurface.append(d1Cecum[idx])
                d2TrackSurface.append(d2Cecum[idx])

                if n1 == 1:
                    nodeWall = []
                    for n3 in range(elementsCountThroughWall, -1, -1):
                        nodeWall.append(idx - elementsCountAround * n3)
                    nodesOnLHS.append(nodeWall)
                if n1 == elementsCountAroundHaustrum:
                    nodeWall = []
                    for n3 in range(elementsCountThroughWall, -1, -1):
                        nodeWall.append(idx + 1 - elementsCountAround * n3)
                    nodesOnRHS.append(nodeWall)
                if segmentIdx and n2 == 0 and n1 > 0 and n1 < elementsCountAroundHaustrum:
                    nodeWall = []
                    for n3 in range(elementsCountThroughWall, -1, -1):
                        nodeWall.append(idx + 1 - elementsCountAround * n3)
                    nodesOnProximal.append(nodeWall)
                if n2 == elementsCountAlongSegment and n1 > 0 and n1 < elementsCountAroundHaustrum:
                    nodeWall = []
                    for n3 in range(elementsCountThroughWall, -1, -1):
                        nodeWall.append(idx + 1 - elementsCountAround * n3)
                    nodesOnDistal.append(nodeWall)

    sideNodes = []
    for n2 in range(elementsCountAlongSegment + 1):
        for n3 in range(elementsCountThroughWall + 1):
            if n2 == 0 and segmentIdx == 0:
                sideNodes.append((n2 + 1) * (n3 + 1))
            elif segmentIdx and n2 == 0:
                sideNodes.append(nodesOnLHS[n2][n3])
                for n in range(len(nodesOnProximal)):
                    sideNodes.append(nodesOnProximal[n][n3])
                sideNodes.append(nodesOnRHS[n2][n3])
            elif n2 == elementsCountAlongSegment:
                sideNodes.append(nodesOnLHS[n2 + (0 if segmentIdx else -1)][n3])
                for n in range(len(nodesOnDistal)):
                    sideNodes.append(nodesOnDistal[n][n3])
                sideNodes.append(nodesOnRHS[n2 + (0 if segmentIdx else -1)][n3])
            else:
                for n1 in range(2):
                    if n1 == 0:
                        sideNodes.append(nodesOnLHS[n2 + (0 if segmentIdx else -1)][n3])
                    else:
                        sideNodes.append(nodesOnRHS[n2 + (0 if segmentIdx else -1)][n3])

    trackSurfaceOstium = TrackSurface(elementsAroundTrackSurface, elementsAlongTrackSurface,
                                      xTrackSurface, d1TrackSurface, d2TrackSurface)

    # # Visualise track surface
    # nodeIdentifier, elementIdentifier = trackSurfaceOstium.generateMesh(region)

    # Find centre position
    # track along ileum path and since cxIleum[1] could be above or below the track surface, we check both side to
    # determine direction to track. At each point, find the nearest position and take the diff between nearest point
    # to the point in line, keep tracking till diff is close to zero.

    xTol = 1.0E-6
    arcStart = 0.0
    arcEnd = networkLayout.arcLengthOfGroupsAlong[1]
    nearestPosition = trackSurfaceOstium.findNearestPosition(cxIleum[0])
    xNearestStart = trackSurfaceOstium.evaluateCoordinates(nearestPosition, derivatives=False)
    distStart = magnitude([cxIleum[0][c] - xNearestStart[c] for c in range(3)])
    nearestPosition = trackSurfaceOstium.findNearestPosition(cxIleum[-1])
    xNearestEnd = trackSurfaceOstium.evaluateCoordinates(nearestPosition, derivatives=False)
    distEnd = magnitude([cxIleum[-1][c] - xNearestEnd[c] for c in range(3)])

    for iter in range(100):
        arcDistance = (arcStart + arcEnd) * 0.5
        x, d1 = interp.getCubicHermiteCurvesPointAtArcDistance(cxIleum, cd1Ileum, arcDistance)[0:2]
        nearestPosition = trackSurfaceOstium.findNearestPosition(x)
        xNearest = trackSurfaceOstium.evaluateCoordinates(nearestPosition, derivatives=False)
        dist = magnitude([x[c] - xNearest[c] for c in range(3)])

        if abs(distStart - distEnd) > xTol:
            if distStart < distEnd:
                arcEnd = arcDistance
                distEnd = dist
            else:
                arcStart = arcDistance
                distStart = dist

        else:
            xCentre, d1Centre, d2Centre = trackSurfaceOstium.evaluateCoordinates(nearestPosition, derivatives=True)
            normAxis = normalize([-d for d in d1])
            eIdx = interp.getNearestPointIndex(cxIleum, xCentre) - 1
            arcLenghtSum = 0.0
            for e in range(eIdx):
                arcLenghtSum += interp.getCubicHermiteArcLength(cxIleum[e], cd1Ileum[e],
                                                                cxIleum[e + 1], cd1Ileum[e + 1])
            xi = (arcDistance - arcLenghtSum)/\
                 interp.getCubicHermiteArcLength(cxIleum[eIdx], cd1Ileum[eIdx], cxIleum[eIdx + 1], cd1Ileum[eIdx + 1])
            d2Centre = interp.interpolateCubicHermite(cd2Ileum[eIdx], cd12Ileum[eIdx], cd2Ileum[eIdx + 1],
                                                      cd12Ileum[eIdx + 1], xi)
            break
    if iter > 98:
        print('Search for ileum entry centre - Max iters reached:', iter)

    ostiumSettings['Number of elements around ostium'] = elementsCountAlongSegment
    elementsCountAroundOstium = ostiumSettings['Number of elements around ostium']
    ostiumSettings['Use linear through ostium wall'] = options['Use linear through wall']
    ostiumSettings['Use linear through vessel wall'] = options['Use linear through wall']

    ileumGroup = AnnotationGroup(region, get_smallintestine_term("ileum"))
    ileumMeshGroup = ileumGroup.getMeshGroup(mesh)
    ileocecalJunctionGroup = AnnotationGroup(region, get_smallintestine_term("ileocecal junction"))
    ileocecalJunctionMeshGroup = ileocecalJunctionGroup.getMeshGroup(mesh)
    smallIntestineGroup = AnnotationGroup(region, get_smallintestine_term("small intestine"))
    smallIntestineMeshGroup = smallIntestineGroup.getMeshGroup(mesh)
    cecumMeshGroup = cecumGroup.getMeshGroup(mesh)
    allAnnotationGroups += [ileumGroup, ileocecalJunctionGroup, smallIntestineGroup]

    ostiumWallAnnotationGroups = []
    if elementsCountThroughWall == 4:
        ileumMucosaGroup = AnnotationGroup(region, get_smallintestine_term("mucosa of ileum"))
        ileumSubmucosaGroup = AnnotationGroup(region, get_smallintestine_term("submucosa of ileum"))
        ileumCircularGroup = AnnotationGroup(region, get_smallintestine_term("circular muscle layer of ileum"))
        ileumLongitudinalGroup = AnnotationGroup(region,
                                                 get_smallintestine_term("longitudinal muscle layer of ileum"))

        ostiumWallAnnotationGroups = [[ileumMucosaGroup, mucosaGroup],
                                      [ileumSubmucosaGroup, submucosaGroup],
                                      [ileumCircularGroup, circularMuscleGroup],
                                      [ileumLongitudinalGroup, longitudinalMuscleGroup]]

        allAnnotationGroups += [ileumMucosaGroup, ileumSubmucosaGroup,
                                ileumCircularGroup, ileumLongitudinalGroup]

    # Points from track surface and vessel end
    xPath = [cxIleum[0], xCentre]
    d1Path = [cd1Ileum[0], [-d for d in normAxis]]
    d2Path = [cd2Ileum[0], d2Centre]
    d3Path = [cd3Ileum[0], [-d for d in d1Centre]]
    d12Path = [cd2Ileum[0], [0.0, 0.0, 0.0]]
    d13Path = [cd3Ileum[0], [0.0, 0.0, 0.0]]

    networkLayoutIleum = CustomNetworkLayout(xPath, d1Path, d2Path, d3Path, d12Path, d13Path)

    nextNodeIdentifier, nextElementIdentifier, (o1_x, o1_d1, o1_d2, o1_d3, o1_NodeId, o1_Positions) = \
        generateOstiumMesh(region, ostiumSettings, trackSurfaceOstium, networkLayoutIleum,
                           startNodeIdentifier=nextNodeIdentifier, startElementIdentifier=nextElementIdentifier,
                           nodeIdProximal=nodeIdProximalIleum, xProximal=xProximalIleum, d1Proximal=d1ProximalIleum,
                           d2Proximal=d2ProximalIleum, d3Proximal=d3ProximalIleum,
                           vesselMeshGroups=[[cecumMeshGroup, smallIntestineMeshGroup, ileumMeshGroup]],
                           ostiumMeshGroups=[cecumMeshGroup, ileocecalJunctionMeshGroup],
                           wallAnnotationGroups=ostiumWallAnnotationGroups, coordinates=None)

    ostiumFaceStartNode = nextNodeIdentifier
    ostiumFaceStartElement = nextElementIdentifier

    # Create location of annulus
    xAnnulusOuter = [[] for x in range(elementsCountAroundOstium)]
    xAnnulusOuterPosition = [[] for x in range(elementsCountAroundOstium)]
    d1AnnulusNorm = []
    d1AnnulusOuter = []
    e1Left = elementsAroundTrackSurface
    e1Right = 0
    e2Top = 0
    e2Bottom = elementsAlongTrackSurface
    sf = magnitude(networkLayoutIleum.cd2Path[-1]) * 0.35
    for n1 in range(elementsCountAroundOstium):
        normD2 = normalize(o1_d2[-1][n1])
        d1AnnulusNorm.append(normD2)
        d1AnnulusOuter.append(set_magnitude(o1_d2[-1][n1], sf))
        x = [o1_x[-1][n1][c] + sf * normD2[c] for c in range(3)]
        nearestPosition = trackSurfaceOstium.findNearestPosition(x)
        e1 = nearestPosition.e1
        e2 = nearestPosition.e2
        if e1 < e1Left:
            e1Left = e1
        if e1 > e1Right:
            e1Right = e1
        if e2 > e2Top:
            e2Top = e2
        if e2 < e2Bottom:
            e2Bottom = e2
        xAnnulusOuterPosition[n1] = nearestPosition
        xAnnulusOuter[n1] = trackSurfaceOstium.evaluateCoordinates(nearestPosition)

    d2AnnulusOuter = []
    for n in range(elementsCountAlongSegment):
        d = findDerivativeBetweenPoints(xAnnulusOuter[n], xAnnulusOuter[(n + 1) % elementsCountAroundOstium])
        d2AnnulusOuter.append(d)
    d2AnnulusOuter = interp.smoothCubicHermiteDerivativesLoop(xAnnulusOuter, d2AnnulusOuter)
    d3Annulus = []
    for n in range(elementsCountAroundOstium):
        d3 = normalize(cross(normalize(d2AnnulusOuter[n]), d1AnnulusNorm[n]))
        d3Annulus.append(d3)
    annulusD2Curvature = interp.getCurvaturesAlongCurve(xAnnulusOuter, d2AnnulusOuter, d3Annulus, loop=True)

    # # Visualise annulus
    # for n1 in range(len(xAnnulusOuter)):
    #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAnnulusOuter[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1AnnulusOuter[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2AnnulusOuter[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    #     nextNodeIdentifier += 1

    # base row counts with 8 elements around ostium
    rowsIncrement = int((elementsCountAlongSegment - 8) / 4)
    rowsAbove = 1 + rowsIncrement
    rowsOstium = 4 + rowsIncrement * 2
    rowsBelow = 3 + rowsIncrement

    # # rowIdx along segment
    startRowIdx = rowsBelow + 1
    endRowIdx = rowsBelow + rowsOstium - 1

    # sample along the midline of ostium
    rotAngle = math.pi
    rotFrame = axis_angle_to_rotation_matrix(d3ApexUnit, rotAngle)
    d1A = [rotFrame[j][0] * d1Cecum[1][0] + rotFrame[j][1] * d1Cecum[1][1] +
           rotFrame[j][2] * d1Cecum[1][2] for j in range(3)]

    rotAngle = math.pi
    rotFrame = axis_angle_to_rotation_matrix(normalize(d3Annulus[0]), rotAngle)
    d2B = [rotFrame[j][0] * d1AnnulusOuter[0][0] + rotFrame[j][1] * d1AnnulusOuter[0][1] +
           rotFrame[j][2] * d1AnnulusOuter[0][2] for j in range(3)]

    # sample along from apex to annulus start
    if segmentIdx == 0:
        xStart = xCecum[1]
    else:
        idx = int(elementsAroundTrackSurface * 0.5)
        xStart = xTrackSurface[idx]

    xPositionA = trackSurfaceOstium.findNearestPosition(xStart)
    xProportionA = trackSurfaceOstium.getProportion(xPositionA)
    derivativeA = trackSurfaceOstium.evaluateCoordinates(xPositionA, derivatives=True)[2]
    derivativeA = set_magnitude(derivativeA, magnitude(d2ApexInner))
    derivativeMagnitudeA = magnitude(derivativeA)

    xPositionB = trackSurfaceOstium.findNearestPosition(xAnnulusOuter[0])
    xProportionB = trackSurfaceOstium.getProportion(xPositionB)

    nx, nd1, nd2, nd3, proportions = \
        trackSurfaceOstium.createHermiteCurvePoints(
            xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1],
            rowsBelow + 1, derivativeStart=derivativeA, derivativeEnd=None)

    pxAlongMidLine, pd2AlongMidLine, pd1AlongMidLine = \
        trackSurfaceOstium.resampleHermiteCurvePointsSmooth(
            nx, nd1, nd2, nd3, proportions, derivativeMagnitudeStart=derivativeMagnitudeA)[0:3]

    annulusD2ZeroMag = magnitude(pd2AlongMidLine[-1])

    for n in range(len(pd1AlongMidLine)):
        pd1AlongMidLine[n] = [-d for d in pd1AlongMidLine[n]]

    xPositionA = trackSurfaceOstium.findNearestPosition(xAnnulusOuter[int(elementsCountAlongSegment * 0.5)])
    xProportionA = trackSurfaceOstium.getProportion(xPositionA)

    xB = xTrackSurface[-elementsCountAroundHalfHaustrum]
    xPositionB = trackSurfaceOstium.findNearestPosition(xB)
    xProportionB = trackSurfaceOstium.getProportion(xPositionB)
    derivativeB = d2TrackSurface[-elementsCountAroundHalfHaustrum]
    derivativeMagnitudeB = magnitude(d2TrackSurface[-elementsCountAroundHalfHaustrum])

    nx, nd1, nd2, nd3, proportions = \
        trackSurfaceOstium.createHermiteCurvePoints(
            xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1],
            rowsAbove + 1, derivativeStart=None, derivativeEnd=derivativeB)

    pxAlongMidLineBottom, pd2AlongMidLineBottom, pd1AlongMidLineBottom = \
        trackSurfaceOstium.resampleHermiteCurvePointsSmooth(
            nx, nd1, nd2, nd3, proportions, derivativeMagnitudeStart=None,
            derivativeMagnitudeEnd=derivativeMagnitudeB)[0:3]

    annulusD2HalfOstiumMag = magnitude(pd2AlongMidLineBottom[0])

    for n in range(len(pd1AlongMidLineBottom)):
        pd1AlongMidLineBottom[n] = [-d for d in pd1AlongMidLineBottom[n]]

    # sample points around colon to ostium
    annulusIdx = 1
    xAroundAlong = []
    d1AroundAlong = []

    sideElements = int(elementsCountAroundHaustrum * 0.5)
    for n2 in range(elementsCountAlongSegment):
        xAround = []
        d1Around = []
        if n2 == 0 and segmentIdx == 0:
            for n1 in range(elementsCountAroundHaustrum + 1):
                xAround.append(xApex)
                d1Around.append(d1A)
        else:
            for n in range(2):
                # LHS
                if n == 0:
                    xPositionA = trackSurfaceOstium.findNearestPosition(
                        xTrackSurface[n2 * (elementsCountAroundHaustrum + 1)])
                    xProportionA = trackSurfaceOstium.getProportion(xPositionA)
                    derivativeA = d1TrackSurface[n2 * (elementsCountAroundHaustrum + 1)]
                    derivativeMagnitudeA = magnitude(derivativeA)

                    if n2 < rowsBelow + 1:
                        xPositionB = trackSurfaceOstium.findNearestPosition(pxAlongMidLine[n2])
                        derivativeB = None
                        derivativeMagnitudeB = None

                    elif n2 >= rowsBelow + rowsOstium:
                        idx = n2 - (rowsBelow + rowsOstium) + 1
                        xPositionB = trackSurfaceOstium.findNearestPosition(pxAlongMidLineBottom[idx])
                        derivativeB = None
                        derivativeMagnitudeB = None

                    else:
                        xPositionB = trackSurfaceOstium.findNearestPosition(xAnnulusOuter[annulusIdx])
                        derivativeB = None
                        derivativeMagnitudeB = None

                    xProportionB = trackSurfaceOstium.getProportion(xPositionB)

                else:  # RHS
                    if n2 < rowsBelow + 1:
                        xPositionA = trackSurfaceOstium.findNearestPosition(pxAlongMidLine[n2])
                        derivativeA = None
                        derivativeMagnitudeA = None

                    elif n2 >= rowsBelow + rowsOstium:
                        idx = n2 - (rowsBelow + rowsOstium) + 1
                        xPositionA = trackSurfaceOstium.findNearestPosition(pxAlongMidLineBottom[idx])
                        derivativeA = None
                        derivativeMagnitudeA = None

                    else:
                        xPositionA = trackSurfaceOstium.findNearestPosition(xAnnulusOuter[-annulusIdx])
                        derivativeA = None
                        derivativeMagnitudeA = None

                    xProportionA = trackSurfaceOstium.getProportion(xPositionA)
                    xPositionB = trackSurfaceOstium.findNearestPosition(
                        xTrackSurface[n2 * (elementsCountAroundHaustrum + 1) + elementsCountAroundHaustrum])
                    xProportionB = trackSurfaceOstium.getProportion(xPositionB)
                    derivativeB = d1TrackSurface[n2 * (elementsCountAroundHaustrum + 1) + elementsCountAroundHaustrum]
                    derivativeMagnitudeB = magnitude(derivativeB)

                nx, nd1, nd2, nd3, proportions = \
                    trackSurfaceOstium.createHermiteCurvePoints(
                        xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1],
                        sideElements +(0 if (n2 < rowsBelow + 1 or n2 >= rowsBelow + rowsOstium) else -1),
                        derivativeStart=derivativeA, derivativeEnd=derivativeB)

                nx, nd1 = \
                    trackSurfaceOstium.resampleHermiteCurvePointsSmooth(
                        nx, nd1, nd2, nd3, proportions, derivativeMagnitudeStart=derivativeMagnitudeA,
                        derivativeMagnitudeEnd=derivativeMagnitudeB)[0:2]

                if n == 0:
                    xAround += nx
                    d1Around += nd1
                else:
                    xAround += nx[(1 if (n2 < rowsBelow + 1 or n2 >= rowsBelow + rowsOstium) else 0):]
                    d1Around += nd1[(1 if (n2 < rowsBelow + 1 or n2 >= rowsBelow + rowsOstium) else 0):]
                if n2 < rowsBelow + 1 or n2 >= rowsBelow + rowsOstium:
                    d1Around = interp.smoothCubicHermiteDerivativesLine(xAround, d1Around, fixStartDerivative=True,
                                                                        fixEndDerivative=True)
            if n2 >= rowsBelow + 1 and n2 < rowsBelow + rowsOstium:
                annulusIdx += 1

        xAroundAlong.append(xAround)
        d1AroundAlong.append(d1Around)

    xAround = []
    d1Around = []
    for n1 in range(elementsCountAroundHaustrum + 1):
        xAround.append(xTrackSurface[n1 + (elementsCountAroundHaustrum + 1) * elementsCountAlongSegment])
        d1Around.append(d1TrackSurface[n1 + (elementsCountAroundHaustrum + 1) * elementsCountAlongSegment])
    xAroundAlong.append(xAround)
    d1AroundAlong.append(d1Around)
    d1AroundAlongOriginal = copy.deepcopy(d1AroundAlong)

    # Calculate d2 along segment
    d2AroundAlong = []
    xAlongAll = []
    d2AlongAll = []

    for n2 in range(len(xAroundAlong)):
        d2Around = []
        for n1 in range(len(xAroundAlong[n2])):
            d2Around.append([0.0, 0.0, 0.0])
        d2AroundAlong.append(d2Around)

    for n1 in range(elementsCountAroundHaustrum + 1):
        nxAlong = []
        nd2Along = []
        nxTop = []
        nd2Top = []
        annulusIdx = 1
        if n1 < elementsCountAroundHalfHaustrum - 2:
            for n2 in range(elementsCountAlongSegment):
                nxAlong.append(xAroundAlong[n2][n1])
                nd2Along.append(findDerivativeBetweenPoints(xAroundAlong[n2][n1], xAroundAlong[n2 + 1][n1]))
            nxAlong.append(xAroundAlong[n2 + 1][n1])
            nd2Along.append(d2TrackSurface[(elementsCountAroundHaustrum + 1) * elementsCountAlongSegment + n1])
            nd2Along = interp.smoothCubicHermiteDerivativesLine(nxAlong, nd2Along, fixStartDerivative=True,
                                                                fixEndDerivative=True)

            xAlongAll.append(nxAlong)
            d2AlongAll.append(nd2Along)

            for n2 in range(len(nd2Along)):
                d2AroundAlong[n2][n1] = nd2Along[n2]

        # Deal with annulus
        elif n1 == elementsCountAroundHalfHaustrum - 2:
            # Smooth from apex to annulus
            for n2 in range(startRowIdx + 1):
                nxAlong.append(xAroundAlong[n2][n1])
                nd2Along.append(findDerivativeBetweenPoints(xAroundAlong[n2][n1], xAroundAlong[n2 + 1][n1]))
            nd2Along = interp.smoothCubicHermiteDerivativesLine(nxAlong, nd2Along, fixStartDerivative=True)
            nd2AlongBottomLHS = interp.smoothCubicHermiteDerivativesLine(nxAlong, nd2Along, fixStartDerivative=True)

            # Make sure d2 at annulus is length of element below
            nd2Along[-1] = set_magnitude(d1AnnulusOuter[1], magnitude(nd2Along[-1]))

            # Make magnitude of annulus d2 between start and end row equivalent to arclength of element on its left
            for m in range(endRowIdx - startRowIdx - 1):
                nxAlong.append(xAnnulusOuter[2 + m])
                n2Idx = m + startRowIdx + 1
                nd2Along.append(set_magnitude(d1AnnulusOuter[2+m], magnitude(d1AroundAlong[n2Idx][n1])))

            # Smooth from annulus to end of cecum
            for n2 in range(endRowIdx, elementsCountAlongSegment):
                nxTop.append(xAroundAlong[n2][n1])
                nd2Top.append(findDerivativeBetweenPoints(xAroundAlong[n2][n1], xAroundAlong[n2 + 1][n1]))
            nxTop.append(xAroundAlong[elementsCountAlongSegment][n1])
            nd2Top.append(d2TrackSurface[(elementsCountAroundHaustrum + 1) * elementsCountAlongSegment + n1])
            nd2Top = interp.smoothCubicHermiteDerivativesLine(nxTop, nd2Top, fixEndDerivative=True)
            nd2AlongTopLHS = interp.smoothCubicHermiteDerivativesLine(nxTop, nd2Top, fixEndDerivative=True)

            # Make sure d2 at annulus is length of element above
            nd2Top[0] = set_magnitude(d1AnnulusOuter[int(elementsCountAroundOstium * 0.5) - 1],
                                            magnitude(nd2Top[0]))

            nxAlong += nxTop
            nd2Along += nd2Top

            xAlongAll.append(nxAlong)
            d2AlongAll.append(nd2Along)

            for n2 in range(len(nd2Along)):
                d2AroundAlong[n2][n1] = nd2Along[n2]
                if n1 == elementsCountAroundHalfHaustrum - 2 and startRowIdx -1 < n2 < endRowIdx + 1:
                    d1AroundAlong[n2][n1] = d2AnnulusOuter[annulusIdx]
                    annulusIdx += 1

        elif n1 == elementsCountAroundHalfHaustrum - 1:
            for n2 in range(len(pxAlongMidLine)):
                d2AroundAlong[n2][n1] = pd2AlongMidLine[n2]
            for n in range(1, len(pxAlongMidLineBottom)):
                d2AroundAlong[n + endRowIdx][n1] = pd2AlongMidLineBottom[n]
            nxAlong = pxAlongMidLine + pxAlongMidLineBottom
            xAlongAll.append(nxAlong)
            d2AlongAll.append(pd2AlongMidLine + pd2AlongMidLineBottom)

        elif n1 == elementsCountAroundHalfHaustrum:
            # Smooth from apex to annulus
            for n2 in range(startRowIdx + 1):
                 nxAlong.append(xAroundAlong[n2][n1 + (0 if n2 < startRowIdx else -1)])
                 nd2Along.append(findDerivativeBetweenPoints(xAroundAlong[n2][n1], xAroundAlong[n2 + 1][n1]))
            nd2Along = interp.smoothCubicHermiteDerivativesLine(nxAlong, nd2Along, fixStartDerivative=True)
            nd2AlongBottomRHS = interp.smoothCubicHermiteDerivativesLine(nxAlong, nd2Along, fixStartDerivative=True)

            # Make sure d2 at annulus is length of element below
            nd2Along[-1] = set_magnitude(d1AnnulusOuter[-1], magnitude(nd2Along[-1]))

            # Make magnitude of annulus d2 between start and end row equivalent to arclength of element on its right
            for m in range(endRowIdx - startRowIdx - 1):
                nxAlong.append(xAnnulusOuter[-(2 + m)])
                n2Idx = m + startRowIdx + 1
                nd2Along.append(set_magnitude(d1AnnulusOuter[-(2 + m)],
                                                    magnitude(d1AroundAlong[n2Idx][n1 - 1])))

            # Smooth from annulus to end of cecum
            for n2 in range(endRowIdx, elementsCountAlongSegment):
                nxTop.append(xAroundAlong[n2][n1 + (0 if n2 > endRowIdx else -1)])
                if n2 > endRowIdx - 1:
                    nxAlongNext = xAroundAlong[n2 + 1][n1]
                else:
                    nxAlongNext = xAroundAlong[n2 + 1][n1 - 1]
                nd2Top.append(findDerivativeBetweenPoints(xAroundAlong[n2][n1 + (0 if n2 > endRowIdx else -1)],
                                                          nxAlongNext))
            nxTop.append(xAroundAlong[elementsCountAlongSegment][n1])
            nd2Top.append(d2TrackSurface[(elementsCountAroundHaustrum + 1) * elementsCountAlongSegment + n1])
            nd2Top = interp.smoothCubicHermiteDerivativesLine(nxTop, nd2Top, fixEndDerivative=True)
            nd2AlongTopRHS = interp.smoothCubicHermiteDerivativesLine(nxTop, nd2Top, fixEndDerivative=True)

            # Make sure d2 at annulus is length of element above
            nd2Top[0] = set_magnitude(d1AnnulusOuter[int(elementsCountAroundOstium * 0.5) + 1],
                                            magnitude(nd2Top[0]))

            nxAlong += nxTop
            nd2Along += nd2Top

            xAlongAll.append(nxAlong)
            d2AlongAll.append(nd2Along)

            for n2 in range(len(nd2Along)):
                if n2 < startRowIdx or n2 > endRowIdx:
                    n1Idx = n1
                else:
                    n1Idx = n1 - 1
                    if n1 == elementsCountAroundHalfHaustrum:
                        d1AroundAlong[n2][n1Idx] = d2AnnulusOuter[-annulusIdx]
                        annulusIdx += 1
                d2AroundAlong[n2][n1Idx] = nd2Along[n2]

        else:
            for n2 in range(elementsCountAlongSegment):
                nxAlong.append(xAroundAlong[n2][n1 + (0 if (n2 < startRowIdx or n2 > endRowIdx) else -1)])
                if n2 < startRowIdx - 1 or n2 > endRowIdx -1:
                    nxAlongNext = xAroundAlong[n2 + 1][n1]
                else:
                    nxAlongNext = xAroundAlong[n2 + 1][n1 - 1]
                nd2Along.append(findDerivativeBetweenPoints(
                    xAroundAlong[n2][n1 + (0 if (n2 < startRowIdx or n2 > endRowIdx) else -1)], nxAlongNext))

            nxAlong.append(xAroundAlong[n2 + 1][n1])
            nd2Along.append(d2TrackSurface[(elementsCountAroundHaustrum + 1) * elementsCountAlongSegment + n1])
            nd2Along = interp.smoothCubicHermiteDerivativesLine(nxAlong, nd2Along, fixStartDerivative=True,
                                                                fixEndDirection=True)

            xAlongAll.append(nxAlong)
            d2AlongAll.append(nd2Along)

            for n2 in range(len(nd2Along)):
                if n2 < startRowIdx or n2 > endRowIdx:
                    n1Idx = n1
                else:
                    n1Idx = n1 - 1
                d2AroundAlong[n2][n1Idx] = nd2Along[n2]

    # Calculate d3
    d3UnitAroundAlong = []
    for n2 in range(len(xAroundAlong)):
        d3Around = []
        for n1 in range(len(xAroundAlong[n2])):
            d3Around.append(normalize(
                cross(normalize(d1AroundAlong[n2][n1]), normalize(d2AroundAlong[n2][n1]))))
        d3UnitAroundAlong.append(d3Around)

    # Calculate curvatures
    # Curvatures around
    d1Curvature = []
    if segmentIdx == 0:
        d1Curvature.append([0.0])

    for n2 in range(1 if segmentIdx == 0 else 0, elementsCountAlongSegment + 1):
        if n2 < startRowIdx or n2 > endRowIdx: #== elementsCountAlongSegment:
            d1Curvature.append(interp.getCurvaturesAlongCurve(xAroundAlong[n2], d1AroundAlongOriginal[n2],
                                                             d3UnitAroundAlong[n2]))
        else:
            d1CurvatureLeft = interp.getCurvaturesAlongCurve(xAroundAlong[n2][:int(0.5 * len(xAroundAlong[n2]))],
                                                     d1AroundAlongOriginal[n2][:int(0.5 * len(xAroundAlong[n2]))],
                                                     d3UnitAroundAlong[n2][:int(0.5 * len(xAroundAlong[n2]))])

            d1CurvatureRight = interp.getCurvaturesAlongCurve(xAroundAlong[n2][int(0.5 * len(xAroundAlong[n2])):],
                                                      d1AroundAlongOriginal[n2][int(0.5 * len(xAroundAlong[n2])):],
                                                      d3UnitAroundAlong[n2][int(0.5 * len(xAroundAlong[n2])):])
            d1Curvature.append(d1CurvatureLeft + d1CurvatureRight)

    # Curvatures along
    d2Curvature = []
    for n2 in range(len(xAroundAlong)):
        d2CurvatureAround = []
        for n1 in range(len(xAroundAlong[n2])):
            d2CurvatureAround.append(0.0)
        d2Curvature.append(d2CurvatureAround)

    for n1 in range(elementsCountAroundHaustrum + 1):
        xAlong = xAlongAll[n1]
        d2Along = d2AlongAll[n1]
        d3UnitAlong = []
        annulusIdx = 1

        if n1 < elementsCountAroundHalfHaustrum - 2:
            for n2 in range(elementsCountAlongSegment):
                d3UnitAlong.append(d3UnitAroundAlong[n2][n1])
            d3UnitAlong.append(d3UnitAroundAlong[n2 + 1][n1])
            d2CurvatureAlong = interp.getCurvaturesAlongCurve(xAlong, d2Along, d3UnitAlong)
            for n2 in range(len(d2CurvatureAlong)):
                d2Curvature[n2][n1] = d2CurvatureAlong[n2]

        elif n1 == elementsCountAroundHalfHaustrum - 2:
            # From apex to annulus
            for n2 in range(elementsCountAlongSegment):
                d3UnitAlong.append(d3UnitAroundAlong[n2][n1])
            d3UnitAlong.append(d3UnitAroundAlong[n2 + 1][n1])
            d2CurvatureAlong = interp.getCurvaturesAlongCurve(xAlong[:startRowIdx + 1], nd2AlongBottomLHS,
                                                                   d3UnitAlong[:startRowIdx + 1])

            # Curvature of nodes along LHS annulus
            for m in range(endRowIdx - startRowIdx - 1):
                d2CurvatureAlong.append(d1Curvature[m + startRowIdx + 1][n1])

            # From annulus to distal end
            d2CurvatureAlong += interp.getCurvaturesAlongCurve(xAlong[endRowIdx:], nd2AlongTopLHS,
                                                             d3UnitAlong[endRowIdx:])

            for n2 in range(len(d2CurvatureAlong)):
                d2Curvature[n2][n1] = d2CurvatureAlong[n2]

            for n2 in range(startRowIdx, endRowIdx + 1):
                d1Curvature[n2][n1] = annulusD2Curvature[annulusIdx]
                annulusIdx += 1

        elif n1 == elementsCountAroundHalfHaustrum - 1:
            # From apex to annulus
            for n2 in range(len(pxAlongMidLine)):
                d3UnitAlong.append(d3UnitAroundAlong[n2][n1])
            d2CurvatureAlong = interp.getCurvaturesAlongCurve(pxAlongMidLine, pd2AlongMidLine, d3UnitAlong)
            d2CurvatureAnnulusZero = d2CurvatureAlong[-1]
            for n2 in range(len(d2CurvatureAlong) - 1):
                d2Curvature[n2][n1] = d2CurvatureAlong[n2]

            d3UnitAlong = []
            for n in range(len(pxAlongMidLineBottom)):
                d3UnitAlong.append(d3UnitAroundAlong[n + endRowIdx][n1])
            d2CurvatureAlong = interp.getCurvaturesAlongCurve(pxAlongMidLineBottom, pd2AlongMidLineBottom, d3UnitAlong)
            d2CurvatureAlongHalfOstium = d2CurvatureAlong[0]
            for n in range(1, len(pd2AlongMidLineBottom)):
                nIdx = n + endRowIdx
                d2Curvature[nIdx][n1] = d2CurvatureAlong[n]

        elif n1 == elementsCountAroundHalfHaustrum:
            # From apex to annulus
            for n2 in range(elementsCountAlongSegment):
                d3UnitAlong.append(d3UnitAroundAlong[n2][n1 + (0 if (n2 < startRowIdx or n2 > endRowIdx) else -1)])
            d3UnitAlong.append(d3UnitAroundAlong[n2 + 1][n1])
            d2CurvatureAlong = interp.getCurvaturesAlongCurve(xAlong[:startRowIdx + 1], nd2AlongBottomRHS,
                                                             d3UnitAlong[:startRowIdx + 1])

            # Curvature of nodes along LHS annulus
            for m in range(endRowIdx - startRowIdx - 1):
                d2CurvatureAlong.append(d1Curvature[m + startRowIdx + 1][n1 - 1])

            # From annulus to distal end
            d2CurvatureAlong += interp.getCurvaturesAlongCurve(xAlong[endRowIdx:], nd2AlongTopRHS,
                                                              d3UnitAlong[endRowIdx:])

            for n2 in range(len(d2CurvatureAlong)):
                if n2 < startRowIdx or n2 > endRowIdx:
                    n1Idx = n1
                else:
                    n1Idx = n1 - 1
                d2Curvature[n2][n1Idx] = d2CurvatureAlong[n2]

            for n2 in range(startRowIdx, endRowIdx + 1):
                d1Curvature[n2][n1] = annulusD2Curvature[-annulusIdx]
                annulusIdx += 1

        else:
            for n2 in range(elementsCountAlongSegment):
                d3UnitAlong.append(d3UnitAroundAlong[n2][n1 + (0 if (n2 < startRowIdx or n2 > endRowIdx) else -1)])
            d3UnitAlong.append(d3UnitAroundAlong[n2 + 1][n1])
            d2CurvatureAlong = interp.getCurvaturesAlongCurve(xAlong, d2Along, d3UnitAlong)

            # Adjust for corners
            if n1 == elementsCountAroundHalfHaustrum:
                d2CurvatureAlong[endRowIdx] = annulusD2Curvature[int(elementsCountAlongSegment * 0.5) + 1]

            for n2 in range(len(d2CurvatureAlong)):
                if n2 < startRowIdx or n2 > endRowIdx:
                    n1Idx = n1
                else:
                    n1Idx = n1 - 1
                d2Curvature[n2][n1Idx] = d2CurvatureAlong[n2]

    # Adjust annulus points
    xAroundAlong[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), xAnnulusOuter[0])
    d1AroundAlong[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d2AnnulusOuter[0])
    d2AroundAlong[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), set_magnitude(d1AnnulusOuter[0],
                                                                                                  annulusD2ZeroMag))
    d3UnitAroundAlong[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d3Annulus[0])
    d1Curvature[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), annulusD2Curvature[0])
    d2Curvature[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d2CurvatureAnnulusZero)

    idx = int(elementsCountAlongSegment * 0.5)
    xAroundAlong[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), xAnnulusOuter[idx])
    d1AroundAlong[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d2AnnulusOuter[idx])
    d2AroundAlong[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), set_magnitude(d1AnnulusOuter[idx],
                                                                                                annulusD2HalfOstiumMag))
    d3UnitAroundAlong[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d3Annulus[idx])
    d1Curvature[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), annulusD2Curvature[idx])
    d2Curvature[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d2CurvatureAlongHalfOstium)

    # Create inner nodes
    sideNodesToDelete = []
    xList = []
    d1List = []
    d2List = []
    d3List = []
    nodeIdx = ostiumFaceStartNode
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

    count = 0
    for n2 in range(len(xAroundAlong)):
        idxThroughWall = []
        for n3 in range(elementsCountThroughWall + 1):
            xi3 = xi3List[n3] if elementsCountThroughWall > 1 else 1.0 / elementsCountThroughWall * n3
            idxAround = []

            for n1 in range(1 if (n2 == 0 and segmentIdx == 0) else len(xAroundAlong[n2])):
                if n1 == 0 or n1 == len(xAroundAlong[n2]) - 1:
                    sideNodesToDelete.append(nodeIdx)
                    newNodeIdx = sideNodes[count] + startNode - 1
                    count += 1

                elif n2 == 0 and segmentIdx and n1 > 0 and n1 < len(xAroundAlong[n2]) - 1:
                    sideNodesToDelete.append(nodeIdx)
                    newNodeIdx = sideNodes[count] + startNode - 1
                    count += 1

                elif n2 == len(xAroundAlong) - 1 and n1 > 0 and n1 < len(xAroundAlong[n2]) - 1:
                    sideNodesToDelete.append(nodeIdx)
                    newNodeIdx = sideNodes[count] + startNode - 1
                    count += 1

                else:
                    newNodeIdx = nodeIdx

                # Coordinates
                norm = d3UnitAroundAlong[n2][n1]
                xOut = xAroundAlong[n2][n1]
                xIn = [xOut[i] - norm[i] * wallThickness for i in range(3)]
                dWall = [wallThickness * c for c in norm]
                x = interp.interpolateCubicHermite(xIn, dWall, xOut, dWall, xi3)
                xList.append(x)

                # d1
                factor = 1.0 + wallThickness * (1.0 - xi3) * d1Curvature[n2][n1]
                d1 = [factor * c for c in d1AroundAlong[n2][n1]]
                d1List.append(d1)

                # d2
                factor = 1.0 + wallThickness * (1.0 - xi3) * d2Curvature[n2][n1]
                d2 = [factor * c for c in d2AroundAlong[n2][n1]]
                d2List.append(d2)

                # d3
                d3 = [c * wallThickness * (thicknessProportions[n3 + 1] if elementsCountThroughWall > 1 else 1.0)
                      for c in norm]
                d3List.append(d3)

                idxAround.append(newNodeIdx)
                nodeIdx += 1
            idxThroughWall.append(idxAround)
        idxMat.append(idxThroughWall)

    for n in range(len(xList)):
        if nextNodeIdentifier in sideNodesToDelete:
            nextNodeIdentifier += 1
            continue
        else:
            node = nodes.createNode(nextNodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xList[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1List[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2List[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3List[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nextNodeIdentifier += 1

    # Create elements
    elementIdxMat = []

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

    radiansPerElementAround = math.pi * 2.0 / elementsCountAround
    elementIdentifier = ostiumFaceStartElement

    for e2 in range(elementsCountAlongSegment):
        elementIdxThroughWall = []
        if segmentIdx == 0 and e2 == 0:  # pole
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(elementsCountAroundHaustrum):
                    va = e1 + startIdxElementsAround
                    vb = (e1 + startIdxElementsAround + 1)
                    eft1 = eftfactory.createEftShellPoleBottom(va * 100, vb * 100)
                    elementtemplateX.defineField(coordinates, -1, eft1)
                    element = mesh.createElement(elementIdentifier, elementtemplateX)
                    bni1 = e3 + 1 + startNode - 1
                    bni2 = idxMat[e2 + 1][e3][e1]
                    bni3 = idxMat[e2 + 1][e3][e1 + 1]
                    bni4 = bni1 + 1
                    bni5 = idxMat[e2 + 1][e3 + 1][e1]
                    bni6 = idxMat[e2 + 1][e3 + 1][e1 + 1]
                    nodeIdentifiers = [bni1, bni2, bni3, bni4, bni5, bni6]
                    element.setNodesByIdentifier(eft1, nodeIdentifiers)
                    # set general linear map coefficients
                    radiansAround = (e1 + 1) * radiansPerElementAround + sectorStartAngle
                    radiansAroundNext = (e1 + 2) * radiansPerElementAround + sectorStartAngle
                    scalefactors = [
                        -1.0,
                        math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                        math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround,
                        math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
                        math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround]
                    element.setScaleFactors(eft1, scalefactors)
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
            for e3 in range(elementsCountThroughWall):
                elementIdxAround = []
                for e1 in range(len(xAroundAlong[e2]) - 1):
                    offset = 0
                    if endRowIdx - startRowIdx > 1:
                        if e2 == startRowIdx and e1 > int(0.5 * len(xAroundAlong[e2])):
                            offset = -1
                        elif e2 == endRowIdx - 1 and e1 >= int(0.5 * len(xAroundAlong[e2])):
                            offset = 1

                    if (startRowIdx <= e2 <= endRowIdx - 1) and 0.5 * len(xAroundAlong[e2]) - 2 < e1 < 0.5 * len(
                            xAroundAlong[e2]):
                        continue
                    else:
                        eft1 = eftStandard
                        scaleFactors = []
                        elementtemplate1 = elementtemplateStandard
                        bni111 = idxMat[e2][e3][e1]
                        bni211 = idxMat[e2][e3][e1 + 1]
                        bni121 = idxMat[e2 + 1][e3][e1 + offset]
                        bni221 = idxMat[e2 + 1][e3][e1 + 1 + offset]
                        bni112 = idxMat[e2][e3 + 1][e1]
                        bni212 = idxMat[e2][e3 + 1][e1 + 1]
                        bni122 = idxMat[e2 + 1][e3 + 1][e1 + offset]
                        bni222 = idxMat[e2 + 1][e3 + 1][e1 + 1 + offset]
                        nodeIdentifiers = [bni111, bni211, bni121, bni221,
                                           bni112, bni212, bni122, bni222]

                        if e2 == startRowIdx - 1:
                            if e1 == elementsCountAroundHalfHaustrum - 3:
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX
                                # print(elementIdentifier) # 331

                            elif e1 == elementsCountAroundHalfHaustrum - 2:  # LHS bottom
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
                                # print(elementIdentifier) # 332
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 == elementsCountAroundHalfHaustrum - 1:  # RHS bottom
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
                                # print(elementIdentifier)  # 333
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX

                            elif e1 == elementsCountAroundHalfHaustrum:
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                       [(Node.VALUE_LABEL_D_DS2, [1])])
                                remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX
                                # print(elementIdentifier) #334

                        elif e2 == startRowIdx:
                            if e1 == elementsCountAroundHalfHaustrum - 3:
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
                                # print(elementIdentifier) # 339

                            elif e1 == elementsCountAroundHalfHaustrum:
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
                                # print(elementIdentifier) # 340

                        elif int(elementsCountAroundOstium*0.25) - 2 > 0 and \
                                startRowIdx + 1 <= e2 < startRowIdx + 1 + \
                                2.0 * (int(elementsCountAroundOstium*0.25) - 2):
                            if e1 == elementsCountAroundHalfHaustrum - 3:
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

                            elif e1 == elementsCountAroundHalfHaustrum - 1:
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

                        elif e2 == endRowIdx - 1:
                            if e1 == elementsCountAroundHalfHaustrum - 3:
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
                                # print(elementIdentifier) # 345

                            elif e1 == elementsCountAroundHalfHaustrum - 1:
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
                                # print(elementIdentifier) # 346

                        elif e2 == endRowIdx:
                            if e1 == elementsCountAroundHalfHaustrum - 3:  # end LHS -1
                                scaleFactors = [-1.0]
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                setEftScaleFactorIds(eft1, [1], [])
                                remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [1])])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX
                                # print('1', elementIdentifier) #351

                            elif e1 == elementsCountAroundHalfHaustrum:  # end RHS + 1
                                eft1 = eftfactory.createEftNoCrossDerivatives()
                                remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1,
                                                       [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                                elementtemplateX.defineField(coordinates, -1, eft1)
                                elementtemplate1 = elementtemplateX
                                # print(elementIdentifier) #354

                        element = mesh.createElement(elementIdentifier, elementtemplate1)
                        element.setNodesByIdentifier(eft1, nodeIdentifiers)
                        if scaleFactors:
                            element.setScaleFactors(eft1, scaleFactors)
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
    endPoints_x = [[None] * elementsCountAroundOstium for n3 in range(elementsCountThroughWall + 1)]
    endPoints_d1 = [[None] * elementsCountAroundOstium for n3 in range(elementsCountThroughWall + 1)]
    endPoints_d2 = [[None] * elementsCountAroundOstium for n3 in range(elementsCountThroughWall + 1)]
    endNode_Id = [[None] * elementsCountAroundOstium for n3 in range(elementsCountThroughWall + 1)]
    endDerivativesMap = [[None] * elementsCountAroundOstium for n3 in range(elementsCountThroughWall + 1)]
    endProportions = []

    for n3 in range(elementsCountThroughWall + 1):
        n1 = 0
        for nAround in range(elementsCountAroundOstium):
            if nAround == 0:
                idx = idxMat[startRowIdx][n3][elementsCountAroundHalfHaustrum - 1]
            elif 0 < nAround < (elementsCountAroundOstium * 0.5):
                idx = idxMat[startRowIdx + n1][n3][elementsCountAroundHalfHaustrum - 2]
                n1 += 1
            elif nAround == int(elementsCountAroundOstium * 0.5):
                n1 -= 1
                idx = idxMat[startRowIdx + n1][n3][elementsCountAroundHalfHaustrum - 1]
            else:
                idx = idxMat[startRowIdx + n1][n3][
                    elementsCountAroundHalfHaustrum - 1 + (
                        1 if (n1 == int(elementsCountAroundOstium * 0.5) - 2 or n1 == 0) else 0)]
                n1 -= 1

            endPoints_x[n3][nAround] = xList[idx - ostiumFaceStartNode]
            endPoints_d1[n3][nAround] = d1List[idx - ostiumFaceStartNode]
            endPoints_d2[n3][nAround] = d2List[idx - ostiumFaceStartNode]
            endNode_Id[n3][nAround] = idx

            if n3 == elementsCountThroughWall:  # outer layer
                endPosition = trackSurfaceOstium.findNearestPosition(endPoints_x[n3][nAround])
                endProportions.append(trackSurfaceOstium.getProportion(endPosition))

    for n3 in range(elementsCountThroughWall + 1):
        for nAround in range(elementsCountAroundOstium):
            endDerivativesMap[n3][nAround] = (None, None, None)

    startProportions = []
    for n in range(elementsCountAroundOstium):
        startProportions.append(trackSurfaceOstium.getProportion(o1_Positions[n]))

    cecumWallAnnotationGroups = []
    if elementsCountThroughWall == 4:
        cecumWallAnnotationGroups = [[mucosaGroup], [submucosaGroup], [circularMuscleGroup],
                                     [longitudinalMuscleGroup]]

    nextNodeIdentifier, nextElementIdentifier = createAnnulusMesh3d(
        nodes, mesh, nextNodeIdentifier, elementIdentifier,
        o1_x, o1_d1, o1_d2, None, o1_NodeId, None,
        endPoints_x, endPoints_d1, endPoints_d2, None, endNode_Id, endDerivativesMap,
        elementsCountRadial=1, meshGroups=[cecumMeshGroup],
        wallAnnotationGroups=cecumWallAnnotationGroups,
        tracksurface=trackSurfaceOstium,
        startProportions=startProportions, endProportions=endProportions,
        rescaleStartDerivatives=True, rescaleEndDerivatives=True, sampleBlend=0.0, fixMinimumStart=True,
        coordinates=coordinates)

    # Delete elements in new haustrum
    mesh_destroy_elements_and_nodes_by_identifiers(mesh, deleteElementIdentifier)

    return allAnnotationGroups, nextNodeIdentifier, nextElementIdentifier, nodesIdDistal, xDistal, d1Distal, \
           d2Distal, d3Distal


class CecumNetworkLayout:
    """
    Generates sampled network layout for cecum scaffold.
    """
    def __init__(self, region, networkLayout, termsAlong=[None], ileumSegmentIdx=0, cecumSegmentIdx=[1,2]):
        """
        :param region: Zinc region to define model in.
        :param networkLayout: Network layout subscaffold from meshtype_1d_networklayout1
        :param termsAlong: Annotation terms along length of network layout
        :param ileumSegmentIdx: Segment index of ileum branch.
        :param cecumSegmentIdx: Segment index of body of cecum.
        """

        # Extract length of each group along cecum from network layout
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
        cecumNodes = []
        lowerCecumNodes = []
        lowerCecumVersions = []
        ileumNodes = []
        ileumVersions = []

        for termName in termsAlong:
            tmpGroup = tmpFieldmodule.findFieldByName(termName).castGroup() if termName else None
            tmpNodeset = tmpGroup.getNodesetGroup(tmpNodes) if tmpGroup else tmpNodes

            if termName == "caecum":
                nodeiterator = tmpNodeset.createNodeiterator()
                node = nodeiterator.next()
                while node.isValid():
                    cecumNodes.append(node.getIdentifier())
                    node = nodeiterator.next()

                for i in range(len(networkSegments[cecumSegmentIdx[1]].getNodeIdentifiers())):
                    if networkSegments[cecumSegmentIdx[1]].getNodeIdentifiers()[i] in cecumNodes:
                        lowerCecumNodes.append(networkSegments[cecumSegmentIdx[1]].getNodeIdentifiers()[i])
                        lowerCecumVersions.append(networkSegments[cecumSegmentIdx[1]].getNodeVersions()[i])

                for i in range(2):
                    cx, cd1, cd2, cd3, cd12, cd13 = get_nodeset_path_ordered_field_parameters(
                        tmpNodeset, tmpCoordinates,
                        [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                         Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                         Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                        (lowerCecumNodes if i else networkSegments[cecumSegmentIdx[0]].getNodeIdentifiers()),
                        (lowerCecumVersions if i else networkSegments[cecumSegmentIdx[0]].getNodeVersions()))

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

            elif termName == "ileum part of cecum":
                for i in range(len(networkSegments[ileumSegmentIdx].getNodeIdentifiers())):
                    if networkSegments[ileumSegmentIdx].getNodeIdentifiers()[i] in cecumNodes:
                        ileumNodes.append(networkSegments[ileumSegmentIdx].getNodeIdentifiers()[i])
                        ileumVersions.append(networkSegments[ileumSegmentIdx].getNodeVersions()[i])

                cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = \
                    get_nodeset_path_ordered_field_parameters(tmpNodes, tmpCoordinates,
                                                              [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                                                               Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                                                               Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                                                              ileumNodes, ileumVersions)

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
        del tmpRegion

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
