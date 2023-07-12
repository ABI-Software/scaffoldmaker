"""
Generates a 3-D cecum mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""

import copy
import math

from cmlibs.maths.vectorops import add, sub
from cmlibs.utils.zinc.field import findOrCreateFieldCoordinates # KM
from cmlibs.zinc.element import Element, Elementbasis
from cmlibs.zinc.field import Field
from cmlibs.zinc.node import Node
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup, mergeAnnotationGroups
from scaffoldmaker.annotation.cecum_terms import get_cecum_term
from scaffoldmaker.annotation.smallintestine_terms import get_smallintestine_term
from scaffoldmaker.meshtypes.meshtype_1d_network_layout1 import MeshType_1d_network_layout1
from scaffoldmaker.meshtypes.meshtype_3d_colonsegment1 import ColonSegmentTubeMeshInnerPoints, \
    getFullProfileFromHalfHaustrum, getTeniaColi, createNodesAndElementsTeniaColi
from scaffoldmaker.meshtypes.meshtype_3d_ostium2 import MeshType_3d_ostium2, generateOstiumMesh
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils.geometry import createEllipsePoints # KM
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import setEftScaleFactorIds, remapEftNodeValueLabel, remapEftNodeValueLabelsVersion
from scaffoldmaker.utils.networkmesh import NetworkMesh, NetworkNode
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition
from scaffoldmaker.utils.zinc_utils import exnode_string_from_nodeset_field_parameters, \
    mesh_destroy_elements_and_nodes_by_identifiers, get_nodeset_path_field_parameters, \
    get_nodeset_path_ordered_field_parameters


class MeshType_3d_cecum1(Scaffold_base):
    '''
    Generates a 3-D cecum mesh with variable numbers
    of elements around, along the central line, and through wall.
    The cecum is created by a function that generates a cecum
    segment and uses tubemesh to map the segment along a central
    line profile. The proximal end of the cecum is closed up with
    an apex plate. An ostium is included to generate the
    ileo-cecal junction.
    '''

    parameterSetStructureStrings = {
        'Human 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-3-5"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                (1, [[7.50,-20.00,0.00], [0.00,8.28,0.00], [-1.50,0.00,0.00], [0.00,0.00,0.00], [0.00,-0.00,1.50], [0.00,0.00,0.00]]),
                (2, [[7.50,-10.86,0.00], [0.00,10.00,0.00], [-4.50,0.00,0.00], [0.00,0.00,0.00], [0.00,-0.00,4.50], [0.00,0.00,0.00]]),
                (3, [[7.50,0.00,0.00], [[7.50,0.00,0.00],[0.00,11.72,0.00]], [[0.00,10.00,0.00],[-4.50,0.00,0.00]], [[1.02,6.79,0.00],[0.00,0.00,0.00]], [[0.00,0.00,10.00],[0.00,-0.00,4.50]], [[0.00,0.00,5.77],[0.00,0.00,0.00]]]),
                (4, [[0.00,0.00,0.00], [7.50,0.00,0.00], [0.00,10.00,0.00], [0.00,0.00,0.00], [0.00,0.00,10.00], [0.00,0.00,0.00]]),
                (5, [[15.00,0.00,0.00], [7.50,0.00,0.00], [0.00,10.00,0.00], [0.00,0.00,0.00], [0.00,0.00,10.00], [0.00,0.00,0.00]])]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-4',
                    'name': get_cecum_term('caecum')[0],
                    'ontId': get_cecum_term('caecum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_smallintestine_term('ileum')[0],
                    'ontId': get_smallintestine_term('ileum')[1]
                }]
        }),
        'Human 2': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-3-5"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                (1, [[9.21,-19.56,6.43], [-2.02,7.27,-6.55], [-1.46,-0.27,0.15], [0.00,0.00,0.00], [-0.06,0.99,1.12], [0.00,0.00,0.00]]),
                (2, [[7.72,-10.79,1.37], [-0.90,10.04,-3.37], [-4.48,-0.36,0.11], [0.00,0.00,0.00], [-0.01,1.43,4.27], [0.00,0.00,0.00]]),
                (3, [[7.50,0.00,0.00], [[7.50,0.00,0.00],[0.45,11.25,0.61]], [[0.00,10.00,0.00],[-4.50,0.18,0.01]], [[1.02,6.79,0.00],[0.00,0.00,0.00]], [[0.00,0.00,10.00],[0.00,-0.24,4.49]], [[0.00,0.00,5.77],[0.00,0.00,0.00]]]),
                (4, [[0.00,0.00,0.00], [7.50,0.00,0.00], [0.00,10.00,0.00], [0.00,0.00,0.00], [0.00,0.00,10.00], [0.00,0.00,0.00]]),
                (5, [[15.00,0.00,0.00], [7.50,0.00,0.00], [0.00,10.00,0.00], [0.00,0.00,0.00], [0.00,0.00,10.00], [0.00,0.00,0.00]])]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-4',
                    'name': get_cecum_term('caecum')[0],
                    'ontId': get_cecum_term('caecum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_smallintestine_term('ileum')[0],
                    'ontId': get_smallintestine_term('ileum')[1]
                }]
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_network_layout1, {
            'scaffoldSettings': {
                "Structure": "1-2-3.2, 4-5-6-3-7"
            },
            'meshEdits': exnode_string_from_nodeset_field_parameters(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D_DS3, Node.VALUE_LABEL_D2_DS1DS3], [
                (1, [[167.00,29.75,51.53], [0.00,-9.63,-16.67], [-4.00,-0.00,0.00], [0.00,0.00,0.00], [0.00,6.64,-3.83], [0.00,0.00,0.00]]),
                (2, [[167.00,17.50,30.31], [0.00,-14.87,-25.77], [-15.00,-0.00,0.00], [0.00,0.00,0.00], [0.00,10.49,-6.05], [0.00,0.00,0.00]]),
                (3, [[167.00,0.00,0.00], [[30.00,0.00,0.00],[0.00,-20.13,-34.85]], [[0.00,35.00,0.00],[-15.00,0.00,0.00]], [[0.00,0.00,0.00],[0.00,0.00,0.00]], [[0.00,0.00,35.00],[0.00,8.49,-4.90]], [[0.00,0.00,0.00],[0.00,0.00,0.00]]]),
                (4, [[0.00,0.00,0.00], [60.00,0.00,0.00], [0.00,35.00,0.00], [0.00,0.00,0.00], [0.00,0.00,35.00], [0.00,0.00,0.00]]),
                (5, [[60.00,0.00,0.00], [60.00,0.00,0.00], [0.00,35.00,0.00], [0.00,0.00,0.00], [0.00,0.00,35.00], [0.00,0.00,0.00]]),
                (6, [[120.00,0.00,0.00], [53.50,0.00,0.00], [0.00,35.00,0.00], [0.00,0.00,0.00], [0.00,0.00,35.00], [0.00,0.00,0.00]]),
                (7, [[180.00,0.00,0.00], [4.00,-0.00,-0.00], [0.00,35.00,0.00], [0.00,0.00,0.00], [0.00,0.00,35.00], [0.00,0.00,0.00]])]),

            'userAnnotationGroups': [
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '3-6',
                    'name': get_cecum_term('caecum')[0],
                    'ontId': get_cecum_term('caecum')[1]
                },
                {
                    '_AnnotationGroup': True,
                    'dimension': 1,
                    'identifierRanges': '1-2',
                    'name': get_smallintestine_term('ileum')[0],
                    'ontId': get_smallintestine_term('ileum')[1]
                }]
        }),
    }
   
    ostiumDefaultScaffoldPackages = {
        'Human 1': ScaffoldPackage(MeshType_3d_ostium2, {
            'scaffoldSettings': {
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium wall thickness': 1.6,
                'Ostium wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Use linear through ostium wall': True,
                'Vessel wall thickness': 0.45,
                'Vessel wall relative thicknesses': [0.55, 0.15, 0.25, 0.05],
                'Use linear through vessel wall': True,
                'Use cross derivatives': False,
                'Refine': False,
                'Refine number of elements around': 4,
                'Refine number of elements along': 4,
                'Refine number of elements through wall': 1
            },
        }),
        'Pig 1': ScaffoldPackage(MeshType_3d_ostium2, {
            'scaffoldSettings': {
                'Number of elements around ostium': 8,
                'Number of elements along': 2,
                'Number of elements through wall': 1,  # not implemented for > 1
                'Unit scale': 1.0,
                'Outlet': False,
                'Ostium wall thickness': 2.0,
                'Use linear through ostium wall': True,
                'Vessel wall thickness': 2.0,
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
        if 'Human 2' in parameterSetName:
            centralPathOption = cls.parameterSetStructureStrings['Human 2']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Human 1']
        elif 'Pig 1' in parameterSetName:
            centralPathOption = cls.parameterSetStructureStrings['Pig 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Pig 1']
        else:
            centralPathOption = cls.parameterSetStructureStrings['Human 1']
            ostiumOption = cls.ostiumDefaultScaffoldPackages['Human 1']

        options = {
            'Central path': copy.deepcopy(centralPathOption),
            'Number of segments': 1,
            'Number of elements around tenia coli': 2,
            'Number of elements around haustrum': 8,
            'Number of elements along segment': 8,
            'Number of elements through wall': 1, # 4 later!
            'Corner inner radius factor': 0.5,
            'Haustrum inner radius factor': 0.4,
            'Segment length end derivative factor': 0.5,
            'Segment length mid derivative factor': 1.0, # 3.0,
            'Number of tenia coli': 3,
            'Start tenia coli width': 2.0, #10.0,
            'Start tenia coli width derivative': 2.0, #0.0,
            'End tenia coli width': 4.0, #10.0,
            'End tenia coli width derivative': 2.0, #0.0,
            'Tenia coli thickness': 0.5, #1.6,
            'Wall thickness': 1.6,
            'Mucosa relative thickness': 0.55,
            'Submucosa relative thickness': 0.15,
            'Circular muscle layer relative thickness': 0.25,
            'Longitudinal muscle layer relative thickness': 0.05,
            'Ileocecal junction': copy.deepcopy(ostiumOption),
            'Use cross derivatives': False,
            'Use linear through wall': True,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }

        if 'Pig 1' in parameterSetName:
            options['Number of segments'] = 5
            options['Haustrum inner radius factor'] = 0.25
            options['Segment length end derivative factor'] = 1.0
            options['Segment length mid derivative factor'] = 4.0
            options['Start tenia coli width'] = 5.0
            options['Start tenia coli width derivative'] = 0.0
            options['End tenia coli width'] = 5.0
            options['End tenia coli width derivative'] = 0.0
            options['Wall thickness'] = 2.0

        cls.updateSubScaffoldOptions(options)
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of segments',
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Number of elements along segment',
            'Number of elements through wall',
            'Corner inner radius factor',
            'Haustrum inner radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Number of tenia coli',
            'Start tenia coli width',
            'Start tenia coli width derivative',
            'End tenia coli width',
            'End tenia coli width derivative',
            'Tenia coli thickness',
            'Wall thickness',
            'Mucosa relative thickness',
            'Submucosa relative thickness',
            'Circular muscle layer relative thickness',
            'Longitudinal muscle layer relative thickness',
            'Ileocecal junction',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [ MeshType_1d_network_layout1 ]
        if optionName == 'Ileocecal junction':
            return [ MeshType_3d_ostium2 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.parameterSetStructureStrings.keys())
        if optionName == 'Ileocecal junction':
            return list(cls.ostiumDefaultScaffoldPackages.keys())
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
                parameterSetName = list(cls.parameterSetStructureStrings.keys())[0]
            return copy.deepcopy(cls.parameterSetStructureStrings[parameterSetName])
        if optionName == 'Ileocecal junction':
            if not parameterSetName:
                parameterSetName = list(cls.ostiumDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.ostiumDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_network_layout1)
        if not options['Ileocecal junction'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Ileocecal junction'):
            options['Ileocecal junction'] = cls.getOptionScaffoldPackage('Ileocecal junction', MeshType_3d_ostium2)
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
            'Number of elements along segment',
            'Corner inner radius factor',
            'Haustrum inner radius factor',
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
            if options['Number of elements through wall'] != (1 or 4):
                options['Number of elements through wall'] = 4
        cls.updateSubScaffoldOptions(options)

    @classmethod
    def updateSubScaffoldOptions(cls, options):
        '''
        Update ostium sub-scaffold options which depend on parent options.
        '''
        wallThickness = options['Wall thickness']
        ostiumOptions = options['Ileocecal junction']
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
        Generate the base tricubic Hermite mesh.
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """

        nextNodeIdentifier = 1
        nextElementIdentifier = 1
        cls.updateSubScaffoldOptions(options)
        geometricCentralPath = options['Central path']
        cecumTermsAlong = ['caecum', 'ileum']
        geometricCentralPath = CecumCentralPath(region, geometricCentralPath, cecumTermsAlong)
        annotationGroups, nextNodeIdentifier, nextElementIdentifier = \
            createCecumMesh3d(region, options, geometricCentralPath, nextNodeIdentifier,
                              nextElementIdentifier)

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

def getApexSegmentForCecum(xInner, d1Inner, d2Inner, elementsCountAroundHalfHaustrum,
                           elementsCountAroundTC, elementsCountAround, elementsCountAlongSegment, tcCount):
    """
    Generates the inner coordinates and derivatives for a cecum segment on the closed end.
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
    xFirstSegment += xInner[elementsCountAround * int(elementsCountAlongSegment * 0.5):] # second half of first regular segment
    d1FirstDirectionVector = vector.normalise(d1Inner[elementsCountAround]) # Store direction vector of first d1 intra-haustral for later
    d2Vector = xInner[elementsCountAround * int(elementsCountAlongSegment * 0.5):
                      elementsCountAround * (int(elementsCountAlongSegment * 0.5) + 1)] # half face of segment - apex
    d2FirstSegment = []
    for c in range(elementsCountAround):
        d2 = [d2Vector[c][0], d2Vector[c][1], 0.0 ]  # project onto x-y plane to get d2 pointing vertically
        d2FirstSegment.append(d2)
    d2FirstSegment += d2Inner[elementsCountAround * int(elementsCountAlongSegment*0.5):]

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
                d1 = matrix.rotateAboutZAxis(d2, math.pi*0.5)
                d1FirstSegmentSampled.append(d1)

        if n2 > 0:
            d1Around = []
            for n1 in range(elementsCountAroundHalfHaustrum):
                v1 = xAround[n1]
                v2 = xAround[n1 + 1]
                d1 = d1FirstDirectionVector if n1 == 0 else [v2[c] - v1[c] for c in range(3)]
                d2 = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c*arcLengthAround for c in vector.normalise(d1)]
                d1Around.append(dx_ds1)
            # Account for d1 of node sitting on half haustrum
            d1 = vector.normalise(
                [xAround[elementsCountAroundHalfHaustrum][c] - xAround[elementsCountAroundHalfHaustrum - 1][c]
                 for c in range(3)])
            dx_ds1 = [c * arcLengthAround for c in d1]
            d1Around.append(dx_ds1)

            d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, d1Around, fixStartDerivative=True)
            d1TCEdge = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC * 0.5)],
                                           vector.magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 - 1)]))
            d1Transition = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 1)],
                                               vector.magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 2)]))
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
    :param d1HaustrumHalfSet:
    :param tcCount:
    :return:
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

    xInnerRaw = []
    d2InnerRaw = []
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
        xInnerRaw.append(xSampled)
        d2InnerRaw.append(d2Sampled)

    # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
    for n2 in range(elementsCountAlong + 1):
        xAround = []
        d2Around = []

        for n1 in range(elementsCountAroundHalfHaustrum + 1):
            x = xInnerRaw[n1][n2]
            d2 = d2InnerRaw[n1][n2]
            xAround.append(x)
            d2Around.append(d2)

        d1InnerAroundList = []
        if n2 == 0:
            d1Corrected = d1ToSample[:elementsCountAroundHalfHaustrum + 1]

        else:
            for n1 in range(elementsCountAroundHalfHaustrum):
                v1 = xAround[n1]
                v2 = xAround[n1 + 1]
                d1 = d1FirstDirectionVector if n1 == 0 else [v2[c] - v1[c] for c in range(3)]
                d2 = [v2[c] - v1[c] for c in range(3)]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c * arcLengthAround for c in vector.normalise(d1)]
                d1InnerAroundList.append(dx_ds1)
            # Account for d1 of node sitting on half haustrum
            d1 = vector.normalise([xAround[elementsCountAroundHalfHaustrum][c] -
                                   xAround[elementsCountAroundHalfHaustrum - 1][c] for c in range(3)])
            dx_ds1 = [c * arcLengthAround for c in d1]
            d1InnerAroundList.append(dx_ds1)

        if d1InnerAroundList:
            d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, d1InnerAroundList, fixStartDerivative=True)
            d1TCEdge = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC * 0.5)],
                                           vector.magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 - 1)]))
            d1Transition = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 1)],
                                               vector.magnitude(d1Smoothed[int(elementsCountAroundTC * 0.5 + 2)]))
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

def getElementIdxOfOstiumBoundary(centrePosition, trackSurfaceOstium, ostiumDiameter):
    """
    Finds the element indices of the boundaries of elements on tracksurface that surround
    the ostium. Indices based on numbering for elements around and along tracksurface.
    Boundary lies on xi=0 of elements on left and bottom boundaries and xi = 1 for right and
    top boundaries.
    :param centrePosition: surface description for centre of ostium.
    :param trackSurfaceOstium: surface description for tracksurface.
    :param ostiumDiameter: Diameter of ostium.
    :return: element indices on the left, right, bottom and top boundaries around tracksurface.
    """

    elementsAroundTrackSurface = trackSurfaceOstium.elementsCount1
    elementsAlongTrackSurface = trackSurfaceOstium.elementsCount2
    ei1 = centrePosition.e1
    ei2 = centrePosition.e2
    xi1 = centrePosition.xi1
    xi2 = centrePosition.xi2
    xCentre, d1Centre, d2Centre = trackSurfaceOstium.evaluateCoordinates(centrePosition, derivatives=True)

    # Left boundary
    leftPositionOfCentreElement = TrackSurfacePosition(ei1, ei2, 0, xi2)
    xLeft, d1Left, _ = trackSurfaceOstium.evaluateCoordinates(leftPositionOfCentreElement, derivatives=True)
    distxLeftToxCentre = interp.computeCubicHermiteArcLength(xLeft, d1Left, xCentre, d1Centre, False)
    remainingLength = ostiumDiameter * 0.5 - distxLeftToxCentre
    xCurrent = xLeft
    d1Current = d1Left

    for n1 in range(ei1, -1, -1):
        if remainingLength > 0.0:
            prevPosition = TrackSurfacePosition(n1-1, ei2, 0, xi2)
            xPrev, d1Prev, _ = trackSurfaceOstium.evaluateCoordinates(prevPosition, derivatives=True)
            distPrevToxCurrent = interp.computeCubicHermiteArcLength(xPrev, d1Prev, xCurrent, d1Current, False)
            remainingLength -= distPrevToxCurrent
            xCurrent = xPrev
            d1Current = d1Prev
        else:
            ei1Left = n1
            break

    # Right boundary
    rightPositionOfCentreElement = TrackSurfacePosition(ei1, ei2, 1.0, xi2)
    xRight, d1Right, _ = trackSurfaceOstium.evaluateCoordinates(rightPositionOfCentreElement, derivatives=True)
    distxCentreToxRight = interp.computeCubicHermiteArcLength(xCentre, d1Centre, xRight, d1Right, False)
    remainingLength = ostiumDiameter * 0.5 - distxCentreToxRight
    xCurrent = xRight
    d1Current = d1Right

    for n1 in range(ei1, elementsAroundTrackSurface):
        if remainingLength > 0.0:
            nextPosition = TrackSurfacePosition(n1+1, ei2, 1.0, xi2)
            xNext, d1Next, _ = trackSurfaceOstium.evaluateCoordinates(nextPosition, derivatives=True)
            distxCurrentToxNext = interp.computeCubicHermiteArcLength(xCurrent, d1Current, xNext, d1Next, False)
            remainingLength -= distxCurrentToxNext
            xCurrent = xNext
            d1Current = d1Next
        else:
            ei1Right = n1
            break

    # Bottom boundary
    bottomPositionOfCentreElement = TrackSurfacePosition(ei1, ei2, xi1, 0)
    xBottom, _, d2Bottom = trackSurfaceOstium.evaluateCoordinates(bottomPositionOfCentreElement, derivatives=True)
    distxBottomToxCentre = interp.computeCubicHermiteArcLength(xBottom, d2Bottom, xCentre, d2Centre, False)
    remainingLength = ostiumDiameter * 0.5 - distxBottomToxCentre
    xCurrent = xBottom
    d2Current = d2Bottom

    for n2 in range(ei2, -1, -1):
        if remainingLength > 0.0:
            prevPosition = TrackSurfacePosition(ei1, n2 - 1, xi1, 0)
            xPrev, _, d2Prev = trackSurfaceOstium.evaluateCoordinates(prevPosition, derivatives=True)
            distPrevToxCurrent = interp.computeCubicHermiteArcLength(xPrev, d2Prev, xCurrent, d2Current, False)
            remainingLength -= distPrevToxCurrent
            xCurrent = xPrev
            d2Current = d2Prev
        else:
            ei2Bottom = n2
            break

    # Top boundary
    topPositionOfCentreElement = TrackSurfacePosition(ei1, ei2, xi1, 1.0)
    xTop, _, d2Top = trackSurfaceOstium.evaluateCoordinates(topPositionOfCentreElement, derivatives=True)
    distxCentreToxTop = interp.computeCubicHermiteArcLength(xCentre, d2Centre, xTop, d2Top, False)
    remainingLength = ostiumDiameter * 0.5 - distxCentreToxTop
    xCurrent = xTop
    d2Current = d2Top

    for n2 in range(ei2, elementsAlongTrackSurface):
        if remainingLength > 0.0:
            nextPosition = TrackSurfacePosition(ei1, n2+1, xi1, 1.0)
            xNext, _, d2Next = trackSurfaceOstium.evaluateCoordinates(nextPosition, derivatives=True)
            distxCurrentToxNext = interp.computeCubicHermiteArcLength(xCurrent, d2Current, xNext, d2Next, False)
            remainingLength -= distxCurrentToxNext
            xCurrent = xNext
            d2Current = d2Next
        else:
            ei2Top = n2
            break

    return ei1Left, ei1Right, ei2Bottom, ei2Top

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

def createCecumMesh3d(region, options, centralPath, nextNodeIdentifier, nextElementIdentifier):
    """

    """
    segmentCount = options['Number of segments']
    startPhase = 0.0
    elementsCountAroundTC = options['Number of elements around tenia coli']
    elementsCountAroundHaustrum = options['Number of elements around haustrum']
    elementsCountAlongSegment = options['Number of elements along segment']
    elementsCountThroughWall = options['Number of elements through wall']
    cornerInnerRadiusFactor = options['Corner inner radius factor']
    haustrumInnerRadiusFactor = options['Haustrum inner radius factor']
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

    ostiumOptions = options['Ileocecal junction']
    ostiumSettings = ostiumOptions.getScaffoldSettings()

    zero = [0.0, 0.0, 0.0]
    firstNodeIdentifier = nextNodeIdentifier
    firstElementIdentifier = nextElementIdentifier
    startNode = nextNodeIdentifier
    startElement = nextElementIdentifier

    fm = region.getFieldmodule()
    coordinates = findOrCreateFieldCoordinates(fm)
    cache = fm.createFieldcache()

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
        annotationGroupsThroughWall = [[mucosaGroup],
                                       [submucosaGroup],
                                       [circularMuscleGroup],
                                       [longitudinalMuscleGroup]]

    # Sample central path along cecum
    # print(len(centralPath.cxGroups))
    cecumLength = centralPath.arcLengthOfGroupsAlong[0]
    cx = centralPath.cxGroups[0]
    cd1 = centralPath.cd1Groups[0]
    cd2 = centralPath.cd2Groups[0]
    cd12 = centralPath.cd12Groups[0]

    cxIleum = centralPath.cxGroups[1]
    cd1Ileum = centralPath.cd1Groups[1]
    cd2Ileum = centralPath.cd2Groups[1]
    cd3Ileum = centralPath.cd3Groups[1]
    cd12Ileum = centralPath.cd12Groups[1]
    cd13Ileum = centralPath.cd13Groups[1]

    # xBranchPt = centralPath.xBranchPt
    d2BranchPt = centralPath.d2BranchPt

    # smoothd1 = interp.smoothCubicHermiteDerivativesLine(cx, cd1, fixAllDirections=True,
    #                                                     magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
    sxCecum, sd1Cecum, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
    sd2Cecum, sd12Cecum = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)

    # Calculate segment length
    segmentLength = cecumLength / segmentCount

    # Generate variation of radius & tc width along length
    innerRadiusAlongCecum = []
    dInnerRadiusAlongCecum = []
    tcWidthAlongCecum = []

    closedProximalEnd = True

    innerRadiusListCP = [vector.magnitude(c) for c in cd2]
    dInnerRadiusListCP = []
    for n in range(len(innerRadiusListCP) - 1):
        dInnerRadiusListCP.append(innerRadiusListCP[n + 1] - innerRadiusListCP[n])
    dInnerRadiusListCP.append(innerRadiusListCP[-1] - innerRadiusListCP[-2])
    innerRadiusAlongElementList, dInnerRadiusAlongElementList = \
        interp.interpolateSampleCubicHermite(innerRadiusListCP, dInnerRadiusListCP, se, sxi, ssf)

    for n2 in range(elementsCountAlongSegment * segmentCount + 1):
        xi = 1 / (elementsCountAlongSegment * segmentCount) * n2

        radius = innerRadiusAlongElementList[n2]
        innerRadiusAlongCecum.append(radius)
        dRadius = dInnerRadiusAlongElementList[n2]
        dInnerRadiusAlongCecum.append(dRadius)
        tcWidth = interp.interpolateCubicHermite([startTCWidth], [startTCWidthDerivative],
                                                 [endTCWidth], [endTCWidthDerivative], xi)[0]
        tcWidthAlongCecum.append(tcWidth)

    haustrumInnerRadiusFactorAlongCecum = [haustrumInnerRadiusFactor] * (elementsCountAlong + 1)

    xToSample = []
    d1ToSample = []
    d2ToSample = []

    elementsCountAroundHalfHaustrum = int((elementsCountAroundTC + elementsCountAroundHaustrum) * 0.5)

    # Create object
    colonSegmentTubeMeshInnerPoints = ColonSegmentTubeMeshInnerPoints(
        region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment,
        tcCount, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor,
        segmentLength, wallThickness, cornerInnerRadiusFactor, haustrumInnerRadiusFactorAlongCecum,
        innerRadiusAlongCecum, dInnerRadiusAlongCecum, tcWidthAlongCecum, startPhase)

    # Create annotation
    cecumGroup = AnnotationGroup(region, get_cecum_term("caecum"))
    annotationGroupsAlong = []
    for i in range(elementsCountAlong):
        annotationGroupsAlong.append([cecumGroup])

    for nSegment in range(segmentCount):
        # Make regular segments
        xInner, d1Inner, d2Inner, transitElementList, segmentAxis, annotationGroupsAround \
            = colonSegmentTubeMeshInnerPoints.getColonSegmentTubeMeshInnerPoints(nSegment)

        # Replace first half of first segment with apex and sample along apex and second half of segment
        if nSegment == 0:
            xFirstSegmentSampled, d1FirstSegmentSampled, d2FirstSegmentSampled, d1FirstDirectionVector = \
                getApexSegmentForCecum(xInner, d1Inner, d2Inner, elementsCountAroundHalfHaustrum,
                                       elementsCountAroundTC, elementsCountAround, elementsCountAlongSegment,
                                       tcCount)

            xToSample += xFirstSegmentSampled
            d1ToSample += d1FirstSegmentSampled
            d2ToSample += d2FirstSegmentSampled
        else:
            xInnerExtrude = []
            for n in range(len(xInner)):
                xInnerExtrude.append([xInner[n][0], xInner[n][1], xInner[n][2] + segmentLength * nSegment])
            xToSample += xInnerExtrude[elementsCountAround:]
            d1ToSample += d1Inner[elementsCountAround:]
            d2ToSample += d2Inner[elementsCountAround:]

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

    # Project reference point for warping onto central path
    sxRefList, sd1RefList, sd2ProjectedListRef, zRefList = \
        tubemesh.getPlaneProjectionOnCentralPath(xToWarp, elementsCountAround, elementsCountAlong,
                                                 cecumLength, sxCecum, sd1Cecum, sd2Cecum, sd12Cecum)

    # Warp points
    xWarpedList, d1WarpedList, d2WarpedList, d3WarpedUnitList = \
        tubemesh.warpSegmentPoints(xToWarp, d1ToWarp, d2ToWarp, segmentAxis, sxRefList, sd1RefList,
                                   sd2ProjectedListRef, elementsCountAround, elementsCountAlong, zRefList)

    # Create coordinates and derivatives
    wallThicknessList = [wallThickness] * (elementsCountAlong + 1)

    xList, d1List, d2List, d3List, curvatureList = \
        tubemesh.getCoordinatesFromInner(xWarpedList, d1WarpedList,d2WarpedList, d3WarpedUnitList,
                                         wallThicknessList, relativeThicknessList, elementsCountAround,
                                         elementsCountAlong, elementsCountThroughWall, transitElementList)

    # Deal with multiple nodes at end point for closed proximal end
    xApexInner = xList[0]
    # arclength between apex point and corresponding point on next face
    mag = interp.getCubicHermiteArcLength(xList[0], d2List[0],
                                          xList[elementsCountAround * (elementsCountThroughWall + 1)],
                                          d2List[elementsCountAround * (elementsCountThroughWall + 1)])
    d2ApexInner = vector.setMagnitude(sd2Cecum[0], mag * 0.5)
    d1ApexInner = vector.crossproduct3(sd1Cecum[0], d2ApexInner)
    d1ApexInner = vector.setMagnitude(d1ApexInner, mag * 0.5)
    d3ApexUnit = vector.normalise(vector.crossproduct3(vector.normalise(d1ApexInner),
                                                       vector.normalise(d2ApexInner)))
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
        tubeTCWidthList = colonSegmentTubeMeshInnerPoints.getTubeTCWidthList()
        xCecum, d1Cecum, d2Cecum, d3Cecum, annotationArrayAround = getTeniaColi(
            region, xCecum, d1Cecum, d2Cecum, d3Cecum, curvatureList, tcCount, elementsCountAroundTC,
            elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
            tubeTCWidthList, tcThickness, sxRefList, annotationGroupsAround, closedProximalEnd)

        nextNodeIdentifier, nextElementIdentifier, allAnnotationGroups = createNodesAndElementsTeniaColi(
            region, xCecum, d1Cecum, d2Cecum, d3Cecum, xFlat, d1Flat, d2Flat, xOrgan, d1Organ, d2Organ, None,
            elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
            tcCount, annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
            closedProximalEnd)

    else:
        nextNodeIdentifier, nextElementIdentifier, allAnnotationGroups = tubemesh.createNodesAndElements(
            region, xCecum, d1Cecum, d2Cecum, d3Cecum, xFlat, d1Flat, d2Flat, xOrgan, d1Organ, d2Organ, None,
            elementsCountAround, elementsCountAlong, elementsCountThroughWall,
            annotationGroupsAround, annotationGroupsAlong, annotationGroupsThroughWall,
            firstNodeIdentifier, firstElementIdentifier, useCubicHermiteThroughWall, useCrossDerivatives,
            closedProximalEnd)

    nodeIdentifierCecum = nextNodeIdentifier
    for n2 in range(len(cx)):
        node = nodes.createNode(nextNodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
        nextNodeIdentifier += 1

    #################
    # Create elements
    #################

    mesh = fm.findMeshByDimension(1)
    cubicHermiteBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    eft = mesh.createElementfieldtemplate(cubicHermiteBasis)
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
    result = elementtemplate.defineField(coordinates, -1, eft)

    elementIdentifier = nextElementIdentifier
    for e in range(len(cx) - 1):
        element = mesh.createElement(elementIdentifier, elementtemplate)
        element.setNodesByIdentifier(eft, [nodeIdentifierCecum + e, nodeIdentifierCecum + 1 + e])
        elementIdentifier = elementIdentifier + 1

    cx = centralPath.cxGroups[1]
    cd1 = centralPath.cd1Groups[1]
    cd2 = centralPath.cd2Groups[1]

    nodeIdentifierIleum = nextNodeIdentifier
    for n2 in range(len(cx)):
        node = nodes.createNode(nextNodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
        nextNodeIdentifier += 1

    for e in range(2):
        element = mesh.createElement(elementIdentifier, elementtemplate)
        element.setNodesByIdentifier(eft, [nodeIdentifierIleum + e, nodeIdentifierIleum + 1 + e])
        elementIdentifier = elementIdentifier + 1

    # Add ostium on track surface
    elementsAroundTrackSurface = elementsCountAroundHaustrum
    elementsAlongTrackSurface = elementsCountAlongSegment

    # Find region where ostium sits
    # angle between d2 of branch point and vector between branch point and 1st point on ileum
    dV = [cxIleum[0][c] - cxIleum[-1][c] for c in range(3)]
    ostiumPositionAngleAround = math.acos(vector.dotproduct(dV, d2BranchPt)/
                                          (vector.magnitude(dV) * vector.magnitude(d2BranchPt)))
    sectorIdx = ostiumPositionAngleAround // (2 * math.pi / tcCount)
    sectorStartAngle = sectorIdx * (2 * math.pi / tcCount)

    startIdxElementsAround = int((elementsCountAroundHaustrum + elementsCountAroundTC) * sectorIdx +
                                 elementsCountAroundTC * 0.5)

    segmentIdx = int(centralPath.arcLengthToBranchPt // segmentLength)

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
                    startIdxElementsAround + int(elementsCountAroundTC * (0.5 if tcThickness > 0.0 else 0)) + \
                    n1 + (elementsCountAround * n3) + (elementsCountAround * elementsCountThroughWall +
                                                       elementsCountAroundTC * (tcCount if tcThickness > 0.0 else 0)) \
                    * n2 + startElement - 1 + segmentIdx * ((elementsCountAround * n3) +
                                                            (elementsCountAround * elementsCountThroughWall +
                                                             elementsCountAroundTC *
                                                             (tcCount if tcThickness > 0.0 else 0)) *
                                                            elementsCountAlongSegment)

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

    # # Visualise track surface
    # for n1 in range(len(xTrackSurface)):
    #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTrackSurface[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d2TrackSurface[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d1TrackSurface[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    #     nextNodeIdentifier += 1

    trackSurfaceOstium = TrackSurface(elementsAroundTrackSurface, elementsAlongTrackSurface,
                                      xTrackSurface, d1TrackSurface, d2TrackSurface)

    # Find centre position
    # track along ileum path and since cxIleum[1] could be above or below the track surface, we check both side to determine direction to track
    # At each point, find the nearest position and take the diff between nearest point to the point in line, keep tracking till diff is close to zero.
    xTol = 1.0E-6
    arcStart = 0.0
    arcEnd = centralPath.arcLengthOfGroupsAlong[1]
    nearestPosition = trackSurfaceOstium.findNearestPosition(cxIleum[0])
    xNearestStart = trackSurfaceOstium.evaluateCoordinates(nearestPosition, derivatives=False)
    distStart = vector.magnitude([cxIleum[0][c] - xNearestStart[c] for c in range(3)])
    nearestPosition = trackSurfaceOstium.findNearestPosition(cxIleum[-1])
    xNearestEnd = trackSurfaceOstium.evaluateCoordinates(nearestPosition, derivatives=False)
    distEnd = vector.magnitude([cxIleum[-1][c] - xNearestEnd[c] for c in range(3)])

    for iter in range(100):
        arcDistance = (arcStart + arcEnd) * 0.5
        x, d1 = interp.getCubicHermiteCurvesPointAtArcDistance(cxIleum, cd1Ileum, arcDistance)[0:2]
        nearestPosition = trackSurfaceOstium.findNearestPosition(x)
        xNearest = trackSurfaceOstium.evaluateCoordinates(nearestPosition, derivatives=False)
        dist = vector.magnitude([x[c] - xNearest[c] for c in range(3)])

        if abs(distStart - distEnd) > xTol:
            if distStart < distEnd:
                arcEnd = arcDistance
                distEnd = dist
            else:
                arcStart = arcDistance
                distStart = dist

        else:
            xCentre, d1Centre, d2Centre = trackSurfaceOstium.evaluateCoordinates(nearestPosition, derivatives=True)
            normAxis = vector.normalise([-d for d in d1])
            eIdx = interp.getNearestPointIndex(cxIleum, xCentre) - 1
            arcLenghtSum = 0.0
            for e in range(eIdx):
                arcLenghtSum += interp.getCubicHermiteArcLength(cxIleum[e], cd1Ileum[e],
                                                                cxIleum[e + 1], cd1Ileum[e + 1])
            xi = (arcDistance - arcLenghtSum)/\
                 interp.getCubicHermiteArcLength(cxIleum[eIdx], cd1Ileum[eIdx], cxIleum[eIdx + 1], cd1Ileum[eIdx + 1])
            d2Centre = interp.interpolateCubicHermite(cd2Ileum[eIdx], cd12Ileum[eIdx], cd2Ileum[eIdx + 1], cd12Ileum[eIdx + 1], xi)
            break
    if iter > 98:
        print('Search for ileum entry centre - Max iters reached:', iter)

    # node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    # cache.setNode(node)
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xNearest)
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, 0.0])
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, 0.0, 0.0])
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    # print('x', nextNodeIdentifier)
    # nextNodeIdentifier += 1

    # node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    # cache.setNode(node)
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxIleum[eIdx])
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, 0.0])
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2Ileum[eIdx])
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    # nextNodeIdentifier += 1
    #
    #
    # node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    # cache.setNode(node)
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxIleum[eIdx + 1])
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, 0.0])
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2Ileum[eIdx + 1])
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    # nextNodeIdentifier += 1
    #
    # node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    # cache.setNode(node)
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xCentre)
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Centre)
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Centre)
    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, normAxis)
    # nextNodeIdentifier += 1

    ostiumSettings['Number of elements around ostium'] = elementsCountAlongSegment
    elementsCountAroundOstium = ostiumSettings['Number of elements around ostium']

    fm = region.getFieldmodule()
    mesh = fm.findMeshByDimension(3)
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)

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

    centralPathIleum = CustomCentralPath(xPath, d1Path, d2Path, d3Path, d12Path, d13Path)

    # for n in range(2):
    #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xPath[n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Path[n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Path[n])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3Path[n])
    #     nextNodeIdentifier += 1

    nextNodeIdentifier, nextElementIdentifier, (o1_x, o1_d1, o1_d2, o1_d3, o1_NodeId, o1_Positions) = \
        generateOstiumMesh(region, ostiumSettings, trackSurfaceOstium, centralPathIleum,
                           startNodeIdentifier=nextNodeIdentifier, startElementIdentifier=nextElementIdentifier,
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
    sf = vector.magnitude(centralPathIleum.cd2Path[-1]) * 0.25
    for n1 in range(elementsCountAroundOstium):
        normD2 = vector.normalise(o1_d2[-1][n1])
        d1AnnulusNorm.append(normD2)
        d1AnnulusOuter.append(vector.setMagnitude(o1_d2[-1][n1], sf))
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
        d3 = vector.normalise(vector.crossproduct3(vector.normalise(d2AnnulusOuter[n]), d1AnnulusNorm[n]))
        d3Annulus.append(d3)
    annulusD2Curvature = findCurvatureAroundLoop(xAnnulusOuter, d2AnnulusOuter, d3Annulus)

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
    rotFrame = matrix.getRotationMatrixFromAxisAngle(d3ApexUnit, rotAngle)
    d1A = [rotFrame[j][0] * d1Cecum[1][0] + rotFrame[j][1] * d1Cecum[1][1] +
           rotFrame[j][2] * d1Cecum[1][2] for j in range(3)]

    rotAngle = math.pi
    rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(d3Annulus[0]), rotAngle)
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

    xPositionB = trackSurfaceOstium.findNearestPosition(xAnnulusOuter[0])
    xProportionB = trackSurfaceOstium.getProportion(xPositionB)
    derivativeB = d2B
    derivativeMagnitudeB = vector.magnitude(derivativeB)

    nx, nd1, nd2, nd3, proportions = \
        trackSurfaceOstium.createHermiteCurvePoints(
            xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1],
            rowsBelow + 1, derivativeStart=derivativeA, derivativeEnd=derivativeB)

    pxAlongMidLine, pd2AlongMidLine, pd1AlongMidLine = \
        trackSurfaceOstium.resampleHermiteCurvePointsSmooth(
            nx, nd1, nd2, nd3, proportions,
            derivativeMagnitudeEnd=derivativeMagnitudeB)[0:3]

    for n in range(len(pd1AlongMidLine)):
        pd1AlongMidLine[n] = [-d for d in pd1AlongMidLine[n]]

    # for n1 in range(len(nx)):
    #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd2[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, 0.0, 0.0])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    #     # print(n1, nextNodeIdentifier)
    #     nextNodeIdentifier += 1

    # Apex half
    # print(startRowIdx)
    # pxAlongMidLine, pd2AlongMidLine = \
    #     interp.sampleCubicHermiteCurvesSmooth(
    #         nx, nd2, rowsBelow + 1,
    #         # derivativeMagnitudeStart=vector.magnitude([nx[1][c] - nx[0][c] for c in range(3)]),
    #         derivativeMagnitudeEnd=vector.magnitude(d2B))[:2]

    # for n1 in range(len(pxAlongMidLine)):
    #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, pxAlongMidLine[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, pd2AlongMidLine[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, pd1AlongMidLine[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    #     nextNodeIdentifier += 1

    # Next haustrum half
    # nx = [xAnnulusOuter[int(elementsCountAlongSegment * 0.5)]]
    # nd2 = [d1AnnulusOuter[int(elementsCountAlongSegment * 0.5)]]
    # nd1 = [d2AnnulusOuter[int(elementsCountAlongSegment * 0.5)]]
    #
    # for n in range(elementsAlongTrackSurface - e2Top):
    #     n1 = e2Top + 1 + n
    #     idx = elementsCountAroundHalfHaustrum + n1 * (elementsCountAroundHaustrum + 1) - 1
    #     nx.append(xTrackSurface[idx])
    #     nd2.append(d2TrackSurface[idx])
    #     nd1.append(d1TrackSurface[idx])
    #
    # pxAlongMidLineBottom, pd2AlongMidLineBottom, pe, pxi, psf = \
    #     interp.sampleCubicHermiteCurvesSmooth(
    #         nx, nd2, rowsAbove + 1,
    #         derivativeMagnitudeStart=vector.magnitude(nd2[0]),
    #         derivativeMagnitudeEnd=vector.magnitude(nd2[-1]))
    #
    # pd1AlongMidLineBottom = nd1

    xPositionA = trackSurfaceOstium.findNearestPosition(xAnnulusOuter[int(elementsCountAlongSegment * 0.5)])
    xProportionA = trackSurfaceOstium.getProportion(xPositionA)
    xA, derivative2A, derivativeA = trackSurfaceOstium.evaluateCoordinates(xPositionA, derivatives=True)
    derivativeMagnitudeA = vector.magnitude(d2AnnulusOuter[int(elementsCountAlongSegment * 0.5)]) # derivativeA)

    xB = xTrackSurface[-elementsCountAroundHalfHaustrum]
    xPositionB = trackSurfaceOstium.findNearestPosition(xB)
    xProportionB = trackSurfaceOstium.getProportion(xPositionB)
    derivativeB = d2TrackSurface[-elementsCountAroundHalfHaustrum]

    nx, nd1, nd2, nd3, proportions = \
        trackSurfaceOstium.createHermiteCurvePoints(
            xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1],
            rowsAbove + 1, derivativeStart=derivativeA, derivativeEnd=derivativeB)

    pxAlongMidLineBottom, pd2AlongMidLineBottom, pd1AlongMidLineBottom = \
        trackSurfaceOstium.resampleHermiteCurvePointsSmooth(
            nx, nd1, nd2, nd3, proportions, derivativeMagnitudeStart=derivativeMagnitudeA,
            derivativeMagnitudeEnd=derivativeMagnitudeA)[0:3]

    for n in range(len(pd1AlongMidLineBottom)):
        pd1AlongMidLineBottom[n] = [-d for d in pd1AlongMidLineBottom[n]]

    # for n1 in range(len(pxAlongMidLineBottom)):
    #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #     cache.setNode(node)
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, pxAlongMidLineBottom[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, pd2AlongMidLineBottom[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, pd1AlongMidLineBottom[n1])
    #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
    #     nextNodeIdentifier += 1

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

                    if n2 < rowsBelow + 1:
                        xPositionB = trackSurfaceOstium.findNearestPosition(pxAlongMidLine[n2])
                        derivativeB = None

                    elif n2 >= rowsBelow + rowsOstium:
                        idx = n2 - (rowsBelow + rowsOstium) + 1
                        xPositionB = trackSurfaceOstium.findNearestPosition(pxAlongMidLineBottom[idx])
                        derivativeB = None

                    else:
                        xPositionB = trackSurfaceOstium.findNearestPosition(xAnnulusOuter[annulusIdx])
                        rotAngle = math.pi
                        rotFrame = matrix.getRotationMatrixFromAxisAngle(vector.normalise(d3Annulus[annulusIdx]),
                                                                         rotAngle)
                        derivativeB = [rotFrame[j][0] * d1AnnulusOuter[annulusIdx][0] +
                                       rotFrame[j][1] * d1AnnulusOuter[annulusIdx][1] +
                                       rotFrame[j][2] * d1AnnulusOuter[annulusIdx][2] for j in range(3)]

                    xProportionB = trackSurfaceOstium.getProportion(xPositionB)

                    # node = nodes.createNode(nextNodeIdentifier, nodetemplate)
                    # cache.setNode(node)
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTrackSurface[n2 * (elementsCountAroundHaustrum + 1)])
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, 0.0])
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, derivativeA)
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
                    # nextNodeIdentifier += 1
                    #
                    # node = nodes.createNode(nextNodeIdentifier, nodetemplate)
                    # cache.setNode(node)
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xB)
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [0.0, 0.0, 0.0])
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, derivativeB)
                    # coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
                    # nextNodeIdentifier += 1

                else:  # RHS
                    if n2 < rowsBelow + 1:
                        xPositionA = trackSurfaceOstium.findNearestPosition(pxAlongMidLine[n2])
                        derivativeA = None

                    elif n2 >= rowsBelow + rowsOstium:
                        idx = n2 - (rowsBelow + rowsOstium) + 1
                        xPositionA = trackSurfaceOstium.findNearestPosition(pxAlongMidLineBottom[idx])
                        derivativeA = None

                    else:
                        xPositionA = trackSurfaceOstium.findNearestPosition(xAnnulusOuter[-annulusIdx])
                        derivativeA = d1AnnulusOuter[-annulusIdx]
                    xProportionA = trackSurfaceOstium.getProportion(xPositionA)

                    xPositionB = trackSurfaceOstium.findNearestPosition(
                        xTrackSurface[n2 * (elementsCountAroundHaustrum + 1) + elementsCountAroundHaustrum])
                    xProportionB = trackSurfaceOstium.getProportion(xPositionB)
                    derivativeB = d1TrackSurface[n2 * (elementsCountAroundHaustrum + 1) + elementsCountAroundHaustrum]

                nx, nd1, nd2, nd3, proportions = \
                    trackSurfaceOstium.createHermiteCurvePoints(
                        xProportionA[0], xProportionA[1], xProportionB[0], xProportionB[1],
                        sideElements +(0 if (n2 < rowsBelow + 1 or n2 >= rowsBelow + rowsOstium) else -1),
                        derivativeStart=derivativeA, derivativeEnd=derivativeB)

                nx, nd1 = \
                    trackSurfaceOstium.resampleHermiteCurvePointsSmooth(
                        nx, nd1, nd2, nd3, proportions)[0:2]

                # for n in range(len(nx)):
                #     node = nodes.createNode(nextNodeIdentifier, nodetemplate)
                #     cache.setNode(node)
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, nx[n])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, nd1[n])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, [0.0, 0.0, 0.0])
                #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, [0.0, 0.0, 0.0])
                #     nextNodeIdentifier += 1

                if n2 == rowsBelow + 1 or n2 == rowsBelow + rowsOstium - 1:
                    mag = 0.5 * (vector.magnitude(d1AnnulusOuter[0]) + vector.magnitude(
                        d2AnnulusOuter[0])) if n2 == startRowIdx else vector.magnitude(d1AnnulusOuter[0])
                    if n == 0:
                        nd1[-1] = vector.setMagnitude(nd1[-1], mag)
                    else:
                        nd1[0] = vector.setMagnitude(nd1[0], mag)
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
        if n1 < elementsCountAroundHalfHaustrum - 1:
            for n2 in range(elementsCountAlongSegment):
                nxAlong.append(xAroundAlong[n2][n1])
                nd2Along.append(findDerivativeBetweenPoints(xAroundAlong[n2][n1], xAroundAlong[n2 + 1][n1]))
            nxAlong.append(xAroundAlong[n2 + 1][n1])
            nd2Along.append(d2TrackSurface[(elementsCountAroundHaustrum + 1) * elementsCountAlongSegment + n1])
            nd2Along = interp.smoothCubicHermiteDerivativesLine(nxAlong, nd2Along, fixStartDerivative=True,
                                                                fixEndDerivative=True)

            # Replace d2 on node along annulus LHS with d2Annulus
            if n1 == elementsCountAroundHalfHaustrum - 2:
                nd2Along[startRowIdx] = d2AnnulusOuter[1]
                nd2Along[endRowIdx] = d2AnnulusOuter[int(elementsCountAroundOstium * 0.5) - 1]

            xAlongAll.append(nxAlong)
            d2AlongAll.append(nd2Along)

            for n2 in range(len(nd2Along)):
                d2AroundAlong[n2][n1] = nd2Along[n2]

        elif n1 == elementsCountAroundHalfHaustrum - 1:
            for n2 in range(len(pxAlongMidLine)):
                d2AroundAlong[n2][n1] = pd2AlongMidLine[n2]
            for n in range(1, len(pxAlongMidLineBottom)):
                d2AroundAlong[n + endRowIdx][n1] = pd2AlongMidLineBottom[n]
            nxAlong = pxAlongMidLine + pxAlongMidLineBottom
            xAlongAll.append(nxAlong)
            d2AlongAll.append(pd2AlongMidLine + pd2AlongMidLineBottom)

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

            # Replace d2 on node along annulus RHS with d2Annulus
            if n1 == elementsCountAroundHalfHaustrum:
                nd2Along[startRowIdx] = [-d for d in d2AnnulusOuter[-1]]
                nd2Along[endRowIdx] = [-d for d in d2AnnulusOuter[int(elementsCountAroundOstium * 0.5) + 1]]

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
            d3Around.append(vector.normalise(
                vector.crossproduct3(vector.normalise(d1AroundAlong[n2][n1]),
                                     vector.normalise(d2AroundAlong[n2][n1]))))
        d3UnitAroundAlong.append(d3Around)

    # for n2 in range(len(xAroundAlong)):
    #     for n1 in range(len(xAroundAlong[n2])):
    #         node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1AroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2AroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3UnitAroundAlong[n2][n1])
    #         nextNodeIdentifier += 1

    # Calculate curvatures
    # Curvatures around
    d1Curvature = []
    if segmentIdx == 0:
        d1Curvature.append([0.0])

    for n2 in range(1 if segmentIdx == 0 else 0, elementsCountAlongSegment + 1):
        if n2 < startRowIdx or n2 == elementsCountAlongSegment:
            d1Curvature.append(findCurvatureAlongLine(xAroundAlong[n2], d1AroundAlong[n2], d3UnitAroundAlong[n2]))
        else:
            d1CurvatureLeft = findCurvatureAlongLine(xAroundAlong[n2][:int(0.5 * len(xAroundAlong[n2]))],
                                                     d1AroundAlong[n2][:int(0.5 * len(xAroundAlong[n2]))],
                                                     d3UnitAroundAlong[n2][:int(0.5 * len(xAroundAlong[n2]))])

            d1CurvatureRight = findCurvatureAlongLine(xAroundAlong[n2][int(0.5 * len(xAroundAlong[n2])):],
                                                      d1AroundAlong[n2][int(0.5 * len(xAroundAlong[n2])):],
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

        if n1 < elementsCountAroundHalfHaustrum - 1:
            for n2 in range(elementsCountAlongSegment):
                d3UnitAlong.append(d3UnitAroundAlong[n2][n1])
            d3UnitAlong.append(d3UnitAroundAlong[n2 + 1][n1])
            d2CurvatureAlong = findCurvatureAlongLine(xAlong, d2Along, d3UnitAlong)
            # Adjust for corners
            if n1 == elementsCountAroundHalfHaustrum - 2:
                d2CurvatureAlong[endRowIdx] = annulusD2Curvature[int(elementsCountAlongSegment * 0.5) - 1]

            for n2 in range(len(d2CurvatureAlong)):
                d2Curvature[n2][n1] = d2CurvatureAlong[n2]

        elif n1 == elementsCountAroundHalfHaustrum - 1:
            xAlong = []
            d2Along = []
            d3UnitAlong = []
            for n2 in range(len(pd2AlongMidLine)):
                xAlong.append(xAroundAlong[n2][n1])
                d2Along.append(d2AroundAlong[n2][n1])
                d3UnitAlong.append(d3UnitAroundAlong[n2][n1])
            d2CurvatureAlong = findCurvatureAlongLine(xAlong, d2Along, d3UnitAlong)
            for n2 in range(len(pd2AlongMidLine)):
                d2Curvature[n2][n1] = d2CurvatureAlong[n2]

            # Calculate curvature for node on edge
            xAlong = []
            d2Along = []
            d3UnitAlong = []
            for n in range(1, len(pd2AlongMidLineBottom)):
                nIdx = n + endRowIdx
                xAlong.append(xAroundAlong[nIdx][n1])
                d2Along.append(d2AroundAlong[nIdx][n1])
                d3UnitAlong.append(d3UnitAroundAlong[nIdx][n1])
            d2CurvatureAlong = findCurvatureAlongLine(xAlong, d2Along, d3UnitAlong)
            for n in range(len(pd2AlongMidLineBottom) - 1):
                nIdx = n + endRowIdx + 1
                d2Curvature[nIdx][n1] = d2CurvatureAlong[n]
        else:
            for n2 in range(elementsCountAlongSegment):
                d3UnitAlong.append(d3UnitAroundAlong[n2][n1 + (0 if (n2 < startRowIdx or n2 > endRowIdx) else -1)])
            d3UnitAlong.append(d3UnitAroundAlong[n2 + 1][n1])
            d2CurvatureAlong = findCurvatureAlongLine(xAlong, d2Along, d3UnitAlong)

            # Adjust for corners
            if n1 == elementsCountAroundHalfHaustrum:
                d2CurvatureAlong[endRowIdx] = annulusD2Curvature[int(elementsCountAlongSegment * 0.5) + 1]

            for n2 in range(len(d2CurvatureAlong)):
                if n2 < startRowIdx or n2 > endRowIdx:
                    n1Idx = n1
                else:
                    n1Idx = n1 - 1
                d2Curvature[n2][n1Idx] = d2CurvatureAlong[n2]

    # Slot in annulus points along the midline
    rotFrame = matrix.getRotationMatrixFromAxisAngle(d3Annulus[0], math.pi)
    rotD2 = [rotFrame[j][0] * d2AnnulusOuter[0][0] + rotFrame[j][1] * d2AnnulusOuter[0][1] +
             rotFrame[j][2] * d2AnnulusOuter[0][2] for j in range(3)]
    xAroundAlong[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), xAnnulusOuter[0])
    d1AroundAlong[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), rotD2)
    d2AroundAlong[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d1AnnulusOuter[0])
    d3UnitAroundAlong[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d3Annulus[0])
    d1Curvature[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), annulusD2Curvature[0])
    d2Curvature[startRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), annulusD2Curvature[0])

    idx = int(elementsCountAlongSegment * 0.5)
    xAroundAlong[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), xAnnulusOuter[idx])
    d1AroundAlong[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d2AnnulusOuter[idx])
    d2AroundAlong[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d1AnnulusOuter[idx])
    d3UnitAroundAlong[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), d3Annulus[idx])
    d1Curvature[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), annulusD2Curvature[idx])
    d2Curvature[endRowIdx].insert(int(elementsCountAroundHaustrum * 0.5), annulusD2Curvature[idx])

    # for n2 in range(len(xAroundAlong)):
    #     for n1 in range(len(xAroundAlong[n2])):
    #         node = nodes.createNode(nextNodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xAroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1AroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2AroundAlong[n2][n1])
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3UnitAroundAlong[n2][n1])
    #         nextNodeIdentifier += 1

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

                        if e2 == startRowIdx - 1 and e1 == elementsCountAroundHalfHaustrum - 2:  # LHS bottom
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])

                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e2 == startRowIdx - 1 and e1 == elementsCountAroundHalfHaustrum - 1:  # RHS bottom
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e2 == startRowIdx - 1 and e1 == elementsCountAroundHalfHaustrum:  # RHS bottom + 1
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [3, 7], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e2 == startRowIdx - 1 and e1 == elementsCountAroundHalfHaustrum - 3:  # LHS Bottom Left - 1
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            remapEftNodeValueLabel(eft1, [4, 8], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e2 == endRowIdx and e1 == elementsCountAroundHalfHaustrum - 2:  # end LHS
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e2 == endRowIdx and e1 == elementsCountAroundHalfHaustrum - 3:  # end LHS -1
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, [1]), (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e2 == endRowIdx and e1 == elementsCountAroundHalfHaustrum - 1:  # end RHS
                            scaleFactors = [-1.0]
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            setEftScaleFactorIds(eft1, [1], [])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS1,
                                                   [(Node.VALUE_LABEL_D_DS2, [1])])
                            remapEftNodeValueLabel(eft1, [2, 6], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])

                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

                        elif e2 == endRowIdx and e1 == elementsCountAroundHalfHaustrum:  # end RHS + 1
                            eft1 = eftfactory.createEftNoCrossDerivatives()
                            remapEftNodeValueLabel(eft1, [1, 5], Node.VALUE_LABEL_D_DS2,
                                                   [(Node.VALUE_LABEL_D_DS1, []), (Node.VALUE_LABEL_D_DS2, [])])
                            elementtemplateX.defineField(coordinates, -1, eft1)
                            elementtemplate1 = elementtemplateX

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
            if nAround == 0:
                endDerivativesMap[n3][nAround] = ((-1, 0, 0), (0, 1, 0), None)
            elif 1 <= nAround < int(elementsCountAroundOstium * 0.5):
                endDerivativesMap[n3][nAround] = ((0, 1, 0), (-1, 0, 0), None)
            elif nAround == int(elementsCountAroundOstium * 0.5):
                endDerivativesMap[n3][nAround] = (None, None, None)
            else:
                endDerivativesMap[n3][nAround] = ((0, -1, 0), (1, 0, 0), None)

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

    return allAnnotationGroups, nextNodeIdentifier, nextElementIdentifier #, nodeIdDistal, xDistal, d1Distal, \
           # d2Distal, d3Distal


class CecumCentralPath:
    """
    Generates sampled central path for cecum scaffold.
    """
    def __init__(self, region, centralPath, termsAlong=[None]):
        """
        :param region: Zinc region to define model in.
        :param centralPath: Central path subscaffold from meshtype_1d_path1
        :param termsAlong: Annotation terms along length of central path
        """
        # Extract length of each group along stomach from central path
        arcLengthOfGroupsAlong = []
        cxGroups = []
        cd1Groups = []
        cd2Groups = []
        cd3Groups = []
        cd12Groups = []
        cd13Groups = []

        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        pathNetworkMesh = centralPath.getConstructionObject()
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

        for termName in termsAlong:
            tmpGroup = tmpFieldmodule.findFieldByName(termName).castGroup() if termName else None
            tmpNodeset = tmpGroup.getNodesetGroup(tmpNodes) if tmpGroup else tmpNodes

            if termName == "caecum":
                for i in range(2):
                    cx, cd1, cd2, cd3, cd12, cd13 = get_nodeset_path_ordered_field_parameters(
                        tmpNodeset, tmpCoordinates,
                        [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                         Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                         Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                        networkSegments[i + 1].getNodeIdentifiers(), networkSegments[i + 1].getNodeVersions())

                    cxGroup += cx[(1 if i else 0):]
                    cd1Group += cd1[(1 if i else 0):]
                    cd2Group += cd2[(1 if i else 0):]
                    cd3Group += cd3[(1 if i else 0):]
                    cd12Group += cd12[(1 if i else 0):]
                    cd13Group += cd13[(1 if i else 0):]

                    if i == 0:
                        xbranchpt = cx[-1]
                        d2branchpt = cd2[-1]
                        arcLengthToBranchPt = 0.0
                        for n in range(len(cx) - 1):
                            arcLengthToBranchPt += interp.getCubicHermiteArcLength(cx[n], cd1[n], cx[n + 1], cd1[n + 1])

            elif termName == "ileum":
                cxGroup, cd1Group, cd2Group, cd3Group, cd12Group, cd13Group = get_nodeset_path_ordered_field_parameters(
                    tmpNodes, tmpCoordinates,
                    [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
                     Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3,
                     Node.VALUE_LABEL_D2_DS1DS2, Node.VALUE_LABEL_D2_DS1DS3],
                    networkSegments[0].getNodeIdentifiers(), networkSegments[0].getNodeVersions())

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
        self.arcLengthToBranchPt = arcLengthToBranchPt

class CustomCentralPath:
    """
    Generates sampled central path for part of central path.
    """
    def __init__(self, cx, cd1, cd2, cd3, cd12, cd13):
        self.cxPath = cx
        self.cd1Path = cd1
        self.cd2Path = cd2
        self.cd3Path = cd3
        self.cd12Path = cd12
        self.cd13Path = cd13
        