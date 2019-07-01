"""
Generates a single 3-D colon segment with a simple
mesentery mesh along a central line, with variable
numbers of elements around, along and through wall,
with variable radius and thickness along.
"""

import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.meshtype_3d_colonsegmentteniacoli1 import sampleHaustrum, getuListFromOuterMidLengthProfile
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import tubemesh

class MeshType_3d_colonsegmentsimplemesentery1(Scaffold_base):
    '''
    Generates a single 3-D colon segment with simple mesentery mesh
    with variable numbers of elements around, along the central line,
    and through wall. The cross-section profile of the colon
    segment is largely cylindrical due to the lack of a tenia coli.
    '''
    @staticmethod
    def getName():
        return '3D Colon Segment Simple Mesentery 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Mouse 1']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around mesenteric zone' : 2,
            'Number of elements around non-mesenteric zone' : 10,
            'Number of elements along segment' : 4,
            'Number of elements through wall' : 1,
            'Inner radius': 1.35,
            'Mesenteric zone width': 0.8,
            'Segment length': 1.5,
            'Wall thickness': 0.05,
            'Use cross derivatives' : False,
            'Use linear through wall' : True,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along segment' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around mesenteric zone',
            'Number of elements around non-mesenteric zone',
            'Number of elements along segment',
            'Number of elements through wall',
            'Inner radius',
            'Mesenteric zone width',
            'Segment length',
            'Wall thickness',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements along segment',
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if options['Number of elements around mesenteric zone'] < 2:
            options['Number of elements around mesenteric zone'] = 2
        if options['Number of elements around non-mesenteric zone'] < 4:
            options['Number of elements around non-mesenteric zone'] = 4
        for key in [
            'Number of elements around mesenteric zone',
            'Number of elements around non-mesenteric zone']:
            if options[key] % 2 > 0:
                options[key] = options[key] + 1
        for key in [
            'Inner radius',
            'Mesenteric zone width',
            'Segment length',
            'Wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Mesenteric zone width'] < 10.0*math.pi/180.0*options['Inner radius']:
           options['Mesenteric zone width'] = round(10.0*math.pi/180.0*options['Inner radius'], 2)
        if options['Mesenteric zone width'] > math.pi*0.5*options['Inner radius']:
           options['Mesenteric zone width'] = round(math.pi*0.5*options['Inner radius'], 2)

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAroundMZ = options['Number of elements around mesenteric zone']
        elementsCountAroundNonMZ = options['Number of elements around non-mesenteric zone']
        elementsCountAround = elementsCountAroundMZ + elementsCountAroundNonMZ
        elementsCountAlongSegment = options['Number of elements along segment']
        elementsCountThroughWall = options['Number of elements through wall']
        radius = options['Inner radius']
        segmentLength = options['Segment length']
        widthMZ = options['Mesenteric zone width']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        segmentCount = 1

        cx = [ [ 0.0, 0.0, 0.0 ], [ segmentLength, 0.0, 0.0 ] ]
        cd1 = [ [ segmentLength, 0.0, 0.0 ], [ segmentLength, 0.0, 0.0 ] ]
        cd2 = [ [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ]
        cd12 = [ [0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]

        # Generate inner surface of a colon segment
        annotationGroups, annotationArray, transitElementList, uList, arcLengthOuterMidLength, xInner, d1Inner, d2Inner, segmentAxis = getColonSegmentInnerPointsNoTeniaColi(region, elementsCountAroundMZ,
           elementsCountAroundNonMZ, elementsCountAlongSegment, widthMZ, radius, segmentLength, wallThickness)

        # Generate tube mesh
        annotationGroups, nextNodeIdentifier, nextElementIdentifier, xList, d1List, d2List, d3List, sx, curvatureAlong, factorList = tubemesh.generatetubemesh(region,
            elementsCountAround, elementsCountAlongSegment, elementsCountThroughWall, segmentCount, cx, cd1, cd2, cd12,
            xInner, d1Inner, d2Inner, wallThickness, segmentAxis, segmentLength, useCrossDerivatives, useCubicHermiteThroughWall,
            annotationGroups, annotationArray, transitElementList, uList, arcLengthOuterMidLength)

        return annotationGroups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along segment']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)
        return meshrefinement.getAnnotationGroups()

def getColonSegmentInnerPointsNoTeniaColi(region, elementsCountAroundMZ, elementsCountAroundNonMZ, elementsCountAlongSegment,
    widthMZ, radius, segmentLength, wallThickness):
    """
    Generates a 3-D colon segment mesh with a simple mesentery
    (no tenia coli) with variable numbers of elements around,
    along the central line, and through wall. The colon segment
    has a cylindrical profile.
    :param elementsCountAroundMZ: Number of elements around mesenteric zone.
    :param elementsCountAroundNonMZ: Number of elements around non-mesenteric zone.
    :param elementsCountAlongSegment: Number of elements along colon segment.
    :param widthMZ: Width of mesenteric zone in flat preparation.
    :param radius: Inner radius of colon segment.
    :param segmentLength: Length of a colon segment.
    :param wallThickness: Thickness of wall.
    :return annotationGroups, annotationArray: annotationArray stores annotation
    names of elements around
    :return transitElementList: stores true if element around is an element that
    transits from tenia coli / mesenteric zone to haustrum / non-mesenteric zone.
    :return uList: List of xi for node around
    : return totalArcLengthOuterMidLength: total arclength of elements on outer
    surface along mid-length of segment.
    :return coordinates, derivatives on inner surface of a colon segment.
    """
    MZGroup = AnnotationGroup(region, 'mesenteric zone', FMANumber = 'FMANumber unknown', lyphID = 'Lyph ID unknown')
    NonMZGroup = AnnotationGroup(region, 'non-mesenteric zone', FMANumber = 'FMANumber unknown', lyphID = 'Lyph ID unknown')
    annotationGroups = [MZGroup, NonMZGroup]
    annotationArray = ['mesenteric zone']*int(elementsCountAroundMZ*0.5) + ['non-mesenteric zone']*elementsCountAroundNonMZ + ['mesenteric zone']*int(elementsCountAroundMZ*0.5)
    transitElementList = [0]*int(elementsCountAroundMZ*0.5) + [1] + [0]*int(elementsCountAroundNonMZ - 2) + [1] + [0]*int(elementsCountAroundMZ*0.5)

    # create nodes
    x = [ 0.0, 0.0, 0.0 ]
    d1 = [ 0.0, 0.0, 0.0 ]
    segmentAxis = [0.0, 0.0, 1.0]

    xMZ = []
    d1MZ = []
    xListBaseLayer = []
    d1ListBaseLayer = []
    xList = []
    d1List = []
    d2List = []
    uList = []

    #Set up profile
    sampleElementOut = 20
    nx, nd1 = createCirclePoints([ 0.0, 0.0, 0.0 ], [ radius, 0.0, 0.0 ], [ 0.0, radius, 0.0 ], sampleElementOut, startRadians = 0.0)

    # Sample half mesenteric zone into equally spaced nodes
    radiansAroundMZ = widthMZ/radius
    radiansPerElementAroundMZ = radiansAroundMZ / elementsCountAroundMZ
    for n1 in range(int(elementsCountAroundMZ*0.5 + 1)):
        radiansAround = radiansPerElementAroundMZ * n1
        cosRadiansAround = math.cos(radiansAround)
        sinRadiansAround = math.sin(radiansAround)
        x = [ radius*cosRadiansAround,
              radius*sinRadiansAround,
              0.0 ]
        d1 = [ -radius*sinRadiansAround*radiansPerElementAroundMZ,
                radius*cosRadiansAround*radiansPerElementAroundMZ,
                0.0]
        xMZ.append(x)
        d1MZ.append(d1)
    arcLengthPerMZ = radius*radiansPerElementAroundMZ

    # Sample half non-mesenteric zone into equally spaced nodes
    halfCircumference = math.pi*radius
    xNonMZ, d1NonMZ, arcLengthPerNonMZ, arcLengthPerTransition = sampleHaustrum(nx, nd1, xMZ[-1], d1MZ[-1], halfCircumference, widthMZ*0.5, elementsCountAroundNonMZ)

    xListBaseLayer = xList + xMZ + xNonMZ[1:]
    d1ListBaseLayer = d1List + d1MZ + d1NonMZ[1:]
    lengthHalfList = len(xListBaseLayer)

    # Reflect to get other half
    for n in range(1,int(lengthHalfList)-1):
        idx =  -n + lengthHalfList - 1
        x = xListBaseLayer[idx]
        d1 = d1ListBaseLayer[idx]
        xReflect = [x[0], -x[1], x[2]]
        d1Reflect = [-d1[0], d1[1], d1[2]]
        xListBaseLayer.append(xReflect)
        d1ListBaseLayer.append(d1Reflect)

    # Generate node along segment length
    lengthPerElementAlong = segmentLength / elementsCountAlongSegment
    d2 = [0.0, 0.0, lengthPerElementAlong]
    for n2 in range(elementsCountAlongSegment + 1):
        for n1 in range(elementsCountAroundNonMZ + elementsCountAroundMZ):
            z = lengthPerElementAlong*n2
            x = [xListBaseLayer[n1][0], xListBaseLayer[n1][1], z]
            d1 = d1ListBaseLayer[n1]
            xList.append(x)
            d1List.append(d1)
            d2List.append(d2)

    # Calculate uList for elements on outer surface along mid-length of segment
    uList, totalArcLengthOuterMidLength = getuListFromOuterMidLengthProfile(xListBaseLayer, d1ListBaseLayer, segmentAxis, wallThickness, transitElementList)

    return annotationGroups, annotationArray, transitElementList, uList, totalArcLengthOuterMidLength, xList, d1List, d2List, segmentAxis
