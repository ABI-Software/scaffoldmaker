"""
Generates a single 3-D colon segment mesh along a central
line, with variable numbers of elements around, along and
through wall, with variable radius and thickness along.
"""

import math
from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector
from scaffoldmaker.utils import zinc_utils
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

class MeshType_3d_colonsegmentteniacoli1(Scaffold_base):
    '''
    Generates a single 3-D colon segment mesh with variable
    numbers of tenia coli, elements around, along the central
    line, and through wall. The cross-section profile of the colon
    segment varies with species and is dependent on the number
    of tenia coli.
    Pig: 2 tenia coli, bow tie profile
    Human (Default): 3 tenia coli, triangular profile with rounded
    corners at the inter-haustral septa, and a clover
    profile in the intra-haustral region.
    '''
    @staticmethod
    def getName():
        return '3D Colon Segment Tenia Coli 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around tenia coli' : 2,
            'Number of elements around haustrum' : 8,
            'Number of elements along segment' : 4,
            'Number of elements through wall' : 1,
            'Inner radius': 1.0,
            'Corner inner radius factor': 0.5,
            'Haustrum inner radius factor': 0.5,
            'Segment length end derivative factor': 0.5,
            'Segment length mid derivative factor': 2.0,
            'Segment length': 1.5,
            'Number of tenia coli': 3,
            'Tenia coli width': 0.2,
            'Tenia coli thickness': 0.03,
            'Wall thickness': 0.02,
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
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Number of elements along segment',
            'Number of elements through wall',
            'Inner radius',
            'Corner inner radius factor',
            'Haustrum inner radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Segment length',
            'Number of tenia coli',
            'Tenia coli width',
            'Tenia coli thickness',
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
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along segment',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around tenia coli',
            'Number of elements along segment']:
            if options[key] < 2:
                options[key] = 2
        if options['Number of elements around haustrum'] < 4:
            options['Number of elements around haustrum'] = 4
        for key in [
            'Number of elements around tenia coli',
            'Number of elements around haustrum']:
            if options[key] % 2 > 0:
                options[key] = options[key] + 1
        for key in [
            'Inner radius',
            'Haustrum inner radius factor',
            'Segment length end derivative factor',
            'Segment length mid derivative factor',
            'Segment length',
            'Tenia coli thickness',
            'Wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Corner inner radius factor'] < 0.1:
            options['Corner inner radius factor'] = 0.1
        for key in [
            'Corner inner radius factor',
            'Segment length end derivative factor']:
            if options[key] > 1.0:
                options[key] = 1.0
        if options['Number of tenia coli'] < 2:
            options['Number of tenia coli'] = 2
        elif options['Number of tenia coli'] > 3:
            options['Number of tenia coli'] = 3
        if options['Tenia coli width'] < 0.2*options['Inner radius']:
            options['Tenia coli width'] = round(0.2*options['Inner radius'], 2)
        if options['Tenia coli width'] > round(math.sqrt(3)*0.5*options['Inner radius'],2):
            options['Tenia coli width'] = round(math.sqrt(3)*0.5*options['Inner radius'],2)

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAroundTC = options['Number of elements around tenia coli']
        elementsCountAroundHaustrum = options['Number of elements around haustrum']
        elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*3
        elementsCountAlongSegment = options['Number of elements along segment']
        elementsCountThroughWall = options['Number of elements through wall']
        radius = options['Inner radius']
        cornerInnerRadiusFactor = options['Corner inner radius factor']
        haustrumInnerRadiusFactor = options['Haustrum inner radius factor']
        segmentLengthEndDerivativeFactor = options['Segment length end derivative factor']
        segmentLengthMidDerivativeFactor = options['Segment length mid derivative factor']
        segmentLength = options['Segment length']
        tcCount = options['Number of tenia coli']
        tcWidth = options['Tenia coli width']
        tcThickness = options['Tenia coli thickness']
        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        haustraSegmentCount = 1

        cx = [ [ 0.0, 0.0, 0.0 ], [ segmentLength, 0.0, 0.0 ] ]
        cd1 = [ [ segmentLength, 0.0, 0.0 ], [ segmentLength, 0.0, 0.0 ] ]
        cd2 = [ [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ]
        cd12 = [ [0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]

        # Generate inner surface of a colon segment
        annotationGroups, annotationArray, transitElementList, uList, arcLengthOuterMidLength, xHaustraInner, d1HaustraInner, d2HaustraInner, haustraSegmentAxis = getColonSegmentInnerPointsTeniaColi(region, elementsCountAroundTC,
            elementsCountAroundHaustrum, elementsCountAlongSegment, tcCount, tcWidth, radius, cornerInnerRadiusFactor, haustrumInnerRadiusFactor,
            segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor, segmentLength, wallThickness)

        # Generate tube mesh
        annotationGroups, nextNodeIdentifier, nextElementIdentifier, xList, d1List, d2List, d3List, sx, curvatureAlong, factorList = tubemesh.generatetubemesh(region,
            elementsCountAround, elementsCountAlongSegment, elementsCountThroughWall, haustraSegmentCount, cx, cd1, cd2, cd12,
            xHaustraInner, d1HaustraInner, d2HaustraInner, wallThickness, haustraSegmentAxis, segmentLength,
            useCrossDerivatives, useCubicHermiteThroughWall, annotationGroups, annotationArray, transitElementList, uList, arcLengthOuterMidLength)

        # Generate tenia coli
        annotationGroupsTC, nextNodeIdentifier, nextElementIdentifier = getTeniaColi(region, nextNodeIdentifier, nextElementIdentifier,
           useCrossDerivatives, useCubicHermiteThroughWall, xList, d1List, d2List, d3List,
           elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment, elementsCountThroughWall, wallThickness,
           tcWidth, tcThickness, sx, curvatureAlong, factorList, uList, arcLengthOuterMidLength, segmentLength)

        annotationGroups += annotationGroupsTC

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

def getColonSegmentInnerPointsTeniaColi(region, elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongSegment, tcCount, tcWidth, radius,
    cornerInnerRadiusFactor, haustrumInnerRadiusFactor, segmentLengthEndDerivativeFactor, segmentLengthMidDerivativeFactor, segmentLength, wallThickness):
    """
    Generates a 3-D colon segment mesh with 3 tenia coli with variable
    numbers of elements around, along the central line, and through wall.
    The colon segment has a triangular profile with rounded corners
    at the inter-haustral septa, and a clover profile in the intra-haustral
    region.
    :param elementsCountAroundTC: Number of elements around each tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlongSegment: Number of elements along colon segment.
    :param tcCount: Number of tenia coli.
    :param tcWidth: Width of tenia coli.
    :param radius: Inner radius defined from center of triangular
    profile to vertex of the triangle.
    :param cornerInnerRadiusFactor: Roundness of triangular corners of
    inter-haustral septa. Factor is multiplied by inner radius
    to get a radius of curvature at the corners.
    :param haustrumInnerRadiusFactor: Factor is multiplied by inner
    radius to obtain radius of intersecting circles in the middle cross-section
    along a haustra segment.
    :param segmentLengthEndDerivativeFactor: Factor is multiplied by segment
    length to scale derivative along the end of a segment length.
    :param segmentLengthMidDerivativeFactor: Factor is multiplied by segment
    length to scale derivative along the mid length of the segment.
    :param segmentLength: Length of a colon segment.
    :param wallThickness: Thickness of wall.
    :return annotationGroups, annotationArray: annotationArray stores annotation
    names of elements around
    :return transitElementList: stores true if element around is an element that
    transits from tenia coli / mesenteric zone to haustrum / non-mesenteric zone.
    :return uList: List of xi for each node around mid-length haustra.
    : return totalArcLengthOuterMidLengthHaustra: Total arclength around outer
    surface on mid-length haustra.
    :return coordinates, derivatives on inner surface of a colon segment.
    """
    annotationGroups = []
    annotationArray = []
    transitElementListHaustrum = [0]*int(elementsCountAroundTC*0.5) + [1] + [0]*int(elementsCountAroundHaustrum - 2) + [1] + [0]*int(elementsCountAroundTC*0.5)
    transitElementList = transitElementListHaustrum*3

    assert tcCount == 3, 'getColonSegmentInnerPointsTeniaColi Only implemented for 3 tenia coli'

    # create nodes
    x = [ 0.0, 0.0, 0.0 ]
    dx_ds1 = [ 0.0, 0.0, 0.0 ]
    dx_ds2 = [ 0.0, 0.0, 0.0 ]
    dx_ds3 = [ 0.0, 0.0, 0.0 ]
    cornerRC = cornerInnerRadiusFactor*radius
    radiansRangeRC = [7*math.pi/4, 0.0, math.pi/4]
    sampleElementOut = 20
    segmentAxis = [0.0, 0.0, 1.0]
    haustrumRadius = (haustrumInnerRadiusFactor + 1)*radius

    xAround = []
    d1Around = []
    xHalfSetInterHaustra = []
    d1HalfSetInterHaustra = []
    nxHaustrum = []
    nd1Haustrum = []
    xHalfSetIntraHaustra = []
    d1HalfSetIntraHaustra = []
    xInnerRaw = []
    dx_ds2InnerRaw = []
    xInnerList = []
    dx_ds2InnerList = []
    xFinal = []
    d1Final = []
    d2Final = []
    uList = []

    # Inter-haustral segment
    # Set up profile
    for n1 in range(3):
        radiansAround = n1*2.0*math.pi / 3.0
        cosRadiansAround = math.cos(radiansAround)
        sinRadiansAround = math.sin(radiansAround)
        xc = [(radius - cornerRC) * cosRadiansAround, (radius - cornerRC) * sinRadiansAround, 0.0]

        for n in range(3):
            radiansRC = radiansAround + radiansRangeRC[n]
            cosRadiansRC = math.cos(radiansRC)
            sinRadiansRC = math.sin(radiansRC)
            x = [xc[0] + cornerRC*cosRadiansRC, xc[1] + cornerRC*sinRadiansRC, 0.0]
            xAround.append(x)
            d1 = [ cornerRC*math.pi/4.0 * -sinRadiansRC, cornerRC*math.pi/4.0 * cosRadiansRC, 0.0]
            d1Around.append(d1)

    xSample = xAround[1:9] +[xAround[0], xAround[1]]
    d1Sample = d1Around[1:9] +[d1Around[0], d1Around[1]]
    sx, sd1, se, sxi, _= interp.sampleCubicHermiteCurves(xSample, d1Sample, sampleElementOut)
    xLoop = sx[:-1]
    d1Loop = interp.smoothCubicHermiteDerivativesLoop(sx[:-1], sd1[:-1])

    # Calculate arc length
    arcLength = 0.0
    arcDistance = interp.getCubicHermiteArcLength(xLoop[0], d1Loop[0], xLoop[1], d1Loop[1])
    arcLength = arcDistance * sampleElementOut
    arcStart = 0.0
    arcEnd = arcLength/6.0

    # Find edge of TC
    arcDistanceTCEdge = findEdgeOfTeniaColi(xLoop, d1Loop, tcWidth, arcStart, arcEnd)

    # Sample TC into equally sized elements
    xTC, d1TC, _ = sampleTeniaColi(xLoop, d1Loop, arcDistanceTCEdge, elementsCountAroundTC)

    # Sample haustrum into equally sized elements
    xHaustrum, d1Haustrum, _, _ = sampleHaustrum(xLoop, d1Loop, xTC[-1], d1TC[-1], arcLength/6.0, arcDistanceTCEdge, elementsCountAroundHaustrum)

    xHalfSetInterHaustra = xHalfSetInterHaustra + xTC + xHaustrum[1:]
    d1HalfSetInterHaustra = d1HalfSetInterHaustra + d1TC + d1Haustrum[1:]

    # Intra-haustral segment
    # Set up profile
    xc = [(radius - cornerRC)* math.cos(0.0), (radius - cornerRC)*math.sin(0.0), 0.0]
    pt1 = [xc[0] + cornerRC*math.cos(0.0), xc[1] + cornerRC*math.sin(0.0), 0.0]
    xTC2 = radius* math.cos(2.0*math.pi/3.0)
    yTC2 = radius* math.sin(2.0*math.pi/3.0)
    originRC = (xTC2*xTC2 + yTC2*yTC2 - haustrumRadius*haustrumRadius) / (2*(-xTC2 - haustrumRadius))
    RC = haustrumRadius - originRC

    # Rotate to find originRC of 1st haustrum
    yTC1 = pt1[1]
    rotOriginRC = [ originRC*math.cos(-2.0/3.0*math.pi), originRC*math.sin(-2.0/3.0*math.pi), 0.0]

    thetaStart = math.asin((yTC1 + rotOriginRC[1]) / RC)
    thetaEnd = math.pi - math.asin((yTC2 + rotOriginRC[1])/ RC)
    thetaHalfHaustrum = (thetaEnd + thetaStart)*0.5
    concaveFactor = 0.15
    thetaConcave = (thetaEnd - thetaStart)*concaveFactor + thetaStart
    thetaCircularStart = (thetaEnd - thetaStart)*(concaveFactor + 0.5)*0.5 + thetaStart

    thetaSet = [thetaStart, thetaConcave]
    xConcave, _ = getCircleXandD1FromRadians(thetaSet, RC, rotOriginRC)
    d1 = [0.0, xConcave[1][1], 0.0]
    nxHaustrum = nxHaustrum + xConcave
    nd1Haustrum = nd1Haustrum + [d1, d1]
    thetaCircularSet = [thetaCircularStart, thetaHalfHaustrum]
    xCircular, d1Circular = getCircleXandD1FromRadians(thetaCircularSet, RC, rotOriginRC)
    nxHaustrum = nxHaustrum + xCircular
    nd1Haustrum = nd1Haustrum + d1Circular
    smoothd1 = interp.smoothCubicHermiteDerivativesLine(nxHaustrum, nd1Haustrum, fixStartDirection = True, fixEndDirection = True)

    # Find max x along path
    sxHaustrum, sd1Haustrum, _, _ , _ = interp.sampleCubicHermiteCurves(nxHaustrum, smoothd1, sampleElementOut)
    xVal = []
    for i in range(len(sxHaustrum)):
        xVal.append(sxHaustrum[i][0])
    maxValue = max(xVal)
    maxIndex = xVal.index(maxValue)
    yAtMaxValue = sxHaustrum[maxIndex][1]

    arcLength = 0.0
    for e in range(len(nxHaustrum)-1):
        arcDistance = interp.getCubicHermiteArcLength(nxHaustrum[e], smoothd1[e], nxHaustrum[e+1], smoothd1[e+1])
        arcLength = arcLength + arcDistance

    arcDistanceAtMaxValue = arcLength / sampleElementOut * maxIndex
    arcStart = 0.0 if yAtMaxValue > tcWidth*0.5 else arcDistanceAtMaxValue
    arcEnd = arcDistanceAtMaxValue if yAtMaxValue > tcWidth*0.5 else arcLength

    # Find edge of TC
    arcDistanceTCEdge = findEdgeOfTeniaColi(nxHaustrum, smoothd1, tcWidth, arcStart, arcEnd)

    # Sample TC into equally sized elements
    xTC, d1TC, arcLengthPerTC = sampleTeniaColi(nxHaustrum, smoothd1, arcDistanceTCEdge, elementsCountAroundTC)

    # Sample haustrum into equally sized elements
    xHaustrum, d1Haustrum, arcLengthPerHaustrum, arcLengthPerTransition = sampleHaustrum(nxHaustrum, smoothd1, xTC[-1], d1TC[-1], arcLength, arcDistanceTCEdge, elementsCountAroundHaustrum)

    xHalfSetIntraHaustra = xHalfSetIntraHaustra + xTC + xHaustrum[1:]
    d1HalfSetIntraHaustra = d1HalfSetIntraHaustra + d1TC + d1Haustrum[1:]

    # Calculate uList for elements lying around outer surface along mid-length haustra
    d2 = []
    for n in range(len(xHalfSetIntraHaustra)):
        d2.append([0.0, 0.0, 0.0])
    xInnerMidLengthHaustra, d1InnerMidLengthHaustra, _ = getFullProfileFromHalfHaustrum(xHalfSetIntraHaustra, d1HalfSetIntraHaustra, d2)
    uList, totalArcLengthOuterMidLengthHaustra = getuListFromOuterMidLengthProfile(xInnerMidLengthHaustra, d1InnerMidLengthHaustra, segmentAxis, wallThickness, transitElementList)

    # Sample arclength of haustra segment
    elementsCountAroundHalfHaustrum = int((elementsCountAroundTC + elementsCountAroundHaustrum)*0.5)

    for n1 in range(elementsCountAroundHalfHaustrum + 1):
        v1 = [xHalfSetInterHaustra[n1][0], xHalfSetInterHaustra[n1][1], 0.0]
        startArcLength = segmentLengthEndDerivativeFactor * segmentLength
        d1 = [ c*startArcLength for c in segmentAxis]
        v2 = [xHalfSetIntraHaustra[n1][0], xHalfSetIntraHaustra[n1][1], segmentLength/2]
        midArcLength = segmentLengthMidDerivativeFactor * segmentLength
        d2 = [c*midArcLength for c in segmentAxis]
        v3 = [xHalfSetInterHaustra[n1][0], xHalfSetInterHaustra[n1][1], segmentLength]
        d3 = [ c*startArcLength for c in segmentAxis]
        nx = [v1, v2, v3]
        nd1 = [d1, d2, d3]
        sx, sd1, se, sxi, _  = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAlongSegment)
        xInnerRaw.append(sx)
        dx_ds2InnerRaw.append(sd1)

    # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
    for n2 in range(elementsCountAlongSegment + 1):
        xAround = []
        unitdx_ds1Around = []
        d2Around = []

        for n1 in range(elementsCountAroundHalfHaustrum+1):
            x = xInnerRaw[n1][n2]
            xInnerList.append(x)
            dx_ds2 = dx_ds2InnerRaw[n1][n2]
            dx_ds2InnerList.append(dx_ds2)
            unitTangent = vector.normalise(dx_ds2)
            # Intra-Haustra segments
            if n1 == 0:
                unitdx_ds1 = vector.normalise(d1HalfSetIntraHaustra[n1])
            else: # points on clover
                if n2 <= int(elementsCountAlongSegment/2): # first half of segmentLength
                    axisRot = vector.crossproduct3(segmentAxis, unitTangent)
                elif n2 > int(elementsCountAlongSegment/2): # second half of segmentLength
                    axisRot = vector.crossproduct3(unitTangent, segmentAxis)
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, math.pi/2)
                rotNormal = [rotFrame[j][0]*unitTangent[0] + rotFrame[j][1]*unitTangent[1] + rotFrame[j][2]*unitTangent[2] for j in range(3)]
                unitdx_ds3 = vector.normalise(rotNormal)
                unitdx_ds1 = vector.crossproduct3(unitTangent, unitdx_ds3)
            xAround.append(x)
            d2Around.append(dx_ds2)
            unitdx_ds1Around.append(unitdx_ds1)

        if n2 > 0 and n2 < elementsCountAlongSegment:
            dx_ds1InnerAroundList = []
            if elementsCountAlongSegment%2 == 0 and n2 == int(elementsCountAlongSegment*0.5):
                dx_ds1InnerAroundList = dx_ds1InnerAroundList + d1HalfSetIntraHaustra
            else:
                for n1 in range(elementsCountAroundHalfHaustrum):
                    v1 = xAround[n1]
                    d1 = unitdx_ds1Around[n1]
                    v2 = xAround[n1+1]
                    d2 = unitdx_ds1Around[n1+1]
                    arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                    dx_ds1 = [c*arcLengthAround for c in d1]
                    dx_ds1InnerAroundList.append(dx_ds1)
                # Account for d1 of node sitting on half haustrum
                dx_ds1 = [c*arcLengthAround for c in unitdx_ds1Around[elementsCountAroundHalfHaustrum]]
                dx_ds1InnerAroundList.append(dx_ds1)
            d1Smoothed = interp.smoothCubicHermiteDerivativesLine(xAround, dx_ds1InnerAroundList, fixStartDerivative = True)
            d1TCEdge = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC*0.5)], vector.magnitude(d1Smoothed[int(elementsCountAroundTC*0.5 - 1)]))
            d1Transition = vector.setMagnitude(d1Smoothed[int(elementsCountAroundTC*0.5 + 1)], vector.magnitude(d1Smoothed[int(elementsCountAroundTC*0.5 + 2)]))
            d1Corrected = []
            d1Corrected = d1Corrected + d1Smoothed[:int(elementsCountAroundTC*0.5)]
            d1Corrected.append(d1TCEdge)
            d1Corrected.append(d1Transition)
            d1Corrected = d1Corrected + d1Smoothed[int(elementsCountAroundTC*0.5 + 2):]
        else:
            d1Corrected = d1HalfSetInterHaustra

        xAlongList, d1AlongList, d2AlongList = getFullProfileFromHalfHaustrum(xAround, d1Corrected, d2Around)
        xFinal = xFinal + xAlongList
        d1Final = d1Final + d1AlongList
        d2Final = d2Final + d2AlongList

    return annotationGroups, annotationArray, transitElementList, uList, totalArcLengthOuterMidLengthHaustra, xFinal, d1Final, d2Final, segmentAxis

def findEdgeOfTeniaColi(nx, nd1, tcWidth, arcStart, arcEnd):
    """
    Locate edge of tenia coli on a cubic hermite interpolated
    curve defined by nodes nx and derivatives nd1.
    :param nx: Coordinates of nodes along curve.
    :param nd1: Derivatives of nodes along curve.
    :param tcWidth: Width of tenia coli.
    :param arcStart: Lower limit of arc distance to search for
    edge of tenia coli.
    :param arcEnd: Upper limit of arc distance to search for
    edge of tenia coli.
    :return: arc distance covered by tenia coli.
    """
    xTol = 1.0E-6
    for iter in range(100):
        arcDistance = (arcStart + arcEnd)*0.5
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        diff = x[1] - tcWidth*0.5
        if abs(diff) > xTol:
            if diff < 0.0:
                arcStart = arcDistance
            else:
                arcEnd = arcDistance
        else:
            arcDistanceTCEdge = arcDistance
            break
    if iter > 99:
        print('Search for TC boundary - Max iters reached:',iter)

    return arcDistanceTCEdge

def sampleTeniaColi(nx, nd1, arcDistanceTCEdge, elementsCountAroundTC):
    """
    Get systematically spaced points and derivatives over the width of
    tenia coli over a cubic hermite interpolated curves defined by nodes
    nx and derivatives nd1 lying along half of a haustrum.
    :param nx: Coordinates of nodes along curve.
    :param nd1: Derivatives of nodes along curve.
    :param arcDistanceTCEdge: Arc distance covered by tenia coli.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :return: coordinates, derivatives of tenia coli, and arclength of
    each tenia coli element.
    """
    xTC = []
    d1TC = []
    arcDistancePerElementTC = arcDistanceTCEdge / (elementsCountAroundTC*0.5)
    for e in range(int(elementsCountAroundTC*0.5)+1):
        arcDistance = arcDistancePerElementTC * e
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        d1Scaled = vector.setMagnitude(d1, arcDistancePerElementTC)
        xTC.append(x)
        d1TC.append(d1Scaled)

    return xTC, d1TC, arcDistancePerElementTC

def sampleHaustrum(nx, nd1, xTCLast, d1TCLast, arcLength, arcDistanceTCEdge, elementsCountAroundHaustrum):
    """
    Get systematically spaced points and derivatives over the width of
    haustrum over a cubic hermite interpolated curves defined by nodes
    nx and derivatives nd1 lying along half of a haustrum.
    :param nx: Coordinates of nodes along curve.
    :param nd1: Derivatives of nodes along curve.
    :param xTCLast: Coordinates of node lying on edge of tenia coli.
    :param d1TCLast: Derivative of node lying on edge of tenia coli.
    :param arcLength: Arc length of curve.
    :param arcDistanceTCEdge: Arc distance covered by tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :return: coordinates, derivatives on haustrum.
    :return: arclength of haustrum element and transition elements.
    """
    xHaustrum = []
    d1Haustrum = []
    elementLengths = []
    length = arcLength - arcDistanceTCEdge
    elementsCountOut = int(elementsCountAroundHaustrum * 0.5)
    addLengthStart = 0.5 * vector.magnitude(d1TCLast)
    lengthFractionStart = 0.5
    proportionStart = 1.0
    elementLengthMid = (length - addLengthStart) / (elementsCountOut - 1.0 + proportionStart*lengthFractionStart )
    elementLengthProportionStart = proportionStart*lengthFractionStart*elementLengthMid
    elementLengthProportionEnd = elementLengthMid

    for e in range(elementsCountOut):
        xi = e/(elementsCountOut - 1)
        elementLengths.append(elementLengthMid)
    elementLengths[0] = addLengthStart + elementLengthProportionStart
    elementLengths[-1] = 0.0 + elementLengthProportionEnd
    haustrumElementLength = elementLengthMid
    transitionElementLength = elementLengths[0]

    arcDistance = arcDistanceTCEdge
    xHaustrum.append(xTCLast)
    d1Scaled = vector.setMagnitude(d1TCLast, elementLengths[0])
    d1Haustrum.append(d1Scaled)

    for e in range(elementsCountOut):
        arcDistance = arcDistance + elementLengths[e]
        x, d1, _, _ = interp.getCubicHermiteCurvesPointAtArcDistance(nx, nd1, arcDistance)
        d1Scaled = vector.setMagnitude(d1, elementLengths[e] if e > 0 else elementLengths[e+1])
        xHaustrum.append(x)
        d1Haustrum.append(d1Scaled)

    return xHaustrum, d1Haustrum, haustrumElementLength, transitionElementLength

def getCircleXandD1FromRadians(thetaSet, radius, origin):
    """
    Gets the coordinates and derivatives along a circular path
    based on an angular range.
    :param thetaSet: Lower and upper limit of theta.
    :param radius: Radius of circle.
    :param origin: Origin of circle.
    :return: coordinates, derivatives on lower and upper
    limit of thetaSet.
    """
    nx = []
    nd1 = []
    dTheta = thetaSet[1] - thetaSet[0]
    for n in range(2):
        theta = thetaSet[n]
        x = [radius*math.cos(theta) - origin[0],
            radius*math.sin(theta) - origin[1],
            0.0]
        d1 = [-radius*math.sin(theta)*dTheta,
            radius*math.cos(theta)*dTheta,
            0.0]
        nx.append(x)
        nd1.append(d1)

    return nx, nd1

def getFullProfileFromHalfHaustrum(xHaustrumHalfSet, d1HaustrumHalfSet, d2HaustrumHalfSet):
    """
    Gets the coordinates and derivatives of the entire profile
    using points from first half of the first sector. The first
    sector starts from the x-axis. The first half set of points
    are reflected across the x-axis followed by rotation to get
    the points in the second half of the first sector. The full
    set of points in the first sector are then rotated to obtain
    points in the other two sectors.
    :param xHaustrumHalfSet: Coordinates of points in first
    half of the first sector.
    :param d1HaustrumHalfSet: Derivatives of points around first
    half of the first sector.
    :param d2HaustrumHalfSet: Derivatives of points along first
    half of the first sector.
    :return: coordinates, derivatives of points over entire profile.
    """
    xHaustrumHalfSet2 = []
    d1HaustrumHalfSet2 = []
    d2HaustrumHalfSet2 = []
    xHaustrum = []
    d1Haustrum = []
    d2Haustrum = []
    xHaustra = []
    d1Haustra = []
    d2Haustra = []

    for n in range(1,len(xHaustrumHalfSet)):
        idx =  -n + len(xHaustrumHalfSet) - 1
        x = xHaustrumHalfSet[idx]
        d1 = d1HaustrumHalfSet[idx]
        xReflect = [x[0], -x[1], x[2]]
        d1Reflect = [d1[0], -d1[1], d1[2]]
        xRot = [xReflect[0]*math.cos(2/3*math.pi) - xReflect[1]*math.sin(2/3*math.pi),
                xReflect[0]*math.sin(2/3*math.pi) + xReflect[1]*math.cos(2/3*math.pi),
                xReflect[2]]
        d1Rot = [-(d1Reflect[0]*math.cos(2/3*math.pi) - d1Reflect[1]*math.sin(2/3*math.pi)),
                -(d1Reflect[0]*math.sin(2/3*math.pi) + d1Reflect[1]*math.cos(2/3*math.pi)),
                -d1Reflect[2]]
        d2 = d2HaustrumHalfSet[idx]
        d2Reflect = [d2[0], -d2[1], d2[2]]
        d2Rot = [(d2Reflect[0]*math.cos(2/3*math.pi) - d2Reflect[1]*math.sin(2/3*math.pi)),
                (d2Reflect[0]*math.sin(2/3*math.pi) + d2Reflect[1]*math.cos(2/3*math.pi)),
                d2Reflect[2]]
        xHaustrumHalfSet2.append(xRot)
        d1HaustrumHalfSet2.append(d1Rot)
        d2HaustrumHalfSet2.append(d2Rot)

    xHaustrum = xHaustrumHalfSet + xHaustrumHalfSet2
    d1Haustrum = d1HaustrumHalfSet + d1HaustrumHalfSet2
    d2Haustrum = d2HaustrumHalfSet + d2HaustrumHalfSet2

    # Rotate to get all 3 sectors
    xHaustra = xHaustra + xHaustrum[:-1]
    d1Haustra = d1Haustra + d1Haustrum[:-1]
    d2Haustra = d2Haustra + d2Haustrum[:-1]

    ang = [ 2/3*math.pi, -2/3*math.pi]
    for i in range(2):
        rotAng = ang[i]
        cosRotAng = math.cos(rotAng)
        sinRotAng = math.sin(rotAng)
        for n in range(len(xHaustrum)- 1):
            x = xHaustrum[n]
            d1 = d1Haustrum[n]
            d2 = d2Haustrum[n]
            x = [ x[0]*cosRotAng - x[1]*sinRotAng, x[0]*sinRotAng + x[1]*cosRotAng, x[2]]
            xHaustra.append(x)
            dx_ds1 = [ d1[0]*cosRotAng - d1[1]*sinRotAng, d1[0]*sinRotAng + d1[1]*cosRotAng, d1[2]]
            d1Haustra.append(dx_ds1)
            dx_ds2 = [ d2[0]*cosRotAng - d2[1]*sinRotAng, d2[0]*sinRotAng + d2[1]*cosRotAng, d2[2]]
            d2Haustra.append(dx_ds2)

    return xHaustra, d1Haustra, d2Haustra

def getuListFromOuterMidLengthProfile(xInnerMidLength, d1InnerMidLength,
    segmentAxis, wallThickness, transitElementList):
    """
    Gets a list of xi for flat and texture coordinates calculated
    from outer arclength of elements around mid-length of a segment.
    :param xInnerMidLength: Coordinates of points on inner surface around
    mid-length of segment.
    :param d1InnerMidLength: Derivatives of points on inner surface around
    mid-length of segment.
    :param segmentAxis: Axis of segment.
    :param wallThickness: Thickness of wall.
    :param transitElementList: stores true if element around is an element that
    transits from tenia coli / mesenteric zone to haustrum / non-mesenteric zone.
    :return uList: List containing xi for each point on outer surface extruded by
    wall thickness from points on inner surface.
    :return totalArcLengthOuterMidLength: Total arclength around outer surface on
    mid-length section of segment.
    """
    unitNormList = []
    d1OuterMidLength = []

    for n in range(len(xInnerMidLength)):
        unitNormList.append(vector.normalise(vector.crossproduct3(d1InnerMidLength[n], segmentAxis)))

    xOuterMidLength, curvatureInnerMidLength = tubemesh.getOuterCoordinatesAndCurvatureFromInner(xInnerMidLength, d1InnerMidLength, unitNormList, wallThickness, 0, len(xInnerMidLength), transitElementList)

    for n1 in range(len(xOuterMidLength)):
        factor = 1.0 + wallThickness * curvatureInnerMidLength[n1]
        dx_ds1 = [ factor*c for c in d1InnerMidLength[n1]]
        d1OuterMidLength.append(dx_ds1)

    arcLengthList = []
    for n1 in range(len(xOuterMidLength)):
        arcLengthPerElement = interp.getCubicHermiteArcLength(xOuterMidLength[n1], d1OuterMidLength[n1], xOuterMidLength[(n1+1) % len(xOuterMidLength)], d1OuterMidLength[(n1+1)% len(xOuterMidLength)])
        arcLengthList.append(arcLengthPerElement)

    # Total arcLength
    totalArcLengthOuterMidLength = 0.0
    for n1 in range(len(arcLengthList)):
        totalArcLengthOuterMidLength += arcLengthList[n1]

    uList = [0.0]
    arcDistance = 0
    for n in range(len(arcLengthList)):
        arcDistance = arcDistance + arcLengthList[n]
        xi = arcDistance / totalArcLengthOuterMidLength
        uList.append(xi)

    return uList, totalArcLengthOuterMidLength

def getTeniaColi(region, nodeIdentifier, elementIdentifier, useCrossDerivatives,
    useCubicHermiteThroughWall, xList, d1List, d2List, d3List,
    elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlong, elementsCountThroughWall,
    wallThickness, tcWidth, tcThickness, sxCentralLine, curvatureAlong, factorList, uList,
    arcLengthOuterMidLength, totalLengthAlong):
    """
    Create equally spaced nodes and elements for tenia coli over the outer
    surface of the haustra. Nodes of the tenia coli is sampled from a cubic
    hermite curve running through the left edge of tenia coli boundary on the
    outer surface of the haustra, a midpoint lying at a distance of tenia
    coli thickness above the midpoint of the boundary of tenia coli, and the
    right edge of tenia coli boundary.
    :param nodeIdentifier, elementIdentifier: First node and element identifier to
    use for tenia coli.
    :param xList, d1List, d2List, d3List: Coordinates and derivatives of nodes on haustra.
    :param elementsCountAroundTC: Number of elements around tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlong: Number of elements along scaffold.
    :param elementsCountThroughWall: Number of elements through wall.
    :param wallThickness: Thickness of wall.
    :param tcWidth: Width of tenia coli.
    :param tcThickness: Thickness of tenia coli at its thickest part.
    :param sxCentralLine: Coordinates sampled from central line.
    :param curvatureAlong: Curvatures along the colon for nodes on inner surface of colon.
    :param factorList: Factors used for scaling d2 to account for curvature along colon.
    :param uList: List of xi for each node around mid-length haustra.
    :param arcLengthOuterMidLength: Total arclength of elements around outer surface along
    mid-length of segment.
    :param totalLengthAlong: Total length of colon along center line.
    :return: annotationGroups, nodeIdentifier, elementIdentifier
    """

    fm = region.getFieldmodule()
    fm.beginChange()
    cache = fm.createFieldcache()
    coordinates = zinc_utils.getOrCreateCoordinateField(fm)

    tlGroup = AnnotationGroup(region, 'tenia libera', FMANumber = 'FMANumber unknown', lyphID = 'Lyph ID unknown')
    tmGroup = AnnotationGroup(region, 'tenia mesocolica', FMANumber = 'FMANumber unknown', lyphID = 'Lyph ID unknown')
    toGroup = AnnotationGroup(region, 'tenia omentalis', FMANumber = 'FMANumber unknown', lyphID = 'Lyph ID unknown')
    annotationGroups = [tlGroup, tmGroup, toGroup]

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

    mesh = fm.findMeshByDimension(3)

    if useCubicHermiteThroughWall:
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    else:
        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    result = elementtemplate.defineField(coordinates, -1, eft)

    # Tenia coli edge elements
    elementtemplate1 = mesh.createElementtemplate()
    elementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    eft1 = eftfactory.createEftWedgeXi1One()
    elementtemplate1.defineField(coordinates, -1, eft1)

    elementtemplate2 = mesh.createElementtemplate()
    elementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    eft2 = eftfactory.createEftWedgeXi1Zero()
    elementtemplate2.defineField(coordinates, -1, eft2)

    # create nodes
    elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*3
    TCEdgeFactor = 1.5
    zero = [0.0, 0.0, 0.0]
    prevNodeIdentifier = nodeIdentifier - 1
    prevElementIdentifier = elementIdentifier - 1

    xTCOuterList = []
    d1TCOuterList = []
    d2TCOuterList = []
    d3TCOuterList = []
    bniList = []
    bniTCList = []
    idxList = []
    wList = []

    set1 = list(range(int(elementsCountAroundTC/2)+1))
    set2 = list(range(set1[-1] + elementsCountAroundHaustrum, set1[-1] + elementsCountAroundHaustrum + elementsCountAroundTC + 1))
    set3 = list(range(set2[-1] + elementsCountAroundHaustrum, set2[-1] + elementsCountAroundHaustrum + elementsCountAroundTC +1))
    set4 = list(range(set3[-1] + elementsCountAroundHaustrum, set3[-1] + elementsCountAroundHaustrum + int(elementsCountAroundTC/2)))
    setTCIdx = set1 + set2 + set3 + set4

    for n2 in range(elementsCountAlong + 1):
        xTCRaw = []
        d1TCRaw = []
        d2TCRaw = []
        d3TCRaw = []
        for N in range(3):
            idxTCMid = elementsCountThroughWall*(elementsCountAlong+1)*elementsCountAround + n2*elementsCountAround + N*(elementsCountAroundTC + elementsCountAroundHaustrum)
            unitNorm = vector.normalise(d3List[idxTCMid])
            xMid = [xList[idxTCMid][i] + unitNorm[i]*tcThickness for i in range(3)]
            d1Mid = d1List[idxTCMid]
            TCStartIdx = idxTCMid - int(elementsCountAroundTC*0.5) if N > 0 else idxTCMid + 3*(elementsCountAroundTC + elementsCountAroundHaustrum) - int(elementsCountAroundTC*0.5)
            TCEndIdx = idxTCMid + int(elementsCountAroundTC*0.5)
            v1 = xList[TCStartIdx]
            v2 = xMid
            d1MidScaled = [c*tcWidth*TCEdgeFactor for c in vector.normalise(d1Mid)]
            v3 = xList[TCEndIdx]
            nx = [v1, v2, v3]
            nd1 = [d1List[TCStartIdx], d1MidScaled, d1List[TCEndIdx]]
            sx, sd1, se, sxi, _  = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAroundTC)
            xTCRaw = xTCRaw + sx[1:-1]
            if elementsCountAroundTC == 2:
                p = [v2[i] - v1[i] for i in range(3)]
                A = vector.dotproduct(unitNorm, p) # A<0 if v2 is higher than v1
                d1 = [c*tcWidth*0.5 for c in vector.normalise(d1Mid)] if A < 0 else d1MidScaled
            d1TCRaw = d1TCRaw + sd1[1:-1] if elementsCountAroundTC > 2 else d1TCRaw + [d1]
            xTCInnerSet = list(range(TCStartIdx+1, TCEndIdx)) if N > 0 else list(range(TCStartIdx + 1, TCStartIdx + int(elementsCountAroundTC * 0.5))) + list(range(idxTCMid, idxTCMid + int(elementsCountAroundTC * 0.5)))

            for n in range(elementsCountAroundTC - 1):
                xTCInner = xList[xTCInnerSet[n]]
                xTCOuter = sx[n + 1]
                d3 = [xTCOuter[i] - xTCInner[i] for i in range(3)]
                d3TCRaw.append(d3)
                node = nodes.findNodeByIdentifier(xTCInnerSet[n]+1)
                cache.setNode(node)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3)

                innerIdx = xTCInnerSet[n] - elementsCountThroughWall*(elementsCountAlong+1)*elementsCountAround
                d2 = d2List[xTCInnerSet[n]]
                factor = factorList[innerIdx]
                d2Unscaled = [ 1.0/factor*c for c in d2]
                curvature = curvatureAlong[innerIdx]
                distance = vector.magnitude([xTCOuter[i] - sxCentralLine[n2][i] for i in range(3)])
                newFactor = 1.0 + curvature*distance
                dx_ds2 = [ newFactor*c for c in d2Unscaled]
                d2TCRaw.append(dx_ds2)

        xTCOuterList = xTCOuterList + xTCRaw[int((elementsCountAroundTC-2)*0.5):] + xTCRaw[:int((elementsCountAroundTC-2)*0.5)]
        d1TCOuterList = d1TCOuterList + d1TCRaw[int((elementsCountAroundTC-2)*0.5):] + d1TCRaw[:int((elementsCountAroundTC-2)*0.5)]
        d2TCOuterList = d2TCOuterList + d2TCRaw[int((elementsCountAroundTC-2)*0.5):] + d2TCRaw[:int((elementsCountAroundTC-2)*0.5)]
        d3TCOuterList = d3TCOuterList + d3TCRaw[int((elementsCountAroundTC-2)*0.5):] + d3TCRaw[:int((elementsCountAroundTC-2)*0.5)]

    for n in range(len(xTCOuterList)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xTCOuterList[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1TCOuterList[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2TCOuterList[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3TCOuterList[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # create elements
    tlMeshGroup = tlGroup.getMeshGroup(mesh)
    tmMeshGroup = tmGroup.getMeshGroup(mesh)
    toMeshGroup = toGroup.getMeshGroup(mesh)

    for N in range(4):
        if N == 0:
            for e1 in range(int(elementsCountAroundTC*0.5)):
                bni = e1 + 1
                bniTC = e1 + 1
                bniList.append(bni)
                bniTCList.append(bniTC)
        elif N > 0 and N < 3:
            for e1 in range(elementsCountAroundTC):
                bni = e1 + int(elementsCountAroundTC*0.5) + (N-1)*(elementsCountAroundTC+1) + 2
                bniList.append(bni)
            bniTCList.append(bniTC + 1)
            for e1 in range(elementsCountAroundTC - 1):
                bniTC = bniTC + 1
                bniTCList.append(bniTC)
        else:
            for e1 in range(int(elementsCountAroundTC*0.5)):
                bni = e1 + int(elementsCountAroundTC*0.5) + (N-1)*(elementsCountAroundTC+1) + 2
                bniList.append(bni)
            if elementsCountAroundTC > 2:
                bniTCList.append(bniTC + 1)
                for e1 in range(int(elementsCountAroundTC*0.5 - 1)):
                    bniTC = bniTC + 1
                    bniTCList.append(bniTC)
            else:
                bniTCList.append(1)

    for e1 in range((elementsCountAroundTC+1)*3):
        idxTC = elementsCountThroughWall*(elementsCountAlong+1)*elementsCountAround + setTCIdx[e1] +1
        idxList.append(idxTC)

    for e2 in range(elementsCountAlong):
        e1 = -1
        A = e2*(elementsCountAroundTC-1)*3 + prevNodeIdentifier
        for n1 in range(int(elementsCountAroundTC*0.5)):
            e1 = e1 + 1
            regNodeIdentifiers = getNodeIdentifierForRegularElement(e1, e2, A, elementsCountAroundTC, elementsCountAround, bniList, bniTCList, idxList)
            nodeIdentifiers = regNodeIdentifiers if n1 < int(elementsCountAroundTC*0.5) - 1 else [regNodeIdentifiers[i] for i in range(4)] + [regNodeIdentifiers[4], regNodeIdentifiers[6]]
            element = mesh.createElement(elementIdentifier, elementtemplate if n1 < int(elementsCountAroundTC*0.5) - 1 else elementtemplate1)
            result = element.setNodesByIdentifier(eft if n1 < int(elementsCountAroundTC*0.5) - 1 else eft1, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1
            toMeshGroup.addElement(element)

        for N in range(2):
            for n1 in range(elementsCountAroundTC):
                e1 = e1 + 1
                regNodeIdentifiers = getNodeIdentifierForRegularElement(e1, e2, A, elementsCountAroundTC, elementsCountAround, bniList, bniTCList, idxList)
                if n1 == 0:
                    nodeIdentifiers = [regNodeIdentifiers[i] for i in range(4)] + [regNodeIdentifiers[4], regNodeIdentifiers[6]]
                    element = mesh.createElement(elementIdentifier, elementtemplate2)
                    result = element.setNodesByIdentifier(eft2, nodeIdentifiers)
                elif n1 > 0 and n1 < elementsCountAroundTC - 1:
                    nodeIdentifiers = regNodeIdentifiers
                    element = mesh.createElement(elementIdentifier, elementtemplate)
                    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                else:
                    nodeIdentifiers = [regNodeIdentifiers[i] for i in range(4)] + [regNodeIdentifiers[4], regNodeIdentifiers[6]]
                    element = mesh.createElement(elementIdentifier, elementtemplate1)
                    result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1
                tlMeshGroup.addElement(element) if N == 0 else tmMeshGroup.addElement(element)

        for n1 in range(int(elementsCountAroundTC*0.5)):
            e1 = e1 + 1
            regNodeIdentifiers = getNodeIdentifierForRegularElement(e1, e2, A, elementsCountAroundTC, elementsCountAround, bniList, bniTCList, idxList)
            nodeIdentifiers = regNodeIdentifiers if n1 > 0 else [regNodeIdentifiers[i] for i in range(4)] + [regNodeIdentifiers[4], regNodeIdentifiers[6]]
            element = mesh.createElement(elementIdentifier, elementtemplate if n1 > 0 else elementtemplate2)
            result = element.setNodesByIdentifier(eft if n1 > 0 else eft2, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1
            toMeshGroup.addElement(element)

    # Define texture coordinates field for tenia coli
    textureCoordinates = zinc_utils.getOrCreateTextureCoordinateField(fm)
    textureNodetemplate1 = nodes.createNodetemplate()
    textureNodetemplate1.defineField(textureCoordinates)
    textureNodetemplate1.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    textureNodetemplate1.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    textureNodetemplate1.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        textureNodetemplate1.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

    textureNodetemplate2 = nodes.createNodetemplate()
    textureNodetemplate2.defineField(textureCoordinates)
    textureNodetemplate2.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_VALUE, 2)
    textureNodetemplate2.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
    textureNodetemplate2.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D_DS2, 2)
    if useCrossDerivatives:
        textureNodetemplate2.setValueNumberOfVersions(textureCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 2)

    bicubichermitelinear = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eftTexture1 = bicubichermitelinear.createEftBasic()
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(textureCoordinates, -1, eftTexture1)

    eftTexture2 = bicubichermitelinear.createEftWedgeXi1One()
    elementtemplate1.defineField(textureCoordinates, -1, eftTexture2)

    eftTexture3 = bicubichermitelinear.createEftWedgeXi1Zero()
    elementtemplate2.defineField(textureCoordinates, -1, eftTexture3)

    eftTexture4 = bicubichermitelinear.createEftOpenTube()
    elementtemplate3 = mesh.createElementtemplate()
    elementtemplate3.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate3.defineField(textureCoordinates, -1, eftTexture4)

    eftTexture5 = bicubichermitelinear.createEftWedgeXi1ZeroOpenTube()
    elementtemplate4 = mesh.createElementtemplate()
    elementtemplate4.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate4.defineField(textureCoordinates, -1, eftTexture5)

    # Calculate texture coordinates and derivatives
    nodeIdentifier = prevNodeIdentifier + 1
    dTC = uList[1] - uList[0]

    for N in range(4):
        idxTCMid = N*(elementsCountAroundTC + elementsCountAroundHaustrum)
        TCStartIdx = idxTCMid - int(elementsCountAroundTC*0.5)
        TCEndIdx = idxTCMid + int(elementsCountAroundTC*0.5)
        v1 = [uList[TCStartIdx], 0.0, 1.0] if N > 0 else [-dTC, 0.0, 1.0]
        v2 = [  uList[idxTCMid], 0.0, 2.0]
        v3 = [  uList[TCEndIdx], 0.0, 1.0] if N < 3 else [ 1 + dTC , 0.0, 1.0]
        d1 = d2 = d3 = [dTC, 0.0, 0.0]
        nx = [v1, v2, v3]
        nd1 = [d1, d2, d3]
        sx, _, _, _, _ = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAroundTC)
        if N > 0 and N < 3:
            w = sx[1:-1]
        elif N == 0:
            w = sx[int(elementsCountAroundTC*0.5):-1]
        else: # N > 2
            w = sx[1:int(elementsCountAroundTC*0.5)+1]
        wList = wList + w

    d2 = [0.0, 1.0 / elementsCountAlong, 0.0]
    for n2 in range(elementsCountAlong + 1):
        for n1 in range((elementsCountAroundTC-1) * 3):
            u = [ wList[n1][0],
                  1.0 / elementsCountAlong * n2,
                  wList[n1][2]]
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            node.merge(textureNodetemplate2 if n1 == 0 else textureNodetemplate1)
            cache.setNode(node)
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, u)
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, [dTC, 0.0, 0.0])
            textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2)
            if useCrossDerivatives:
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            if n1 == 0:
                u = [ wList[-1][0],
                      1.0 / elementsCountAlong * n2,
                      wList[-1][2]]
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2, u)
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, [dTC, 0.0, 0.0])
                textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2, d2)
                if useCrossDerivatives:
                    textureCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
            nodeIdentifier = nodeIdentifier + 1

    # Define flat coordinates field for tenia coli
    flatCoordinates = zinc_utils.getOrCreateFlatCoordinateField(fm)
    flatNodetemplate1 = nodes.createNodetemplate()
    flatNodetemplate1.defineField(flatCoordinates)
    flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    if useCrossDerivatives:
        flatNodetemplate1.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

    flatNodetemplate2 = nodes.createNodetemplate()
    flatNodetemplate2.defineField(flatCoordinates)
    flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_VALUE, 2)
    flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS1, 2)
    flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D_DS2, 2)
    if useCrossDerivatives:
        flatNodetemplate2.setValueNumberOfVersions(flatCoordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 2)

    flatElementtemplate = mesh.createElementtemplate()
    flatElementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    flatElementtemplate.defineField(flatCoordinates, -1, eftTexture1)

    flatElementtemplate1 = mesh.createElementtemplate()
    flatElementtemplate1.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    flatElementtemplate1.defineField(flatCoordinates, -1, eftTexture2)

    flatElementtemplate2 = mesh.createElementtemplate()
    flatElementtemplate2.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    flatElementtemplate2.defineField(flatCoordinates, -1, eftTexture3)

    flatElementtemplate3 = mesh.createElementtemplate()
    flatElementtemplate3.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    flatElementtemplate3.defineField(flatCoordinates, -1, eftTexture4)

    flatElementtemplate4 = mesh.createElementtemplate()
    flatElementtemplate4.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    flatElementtemplate4.defineField(flatCoordinates, -1, eftTexture5)

    # Calculate flat coordinates and derivatives
    nodeIdentifier = prevNodeIdentifier + 1
    dTC = (uList[1] - uList[0])*arcLengthOuterMidLength
    wList = []
    dwList = []

    for N in range(4):
        idxTCMid = N*(elementsCountAroundTC + elementsCountAroundHaustrum)
        TCStartIdx = idxTCMid - int(elementsCountAroundTC*0.5)
        TCEndIdx = idxTCMid + int(elementsCountAroundTC*0.5)
        v1 = [uList[TCStartIdx]*arcLengthOuterMidLength, 0.0, wallThickness] if N > 0 else [-dTC, 0.0, wallThickness]
        v2 = [  uList[idxTCMid]*arcLengthOuterMidLength, 0.0, wallThickness + tcThickness]
        v3 = [  uList[TCEndIdx]*arcLengthOuterMidLength, 0.0, wallThickness] if N < 3 else [ arcLengthOuterMidLength + dTC , 0.0, wallThickness]
        d1 = d3 = [dTC, 0.0, 0.0]
        d2 = [c*10.0 for c in [dTC, 0.0, 0.0]]
        nx = [v1, v2, v3]
        nd1 = [d1, d2, d3]
        sx, sd1, _, _, _ = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAroundTC)
        if N > 0 and N < 3:
            w = sx[1:-1]
            dw = sd1[1:-1]
        elif N == 0:
            w = sx[int(elementsCountAroundTC*0.5):-1]
            dw = sd1[int(elementsCountAroundTC*0.5):-1]
        else: # N > 2
            w = sx[1:int(elementsCountAroundTC*0.5)+1]
            dw = sd1[1:int(elementsCountAroundTC*0.5)+1]
        wList = wList + w
        dwList = dwList + dw

    dxds2 = [0.0, totalLengthAlong / elementsCountAlong, 0.0]

    for n2 in range(elementsCountAlong + 1):
        for n1 in range((elementsCountAroundTC-1) * 3):
            x = [ wList[n1][0],
                  totalLengthAlong / elementsCountAlong * n2,
                  wList[n1][2]]
            node = nodes.findNodeByIdentifier(nodeIdentifier)
            node.merge(flatNodetemplate2 if n1 == 0 else flatNodetemplate1)
            cache.setNode(node)
            flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dwList[n1])
            flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dxds2)
            if useCrossDerivatives:
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            if n1 == 0:
                x = [ wList[-1][0],
                      totalLengthAlong / elementsCountAlong * n2,
                      wList[-1][2]]
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 2, x)
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 2, dwList[n1])
                flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 2, dxds2)
                if useCrossDerivatives:
                    flatCoordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 2, zero)
            nodeIdentifier = nodeIdentifier + 1

    # Define texture coordinates field & flat coordinates field over elements
    elementIdentifier = prevElementIdentifier + 1
    for e2 in range(elementsCountAlong):
        e1 = -1
        A = e2*(elementsCountAroundTC-1)*3 + prevNodeIdentifier

        for n1 in range(int(elementsCountAroundTC*0.5)):
            e1 = e1 + 1
            regNodeIdentifiers = getNodeIdentifierForRegularElement(e1, e2, A, elementsCountAroundTC, elementsCountAround, bniList, bniTCList, idxList)
            nodeIdentifiers = regNodeIdentifiers if n1 < int(elementsCountAroundTC*0.5) - 1 else [regNodeIdentifiers[i] for i in range(4)] + [regNodeIdentifiers[4], regNodeIdentifiers[6]]
            element = mesh.findElementByIdentifier(elementIdentifier)
            element.merge(elementtemplate if n1 < int(elementsCountAroundTC*0.5) - 1 else elementtemplate1)
            element.merge(flatElementtemplate if n1 < int(elementsCountAroundTC*0.5) - 1 else flatElementtemplate1)
            element.setNodesByIdentifier(eftTexture1 if n1 < int(elementsCountAroundTC*0.5) - 1 else eftTexture2, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

        for N in range(2):
            for n1 in range(elementsCountAroundTC):
                e1 = e1 + 1
                regNodeIdentifiers = getNodeIdentifierForRegularElement(e1, e2, A, elementsCountAroundTC, elementsCountAround, bniList, bniTCList, idxList)
                element = mesh.findElementByIdentifier(elementIdentifier)
                if n1 == 0:
                    nodeIdentifiers = [regNodeIdentifiers[i] for i in range(4)] + [regNodeIdentifiers[4], regNodeIdentifiers[6]]
                    element.merge(elementtemplate2)
                    element.merge(flatElementtemplate2)
                    element.setNodesByIdentifier(eftTexture3, nodeIdentifiers)
                elif n1 > 0 and n1 < elementsCountAroundTC - 1:
                    nodeIdentifiers = regNodeIdentifiers
                    element.merge(elementtemplate)
                    element.merge(flatElementtemplate)
                    element.setNodesByIdentifier(eftTexture1, nodeIdentifiers)
                else:
                    nodeIdentifiers = [regNodeIdentifiers[i] for i in range(4)] + [regNodeIdentifiers[4], regNodeIdentifiers[6]]
                    element.merge(elementtemplate1)
                    element.merge(flatElementtemplate1)
                    element.setNodesByIdentifier(eftTexture2, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

        for n1 in range(int(elementsCountAroundTC*0.5)):
            onOpening = (n1 == int(elementsCountAroundTC*0.5 - 1))
            e1 = e1 + 1
            regNodeIdentifiers = getNodeIdentifierForRegularElement(e1, e2, A, elementsCountAroundTC, elementsCountAround, bniList, bniTCList, idxList)
            nodeIdentifiers = regNodeIdentifiers if n1 > 0 else [regNodeIdentifiers[i] for i in range(4)] + [regNodeIdentifiers[4], regNodeIdentifiers[6]]
            element = mesh.findElementByIdentifier(elementIdentifier)
            if n1 > 0:
                element.merge(elementtemplate3 if onOpening else elementtemplate)
                element.merge(flatElementtemplate3 if onOpening else flatElementtemplate)
                element.setNodesByIdentifier(eftTexture4 if onOpening else eftTexture1, nodeIdentifiers)
            else:
                element.merge(elementtemplate4 if onOpening else elementtemplate2)
                element.merge(flatElementtemplate4 if onOpening else flatElementtemplate2)
                result = element.setNodesByIdentifier(eftTexture5 if onOpening else eftTexture3, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

    fm.endChange()

    return annotationGroups, nodeIdentifier, elementIdentifier

def getNodeIdentifierForRegularElement(e1, e2, A, elementsCountAroundTC, elementsCountAround, bniList, bniTCList, idxList):
    """
    Get node identifiers to create regular elements for tenia coli.
    :param e1: Element count iterator around tenia coli
    :param e2: Element count iterator along tenia coli
    :param A: Element offset
    :param elementsCountAroundTC: Number of elements around tenia coli
    :param elementsCountAround: Number of elements around scaffold
    :param bniList: Base node indices for nodes lying on the outer
    surface of haustra over the boundary of tenia coli
    :param bniTCList: Base node index for tenia coli nodes
    :param idxList: List of global numbering of nodes lying on the outer
    surface of haustra over the boundary of tenia coli
    :return: nodeIdentifiers
    """
    bni111 = idxList[bniList[e1]-1] + e2*elementsCountAround
    bni121 = idxList[bniList[e1]%((elementsCountAroundTC+1)*3)] + e2*elementsCountAround
    bni211 = bni111 + elementsCountAround
    bni221 = bni121 + elementsCountAround
    bni112 = bniTCList[e1] + A
    bni122 = (bniTCList[e1]+1)%((elementsCountAroundTC - 1)*3) + A if (bniTCList[e1]+1)%((elementsCountAroundTC - 1)*3)>0 else bniTCList[e1]+1 + A
    bni212 = bni112 + (elementsCountAroundTC-1)*3
    bni222 = bni122 + (elementsCountAroundTC-1)*3
    nodeIdentifiers = [ bni111, bni121, bni211, bni221, bni112, bni122, bni212, bni222]

    return nodeIdentifiers
