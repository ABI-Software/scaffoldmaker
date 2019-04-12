"""
Generates a single 3-D haustra segment mesh along a central
line, with variable numbers of elements around, along and
through wall, with variable radius and thickness along.
"""

import math
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils import vector

class MeshType_3d_haustra1(Scaffold_base):
    '''
    Generates a single 3-D haustra segment mesh with variable
    numbers of elements around, along the central line, and
    through wall. The haustra segment has a triangular profile
    with rounded corners at the inter-haustral septa, and a
    clover profile in the intra-haustral region.
    '''
    @staticmethod
    def getName():
        return '3D Haustra 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of elements around tenia coli' : 2,
            'Number of elements around haustrum' : 4,
            'Number of elements along haustrum' : 3,
            'Number of elements through wall' : 1,
            'Inner width of tenia coli': 0.2,
            'Inner radius': 0.5,
            'Corner inner radius factor': 0.5,
            'Haustrum inner radius factor': 0.5,
            'Haustrum length end derivative factor': 0.5,
            'Haustrum length mid derivative factor': 1.0,
            'Wall thickness': 0.01,
            'Haustrum length': 1.0,
            'Use cross derivatives' : False,
            'Use linear through wall' : True,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along haustrum' : 1,
            'Refine number of elements through wall' : 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around tenia coli',
            'Number of elements around haustrum',
            'Number of elements along haustrum',
            'Number of elements through wall',
            'Inner width of tenia coli',
            'Inner radius',
            'Corner inner radius factor',
            'Haustrum inner radius factor',
            'Haustrum length end derivative factor',
            'Haustrum length mid derivative factor',
            'Wall thickness',
            'Haustrum length',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along haustrum',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along haustrum',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        for key in [
            'Number of elements around tenia coli',
            'Number of elements along haustrum']:
            if options[key] < 2:
                options[key] = 2
        if options['Number of elements around haustrum'] < 4:
            options['Number of elements around haustrum'] = 4
        for key in [
            'Number of elements around tenia coli',
            'Number of elements around haustrum']:
            if options[key] % 2 > 0:
                options[key] = options[key] + 1
        if options['Inner width of tenia coli'] < 0.1*options['Inner radius']:
            options['Inner width of tenia coli'] = 0.1*options['Inner radius']
        if options['Inner width of tenia coli'] > round(math.sqrt(3)*0.5*options['Inner radius'],2):
            options['Inner width of tenia coli'] = round(math.sqrt(3)*0.5*options['Inner radius'],2)
        for key in [
            'Inner radius',
            'Haustrum inner radius factor',
            'Haustrum length end derivative factor',
            'Haustrum length mid derivative factor',
            'Wall thickness',
            'Haustrum length']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Corner inner radius factor'] < 0.1:
            options['Corner inner radius factor'] = 0.1
        for key in [
            'Corner inner radius factor',
            'Haustrum length end derivative factor']:
            if options[key] > 1.0:
                options[key] = 1.0

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
        elementsCountAlongHaustrum = options['Number of elements along haustrum']
        elementsCountThroughWall = options['Number of elements through wall']
        widthTC = options['Inner width of tenia coli']
        radius = options['Inner radius']
        cornerInnerRadiusFactor = options['Corner inner radius factor']
        haustrumInnerRadiusFactor = options['Haustrum inner radius factor']
        haustrumLengthEndDerivativeFactor = options['Haustrum length end derivative factor']
        haustrumLengthMidDerivativeFactor = options['Haustrum length mid derivative factor']
        wallThickness = options['Wall thickness']
        haustrumLength = options['Haustrum length']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = elementsCountAlongHaustrum
        haustraSegmentCount = 1

        cx = [ [ 0.0, 0.0, 0.0 ], [ haustrumLength, 0.0, 0.0 ] ]
        cd1 = [ [ haustrumLength, 0.0, 0.0 ], [ haustrumLength, 0.0, 0.0 ] ]

        # Generate inner surface of a haustra segment
        xHaustraInner, d1HaustraInner, d2HaustraInner, haustraSegmentAxis = getColonHaustraSegmentInnerPoints(elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongHaustrum, widthTC, radius, cornerInnerRadiusFactor,
            haustrumInnerRadiusFactor, haustrumLengthEndDerivativeFactor, haustrumLengthMidDerivativeFactor, haustrumLength)

        # Generate tube mesh
        annotationGroups, nextNodeIdentifier, nextElementIdentifier = tubemesh.generatetubemesh(region, elementsCountAround, elementsCountAlongHaustrum, elementsCountThroughWall, haustraSegmentCount,
            cx, cd1, xHaustraInner, d1HaustraInner, d2HaustraInner, wallThickness, haustraSegmentAxis, haustrumLength, useCrossDerivatives, useCubicHermiteThroughWall)

        return annotationGroups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)
        return meshrefinement.getAnnotationGroups()

def getColonHaustraSegmentInnerPoints(elementsCountAroundTC, elementsCountAroundHaustrum, elementsCountAlongHaustrum, widthTC, radius, cornerInnerRadiusFactor,
    haustrumInnerRadiusFactor, haustrumLengthEndDerivativeFactor, haustrumLengthMidDerivativeFactor, haustrumLength):
    """
    Generates a 3-D haustra segment mesh with variable numbers
    of elements around, along the central line, and through wall.
    The haustra segment has a triangular profile with rounded corners
    at the inter-haustral septa, and a clover profile in the intra-haustral
    region.
    :param elementsCountAroundTC: Number of elements around each tenia coli.
    :param elementsCountAroundHaustrum: Number of elements around haustrum.
    :param elementsCountAlongHaustrum: Number of elements along haustrum.
    :param widthTC: Width of tenia coli.
    :param radius: Inner radius defined from center of triangular
    profile to vertex of the triangle.
    :param cornerInnerRadiusFactor: Roundness of triangular corners of
    interhaustral septa. Factor is multiplied by inner radius
    to get a radius of curvature at the corners.
    :param haustrumInnerRadiusFactor: Factor is multiplied by inner
    radius to obtain radius of intersecting circles in the middle cross-section
    along a haustra segment.
    :param haustrumLengthEndDerivativeFactor: Factor is multiplied by haustrum
    length to scale derivative along the end of a haustrum length.
    :param haustrumLengthMidDerivativeFactor: Factor is multiplied by haustrum
    length to scale derivative along the mid length of the haustrum.
    :param haustrumLength: Length of a haustrum.
    :return: coordinates, derivatives on inner surface of haustra segment.
    """

    # create nodes
    zero = [0.0, 0.0, 0.0] # Delete later if not in use
    x = [ 0.0, 0.0, 0.0 ]
    dx_ds1 = [ 0.0, 0.0, 0.0 ]
    dx_ds2 = [ 0.0, 0.0, 0.0 ]
    dx_ds3 = [ 0.0, 0.0, 0.0 ]
    cornerRC = cornerInnerRadiusFactor*radius
    unitZ = [0.0, 0.0, 1.0]
    haustrumRadius = (haustrumInnerRadiusFactor + 1)*radius

    xAround = []
    dx_ds1InnerList = []
    xTC = []
    d1TC = []
    xHalfSet = []
    d1HalfSet = []
    nxHaustrum = []
    nd1Haustrum = []
    xInnerRaw = []
    dx_ds2InnerRaw = []
    xInnerList = []
    dx_ds2InnerList = []

    # Inter-haustral segment
    # Calculate boundary of tenia coli and haustrum
    xc = [(radius - cornerRC)* math.cos(0.0), (radius - cornerRC)*math.sin(0.0), 0.0]
    pt1 = [xc[0] + cornerRC*math.cos(0.0), xc[1] + cornerRC*math.sin(0.0), 0.0]
    pt2 = [xc[0] + cornerRC*math.cos(math.pi/4), xc[1] + cornerRC*math.sin(math.pi/4), 0.0]
    pt3 = [(radius - cornerRC) * math.cos(2*math.pi/3) + cornerRC*math.cos(2*math.pi/3 + math.pi*7/4),
           (radius - cornerRC) * math.sin(2*math.pi/3) + cornerRC*math.sin(2*math.pi/3 + math.pi*7/4), 
           0.0]
    xMid = [(pt2[n] + pt3[n])*0.5 for n in range(3)]

    if widthTC/2 < cornerRC*math.sin(math.pi/4):
        TCTheta = math.asin(widthTC/(2*cornerRC))
        thetaSet = [-TCTheta, TCTheta]
        for n in range(2):
            theta = thetaSet[n]
            cosTheta = math.cos(theta)
            sinTheta = math.sin(theta)
            x = [xc[0] + cornerRC*cosTheta, xc[1] + cornerRC*sinTheta, 0.0]
            d1 = [-cornerRC*sinTheta*TCTheta*2.0, cornerRC*cosTheta*TCTheta*2.0, 0.0]
            xTC.append(x)
            d1TC.append(d1)
        thetaDiff = math.pi/4 - TCTheta
        d1Pt2 = [-cornerRC*math.sin(math.pi/4)*thetaDiff, cornerRC*math.cos(math.pi/4)*thetaDiff, 0.0]
        d1HaustrumEnd = [xMid[c] - pt2[c] for c in range(3)]
    else:
        widthOutsideCircle = widthTC/2 - cornerRC*math.sin(math.pi/4)
        m = (pt2[1] - pt3[1])/(pt2[0] - pt3[0])
        c = pt2[1] - m*pt2[0]
        xA = [(widthOutsideCircle + pt2[1] - c)/m, widthOutsideCircle + pt2[1], 0.0]
        d1A = [xA[c] - pt2[c] for c in range(3)]
        xB = [xA[0], -xA[1], 0.0]
        reflectPt2 = [pt2[0], -pt2[1], 0.0]
        d1B = [reflectPt2[c] - xB[c] for c in range(3)]
        d1Pt1 = [-cornerRC*math.sin(0.0), cornerRC*math.cos(0.0), 0.0]
        xTC = [xB, pt1, xA]
        d1TC = [d1B, d1Pt1, d1A]
        d1HaustrumEnd = [xMid[c] - xA[c] for c in range(3)]
    sxTC, sd1TC, se, sxi, _= interp.sampleCubicHermiteCurves(xTC, d1TC, elementsCountAroundTC)

    xHaustrumStart = sxTC[-1]
    d1HaustrumStart = sd1TC[-1]
    cosRotAng = math.cos(2*math.pi/3)
    sinRotAng = math.sin(2*math.pi/3)
    xHaustrumEnd = xMid
    nx = [xHaustrumStart, xHaustrumEnd] if widthTC/2 > cornerRC*math.sin(math.pi/4) else [xHaustrumStart, pt2, xHaustrumEnd]
    nd1 = [d1HaustrumStart, d1HaustrumEnd] if widthTC/2 > cornerRC*math.sin(math.pi/4) else [d1HaustrumStart, d1Pt2, d1HaustrumEnd]
    sxHaustrum, sd1Haustrum, se, sxi, _= interp.sampleCubicHermiteCurves(nx, nd1, int(elementsCountAroundHaustrum/2),
                                         addLengthStart = 0.5*vector.magnitude(nd1[0]), lengthFractionStart = 0.5)

    xHalfSet = xHalfSet + sxTC[int(elementsCountAroundTC/2):-1] + sxHaustrum
    d1HalfSet = d1HalfSet + sd1TC[int(elementsCountAroundTC/2):-1] + sd1Haustrum
    xInner, d1Inner = getFullProfileFromHalfHaustrum(xHalfSet, d1HalfSet)

    # Haustra segment
    # Re-initialise for haustra segment
    xTC = []
    d1TC = []
    xHalfSet = []
    d1HalfSet = []

    xTC2 = radius* math.cos(2*math.pi/3)
    yTC2 = radius* math.sin(2*math.pi/3)
    originRC = (xTC2*xTC2 + yTC2*yTC2 - haustrumRadius*haustrumRadius) / (2*(-xTC2 - haustrumRadius))
    RC = haustrumRadius - originRC

    # Rotate to find originRC of 1st haustrum
    yTC1 = pt1[1]
    rotOriginRC = [ originRC*math.cos(-2/3*math.pi), originRC*math.sin(-2/3*math.pi), 0.0]

    # Teniae coli boundary on 1st haustrum
    thetaStart = math.asin((yTC1 + rotOriginRC[1]) / RC)
    thetaEnd = math.pi - math.asin((yTC2 + rotOriginRC[1])/ RC)
    thetaHalfHaustrum = (thetaEnd - thetaStart)/2 + thetaStart
    thetaStartToHalfTC = math.asin(widthTC/(2*RC)+ math.sin(thetaStart))
    thetaSet1 = [thetaStart, thetaStartToHalfTC]
    thetaDiff = (thetaHalfHaustrum + thetaStartToHalfTC)*0.5
    thetaSet2 = [thetaDiff, thetaHalfHaustrum]

    nxTC, nd1TC = getCircleXandD1FromRadians(thetaSet1, RC, rotOriginRC)
    xHaustrumPt, d1HaustrumPt = getCircleXandD1FromRadians(thetaSet2, RC, rotOriginRC)

    if int(elementsCountAroundTC/2) > 1:
        sxTC, sd1TC, _, _ , _ = interp.sampleCubicHermiteCurves(nxTC, nd1TC, int(elementsCountAroundTC/2))
        xTC = sxTC[:-1]
        d1TC = sd1TC[:-1]
        xHaustrumStart = sxTC[-1]
        d1HaustrumStart = sd1TC[-1]
    else:
        xTC = nxTC[:-1]
        d1TC = nd1TC[:-1]
        xHaustrumStart = nxTC[-1]
        d1HaustrumStart = nd1TC[-1]

    nxHaustrum.append(xHaustrumStart)
    nxHaustrum = nxHaustrum + xHaustrumPt
    nd1Haustrum.append(d1HaustrumStart)
    nd1Haustrum = nd1Haustrum + d1HaustrumPt
    sxHaustrum, sd1Haustrum, _, _ , _ = interp.sampleCubicHermiteCurves(nxHaustrum, nd1Haustrum, int(elementsCountAroundHaustrum/2),
                                                                        addLengthStart = 0.5*vector.magnitude(d1HaustrumStart), lengthFractionStart = 0.5)
    xHalfSet = xHalfSet + xTC + sxHaustrum
    d1HalfSet = d1HalfSet + d1TC + sd1Haustrum
    xHaustraInner, d1HaustraInner = getFullProfileFromHalfHaustrum(xHalfSet, d1HalfSet)

    # Sample arclength of haustra segment
    elementsCountAround = (elementsCountAroundTC + elementsCountAroundHaustrum)*3

    for n1 in range(elementsCountAround):
        if n1%(elementsCountAround/3) > 0.0:
            v1 = [xInner[n1][0], xInner[n1][1], 0.0]
            startArcLength = haustrumLengthEndDerivativeFactor * haustrumLength
            d1 = [ c*startArcLength for c in unitZ]
            v2 = [xHaustraInner[n1][0], xHaustraInner[n1][1], haustrumLength/2]
            midArcLength = haustrumLengthMidDerivativeFactor * haustrumLength
            d2 = [ c*midArcLength for c in unitZ]
            v3 = [xInner[n1][0], xInner[n1][1], haustrumLength]
            d3 = [ c*startArcLength for c in unitZ]
        else:
            v1 = [xInner[n1][0], xInner[n1][1], 0.0]
            v2 = [xInner[n1][0], xInner[n1][1], haustrumLength/2]
            v3 = [xInner[n1][0], xInner[n1][1], haustrumLength]
            d3 = d2 = d1 = [c* haustrumLength/3 for c in unitZ]
        nx = [v1, v2, v3]
        nd1 = [d1, d2, d3]
        sx, sd1, se, sxi, _  = interp.sampleCubicHermiteCurves(nx, nd1, elementsCountAlongHaustrum)
        xInnerRaw.append(sx)
        dx_ds2InnerRaw.append(sd1)

    # Re-arrange sample order & calculate dx_ds1 and dx_ds3 from dx_ds2
    dx_ds1InnerList = dx_ds1InnerList + d1Inner
    for n2 in range(elementsCountAlongHaustrum + 1):
        xAround = []
        unitdx_ds1Around = []
        for n1 in range(elementsCountAround):
            x = xInnerRaw[n1][n2]
            xInnerList.append(x)
            dx_ds2 = dx_ds2InnerRaw[n1][n2]
            dx_ds2InnerList.append(dx_ds2)
            unitTangent = vector.normalise(dx_ds2)
            # Inter-haustra
            if n2 == 0 or n2 > elementsCountAlongHaustrum - 1:
                dx_ds1 = d1Inner[n1]
                unitdx_ds3 = vector.crossproduct3(vector.normalise(dx_ds1), unitTangent)
            else:
                # Intra-Haustra
                if elementsCountAlongHaustrum == 2:
                    unitdx_ds1 = vector.normalise(d1HaustraInner[n1])
                else:
                    if n1%(elementsCountAround/3) == 0: # intersection points
                        unitdx_ds1 = vector.normalise(d1HaustraInner[n1])
                    else: # points on clover
                        if elementsCountAlongHaustrum > 3:
                            if n2 < int(elementsCountAlongHaustrum/2): # first half of haustrumLength
                                axisRot = vector.crossproduct3(unitZ, unitTangent)
                            elif n2 > int(elementsCountAlongHaustrum/2): # second half of haustrumLength
                                axisRot = vector.crossproduct3(unitTangent, unitZ)
                        elif elementsCountAlongHaustrum == 3: # 3 elementsAlongHaustrum
                            axisRot = vector.crossproduct3(unitTangent, unitZ)

                        rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, math.pi/2)
                        rotNormal = [rotFrame[j][0]*unitTangent[0] + rotFrame[j][1]*unitTangent[1] + rotFrame[j][2]*unitTangent[2] for j in range(3)]
                        unitdx_ds3 = vector.normalise(rotNormal)
                        unitdx_ds1 = vector.crossproduct3(unitTangent, unitdx_ds3)
                xAround.append(x)
                unitdx_ds1Around.append(unitdx_ds1)

        if n2 > 0 and n2 < elementsCountAlongHaustrum:
            dx_ds1InnerAroundList = []
            for n1 in range(elementsCountAround):
                v1 = xAround[n1]
                d1 = unitdx_ds1Around[n1]
                v2 = xAround[(n1+1)%elementsCountAround]
                d2 = unitdx_ds1Around[(n1+1)%elementsCountAround]
                arcLengthAround = interp.computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c*arcLengthAround for c in d1]
                dx_ds1InnerAroundList.append(dx_ds1)
            d1Smoothed = interp.smoothCubicHermiteDerivativesLoop(xAround, dx_ds1InnerAroundList)
            dx_ds1InnerList = dx_ds1InnerList + d1Smoothed

    dx_ds1InnerList = dx_ds1InnerList + d1Inner

    return xInnerList, dx_ds1InnerList, dx_ds2InnerList, unitZ

def getFullProfileFromHalfHaustrum(xHaustrumHalfSet, d1HaustrumHalfSet):
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
    :param d1HaustrumHalfSet: Derivatives of points in first
    half of the first sector.
    :return: coordinates, derivatives of points over entire profile.
    """
    xHaustrumHalfSet2 = []
    d1HaustrumHalfSet2 = []
    xHaustrum = []
    d1Haustrum = []
    xHaustra = []
    d1Haustra = []

    for n in range(1,len(xHaustrumHalfSet)):
        idx =  -n + len(xHaustrumHalfSet) - 1
        x = xHaustrumHalfSet[idx]
        d1 = d1HaustrumHalfSet[idx]
        xReflect = [x[0], -x[1], 0.0]
        d1Reflect = [d1[0], -d1[1], 0.0]
        xRot = [xReflect[0]*math.cos(2/3*math.pi) - xReflect[1]*math.sin(2/3*math.pi),
                xReflect[0]*math.sin(2/3*math.pi) + xReflect[1]*math.cos(2/3*math.pi),
                0.0]
        d1Rot = [-(d1Reflect[0]*math.cos(2/3*math.pi) - d1Reflect[1]*math.sin(2/3*math.pi)),
                -(d1Reflect[0]*math.sin(2/3*math.pi) + d1Reflect[1]*math.cos(2/3*math.pi)),
                0.0]
        xHaustrumHalfSet2.append(xRot)
        d1HaustrumHalfSet2.append(d1Rot)

    xHaustrum = xHaustrumHalfSet + xHaustrumHalfSet2
    d1Haustrum = d1HaustrumHalfSet + d1HaustrumHalfSet2

    # Rotate to get all 3 sectors
    xHaustra = xHaustra + xHaustrum[:-1]
    d1Haustra = d1Haustra + d1Haustrum[:-1]
    ang = [ 2/3*math.pi, -2/3*math.pi]
    for i in range(2):
        rotAng = ang[i]
        cosRotAng = math.cos(rotAng)
        sinRotAng = math.sin(rotAng)
        for n in range(len(xHaustrum)- 1):
            x = xHaustrum[n]
            d1 = d1Haustrum[n]
            x = [ x[0]*cosRotAng - x[1]*sinRotAng, x[0]*sinRotAng + x[1]*cosRotAng, 0.0]
            xHaustra.append(x)
            dx_ds1 = [ d1[0]*cosRotAng - d1[1]*sinRotAng, d1[0]*sinRotAng + d1[1]*cosRotAng, 0.0]
            d1Haustra.append(dx_ds1)

    return xHaustra, d1Haustra

def getCircleXandD1FromRadians(thetaSet, radius, origin):
    """
    Gets the coordinates and derivatives along a circular path
    based on the angular range.
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
