"""
Generates a 3-D colon mesh along the central line, with variable
numbers of elements around, along and through wall, with
variable radius and thickness along.
"""
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tubemesh2 import *

class MeshType_3d_colon1(Scaffold_base):
    '''
    Generates a 3-D colon mesh with variable numbers
    of elements around, along the central line, and through wall.
    The colon is created by a function that generates haustra units
    and uses tubemesh to map the units along a central line profile.
    '''
    @staticmethod
    def getName():
        return '3D Colon 1'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Section 1']

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = {
            'Number of elements around': 15,
            'Number of elements along haustrum': 4,
            'Number of elements through wall': 1,
            'Number of haustra': 30,
            'Inner radius': 1.0,
            'Corner radius fraction': 0.5,
            'Haustra radius fraction': 0.5,
            'Interhaustra fold factor': 0.5,
            'Haustra curvature factor': 2.0,
            'Wall thickness': 0.05,
            'Tube type': 2,
            'Use cross derivatives': False,
            'Use linear through wall': False,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements along': 1,
            'Refine number of elements through wall': 1
        }
        if 'Human 1' in parameterSetName:
            options['Tube type'] = 3
        elif 'Section 1' in parameterSetName:
            options['Tube type'] = 1
            options['Number of haustra'] = 3
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around',
            'Number of elements along haustrum',
            'Number of elements through wall',
            'Number of haustra',
            'Inner radius',
            'Corner radius fraction',
            'Haustra radius fraction',
            'Interhaustra fold factor',
            'Haustra curvature factor',
            'Wall thickness',
            'Tube type',
            'Use cross derivatives',
            'Use linear through wall',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall'
        ]

    def checkOptions(options):
        for key in [
            'Number of elements through wall',
            'Number of haustra',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if (options['Number of elements around'] < 9) :
            options['Number of elements around'] = 9
        if (options['Number of elements around'] % 3 != 0) :
            options['Number of elements around'] = (options['Number of elements around']//3)*3
        for key in [
            'Inner radius',
            'Haustra radius fraction',
            'Interhaustra fold factor',
            'Haustra curvature factor',
            'Wall thickness']:
            if options[key] < 0.0:
                options[key] = 0.0
        if options['Corner radius fraction'] < 0.1:
            options['Corner radius fraction'] = 0.1
        for key in [
            'Corner radius fraction',
            'Interhaustra fold factor']:
            if options[key] > 1.0:
                options[key] = 1.0
        if (options['Tube type'] < 1 or options['Tube type'] > 3 ) :
            options['Tube type'] = 1

    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountAround = options['Number of elements around']
        elementsCountAlongHaustrum = options['Number of elements along haustrum']
        elementsCountThroughWall = options['Number of elements through wall']
        haustraCount = options['Number of haustra']
        radius = options['Inner radius']
        cornerRadiusFraction = options['Corner radius fraction']
        haustraRadiusFraction = options['Haustra radius fraction']
        foldFactor = options['Interhaustra fold factor']
        haustraCurvatureFactor = options['Haustra curvature factor']
        wallThickness = options['Wall thickness']
        tubeType = options['Tube type']
        useCrossDerivatives = options['Use cross derivatives']
        useCubicHermiteThroughWall = not(options['Use linear through wall'])
        elementsCountAlong = int(elementsCountAlongHaustrum*haustraCount)

        if tubeType == 1: # Straight tube
            cx = [[-4.0, 1.0, 3.0], [ 1.0, 2.0, 0.0 ] ]
            cd1 = [[ 5.0, 1.0, -3.0 ], [ 5.0, 1.0, -3.0 ]]
        elif tubeType == 2: # Human colon in x-y plane
            cx = [ [ 0.0, 0.0, 0.0], [0.0, 10.0, 0.0], [5.0, 9.0, 0.0], [ 10.0, 10.0, 0.0 ], [ 10.0, -2.0, 0.0], [ 7.0, -4.0, 0.0] ]
            cd1 = [ [ 0.0, 10.0, 0.0 ], [ 5.0, 5.0, 0.0 ], [5.0, 0.0, 0.0], [ 5.0, -5.0, 0.0 ], [ -3.0, -5.0, 0.0 ], [ -3.0, 0.0, 0.0 ]]
        elif tubeType == 3: # Human colon in 3D
            cx = [ [ 0.0, 0.0, 0.0], [0.0, 10.0, 3.0], [5.0, 9.0, 0.0], [ 10.0, 10.0, 2.0 ], [15.0, 15.0, 7.0], [ 20.0, -2.0, 0.0], [ 10.0, -4.0, -0.0] ]
            cd1 = [ [ 0.0, 10.0, 3.0 ], [ 5.0, 5.0, 0.0 ], [5.0, 0.0, 0.0], [ 10.0, -5.0, 0.0 ], [12.0, 12.0, 0.0], [ 5.0, -12.0, -5.0 ], [ -8.0, 0.0, 0.0 ]]

        # find arclength of colon
        length = 0.0
        elementsCountIn = len(cx) - 1
        for e in range(elementsCountIn):
            arcLength = computeCubicHermiteArcLength(cx[e], cd1[e], cx[e + 1], cd1[e + 1], rescaleDerivatives = True)
            length += arcLength
        haustraLength = length / haustraCount

        # Generate outer surface of a haustra unit
        xHaustraOuter, d1HaustraOuter, d2HaustraOuter, haustraAxis = getUnitHaustraOuter(elementsCountAround, elementsCountAlongHaustrum, radius, cornerRadiusFraction,
            haustraRadiusFraction, foldFactor, haustraCurvatureFactor, wallThickness, haustraLength)

        # Generate tube mesh
        nextNodeIdentifier, nextElementIdentifier = generatetubemesh(region, elementsCountAround, elementsCountAlongHaustrum, elementsCountThroughWall, haustraCount,
            cx, cd1, xHaustraOuter, d1HaustraOuter, d2HaustraOuter, wallThickness, haustraAxis, haustraLength, useCrossDerivatives, useCubicHermiteThroughWall)

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
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong, refineElementsCountThroughWall)

def getUnitHaustraOuter(elementsCountAround, elementsCountAlongHaustrum, radius, cornerRadiusFraction,
        haustraRadiusFraction, foldFactor, haustraCurvatureFactor, wallThickness, haustraLength):
    """
    Generates a 3-D haustra unit mesh with variable numbers
    of elements around, along the central line, and through wall.
    The haustra has a triangular profile with rounded corners at the
    interhaustral septa, and a clover profile outside the septa.
    :param elementsCountAround: Number of elements around haustra.
    :param elementsCountAlongHaustrum: Number of elements along haustrum.
    :param radius: Inner radius of haustra.
    :param cornerRadiusFraction: Roundness of triangular corners of
    interhaustral septa. Factor is multiplied by inner haustra radius
    to get a radius of curvature at the corners.
    :param haustraRadiusFraction: Factor is multiplied by inner haustra
    radius to obtain radius of intersecting circles in the middle cross-section
    along a haustra.
    :param foldFactor: Factor is multiplied by haustra length to scale derivative
    along the haustra at the interhaustral septa end.
    :param haustraCurvatureFactor: Factor is multiplied by haustra length to 
    scale derivative along the haustra at the middle section of the haustra.
    :param wallThickness: Thickness of haustra through wall.
    :param haustraLength: Length of a haustra unit.
    :return: coordinates, derivatives on outer surface of haustra unit.
    """

    # create nodes
    x = [ 0.0, 0.0, 0.0 ]
    dx_ds1 = [ 0.0, 0.0, 0.0 ]
    dx_ds2 = [ 0.0, 0.0, 0.0 ]
    dx_ds3 = [ 0.0, 0.0, 0.0 ]
    radiansRangeRC = [7*math.pi/4, 0.0, math.pi/4]
    cornerRC = cornerRadiusFraction*radius
    unitZ = [0.0, 0.0, 1.0]
    haustraRadius = (haustraRadiusFraction + 1)*radius
    xAround = []
    d1Around = []
    xInner = []
    dx_ds1InnerList = []
    xHaustraSide = []
    xHaustraInner = []
    d1InnerHaustraRaw = []
    xInnerRaw = []
    dx_ds2InnerRaw = []
    xInnerList = []
    dx_ds2InnerList = []
    dx_ds3InnerUnitList = []

    # Pre-calculate node locations and derivatives on inner triangle
    for n1 in range(3):
        radiansAround = n1*2*math.pi / 3
        cosRadiansAround = math.cos(radiansAround)
        sinRadiansAround = math.sin(radiansAround)
        xc = [(radius - cornerRC) * cosRadiansAround, (radius - cornerRC) * sinRadiansAround, 0.0]

        for n in range(3):
            radiansRC = radiansAround + radiansRangeRC[n]
            cosRadiansRC = math.cos(radiansRC)
            sinRadiansRC = math.sin(radiansRC)
            x = [xc[0] + cornerRC*cosRadiansRC, xc[1] + cornerRC*sinRadiansRC, 0.0]
            xAround.append(x)
            d1 = [ cornerRC*math.pi/4 * -sinRadiansRC, cornerRC*math.pi/4 * cosRadiansRC, 0.0]
            d1Around.append(d1)

    xSample = xAround[1:9]
    xSample.append(xAround[0])
    xSample.append(xAround[1])
    d1Sample = d1Around[1:9]
    d1Sample.append(d1Around[0])
    d1Sample.append(d1Around[1])
    sx, sd1, se, sxi, _= sampleCubicHermiteCurves(xSample, d1Sample, elementsCountAround)
    xInner = xInner + sx[0:-1]
    d1Inner = smoothCubicHermiteDerivativesLoop(sx[0:-1], sd1[0:-1])

    # Pre-calculate node locations and derivatives on haustra inner cross-section
    elementsCountAroundSide = int(elementsCountAround/3)
    Ax = xInner[elementsCountAroundSide][0]
    Ay = xInner[elementsCountAroundSide][1]
    originRC = (Ax*Ax + Ay*Ay - haustraRadius*haustraRadius) / (2*(-Ax - haustraRadius))
    RC = haustraRadius - originRC

    if originRC > -Ax:
        startTheta = math.asin(Ay/RC)
        thetaRC = (math.pi - startTheta)*2
    else:
        startTheta = math.pi - math.asin(Ay/RC)
        thetaRC = math.asin(Ay/RC)*2
    thetaPerElementAround = thetaRC/(elementsCountAround/3)

    for n in range(elementsCountAroundSide + 1):
        theta = startTheta + thetaPerElementAround * n
        x = [RC*math.cos(theta) - originRC,
            RC*math.sin(theta),
            0.0]
        xHaustraSide.append(x)

    ang = [-2/3*math.pi, 0.0, 2/3*math.pi]

    for i in range(3):
        rotAng = ang[i]
        cosRotAng = math.cos(rotAng)
        sinRotAng = math.sin(rotAng)
        for n in range(elementsCountAroundSide):
            theta = startTheta + thetaPerElementAround * n
            cosTheta = math.cos(theta)
            sinTheta = math.sin(theta)
            x = [ (RC*cosTheta - originRC)*cosRotAng - RC*sinTheta*sinRotAng,
                  (RC*cosTheta - originRC)*sinRotAng + RC*sinTheta*cosRotAng,
                  0.0]
            xHaustraInner.append(x)
            dx_ds1 = [(-RC*sinTheta*cosRotAng - RC*cosTheta*sinRotAng)*thetaPerElementAround,
                  (-RC*sinTheta*sinRotAng + RC*cosTheta*cosRotAng)*thetaPerElementAround,
                  0.0]
            d1InnerHaustraRaw.append(dx_ds1)
    d1InnerHaustra = smoothCubicHermiteDerivativesLoop(xHaustraInner, d1InnerHaustraRaw)

    # Sample arclength of haustra
    for n1 in range(elementsCountAround):
        if n1%(elementsCountAround/3) > 0.0:
            v1 = [xInner[n1][0], xInner[n1][1], 0.0]
            startArcLength = foldFactor * haustraLength
            d1 = [ c*startArcLength for c in unitZ]
            v2 = [xHaustraInner[n1][0], xHaustraInner[n1][1], haustraLength/2]
            midArcLength = haustraCurvatureFactor * haustraLength
            d2 = [ c*midArcLength for c in unitZ]
            v3 = [xInner[n1][0], xInner[n1][1], haustraLength]
            d3 = [ c*startArcLength for c in unitZ]
        else:
            v1 = [xInner[n1][0], xInner[n1][1], 0.0]
            v2 = [xInner[n1][0], xInner[n1][1], haustraLength/2]
            v3 = [xInner[n1][0], xInner[n1][1], haustraLength]
            d3 = d2 = d1 = [c* haustraLength/3 for c in unitZ]
        nx = [v1, v2, v3]
        nd1 = [d1, d2, d3]
        sx, sd1, se, sxi, _  = sampleCubicHermiteCurves(nx, nd1, elementsCountAlongHaustrum)
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
            unitTangent = normalise(dx_ds2)
            # Interhaustra
            if n2 == 0 or n2 > elementsCountAlongHaustrum - 1:
                dx_ds1 = d1Inner[n1]
                unitdx_ds3 = crossproduct3(normalise(dx_ds1), unitTangent)
            else:
                # Haustra
                if elementsCountAlongHaustrum == 2:
                    unitdx_ds1 = normalise(d1InnerHaustra[n1])
                else:
                    if n1%(elementsCountAround/3) == 0: # intersection points
                        unitdx_ds1 = normalise(d1InnerHaustra[n1])
                    else: # points on clover
                        if elementsCountAlongHaustrum > 3:
                            if n2 < int(elementsCountAlongHaustrum/2): # first half of haustraLength
                                axisRot = crossproduct3(unitZ, unitTangent)
                            elif n2 > int(elementsCountAlongHaustrum/2): # second half of haustraLength
                                axisRot = crossproduct3(unitTangent, unitZ)
                        elif elementsCountAlongHaustrum == 3: # 3 elementsAlongHaustrum
                            axisRot = crossproduct3(unitTangent, unitZ)

                        rotFrame = rotationMatrixAboutAxis(axisRot, math.pi/2)
                        rotNormal = [rotFrame[j][0]*unitTangent[0] + rotFrame[j][1]*unitTangent[1] + rotFrame[j][2]*unitTangent[2] for j in range(3)]
                        unitdx_ds3 = normalise(rotNormal)
                        unitdx_ds1 = crossproduct3(unitTangent, unitdx_ds3)
                xAround.append(x)
                unitdx_ds1Around.append(unitdx_ds1)

        if n2 > 0 and n2 < elementsCountAlongHaustrum:
            dx_ds1InnerAroundList = []
            for n1 in range(elementsCountAround):
                v1 = xAround[n1]
                d1 = unitdx_ds1Around[n1]
                v2 = xAround[(n1+1)%elementsCountAround]
                d2 = unitdx_ds1Around[(n1+1)%elementsCountAround]
                arcLengthAround = computeCubicHermiteArcLength(v1, d1, v2, d2, True)
                dx_ds1 = [c*arcLengthAround for c in d1]
                dx_ds1InnerAroundList.append(dx_ds1)
            d1Smoothed = smoothCubicHermiteDerivativesLoop(xAround, dx_ds1InnerAroundList)
            dx_ds1InnerList = dx_ds1InnerList + d1Smoothed

    dx_ds1InnerList = dx_ds1InnerList + d1Inner

    for n in range(len(xInnerList)):
        dx_ds3 = crossproduct3(normalise(dx_ds1InnerList[n]), normalise(dx_ds2InnerList[n]))
        unitdx_ds3 = normalise(dx_ds3)
        dx_ds3InnerUnitList.append(unitdx_ds3)

    # Pre-calculate node locations and derivatives on outer boundary
    xOuter, d1Outer, d2Outer = getOuterCoordinatesAndDerivativesFromInner(xInnerList, dx_ds1InnerList, dx_ds2InnerList, dx_ds3InnerUnitList, wallThickness, elementsCountAlongHaustrum, elementsCountAround)

    return xOuter, d1Outer, d2Outer, unitZ

def getOuterCoordinatesAndDerivativesFromInner(xInner, d1Inner, d2Inner, d3Inner, wallThickness, elementsCountAlongHaustrum, elementsCountAround):
    """
    Generates coordinates and derivatives of outer surface from
    coordinates and derivatives of inner surface using wall thickness
    and normals.
    param xInner: Coordinates on inner surface
    param d1Inner: Derivatives on inner surface around haustra
    param d2Inner: Derivatives on inner surface along haustra
    param d3Inner: Derivatives on inner surface through wall
    param thickness: Thickness of wall
    param elementsCountAlongHaustrum: Number of elements along haustrum
    param elementsCountAround: Number of elements around haustra
    return xOuter: Coordinates on outer surface
    return nd1Outer: Derivatives on outer surface around haustra
    return nd2Outer: Derivatives on outer surface along haustra
    """
    xOuter = []
    nd1Outer = []
    nd2Outer = []

    for n2 in range(elementsCountAlongHaustrum + 1):
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            x = [xInner[n][i] + d3Inner[n][i]*wallThickness for i in range(3)]
            norm = d3Inner[n]
            # d1
            prevIdx = n-1 if (n1 != 0) else (n2+1)*elementsCountAround - 1
            nextIdx = n+1 if (n1 < elementsCountAround-1) else n2*elementsCountAround
            curvatureAround = 0.5*(
                getCubicHermiteCurvature(xInner[prevIdx], d1Inner[prevIdx], xInner[n], d1Inner[n], norm, 1.0) +
                getCubicHermiteCurvature(xInner[n], d1Inner[n], xInner[nextIdx], d1Inner[nextIdx], norm, 0.0))
            factor = 1.0 - curvatureAround*wallThickness
            nd1 = [ factor*c for c in d1Inner[n]]
            # d2
            if n2 > 0 and n2 < elementsCountAlongHaustrum:
                prevIdx = (n2-1)*elementsCountAround + n1
                nextIdx = (n2+1)*elementsCountAround + n1
            elif n2 == 0:
                prevIdx = (elementsCountAlongHaustrum-1)*elementsCountAround + n1
                nextIdx = (n2+1)*elementsCountAround + n1
            else:
                prevIdx = (n2-1)*elementsCountAround + n1
                nextIdx = elementsCountAround + n1
            curvatureAlong = 0.5*(
                getCubicHermiteCurvature(xInner[prevIdx], d2Inner[prevIdx], xInner[n], d2Inner[n], norm, 1.0)+
                getCubicHermiteCurvature(xInner[n], d2Inner[n], xInner[nextIdx], d2Inner[nextIdx], norm, 0.0))
            factor = 1.0 - curvatureAlong*wallThickness
            nd2 = [ factor*c for c in d2Inner[n]]

            xOuter.append(x)
            nd1Outer.append(nd1)
            nd2Outer.append(nd2)

    return xOuter, nd1Outer, nd2Outer

def rotationMatrixAboutAxis(rotAxis, theta):
    """
    Generate the rotation matrix for rotation about an axis.
    :param rotAxis: axis of rotation
    :param theta: angle of rotation
    :return: rotation matrix
    """
    cosTheta = math.cos(theta)
    sinTheta = math.sin(theta)
    C = 1 - cosTheta
    rotMatrix = ([[rotAxis[0]*rotAxis[0]*C + cosTheta, rotAxis[0]*rotAxis[1]*C - rotAxis[2]*sinTheta, rotAxis[0]*rotAxis[2]*C + rotAxis[1]*sinTheta],
        [rotAxis[1]*rotAxis[0]*C + rotAxis[2]*sinTheta, rotAxis[1]*rotAxis[1]*C + cosTheta, rotAxis[1]*rotAxis[2]*C - rotAxis[0]*sinTheta],
        [rotAxis[2]*rotAxis[0]*C - rotAxis[1]*sinTheta, rotAxis[2]*rotAxis[1]*C + rotAxis[0]*sinTheta, rotAxis[2]*rotAxis[2]*C + cosTheta]])

    return rotMatrix